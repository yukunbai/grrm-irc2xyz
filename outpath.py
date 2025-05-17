#!/usr/bin/env python3
"""
GRRM to XYZ Trajectory Converter

Reads a GRRM log (e.g., 1.log) and writes two multi-frame XYZ (.xyz) files:
 - forward IRC (TS → products)
 - backward IRC (TS → reactants)
Suitable for visualization in molecular viewers.
"""
import re
import argparse
import sys
import os


def parse_grrm(logfile):
    with open(logfile) as f:
        lines = f.read().splitlines()

    # Initial (TS) structure
    ts_atoms, ts_energy = [], None
    for i, line in enumerate(lines):
        if "INITIAL STRUCTURE" in line:
            for l in lines[i+1:]:
                parts = l.split()
                if not parts:
                    break
                if parts[0] == "ENERGY":
                    ts_energy = float(parts[2]); break
                if re.match(r'^[A-Z][a-z]?$', parts[0]):
                    ts_atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            break

    def collect(start_marker, end_marker):
        steps = []
        j = 0
        while j < len(lines) and start_marker not in lines[j]: j += 1
        j += 1
        while j < len(lines) and end_marker not in lines[j]:
            if lines[j].strip().startswith("# STEP"):
                j += 1
                atoms = []
                while j < len(lines) and "ENERGY" not in lines[j]:
                    tok = lines[j].split()
                    if re.match(r'^[A-Z][a-z]?$', tok[0]):
                        atoms.append((tok[0], float(tok[1]), float(tok[2]), float(tok[3])))
                    j += 1
                energy = float(lines[j].split()[2])
                steps.append((energy, atoms))
            j += 1
        return steps

    fwd = collect("IRC FOLLOWING (FORWARD)", "EQ EXIST WITHIN STEPSIZE")
    bwd = collect("IRC FOLLOWING (BACKWARD)", "Energy profile along IRC")
    return ts_atoms, ts_energy, fwd, bwd


def write_xyz(frames, outfile):
    """
    Write a sequence of frames to an XYZ file.
    """
    with open(outfile, 'w') as f:
        for energy, atoms, label in frames:
            f.write(f"{len(atoms)}\n")
            f.write(f"Energy={energy:.6f} Label={label}\n")
            for el, x, y, z in atoms:
                f.write(f"{el} {x:.6f} {y:.6f} {z:.6f}\n")
    print(f"Written XYZ trajectory: {outfile}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert GRRM log to separate forward/backward .xyz trajectories"
    )
    parser.add_argument('logfile', nargs='?', default='1.log', help='GRRM log file')
    parser.add_argument('fwd_out', nargs='?', default='irc_forward.xyz', help='Forward IRC .xyz file')
    parser.add_argument('bwd_out', nargs='?', default='irc_backward.xyz', help='Backward IRC .xyz file')
    args = parser.parse_args()

    if not os.path.exists(args.logfile):
        print(f"Error: {args.logfile} not found.")
        sys.exit(1)

    ts_atoms, ts_energy, fwd_steps, bwd_steps = parse_grrm(args.logfile)

    # Build frame lists
    forward_frames = [(ts_energy, ts_atoms, 'TS')] + [(e, atoms, 'FWD') for e, atoms in fwd_steps]
    backward_frames = [(ts_energy, ts_atoms, 'TS')] + [(e, atoms, 'BWD') for e, atoms in reversed(bwd_steps)]

    write_xyz(forward_frames, args.fwd_out)
    write_xyz(backward_frames, args.bwd_out)

if __name__ == '__main__':
    main()
