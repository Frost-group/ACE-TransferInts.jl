import os
import shutil
import subprocess
import tempfile

import ase
from ase.io import read, write

from dipro_full_parser2 import parse_all, run_xtb

def process_dimer(dimer, tmpdir):
    xyz_path = os.path.join(tmpdir, "dimer.xyz")
    write(xyz_path, dimer)

    # Run xtb with dipro
    run_xtb(xyz_path, output_file=os.path.join(tmpdir, "dipro_full.out"))

    # Parse output
    results = parse_all(os.path.join(tmpdir, "dipro_full.out"))

    # Extract HOMO-HOMO J and J_eff (hole transport)
    homo_homo = results.get("HOMO–HOMO")
    if homo_homo is None:
        raise ValueError("HOMO-HOMO block not found in DIPRO output.")

    J = homo_homo["J"]  # Coupling constant
    J_eff = homo_homo["J_eff_calc"]  # Effective coupling

    return J, J_eff

def main(xyz_infile, extxyz_infile, extxyz_outfile):
    # Load all dimers
    dimers = read(xyz_infile, index=":")

    # Load extxyz file
    frames = read(extxyz_infile, index=":")

    if len(dimers) != len(frames):
        raise ValueError("Number of dimers and frames in extxyz do not match.")

    with tempfile.TemporaryDirectory() as tmpdir:
        for i, dimer in enumerate(dimers):
            print(f"Processing dimer {i+1}/{len(dimers)}...")
            J, J_eff = process_dimer(dimer, tmpdir)

            # Attach results to corresponding frame
            frames[i].info["energy"] = J
            frames[i].info["energy_eff"] = J_eff

    # Write updated frames to new extxyz file
    write(extxyz_outfile, frames)

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 4:
        print("Usage: python script.py dimer.xyz input.extxyz output.extxyz")
        sys.exit(1)

    xyz_infile = sys.argv[1]
    extxyz_infile = sys.argv[2]
    extxyz_outfile = sys.argv[3]

    main(xyz_infile, extxyz_infile, extxyz_outfile)
