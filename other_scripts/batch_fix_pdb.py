#!/usr/bin/env python
# batch_fix_pdb.py

import glob
import io # Required for StringIO
from pdbfixer import PDBFixer
from openmm.app import PDBFile

for pdb_path in glob.glob("./complex/*.pdb"):
    print(f"Processing {pdb_path} ...")
    try:
        fixer = PDBFixer(filename=pdb_path)

        # Initialize missingResidues attribute
        fixer.findMissingResidues()
        # Mark missing atoms
        fixer.findMissingAtoms()

        fixer.addMissingAtoms()

        # Do not call addMissingResidues(), do not fill missing residues
        # Do not call addMissingHydrogens(), do not fill missing hydrogens

        with open(pdb_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

        print(f"  → Successfully processed and overwrote {pdb_path}")
    except ValueError as ve:
        if "Misaligned residue name" in str(ve):
            print(f"  Skipping {pdb_path} due to misaligned residue name error: {ve}")
        else:
            print(f"  Skipping {pdb_path} due to a ValueError: {ve}")
    except Exception as e:
        print(f"  Skipping {pdb_path} due to an unexpected error: {e}")

