#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generic pocket structure generation script
Supports different receptor-ligand combinations:
- receptor: protein, nucleic
- ligand: protein, nucleic, small_molecule

Based on ligand_cutoff.py
"""

import os
import sys
import argparse
import pathlib
from Bio import PDB

script_dir = pathlib.Path(__file__).resolve().parent
src_dir = script_dir / "src"
sys.path.insert(0, str(src_dir))

from data.ligand_cutoff import get_motif_center_pos, classify_chain


def create_parser():
    """Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='Generate a pocket structure with pocket residues'
    )
    
    parser.add_argument(
        '--input_pdb',
        type=str,
        required=True,
        help='Input PDB file path'
    )
    
    parser.add_argument(
        '--output_pdb',
        type=str,
        help='Output PDB file path (if not specified, will be generated in the same directory as the input file)'
    )
    
    parser.add_argument(
        '--receptor_type',
        type=str,
        choices=['protein', 'nucleic'],
        default='protein',
        help='Receptor type: protein (protein) or nucleic (nucleic)'
    )
    
    parser.add_argument(
        '--ligand_type',
        type=str,
        choices=['protein', 'nucleic', 'small_molecule'],
        default='nucleic',
        help='Ligand type: protein (protein), nucleic (nucleic) or small_molecule (small molecule)'
    )
    
    parser.add_argument(
        '--hotspot_cutoff',
        type=float,
        default=8.0,
        help='Hotspot residue cutoff distance (Å), default 8.0'
    )
    
    parser.add_argument(
        '--pocket_cutoff',
        type=float,
        default=10.0,
        help='Pocket residue cutoff distance (Å), default 10.0'
    )
    
    parser.add_argument(
        '--receptor_chains',
        type=str,
        nargs='+',
        help='Explicitly specify receptor chain IDs, e.g.: --receptor_chains A B'
    )
    
    parser.add_argument(
        '--ligand_chains',
        type=str,
        nargs='+',
        help='Explicitly specify ligand chain IDs, e.g.: --ligand_chains X Y'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Show detailed output'
    )
    
    parser.add_argument(
        '--list_chains',
        action='store_true',
        help='List chains in PDB file and their types, but do not generate pocket structure'
    )
    
    parser.add_argument(
        '--save_center_info',
        action='store_true',
        help='Save ligand center position information to a text file'
    )
    
    return parser


def list_chains_in_pdb(pdb_file, verbose=False):
    """
    List chains in PDB file and their automatically recognized types
    
    Args:
        pdb_file: PDB file path
        verbose: Whether to show detailed information
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    print(f"\nPDB file: {pdb_file}")
    print("Chain information:")
    print("-" * 60)
    print(f"{'Chain ID':<6} {'Type':<15} {'Residue count':<8} {'Classification source'}")
    print("-" * 60)
    
    for model in structure:
        sorted_chains = sorted(model.get_chains(), key=lambda c: c.id)
        for chain in sorted_chains:
            residues = list(chain.get_residues())
            if not residues:
                continue
                
            chain_type, classification_source = classify_chain(chain, verbose=False)
            
            standard_residues = [r for r in residues if r.id[0] == ' ']
            hetero_residues = [r for r in residues if r.id[0] != ' ']
            
            residue_count = f"{len(standard_residues)}"
            if hetero_residues:
                residue_count += f"+{len(hetero_residues)}H"
            
            print(f"{chain.id:<6} {chain_type:<15} {residue_count:<8} {classification_source}")
            
            if verbose:
                sample_residues = residues[:5]
                residue_names = [r.resname for r in sample_residues]
                if len(residues) > 5:
                    residue_names.append("...")
                print(f"       Residue example: {', '.join(residue_names)}")
    
    print("-" * 60)


def generate_pocket_structure(input_pdb, output_pdb, receptor_type, ligand_type, 
                            hotspot_cutoff, pocket_cutoff, receptor_chains=None, 
                            ligand_chains=None, verbose=False, save_center_info=False):
    """
    Generate pocket structure
    
    Args:
        input_pdb: Input PDB file path
        output_pdb: Output PDB file path
        receptor_type: Receptor type
        ligand_type: Ligand type
        hotspot_cutoff: Hotspot cutoff distance
        pocket_cutoff: Pocket cutoff distance
        receptor_chains: Receptor chain ID list
        ligand_chains: Ligand chain ID list
        verbose: Whether to show detailed output
        save_center_info: Whether to save center position information
    """
    if verbose:
        print(f"Input file: {input_pdb}")
        print(f"Output file: {output_pdb}")
        print(f"Receptor type: {receptor_type}")
        print(f"Ligand type: {ligand_type}")
        print(f"Hotspot cutoff distance: {hotspot_cutoff} Å")
        print(f"Pocket cutoff distance: {pocket_cutoff} Å")
        
        if receptor_chains:
            print(f"Specified receptor chains: {receptor_chains}")
        if ligand_chains:
            print(f"Specified ligand chains: {ligand_chains}")
    
    try:
        # Generate pocket structure
        pocket_struct, center_pos, raw_seq_data = get_motif_center_pos(
            infile=input_pdb,
            receptor_type=receptor_type,
            ligand_type=ligand_type,
            hotspot_cutoff=hotspot_cutoff,
            pocket_cutoff=pocket_cutoff,
            verbose=verbose,
            explicit_receptor_chain_ids=receptor_chains,
            explicit_ligand_chain_ids=ligand_chains
        )
        
        # Save structure
        io = PDB.PDBIO()
        io.set_structure(pocket_struct)
        io.save(output_pdb)
        
        if verbose:
            print(f"Pocket structure saved to: {output_pdb}")
            print(f"Ligand center position: [{center_pos[0]:.3f}, {center_pos[1]:.3f}, {center_pos[2]:.3f}]")
            
            # Count residues in the output structure
            total_residues = 0
            total_chains = 0
            chain_info = {}
            for model in pocket_struct:
                for chain in model:
                    total_chains += 1
                    residues = list(chain.get_residues())
                    total_residues += len(residues)
                    chain_info[chain.id] = len(residues)
                    if verbose:
                        print(f"  Chain {chain.id}: {len(residues)} residues")
            
            print(f"Output structure contains: {total_chains} chains, {total_residues} residues")
            
            # Show chain classification information
            if raw_seq_data:
                receptor_chains_found = raw_seq_data.get("receptor_chains", [])
                ligand_chains_found = raw_seq_data.get("ligand_chains", [])
                print(f"Receptor chains: {receptor_chains_found}")
                print(f"Ligand chains: {ligand_chains_found}")
        
        # Save center position information
        if save_center_info:
            center_info_file = output_pdb.replace('.pdb', '_center_info.txt')
            with open(center_info_file, 'w') as f:
                f.write(f"# Pocket structure center position information\n")
                f.write(f"# Input file: {input_pdb}\n")
                f.write(f"# Receptor type: {receptor_type}\n")
                f.write(f"# Ligand type: {ligand_type}\n")
                f.write(f"# Hotspot cutoff distance: {hotspot_cutoff} Å\n")
                f.write(f"# Pocket cutoff distance: {pocket_cutoff} Å\n")
                f.write(f"\n")
                f.write(f"Ligand center position (x, y, z): {center_pos[0]:.6f}, {center_pos[1]:.6f}, {center_pos[2]:.6f}\n")
                
                if raw_seq_data:
                    f.write(f"\n# Chain information\n")
                    receptor_chains_found = raw_seq_data.get("receptor_chains", [])
                    ligand_chains_found = raw_seq_data.get("ligand_chains", [])
                    f.write(f"Receptor chains: {', '.join(receptor_chains_found)}\n")
                    f.write(f"Ligand chains: {', '.join(ligand_chains_found)}\n")
            
            if verbose:
                print(f"Center position information saved to: {center_info_file}")
        
        return True, center_pos, raw_seq_data
        
    except Exception as e:
        print(f"Error generating pocket structure: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return False, None, None


def main():
    parser = create_parser()
    args = parser.parse_args()
    
    # Check input file
    if not os.path.exists(args.input_pdb):
        print(f"Error: Input PDB file does not exist: {args.input_pdb}")
        return 1
    
    # If only listing chain information
    if args.list_chains:
        list_chains_in_pdb(args.input_pdb, args.verbose)
        return 0
    
    # Generate output file path
    if args.output_pdb:
        output_pdb = args.output_pdb
    else:
        # Generate output file in the same directory as the input file
        input_path = pathlib.Path(args.input_pdb)
        suffix = f"_pocket_{args.receptor_type}_{args.ligand_type}"
        output_pdb = str(input_path.parent / f"{input_path.stem}{suffix}.pdb")
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_pdb)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate pocket structure
    success, _, _ = generate_pocket_structure(
        input_pdb=args.input_pdb,
        output_pdb=output_pdb,
        receptor_type=args.receptor_type,
        ligand_type=args.ligand_type,
        hotspot_cutoff=args.hotspot_cutoff,
        pocket_cutoff=args.pocket_cutoff,
        receptor_chains=args.receptor_chains,
        ligand_chains=args.ligand_chains,
        verbose=args.verbose,
        save_center_info=args.save_center_info
    )
    
    if success:
        print(f"Pocket structure generated successfully: {output_pdb}")
        return 0
    else:
        print("Pocket structure generation failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
