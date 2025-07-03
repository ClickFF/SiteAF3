#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generic hotspot structure generation script
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

from data.ligand_cutoff import get_hotspot_complex_struct, classify_chain


def create_parser():
    """Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='Generate a hotspot structure with hotspot residues'
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
        help='List chains in PDB file and their types, but do not generate hotspot structure'
    )
    
    return parser


def list_chains_in_pdb(pdb_file, verbose=False, return_info=False):
    """
    List chains in PDB file and their automatically recognized types
    
    Args:
        pdb_file: PDB file path
        verbose: Whether to show detailed information
        return_info: Whether to return chain information dictionary instead of just printing
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    chain_info_dict = {}

    if not return_info:
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
                chain_info_dict[chain.id] = {"type": "empty", "residue_count_str": "0", "source": "no_residues"}
                if not return_info and verbose:
                    print(f"{chain.id:<6} {'empty':<15} {'0':<8} {'no_residues'}")
                continue
                
            chain_type, classification_source = classify_chain(chain, verbose=False) # verbose for classify_chain kept false internally for now
            
            standard_residues = [r for r in residues if r.id[0] == ' ']
            hetero_residues = [r for r in residues if r.id[0] != ' ']
            
            residue_count_str = f"{len(standard_residues)}"
            if hetero_residues:
                residue_count_str += f"+{len(hetero_residues)}H"

            chain_info_dict[chain.id] = {
                "type": chain_type,
                "residue_count_str": residue_count_str,
                "source": classification_source
            }
            
            if not return_info:
                print(f"{chain.id:<6} {chain_type:<15} {residue_count_str:<8} {classification_source}")
                if verbose:
                    sample_residues = residues[:5]
                    residue_names = [r.resname for r in sample_residues]
                    if len(residues) > 5:
                        residue_names.append("...")
                    print(f"       Residue example: {', '.join(residue_names)}")
    
    if not return_info:
        print("-" * 60)
    
    if return_info:
        return chain_info_dict


def generate_hotspot_structure(input_pdb, output_pdb, receptor_type, ligand_type, 
                             hotspot_cutoff, receptor_chains=None, ligand_chains=None, 
                             verbose=False):
    """
    Generate hotspot structure
    
    Args:
        input_pdb: Input PDB file path
        output_pdb: Output PDB file path
        receptor_type: Receptor type
        ligand_type: Ligand type
        hotspot_cutoff: Hotspot cutoff distance
        receptor_chains: Receptor chain ID list
        ligand_chains: Ligand chain ID list
        verbose: Whether to show detailed output
    """
    if verbose:
        print(f"Input file: {input_pdb}")
        print(f"Output file: {output_pdb}")
        print(f"Receptor type: {receptor_type}")
        print(f"Ligand type: {ligand_type}")
        print(f"Hotspot cutoff distance: {hotspot_cutoff} Å")
        
        if receptor_chains:
            print(f"Specified receptor chains: {receptor_chains}")
        if ligand_chains:
            print(f"Specified ligand chains: {ligand_chains}")
    
    try:
        # Generate hotspot complex structure
        hotspot_struct = get_hotspot_complex_struct(
            infile=input_pdb,
            receptor_type=receptor_type,
            ligand_type=ligand_type,
            hotspot_cutoff=hotspot_cutoff,
            verbose=verbose,
            explicit_receptor_chain_ids=receptor_chains,
            explicit_ligand_chain_ids=ligand_chains
        )
        
        # Save structure
        io = PDB.PDBIO()
        io.set_structure(hotspot_struct)
        io.save(output_pdb)
        
        if verbose:
            print(f"Hotspot structure saved to: {output_pdb}")
            
            # Count residues in the output structure
            total_residues = 0
            total_chains = 0
            for model in hotspot_struct:
                for chain in model:
                    total_chains += 1
                    residues = list(chain.get_residues())
                    total_residues += len(residues)
                    if verbose:
                        print(f"  Chain {chain.id}: {len(residues)} residues")
            
            print(f"Output structure contains: {total_chains} chains, {total_residues} residues")
        
        return True
        
    except Exception as e:
        print(f"Error generating hotspot structure: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return False


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
        suffix = f"_hotspot_{args.receptor_type}_{args.ligand_type}"
        output_pdb = str(input_path.parent / f"{input_path.stem}{suffix}.pdb")
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_pdb)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate hotspot structure
    success = generate_hotspot_structure(
        input_pdb=args.input_pdb,
        output_pdb=output_pdb,
        receptor_type=args.receptor_type,
        ligand_type=args.ligand_type,
        hotspot_cutoff=args.hotspot_cutoff,
        receptor_chains=args.receptor_chains,
        ligand_chains=args.ligand_chains,
        verbose=args.verbose
    )
    
    if success:
        print(f"Hotspot structure generated successfully: {output_pdb}")
        return 0
    else:
        print("Hotspot structure generation failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())