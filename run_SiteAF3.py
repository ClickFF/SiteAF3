#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Main script to run the SiteAF3 pipeline.
Reads a JSON configuration file, generates embeddings for a protein-RNA complex,
and then predicts the 3D structure using a conditional diffusion model.
"""

import argparse
import json
import os
import pathlib
import sys
import pickle
import datetime
import time
import numpy as np


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR / "src"
sys.path.insert(0, str(SRC_DIR))

from src.embeddings.embed_cond import Embedder as EmbeddingGenerator
from src.diffusion.run_cond_Diff import ModelRunner as PredictionModelRunner, predict_structure_from_embedding
import jax
from alphafold3.model import model as af3_model

def _extract_smiles_from_ccd(ccd_file_path, ccd_code):
    """
    From CCD file, extract the SMILES string of the 
    specified CCD code and calculate the number of heavy atoms
    
    Args:
        ccd_file_path: Path to the CCD file
        ccd_code: The CCD code to find
        
    Returns:
        tuple: (SMILES string, number of heavy atoms)
        
    Raises:
        ValueError: If the specified CCD code or SMILES is not found
    """
    try:
        with open(ccd_file_path, 'r') as f:
            content = f.read()
        
        # Find the SMILES descriptor
        lines = content.split('\n')
        in_descriptor_loop = False
        smiles_found = None
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            if line.startswith('loop_') and i < len(lines) - 1:
                next_line = lines[i + 1].strip()
                if '_pdbx_chem_comp_descriptor' in next_line:
                    in_descriptor_loop = True
                    continue
            
            if in_descriptor_loop:
                if line.startswith('_pdbx_chem_comp_descriptor'):
                    continue
                
                # If it's an empty line or a new loop, end the current loop
                if not line or line.startswith('#') or line.startswith('loop_') or line.startswith('_'):
                    if not line.startswith('_pdbx_chem_comp_descriptor'):
                        in_descriptor_loop = False
                    continue
                
                # Parse the data line
                parts = line.split()
                if len(parts) >= 5:  # comp_id, type, program, version, descriptor
                    comp_id = parts[0]
                    descriptor_type = parts[1]
                    descriptor_value = ' '.join(parts[4:])  # SMILES may contain spaces (inside quotes)
                    
                    if descriptor_value.startswith('"') and descriptor_value.endswith('"'):
                        descriptor_value = descriptor_value[1:-1]
                    
                    # Check if it's the CCD code and SMILES type we're looking for
                    if comp_id == ccd_code and descriptor_type.upper() == 'SMILES':
                        smiles_found = descriptor_value
                        break
        
        if smiles_found:
            # Calculate the number of heavy atoms
            heavy_atoms_count = 0
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(smiles_found)
                if mol is not None:
                    heavy_atoms_count = mol.GetNumHeavyAtoms()
                else: # If the SMILES parsing fails, use a simple estimate
                    heavy_atoms_count = len([c for c in smiles_found if c.isupper() and c not in ['H']])
            except ImportError:
                # If RDKit is not installed, use a simple estimate: the number of uppercase letters (except H)
                heavy_atoms_count = len([c for c in smiles_found if c.isupper() and c not in ['H']])
            except Exception:
                raise ValueError(f"Error parsing SMILES string: {smiles_found}")
            return smiles_found, heavy_atoms_count
        else:
            raise ValueError(f"SMILES descriptor not found for CCD code {ccd_code}")
            
    except Exception as e:
        raise ValueError(f"Error parsing CCD file: {e}")

def create_arg_parser():
    """Creates and returns the argument parser."""
    parser = argparse.ArgumentParser(
        description="Run the SiteAF3 pipeline for protein-RNA complex structure prediction."
    )
    parser.add_argument(
        "--config_file",
        type=str,
        required=True,
        help="Path to the JSON configuration file for the case.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save the results.",
    )
    parser.add_argument(
        "--model_weights_dir",
        type=str,
        default="/work/hat170/aptamer/alphafold3/weight", # Default from user's other scripts
        help="Path to the AlphaFold3 model weights directory.",
    )
    parser.add_argument(
        "--db_dir",
        type=str,
        default="/work/hat170/aptamer/alphafold3/database/", # Default from user's other scripts
        help="Path to the MSA databases (used if use_af3_msa is enabled for embedding).",
    )
    parser.add_argument(
        "--receptor_type",
        type=str,
        choices=['protein', 'nucleic', 'auto'],
        default='auto',
        help="Receptor type: protein, nucleic, or auto (automatically detect). Default: auto"
    )
    parser.add_argument(
        "--ligand_type",
        type=str,
        choices=['protein', 'nucleic', 'small_molecule', 'auto'],
        default='auto',
        help="Ligand type: protein, nucleic, small_molecule, or auto (automatically detect). Default: auto"
    )
    parser.add_argument(
        "--hotspot_cutoff",
        type=float,
        default=8.0,
        help="Cutoff distance (Angstroms) for defining hotspot residues during embedding."
    )
    parser.add_argument(
        "--pocket_cutoff",
        type=float,
        default=10.0,
        help="Cutoff distance (Angstroms) for defining pocket residues during embedding."
    )
    parser.add_argument(
        "--use_af3_msa_for_embedding",
        action="store_true",
        default=False, # Default to False, enable explicitly
        help="Enable AlphaFold3's MSA generation pipeline during embedding. Can be combined with custom MSA options. Default is False.",
    )
    parser.add_argument(
        "--use_pocket_msa_for_embedding",
        action="store_true",
        default=False, # Changed default to False
        help="Enable pocket-based MSA generation for protein chains during embedding. Can be combined with AF3 MSA. Default is False.",
    )
    parser.add_argument(
        "--use_hotspot_msa_for_embedding",
        action="store_true",
        default=False, # Changed default to False
        help="Enable hotspot-based MSA generation for protein chains during embedding. Can be combined with AF3 MSA. Default is False.",
    )
    parser.add_argument(
        "--use_pocket_diffusion_for_prediction",
        action="store_true",
        default=True, # Defaulting to True for SiteAF3 typical use case
        help="Use pocket diffusion logic (fixed protein, diffusing ligand) for structure prediction. Default is True.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output."
    )
    return parser

def main(args):
    """Main execution function."""
    start_pipeline_time = time.time()

    # 1. Load and parse configuration file
    if not os.path.exists(args.config_file):
        print(f"Error: Configuration file not found at {args.config_file}")
        sys.exit(1)
    
    config_file_path = pathlib.Path(args.config_file)
    with open(config_file_path, "r") as f:
        try:
            config_data = json.load(f)
        except json.JSONDecodeError as e:
            print(f"Error decoding JSON config file {args.config_file}: {e}")
            sys.exit(1)

    if args.verbose:
        print(f"Successfully loaded configuration from: {args.config_file}")

    # Extract necessary parameters from config
    case_name = config_data.get("name", config_file_path.stem)
    model_seeds = config_data.get("modelSeeds", [42]) # Default to [42] if not specified
    if not isinstance(model_seeds, list) or not model_seeds:
        print("Warning: modelSeeds in config is not a valid list or is empty. Defaulting to [42].")
        model_seeds = [42]
    if args.verbose:
        print(f"Using model seeds: {model_seeds}")
    
    # Receptor PDB path
    if not config_data.get("receptor") or not config_data["receptor"][0].get("rec_struct_path"):
        print("Error: Receptor PDB path ('rec_struct_path') not found in config.")
        sys.exit(1)
    pdb_path_str = config_data["receptor"][0]["rec_struct_path"]
    
    # Parse receptor configuration
    receptor_config = config_data["receptor"][0]
    receptor_chain_ids = receptor_config.get("fixed_chain_id", [])
    if not isinstance(receptor_chain_ids, list):
        receptor_chain_ids = [str(receptor_chain_ids)] if receptor_chain_ids else []
    
    # Determine receptor type from config or default to protein
    receptor_type_from_config = 'protein'  # Default assumption for receptor
    if args.verbose and receptor_chain_ids:
        print(f"Parsed receptor: type={receptor_type_from_config}, chain IDs={receptor_chain_ids}")

    # Parse ligand configuration
    ligand_specs = {}
    ligand_chain_ids = []
    ligand_type_from_config = None
    
    if config_data.get("ligand"):
        for lig_entry in config_data["ligand"]:
            # Handle protein ligands
            if lig_entry.get("protein"):
                protein_info = lig_entry["protein"]
                chain_id = protein_info.get("id")
                sequence = protein_info.get("sequence")
                
                if chain_id and sequence:
                    ligand_specs[chain_id] = {
                        "sequence": sequence,
                        "type": "protein"
                    }
                    ligand_chain_ids.append(chain_id)
                    ligand_type_from_config = 'protein'
                    if args.verbose:
                        print(f"Parsed protein ligand: chain ID={chain_id}, sequence length={len(sequence)}")
                else:
                    print(f"Warning: Protein ligand missing 'id' or 'sequence'. Skipping.")
                continue
            
            # Handle RNA/DNA ligands
            if lig_entry.get("rna") or lig_entry.get("dna"):
                nucleic_type = "rna" if lig_entry.get("rna") else "dna"
                nucleic_info = lig_entry[nucleic_type]
                chain_id = nucleic_info.get("id")
                
                # Check for direct sequence or sequence file
                sequence = nucleic_info.get("sequence")
                if not sequence and nucleic_info.get("seq_path"):
                    seq_file_path_str = nucleic_info["seq_path"]
                    seq_file_path = pathlib.Path(seq_file_path_str)
                    if not seq_file_path.is_absolute():
                        seq_file_path = config_file_path.parent / seq_file_path
                    
                    if seq_file_path.exists():
                        with open(seq_file_path, "r") as seq_file:
                            lines = seq_file.readlines()
                            if lines[0].startswith(">"):  # FASTA format
                                sequence = "".join(line.strip() for line in lines[1:] if line.strip())
                            else:  # Raw sequence
                                sequence = "".join(line.strip() for line in lines if line.strip())
                    else:
                        print(f"Warning: Sequence file not found at {seq_file_path}. Skipping.")
                        continue
                
                if chain_id and sequence:
                    ligand_specs[chain_id] = {
                        "sequence": sequence,
                        "type": nucleic_type
                    }
                    ligand_chain_ids.append(chain_id)
                    ligand_type_from_config = 'nucleic'
                    if args.verbose:
                        print(f"Parsed {nucleic_type.upper()} ligand: chain ID={chain_id}, sequence length={len(sequence)}")
                else:
                    print(f"Warning: {nucleic_type.upper()} ligand missing 'id' or 'sequence'. Skipping.")
                continue
            
            # Handle small molecule ligands
            if lig_entry.get("small_molecule"):
                small_mol_info = lig_entry["small_molecule"]
                chain_id = small_mol_info.get("id")
                smiles = small_mol_info.get("smiles")
                ccd_codes = small_mol_info.get("ccdCodes")
                
                # If there's no SMILES but there are ccdCodes, try to get the SMILES from the CCD file
                if not smiles and ccd_codes and config_data.get("userCCDPath"):
                    if isinstance(ccd_codes, list) and len(ccd_codes) > 0:
                        ccd_code = ccd_codes[0]  # Use the first CCD code
                        user_ccd_path = config_data["userCCDPath"]
                        
                        # If the path is not absolute, relative to the config file directory
                        if not os.path.isabs(user_ccd_path):
                            user_ccd_path = os.path.join(config_file_path.parent, user_ccd_path)
                        
                        if os.path.exists(user_ccd_path):
                            try:
                                smiles, heavy_atoms_count = _extract_smiles_from_ccd(user_ccd_path, ccd_code)
                                if args.verbose:
                                    print(f"Extracted SMILES: {smiles}, number of heavy atoms: {heavy_atoms_count} from CCD file {user_ccd_path} for CCD code {ccd_code}")
                                
                                if not hasattr(args, '_ccd_heavy_atoms_info'):
                                    args._ccd_heavy_atoms_info = {}
                                args._ccd_heavy_atoms_info[ccd_code] = heavy_atoms_count
                                
                            except Exception as e:
                                print(f"Warning: Failed to extract SMILES from CCD file: {e}")
                        else:
                            print(f"Warning: CCD file not found: {user_ccd_path}")
                
                if chain_id and smiles:
                    ligand_specs[chain_id] = {
                        "smiles": smiles,
                        "type": "small_molecule"
                    }
                    if ccd_codes:
                        ligand_specs[chain_id]["ccd_codes"] = ccd_codes
                    
                    ligand_chain_ids.append(chain_id)
                    ligand_type_from_config = 'small_molecule'
                    if args.verbose:
                        ccd_info = f", CCD codes={ccd_codes}" if ccd_codes else ""
                        print(f"Parsed small molecule ligand: chain ID={chain_id}, SMILES={smiles}{ccd_info}")
                else:
                    missing_info = []
                    if not chain_id:
                        missing_info.append("'id'")
                    if not smiles:
                        missing_info.append("'smiles' or 'ccdCodes' with valid userCCDPath")
                    print(f"Warning: Small molecule ligand missing {' and '.join(missing_info)}. Skipping.")
                continue
            
            # Handle legacy format with explicit type
            if lig_entry.get("type"):
                lig_type = lig_entry["type"].lower()
                if lig_type in ["rna", "dna"] and lig_entry.get("seq_path"):
                    seq_file_path_str = lig_entry["seq_path"]
                    seq_file_path = pathlib.Path(seq_file_path_str)
                    if not seq_file_path.is_absolute():
                        seq_file_path = config_file_path.parent / seq_file_path
                    
                    if seq_file_path.exists():
                        with open(seq_file_path, "r") as seq_file:
                            lines = seq_file.readlines()
                            sequence = ""
                            if lines[0].startswith(">"):  # FASTA format
                                sequence = "".join(line.strip() for line in lines[1:] if line.strip())
                            else:  # Raw sequence
                                sequence = "".join(line.strip() for line in lines if line.strip())
                        
                        if sequence and lig_entry.get("chain_id"):
                            chain_id = lig_entry["chain_id"]
                            ligand_specs[chain_id] = {
                                "sequence": sequence,
                                "type": lig_type
                            }
                            ligand_chain_ids.append(chain_id)
                            ligand_type_from_config = 'nucleic'
                            if args.verbose:
                                print(f"Parsed {lig_type.upper()} ligand (legacy): chain ID={chain_id}, sequence length={len(sequence)}")
    
    if args.verbose:
        print(f"Final molecule type configuration:")
        print(f"  Receptor type: {receptor_type_from_config}, chain IDs: {receptor_chain_ids}")
        print(f"  Ligand type: {ligand_type_from_config}, chain IDs: {ligand_chain_ids}")
        print(f"  Ligand sequence specifications: {list(ligand_specs.keys())}")

    # Read predefined hotspot and pocket paths from config
    hotspot_pdb_path = None
    pocket_pdb_path = None
    if config_data.get("receptor") and config_data["receptor"][0].get("hotspot_path"):
        hotspot_path_str = config_data["receptor"][0]["hotspot_path"]
        hotspot_pdb_path = pathlib.Path(hotspot_path_str)
        if not hotspot_pdb_path.is_absolute():
            hotspot_pdb_path = config_file_path.parent / hotspot_pdb_path
        if hotspot_pdb_path.exists():
            if args.verbose:
                print(f"Using predefined hotspot file: {hotspot_pdb_path}")
        else:
            print(f"Warning: Hotspot PDB file not found at {hotspot_pdb_path}, will use dynamic calculation")
            hotspot_pdb_path = None
    
    if config_data.get("receptor") and config_data["receptor"][0].get("pocket_path"):
        pocket_path_str = config_data["receptor"][0]["pocket_path"]
        pocket_pdb_path = pathlib.Path(pocket_path_str)
        if not pocket_pdb_path.is_absolute():
            pocket_pdb_path = config_file_path.parent / pocket_pdb_path
        if pocket_pdb_path.exists():
            if args.verbose:
                print(f"Using predefined pocket file: {pocket_pdb_path}")
        else:
            print(f"Warning: Pocket PDB file not found at {pocket_pdb_path}, will use dynamic calculation")
            pocket_pdb_path = None

    # Resolve PDB path (if relative, assume relative to config file's directory)
    pdb_path = pathlib.Path(pdb_path_str)
    if not pdb_path.is_absolute():
        pdb_path = config_file_path.parent / pdb_path
    if not pdb_path.exists():
        print(f"Error: Receptor PDB file not found at {pdb_path}")
        sys.exit(1)

    # Auto-detect receptor and ligand types if needed
    receptor_type_final = args.receptor_type
    ligand_type_final = args.ligand_type
    
    if receptor_type_final == 'auto' or ligand_type_final == 'auto':
        if args.verbose:
            print("Auto-detecting molecule types...")
        
        # Import the classification function
        from data.ligand_cutoff import classify_chain
        from Bio import PDB
        
        # Parse PDB to analyze chains
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("temp", pdb_path)[0]
        
        detected_protein_chains = []
        detected_nucleic_chains = []
        detected_small_mol_chains = []
        
        for chain in structure.get_chains():
            chain_type, classification_source = classify_chain(chain, verbose=False)
            if chain_type == 'protein':
                detected_protein_chains.append(chain.id)
            elif chain_type == 'nucleic':
                detected_nucleic_chains.append(chain.id)
            elif chain_type == 'small_molecule':
                detected_small_mol_chains.append(chain.id)
            
            if args.verbose:
                print(f"链 {chain.id}: {chain_type} ({classification_source})")
        
        # Auto-detect receptor type
        if receptor_type_final == 'auto':
            if receptor_chain_ids:
                # Use explicitly specified protein chains
                receptor_type_final = 'protein'
                if args.verbose:
                    print(f"Auto-detected: receptor type is protein (based on explicitly specified chains: {receptor_chain_ids})")
            elif detected_protein_chains:
                receptor_type_final = 'protein'
                if args.verbose:
                    print(f"Auto-detected: receptor type is protein (detected protein chains: {detected_protein_chains})")
            elif detected_nucleic_chains:
                receptor_type_final = 'nucleic'
                if args.verbose:
                    print(f"Auto-detected: receptor type is nucleic (detected nucleic chains: {detected_nucleic_chains})")
            else:
                print("Warning: cannot auto-detect receptor type, defaulting to protein")
                receptor_type_final = 'protein'
        
        # Auto-detect ligand type
        if ligand_type_final == 'auto':
            if ligand_chain_ids:
                # Check the type from ligand_specs
                ligand_types_from_specs = [ligand_specs[lid]["type"] for lid in ligand_chain_ids]
                unique_types = list(set(ligand_types_from_specs))
                if len(unique_types) == 1:
                    if unique_types[0] in ['rna', 'dna']:
                        ligand_type_final = 'nucleic'
                    else:
                        ligand_type_final = unique_types[0]
                    if args.verbose:
                        print(f"Auto-detected: ligand type is {ligand_type_final} (based on definitions in the config file)")
                else:
                    print(f"Warning: detected multiple ligand types {unique_types}, using the first: {unique_types[0]}")
                    ligand_type_final = 'nucleic' if unique_types[0] in ['rna', 'dna'] else unique_types[0]
            elif detected_small_mol_chains:
                ligand_type_final = 'small_molecule'
                if args.verbose:
                    print(f"Auto-detected: ligand type is small molecule (detected small molecule chains: {detected_small_mol_chains})")
            elif detected_nucleic_chains:
                ligand_type_final = 'nucleic'
                if args.verbose:
                    print(f"Auto-detected: ligand type is nucleic (detected nucleic chains: {detected_nucleic_chains})")
            elif detected_protein_chains and receptor_type_final != 'protein':
                ligand_type_final = 'protein'
                if args.verbose:
                    print(f"Auto-detected: ligand type is protein (detected protein chains but receptor is not protein)")
            else:
                print("Warning: cannot auto-detect ligand type, defaulting to nucleic")
                ligand_type_final = 'nucleic'
    else:
        # Use types from configuration if available, otherwise use command line args
        if receptor_type_from_config:
            receptor_type_final = receptor_type_from_config
        if ligand_type_from_config:
            ligand_type_final = ligand_type_from_config
    
    if args.verbose:
        print(f"Final molecule types: receptor={receptor_type_final}, ligand={ligand_type_final}")

    # Create the main timestamped output directory for the whole case
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    case_specific_output_dir = pathlib.Path(args.output_dir) / f"{case_name}_{timestamp}"
    os.makedirs(case_specific_output_dir, exist_ok=True)
    if args.verbose:
        print(f"All results for this case will be saved under: {case_specific_output_dir}")

    # Initialize PredictionModelRunner once, as it's not seed-dependent for its parameters/config
    if args.verbose:
        print("Initializing prediction model runner (once)...")
    try:
        prediction_devices = jax.devices()
        if not prediction_devices:
            print("Error: No JAX devices found for prediction.")
            sys.exit(1)
        prediction_device = prediction_devices[0]
        if args.verbose:
            print(f"Using JAX device for prediction: {prediction_device}")
    except Exception as e:
        print(f"Error initializing JAX for prediction: {e}")
        sys.exit(1)

    prediction_model_config = af3_model.Model.Config()
    pred_runner = PredictionModelRunner(
        config=prediction_model_config, 
        device=prediction_device,
        model_dir=pathlib.Path(args.model_weights_dir),
        verbose=args.verbose,
        use_precomputed_embed=True, 
        use_pocket_diffusion=args.use_pocket_diffusion_for_prediction
    )

    # --- Iterate over each seed specified in the config --- 
    for seed_value in model_seeds:
        if args.verbose:
            print(f"\n--- Processing for SEED: {seed_value} ---")

        # Create a subdirectory for this specific seed's results
        seed_specific_run_dir = case_specific_output_dir / f"seed_{seed_value}"
        os.makedirs(seed_specific_run_dir, exist_ok=True)
        if args.verbose:
            print(f"Seed-specific results will be in: {seed_specific_run_dir}")

        # 3. Initialize Embedder for this seed (or can be outside if it doesn't hold seed-specific state)
        # For simplicity and clarity, let's keep it inside if any internal part of embedder might cache seed-based results, though unlikely for its current setup.
        if args.verbose:
            print(f"Initializing embedding generator for seed {seed_value}...")
        embedder = EmbeddingGenerator(
            weight_dir=pathlib.Path(args.model_weights_dir), 
            verbose=args.verbose
        )

        if hasattr(args, '_ccd_heavy_atoms_info'):
            embedder._current_ccd_heavy_atoms_info = args._ccd_heavy_atoms_info
            if args.verbose:
                print(f"Passing CCD heavy atoms info to the embedder: {args._ccd_heavy_atoms_info}")

        # 4. Generate embeddings for the current seed
        if args.verbose:
            print(f"Generating structural embeddings for {pdb_path} (seed: {seed_value})...")
            print(f"  Ligand specifications: {list(ligand_specs.keys())}")
            print(f"  Use AF3 MSA for embedding: {args.use_af3_msa_for_embedding}")
            print(f"  Use Pocket MSA for embedding: {args.use_pocket_msa_for_embedding}")
            print(f"  Use Hotspot MSA for embedding: {args.use_hotspot_msa_for_embedding}")

        embedding_start_time = time.time()
        embedding_data = None # Ensure it's defined in this scope
        try:
            # Select the appropriate method based on the detected molecule types
            if receptor_type_final in ['protein', 'nucleic'] and ligand_type_final in ['protein', 'nucleic', 'small_molecule']:
                if args.verbose:
                    print(f"Using generic embedding method (receptor: {receptor_type_final}, ligand: {ligand_type_final})")
                
                # Accurately map receptor and ligand chain IDs
                final_receptor_chain_ids = receptor_chain_ids if receptor_type_final == 'protein' else None
                final_ligand_chain_ids = ligand_chain_ids if ligand_type_final in ['protein', 'nucleic', 'small_molecule'] else None
                
                # For auto-detected cases, use the detected chains
                if not final_receptor_chain_ids and receptor_type_final == 'protein':
                    final_receptor_chain_ids = detected_protein_chains
                elif not final_receptor_chain_ids and receptor_type_final == 'nucleic':
                    final_receptor_chain_ids = detected_nucleic_chains
                
                if not final_ligand_chain_ids and ligand_type_final == 'nucleic':
                    final_ligand_chain_ids = detected_nucleic_chains
                elif not final_ligand_chain_ids and ligand_type_final == 'small_molecule':
                    final_ligand_chain_ids = detected_small_mol_chains
                elif not final_ligand_chain_ids and ligand_type_final == 'protein':
                    final_ligand_chain_ids = detected_protein_chains
                
                if args.verbose:
                    print(f"  Final receptor chain IDs: {final_receptor_chain_ids}")
                    print(f"  Final ligand chain IDs: {final_ligand_chain_ids}")
                
                # Use the generic struct_emb_af3 method
                embedding_data = embedder.struct_emb_af3(
                    pdb_file=str(pdb_path),
                    db_dir=args.db_dir if args.use_af3_msa_for_embedding else None,
                    msa_dir=None, 
                    use_af3_msa=args.use_af3_msa_for_embedding,
                    use_pocket_msa=args.use_pocket_msa_for_embedding,
                    use_hotspot_msa=args.use_hotspot_msa_for_embedding,
                    hotspot_cutoff=args.hotspot_cutoff,
                    pocket_cutoff=args.pocket_cutoff,
                    ligand_sequences_override=None,
                    ligand_specs=ligand_specs,
                    explicit_protein_chain_ids=final_receptor_chain_ids,
                    explicit_ligand_chain_ids=final_ligand_chain_ids,
                    rng_seed_override=seed_value,
                    predefined_hotspot_pdb=str(hotspot_pdb_path) if hotspot_pdb_path else None,
                    predefined_pocket_pdb=str(pocket_pdb_path) if pocket_pdb_path else None,
                    receptor_type=receptor_type_final,
                    ligand_type=ligand_type_final
                )
            else:
                # Use the original method as a fallback
                if args.verbose:
                    print(f"Using original embedding method (protein-nucleic complex)")
            embedding_data = embedder.struct_emb_af3(
                pdb_file=str(pdb_path),
                db_dir=args.db_dir if args.use_af3_msa_for_embedding else None,
                msa_dir=None, 
                use_af3_msa=args.use_af3_msa_for_embedding,
                use_pocket_msa=args.use_pocket_msa_for_embedding,
                use_hotspot_msa=args.use_hotspot_msa_for_embedding,
                hotspot_cutoff=args.hotspot_cutoff,
                pocket_cutoff=args.pocket_cutoff,
                    ligand_sequences_override=ligand_specs,
                ligand_specs=ligand_specs,
                    explicit_protein_chain_ids=receptor_chain_ids,
                    explicit_ligand_chain_ids=ligand_chain_ids,
                rng_seed_override=seed_value,
                predefined_hotspot_pdb=str(hotspot_pdb_path) if hotspot_pdb_path else None,
                    predefined_pocket_pdb=str(pocket_pdb_path) if pocket_pdb_path else None,
                    receptor_type=receptor_type_final,
                    ligand_type=ligand_type_final
            )
        except Exception as e:
            print(f"Error during embedding generation for seed {seed_value}: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
            continue # Skip to the next seed if embedding fails for this one
        
        embedding_time = time.time() - embedding_start_time
        if args.verbose:
            print(f"Embedding generation completed for seed {seed_value} in {embedding_time:.2f} seconds.")

        # 5. Predict structure using the generated embeddings for the current seed
        if args.verbose:
            print(f"Starting structure prediction for seed {seed_value}. Output directory: {seed_specific_run_dir}")
            print(f"  Use pocket diffusion for prediction: {args.use_pocket_diffusion_for_prediction}")

        prediction_start_time = time.time()
        try:
            predict_structure_from_embedding(
                embedding_data=embedding_data, 
                model_runner=pred_runner, # Use the runner initialized outside the loop
                output_dir=seed_specific_run_dir, 
                verbose=args.verbose,
                rng_seed=seed_value
            )
        except Exception as e:
            print(f"Error during structure prediction for seed {seed_value}: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
            continue # Skip to the next seed if prediction fails
        
        prediction_time = time.time() - prediction_start_time
        if args.verbose:
            print(f"Structure prediction completed for seed {seed_value} in {prediction_time:.2f} seconds.")
    
    # End of seed loop

    total_pipeline_time = time.time() - start_pipeline_time
    print(f"\nSiteAF3 pipeline finished for {case_name}.")
    print(f"Total execution time: {total_pipeline_time:.2f} seconds.")
    print(f"All outputs saved in directory: {case_specific_output_dir}")

    return 0

if __name__ == "__main__":
    parser = create_arg_parser()
    args = parser.parse_args()
    
    if not pathlib.Path(args.model_weights_dir).exists():
        print(f"Error: Model weights directory not found at {args.model_weights_dir}")
        sys.exit(1)
    
    sys.exit(main(args))
