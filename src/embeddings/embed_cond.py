"""
This module is used to process AlphaFold3 embeddings, and can be run independently.
Use AlphaFold3 to generate embeddings.
"""
import os
import argparse
import pathlib
import copy
import pickle
from tqdm import tqdm
import numpy as np
import pandas as pd
import json
import time
import datetime
import mdtraj as md
from Bio import PDB
from pathlib import Path
import sys
from Bio import SeqIO # Added for FASTA parsing

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Import AlphaFold3 related modules
import jax
import jax.numpy as jnp
import haiku as hk
from alphafold3 import structure
from alphafold3.common import folding_input
from alphafold3.model import model
from alphafold3.model import params as af3_params
from alphafold3.model import features as af3_features
from alphafold3.model.components import utils as af3_utils
from alphafold3.model.network import featurization
from alphafold3.data import featurisation
from alphafold3.constants import chemical_components
from alphafold3.constants import atom_types

# Import data processing modules
from data import utils as du
from data import errors
from data import parsers
from data import residue_constants
from data.ligand_cutoff import (get_seq, resid_unique, 
                              get_motif_center_pos, classify_chain, get_hotspot_complex_struct,
                              get_representative_atom, is_small_molecule)

# Import MSA masking functions
from embeddings.masking_af3_msa import remove_duplicate_msa_sequences
try:
    from embeddings.af3_msa_masking import apply_pocket_mask_to_msa_af3
    AF3_MSA_MASKING_AVAILABLE = True
except ImportError:
    from embeddings.masking_af3_msa import apply_pocket_mask_to_msa
    AF3_MSA_MASKING_AVAILABLE = False
    print("Warning: AF3 MSA masking not available, using simple method")

MAX_ATOMS_PER_TOKEN = 24

# Helper function to parse PDB atoms (simplified)
def parse_pdb_atoms(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_file)
    atom_coords = {}
    for model in structure:
        # Sort chains by chain.id before iterating
        sorted_chains = sorted(model.get_chains(), key=lambda chain: chain.id)
        for chain in sorted_chains:
            for residue in chain:
                res_id_tuple = (chain.id, residue.id[1])
                if residue.id[0] == ' ': # Standard residues
                    atom_coords[res_id_tuple] = {}
                    for atom in residue:
                        # Normalize atom name
                        atom_name = atom.get_name().strip()
                        atom_coords[res_id_tuple][atom_name] = atom.get_coord()
    return atom_coords

class Embedder():
    """
    AlphaFold3 Embedder class
    """
    def __init__(self, 
                 weight_dir=None, 
                 device=None, 
                 verbose=False):
        """
        Initialize Embedder
        
        Args:
            weight_dir: AlphaFold3 weight directory
            device: Device to use
            verbose: Whether to display detailed output
        """
        self.verbose = verbose
        
        # Set weight directory
        if weight_dir is None:
            raise ValueError("AlphaFold3 weight directory is required")
        else:
            self.weight_dir = Path(weight_dir)
            if not self.weight_dir.exists():
                raise FileNotFoundError(f"AlphaFold3 weight directory does not exist: {self.weight_dir}")
        
        # Set device
        if device is None:
            devices = jax.devices()
            if devices:
                self.device = devices[0]
            else:
                raise RuntimeError("No devices found")
        else:
            self.device = device
            
        # Initialize model parameters
        try:
            self.model_params = af3_params.get_model_haiku_params(
                model_dir=self.weight_dir,
            )
            if self.verbose:
                print(f"Successfully loaded AlphaFold3 model parameters from {self.weight_dir}")
        except Exception as e:
            raise RuntimeError(f"Failed to load AlphaFold3 model parameters: {e}")
            
        # Initialize model configuration
        self.model_config = model.Model.Config()
        self.model_config.return_embeddings = True  # Ensure return embeddings
        
        # Create forward propagation function
        self._init_forward_fn()
    
    def _init_forward_fn(self):
        """Initialize model forward propagation function"""
        @hk.transform
        def forward_fn(batch):
            with af3_utils.bfloat16_context():
                af3_model = model.Model(self.model_config)
                return af3_model(batch)
                
        self._forward_fn_raw = forward_fn
        
    def _forward_fn(self, key, features_batch):
        """
        Call model forward propagation to get embeddings
        
        Args:
            key: Random seed
            features_batch: Features batch
            
        Returns:
            Model output result
        """
        batch = jax.device_put(
            jax.tree_util.tree_map(
                jnp.asarray, af3_utils.remove_invalidly_typed_feats(features_batch)
            ),
            self.device,
        )
        result = jax.jit(self._forward_fn_raw.apply, device=self.device)(
            self.model_params, key, batch
        )
        
        return result
    
    def _is_nucleic(self, residue):
        """
        Check if residue is a nucleic acid residue
        
        Args:
            residue: Biopython residue object
            
        Returns:
            bool: True if residue is a nucleic acid residue
        """
        if residue.id[0] != ' ':  # Non-standard residue
            return False
            
        # Determine based on residue name
        resname = residue.resname
        return resname in residue_constants.resname_to_nucleic_acid_1letter
    
    def _get_representative_atom(self, residue):
        """
        Get the representative atom coordinates of the residue, for center position calculation
        
        Args:
            residue: Biopython residue object
            
        Returns:
            np.ndarray: Representative atom coordinates
        """
        # For nucleic acids, use P atom first, then C4', and finally any atom
        if self._is_nucleic(residue):
            atom_names = ['P', "C4'", 'C4*']
        else:  # For proteins, use CA atom first
            atom_names = ['CA']
            
        # Try to find the atom based on priority
        for atom_name in atom_names:
            if atom_name in residue:
                return residue[atom_name].get_coord()
                
        # If no atom is found, return the coordinates of the first atom
        for atom in residue:
            return atom.get_coord()
            
        # If residue has no atoms, raise an error
        raise ValueError(f"Residue {residue.id} has no atoms")
 
    def _get_protein_mask(self, fold_input):
        """
        Generate a protein mask, shape is (token_num,)
        Protein residues are 1, others are 0
        """
        protein_mask = []
        protein_indices = [] # Store indices of protein residues
        other_indices = [] # Store indices of other residues
        current_idx = 0
        for chain in fold_input.chains:
            is_protein_chain = isinstance(chain, folding_input.ProteinChain)
            
            # Determine the number of tokens in the chain
            if isinstance(chain, folding_input.Ligand):
                # Small molecule: number of tokens is the number of atoms
                if hasattr(chain, 'smiles') and chain.smiles:
                    # Use SMILES to calculate number of atoms
                    try:
                        from rdkit import Chem
                        mol = Chem.MolFromSmiles(chain.smiles)
                        if mol is not None:
                            chain_length = mol.GetNumAtoms()
                        else:
                            raise ValueError("Failed to process small molecule from SMILES")
                    except ImportError:
                        raise ImportError("RDKit is required for small molecule processing")
                    except Exception:
                        raise ValueError("Failed to process small molecule from SMILES")
                elif hasattr(chain, 'ccd_ids') and chain.ccd_ids:
                    # For CCD-defined small molecules, use actual atoms in PDB instead of theoretical heavy atoms
                    # This avoids the mismatch between CCD heavy atom count and actual PDB atoms
                    actual_atoms_in_pdb = 0
                    
                    # Count actual atoms in PDB for this chain
                    for model in fold_input.chains:
                        if hasattr(model, 'id') and model.id == chain.id:
                            # This is the current chain we're checking
                            continue
                    
                    # Use CCD heavy atoms as fallback if PDB parsing fails
                    if hasattr(self, '_current_ccd_heavy_atoms_info'):
                        ccd_heavy_atoms_info = getattr(self, '_current_ccd_heavy_atoms_info', {})
                        total_heavy_atoms = 0
                        for ccd_id in chain.ccd_ids:
                            if ccd_id in ccd_heavy_atoms_info:
                                total_heavy_atoms += ccd_heavy_atoms_info[ccd_id]
                        chain_length = total_heavy_atoms
                        if self.verbose:
                            print(f"Using CCD heavy atoms count ({total_heavy_atoms}) for small molecule chain {chain.id}")
                    else:
                        raise ValueError(f"No cached CCD heavy atoms info, cannot calculate number of heavy atoms for small molecule {chain.id}")
                else:
                    raise ValueError(f"Small molecule chain {chain.id} neither has SMILES nor CCD codes")
            else:
                # Protein/nucleic acid: number of tokens is the number of residues
                chain_length = len(chain.sequence)
            
            for _ in range(chain_length):
                if is_protein_chain:
                    protein_mask.append(1)
                    protein_indices.append(current_idx)
                else:
                    protein_mask.append(0)
                    other_indices.append(current_idx)
                current_idx += 1
        protein_mask = np.array(protein_mask)
        if self.verbose:
            num_protein_residues = np.sum(protein_mask)
            print(f"protein_mask shape: {protein_mask.shape}, number of protein residues: {num_protein_residues}")
        return protein_mask, protein_indices, other_indices
    
    def _get_pocket_mask(self, pdb_file, receptor_type='protein', ligand_type='nucleic', hot_spot_cutoff=8.0, pocket_cutoff=10.0, 
                        explicit_receptor_chain_ids: list = None, explicit_ligand_chain_ids: list = None, 
                        predefined_pocket_pdb: str = None):
        """
        Universal pocket mask generation function, supports different types of receptor and ligand combinations
        
        Args:
            pdb_file: PDB file path
            receptor_type: Receptor type ('protein', 'nucleic')
            ligand_type: Ligand type ('protein', 'nucleic', 'small_molecule')
            hot_spot_cutoff: Hotspot residue cutoff distance (Å)
            pocket_cutoff: Pocket residue cutoff distance (Å)
            explicit_receptor_chain_ids: Optional, explicitly specify the receptor chain ID list
            explicit_ligand_chain_ids: Optional, explicitly specify the ligand chain ID list
            predefined_pocket_pdb: Optional, predefined pocket PDB file path
            
        Returns:
            tuple: (pocket_mask, center_pos, raw_seq_data)
        """
        if predefined_pocket_pdb:
            if self.verbose:
                print(f"Using predefined pocket PDB file: {predefined_pocket_pdb}")
            # Set current chain information for _get_pocket_mask_from_pdb_file
            self._current_explicit_receptor_chain_ids = explicit_receptor_chain_ids
            self._current_explicit_ligand_chain_ids = explicit_ligand_chain_ids
            self._current_ligand_type = ligand_type
            return self._get_pocket_mask_from_pdb_file(pdb_file, predefined_pocket_pdb)
        
        try:
            if self.verbose:
                print(f"Using generic method to calculate pocket mask")
                print(f"  Receptor type: {receptor_type}")
                print(f"  Ligand type: {ligand_type}")
                print(f"  Hotspot cutoff: {hot_spot_cutoff}Å")
                print(f"  Pocket cutoff: {pocket_cutoff}Å")
            
            processed_struct, center_pos, raw_seq_data = get_motif_center_pos(
                infile=pdb_file,
                receptor_type=receptor_type,
                ligand_type=ligand_type,
                hotspot_cutoff=hot_spot_cutoff,
                pocket_cutoff=pocket_cutoff,
                verbose=self.verbose,
                explicit_receptor_chain_ids=explicit_receptor_chain_ids,
                explicit_ligand_chain_ids=explicit_ligand_chain_ids
            )
            
            if self.verbose:
                print(f"Pocket calculation completed, center position: [{center_pos[0]:.3f}, {center_pos[1]:.3f}, {center_pos[2]:.3f}]")
            
            # Generate mask - based on receptor chain residues in the original PDB structure
            parser = PDB.PDBParser(QUIET=True)
            original_structure = parser.get_structure("original", pdb_file)[0]
            processed_structure = processed_struct[0] if hasattr(processed_struct, '__getitem__') else processed_struct
            
            # Get receptor chain residues in the original PDB structure (sorted by chain ID)
            # Ensure only chains with residue index mapping are included
            original_protein_residues = []
            for chain in sorted(original_structure.get_chains(), key=lambda c: c.id):
                # Only consider chains specified in explicit_receptor_chain_ids
                if explicit_receptor_chain_ids and chain.id in explicit_receptor_chain_ids:
                    for residue in chain.get_residues():
                        if residue.id[0] == ' ':  # Standard residue
                            original_protein_residues.append(residue)
            
            processed_residue_ids = set()
            for chain in processed_structure.get_chains():
                for residue in chain.get_residues():
                    if residue.id[0] == ' ':
                        residue_id = f"{residue.id[1]}_{chain.id}"
                        processed_residue_ids.add(residue_id)
            
            # Create mask
            pocket_mask = []
            for residue in original_protein_residues:
                residue_id = f"{residue.id[1]}_{residue.parent.id}"
                pocket_mask.append(residue_id in processed_residue_ids)
            
            pocket_mask = np.array(pocket_mask, dtype=bool)
            
            if self.verbose:
                total_residues = len(pocket_mask)
                pocket_count = np.sum(pocket_mask)
                print(f"Generated pocket mask: {pocket_count}/{total_residues} protein residues in pocket ({pocket_count/total_residues*100:.1f}%)")
            
            return pocket_mask, center_pos, raw_seq_data
            
        except Exception as e:
            raise ValueError(f"Failed to generate pocket mask: {e}")

    def _get_hotspot_mask(self, pdb_file, receptor_type='protein', ligand_type='nucleic', hotspot_cutoff=8.0, 
                         explicit_receptor_chain_ids: list = None, explicit_ligand_chain_ids: list = None, 
                         predefined_hotspot_pdb: str = None):
        """
        Universal hotspot mask generation function, supports different types of receptor and ligand combinations
        
        Args:
            pdb_file: PDB file path
            receptor_type: Receptor type ('protein', 'nucleic')
            ligand_type: Ligand type ('protein', 'nucleic', 'small_molecule')
            hotspot_cutoff: Hotspot residue cutoff distance (Å)
            explicit_receptor_chain_ids: Optional, explicitly specify the receptor chain ID list
            explicit_ligand_chain_ids: Optional, explicitly specify the ligand chain ID list
            predefined_hotspot_pdb: Optional, predefined hotspot PDB file path
            
        Returns:
            np.ndarray: Boolean hotspot mask array
        """
        if predefined_hotspot_pdb:
            if self.verbose:
                print(f"Using predefined hotspot PDB file: {predefined_hotspot_pdb}")
            # Set current chain information for _get_hotspot_mask_from_pdb_file
            self._current_explicit_receptor_chain_ids = explicit_receptor_chain_ids
            self._current_explicit_ligand_chain_ids = explicit_ligand_chain_ids
            self._current_ligand_type = ligand_type
            return self._get_hotspot_mask_from_pdb_file(pdb_file, predefined_hotspot_pdb)
        
        try:
            if self.verbose:
                print(f"Using generic method to calculate hotspot mask")
                print(f"  Receptor type: {receptor_type}")
                print(f"  Ligand type: {ligand_type}")
                print(f"  Hotspot cutoff: {hotspot_cutoff}Å")
            
            hotspot_struct = get_hotspot_complex_struct(
                infile=pdb_file,
                receptor_type=receptor_type,
                ligand_type=ligand_type,
                hotspot_cutoff=hotspot_cutoff,
                verbose=self.verbose,
                explicit_receptor_chain_ids=explicit_receptor_chain_ids,
                explicit_ligand_chain_ids=explicit_ligand_chain_ids
            )
            
            if self.verbose:
                print("hotspot calculation completed")
            
            # Generate mask - based on receptor chain residues in the original PDB structure
            parser = PDB.PDBParser(QUIET=True)
            original_structure = parser.get_structure("original", pdb_file)[0]
            hotspot_structure = hotspot_struct[0] if hasattr(hotspot_struct, '__getitem__') else hotspot_struct
            
            # Get receptor chain residues in the original PDB structure (sorted by chain ID)
            # Ensure only chains with residue index mapping are included
            original_protein_residues = []
            for chain in sorted(original_structure.get_chains(), key=lambda c: c.id):
                # Only consider chains specified in explicit_receptor_chain_ids
                if explicit_receptor_chain_ids and chain.id in explicit_receptor_chain_ids:
                    for residue in chain.get_residues():
                        if residue.id[0] == ' ':  # Standard residue
                            original_protein_residues.append(residue)
            
            hotspot_residue_ids = set()
            for chain in hotspot_structure.get_chains():
                for residue in chain.get_residues():
                    if residue.id[0] == ' ':
                        residue_id = f"{residue.id[1]}_{chain.id}"
                        hotspot_residue_ids.add(residue_id)
            
            hotspot_mask = []
            for residue in original_protein_residues:
                residue_id = f"{residue.id[1]}_{residue.parent.id}"
                hotspot_mask.append(residue_id in hotspot_residue_ids)
            
            hotspot_mask = np.array(hotspot_mask, dtype=bool)
            
            if self.verbose:
                total_residues = len(hotspot_mask)
                hotspot_count = np.sum(hotspot_mask)
                print(f"Generated hotspot mask: {hotspot_count}/{total_residues} protein residues are hotspots ({hotspot_count/total_residues*100:.1f}%)")
            
            return hotspot_mask
            
        except Exception as e:
            raise ValueError(f"Failed to generate hotspot mask: {e}")

    # Backward compatibility wrapper functions
    def _get_pocket_mask_by_type(self, pdb_file, receptor_type, ligand_type, hot_spot_cutoff=8.0, pocket_cutoff=10.0, 
                                explicit_receptor_chain_ids: list = None, explicit_ligand_chain_ids: list = None, 
                                predefined_pocket_pdb: str = None):
        return self._get_pocket_mask(pdb_file, receptor_type, ligand_type, hot_spot_cutoff, pocket_cutoff,
                                   explicit_receptor_chain_ids, explicit_ligand_chain_ids, predefined_pocket_pdb)

    def _get_hotspot_mask_by_type(self, pdb_file, receptor_type, ligand_type, hotspot_cutoff=8.0, 
                                 explicit_receptor_chain_ids: list = None, explicit_ligand_chain_ids: list = None, 
                                 predefined_hotspot_pdb: str = None):
        return self._get_hotspot_mask(pdb_file, receptor_type, ligand_type, hotspot_cutoff,
                                    explicit_receptor_chain_ids, explicit_ligand_chain_ids, predefined_hotspot_pdb)
    
    def _get_pocket_mask_legacy(self, pdb_file, hot_spot_cutoff=8.0, pocket_cutoff=10.0,
                               explicit_receptor_chain_ids=None, explicit_ligand_chain_ids=None,
                               predefined_pocket_pdb=None):
        return self._get_pocket_mask(
            pdb_file, 'protein', 'nucleic', hot_spot_cutoff, pocket_cutoff,
            explicit_receptor_chain_ids, explicit_ligand_chain_ids, predefined_pocket_pdb
        )
    
    def _get_hotspot_mask_legacy(self, pdb_file, hotspot_cutoff=8.0,
                                explicit_receptor_chain_ids=None, explicit_ligand_chain_ids=None,
                                predefined_hotspot_pdb=None):
        return self._get_hotspot_mask(
            pdb_file, 'protein', 'nucleic', hotspot_cutoff,
            explicit_receptor_chain_ids, explicit_ligand_chain_ids, predefined_hotspot_pdb
        )
    
    def _get_pocket_mask_from_pdb_file(self, original_pdb_file, pocket_pdb_file):
        """
        Calculate mask from predefined pocket PDB file
        
        Args:
            original_pdb_file: Original PDB file path
            pocket_pdb_file: Pocket PDB file path
            
        Returns:
            pocket_mask: Pocket mask array
            center_pos: Center position of the complex
            raw_seq_data: Sequence data dictionary (for compatibility)
        """
        try:
            # Parse original PDB file and pocket PDB file
            parser = PDB.PDBParser(QUIET=True)
            original_structure = parser.get_structure('original', original_pdb_file)[0]
            pocket_structure = parser.get_structure('pocket', pocket_pdb_file)[0]
            
            # Get specified receptor chain IDs
            receptor_chain_ids = getattr(self, '_current_explicit_receptor_chain_ids', [])
            if not receptor_chain_ids:
                if self.verbose:
                    print("Warning: No receptor chain IDs specified, will process all protein chains")
                # If no specified, process all protein chains
                receptor_chain_ids = []
                for chain in sorted(original_structure.get_chains(), key=lambda c: c.id):
                    from data.ligand_cutoff import classify_chain
                    chain_type, _ = classify_chain(chain, verbose=False)
                    if chain_type == 'protein':
                        receptor_chain_ids.append(chain.id)
            
            if self.verbose:
                print(f"Processing receptor chains: {receptor_chain_ids}")
            
            # Get all pocket residues by chain ID from pocket PDB
            pocket_residues_by_chain = {}
            for chain in pocket_structure.get_chains():
                chain_id = chain.id
                if chain_id not in pocket_residues_by_chain:
                    pocket_residues_by_chain[chain_id] = []
                for residue in chain.get_residues():
                    if residue.id[0] == ' ':  
                        pocket_residues_by_chain[chain_id].append(residue.id[1])
            
            if self.verbose:
                for chain_id, res_ids in pocket_residues_by_chain.items():
                    print(f"Number of residues in chain {chain_id} in pocket PDB: {len(res_ids)}")
            
            # Build mask based on chain order in original PDB
            all_receptor_residues = []
            pocket_mask_values = []
            
            for chain in sorted(original_structure.get_chains(), key=lambda c: c.id):
                chain_id = chain.id
                if chain_id in receptor_chain_ids:
                    chain_residues = [res for res in chain.get_residues() if res.id[0] == ' ']
                    pocket_res_ids = set(pocket_residues_by_chain.get(chain_id, []))
                    
                    for residue in chain_residues:
                        res_id = residue.id[1]
                        all_receptor_residues.append((chain_id, res_id))
                        # Check if the residue is in the pocket
                        pocket_mask_values.append(1 if res_id in pocket_res_ids else 0)
            
            pocket_mask = np.array(pocket_mask_values, dtype=float)
            
            if self.verbose:
                print(f"Total number of receptor residues in original PDB: {len(all_receptor_residues)}")
                pocket_count = int(np.sum(pocket_mask))
                print(f"Read {pocket_count}/{len(all_receptor_residues)} receptor residues from predefined pocket PDB")
            
            # Calculate center position (using the geometric center of pocket residues)
            center_coords = []
            for chain in pocket_structure.get_chains():
                for residue in chain.get_residues():
                    if residue.id[0] == ' ':
                        try:
                            coord = self._get_representative_atom(residue)
                            center_coords.append(coord)
                        except ValueError:
                            continue
            
            if center_coords:
                center_pos = np.mean(center_coords, axis=0)
                print(f"Pocket center position: {center_pos}")
            else:
                center_pos = np.zeros(3)
                print("Warning: No pocket residues found, using zero center position")
            
            raw_seq_data = {
                'receptor_chains': set(receptor_chain_ids),
                'ligand_chains': set()
            }
                
            return pocket_mask, center_pos, raw_seq_data
            
        except Exception as e:
            raise ValueError(f"Failed to generate pocket mask from predefined pocket PDB file: {e}")
    
    def _get_hotspot_mask_from_pdb_file(self, original_pdb_file, hotspot_pdb_file):
        """
        Calculate mask from predefined hotspot PDB file
        
        Args:
            original_pdb_file: Original PDB file path
            hotspot_pdb_file: Hotspot PDB file path
            
        Returns:
            hotspot_mask: Hotspot mask array
        """
        try:
            # Parse original PDB file and hotspot PDB file
            parser = PDB.PDBParser(QUIET=True)
            original_structure = parser.get_structure('original', original_pdb_file)[0]
            hotspot_structure = parser.get_structure('hotspot', hotspot_pdb_file)[0]
            
            # Get specified receptor chain IDs
            receptor_chain_ids = getattr(self, '_current_explicit_receptor_chain_ids', [])
            if not receptor_chain_ids:
                if self.verbose:
                    print("Warning: No receptor chain IDs specified, will process all protein chains")
                # If no specified, process all protein chains
                receptor_chain_ids = []
                for chain in sorted(original_structure.get_chains(), key=lambda c: c.id):
                    from data.ligand_cutoff import classify_chain
                    chain_type, _ = classify_chain(chain, verbose=False)
                    if chain_type == 'protein':
                        receptor_chain_ids.append(chain.id)
            
            if self.verbose:
                print(f"Processing receptor chains: {receptor_chain_ids}")
            
            # Get all hotspot residues by chain ID from hotspot PDB
            hotspot_residues_by_chain = {}
            for chain in hotspot_structure.get_chains():
                chain_id = chain.id
                if chain_id not in hotspot_residues_by_chain:
                    hotspot_residues_by_chain[chain_id] = []
                for residue in chain.get_residues():
                    if residue.id[0] == ' ': 
                        hotspot_residues_by_chain[chain_id].append(residue.id[1])
            
            if self.verbose:
                for chain_id, res_ids in hotspot_residues_by_chain.items():
                    print(f"Number of residues in chain {chain_id} in hotspot PDB: {len(res_ids)}")
            
            # Build mask based on chain order in original PDB
            all_receptor_residues = []
            hotspot_mask_values = []
            
            for chain in sorted(original_structure.get_chains(), key=lambda c: c.id):
                chain_id = chain.id
                if chain_id in receptor_chain_ids:
                    chain_residues = [res for res in chain.get_residues() if res.id[0] == ' ']
                    hotspot_res_ids = set(hotspot_residues_by_chain.get(chain_id, []))
                    
                    for residue in chain_residues:
                        res_id = residue.id[1]
                        all_receptor_residues.append((chain_id, res_id))
                        # Check if the residue is a hotspot
                        hotspot_mask_values.append(1 if res_id in hotspot_res_ids else 0)
            
            hotspot_mask = np.array(hotspot_mask_values, dtype=float)
            
            if self.verbose:
                print(f"Total number of receptor residues in original PDB: {len(all_receptor_residues)}")
                hotspot_count = int(np.sum(hotspot_mask))
                print(f"Read {hotspot_count}/{len(all_receptor_residues)} receptor residues from predefined hotspot PDB")
            
            return hotspot_mask
            
        except Exception as e:
            raise ValueError(f"Failed to generate hotspot mask from predefined hotspot PDB file: {e}")
    
    def seq_emb_af3(self, pdb_file, chain_id=None, msa_dir=None, db_dir=None, use_af3_msa=True, use_pocket_msa=False, use_hotspot_msa=False, use_pocket_masked_af3_msa=False, hotspot_cutoff=8.0, pocket_cutoff=10.0, ligand_sequences_override=None, ligand_specs: dict = None, explicit_protein_chain_ids: list = None, explicit_ligand_chain_ids: list = None, rng_seed_override: int = None, predefined_hotspot_pdb: str = None, predefined_pocket_pdb: str = None, receptor_type: str = 'protein', ligand_type: str = 'nucleic'):
        """
        Read PDB file, extract sequence information
        Use pre-trained AlphaFold3 to generate sequence embeddings
        
        Args:
            pdb_file: PDB file path
            chain_id: Select specific chain, if None, process all chains
            msa_dir: MSA file directory, if None, use default directory
            db_dir: MSA database directory, if None, use default directory
            use_af3_msa: Whether to use AlphaFold3's MSA generation process
            use_pocket_msa: Whether to generate unpaired_msa based on pocket for protein chains
            use_hotspot_msa: Whether to generate unpaired_msa based on hotspot for protein chains
            use_pocket_masked_af3_msa: Whether to use AF3 MSA then apply pocket mask to receptor chains only
            hotspot_cutoff: Hotspot residue cutoff distance
            pocket_cutoff: Pocket residue cutoff distance
            ligand_sequences_override: (DEPRECATED, use ligand_specs) Optional dictionary, used to directly provide mapping from ligand chain ID to sequence.
            ligand_specs: Dictionary containing ligand chain ID, sequence, and type ('rna' or 'dna')
            explicit_protein_chain_ids: Optional, explicitly specify protein chain ID list.
            explicit_ligand_chain_ids: Optional, explicitly specify ligand chain ID list.
            rng_seed_override: Optional integer, used to override default random seed in folding_input.
            predefined_hotspot_pdb: Optional predefined hotspot PDB file path
            predefined_pocket_pdb: Optional predefined pocket PDB file path
            receptor_type: Receptor type ('protein' or 'nucleic')
            ligand_type: Ligand type ('protein' or 'nucleic')
        Returns:
            {
                'embeddings': {
                    'single': np.ndarray,  # sequence embedding (token_num, 384)
                    'pair': np.ndarray,    # pair embedding (token_num, token_num, 128)
                    'target_feat': np.ndarray,  # target feature embedding (token_num, 447)
                },
                'seq_info': dict,
                'input_features': dict,
                'pocket_mask': np.ndarray,  # pocket mask (token_num,)
                'center_pos': np.ndarray,  # center position of the complex (3,)
                'protein_mask': np.ndarray,  # protein mask (token_num,)
            }
        """
        if self.verbose:
            print(f"Loading structure from {pdb_file}")
        
        # First ensure model parameters are loaded
        if self.model_params is None:
            raise ValueError("AlphaFold3 model parameters not loaded, cannot generate embeddings")
        
        # Convert PDB file to folding_input (process chains and sequences)
        fold_input = self._pdb_to_fold_input(
            pdb_file, 
            chain_id, 
            use_af3_msa, 
            use_pocket_msa,
            use_hotspot_msa,
            use_pocket_masked_af3_msa,
            hotspot_cutoff,
            pocket_cutoff,
            ligand_sequences_override=ligand_sequences_override,
            ligand_specs=ligand_specs,
            explicit_protein_chain_ids=explicit_protein_chain_ids,
            explicit_ligand_chain_ids=explicit_ligand_chain_ids,
            rng_seed_override=rng_seed_override,
            predefined_hotspot_pdb=predefined_hotspot_pdb,
            predefined_pocket_pdb=predefined_pocket_pdb,
            receptor_type=receptor_type,
            ligand_type=ligand_type
        )
        
        # If using AlphaFold3's MSA process (including pocket masked mode) and no MSA directory provided, need to automatically calculate MSA
        if (use_af3_msa or use_pocket_masked_af3_msa) and msa_dir is None:
            
            # Create data pipeline configuration
            data_pipeline_config = self._create_data_pipeline_config(db_dir=db_dir)
            
            print("Using AlphaFold3 data pipeline to generate MSA and template")
            fold_input = self._run_data_pipeline(fold_input, data_pipeline_config)
            
            # If using pocket masked AF3 MSA mode, apply pocket mask
            if use_pocket_masked_af3_msa:
                if self.verbose:
                    print("Applying pocket mask to AF3 MSA...")
                fold_input = self._apply_pocket_mask_to_af3_msa(
                    fold_input, pdb_file, hotspot_cutoff, pocket_cutoff,
                    explicit_protein_chain_ids, explicit_ligand_chain_ids, 
                    predefined_pocket_pdb, receptor_type, ligand_type
                )
            
            # Check if there are custom MSA data to merge
            if hasattr(self, '_current_custom_msa_data') and self._current_custom_msa_data:
                if self.verbose:
                    print("Detected custom MSA data, merging...")
                fold_input = self._merge_custom_msa_with_af3(fold_input, self._current_custom_msa_data)
                # Clean up temporary data
                self._current_custom_msa_data = None
            
            if self.verbose:
                print("AlphaFold3 MSA and template generation completed")
                self._print_msa_stats(fold_input)
        elif not use_af3_msa and not use_pocket_masked_af3_msa and self.verbose:
            print("Skipping AlphaFold3 MSA generation step")
        
        if self.verbose:
            print("Featurizing input...")
        
        # Featurize input data
        batches = featurisation.featurise_input(
            fold_input=fold_input,
            buckets=None,
            ccd=chemical_components.cached_ccd(),
            verbose=self.verbose,
            conformer_max_iterations=None
        )
        features_batch = batches[0]
        
        print("Running model to generate embeddings...")
        
        key_seed_for_jax = fold_input.rng_seeds[0]
        print(f"Using JAX PRNGKey seed: {key_seed_for_jax} for embedding model call.")
        key = jax.random.PRNGKey(key_seed_for_jax) 
        result = self._forward_fn(key, features_batch)
        
        single_emb = np.array(result['single_embeddings'])
        pair_emb = np.array(result['pair_embeddings'])
        target_feat = np.array(result['target_feat'])
        seq_info = self._extract_sequence_info(fold_input)
        
        # Add receptor chain information
        if explicit_protein_chain_ids:
            for chain_info in seq_info['chains']:
                chain_id = chain_info['chain_id']
                chain_info['is_receptor'] = (chain_id in explicit_protein_chain_ids)
                if self.verbose and chain_info['is_receptor']:
                    print(f"Marking chain {chain_id} as receptor chain")
        else:
            # If no receptor chain IDs are explicitly specified, use heuristic method
            protein_chains = [c for c in seq_info['chains'] if c['chain_type'] == 'protein']
            if len(protein_chains) == 1:
                # Only one protein chain, mark as receptor
                protein_chains[0]['is_receptor'] = True
                print(f"Marking chain {protein_chains[0]['chain_id']} as receptor chain (only one protein chain)")
            elif len(protein_chains) > 1:
                # Multiple protein chains, use the longest as receptor
                longest_chain = max(protein_chains, key=lambda c: len(c['sequence']))
                for chain_info in seq_info['chains']:
                    chain_info['is_receptor'] = (chain_info['chain_id'] == longest_chain['chain_id'])
                print(f"Automatically marking chain {longest_chain['chain_id']} as receptor chain (longest protein chain)")
            else:
                # No protein chains, all chains are not receptors
                for chain_info in seq_info['chains']:
                    chain_info['is_receptor'] = False
        
        # Calculate protein mask
        protein_mask_full, protein_indices, _ = self._get_protein_mask(fold_input)
        
        # Get pocket and hotspot masks from original PDB (these masks are based on original PDB structure atoms)
        # They may not match the final fold_input (if FASTA replacement is used)
        # We will only use the parts corresponding to protein in these masks
        center_pos = np.zeros(3) # Initialize center_pos
        original_pocket_mask_pdb, center_pos_cutoff, _ = self._get_pocket_mask(pdb_file, receptor_type, ligand_type, hotspot_cutoff, pocket_cutoff, explicit_receptor_chain_ids=explicit_protein_chain_ids, explicit_ligand_chain_ids=explicit_ligand_chain_ids, predefined_pocket_pdb=predefined_pocket_pdb)
        if center_pos_cutoff is not None: center_pos = center_pos_cutoff

        original_hotspot_mask_pdb = self._get_hotspot_mask(pdb_file, receptor_type, ligand_type, hotspot_cutoff, explicit_receptor_chain_ids=explicit_protein_chain_ids, explicit_ligand_chain_ids=explicit_ligand_chain_ids, predefined_hotspot_pdb=predefined_hotspot_pdb)

        final_num_tokens = single_emb.shape[0]
        final_pocket_mask = np.zeros(final_num_tokens)
        final_hotspot_mask = np.zeros(final_num_tokens)

        # Check if we need to recalculate masks (when sequence replacement causes length mismatch)
        need_recalculate_masks = False
        if original_pocket_mask_pdb is not None and original_hotspot_mask_pdb is not None:
            # Create a temporary fold_input, based only on PDB, to align original masks and protein indices
            temp_fold_input_for_pdb_masks = self._pdb_to_fold_input(
                pdb_file, None, False, False, False, False, 0, 0, 
                ligand_sequences_override=None, ligand_specs=None, 
                explicit_protein_chain_ids=None, explicit_ligand_chain_ids=None, 
                rng_seed_override=None, predefined_hotspot_pdb=None, predefined_pocket_pdb=None,
                receptor_type=receptor_type, ligand_type=ligand_type
            )
            _, protein_indices_pdb, _ = self._get_protein_mask(temp_fold_input_for_pdb_masks)
            num_res_temp_fold_input = sum(len(c.sequence) for c in temp_fold_input_for_pdb_masks.chains)

            # Check if lengths match
            if (len(original_pocket_mask_pdb) != num_res_temp_fold_input or 
                len(protein_indices_pdb) != len(protein_indices)):
                need_recalculate_masks = True
                if self.verbose:
                    print(f"Detected sequence length change, need to recalculate masks:")
                    print(f"  Original PDB mask length: {len(original_pocket_mask_pdb)}")
                    print(f"  PDB token number: {num_res_temp_fold_input}")
                    print(f"  Original protein index number: {len(protein_indices_pdb)}")
                    print(f"  Final protein index number: {len(protein_indices)}")
            else:
                # Length matches, use original masks
                protein_pocket_mask_values = original_pocket_mask_pdb[protein_indices_pdb]
                protein_hotspot_mask_values = original_hotspot_mask_pdb[protein_indices_pdb]
                
                if len(protein_pocket_mask_values) == len(protein_indices):
                    final_pocket_mask[protein_indices] = protein_pocket_mask_values
                    final_hotspot_mask[protein_indices] = protein_hotspot_mask_values
                else:
                    need_recalculate_masks = True
        else:
            need_recalculate_masks = True
            print("Original mask calculation failed, need to recalculate")

        # If we need to recalculate masks
        if need_recalculate_masks:          
            try:
                # Recalculate pocket mask, using final fold_input information
                recalc_pocket_mask, recalc_center_pos, _ = self._get_pocket_mask(
                    pdb_file, receptor_type, ligand_type, hotspot_cutoff, pocket_cutoff,
                    explicit_receptor_chain_ids=explicit_protein_chain_ids, 
                    explicit_ligand_chain_ids=explicit_ligand_chain_ids, 
                    predefined_pocket_pdb=predefined_pocket_pdb
                )
                
                recalc_hotspot_mask = self._get_hotspot_mask(
                    pdb_file, receptor_type, ligand_type, hotspot_cutoff,
                    explicit_receptor_chain_ids=explicit_protein_chain_ids, 
                    explicit_ligand_chain_ids=explicit_ligand_chain_ids, 
                    predefined_hotspot_pdb=predefined_hotspot_pdb
                )
                
                if recalc_pocket_mask is not None and recalc_hotspot_mask is not None:
                    print(f"Recalculated mask lengths: pocket={len(recalc_pocket_mask)}, hotspot={len(recalc_hotspot_mask)}")
                    print(f"Mapping masks to final token sequence...")
                    # Recalculated masks are based on receptor chain residues, need to correctly map to final token sequence
                    # The order of masks should be consistent with the receptor chain residues in the original PDB
                 
                    # Parse original PDB to get receptor chain residues order
                    from Bio import PDB
                    parser = PDB.PDBParser(QUIET=True)
                    structure = parser.get_structure("temp", pdb_file)[0]
                    
                    # Build receptor chain residues list according to the order of chains in the original PDB
                    pdb_receptor_residues = []
                    for chain in sorted(structure.get_chains(), key=lambda c: c.id):
                        if explicit_protein_chain_ids and chain.id in explicit_protein_chain_ids:
                            chain_residues = [res for res in chain.get_residues() if res.id[0] == ' ']
                            for res in chain_residues:
                                pdb_receptor_residues.append((chain.id, res.id[1]))
                    
                    if self.verbose:
                        print(f"Total receptor residues in PDB: {len(pdb_receptor_residues)}")
                        print(f"Recalculated mask lengths: pocket={len(recalc_pocket_mask)}, hotspot={len(recalc_hotspot_mask)}")
                    
                    # Check if mask lengths match with PDB receptor residues number
                    if (len(recalc_pocket_mask) == len(pdb_receptor_residues) and 
                        len(recalc_hotspot_mask) == len(pdb_receptor_residues)):
                        
                        # Create mapping from PDB residues to mask indices
                        pdb_res_to_mask_idx = {}
                        for i, (chain_id, res_id) in enumerate(pdb_receptor_residues):
                            pdb_res_to_mask_idx[(chain_id, res_id)] = i
                        
                        # Iterate over final sequence, set masks for receptor chain tokens
                        current_token_idx = 0
                        for chain_info in seq_info['chains']:
                            chain_id = chain_info['chain_id']
                            chain_type = chain_info['chain_type']
                            chain_length = len(chain_info['sequence'])
                            
                            if (chain_type == 'protein' and explicit_protein_chain_ids and 
                                chain_id in explicit_protein_chain_ids):
                                # This is receptor chain
                                if self.verbose:
                                    print(f"Processing receptor chain {chain_id}: token range {current_token_idx} - {current_token_idx + chain_length - 1}")
                                
                                if chain_id in structure:
                                    pdb_chain = structure[chain_id]
                                    pdb_residues = [res for res in pdb_chain.get_residues() if res.id[0] == ' ']
                                    
                                    # Set masks for each token of this chain
                                    for i in range(chain_length):
                                        if i < len(pdb_residues):
                                            pdb_res_id = pdb_residues[i].id[1]
                                            mask_key = (chain_id, pdb_res_id)
                                            
                                            if mask_key in pdb_res_to_mask_idx:
                                                mask_idx = pdb_res_to_mask_idx[mask_key]
                                                final_pocket_mask[current_token_idx + i] = recalc_pocket_mask[mask_idx]
                                                final_hotspot_mask[current_token_idx + i] = recalc_hotspot_mask[mask_idx]
                                    
                                    if self.verbose:
                                        chain_pocket_count = np.sum(final_pocket_mask[current_token_idx:current_token_idx + chain_length])
                                        chain_hotspot_count = np.sum(final_hotspot_mask[current_token_idx:current_token_idx + chain_length])
                                        print(f"  chain {chain_id}: pocket={chain_pocket_count}, hotspot={chain_hotspot_count}")
                            
                            current_token_idx += chain_length
                    else:
                        print(f"Warning: mask length does not match with PDB receptor residues number, cannot map masks")
                    
                    pocket_count = np.sum(final_pocket_mask)
                    hotspot_count = np.sum(final_hotspot_mask)
                    print(f"Recalculated masks applied: pocket={pocket_count}, hotspot={hotspot_count}")
                else:
                    print(f"Recalculated masks failed, using zero masks")

                if recalc_center_pos is not None:
                    center_pos = recalc_center_pos
                    
            except Exception as e:
                print(f"Recalculated masks failed: {e}")
                import traceback
                traceback.print_exc()

        if center_pos is None: # Should not happen if _get_pocket_mask succeeded for coords
            center_pos = np.zeros(3)
            print("Warning: center position not calculated, using origin")
        
        # Build results dictionary
        results = {
            'embeddings': {
                'single': single_emb,
                'pair': pair_emb,
                'target_feat': target_feat,
            },
            'seq_info': seq_info,
            'input_features': features_batch,
            'center_pos': center_pos,
            'protein_mask': protein_mask_full, # This is the mask aligned with final fold_input
            'pocket_mask': final_pocket_mask,   # This is the adjusted mask
            'hotspot_mask': final_hotspot_mask, # This is the adjusted mask
        }
        
        return results
    
    def _create_data_pipeline_config(self, db_dir=None):
        """
        Create AlphaFold3 data pipeline configuration
        
        Args:
            db_dir: User-specified database directory path, if None then use default path
            
        Returns:
            pipeline.DataPipelineConfig: 数据管道配置
        """
        from alphafold3.data import pipeline
        import datetime
        import shutil
        
        # Search for executable paths
        jackhmmer_binary_path = shutil.which('jackhmmer')
        nhmmer_binary_path = shutil.which('nhmmer')
        hmmalign_binary_path = shutil.which('hmmalign')
        hmmsearch_binary_path = shutil.which('hmmsearch')
        hmmbuild_binary_path = shutil.which('hmmbuild')
        
        # Check if necessary tools are available
        if not all([jackhmmer_binary_path, nhmmer_binary_path, 
                   hmmalign_binary_path, hmmsearch_binary_path, hmmbuild_binary_path]):
            raise RuntimeError(
                f"Missing necessary MSA tools. Please ensure the following tools are in PATH: "
                f"jackhmmer{'(missing)' if not jackhmmer_binary_path else ''}, "
                f"nhmmer{'(missing)' if not nhmmer_binary_path else ''}, "
                f"hmmalign{'(missing)' if not hmmalign_binary_path else ''}, "
                f"hmmsearch{'(missing)' if not hmmsearch_binary_path else ''}, "
                f"hmmbuild{'(missing)' if not hmmbuild_binary_path else ''}"
            )
        
        # Use user-specified database directory or default directory
        import os
        from pathlib import Path
        
        if db_dir is not None:
            db_dir = Path(db_dir).expanduser()
        else:
            home_dir = Path(os.environ.get('HOME'))
            db_dir = home_dir / 'public_databases'
        
        # Check if database directory exists
        if not db_dir.exists():
            raise RuntimeError(f"Database directory does not exist: {db_dir}")
            
        if self.verbose:
            print(f"Using MSA database directory: {db_dir}")
        
        # Create data pipeline configuration
        data_pipeline_config = pipeline.DataPipelineConfig(
            jackhmmer_binary_path=jackhmmer_binary_path,
            nhmmer_binary_path=nhmmer_binary_path,
            hmmalign_binary_path=hmmalign_binary_path,
            hmmsearch_binary_path=hmmsearch_binary_path,
            hmmbuild_binary_path=hmmbuild_binary_path,
            small_bfd_database_path=str(db_dir / 'bfd-first_non_consensus_sequences.fasta'),
            mgnify_database_path=str(db_dir / 'mgy_clusters_2022_05.fa'),
            uniprot_cluster_annot_database_path=str(db_dir / 'uniprot_all_2021_04.fa'),
            uniref90_database_path=str(db_dir / 'uniref90_2022_05.fa'),
            ntrna_database_path=str(db_dir / 'nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta'),
            rfam_database_path=str(db_dir / 'rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta'),
            rna_central_database_path=str(db_dir / 'rnacentral_active_seq_id_90_cov_80_linclust.fasta'),
            pdb_database_path=str(db_dir / 'mmcif_files'),
            seqres_database_path=str(db_dir / 'pdb_seqres_2022_09_28.fasta'),
            jackhmmer_n_cpu=16,
            nhmmer_n_cpu=16,
            max_template_date=datetime.date.today()
        )
        
        return data_pipeline_config
    
    def _run_data_pipeline(self, fold_input, data_pipeline_config):
        """
        Run AlphaFold3 data pipeline, calculate MSA and templates
        
        Args:
            fold_input: folding_input.Input object
            data_pipeline_config: data pipeline configuration
            
        Returns:
            folding_input.Input: input object with MSA and templates
        """
        from alphafold3.data import pipeline
        
        data_pipeline = pipeline.DataPipeline(data_pipeline_config)
        processed_fold_input = data_pipeline.process(fold_input)
        
        return processed_fold_input
    
    def _print_msa_stats(self, fold_input):
        """
        Print MSA statistics
        
        Args:
            fold_input: folding_input.Input object
        """
        print("MSA statistics:")
        for chain in fold_input.chains:
            if isinstance(chain, folding_input.ProteinChain):
                # unpaired_msa can be str or None
                unpaired_depth = chain.unpaired_msa.count('>') if isinstance(chain.unpaired_msa, str) else 0
                paired_depth = chain.paired_msa.count('>') if isinstance(chain.paired_msa, str) else 0
                print(f"  Protein chain {chain.id}: {len(chain.sequence)} residues")
                print(f"    Unpaired MSA depth: {unpaired_depth}")
                print(f"    Paired MSA depth: {paired_depth}")
                print(f"    Template number: {len(chain.templates) if chain.templates else 0}")
            elif isinstance(chain, folding_input.RnaChain):
                unpaired_depth = chain.unpaired_msa.count('>') if isinstance(chain.unpaired_msa, str) else 0
                print(f"  RNA chain {chain.id}: {len(chain.sequence)} residues")
                print(f"    MSA depth: {unpaired_depth}")
    
    def _pdb_to_fold_input(self, pdb_file, chain_id=None, use_af3_msa=True, use_pocket_msa=False, use_hotspot_msa=False, use_pocket_masked_af3_msa=False, hotspot_cutoff=8.0, pocket_cutoff=10.0, ligand_sequences_override=None, ligand_specs: dict = None, explicit_protein_chain_ids: list = None, explicit_ligand_chain_ids: list = None, rng_seed_override: int = None, predefined_hotspot_pdb: str = None, predefined_pocket_pdb: str = None, receptor_type: str = 'protein', ligand_type: str = 'nucleic'):
        """
        Convert PDB file to AlphaFold3 input format
        
        Args:
            pdb_file: PDB file path
            chain_id: Select specific chain, if None then process all chains
            use_af3_msa: Whether to prepare AF3 MSA for chain (set to None, then pipeline will fill it).
                          If True, and unpaired_msa is None in this function, AF3's pipeline will try to fill it.
            use_pocket_msa: Whether to generate unpaired_msa for protein chain based on pocket.
            use_hotspot_msa: Whether to generate unpaired_msa for protein chain based on hotspot.
            use_pocket_masked_af3_msa: Whether to use AF3 MSA then apply pocket mask to receptor chains only.
            hotspot_cutoff: Hotspot residue cutoff distance
            pocket_cutoff: Pocket residue cutoff distance
            ligand_sequences_override: (DEPRECATED, use ligand_specs) Optional dictionary.
            ligand_specs: Dictionary containing ligand chain ID, sequence and type ('rna' or 'dna').
            explicit_protein_chain_ids: Optional, explicitly specify protein chain ID list.
            explicit_ligand_chain_ids: Optional, explicitly specify ligand chain ID list.
            rng_seed_override: Optional integer, used to override default random seed in folding_input.
            predefined_hotspot_pdb: Optional predefined hotspot PDB file path
            predefined_pocket_pdb: Optional predefined pocket PDB file path
            receptor_type: Receptor type ('protein' or 'nucleic')
            ligand_type: Ligand type ('protein' or 'nucleic')
            
        Returns:
            folding_input.Input object
        """
        # Prioritize ligand_specs if available
        current_ligand_specs = ligand_specs if ligand_specs is not None else {}
        if not current_ligand_specs and ligand_sequences_override:
            if self.verbose: print("Warning: ligand_sequences_override is deprecated. Converting to ligand_specs.")
            for chain_id_override, seq_override in ligand_sequences_override.items():
                current_ligand_specs[chain_id_override] = {"sequence": seq_override, "type": "rna"} # Assume RNA if only old dict given

        # fasta_ligand_seqs = ligand_sequences_override if ligand_sequences_override is not None else {}

        flat_pocket_mask_global = None
        flat_hotspot_mask_global = None
        residue_indices_map_global = None
        custom_msa_data = {}  # Save custom MSA data for later merging
        
        # If any custom MSA features are needed, try to generate global mapping
        attempt_custom_msa_map = use_pocket_msa or use_hotspot_msa
        
        # Dynamically detect molecule type for mask calculation
        auto_receptor_type = 'protein'  # Default
        auto_ligand_type = 'nucleic'    # Default
        
        if attempt_custom_msa_map:
            local_parser = PDB.PDBParser(QUIET=True)
            try:
                local_structure_root = local_parser.get_structure("temp_orig_for_map", pdb_file)
                local_structure = local_structure_root[0]
                local_all_residues = []
                local_sorted_chains = sorted(local_structure.get_chains(), key=lambda c: c.id)
                
                # Analyze chain types to determine receptor and ligand types
                protein_chains = []
                nucleic_chains = []
                small_mol_chains = []
                
                for ch_local in local_sorted_chains:
                    res_list_local = [r for r in ch_local.get_residues() if r.id[0] == ' ']
                    
                    # Classify chain types
                    chain_type, _ = classify_chain(ch_local, verbose=False)
                    if chain_type == 'protein':
                        protein_chains.append(ch_local.id)
                    elif chain_type == 'nucleic':
                        nucleic_chains.append(ch_local.id)
                    elif chain_type == 'small_molecule':
                        small_mol_chains.append(ch_local.id)
                    
                    # Add residues to mapping only for receptor chains (for mask calculation)
                    should_include_in_mapping = False
                    if explicit_protein_chain_ids and ch_local.id in explicit_protein_chain_ids:
                        should_include_in_mapping = True
                    elif not explicit_protein_chain_ids and chain_type == 'protein':
                        # If no explicit receptor chain is specified, include all protein chains
                        should_include_in_mapping = True
                    else:
                        if self.verbose:
                            print(f"Chain {ch_local.id} will not be included in mapping (explicit_protein_chain_ids={explicit_protein_chain_ids}, chain_type={chain_type})")
                    
                    if should_include_in_mapping:
                        local_all_residues.extend([(ch_local.id, r.id[1]) for r in res_list_local])
                        if self.verbose:
                            print(f"Added {len(res_list_local)} residues from receptor chain {ch_local.id} to mapping")
                
                # Determine molecule types based on explicitly specified chains
                if explicit_protein_chain_ids:
                    # If explicitly specified protein chains, they are receptors
                    # Determine ligand type: check other chains not in receptor chains
                    auto_receptor_type = 'protein'
                    remaining_chain_ids = [ch.id for ch in local_sorted_chains if ch.id not in explicit_protein_chain_ids]
                    
                    if explicit_ligand_chain_ids:
                        # If also explicitly specified ligand chains, check their types
                        ligand_chain_types = []
                        for lig_chain_id in explicit_ligand_chain_ids:
                            for ch_local in local_sorted_chains:
                                if ch_local.id == lig_chain_id:
                                    chain_type, _ = classify_chain(ch_local, verbose=False)
                                    ligand_chain_types.append(chain_type)
                                    break
                        
                        # Determine ligand type
                        if 'protein' in ligand_chain_types:
                            auto_ligand_type = 'protein'
                        elif 'nucleic' in ligand_chain_types:
                            auto_ligand_type = 'nucleic'
                        elif 'small_molecule' in ligand_chain_types:
                            auto_ligand_type = 'small_molecule'
                        else:
                            auto_ligand_type = 'nucleic'
                    else:
                        # If no explicitly specified ligand chains, infer based on remaining chain types
                        remaining_chain_types = []
                        for ch_id in remaining_chain_ids:
                            for ch_local in local_sorted_chains:
                                if ch_local.id == ch_id:
                                    chain_type, _ = classify_chain(ch_local, verbose=False)
                                    remaining_chain_types.append(chain_type)
                                    break
                        
                        # Determine ligand type based on remaining chain types
                        if 'nucleic' in remaining_chain_types:
                            auto_ligand_type = 'nucleic'
                        elif 'small_molecule' in remaining_chain_types:
                            auto_ligand_type = 'small_molecule'
                        elif 'protein' in remaining_chain_types:
                            auto_ligand_type = 'protein'
                        else:
                            auto_ligand_type = 'nucleic'
                else:
                    # 自动检测：蛋白质优先作为受体
                    if protein_chains:
                        auto_receptor_type = 'protein'
                        if nucleic_chains:
                            auto_ligand_type = 'nucleic'
                        elif small_mol_chains:
                            auto_ligand_type = 'small_molecule'
                        elif len(protein_chains) > 1:
                            # 多条蛋白质链，可能是蛋白-蛋白复合物
                            auto_ligand_type = 'protein'
                        else:
                            auto_ligand_type = 'nucleic'  # 默认
                    elif nucleic_chains and small_mol_chains:
                        auto_receptor_type = 'nucleic'
                        auto_ligand_type = 'small_molecule'
                    elif nucleic_chains and len(nucleic_chains) > 1:
                        # 多条核酸链，假设一些是受体一些是配体
                        auto_receptor_type = 'nucleic'
                        auto_ligand_type = 'nucleic'
                    else:
                        auto_receptor_type = 'protein'  # 默认
                        auto_ligand_type = 'nucleic'   # 默认
                
                print(f"Automatically detected molecule types for mask calculation: receptor={auto_receptor_type}, ligand={auto_ligand_type}")
                
                # Create mapping assuming mask length matches local_all_residues length
                # Will be validated later based on actual mask length
                if local_all_residues: # Ensure at least one residue
                    residue_indices_map_global = {res_tuple: i for i, res_tuple in enumerate(local_all_residues)}
                    if self.verbose:
                        print(f"Generated preliminary global residue-to-index mapping for custom MSA features (total {len(local_all_residues)} standard residues).")
                else:
                    if self.verbose:
                        print("Warning: No standard residues found in PDB file, cannot generate mapping for custom MSA.")
                    attempt_custom_msa_map = False # Cannot generate mapping, disable custom MSA

            except Exception as e:
                if self.verbose:
                    print(f"Warning: Error occurred while building local mapping for custom MSA: {e}. Custom MSA will be disabled.")
                attempt_custom_msa_map = False # If error, disable custom MSA

        # If mapping succeeded (or attempted) and pocket mask is needed
        if attempt_custom_msa_map and use_pocket_msa:
            # Ensure mask calculation uses the same receptor chain definition as mapping
            receptor_chains_for_mask = explicit_protein_chain_ids if explicit_protein_chain_ids else [ch.id for ch in local_sorted_chains if classify_chain(ch, verbose=False)[0] == 'protein']
            
            temp_pocket_mask_flat, _, _ = self._get_pocket_mask(
                pdb_file, auto_receptor_type, auto_ligand_type, 
                hotspot_cutoff, pocket_cutoff,
                explicit_receptor_chain_ids=receptor_chains_for_mask, 
                explicit_ligand_chain_ids=explicit_ligand_chain_ids, 
                predefined_pocket_pdb=predefined_pocket_pdb
            )
            
            if temp_pocket_mask_flat is not None and residue_indices_map_global and \
               len(temp_pocket_mask_flat) == len(residue_indices_map_global):
                flat_pocket_mask_global = temp_pocket_mask_flat
                if self.verbose:
                    print(f"Generated global pocket mask, length: {len(flat_pocket_mask_global)}")
                    pocket_count = sum(flat_pocket_mask_global)
                    print(f"Number of 1s in global pocket mask: {pocket_count}/{len(flat_pocket_mask_global)}")
                    print(f"First 10 values of global pocket mask: {flat_pocket_mask_global[:10]}")
            else:
                if self.verbose:
                    print(f"Pocket mask does not match residue index mapping length, will be recalculated later")
                    if temp_pocket_mask_flat is not None:
                        print(f"  Pocket mask length: {len(temp_pocket_mask_flat)}")
                    else:
                        print(f"  Pocket mask is None")
                    if residue_indices_map_global:
                        print(f"  Residue index mapping length: {len(residue_indices_map_global)}")
                    else:
                        print(f"  Residue index mapping is None")
                flat_pocket_mask_global = None # Disable, will be recalculated later

        # Hotspot mask
        if attempt_custom_msa_map and use_hotspot_msa:
            # Ensure mask calculation uses the same receptor chain definition as mapping
            receptor_chains_for_mask = explicit_protein_chain_ids if explicit_protein_chain_ids else [ch.id for ch in local_sorted_chains if classify_chain(ch, verbose=False)[0] == 'protein']
            
            temp_hotspot_mask_flat = self._get_hotspot_mask(
                pdb_file, auto_receptor_type, auto_ligand_type, hotspot_cutoff,
                explicit_receptor_chain_ids=receptor_chains_for_mask, 
                explicit_ligand_chain_ids=explicit_ligand_chain_ids, 
                predefined_hotspot_pdb=predefined_hotspot_pdb
            )
            if temp_hotspot_mask_flat is not None and residue_indices_map_global and \
               len(temp_hotspot_mask_flat) == len(residue_indices_map_global):
                flat_hotspot_mask_global = temp_hotspot_mask_flat
                if self.verbose:
                    print(f"Generated global hotspot mask, length: {len(flat_hotspot_mask_global)}")
                    hotspot_count = sum(flat_hotspot_mask_global)
                    print(f"Number of 1s in global hotspot mask: {hotspot_count}/{len(flat_hotspot_mask_global)}")
                    print(f"First 10 values of global hotspot mask: {flat_hotspot_mask_global[:10]}")
            else:
                if self.verbose:
                    print(f"Hotspot mask does not match residue index mapping length, will be recalculated later")
                    if temp_hotspot_mask_flat is not None:
                        print(f"  Hotspot mask length: {len(temp_hotspot_mask_flat)}")
                    else:
                        print(f"  Hotspot mask is None")
                    if residue_indices_map_global:
                        print(f"  Residue index mapping length: {len(residue_indices_map_global)}")
                    else:
                        print(f"  Residue index mapping is None")
                flat_hotspot_mask_global = None # Disable, will be recalculated later

        pdb_parser = PDB.PDBParser(QUIET=True)
        structure = pdb_parser.get_structure("structure", pdb_file)
        chains = []
        processed_pdb_chain_ids = set() # Keep track of PDB chains already processed
        
        for model_pdb in structure:
            sorted_chains_pdb = sorted(model_pdb.get_chains(), key=lambda chain_obj: chain_obj.id) 
            for chain_pdb_obj in sorted_chains_pdb: 
                current_chain_id = chain_pdb_obj.id
                processed_pdb_chain_ids.add(current_chain_id)

                if chain_id is not None and current_chain_id != chain_id:
                    continue
                
                # 1. Determine chain type and sequence (priority: ligand_specs > PDB analysis)
                chain_info = self._determine_chain_type_and_sequence(
                    chain_pdb_obj, current_ligand_specs, self.verbose
                )
                if not chain_info:
                    continue  # Skip invalid chain
                    
                # 2. Process masks (only for protein chains)
                mask_info = self._process_chain_masks(
                    chain_pdb_obj, chain_info, 
                    use_pocket_msa, use_hotspot_msa,
                    residue_indices_map_global, 
                    flat_pocket_mask_global, flat_hotspot_mask_global,
                    attempt_custom_msa_map
                )
                
                # 3. Generate MSA and create chain object
                chain_obj = self._create_chain_with_msa(
                    chain_info, mask_info, 
                    use_af3_msa, use_pocket_msa, use_hotspot_msa,
                    custom_msa_data, current_ligand_specs, use_pocket_masked_af3_msa
                )
                
                if chain_obj:
                    chains.append(chain_obj)
        
        # Add any ligands from ligand_specs that were not in the PDB file
        for spec_chain_id, spec_info in current_ligand_specs.items():
            if spec_chain_id not in processed_pdb_chain_ids:
                if chain_id is not None and spec_chain_id != chain_id:
                    continue

                if spec_info["type"].lower() == "small_molecule":
                    # Special handling for small molecules - don't create chain_info with sequence
                    chain_obj = self._create_small_molecule_chain_from_specs(spec_chain_id, spec_info)
                    
                    # Calculate token count for small molecules based on heavy atoms
                    token_count = "unknown"
                    if spec_info.get("ccd_codes") and hasattr(self, '_current_ccd_heavy_atoms_info'):
                        ccd_heavy_atoms_info = getattr(self, '_current_ccd_heavy_atoms_info', {})
                        total_heavy_atoms = 0
                        for ccd_id in spec_info["ccd_codes"]:
                            if ccd_id in ccd_heavy_atoms_info:
                                total_heavy_atoms += ccd_heavy_atoms_info[ccd_id]
                        token_count = total_heavy_atoms
                    elif spec_info.get("smiles"):
                        try:
                            from rdkit import Chem
                            mol = Chem.MolFromSmiles(spec_info["smiles"])
                            if mol is not None:
                                token_count = mol.GetNumAtoms()
                        except:
                            pass
                    
                    if self.verbose:
                        print(f"Added new chain from ligand_specs: ID={spec_chain_id}, type=SMALL_MOLECULE, tokens={token_count}")
                else:
                    # Create virtual chain info for protein/nucleic chains
                    new_chain_info = {
                        'chain_id': spec_chain_id,
                        'sequence': spec_info.get("sequence", ""),
                        'chain_type': spec_info["type"].lower(),
                        'is_from_ligand_specs': True
                    }
                    
                    if not new_chain_info['sequence']:
                        if self.verbose:
                            print(f"Warning: Chain {spec_chain_id} in ligand_specs has empty sequence, skipping.")
                        continue
                    
                    if self.verbose:
                        print(f"Added new chain from ligand_specs: ID={spec_chain_id}, type={new_chain_info['chain_type'].upper()}, length={len(new_chain_info['sequence'])}")

                    # Empty mask info (new chains don't need masks)
                    empty_mask_info = {
                        'pocket_mask_values': [],
                        'hotspot_mask_values': []
                    }
                    
                    chain_obj = self._create_chain_with_msa(
                        new_chain_info, empty_mask_info,
                        use_af3_msa, use_pocket_msa, use_hotspot_msa,
                        custom_msa_data, current_ligand_specs, use_pocket_masked_af3_msa
                    )
                
                if chain_obj:
                    chains.append(chain_obj)

        if not chains:
            raise ValueError(f"No valid chains found in {pdb_file} or ligand_specs after type classification, or all chain sequences are empty.")
        
        pdb_name = os.path.basename(pdb_file).replace(".pdb", "")
        
        # Determine the RNG seeds for folding_input
        final_rng_seeds = [42] # Default if nothing else is provided
        if rng_seed_override is not None:
            final_rng_seeds = [rng_seed_override]
            if self.verbose:
                print(f"folding_input will use rng_seed_override: {final_rng_seeds}")
        
        fold_input = folding_input.Input(
            name=pdb_name, chains=chains, rng_seeds=final_rng_seeds)
        self._current_custom_msa_data = custom_msa_data
        
        return fold_input
    
    def _extract_sequence_info(self, fold_input):
        """
        Extract sequence information from fold_input
        
        Args:
            fold_input: folding_input.Input object
            
        Returns:
            Dictionary containing sequence information
        """
        seq_info = {
            'name': fold_input.name,
            'chains': []
        }
        
        current_total_residues = 0
        for chain in fold_input.chains:
            chain_type = "unknown"
            sequence = ""
            chain_length = 0
            
            if isinstance(chain, folding_input.ProteinChain):
                chain_type = "protein"
                sequence = chain.sequence
                chain_length = len(chain.sequence)
            elif isinstance(chain, folding_input.RnaChain):
                chain_type = "rna"
                sequence = chain.sequence
                chain_length = len(chain.sequence)
            elif isinstance(chain, folding_input.DnaChain):
                chain_type = "dna"
                sequence = chain.sequence
                chain_length = len(chain.sequence)
            elif isinstance(chain, folding_input.Ligand):
                chain_type = "ligand"
                # For small molecules, determine sequence representation and calculate heavy atoms
                if chain.smiles is not None:
                    # Use SMILES string as "sequence", The number of tokens for small 
                    # molecules is the number of heavy atoms, need to calculate from SMILES
                    sequence = chain.smiles
                    try:
                        from rdkit import Chem
                        mol = Chem.MolFromSmiles(chain.smiles)
                        if mol is not None:
                            chain_length = mol.GetNumHeavyAtoms()
                        else:
                            if self.verbose:
                                print(f"Warning: Failed to parse SMILES for small molecule chain {chain.id}: {chain.smiles}")
                            chain_length = 1  # Default at least 1 token
                    except ImportError:
                        if self.verbose:
                            print(f"Warning: rdkit is not installed, cannot calculate number of heavy atoms for small molecule. Using SMILES length as approximation.")
                        chain_length = len(chain.smiles) if chain.smiles else 1
                    except Exception as e:
                        if self.verbose:
                            print(f"Warning: Error calculating number of heavy atoms for small molecule {chain.id}: {e}")
                        chain_length = 1  # Default at least 1 token
                elif chain.ccd_ids is not None:
                    # For CCD codes, calculate heavy atoms first
                    chain_length = 0
                    total_heavy_atoms = 0
                    
                    # Try to get heavy atoms info from current CCD extraction process
                    # Check if heavy atoms info has been calculated in run_SiteAF3.py
                    if hasattr(self, '_current_ccd_heavy_atoms_info'):
                        ccd_heavy_atoms_info = getattr(self, '_current_ccd_heavy_atoms_info', {})
                        for ccd_id in chain.ccd_ids:
                            if ccd_id in ccd_heavy_atoms_info:
                                total_heavy_atoms += ccd_heavy_atoms_info[ccd_id]
                        chain_length = total_heavy_atoms
                        # Create a placeholder sequence string that matches the token count
                        # Use 'X' for each heavy atom to maintain correct length
                        sequence = 'X' * chain_length
                        if self.verbose:
                            print(f"Small molecule chain {chain.id}: CCD codes {chain.ccd_ids}, heavy atoms: {total_heavy_atoms}, sequence placeholder: {sequence[:10]}...")
                    else:
                        raise ValueError(f"No cached CCD heavy atoms info, cannot calculate number of heavy atoms for small molecule {chain.id}")
                else:
                    # Neither SMILES nor CCD codes, this should not happen
                    sequence = ""
                    chain_length = 1
                    if self.verbose:
                        print(f"Warning: Small molecule chain {chain.id} neither has SMILES nor CCD codes.")
            
            seq_info['chains'].append({
                'chain_id': chain.id,
                'sequence': sequence,
                'chain_type': chain_type,
                'start_residue_index': current_total_residues, # inclusive
                'end_residue_index': current_total_residues + chain_length - 1 # inclusive
            })
            current_total_residues += chain_length
        seq_info['total_residues'] = current_total_residues
            
        return seq_info

    def struct_emb_af3(self, pdb_file, db_dir=None, msa_dir=None, use_af3_msa=True, use_pocket_msa=False, use_hotspot_msa=False, use_pocket_masked_af3_msa=False, hotspot_cutoff=8.0, pocket_cutoff=10.0, ligand_sequences_override=None, ligand_specs: dict = None, explicit_protein_chain_ids: list = None, explicit_ligand_chain_ids: list = None, rng_seed_override: int = None, predefined_hotspot_pdb: str = None, predefined_pocket_pdb: str = None, receptor_type: str = 'protein', ligand_type: str = 'nucleic'):
        """
        Generate AlphaFold 3 compatible embeddings, including structure information and atom masks.

        Args:
            pdb_file: Input PDB file path
            db_dir: Database directory (for MSA generation)
            msa_dir: Precomputed MSA file directory
            use_af3_msa: Whether to use AlphaFold3 MSA features
            use_pocket_msa: Whether to generate unpaired_msa based on pockets for protein chains
            use_hotspot_msa: Whether to generate unpaired_msa based on hotspots for protein chains
            use_pocket_masked_af3_msa: Whether to use AF3 MSA then apply pocket mask to receptor chains only
            hotspot_cutoff: Hotspot residue cutoff distance (Å)
            pocket_cutoff: Pocket residue cutoff distance (Å)
            ligand_sequences_override: (DEPRECATED) Optional dictionary.
            ligand_specs: Dictionary containing ligand chain IDs, sequences, and types ('rna' or 'dna')
            explicit_protein_chain_ids: Optional, explicitly specify protein chain IDs list.
            explicit_ligand_chain_ids: Optional, explicitly specify ligand chain IDs list.
            rng_seed_override: Optional integer, to override default random seeds in folding_input, and used for JAX PRNGKey.
            predefined_hotspot_pdb: Optional predefined hotspot PDB file path
            predefined_pocket_pdb: Optional predefined pocket PDB file path

        Returns:
            A dictionary containing embeddings and related information
            'embeddings': {'single': ..., 'pair': ..., 'structure_atom_coords': ..., 'structure_atom_mask': ...}
            'seq_info': ...
        """
        print(f"Start generating structure embeddings: {pdb_file}")

        # Parse PDB to get structure
        parser = PDB.PDBParser(QUIET=True)
        try:
            structure = parser.get_structure("molecule", pdb_file)
        except Exception as e:
            raise ValueError(f"Failed to parse PDB file {pdb_file}: {e}")

        # Parse all atom coordinates
        try:
            all_atom_coords = parse_pdb_atoms(pdb_file)
        except Exception as e:
            raise ValueError(f"Error parsing PDB atom coordinates: {e}")

        # Use seq_emb_af3 to generate basic sequence embeddings
        # Note: Here we pass chain_id=None to ensure all chains are included in the base embeddings
        base_embeddings = self.seq_emb_af3(
            pdb_file, 
            chain_id=None, 
            msa_dir=msa_dir, 
            db_dir=db_dir, 
            use_af3_msa=use_af3_msa,
            use_pocket_msa=use_pocket_msa,
            use_hotspot_msa=use_hotspot_msa,
            use_pocket_masked_af3_msa=use_pocket_masked_af3_msa,
            hotspot_cutoff=hotspot_cutoff,
            pocket_cutoff=pocket_cutoff,
            ligand_sequences_override=ligand_sequences_override,
            ligand_specs=ligand_specs,
            explicit_protein_chain_ids=explicit_protein_chain_ids,
            explicit_ligand_chain_ids=explicit_ligand_chain_ids,
            rng_seed_override=rng_seed_override,
            predefined_hotspot_pdb=predefined_hotspot_pdb,
            predefined_pocket_pdb=predefined_pocket_pdb,
            receptor_type=receptor_type,
            ligand_type=ligand_type
        )
        num_tokens = base_embeddings['embeddings']['single'].shape[0]
        
        # --- Create mapping from PDB residue ID to internal index ---
        pdb_residue_map = {} # (pdb_chain, pdb_res_id) -> (internal_chain_id, chain_idx)
        pdb_to_embedding = {} # (internal_chain_id, chain_idx) -> embed_idx
        embedding_to_pdb = {} # embed_idx -> (internal_chain_id, chain_idx)

        current_idx = 0
        for chain_info in base_embeddings['seq_info']['chains']:
            internal_chain_id = chain_info['chain_id']
            chain_type = chain_info['chain_type']
            chain_length = chain_info['end_residue_index'] - chain_info['start_residue_index'] + 1
            
            # Find the corresponding chain in PDB
            pdb_chain = None
            for model in structure:
                if internal_chain_id in model:
                    pdb_chain = model[internal_chain_id]
                    break
            
            if pdb_chain:
                if chain_type == "ligand":
                    # Get all atoms (including HETATM)
                    atom_list = []
                    for residue in pdb_chain.get_residues():
                        for atom in residue.get_atoms():
                            atom_list.append(atom)
                    
                    if len(atom_list) != chain_length and self.verbose:
                        print(f"Warning: Calculated number of atoms ({chain_length}) for small molecule chain {internal_chain_id} does not match the number of atoms in PDB ({len(atom_list)}).")
                    
                    # Create mapping for each token (use calculated chain_length to match embeddings)
                    for i in range(chain_length):
                        if i < len(atom_list):
                            atom = atom_list[i]
                            # For small molecules, use the unique identifier of the atom
                            atom_key = (internal_chain_id, f"atom_{i}")
                            pdb_residue_map[atom_key] = (internal_chain_id, i)
                        
                        pdb_to_embedding[(internal_chain_id, i)] = current_idx
                        embedding_to_pdb[current_idx] = (internal_chain_id, i)
                        current_idx += 1
                else:
                    # Protein/nucleic acid processing: each residue corresponds to a token
                    res_list_pdb = [res for res in pdb_chain if res.id[0] == ' '] # Get standard residue list in PDB
                    
                    # Ensure the number of standard residues in PDB matches the sequence length
                    if len(res_list_pdb) != chain_length:
                         print(f"Warning: Sequence length ({chain_length}) for chain {internal_chain_id} does not match the number of standard residues ({len(res_list_pdb)}) in PDB. Mapping may be inaccurate.")
                    
                    for i in range(chain_length): # i is the index in the sequence/embedding (0-based)
                        if i < len(res_list_pdb):
                             pdb_res = res_list_pdb[i]
                             pdb_res_id = pdb_res.id[1]
                             # Store mapping: (PDB chain ID, PDB residue ID) -> (internal chain ID, chain index)
                             pdb_residue_map[(internal_chain_id, pdb_res_id)] = (internal_chain_id, i)

                        # Store mapping: (internal chain ID, chain index) -> embedding index
                        pdb_to_embedding[(internal_chain_id, i)] = current_idx
                        # Store reverse mapping: embedding index -> (internal chain ID, chain index)
                        embedding_to_pdb[current_idx] = (internal_chain_id, i)
                        current_idx += 1
            else:
                 print(f"Warning: Chain {internal_chain_id} not found in PDB structure (exists in sequence info). Skipping mapping for this chain.")
                 # Still need to assign embedding index for this chain's token
                 for i in range(chain_length):
                     embedding_to_pdb[current_idx] = (internal_chain_id, i)
                     current_idx += 1

        # --- Mapping creation ends ---

        # Create structure embedding and atom mask
        structure_embedding = np.zeros((num_tokens, MAX_ATOMS_PER_TOKEN, 3))
        atom_mask = np.zeros((num_tokens, MAX_ATOMS_PER_TOKEN), dtype=np.float32) # Atom mask

        # Single loop to process atom_mask and structure_embedding
        for embed_idx in range(num_tokens):
            # 1. Get residue information
            if embed_idx not in embedding_to_pdb:
                if self.verbose:
                    print(f"Warning: embed_idx {embed_idx} not found in embedding_to_pdb mapping. Skipping this token.")
                continue

            internal_chain_id, residue_idx_in_chain = embedding_to_pdb[embed_idx]
            chain_info_found = None
            for ch_info in base_embeddings['seq_info']['chains']:
                if ch_info['chain_id'] == internal_chain_id:
                    chain_info_found = ch_info
                    break
            
            if not chain_info_found:
                if self.verbose:
                    print(f"Warning: Chain {internal_chain_id} not found in seq_info (embed_idx {embed_idx}). Skipping.")
                continue

            if residue_idx_in_chain >= len(chain_info_found['sequence']):
                if self.verbose:
                    print(f"Warning: residue_idx_in_chain {residue_idx_in_chain} exceeds the sequence length of chain {internal_chain_id} ({len(chain_info_found['sequence'])}). Skipping.")
                continue
                
            residue_char = chain_info_found['sequence'][residue_idx_in_chain]
            chain_type = chain_info_found['chain_type']

            # 2. Determine canonical atom list and current token's atom order dictionary
            canonical_atom_names_for_residue = []
            current_atom_order_for_token = {} 

            try:
                if chain_type == "protein":
                    res_name_3letter = residue_constants.restype_1to3.get(residue_char)
                    if res_name_3letter:
                        if res_name_3letter == "UNK":
                            pass # No standard atoms for UNK protein
                        elif res_name_3letter in atom_types.RESIDUE_ATOMS:
                            canonical_atom_names_for_residue = atom_types.RESIDUE_ATOMS[res_name_3letter]
                    current_atom_order_for_token = atom_types.ATOM37_ORDER
                
                elif chain_type == "rna":
                    if residue_char in atom_types.DENSE_ATOM: # Keys like 'A', 'U', 'C', 'G'
                        canonical_atom_names_for_residue = atom_types.DENSE_ATOM[residue_char]
                    elif residue_char == 'N':
                        pass # No standard atoms for UNK RNA (N)
                    current_atom_order_for_token = atom_types.ATOM29_ORDER

                elif chain_type == "dna":
                    dna_key = f"D{residue_char}" # e.g. 'DA', 'DT'
                    if residue_char == 'N':
                        pass # No standard atoms for UNK DNA (N)
                    elif dna_key in atom_types.DENSE_ATOM:
                        canonical_atom_names_for_residue = atom_types.DENSE_ATOM[dna_key]
                    current_atom_order_for_token = atom_types.ATOM29_ORDER
                
                elif chain_type == "ligand":
                    # For heavy atoms calculated from CCD/SMILES, all tokens should have mask, regardless of whether they exist in PDB
                    found_atom = None
                    for model in structure:
                        if internal_chain_id in model:
                            chain_pdb = model[internal_chain_id]
                            atom_idx = 0
                            for residue in chain_pdb.get_residues():
                                for atom in residue.get_atoms():
                                    if atom_idx == residue_idx_in_chain:
                                        found_atom = atom
                                        break
                                    atom_idx += 1
                                if found_atom:
                                    break
                            break
                    
                    if found_atom:
                        # For atoms that exist in PDB, use actual atom name
                        atom_name = found_atom.get_name().strip()
                        canonical_atom_names_for_residue = [atom_name]
                        current_atom_order_for_token = {atom_name: 0}
                        if self.verbose:
                            print(f"Small molecule atom {internal_chain_id}:{atom_name} (element: {found_atom.element.strip() if hasattr(found_atom, 'element') and found_atom.element else atom_name[0]})")
                    else:
                        # For atoms that do not exist in PDB, also need to set mask so that diffusion process can predict
                        generic_atom_name = f"ATOM_{residue_idx_in_chain}"
                        canonical_atom_names_for_residue = [generic_atom_name]
                        current_atom_order_for_token = {generic_atom_name: 0}
                        if self.verbose:
                            print(f"Small molecule chain {internal_chain_id} has no atom at residue {residue_idx_in_chain} in PDB, but still set mask for diffusion prediction")
                
                elif self.verbose:
                     print(f"Info: Chain type {chain_type} (residue {residue_char}, embed_idx {embed_idx}) has no predefined atom order. Mask and coordinates will not be filled.")
            
            except KeyError as e:
                if self.verbose:
                    print(f"Warning: KeyError occurred while finding canonical atom name: {e} (residue: {residue_char}, type: {chain_type}, embed_idx: {embed_idx})")
            except Exception as e:
                if self.verbose:
                    print(f"Warning: Unknown error occurred while determining canonical atom: {e} (residue: {residue_char}, type: {chain_type}, embed_idx: {embed_idx})")

            # 3. Fill atom_mask
            if current_atom_order_for_token: # Proceed if an atom order is defined for this chain type
                if not canonical_atom_names_for_residue and self.verbose and residue_char not in ['X', 'N']:
                     # This case might happen if residue_char is valid but not in RESIDUE_ATOMS/DENSE_ATOM for some reason
                     print(f"Warning: Residue {residue_char} (type: {chain_type}, embed_idx {embed_idx}) has no canonical atom list, but atom order is defined. Mask may be incomplete.")               
                
                # For protein chain's C-terminal residue, add OXT atom
                final_canonical_atoms = list(canonical_atom_names_for_residue)
                if (chain_type == "protein" and 
                    residue_idx_in_chain == len(chain_info_found['sequence']) - 1 and  # Last residue of the chain
                    'OXT' in current_atom_order_for_token and  # OXT in atom order
                    'OXT' not in final_canonical_atoms):  # Avoid duplicate addition
                    final_canonical_atoms.append('OXT')
                    if self.verbose:
                        print(f"For protein chain {internal_chain_id} C-terminal residue {residue_char} (embed_idx {embed_idx}), add OXT atom")
                
                # Fill atom_mask for this token
                pos = 0
                for atom_name_canonical in final_canonical_atoms:
                    atom_idx_canonical = current_atom_order_for_token.get(atom_name_canonical)
                    if atom_idx_canonical is not None and pos < MAX_ATOMS_PER_TOKEN: 
                        atom_mask[embed_idx, pos] = 1.0
                        pos += 1
                        
            # 4. Fill structure_embedding (using the same current_atom_order_for_token)
            if current_atom_order_for_token:
                if chain_type == "ligand":
                    # For small molecules, directly get atom coordinates from PDB
                    found_atom_coords = None
                    found_atom = None
                    for model in structure:
                        if internal_chain_id in model:
                            chain_pdb = model[internal_chain_id]
                            atom_idx = 0
                            for residue in chain_pdb.get_residues():
                                for atom in residue.get_atoms():
                                    if atom_idx == residue_idx_in_chain:
                                        found_atom_coords = atom.get_coord()
                                        found_atom = atom
                                        break
                                    atom_idx += 1
                                if found_atom_coords is not None:
                                    break
                            break
                    
                    if found_atom_coords is not None:
                        # Fill small molecule atom coordinates to the first position
                        structure_embedding[embed_idx, 0, :] = found_atom_coords
                        if self.verbose:
                            print(f"Small molecule atom coordinates {internal_chain_id} token {embed_idx}: {found_atom_coords}")
                    else:
                        # For CCD-defined atoms that don't exist in PDB, set default coordinates
                        # This allows the diffusion process to predict their positions
                        default_coords = np.array([0.0, 0.0, 0.0])
                        structure_embedding[embed_idx, 0, :] = default_coords
                        if self.verbose:
                            print(f"Small molecule chain {internal_chain_id} token {embed_idx}: using default coordinates (will be predicted by diffusion)")
                else:
                    # Protein/nucleic acid: process as usual
                    found_pdb_key = None
                    for pdb_key_iter, af3_val_iter in pdb_residue_map.items():
                        if af3_val_iter == (internal_chain_id, residue_idx_in_chain):
                            found_pdb_key = pdb_key_iter
                            break
                    
                    if found_pdb_key:
                        residue_pdb_atoms = all_atom_coords.get(found_pdb_key)
                        if residue_pdb_atoms:
                            pos = 0
                            for atom_name_pdb, coords in residue_pdb_atoms.items():
                                atom_idx_coords = current_atom_order_for_token.get(atom_name_pdb)
                                if atom_idx_coords is not None and pos < MAX_ATOMS_PER_TOKEN: 
                                    structure_embedding[embed_idx, pos, :] = coords
                                    pos += 1

        # Build result dictionary
        result = copy.deepcopy(base_embeddings)
        result['embeddings']['structure_atom_coords'] = structure_embedding
        result['embeddings']['structure_atom_mask'] = atom_mask

        if self.verbose:
            print(f"Structure embedding generation completed.")
            print(f"Structure coordinate embedding shape: {structure_embedding.shape}")
            print(f"Atom mask shape: {atom_mask.shape}")
            print(f"Number of 1s in atom mask: {int(np.sum(atom_mask))}")

        return result

    def _merge_custom_msa_with_af3(self, fold_input, custom_msa_data):
        """
        After AF3 pipeline processing, merge custom MSA into the result
        
        Args:
            fold_input: folding_input.Input object after AF3 pipeline processing
            custom_msa_data: Custom MSA data dictionary, format: {chain_id: custom_msa_string}
            
        Returns:
            folding_input.Input: Input object with custom MSA merged
        """
        if not custom_msa_data:
            return fold_input
            
        print("Merge custom MSA with AF3-generated MSA...")
            
        # Create new chain list
        new_chains = []
        
        for chain in fold_input.chains:
            if isinstance(chain, folding_input.ProteinChain) and chain.id in custom_msa_data:
                custom_msa = custom_msa_data[chain.id]
                
                # Merge unpaired_msa
                if chain.unpaired_msa and custom_msa:
                    merged_unpaired_msa = custom_msa + "\n" + chain.unpaired_msa
                elif custom_msa:
                    merged_unpaired_msa = custom_msa
                else: # Only AF3 MSA
                    merged_unpaired_msa = chain.unpaired_msa
                
                if self.verbose:
                    custom_count = custom_msa.count('>') if custom_msa else 0
                    af3_count = chain.unpaired_msa.count('>') if chain.unpaired_msa else 0
                    total_count = merged_unpaired_msa.count('>') if merged_unpaired_msa else 0
                    print(f"Chain {chain.id}: custom MSA={custom_count} lines, AF3 MSA={af3_count} lines, merged={total_count} lines")
                
                # Create new protein chain object
                new_chain = folding_input.ProteinChain(
                    id=chain.id,
                    sequence=chain.sequence,
                    ptms=chain.ptms,
                    unpaired_msa=merged_unpaired_msa,
                    paired_msa=chain.paired_msa,
                    templates=chain.templates
                )
                new_chains.append(new_chain)
            else:
                # Non-protein chain or chain without custom MSA, keep as is
                new_chains.append(chain)
        
        # Create new fold_input object
        new_fold_input = folding_input.Input(
            name=fold_input.name,
            chains=new_chains,
            rng_seeds=fold_input.rng_seeds
        )
        
        return new_fold_input

    def _apply_pocket_mask_to_af3_msa(self, fold_input, pdb_file, hotspot_cutoff, pocket_cutoff,
                                    explicit_protein_chain_ids, explicit_ligand_chain_ids,
                                    predefined_pocket_pdb, receptor_type, ligand_type):
        """
        Apply pocket mask to AF3 generated MSA for receptor protein chains only.
        The first sequence (query sequence) in each MSA is kept unchanged.
        Only subsequent sequences are masked with pocket residues.
        Non-receptor protein chains (ligands) will keep their original AF3 MSA unchanged.
        
        Both unpaired and paired MSAs are processed. Handles sequence alignment mismatches
        by truncating or padding sequences as needed.
        
        Args:
            fold_input: folding_input.Input object with AF3 generated MSA
            pdb_file: PDB file path
            hotspot_cutoff: Hotspot residue cutoff distance
            pocket_cutoff: Pocket residue cutoff distance
            explicit_protein_chain_ids: Receptor chain IDs (only these chains will have MSA masked)
            explicit_ligand_chain_ids: Ligand chain IDs
            predefined_pocket_pdb: Predefined pocket PDB path
            receptor_type: Receptor type
            ligand_type: Ligand type
            
        Returns:
            folding_input.Input: Input object with pocket masked MSA for receptor chains only
        """
        if self.verbose:
            print("Generating pocket masks for AF3 MSA masking...")
        
        # Check if receptor chain IDs are provided
        if not explicit_protein_chain_ids:
            if self.verbose:
                print("Warning: No explicit receptor chain IDs provided for pocket masked AF3 MSA. Skipping pocket masking.")
            return fold_input
        
        try:
            # Get pocket masks for protein chains
            pocket_masks_by_chain = {}
            
            # Generate pocket mask using the same logic as other MSA modes
            original_pocket_mask_pdb, _, _ = self._get_pocket_mask(
                pdb_file, receptor_type, ligand_type, hotspot_cutoff, pocket_cutoff,
                explicit_receptor_chain_ids=explicit_protein_chain_ids,
                explicit_ligand_chain_ids=explicit_ligand_chain_ids,
                predefined_pocket_pdb=predefined_pocket_pdb
            )
            
            if original_pocket_mask_pdb is None:
                if self.verbose:
                    print("Warning: Could not generate pocket mask, returning original fold_input")
                return fold_input
            
            # Parse PDB to get chain residue mapping
            from Bio import PDB
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("temp", pdb_file)[0]
            
            # Build mapping from chain residues to pocket mask indices
            pdb_receptor_residues = []
            for chain in sorted(structure.get_chains(), key=lambda c: c.id):
                if explicit_protein_chain_ids and chain.id in explicit_protein_chain_ids:
                    chain_residues = [res for res in chain.get_residues() if res.id[0] == ' ']
                    for res in chain_residues:
                        pdb_receptor_residues.append((chain.id, res.id[1]))
            
            if len(original_pocket_mask_pdb) != len(pdb_receptor_residues):
                if self.verbose:
                    print(f"Warning: Pocket mask length ({len(original_pocket_mask_pdb)}) doesn't match receptor residues ({len(pdb_receptor_residues)})")
                return fold_input
            
            # Create chain-specific pocket masks
            current_mask_idx = 0
            for chain in sorted(structure.get_chains(), key=lambda c: c.id):
                if explicit_protein_chain_ids and chain.id in explicit_protein_chain_ids:
                    chain_residues = [res for res in chain.get_residues() if res.id[0] == ' ']
                    chain_mask = []
                    for _ in chain_residues:
                        if current_mask_idx < len(original_pocket_mask_pdb):
                            chain_mask.append(original_pocket_mask_pdb[current_mask_idx])
                            current_mask_idx += 1
                    pocket_masks_by_chain[chain.id] = chain_mask
                    
                    if self.verbose:
                        mask_count = sum(chain_mask)
                        print(f"Chain {chain.id}: pocket mask with {mask_count}/{len(chain_mask)} residues marked")
            
            # Apply pocket masks to receptor protein chains only in fold_input
            new_chains = []
            for chain in fold_input.chains:
                if (isinstance(chain, folding_input.ProteinChain) and 
                    chain.id in pocket_masks_by_chain and
                    explicit_protein_chain_ids and 
                    chain.id in explicit_protein_chain_ids):
                    
                    # This is a receptor chain, apply pocket mask
                    pocket_mask = pocket_masks_by_chain[chain.id]
                    
                    # Apply mask to unpaired MSA using AF3's native method
                    masked_unpaired_msa = chain.unpaired_msa
                    if chain.unpaired_msa:
                        if AF3_MSA_MASKING_AVAILABLE:
                            masked_unpaired_msa = apply_pocket_mask_to_msa_af3(
                                chain.unpaired_msa, pocket_mask, 'protein', verbose=self.verbose
                            )
                        else:
                            masked_unpaired_msa = apply_pocket_mask_to_msa(
                                chain.unpaired_msa, pocket_mask, verbose=self.verbose
                            )
                    
                    # Apply mask to paired MSA using AF3's native method
                    masked_paired_msa = chain.paired_msa
                    if chain.paired_msa:
                        if AF3_MSA_MASKING_AVAILABLE:
                            masked_paired_msa = apply_pocket_mask_to_msa_af3(
                                chain.paired_msa, pocket_mask, 'protein', verbose=self.verbose
                            )
                        else:
                            masked_paired_msa = apply_pocket_mask_to_msa(
                                chain.paired_msa, pocket_mask, verbose=self.verbose
                            )
                    
                    # Create new protein chain with masked MSAs
                    new_chain = folding_input.ProteinChain(
                        id=chain.id,
                        sequence=chain.sequence,
                        ptms=chain.ptms,
                        unpaired_msa=masked_unpaired_msa,
                        paired_msa=masked_paired_msa,
                        templates=chain.templates
                    )
                    new_chains.append(new_chain)
                    
                    if self.verbose:
                        unpaired_original_count = chain.unpaired_msa.count('>') if chain.unpaired_msa else 0
                        unpaired_masked_count = masked_unpaired_msa.count('>') if masked_unpaired_msa else 0
                        paired_original_count = chain.paired_msa.count('>') if chain.paired_msa else 0
                        paired_masked_count = masked_paired_msa.count('>') if masked_paired_msa else 0
                        print(f"Receptor chain {chain.id}: Unpaired MSA masked from {unpaired_original_count} to {unpaired_masked_count} sequences")
                        print(f"Receptor chain {chain.id}: Paired MSA masked from {paired_original_count} to {paired_masked_count} sequences")
                else:
                    # Non-receptor chain, ligand chain, or chain without pocket mask, keep original MSA
                    new_chains.append(chain)
                    if (isinstance(chain, folding_input.ProteinChain) and 
                        explicit_protein_chain_ids and 
                        chain.id not in explicit_protein_chain_ids and 
                        self.verbose):
                        print(f"Non-receptor protein chain {chain.id}: keeping original AF3 MSA")
            
            # Create new fold_input object
            new_fold_input = folding_input.Input(
                name=fold_input.name,
                chains=new_chains,
                rng_seeds=fold_input.rng_seeds
            )
            
            return new_fold_input
            
        except Exception as e:
            if self.verbose:
                print(f"Error applying pocket mask to AF3 MSA: {e}")
                import traceback
                traceback.print_exc()
            return fold_input

    def _determine_chain_type_and_sequence(self, chain_pdb_obj, current_ligand_specs, verbose=False):
        """
        Determine chain type and sequence, priority: ligand_specs > PDB analysis
        
        Args:
            chain_pdb_obj: PDB chain object
            current_ligand_specs: ligand_specs configuration dictionary
            verbose: Whether to output detailed information
            
        Returns:
            dict: Dictionary containing chain information, if invalid, return None
                - chain_id: Chain ID
                - sequence: Sequence string
                - chain_type: Chain type ('protein', 'rna', 'dna', 'small_molecule')
                - is_from_ligand_specs: Whether from ligand_specs configuration
        """
        current_chain_id = chain_pdb_obj.id
        
        # Check ligand_specs configuration first
        if current_chain_id in current_ligand_specs:
            spec = current_ligand_specs[current_chain_id]
            lig_type_from_spec = spec["type"].lower()
            
            if lig_type_from_spec == "small_molecule":
                sequence_to_use = spec.get("smiles", "")
                chain_type_to_use = "small_molecule"
                if verbose:
                    print(f"Chain {current_chain_id} is a small molecule ligand (SMILES: {sequence_to_use})")
            else:
                sequence_to_use = spec["sequence"]
                if lig_type_from_spec == "rna":
                    chain_type_to_use = "rna"
                elif lig_type_from_spec == "dna":
                    chain_type_to_use = "dna"
                elif lig_type_from_spec == "protein":
                    chain_type_to_use = "protein"
                else:
                    if verbose:
                        print(f"Warning: Unknown ligand type '{spec['type']}' for chain {current_chain_id}.")
                    chain_type_to_use = "rna"
                
                if verbose:
                    print(f"Chain {current_chain_id} sequence and type are overridden by ligand_specs "
                          f"(type: {chain_type_to_use.upper()}, length: {len(sequence_to_use)}): "
                          f"{sequence_to_use[:20]}...")
            
            return {
                'chain_id': current_chain_id,
                'sequence': sequence_to_use,
                'chain_type': chain_type_to_use,
                'is_from_ligand_specs': True
            }
        
        # If not in ligand_specs, analyze based on PDB content
        chain_type_classification, _ = classify_chain(chain_pdb_obj, verbose=verbose)
        is_nucleic_from_pdb = (chain_type_classification == 'nucleic')
        is_protein_from_pdb = (chain_type_classification == 'protein')
        
        # Build sequence
        current_chain_sequence_chars = []
        for residue in chain_pdb_obj:
            if residue.id[0] == ' ':
                resname = residue.resname
                if is_nucleic_from_pdb:
                    res_code = residue_constants.resname_to_nucleic_acid_1letter.get(resname, "X")
                else:
                    res_code = residue_constants.restype_3to1.get(resname, "X")
                current_chain_sequence_chars.append(res_code)
        
        sequence_to_use = "".join(current_chain_sequence_chars)
        
        # Determine chain type
        if is_protein_from_pdb:
            chain_type_to_use = "protein"
        elif is_nucleic_from_pdb:
            is_rna_type = False
            for residue in chain_pdb_obj:
                if residue.id[0] == ' ':
                    if residue.resname in residue_constants.resname_to_nucleic_acid_1letter:
                        if residue_constants.resname_to_nucleic_acid_1letter[residue.resname] == 'U':
                            is_rna_type = True
                            break
                    if "O2'" in residue:
                        is_rna_type = True
                        break
            chain_type_to_use = "rna" if is_rna_type else "dna"
        else:
            if verbose:
                print(f"Warning: Unable to determine chain type for {current_chain_id} from PDB. Skipping this chain.")
            return None
        
        if not sequence_to_use:
            if verbose:
                print(f"Warning: Chain {current_chain_id} has empty sequence or no standard residue. Skipping.")
            return None
        
        return {
            'chain_id': current_chain_id,
            'sequence': sequence_to_use,
            'chain_type': chain_type_to_use,
            'is_from_ligand_specs': False
        }

    def _process_chain_masks(self, chain_pdb_obj, chain_info, use_pocket_msa, use_hotspot_msa,
                           residue_indices_map_global, flat_pocket_mask_global, flat_hotspot_mask_global,
                           attempt_custom_msa_map):
        """
        Process chain's pocket and hotspot masks
        
        Args:
            chain_pdb_obj: PDB chain object
            chain_info: Chain information dictionary
            use_pocket_msa: Whether to use pocket MSA
            use_hotspot_msa: Whether to use hotspot MSA
            residue_indices_map_global: Global residue index mapping
            flat_pocket_mask_global: Global pocket mask
            flat_hotspot_mask_global: Global hotspot mask
            attempt_custom_msa_map: Whether to attempt custom MSA mapping
            
        Returns:
            dict: Dictionary containing mask information
                - pocket_mask_values: Pocket mask value list
                - hotspot_mask_values: Hotspot mask value list
        """
        current_chain_pocket_mask_values = []
        current_chain_hotspot_mask_values = []
        
        # Only process masks for protein chains from PDB
        if (chain_info['chain_type'] == 'protein' and 
            not chain_info['is_from_ligand_specs'] and 
            residue_indices_map_global):
            
            for residue in chain_pdb_obj:
                if residue.id[0] == ' ':
                    pdb_res_key = (chain_pdb_obj.id, residue.id[1])
                    
                    if pdb_res_key in residue_indices_map_global:
                        idx_in_flat = residue_indices_map_global[pdb_res_key]
                        
                        if use_pocket_msa and flat_pocket_mask_global is not None:
                            current_chain_pocket_mask_values.append(flat_pocket_mask_global[idx_in_flat])
                        elif use_pocket_msa:
                            current_chain_pocket_mask_values.append(0)
                        
                        if use_hotspot_msa and flat_hotspot_mask_global is not None:
                            current_chain_hotspot_mask_values.append(flat_hotspot_mask_global[idx_in_flat])
                        elif use_hotspot_msa:
                            current_chain_hotspot_mask_values.append(0)
                    
                    elif attempt_custom_msa_map:
                        if self.verbose:
                            print(f"Warning: Custom MSA - residue {residue.id[1]} of protein chain {chain_pdb_obj.id} "
                                  f"not found in global mask mapping.")
                        if use_pocket_msa:
                            current_chain_pocket_mask_values.append(0)
                        if use_hotspot_msa:
                            current_chain_hotspot_mask_values.append(0)
        
        return {
            'pocket_mask_values': current_chain_pocket_mask_values,
            'hotspot_mask_values': current_chain_hotspot_mask_values
        }

    def _create_chain_with_msa(self, chain_info, mask_info, use_af3_msa, use_pocket_msa, use_hotspot_msa, custom_msa_data, current_ligand_specs=None, use_pocket_masked_af3_msa=False):
        """
        Create folding_input chain object based on chain information and mask information
        
        Args:
            chain_info: Chain information dictionary
            mask_info: Mask information dictionary
            use_af3_msa: Whether to use AF3 MSA
            use_pocket_msa: Whether to use pocket MSA
            use_hotspot_msa: Whether to use hotspot MSA
            custom_msa_data: Custom MSA data dictionary
            current_ligand_specs: ligand_specs configuration dictionary
            use_pocket_masked_af3_msa: Whether to use pocket masked AF3 MSA

        Returns:
            folding_input chain object, if creation fails, return None
        """
        chain_id = chain_info['chain_id']
        chain_type = chain_info['chain_type']
        
        if chain_type in ["rna", "dna"]:
            return self._create_nucleic_chain(chain_info, use_af3_msa or use_pocket_masked_af3_msa, current_ligand_specs)
        elif chain_type == "protein":
            return self._create_protein_chain(chain_info, mask_info, use_af3_msa, use_pocket_msa, use_hotspot_msa, custom_msa_data, use_pocket_masked_af3_msa)
        elif chain_type == "small_molecule":
            return self._create_small_molecule_chain(chain_info)
        else:
            if self.verbose:
                print(f"Warning: Unknown chain type '{chain_type}' for chain {chain_id}. Skipping this chain.")
            return None

    def _create_nucleic_chain(self, chain_info, use_af3_msa, current_ligand_specs=None):
        """Create nucleic acid chain object"""
        chain_id = chain_info['chain_id']
        sequence = chain_info['sequence']
        chain_type = chain_info['chain_type']
        is_from_ligand_specs = chain_info['is_from_ligand_specs']
        
        msa_to_pass = ""
        
        if is_from_ligand_specs and current_ligand_specs:
            json_msa_preference = current_ligand_specs.get(chain_id, {}).get(
                "unpaired_msa_preference", "KEY_WAS_LEGITIMATELY_MISSING"
            )
            
            if json_msa_preference == "" or json_msa_preference is None:
                msa_to_pass = ""
                if self.verbose:
                    print(f"Ligand chain {chain_id}: JSON specifies unpaired_msa as empty/null. MSA search will be skipped.")
            elif json_msa_preference == "KEY_WAS_LEGITIMATELY_MISSING":
                if use_af3_msa:
                    msa_to_pass = None
                    if self.verbose:
                        print(f"Ligand chain {chain_id}: JSON has no unpaired_msa key. Global MSA is enabled, AF3 will search MSA.")
                else:
                    msa_to_pass = ""
                    if self.verbose:
                        print(f"Ligand chain {chain_id}: JSON has no unpaired_msa key. Global MSA is disabled, MSA search will be skipped.")
            else:
                msa_to_pass = json_msa_preference
                if self.verbose:
                    print(f"Ligand chain {chain_id}: JSON provides a specific unpaired_msa string.")
        else:
            # PDB chain, based on global settings
            msa_to_pass = None if use_af3_msa else ""
            if self.verbose:
                print(f"PDB chain {chain_id} (type {chain_type}): MSA behavior based on global use_af3_msa ({use_af3_msa}).")
        
        # Create chain object
        if chain_type == "rna":
            return folding_input.RnaChain(
                id=chain_id, sequence=sequence, modifications=[], 
                unpaired_msa=msa_to_pass
            )
        else:  # DNA
            return folding_input.DnaChain(
                id=chain_id, sequence=sequence, modifications=[]
            )

    def _create_protein_chain(self, chain_info, mask_info, use_af3_msa, use_pocket_msa, use_hotspot_msa, custom_msa_data, use_pocket_masked_af3_msa=False):
        """Create protein chain object"""
        chain_id = chain_info['chain_id']
        sequence = chain_info['sequence']
        pocket_mask_values = mask_info['pocket_mask_values']
        hotspot_mask_values = mask_info['hotspot_mask_values']
        
        # Generate custom MSA sequence
        final_unpaired_msa_str = None
        pocket_msa_seq_generated = False
        hotspot_msa_seq_generated = False
        
        # Try to generate pocket MSA part
        if (use_pocket_msa and pocket_mask_values and 
            len(sequence) == len(pocket_mask_values)):
            pocket_msa_seq_parts = [
                s_char if mask_val == 1 else '-'
                for s_char, mask_val in zip(sequence, pocket_mask_values)
            ]
            pocket_msa_final_seq = "".join(pocket_msa_seq_parts)
            pocket_msa_seq_generated = True
            
            if self.verbose:
                print(f"Generated pocket MSA sequence part for protein chain {chain_id}.")
                print(f"Pocket MSA sequence first 20 characters: {pocket_msa_final_seq[:20]}...")
        elif use_pocket_msa and self.verbose:
            print(f"Warning: Failed to generate pocket MSA sequence part for protein chain {chain_id}.")
        
        # Try to generate hotspot MSA part
        if (use_hotspot_msa and hotspot_mask_values and 
            len(sequence) == len(hotspot_mask_values)):
            hotspot_msa_seq_parts = [
                s_char if mask_val == 1 else '-'
                for s_char, mask_val in zip(sequence, hotspot_mask_values)
            ]
            hotspot_msa_final_seq = "".join(hotspot_msa_seq_parts)
            hotspot_msa_seq_generated = True
            
            if self.verbose:
                print(f"Generated hotspot MSA sequence part for protein chain {chain_id}.")
                print(f"Hotspot MSA sequence first 20 characters: {hotspot_msa_final_seq[:20]}...")
        elif use_hotspot_msa and self.verbose:
            print(f"Warning: Failed to generate hotspot MSA sequence part for protein chain {chain_id}.")
        
        # Build MSA
        msa_builder = [f">seq1\n{sequence}"]
        current_seq_idx = 2
        
        if pocket_msa_seq_generated:
            msa_builder.append(f">seq{current_seq_idx}\n{pocket_msa_final_seq}")
            current_seq_idx += 1
        
        if hotspot_msa_seq_generated:
            msa_builder.append(f">seq{current_seq_idx}\n{hotspot_msa_final_seq}")
            current_seq_idx += 1
        
        # Decide final MSA
        if len(msa_builder) > 1:
            final_unpaired_msa_str = "\n".join(msa_builder)
        elif len(msa_builder) == 1 and (use_pocket_msa or use_hotspot_msa) and not use_af3_msa:
            final_unpaired_msa_str = msa_builder[0]
        
        # Set MSA strategy
        current_unpaired_msa = final_unpaired_msa_str
        current_paired_msa = None
        current_templates = []
        
        if current_unpaired_msa is not None:
            if use_af3_msa or use_pocket_masked_af3_msa:
                custom_msa_data[chain_id] = current_unpaired_msa
                if self.verbose:
                    print(f"Compatible mode: custom MSA for protein chain {chain_id} has been saved, will be merged after AF3 pipeline.")
                current_unpaired_msa = None
                current_paired_msa = None
                current_templates = None
            else:
                current_paired_msa = ""
                current_templates = []
                if self.verbose:
                    print(f"Custom MSA mode: protein chain {chain_id} only uses custom MSA, no AF3 pipeline.")
        else:
            if use_af3_msa or use_pocket_masked_af3_msa:
                current_unpaired_msa = None
                current_paired_msa = None
                current_templates = None
                if self.verbose:
                    if use_pocket_masked_af3_msa:
                        print(f"Pocket masked AF3 MSA mode: protein chain {chain_id} will be generated by AF3 pipeline then pocket masked.")
                    else:
                        print(f"AF3 MSA mode: protein chain {chain_id} will be generated by AF3 pipeline.")
            else:
                current_unpaired_msa = ""
                current_paired_msa = ""
                current_templates = []
                if self.verbose:
                    print(f"No MSA mode: protein chain {chain_id} does not use any MSA.")
        
        return folding_input.ProteinChain(
            id=chain_id, sequence=sequence, ptms=[],
            paired_msa=current_paired_msa, 
            unpaired_msa=current_unpaired_msa,
            templates=current_templates
        )

    def _create_small_molecule_chain(self, chain_info):
        """Create small molecule chain object"""
        chain_id = chain_info['chain_id']
        smiles = chain_info['sequence']
        
        if smiles:
            if self.verbose:
                print(f"Add small molecule chain {chain_id} (SMILES: {smiles})")
            return folding_input.Ligand(id=chain_id, smiles=smiles)
        else:
            if self.verbose:
                print(f"Warning: Small molecule chain {chain_id} is missing SMILES string.")
            return None

    def _create_small_molecule_chain_from_specs(self, spec_chain_id, spec_info):
        """
        Create small molecule chain object from ligand_specs configuration
        
        Args:
            spec_chain_id: Chain ID
            spec_info: Configuration information in ligand_specs
            
        Returns:
            folding_input.Ligand object, if creation fails, return None
        """
        spec_smiles = spec_info.get("smiles")
        spec_ccd_codes = spec_info.get("ccd_codes")
        
        if spec_ccd_codes:
            if self.verbose:
                print(f"Add small molecule chain from ligand_specs: ID={spec_chain_id}, CCD codes={spec_ccd_codes}")
            return folding_input.Ligand(id=spec_chain_id, ccd_ids=spec_ccd_codes)
        elif spec_smiles:
            if self.verbose:
                print(f"Add small molecule chain from ligand_specs: ID={spec_chain_id}, SMILES={spec_smiles}")
            return folding_input.Ligand(id=spec_chain_id, smiles=spec_smiles)
        else:
            if self.verbose:
                print(f"Warning: Small molecule chain {spec_chain_id} in ligand_specs is missing CCD codes or SMILES string. Skipping.")
            return None


# def test_seq_emb_af3(pdb_file, weight_dir=None, verbose=True):
#     """
#     测试函数：生成序列嵌入并打印统计信息
    
#     Args:
#         pdb_file: PDB文件路径
#         weight_dir: 权重目录
#         verbose: 是否显示详细输出
        
#     Returns:
#         bool: 成功为True，失败为False
#     """
#     try:
#         # 初始化嵌入器
#         embedder = Embedder(weight_dir=weight_dir, verbose=verbose)
        
#         # 记录开始时间
#         start_time = time.time()
        
#         # 生成嵌入
#         data = embedder.seq_emb_af3(pdb_file)
        
#         # 计算运行时间
#         end_time = time.time()
#         elapsed = end_time - start_time
        
#         # 打印统计信息
#         if verbose:
#             print(f"生成嵌入用时: {elapsed:.2f}秒")
            
#             # 嵌入统计信息
#             if 'embeddings' in data and 'single' in data['embeddings']:
#                 emb = data['embeddings']['single']
#                 pair_emb = data['embeddings'].get('pair')
                
#                 print(f"嵌入形状: {emb.shape}")
#                 if pair_emb is not None:
#                     print(f"配对嵌入形状: {pair_emb.shape}")
                
#                 print(f"均值: {np.mean(emb):.4f}")
#                 print(f"标准差: {np.std(emb):.4f}")
#                 print(f"最小值: {np.min(emb):.4f}")
#                 print(f"最大值: {np.max(emb):.4f}")
                
#             # 序列信息
#             if 'seq_info' in data:
#                 seq_info = data['seq_info']
#                 print(f"序列名称: {seq_info['name']}")
#                 for i, chain in enumerate(seq_info['chains']):
#                     print(f"链 {i+1}: {chain['chain_id']}, 类型: {chain['chain_type']}, 长度: {len(chain['sequence'])}")
        
#         # 测试特定链的嵌入
#         if 'seq_info' in data and 'chains' in data['seq_info'] and len(data['seq_info']['chains']) > 0:
#             chain_id = data['seq_info']['chains'][0]['chain_id']
#             try:
#                 chain_data = embedder.seq_emb_af3(pdb_file, chain_id=chain_id)
#                 if verbose:
#                     print(f"成功生成链 {chain_id} 的嵌入 (形状: {chain_data['embeddings']['single'].shape})")
#             except Exception as e:
#                 print(f"生成链 {chain_id} 的嵌入时出错: {e}")
#                 import traceback
#                 traceback.print_exc()
#                 return False
        
#         return True
#     except Exception as e:
#         print(f"嵌入生成失败: {e}")
#         import traceback
#         traceback.print_exc()
#         return False


# def test_struct_emb_af3(pdb_file=None, weight_dir=None, verbose=True):
#     """
#     测试结构嵌入生成
    
#     Args:
#         pdb_file: PDB文件路径，默认为None，使用测试PDB
#         weight_dir: 权重目录路径
#         verbose: 是否显示详细输出
        
#     Returns:
#         bool: 成功为True，失败为False
#     """
#     try:
#         # 设置默认测试PDB文件
#         if pdb_file is None:
#             pdb_file = "/work/hat170/aptamer/test_input/1A1V.pdb"
            
#         if not os.path.exists(pdb_file):
#             print(f"错误: 测试PDB文件 {pdb_file} 不存在")
#             return False
        
#         # 初始化嵌入器
#         embedder = Embedder(weight_dir=weight_dir, verbose=verbose)
        
#         # 记录开始时间
#         start_time = time.time()
        
#         # 生成结构嵌入
#         data = embedder.struct_emb_af3(
#             pdb_file=pdb_file,
#             db_dir=None,
#             msa_dir=None,
#             use_af3_msa=False, # 测试时通常禁用 MSA 以加快速度
#             use_pocket_msa=True, # 测试启用口袋MSA
#             use_hotspot_msa=True, # 测试启用热点MSA
#             hotspot_cutoff=8.0,
#             pocket_cutoff=10.0
#         )
        
#         # 计算运行时间
#         elapsed = time.time() - start_time
        
#         # 打印统计信息
#         if verbose:
#             print(f"\n--- 结构嵌入生成用时: {elapsed:.2f}秒 ---")
            
#             # 基本嵌入信息
#             if 'embeddings' in data:
#                 if 'single' in data['embeddings']:
#                     print(f"序列嵌入形状: {data['embeddings']['single'].shape}")
#                 if 'pair' in data['embeddings']:
#                     print(f"配对嵌入形状: {data['embeddings']['pair'].shape}")
#                 if 'structure_atom_coords' in data['embeddings']:
#                     print(f"结构坐标嵌入形状: {data['embeddings']['structure_atom_coords'].shape}")
#                 if 'structure_atom_mask' in data['embeddings']:
#                     print(f"原子掩码形状: {data['embeddings']['structure_atom_mask'].shape}")
#                     print(f"原子掩码中为1的数量: {int(np.sum(data['embeddings']['structure_atom_mask']))}")
            
#             # 口袋掩码信息
#             if 'pocket_mask' in data:
#                 print(f"\n--- 口袋掩码信息 ---")
#                 pocket_mask = data['pocket_mask']
#                 pocket_count = np.sum(pocket_mask)
#                 total_count = len(pocket_mask)
#                 print(f"口袋掩码中标记为1的数量: {pocket_count} / {total_count} ({pocket_count/total_count*100:.1f}%)")
            
#             # 序列信息
#             if 'seq_info' in data:
#                 print(f"\n--- 序列信息 ---")
#                 for i, chain in enumerate(data['seq_info']['chains']):
#                     print(f"链 {i+1}: ID={chain['chain_id']}, 类型={chain['chain_type']}, 长度={len(chain['sequence'])}")
#                     # 显示序列前10个字符
#                     print(f"  序列前缀: {chain['sequence'][:10]}...")
        
#         # 保存结果用于后续分析
#         output_dir = "/work/hat170/aptamer/test_output"
#         os.makedirs(output_dir, exist_ok=True)
#         # 文件名反映是结构嵌入
#         output_filename = f"{os.path.basename(pdb_file).replace('.pdb', '')}_struct_emb.pkl"
#         output_path = os.path.join(output_dir, output_filename)
        
#         with open(output_path, 'wb') as f:
#             pickle.dump(data, f)
        
#         if verbose:
#             print(f"\n结构嵌入已保存到: {output_path}")
        
#         return True
#     except Exception as e:
#         print(f"结构嵌入生成失败: {e}")
#         import traceback
#         traceback.print_exc()
#         return False


# def create_parser():
#     """创建命令行参数解析器"""
#     parser = argparse.ArgumentParser(description='使用AlphaFold3生成序列或结构嵌入')
#     parser.add_argument('--pdb_file', type=str, help='要处理的PDB文件路径')
#     parser.add_argument('--output_dir', type=str, default='/work/hat170/aptamer/test_output', help='输出目录路径')
#     parser.add_argument('--weight_dir', type=str, default='/work/hat170/aptamer/alphafold3/weight', 
#                         help='AlphaFold3权重目录路径')
#     parser.add_argument('--hotspot_cutoff', type=float, default=8.0, help='热点残基的截断距离')
#     parser.add_argument('--pocket_cutoff', type=float, default=10.0, help='口袋残基的截断距离')
#     parser.add_argument('--generate_struct_emb', action='store_true', 
#                         help='生成结构嵌入（包含坐标和原子掩码），而不是仅生成序列嵌入')
#     parser.add_argument('--verbose', action='store_true', help='显示详细输出')
#     parser.add_argument('--db_dir', type=str, default="/work/hat170/aptamer/alphafold3/database/", help='MSA数据库目录路径')
#     parser.add_argument('--msa_dir', type=str, default=None, help='预计算的MSA文件目录路径')
#     parser.add_argument('--use_af3_msa', action='store_true', default=False, 
#                         help='是否使用AlphaFold3的MSA生成流程（例如，运行jackhmmer等），默认为False。不应与 --use_pocket_msa 或 --use_hotspot_msa 同时使用。')
#     parser.add_argument('--use_pocket_msa', action='store_true', default=False,
#                         help='是否为蛋白质链生成基于口袋的unpaired MSA，默认为False。不应与 --use_af3_msa 同时使用。')
#     parser.add_argument('--use_hotspot_msa', action='store_true', default=False,
#                         help='是否为蛋白质链生成基于热点的unpaired MSA，默认为False。不应与 --use_af3_msa 同时使用。')
    
#     return parser


# def main():
#     """
#     主函数：处理命令行参数并执行嵌入生成
#     """
#     parser = create_parser()
#     args = parser.parse_args()
    
#     # 验证参数
#     if not args.pdb_file:
#         parser.error("必须提供--pdb_file参数")
    
#     # 移除原有的互斥检查，允许MSA选项同时使用
#     # if args.use_af3_msa and (args.use_pocket_msa or args.use_hotspot_msa):
#     #     parser.error("--use_af3_msa 不能与 --use_pocket_msa 或 --use_hotspot_msa 同时指定。请选择一种MSA策略。")
    
#     # 添加新的兼容性提示
#     if args.use_af3_msa and (args.use_pocket_msa or args.use_hotspot_msa):
#         if args.verbose:
#             print("注意: 同时启用了AF3 MSA和自定义MSA。将先生成自定义MSA，然后由AF3管道补充更多序列。")
    
#     # 创建输出目录
#     os.makedirs(args.output_dir, exist_ok=True)
    
#     # 创建嵌入器
#     try:
#         embedder = Embedder(weight_dir=args.weight_dir, verbose=args.verbose)
#     except Exception as e:
#         print(f"初始化嵌入器失败: {e}")
#         import traceback
#         traceback.print_exc()
#         return 1
        
#     # 处理PDB文件
#     if not os.path.exists(args.pdb_file):
#         # 使用 FileNotFoundError 更合适
#         raise FileNotFoundError(f"PDB文件不存在: {args.pdb_file}")
        
#     pdb_name = os.path.basename(args.pdb_file).replace(".pdb", "")
    
#     try:
#         # 根据模式选择生成序列嵌入或结构嵌入
#         start_time = time.time()
        
#         if args.generate_struct_emb:
#             # 生成结构嵌入
#             output_filename = f"{pdb_name}_struct_emb.pkl"
#             output_path = os.path.join(args.output_dir, output_filename)
#             data = embedder.struct_emb_af3(
#                 pdb_file=args.pdb_file,
#                 db_dir=args.db_dir,
#                 msa_dir=args.msa_dir,
#                 use_af3_msa=args.use_af3_msa,
#                 use_pocket_msa=args.use_pocket_msa,
#                 use_hotspot_msa=args.use_hotspot_msa,
#                 hotspot_cutoff=args.hotspot_cutoff,
#                 pocket_cutoff=args.pocket_cutoff
#             )
        
#         else:
#             # 生成序列嵌入
#             output_filename = f"{pdb_name}_seq_emb.pkl" # 明确是序列嵌入
#             output_path = os.path.join(args.output_dir, output_filename)
#             data = embedder.seq_emb_af3(
#                 pdb_file=args.pdb_file, 
#                 use_af3_msa=args.use_af3_msa,
#                 use_pocket_msa=args.use_pocket_msa,
#                 use_hotspot_msa=args.use_hotspot_msa,
#                 hotspot_cutoff=args.hotspot_cutoff,
#                 pocket_cutoff=args.pocket_cutoff
#             )
        
#         elapsed = time.time() - start_time
        
#         # 保存结果
#         with open(output_path, 'wb') as f:
#             pickle.dump(data, f)
            
#         if args.verbose:
#             print(f"嵌入已保存到 {output_path}")
#             print(f"处理时间: {elapsed:.2f}秒")
            
#     except Exception as e:
#         print(f"处理PDB文件 '{args.pdb_file}' 时出错: {e}")
#         import traceback
#         traceback.print_exc()
#         # 考虑返回非零退出码或继续处理其他文件（如果适用）
#         return 1 # 返回错误码
    
#     return 0 # 成功完成


# if __name__ == "__main__":
#     # 使用 sys.exit 更标准
#     sys.exit(main())
