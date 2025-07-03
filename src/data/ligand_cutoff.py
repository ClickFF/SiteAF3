"""
Script to preprocess PDB datasets.
NOTE: This script will overwrite the raw data.
"""
import pyrootutils

# See: https://github.com/ashleve/pyrootutils
root = pyrootutils.setup_root(
    search_from=__file__,
    indicator=[".git"],
    pythonpath=True,
    dotenv=True
)

import os
import copy
import argparse
import dataclasses
import pandas as pd
import numpy as np
import mdtraj as md
import time
import multiprocessing as mp
import functools as fn
from tqdm import tqdm
from Bio import PDB
from Bio.PDB import Select
import pathlib

import torch
from torch import nn as nn
from torch.utils.data import DataLoader

from src.data import utils as du
from src.data import errors
from src.data import parsers
from src.data import residue_constants


def create_parser():
    parser = argparse.ArgumentParser(
        description='PDB processing script.'
    )

    parser.add_argument(
        '--pdb_dir',
        help='Path to directory with PDB files.',
        type=str
    )

    parser.add_argument(
        '--num_processes',
        help='Number of processes.',
        type=int,
        default=16
    )

    parser.add_argument(
        '--write_dir',
        help='Path to write results to.',
        type=str,
        default=os.path.join(os.getcwd(), 'processed_pdbs')
    )

    parser.add_argument(
        '--hotspot_cutoff',
        help='Cutoff for hotspot residues.',
        type=int,
        default=8
    )

    parser.add_argument(
        '--pocket_cutoff',
        help='Cutoff for pocket residues.',
        type=int,
        default=10
    )

    parser.add_argument(
        '--verbose',
        help='Whether to log everything.',
        action='store_true'
    )

    return parser


def resid_unique(res):
    if type(res) == str:
        res_ = res.split()
        return f'{res_[1]}_{res_[2]}'
    return f'{res.id[1]}_{res.get_parent().id}'


def get_seq(entity, lig_chains):
    """
    Extract sequence information from entity, distinguishing ligand and receptor chains
    
    Args:
        entity: PDB structure entity
        lig_chains: List of ligand chain IDs
        
    Returns:
        aatype_full: Amino acid/nucleotide types for all residues
        mask_full: Mask marking ligand (1) and receptor (0)
        chain_id_full: Chain IDs for all residues
    """
    aatype_lig, aatype_rec = [], []
    chain_id_lig, chain_id_rec = [], []
    
    for res in entity.get_residues():
        res_name = residue_constants.substitute_non_standard_restype.get(res.resname, res.resname)
        try:
            float(res_name)
            raise errors.DataError(f"Residue name should not be a number: {res.resname}")
        except ValueError:
            pass
        
        res_shortname = residue_constants.restype_3to1.get(res_name, '<unk>')
        
        if res.parent.id in lig_chains:
            aatype_lig.append(res_shortname)
            chain_id_lig.append(res.parent.id)
        else:
            aatype_rec.append(res_shortname)
            chain_id_rec.append(res.parent.id)
            
    aatype_full = aatype_lig + aatype_rec
    chain_id_full = chain_id_lig + chain_id_rec
    mask_full = np.concatenate([np.ones(len(aatype_lig)), np.zeros(len(aatype_rec))])
    
    assert len(aatype_full) == len(mask_full) == len(chain_id_full), "Shape mismatch when processing raw data!"
    return aatype_full, mask_full, chain_id_full


class resSelector(Select):
    def __init__(self, res_ids):
        self.res_ids = res_ids
    def accept_residue(self, residue):
        resid = resid_unique(residue)
        if resid not in self.res_ids:
            return False
        else:
            return True


def is_nucleic(residue):
    """Determine if residue is a nucleic acid residue (DNA or RNA)"""
    standard_na_residues = {'A', 'C', 'G', 'U', 'T', 'DA', 'DC', 'DG', 'DT', 'DU'}
    return residue.resname in standard_na_residues or residue.resname.startswith('D')


def is_small_molecule(residue):
    """Determine if residue is a small molecule ligand"""
    excluded_molecules = {
        'HOH', 'WAT', 'H2O', 'TIP', 'TIP3', 'TIP4', 'SPC',
        'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU',
        'SO4', 'PO4', 'NO3', 'BR', 'I', 'F',
        'ACE', 'NME',
        'MSE', 'SEP', 'TPO', 'PTR',
        'EDO', 'PEG', 'GOL', 'ETG'
    }
    
    resname = residue.resname.strip()
    
    if resname in excluded_molecules:
        return False
    
    if resname in residue_constants.restype_3to1:
        return False
        
    if resname in residue_constants.resname_to_nucleic_acid_1letter:
        return False
    
    if is_nucleic(residue):
        return False
    
    num_atoms = len(list(residue.get_atoms()))
    
    if num_atoms < 3:
        return False
        
    if num_atoms > 200:
        return False
    
    if residue.id[0] != ' ':
        return True
    
    atom_types = {atom.element for atom in residue.get_atoms() if atom.element}
    common_protein_atoms = {'C', 'N', 'O', 'S'}
    
    if atom_types - common_protein_atoms:
        return True
    
    return False


def classify_chain(chain, verbose=False):
    """
    Classify chain as protein, nucleic acid, or small molecule using robust detection mechanism
    
    Args:
        chain: PDB chain object
        verbose: Whether to output detailed information
        
    Returns:
        chain_type: 'protein', 'nucleic', 'small_molecule', or 'unknown'
        classification_source: Classification basis
    """
    standard_residues = [res for res in chain if res.id[0] == ' ']
    hetero_residues = [res for res in chain if res.id[0] != ' ']
    
    if not standard_residues and not hetero_residues:
        return 'unknown', 'No residues found'
    
    if not standard_residues and hetero_residues:
        small_mol_count = sum(1 for res in hetero_residues if is_small_molecule(res))
        if small_mol_count > 0:
            return 'small_molecule', f'HETATM only, {small_mol_count} small molecules'
        else:
            return 'unknown', 'HETATM only but no recognizable small molecules'
    
    is_nucleic_result, is_protein_result, classification_source = classify_chain_legacy(chain, verbose)
    
    if is_nucleic_result:
        return 'nucleic', classification_source
    elif is_protein_result:
        return 'protein', classification_source
    else:
        all_residues = standard_residues + hetero_residues
        small_mol_count = sum(1 for res in all_residues if is_small_molecule(res))
        if small_mol_count / len(all_residues) > 0.5:
            return 'small_molecule', f'Majority small molecules ({small_mol_count}/{len(all_residues)})'
        else:
            return 'unknown', 'Could not determine chain type'


def classify_chain_legacy(chain, verbose=False):
    """
    Classify chain as protein or nucleic acid using robust detection mechanism (legacy function for backward compatibility)
    
    Args:
        chain: PDB chain object
        verbose: Whether to output detailed information
        
    Returns:
        is_nucleic: Whether it is a nucleic acid chain
        is_protein: Whether it is a protein chain
        classification_source: Classification basis
    """
    is_nucleic = False
    is_protein = False
    classification_source = "Unknown"
    standard_residues_in_chain = [res for res in chain if res.id[0] == ' ']
    num_standard_residues = len(standard_residues_in_chain)

    if num_standard_residues > 0:
        nucleic_evidence_count = 0
        protein_evidence_count = 0
        atom_check_inconclusive = False
        
        for residue in standard_residues_in_chain[:5]:
            resname = residue.resname
            try:
                has_P = residue.has_id('P')
                has_N9 = residue.has_id('N9')
                has_N1_pyrimidine = residue.has_id('N1') and not residue.has_id('N') and resname in ['U', 'C', 'T']
                has_CA = residue.has_id('CA')

                is_potential_nucleic = (has_P or has_N9 or has_N1_pyrimidine) and not has_CA
                is_potential_protein = has_CA and not (has_P or has_N9 or has_N1_pyrimidine)

                if is_potential_nucleic:
                    nucleic_evidence_count += 1
                if is_potential_protein:
                    protein_evidence_count += 1
                
                if is_potential_nucleic and is_potential_protein:
                    atom_check_inconclusive = True
                    if verbose:
                        print(f"Atom check on {resname} ({residue.id}): Ambiguous (both Nuc and Prot features).")
                    break
            except Exception as e:
                if verbose:
                    print(f"Warning: Error during atom check for residue {residue.id} in chain {chain.id}: {e}. Skipping residue atom check.")

        if not atom_check_inconclusive:
            if nucleic_evidence_count > 0 and protein_evidence_count == 0:
                is_nucleic = True
                classification_source = "Atom Check (Nucleic Evidence Only)"
            elif protein_evidence_count > 0 and nucleic_evidence_count == 0:
                is_protein = True
                classification_source = "Atom Check (Protein Evidence Only)"
            
        if not is_nucleic and not is_protein:
            standard_amino_count = sum(1 for res in standard_residues_in_chain if res.resname in residue_constants.restype_3to1)
            standard_nucleic_count = sum(1 for res in standard_residues_in_chain if res.resname in residue_constants.resname_to_nucleic_acid_1letter)
            
            if num_standard_residues > 0:
                nucleic_ratio = standard_nucleic_count / num_standard_residues
                amino_ratio = standard_amino_count / num_standard_residues

                if nucleic_ratio > 0.7:
                    is_nucleic = True
                    classification_source = "Ratio Check (>70% Nucleic)"
                elif amino_ratio > 0.7:
                    is_protein = True
                    classification_source = "Ratio Check (>70% Protein)"
            
        if not is_nucleic and not is_protein:
            if 'standard_amino_count' not in locals():
                standard_amino_count = sum(1 for res in standard_residues_in_chain if res.resname in residue_constants.restype_3to1)
                standard_nucleic_count = sum(1 for res in standard_residues_in_chain if res.resname in residue_constants.resname_to_nucleic_acid_1letter)
                
            if standard_nucleic_count > standard_amino_count:
                is_nucleic = True
                classification_source = "Count Check (Nucleic > Protein)"
            elif standard_amino_count > standard_nucleic_count:
                is_protein = True
                classification_source = "Count Check (Protein > Nucleic)"
                
        if not is_nucleic and not is_protein:
            is_protein = True
            classification_source = "Default (Protein)"
            if verbose:
                print(f"Warning: Chain {chain.id} classification highly ambiguous. Defaulting to Protein.")
                
    else:
        if verbose:
            print(f"Warning: Chain {chain.id} has no standard residues. Cannot determine type.")
        classification_source = "No Standard Residues"
         
    if verbose and classification_source != "Unknown" and classification_source != "No Standard Residues":
        type_str = "Nucleic" if is_nucleic else ("Protein" if is_protein else "Undetermined")
        print(f"Chain {chain.id}: Classified as {type_str} based on: {classification_source}")
        
    return is_nucleic, is_protein, classification_source


def get_representative_atom(residue, molecule_type=None):
    """
    Get representative atom coordinates for residue based on molecule type
    
    Args:
        residue: Biopython residue object
        molecule_type: Molecule type ('protein', 'nucleic', 'small_molecule'), optional
        
    Returns:
        np.ndarray: Coordinates of representative atom
    """
    if molecule_type is None:
        if is_nucleic(residue):
            molecule_type = 'nucleic'
        elif residue.resname in residue_constants.restype_3to1:
            molecule_type = 'protein'
        elif is_small_molecule(residue):
            molecule_type = 'small_molecule'
    else:
            molecule_type = 'protein'
    
    if molecule_type == 'nucleic':
        for atom_name in ['P', "C4'", 'C4*', "C3'", "C5'", "C2'", "C1'"]:
            if atom_name in residue:
                return residue[atom_name].get_coord()
    elif molecule_type == 'protein':
        if "CA" in residue:
            return residue["CA"].get_coord()
    elif molecule_type == 'small_molecule':
        coords = []
        for atom in residue.get_atoms():
            if not atom.get_name().startswith('H'):
                coords.append(atom.get_coord())
        if coords:
            return np.mean(coords, axis=0)
        
    for atom in residue.get_atoms():
        return atom.get_coord()
    
    return np.array([0.0, 0.0, 0.0])


def get_seq(entity, lig_chains, ligand_type=None):
    """
    Get sequence information based on ligand type
    
    Args:
        entity: PDB structure entity
        lig_chains: List of ligand chain IDs
        ligand_type: Ligand type, optional
        
    Returns:
        aatype_full: Types of all residues
        mask_full: Mask marking ligand (1) and receptor (0)
        chain_id_full: Chain IDs of all residues
    """
    aatype_lig, aatype_rec = [], []
    chain_id_lig, chain_id_rec = [], []
    
    for res in entity.get_residues():
        if ligand_type == 'small_molecule' and res.parent.id in lig_chains:
            res_name = res.resname
            res_shortname = res_name
            aatype_lig.append(res_shortname)
            chain_id_lig.append(res.parent.id)
        elif res.id[0] == ' ':
            res_name = residue_constants.substitute_non_standard_restype.get(res.resname, res.resname)
            try:
                float(res_name)
                raise errors.DataError(f"Residue name should not be a number: {res.resname}")
            except ValueError:
                pass
            
            res_shortname = residue_constants.restype_3to1.get(res_name, '<unk>')
            
            if res.parent.id in lig_chains:
                aatype_lig.append(res_shortname)
                chain_id_lig.append(res.parent.id)
            else:
                aatype_rec.append(res_shortname)
                chain_id_rec.append(res.parent.id)
            
    aatype_full = aatype_lig + aatype_rec
    chain_id_full = chain_id_lig + chain_id_rec
    mask_full = np.concatenate([np.ones(len(aatype_lig)), np.zeros(len(aatype_rec))])
    
    assert len(aatype_full) == len(mask_full) == len(chain_id_full), "Shape mismatch when processing raw data!"
    return aatype_full, mask_full, chain_id_full


def get_motif_center_pos(infile: str, receptor_type: str = 'protein', ligand_type: str = 'nucleic', 
                        hotspot_cutoff=8, pocket_cutoff=10, verbose=False, 
                        explicit_receptor_chain_ids: list = None, 
                        explicit_ligand_chain_ids: list = None):
    """
    Calculate complex center position supporting different receptor and ligand type combinations
    
    Args:
        infile: PDB file path
        receptor_type: Receptor type ('protein', 'nucleic')
        ligand_type: Ligand type ('protein', 'nucleic', 'small_molecule')
        hotspot_cutoff: Cutoff distance for hotspot residues (Angstrom)
        pocket_cutoff: Cutoff distance for pocket residues (Angstrom)
        verbose: Whether to output detailed information
        explicit_receptor_chain_ids: Optional, explicitly specify receptor chain ID list
        explicit_ligand_chain_ids: Optional, explicitly specify ligand chain ID list
        
    Returns:
        struct: Processed structure
        center_pos: Calculated center position
        raw_seq_data: Sequence data dictionary
    """
    p = PDB.PDBParser(QUIET=1)
    struct_orig = p.get_structure('original_structure', infile)[0]
    struct_to_modify = p.get_structure('structure_to_modify', infile)[0]
    
    receptor_chains = []
    ligand_chains = []
    
    if explicit_receptor_chain_ids is not None and explicit_ligand_chain_ids is not None:
        receptor_chains = list(explicit_receptor_chain_ids)
        ligand_chains = list(explicit_ligand_chain_ids)
        if verbose:
            print(f"Using explicitly specified receptor chains: {receptor_chains}")
            print(f"Using explicitly specified ligand chains: {ligand_chains}")
    else:
        if verbose:
            print(f"Auto-detecting receptor type: {receptor_type}, ligand type: {ligand_type}")
        
        for chain in struct_to_modify.get_chains():
            chain_residues = list(chain.get_residues())
            if not chain_residues:
                continue
                
            chain_type, classification_source = classify_chain(chain, verbose=False)
            
            if chain_type == receptor_type:
                receptor_chains.append(chain.id)
            elif chain_type == ligand_type:
                ligand_chains.append(chain.id)
            else:
                if verbose:
                    print(f"Warning: Chain {chain.id} classified as {chain_type}, does not match specified receptor or ligand type")
    
    if not receptor_chains:
        raise errors.DataError(f"No {receptor_type} chains found in {os.path.basename(infile)}")
    if not ligand_chains:
        raise errors.DataError(f"No {ligand_type} chains found in {os.path.basename(infile)}")
    
    seq_full, mask_full, chain_id_full = get_seq(struct_orig, ligand_chains, ligand_type)
    
    lig_residues = []
    for res in struct_orig.get_residues():
        if res.parent.id in ligand_chains:
            if ligand_type == 'small_molecule' or res.id[0] == ' ':
                lig_residues.append(res)
    
    rec_residues = []
    for res in struct_orig.get_residues():
        if res.parent.id in receptor_chains:
            if res.id[0] == ' ':
                rec_residues.append(res)
    
    ligand_repr_coords = [get_representative_atom(res, ligand_type) for res in lig_residues]
    
    motif_res_ids = [resid_unique(res) for res in lig_residues]
    
    hotspot_coords = []
    
    for lig_coord in ligand_repr_coords:
        for rec_res in rec_residues:
            rec_coord = get_representative_atom(rec_res, receptor_type)
            if np.linalg.norm(rec_coord - lig_coord) <= hotspot_cutoff:
                hotspot_coords.append(rec_coord)
                res_id = resid_unique(rec_res)
                if res_id not in motif_res_ids:
                    motif_res_ids.append(res_id)
                    idx = len(lig_residues) + rec_residues.index(rec_res)
                    if idx < len(mask_full):
                        mask_full[idx] = 1
    
    for hotspot_coord in hotspot_coords:
        for rec_res in rec_residues:
            res_id = resid_unique(rec_res)
            if res_id in motif_res_ids:
                continue
                
            rec_coord = get_representative_atom(rec_res, receptor_type)
            if np.linalg.norm(rec_coord - hotspot_coord) <= pocket_cutoff:
                if res_id not in motif_res_ids:
                    motif_res_ids.append(res_id)
                    idx = len(lig_residues) + rec_residues.index(rec_res)
                    if idx < len(mask_full):
                        mask_full[idx] = 1
    
    center_pos = np.sum(np.array(ligand_repr_coords), axis=0) / len(ligand_repr_coords)
    
    new_model = PDB.Model.Model(0)
    for ch_id in struct_to_modify.child_dict.keys():
        orig_chain = struct_to_modify[ch_id]
        new_chain = PDB.Chain.Chain(orig_chain.id)
        has_residues_in_chain = False
        for res in orig_chain:
            if resid_unique(res) in motif_res_ids:
                new_chain.add(copy.deepcopy(res))
                has_residues_in_chain = True
        if has_residues_in_chain:
            new_model.add(new_chain)
    
    final_struct = PDB.Structure.Structure("motif_structure")
    final_struct.add(new_model)

    raw_seq_data = {
        "receptor_chains": receptor_chains,
        "ligand_chains": ligand_chains,
        "receptor_type": receptor_type,
        "ligand_type": ligand_type
    }
    
    seq_full_orig, mask_full_orig, chain_id_full_orig = get_seq(struct_orig, ligand_chains, ligand_type)
    for i, chain_id in enumerate(chain_id_full_orig):
        if chain_id not in raw_seq_data:
            raw_seq_data[chain_id] = {"seq": "", "mask": []}
        raw_seq_data[chain_id]["seq"] += seq_full_orig[i]
        raw_seq_data[chain_id]["mask"].append(mask_full_orig[i])
    
    return final_struct, center_pos, raw_seq_data


def get_hotspot_complex_struct(infile: str, receptor_type: str = 'protein', ligand_type: str = 'nucleic',
                              hotspot_cutoff=8, verbose=False, 
                              explicit_receptor_chain_ids: list = None, 
                              explicit_ligand_chain_ids: list = None):
    """
    Calculate hotspot complex structure supporting different receptor and ligand type combinations
    
    Args:
        infile: PDB file path
        receptor_type: Receptor type ('protein', 'nucleic')
        ligand_type: Ligand type ('protein', 'nucleic', 'small_molecule')
        hotspot_cutoff: Cutoff distance for hotspot residues (Angstrom)
        verbose: Whether to output detailed information
        explicit_receptor_chain_ids: Optional, explicitly specify receptor chain ID list
        explicit_ligand_chain_ids: Optional, explicitly specify ligand chain ID list
        
    Returns:
        struct: Processed structure containing only ligand and hotspot residues
    """
    p = PDB.PDBParser(QUIET=1)
    struct_to_modify = p.get_structure('structure_to_modify_hotspot', infile)[0]
    
    receptor_chains = []
    ligand_chains = []
    
    if explicit_receptor_chain_ids is not None and explicit_ligand_chain_ids is not None:
        receptor_chains = list(explicit_receptor_chain_ids)
        ligand_chains = list(explicit_ligand_chain_ids)
        if verbose:
            print(f"(Hotspot) Using explicitly specified receptor chains: {receptor_chains}")
            print(f"(Hotspot) Using explicitly specified ligand chains: {ligand_chains}")
    else:
        if verbose:
            print(f"(Hotspot) Auto-detecting receptor type: {receptor_type}, ligand type: {ligand_type}")
        
        for chain in struct_to_modify.get_chains():
            chain_residues = list(chain.get_residues())
            if not chain_residues:
                continue
                
            chain_type, classification_source = classify_chain(chain, verbose=False)
            
            if chain_type == receptor_type:
                receptor_chains.append(chain.id)
            elif chain_type == ligand_type:
                ligand_chains.append(chain.id)

    if not ligand_chains:
        raise errors.DataError(f"(Hotspot) No {ligand_type} chains found in {os.path.basename(infile)}")

    # 收集配体和受体残基
    lig_residues = []
    for res in struct_to_modify.get_residues():
        if res.parent.id in ligand_chains:
            if ligand_type == 'small_molecule' or res.id[0] == ' ':
                lig_residues.append(res)
    
    rec_residues = []
    for res in struct_to_modify.get_residues():
        if res.parent.id in receptor_chains:
            if res.id[0] == ' ':
                rec_residues.append(res)
    
    ligand_repr_coords = [get_representative_atom(res, ligand_type) for res in lig_residues]
    
    hotspot_res_ids = [resid_unique(res) for res in lig_residues]  # 包含所有配体残基
    
    # 计算热点残基
    for lig_coord in ligand_repr_coords:
        for rec_res in rec_residues:
            rec_coord = get_representative_atom(rec_res, receptor_type)
            if np.linalg.norm(rec_coord - lig_coord) <= hotspot_cutoff:
                res_id = resid_unique(rec_res)
                if res_id not in hotspot_res_ids:
                    hotspot_res_ids.append(res_id)
    
    # 创建新结构
    new_model_hotspot = PDB.Model.Model(0)
    for ch_id in struct_to_modify.child_dict.keys():
        orig_chain = struct_to_modify[ch_id]
        new_chain = PDB.Chain.Chain(orig_chain.id)
        has_residues_in_chain = False
        for res in orig_chain:
            if resid_unique(res) in hotspot_res_ids:
                new_chain.add(copy.deepcopy(res))
                has_residues_in_chain = True
        if has_residues_in_chain:
            new_model_hotspot.add(new_chain)
    
    final_hotspot_struct = PDB.Structure.Structure("hotspot_structure")
    final_hotspot_struct.add(new_model_hotspot)
        
    return final_hotspot_struct


# 保持向后兼容性的包装函数
def get_motif_center_pos_legacy(infile:str, hotspot_cutoff=8, pocket_cutoff=10, verbose=False, 
                                explicit_protein_chain_ids: list = None, 
                                explicit_ligand_chain_ids: list = None):
    """
    Calculate protein-nucleic acid complex center position, identify hotspot and pocket residues
    (backward compatibility function)
    
    Args:
        infile: PDB file path
        hotspot_cutoff: Cutoff distance for hotspot residues (Angstrom)
        pocket_cutoff: Cutoff distance for pocket residues (Angstrom)
        verbose: Whether to output detailed information
        explicit_protein_chain_ids: Optional, explicitly specify protein chain ID list
        explicit_ligand_chain_ids: Optional, explicitly specify ligand chain ID list
        
    Returns:
        struct: Processed structure
        center_pos: Calculated center position
        raw_seq_data: Sequence data dictionary
    """
    return get_motif_center_pos(
        infile=infile, 
        receptor_type='protein', 
        ligand_type='nucleic',
        hotspot_cutoff=hotspot_cutoff,
        pocket_cutoff=pocket_cutoff,
        verbose=verbose,
        explicit_receptor_chain_ids=explicit_protein_chain_ids,
        explicit_ligand_chain_ids=explicit_ligand_chain_ids
    )


def get_hotspot_complex_struct_legacy(infile:str, hotspot_cutoff=8, verbose=False, 
                                     explicit_protein_chain_ids: list = None, 
                                     explicit_ligand_chain_ids: list = None):
    """
    Calculate protein-nucleic acid complex structure containing only hotspot residues and ligand
    (backward compatibility function)
    
    Args:
        infile: PDB file path
        hotspot_cutoff: Cutoff distance for hotspot residues (Angstrom)
        verbose: Whether to output detailed information
        explicit_protein_chain_ids: Optional, explicitly specify protein chain ID list
        explicit_ligand_chain_ids: Optional, explicitly specify ligand chain ID list
        
    Returns:
        struct: Processed structure containing only ligand and hotspot residues
    """
    return get_hotspot_complex_struct(
        infile=infile,
        receptor_type='protein',
        ligand_type='nucleic',
        hotspot_cutoff=hotspot_cutoff,
        verbose=verbose,
        explicit_receptor_chain_ids=explicit_protein_chain_ids,
        explicit_ligand_chain_ids=explicit_ligand_chain_ids
    )


def process_file(file_path:str, write_dir:str, hotspot_cutoff:int=8, pocket_cutoff:int=10):
    """
    Process protein-nucleic acid complex files and generate usable pickle data
    
    Args:
        file_path: File path to read
        write_dir: Directory to write results to
        hotspot_cutoff: Cutoff distance for hotspot residues
        pocket_cutoff: Cutoff distance for pocket residues
        
    Returns:
        Saves processed protein to pickle and returns metadata
    
    Raises:
        DataError: If known filtering rules are triggered
    """
    metadata = {}
    pdb_name = os.path.basename(file_path).replace('.pdb', '')
    metadata['pdb_name'] = pdb_name
    processed_path = os.path.join(write_dir, f'{pdb_name}.pkl')
    metadata['processed_path'] = os.path.abspath(processed_path)
    
    try:
        structure, center_pos, raw_seq_data = get_motif_center_pos(
            file_path, 
            hotspot_cutoff=hotspot_cutoff, 
            pocket_cutoff=pocket_cutoff
        )
        raw_seq_dict = {pdb_name: raw_seq_data}
    except Exception as e:
        print(f'Failed to parse {pdb_name} with error {e}')
        return None, None
    
    # 提取所有链
    struct_chains = {}
    ligand_chains = raw_seq_data["ligand_chains"]
    receptor_chains = raw_seq_data["receptor_chains"]
    
    # 添加所有配体链
    for lig_chain in ligand_chains:
        try:
            struct_chains[lig_chain.upper()] = structure[lig_chain]
        except:
            print(f"Warning: Could not find ligand chain {lig_chain} in processed structure")
    
    for rec_chain in receptor_chains:
        try:
            struct_chains[rec_chain.upper()] = structure[rec_chain]
        except:
            print(f"Warning: Could not find receptor chain {rec_chain} in processed structure")
    
    com_center = center_pos
    metadata['num_chains'] = len(struct_chains)
    metadata['num_ligand_chains'] = len(ligand_chains)
    metadata['num_receptor_chains'] = len(receptor_chains)
    
    struct_feats = []
    all_seqs = set()
    complex_length = 0
    for chain_id, chain in struct_chains.items():
        complex_length += len([i for i in chain.get_residues()])
    
    chain_masks = {}
    res_count = 0
    for chain_id_str, chain in struct_chains.items():
        chain_id = du.chain_str_to_int(chain_id_str)
        chain_prot = parsers.process_chain(chain, chain_id)
        chain_dict = dataclasses.asdict(chain_prot)
        chain_dict = du.parse_chain_feats(chain_dict, center_pos=com_center)
        all_seqs.add(tuple(chain_dict['aatype']))
        struct_feats.append(chain_dict)
        
        chain_mask = np.zeros(complex_length)
        chain_mask[res_count: res_count + len(chain_dict['aatype'])] = 1
        chain_masks[chain_id_str] = chain_mask
        res_count += len(chain_dict['aatype'])
    
    if len(all_seqs) == 1:
        metadata['quaternary_category'] = 'homomer'
    else:
        metadata['quaternary_category'] = 'heteromer'
    
    complex_feats = du.concat_np_features(struct_feats, False)
    complex_feats['center_pos'] = center_pos
    
    complex_aatype = complex_feats['aatype']
    metadata['seq_len'] = len(complex_aatype)
    modeled_idx = np.where(complex_aatype != 20)[0]
    if np.sum(complex_aatype != 20) == 0:
        raise errors.LengthError('No modeled residues')
    
    metadata['modeled_seq_len'] = len(modeled_idx)
    complex_feats['modeled_idx'] = modeled_idx
    
    complex_feats['ligand_mask'] = np.zeros(complex_length)
    res_count = 0
    for chain_id_str, chain in struct_chains.items():
        if chain_id_str in ligand_chains:
            chain_len = len([i for i in chain.get_residues()])
            complex_feats['ligand_mask'][res_count:res_count + chain_len] = 1
        res_count += len([i for i in chain.get_residues()])
    
    try:
        traj = md.load(file_path)
        pdb_ss = md.compute_dssp(traj, simplified=True)
        pdb_dg = md.compute_rg(traj)
    except Exception as e:
        raise errors.DataError(f'Mdtraj failed with error {e}')
    
    chain_dict['ss'] = pdb_ss[0]
    metadata['coil_percent'] = np.sum(pdb_ss == 'C') / metadata['modeled_seq_len']
    metadata['helix_percent'] = np.sum(pdb_ss == 'H') / metadata['modeled_seq_len']
    metadata['strand_percent'] = np.sum(pdb_ss == 'E') / metadata['modeled_seq_len']
    
    metadata['radius_gyration'] = pdb_dg[0]
    
    du.write_pkl(processed_path, complex_feats)
    
    return metadata, raw_seq_dict


def process_serially(all_paths, write_dir, hotspot_cutoff=8, pocket_cutoff=10):
    all_metadata = []
    all_raw_data = {}

    for i, file_path in enumerate(all_paths):
        try:
            start_time = time.time()
            metadata, raw_seq_data = process_file(
                file_path,
                write_dir,
                hotspot_cutoff=hotspot_cutoff,
                pocket_cutoff=pocket_cutoff
            )
            elapsed_time = time.time() - start_time
            print(f'Finished {file_path} in {elapsed_time:2.2f}s')
            all_metadata.append(metadata)
            all_raw_data.update(raw_seq_data)
        except errors.DataError as e:
            print(f'Failed {file_path}: {e}')

    return all_metadata, all_raw_data


def process_fn(
        file_path,
        verbose=None,
        write_dir=None,
        hotspot_cutoff=8,
        pocket_cutoff=10        
        ):
    try:
        start_time = time.time()
        metadata, raw_seq_data = process_file(
            file_path,
            write_dir,
            hotspot_cutoff=hotspot_cutoff,
            pocket_cutoff=pocket_cutoff
        )
        elapsed_time = time.time() - start_time
        if verbose:
            print(f'Finished {file_path} in {elapsed_time:2.2f}s')
        return metadata, raw_seq_data
    except errors.DataError as e:
        if verbose:
            print(f'Failed {file_path}: {e}')
        return None, None


def main(args):
    pdb_dir = args.pdb_dir
    for file_name in os.listdir(pdb_dir):
        if file_name.endswith('_processed.pdb'):
            file_path = os.path.join(pdb_dir, file_name)
            try:
                os.remove(file_path)
                print(f'Removed file: {file_path}')
            except OSError as e:
                print(f'Error while deleting file {file_path}: {e}')
                
    all_file_paths = [
        os.path.join(pdb_dir, x)
        for x in os.listdir(args.pdb_dir) if '.pdb' in x]
    total_num_paths = len(all_file_paths)
    write_dir = args.write_dir
    if not os.path.exists(write_dir):
        os.makedirs(write_dir)
    metadata_file_name = 'metadata.csv'
    metadata_path = os.path.join(write_dir, metadata_file_name)
    print(f'Files will be written to {write_dir}')

    # Process each pdb file
    if args.num_processes == 1:
        all_metadata, all_raw_data = process_serially(
            all_file_paths,
            write_dir,
            hotspot_cutoff=args.hotspot_cutoff,
            pocket_cutoff=args.pocket_cutoff
        )
    else:
        _process_fn = fn.partial(
            process_fn,
            verbose=args.verbose,
            write_dir=write_dir,
            hotspot_cutoff=args.hotspot_cutoff,
            pocket_cutoff=args.pocket_cutoff
        )
        with mp.Pool(processes=args.num_processes) as pool:
            results = pool.map(_process_fn, all_file_paths)
        all_metadata, all_raw_data = [], {}
        for metadata, raw_seq_data in results:
            if metadata is not None:
                all_metadata.append(metadata)
            if raw_seq_data is not None:
                all_raw_data.update(raw_seq_data)

    metadata_df = pd.DataFrame(all_metadata)
    metadata_df.to_csv(metadata_path, index=False)
    succeeded = len(all_metadata)
    print(f'Finished processing {succeeded}/{total_num_paths} files.')


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)

 