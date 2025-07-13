"""
MSA masking using AlphaFold3's native alignment processing methods.

This module uses AF3's extract_msa_features and related functions to properly
handle sequence alignments before applying pocket masks.
"""

import numpy as np
from typing import Sequence, Tuple
import string

# Import AF3's MSA processing functions
try:
    from alphafold3.data.msa_features import extract_msa_features
    from alphafold3.constants import mmcif_names
    from alphafold3.data import parsers
    AF3_AVAILABLE = True
except ImportError:
    AF3_AVAILABLE = False
    print("Warning: AlphaFold3 modules not available, falling back to simple masking")

# Character mappings from AF3 for reverse conversion
_ID_TO_PROTEIN = {
    0: 'A', 1: 'R', 2: 'N', 3: 'D', 4: 'C', 5: 'Q', 6: 'E', 7: 'G',
    8: 'H', 9: 'I', 10: 'L', 11: 'K', 12: 'M', 13: 'F', 14: 'P', 15: 'S',
    16: 'T', 17: 'W', 18: 'Y', 19: 'V', 20: 'X', 21: '-'
}

_ID_TO_RNA = {
    21: '-', 22: 'A', 23: 'G', 24: 'C', 25: 'U', 30: 'N'
}

_ID_TO_DNA = {
    21: '-', 26: 'A', 27: 'G', 28: 'C', 29: 'T', 30: 'N'
}

def get_chain_poly_type(chain_type: str) -> str:
    """Convert our chain type to AF3's polymer type constants."""
    if chain_type.lower() == 'protein':
        return mmcif_names.PROTEIN_CHAIN
    elif chain_type.lower() == 'rna':
        return mmcif_names.RNA_CHAIN
    elif chain_type.lower() == 'dna':
        return mmcif_names.DNA_CHAIN
    else:
        raise ValueError(f"Unknown chain type: {chain_type}")

def get_id_to_char_map(chain_type: str) -> dict:
    """Get the ID to character mapping for the given chain type."""
    if chain_type.lower() == 'protein':
        return _ID_TO_PROTEIN
    elif chain_type.lower() == 'rna':
        return _ID_TO_RNA
    elif chain_type.lower() == 'dna':
        return _ID_TO_DNA
    else:
        raise ValueError(f"Unknown chain type: {chain_type}")

def msa_arrays_to_sequences(msa_array: np.ndarray, deletions_array: np.ndarray, 
                           chain_type: str) -> Sequence[str]:
    """
    Convert AF3's MSA arrays back to sequence strings.
    
    Args:
        msa_array: MSA array from AF3's extract_msa_features
        deletions_array: Deletions array from AF3's extract_msa_features
        chain_type: Chain type ('protein', 'rna', 'dna')
        
    Returns:
        List of sequence strings
    """
    id_to_char = get_id_to_char_map(chain_type)
    sequences = []
    
    for row_idx in range(msa_array.shape[0]):
        sequence_parts = []
        for col_idx in range(msa_array.shape[1]):
            # Add deletions (insertions from original sequence)
            num_deletions = deletions_array[row_idx, col_idx]
            # Use lowercase residues for insertions (AF3 standard)
            if chain_type.lower() == 'protein':
                insertion_chars = ['a'] * num_deletions  # Use 'a' for protein insertions
            else:
                insertion_chars = ['a'] * num_deletions  # Use 'a' for nucleic acid insertions
            sequence_parts.extend(insertion_chars)
            
            # Add the actual residue
            residue_id = msa_array[row_idx, col_idx]
            residue_char = id_to_char.get(residue_id, 'X')
            sequence_parts.append(residue_char)
        
        sequences.append(''.join(sequence_parts))
    
    return sequences

def apply_mask_to_msa_array(msa_array: np.ndarray, deletions_array: np.ndarray,
                           pocket_mask: Sequence[int], chain_type: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply pocket mask to AF3's MSA arrays.
    
    Args:
        msa_array: MSA array from AF3's extract_msa_features
        deletions_array: Deletions array from AF3's extract_msa_features  
        pocket_mask: Pocket mask (1 for pocket positions, 0 for others)
        chain_type: Chain type ('protein', 'rna', 'dna')
        
    Returns:
        Tuple of (masked_msa_array, masked_deletions_array)
    """
    # Handle length mismatch
    if msa_array.shape[1] != len(pocket_mask):
        if msa_array.shape[1] > len(pocket_mask):
            # Truncate to mask length
            msa_array = msa_array[:, :len(pocket_mask)]
            deletions_array = deletions_array[:, :len(pocket_mask)]
        else:
            # Pad mask to MSA length  
            pocket_mask = list(pocket_mask) + [0] * (msa_array.shape[1] - len(pocket_mask))
    
    # Convert mask to numpy array
    mask_array = np.array(pocket_mask, dtype=bool)
    
    # Create masked arrays
    masked_msa_array = msa_array.copy()
    masked_deletions_array = deletions_array.copy()
    
    # Get gap ID for this chain type
    if chain_type.lower() == 'protein':
        gap_id = 21  # '-' for protein
    else:  # RNA/DNA
        gap_id = 21  # '-' for nucleic acids
    
    # Apply mask (keep first sequence unchanged)
    for row_idx in range(1, msa_array.shape[0]):  # Skip first sequence (query)
        for col_idx in range(msa_array.shape[1]):
            if not mask_array[col_idx]:  # Not a pocket position
                masked_msa_array[row_idx, col_idx] = gap_id
                masked_deletions_array[row_idx, col_idx] = 0  # Remove insertions for masked positions
    
    return masked_msa_array, masked_deletions_array

def apply_pocket_mask_to_msa_af3(msa_string: str, pocket_mask: Sequence[int], 
                                chain_type: str = 'protein', verbose: bool = False) -> str:
    """
    Apply pocket mask to MSA using AlphaFold3's native alignment processing.
    
    Args:
        msa_string: MSA in FASTA format
        pocket_mask: Pocket mask array (1 for pocket positions, 0 for others)
        chain_type: Chain type ('protein', 'rna', 'dna')
        verbose: Whether to print debug information
        
    Returns:
        Masked MSA string
    """
    if not AF3_AVAILABLE:
        # Fallback to simple method if AF3 not available
        from .masking_af3_msa import apply_pocket_mask_to_msa
        return apply_pocket_mask_to_msa(msa_string, pocket_mask, verbose)
    
    if not msa_string or not pocket_mask:
        return msa_string
    
    try:
        # Step 1: Parse FASTA to get sequences and descriptions
        sequences, descriptions = parsers.parse_fasta(msa_string)
        
        if not sequences:
            return msa_string
        
        if verbose:
            print(f"Parsed {len(sequences)} sequences from MSA")
            print(f"First sequence: {sequences[0][:50]}...")
        
        # Step 2: Use AF3's extract_msa_features to get standardized MSA
        poly_type = get_chain_poly_type(chain_type)
        msa_array, deletions_array = extract_msa_features(sequences, poly_type)
        
        if verbose:
            print(f"AF3 MSA array shape: {msa_array.shape}")
            print(f"AF3 deletions array shape: {deletions_array.shape}")
            print(f"Pocket mask length: {len(pocket_mask)}")
        
        # Step 3: Apply pocket mask to the standardized MSA
        masked_msa_array, masked_deletions_array = apply_mask_to_msa_array(
            msa_array, deletions_array, pocket_mask, chain_type
        )
        
        # Step 4: Convert back to sequence strings
        masked_sequences = msa_arrays_to_sequences(masked_msa_array, masked_deletions_array, chain_type)
        
        if verbose:
            print(f"Generated {len(masked_sequences)} masked sequences")
        
        # Step 5: Remove duplicates (using AF3's approach - ignore lowercase)
        deletion_table = str.maketrans('', '', string.ascii_lowercase)
        unique_sequences = []
        unique_descriptions = []
        seen_sequences = set()
        
        for i, (seq, desc) in enumerate(zip(masked_sequences, descriptions)):
            if i == 0:
                # Always keep the first sequence (query)
                unique_sequences.append(sequences[0])  # Use original query sequence
                unique_descriptions.append(desc)
                seq_no_insertions = sequences[0].translate(deletion_table)
                seen_sequences.add(seq_no_insertions)
                if verbose:
                    print(f"Kept query sequence unchanged: {sequences[0][:50]}...")
            else:
                # Remove lowercase and check for duplicates
                seq_no_insertions = seq.translate(deletion_table)
                if seq_no_insertions not in seen_sequences:
                    seen_sequences.add(seq_no_insertions)
                    unique_sequences.append(seq)
                    unique_descriptions.append(desc)
        
        # Step 6: Reconstruct FASTA format
        result_lines = []
        for seq, desc in zip(unique_sequences, unique_descriptions):
            result_lines.append(f">{desc}")
            result_lines.append(seq)
        
        result = '\n'.join(result_lines)
        
        if verbose:
            original_count = len(sequences)
            final_count = len(unique_sequences)
            print(f"AF3-style MSA masking: {original_count} -> {final_count} sequences")
        
        return result
        
    except Exception as e:
        if verbose:
            print(f"Error in AF3 MSA processing: {e}")
            print("Falling back to simple masking method")
        
        # Fallback to simple method if AF3 processing fails
        from .masking_af3_msa import apply_pocket_mask_to_msa
        return apply_pocket_mask_to_msa(msa_string, pocket_mask, verbose)

# For backward compatibility
def apply_pocket_mask_to_msa(msa_string: str, pocket_mask: Sequence[int], verbose: bool = False) -> str:
    """Backward compatibility wrapper that uses AF3 method by default."""
    return apply_pocket_mask_to_msa_af3(msa_string, pocket_mask, 'protein', verbose) 