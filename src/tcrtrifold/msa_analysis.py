"""
MSA analysis utilities for calculating Neff and DCA interface signal.

This module provides functions to:
1. Calculate the number of effective sequences (Neff) from MSAs
2. Perform direct coupling analysis (DCA) to measure interface signal
3. Process both paired and unpaired MSAs for p:MHC and triad data
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import polars as pl
from .utils import hash_sequence


def parse_msa(msa_string: str) -> Dict[str, str]:
    """
    Parse FASTA format MSA string into a dictionary.
    
    Args:
        msa_string (str): MSA in FASTA format
        
    Returns:
        Dict[str, str]: Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    for line in msa_string.strip().split('\n'):
        if line.startswith('>'):
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line.strip())
    
    if current_id is not None:
        sequences[current_id] = ''.join(current_seq)
    
    return sequences


def calculate_neff(msa_sequences: Dict[str, str], threshold: float = 0.8) -> float:
    """
    Calculate the number of effective sequences (Neff) from an MSA.
    
    Uses the approach described in the AlphaFold3 paper: counting non-redundant
    sequences with a threshold of 80% sequence identity.
    
    Args:
        msa_sequences (Dict[str, str]): Dictionary mapping sequence IDs to sequences
        threshold (float): Sequence identity threshold (default 0.8)
        
    Returns:
        float: Number of effective sequences (Neff)
    """
    if not msa_sequences:
        return 0.0
    
    sequences = list(msa_sequences.values())
    n_sequences = len(sequences)
    
    if n_sequences == 1:
        return 1.0
    
    # Calculate sequence weights using the Neff scheme
    weights = np.ones(n_sequences)
    
    for i in range(n_sequences):
        seq_i = sequences[i]
        n_similar = 0
        
        for j in range(n_sequences):
            seq_j = sequences[j]
            
            # Calculate sequence identity
            # Handle sequences of different lengths by using alignment length
            min_len = min(len(seq_i), len(seq_j))
            if min_len == 0:
                identity = 0.0
            else:
                # Count matches at each position, ignoring gaps
                matches = 0
                valid_positions = 0
                
                for k in range(min_len):
                    if seq_i[k] != '-' and seq_j[k] != '-':
                        valid_positions += 1
                        if seq_i[k] == seq_j[k]:
                            matches += 1
                
                identity = matches / valid_positions if valid_positions > 0 else 0.0
            
            if identity >= threshold:
                n_similar += 1
        
        weights[i] = 1.0 / n_similar if n_similar > 0 else 1.0
    
    return np.sum(weights)


def calculate_sequence_identity_matrix(msa_sequences: Dict[str, str]) -> np.ndarray:
    """
    Calculate pairwise sequence identity matrix for MSA sequences.
    
    Args:
        msa_sequences (Dict[str, str]): Dictionary mapping sequence IDs to sequences
        
    Returns:
        np.ndarray: Symmetric matrix of pairwise sequence identities
    """
    sequences = list(msa_sequences.values())
    n_sequences = len(sequences)
    identity_matrix = np.zeros((n_sequences, n_sequences))
    
    for i in range(n_sequences):
        for j in range(i, n_sequences):
            seq_i, seq_j = sequences[i], sequences[j]
            
            # Calculate sequence identity
            min_len = min(len(seq_i), len(seq_j))
            if min_len == 0:
                identity = 0.0
            else:
                matches = 0
                valid_positions = 0
                
                for k in range(min_len):
                    if seq_i[k] != '-' and seq_j[k] != '-':
                        valid_positions += 1
                        if seq_i[k] == seq_j[k]:
                            matches += 1
                
                identity = matches / valid_positions if valid_positions > 0 else 0.0
            
            identity_matrix[i, j] = identity
            identity_matrix[j, i] = identity
    
    return identity_matrix


def simple_dca_score(msa_sequences: Dict[str, str], query_positions: List[int]) -> float:
    """
    Calculate a simplified DCA-like score for measuring interface signal.
    
    This is a simplified version that measures covariation at interface positions.
    A full DCA implementation would require more sophisticated statistical modeling.
    
    Args:
        msa_sequences (Dict[str, str]): Dictionary mapping sequence IDs to sequences
        query_positions (List[int]): List of positions to analyze for interface signal
        
    Returns:
        float: Simplified DCA score representing interface signal strength
    """
    if not msa_sequences or not query_positions:
        return 0.0
    
    sequences = list(msa_sequences.values())
    n_sequences = len(sequences)
    
    if n_sequences < 2:
        return 0.0
    
    # Get the query sequence (first sequence)
    query_seq = sequences[0]
    query_len = len(query_seq)
    
    # Filter positions that are within sequence bounds
    valid_positions = [pos for pos in query_positions if 0 <= pos < query_len]
    
    if len(valid_positions) < 2:
        return 0.0
    
    # Calculate conservation at interface positions
    conservation_scores = []
    
    for pos in valid_positions:
        # Count amino acid frequencies at this position
        aa_counts = {}
        valid_seqs = 0
        
        for seq in sequences:
            if pos < len(seq) and seq[pos] != '-':
                aa = seq[pos]
                aa_counts[aa] = aa_counts.get(aa, 0) + 1
                valid_seqs += 1
        
        if valid_seqs > 0:
            # Calculate entropy as measure of conservation
            frequencies = np.array(list(aa_counts.values())) / valid_seqs
            entropy = -np.sum(frequencies * np.log2(frequencies + 1e-10))
            conservation = 1.0 - entropy / np.log2(20)  # Normalize by max entropy
            conservation_scores.append(max(0.0, conservation))
    
    return np.mean(conservation_scores) if conservation_scores else 0.0


def load_msa_from_file(msa_file_path: Union[str, Path]) -> Dict[str, Union[str, Dict[str, str]]]:
    """
    Load MSA data from JSON file.
    
    Args:
        msa_file_path (Union[str, Path]): Path to MSA JSON file
        
    Returns:
        Dict: MSA data with parsed sequences
    """
    with open(msa_file_path, 'r') as f:
        msa_data = json.load(f)
    
    result = {
        'sequence': msa_data['sequence'],
        'unpairedMsa': parse_msa(msa_data['unpairedMsa']),
        'pairedMsa': parse_msa(msa_data['pairedMsa']),
    }
    
    return result


def get_msa_file_path(sequence: str, msa_dir: Union[str, Path], seq_type: str = 'mhc') -> Optional[Path]:
    """
    Get MSA file path for a given sequence using hash-based lookup.
    
    Args:
        sequence (str): Protein sequence
        msa_dir (Union[str, Path]): Directory containing MSA files
        seq_type (str): Type of sequence ('mhc' or 'tcr')
        
    Returns:
        Optional[Path]: Path to MSA file if it exists, None otherwise
    """
    seq_hash = hash_sequence(sequence, hash_type="sha256")
    msa_file = Path(msa_dir) / seq_type / f"{seq_hash}.json"
    
    return msa_file if msa_file.exists() else None


def analyze_msa_metrics(
    sequence: str,
    msa_dir: Union[str, Path],
    seq_type: str = 'mhc',
    interface_positions: Optional[List[int]] = None
) -> Dict[str, float]:
    """
    Analyze MSA metrics for a given sequence.
    
    Args:
        sequence (str): Protein sequence
        msa_dir (Union[str, Path]): Directory containing MSA files
        seq_type (str): Type of sequence ('mhc' or 'tcr')
        interface_positions (Optional[List[int]]): Positions to analyze for interface signal
        
    Returns:
        Dict[str, float]: Dictionary containing MSA metrics
    """
    msa_file_path = get_msa_file_path(sequence, msa_dir, seq_type)
    
    if msa_file_path is None:
        return {
            'unpaired_neff': 0.0,
            'paired_neff': 0.0,
            'unpaired_dca_score': 0.0,
            'paired_dca_score': 0.0,
        }
    
    try:
        msa_data = load_msa_from_file(msa_file_path)
        
        # Calculate Neff for both MSAs
        unpaired_neff = calculate_neff(msa_data['unpairedMsa'])
        paired_neff = calculate_neff(msa_data['pairedMsa'])
        
        # Calculate DCA scores if interface positions are provided
        unpaired_dca = 0.0
        paired_dca = 0.0
        
        if interface_positions:
            unpaired_dca = simple_dca_score(msa_data['unpairedMsa'], interface_positions)
            paired_dca = simple_dca_score(msa_data['pairedMsa'], interface_positions)
        
        return {
            'unpaired_neff': unpaired_neff,
            'paired_neff': paired_neff,
            'unpaired_dca_score': unpaired_dca,
            'paired_dca_score': paired_dca,
        }
    
    except Exception as e:
        print(f"Error processing MSA file {msa_file_path}: {e}")
        return {
            'unpaired_neff': 0.0,
            'paired_neff': 0.0,
            'unpaired_dca_score': 0.0,
            'paired_dca_score': 0.0,
        }


def estimate_interface_positions(sequence_length: int, interface_ratio: float = 0.3) -> List[int]:
    """
    Estimate interface positions for a sequence.
    
    This is a simplified approach - in practice, interface positions would be
    determined from structural data or domain knowledge.
    
    Args:
        sequence_length (int): Length of the sequence
        interface_ratio (float): Fraction of sequence assumed to be interface
        
    Returns:
        List[int]: Estimated interface positions
    """
    n_interface = max(1, int(sequence_length * interface_ratio))
    
    # Assume interface positions are in the middle region
    start_pos = max(0, sequence_length // 2 - n_interface // 2)
    end_pos = min(sequence_length, start_pos + n_interface)
    
    return list(range(start_pos, end_pos))