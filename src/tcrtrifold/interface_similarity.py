"""
Module for calculating interface similarity between complex sequences and their MSA alignments.
"""

import json
import numpy as np
import polars as pl
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from difflib import SequenceMatcher

from tcrtrifold.utils import hash_sequence
from tcrtrifold.mhc import _hla_slice_to_topological, _h2_slice_to_topological
from tcrtrifold.tcr import tcr_by_imgt_region


def parse_fasta_sequences(fasta_string: str) -> List[Tuple[str, str]]:
    """
    Parse FASTA format string into list of (header, sequence) tuples.

    Args:
        fasta_string: FASTA format string

    Returns:
        List of (header, sequence) tuples
    """
    sequences = []
    current_seq = ""
    current_header = ""

    for line in fasta_string.split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                sequences.append((current_header, current_seq))
            current_header = line
            current_seq = ""
        else:
            current_seq += line

    if current_seq:
        sequences.append((current_header, current_seq))

    return sequences


def extract_mhc_interface_region(sequence: str, mhc_name: str = None) -> str:
    """
    Extract the interface region (extracellular domain) of an MHC sequence.

    Args:
        sequence: MHC sequence
        mhc_name: MHC allele name (optional, for better domain extraction)

    Returns:
        Interface region sequence
    """
    if mhc_name:
        try:
            return _hla_slice_to_topological(mhc_name, sequence)
        except ValueError:
            # Try H2 if HLA fails
            try:
                return _h2_slice_to_topological(mhc_name, sequence, chain=None)
            except:
                pass

    # Default fallback: assume it's already the interface region
    return sequence


def extract_tcr_interface_region(
    sequence: str, chain: str, species: str = "human"
) -> str:
    """
    Extract the interface region (CDR regions) of a TCR sequence.

    Args:
        sequence: TCR sequence
        chain: TCR chain type ("alpha" or "beta")
        species: Species (default "human")

    Returns:
        Concatenated CDR regions (CDR1 + CDR2 + CDR3)
    """
    try:
        resindices = np.arange(len(sequence))
        regions = tcr_by_imgt_region(sequence, resindices, chain, species)

        # Concatenate CDR regions
        cdr_regions = []
        for region in ["cdr_1", "cdr_2", "cdr_3"]:
            if region in regions:
                cdr_seq = "".join(np.array(list(sequence))[regions[region]])
                cdr_regions.append(cdr_seq)

        return "".join(cdr_regions)
    except Exception:
        # Fallback: return full sequence if CDR extraction fails
        return sequence


def load_msa_data(msa_dir: Path, sequence: str) -> Dict:
    """
    Load MSA data for a given sequence.

    Args:
        msa_dir: Directory containing MSA files
        sequence: Sequence to look up

    Returns:
        MSA data dictionary
    """
    seq_hash = hash_sequence(sequence, "sha256")
    msa_file = msa_dir / f"{seq_hash}.json"

    if not msa_file.exists():
        raise FileNotFoundError(f"MSA file not found: {msa_file}")

    with open(msa_file, "r") as f:
        return json.load(f)


def calculate_sequence_similarity(seq1: str, seq2: str) -> float:
    """
    Calculate sequence similarity using SequenceMatcher.

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        Similarity ratio (0.0 to 1.0)
    """
    if not seq1 or not seq2:
        return 0.0

    # Remove gaps
    seq1_clean = seq1.replace("-", "")
    seq2_clean = seq2.replace("-", "")

    if not seq1_clean or not seq2_clean:
        return 0.0

    matcher = SequenceMatcher(None, seq1_clean, seq2_clean)
    return matcher.ratio()


def extract_interface_similarities_from_msa(
    msa_data: Dict, query_interface: str, max_sequences: int = 100
) -> Dict[str, float]:
    """
    Extract interface similarities from MSA data.

    Args:
        msa_data: MSA data dictionary
        query_interface: Query interface sequence
        max_sequences: Maximum number of sequences to compare

    Returns:
        Dictionary with similarity metrics
    """
    paired_msa = msa_data.get("pairedMsa", "")
    if not paired_msa:
        return {"mean_similarity": 0.0, "max_similarity": 0.0, "num_sequences": 0}

    # Parse sequences from pairedMsa
    sequences = parse_fasta_sequences(paired_msa)

    # Skip query sequence (first one)
    msa_sequences = sequences[1 : max_sequences + 1]

    similarities = []
    for header, seq in msa_sequences:
        # Extract interface region from MSA sequence
        # For now, assume the MSA sequences have the same structure as query
        msa_interface = seq.replace("-", "")  # Remove gaps

        if msa_interface:
            similarity = calculate_sequence_similarity(query_interface, msa_interface)
            similarities.append(similarity)

    if not similarities:
        return {"mean_similarity": 0.0, "max_similarity": 0.0, "num_sequences": 0}

    return {
        "mean_similarity": np.mean(similarities),
        "max_similarity": np.max(similarities),
        "num_sequences": len(similarities),
    }


def calculate_pmhc_interface_similarity(
    row: Dict, msa_dir: Path, mhc_subdir: str = "mhc"
) -> Dict:
    """
    Calculate interface similarity for a p:MHC complex.

    Args:
        row: Row from parquet file with sequence data
        msa_dir: Base MSA directory
        mhc_subdir: MHC subdirectory name

    Returns:
        Dictionary with similarity metrics
    """
    result = {}

    # Extract MHC interface similarities
    for mhc_idx in [1, 2]:
        mhc_seq_col = f"mhc_{mhc_idx}_seq"
        mhc_name_col = f"mhc_{mhc_idx}_name"

        if mhc_seq_col in row and row[mhc_seq_col] is not None:
            mhc_seq = row[mhc_seq_col]
            mhc_name = row.get(mhc_name_col, None)

            try:
                # Load MSA data
                msa_data = load_msa_data(msa_dir / mhc_subdir, mhc_seq)

                # Extract interface region
                query_interface = extract_mhc_interface_region(mhc_seq, mhc_name)

                # Calculate similarities
                similarities = extract_interface_similarities_from_msa(
                    msa_data, query_interface
                )

                # Add to result with prefix
                for key, value in similarities.items():
                    result[f"mhc_{mhc_idx}_interface_{key}"] = value

            except Exception as e:
                # Set default values if MSA loading fails
                for key in ["mean_similarity", "max_similarity", "num_sequences"]:
                    result[f"mhc_{mhc_idx}_interface_{key}"] = (
                        0.0 if key != "num_sequences" else 0
                    )

    return result


def calculate_triad_interface_similarity(
    row: Dict, msa_dir: Path, mhc_subdir: str = "mhc", tcr_subdir: str = "tcr"
) -> Dict:
    """
    Calculate interface similarity for a TCR:p:MHC triad.

    Args:
        row: Row from parquet file with sequence data
        msa_dir: Base MSA directory
        mhc_subdir: MHC subdirectory name
        tcr_subdir: TCR subdirectory name

    Returns:
        Dictionary with similarity metrics
    """
    result = {}

    # First calculate p:MHC similarities
    pmhc_similarities = calculate_pmhc_interface_similarity(row, msa_dir, mhc_subdir)
    result.update(pmhc_similarities)

    # Extract TCR interface similarities
    tcr_chains = {"1": "alpha", "2": "beta"}  # Default mapping

    for tcr_idx in [1, 2]:
        tcr_seq_col = f"tcr_{tcr_idx}_seq"
        tcr_chain_col = f"tcr_{tcr_idx}_chain"
        tcr_species_col = f"tcr_{tcr_idx}_species"

        if tcr_seq_col in row and row[tcr_seq_col] is not None:
            tcr_seq = row[tcr_seq_col]
            tcr_chain = row.get(tcr_chain_col, tcr_chains[str(tcr_idx)])
            tcr_species = row.get(tcr_species_col, "human")

            try:
                # Load MSA data
                msa_data = load_msa_data(msa_dir / tcr_subdir, tcr_seq)

                # Extract interface region
                query_interface = extract_tcr_interface_region(
                    tcr_seq, tcr_chain, tcr_species
                )

                # Calculate similarities
                similarities = extract_interface_similarities_from_msa(
                    msa_data, query_interface
                )

                # Add to result with prefix
                for key, value in similarities.items():
                    result[f"tcr_{tcr_idx}_interface_{key}"] = value

            except Exception as e:
                # Set default values if MSA loading fails
                for key in ["mean_similarity", "max_similarity", "num_sequences"]:
                    result[f"tcr_{tcr_idx}_interface_{key}"] = (
                        0.0 if key != "num_sequences" else 0
                    )

    return result


def add_interface_similarity_features(
    df: pl.DataFrame,
    msa_dir: Path,
    complex_type: str = "pmhc",
    mhc_subdir: str = "mhc",
    tcr_subdir: str = "tcr",
) -> pl.DataFrame:
    """
    Add interface similarity features to a DataFrame.

    Args:
        df: Input DataFrame
        msa_dir: Base MSA directory
        complex_type: Type of complex ("pmhc" or "triad")
        mhc_subdir: MHC subdirectory name
        tcr_subdir: TCR subdirectory name

    Returns:
        DataFrame with added interface similarity features
    """
    # Convert to list of dictionaries for processing
    rows = df.to_dicts()

    # Process each row
    for row in rows:
        if complex_type == "pmhc":
            similarities = calculate_pmhc_interface_similarity(row, msa_dir, mhc_subdir)
        elif complex_type == "triad":
            similarities = calculate_triad_interface_similarity(
                row, msa_dir, mhc_subdir, tcr_subdir
            )
        else:
            raise ValueError(f"Unknown complex type: {complex_type}")

        # Add similarities to row
        row.update(similarities)

    # Convert back to DataFrame
    return pl.DataFrame(rows)
