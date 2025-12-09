import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import os


def read_profile(file_path):
    """
    Read profile CSV file (space-separated) and return a (6, L) numpy array of probabilities.
    """
    with open(file_path) as f:
        lines = f.readlines()
    if len(lines) < 7:
        raise RuntimeError(f"Profile file {file_path} has insufficient lines: {len(lines)}")
    probs = []
    for line in lines[1:7]:
        parts = line.strip().split()
        # parts[0] is type name, the rest are numeric values
        try:
            row = [float(x) for x in parts[1:]]
        except ValueError as e:
            raise RuntimeError(f"Failed to parse floats in {file_path}: {e}")
        probs.append(row)
    return np.array(probs)

def read_alignment(file_path):
    """
    Read a FASTA-formatted alignment and return a dict of id -> aligned sequence (with gaps).
    """
    seqs = {}
    current_id = None
    current_seq = []
    with open(file_path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                # Save previous
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_id is not None:
        seqs[current_id] = ''.join(current_seq)
    return seqs



def map_alignment_to_original(seqs):
    """
    Build a mapping from alignment columns to original sequence indices for each sequence.
    Returns dict: {id: list of length alignment_length, with original 0-based index or None for gaps}.
    """
    mapping = {}
    aln_len = len(next(iter(seqs.values())))
    for sid, seq in seqs.items():
        orig_idx = 0
        map_list = []
        for char in seq:
            if char == '-':
                map_list.append(None)
            else:
                map_list.append(orig_idx)
                orig_idx += 1
        if len(map_list) != aln_len:
            raise RuntimeError(f"Alignment length mismatch for {sid}: expected {aln_len}, got {len(map_list)}")
        mapping[sid] = map_list
    return mapping

def get_non_gap_positions(seqs):
    """
    Return list of alignment indices where no sequence has a gap at that position.
    """
    aln_len = len(next(iter(seqs.values())))
    non_gap = []
    for i in range(aln_len):
        if all(seq[i] != '-' for seq in seqs.values()):
            non_gap.append(i)
    return non_gap