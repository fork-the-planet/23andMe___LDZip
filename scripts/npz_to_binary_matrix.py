#!/usr/bin/env python3
"""
Convert sparse COO format .npz LD matrix to dense binary float32 matrix (PLINK-style).

IMPORTANT: This script is designed specifically for Alkes Price lab NPZ file format
(sparse COO with 'row', 'col', 'data', 'shape' keys). It may not work with other
NPZ formats without modification.

Usage:
    python npz_to_binary_matrix.py input.npz variants.gz output.bin

Arguments:
    input.npz    : Input NPZ file containing sparse LD matrix (Alkes Price format)
    variants.gz  : Variant list file (rsid, chrom, pos, ref, alt)
    output.bin   : Output binary matrix file
"""

import sys
import os
import gzip
import time
import numpy as np
from scipy.sparse import coo_matrix

def create_vars_file(variants_path, output_path, num_variants, keep_indices=None):
    """
    Create a .vars file from the variant list (PLINK format).
    Appends :REF:ALT to duplicate RSIDs (multi-allelic variants) to make them unique.

    Args:
        variants_path: Path to .gz file with variant info
        output_path: Path to output .vars file
        num_variants: Expected number of variants
        keep_indices: Optional list of indices to keep (after duplicate removal)
    """
    print(f"Creating .vars file from {os.path.basename(variants_path)}...", flush=True)

    with gzip.open(variants_path, 'rt') as f:
        # Skip header
        f.readline()

        variants = []
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                rsid, chrom, pos, a1, a2 = parts[:5]
                variants.append((rsid, chrom, pos, a1, a2))

    if len(variants) != num_variants:
        print(f"  WARNING: variants in .gz ({len(variants)}) != matrix dim ({num_variants})", flush=True)

    # Filter to kept indices if provided
    if keep_indices is not None:
        variants = [variants[i] for i in keep_indices]

    # Detect duplicate RSIDs (multi-allelic: same rsID, different alleles) and make them unique by appending :REF:ALT
    rsid_counts = {}
    for rsid, chrom, pos, a1, a2 in variants:
        rsid_counts[rsid] = rsid_counts.get(rsid, 0) + 1

    duplicates = {rsid for rsid, count in rsid_counts.items() if count > 1}

    if duplicates:
        print(f"  Found {len(duplicates)} multi-allelic RSIDs, appending :REF:ALT to make unique", flush=True)

    # Write .vars file with header: #CHROM POS ID REF ALT
    with open(output_path, 'w') as f:
        f.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for rsid, chrom, pos, a1, a2 in variants:
            # If this RSID appears multiple times (multi-allelic), append :REF:ALT to make it unique
            if rsid in duplicates:
                unique_id = f"{rsid}:{a1}:{a2}"
            else:
                unique_id = rsid
            f.write(f"{chrom}\t{pos}\t{unique_id}\t{a1}\t{a2}\n")

    print(f"  Wrote {len(variants)} variants", flush=True)

def convert_npz_to_binary(npz_path, output_path, variants_path=None):
    """
    Load sparse .npz LD matrix and save as dense binary float32.
    Removes exact duplicate variants (same chr:pos:ref:alt) and corresponding matrix rows/cols.

    Args:
        npz_path: Path to input .npz file (COO sparse format)
        output_path: Path to output binary file
        variants_path: Path to .gz file with variant info
    """
    start_time = time.time()

    # Load sparse matrix data
    print(f"[1/5] Loading {npz_path}...", flush=True)
    data = np.load(npz_path)

    row = data['row']
    col = data['col']
    values = data['data']
    shape = tuple(data['shape'])

    print(f"  {shape[0]} x {shape[1]} matrix, {len(values):,} non-zero entries ({time.time()-start_time:.1f}s)", flush=True)

    # Load variants to detect exact duplicates and boundary variants
    keep_indices = None
    if variants_path:
        print(f"[2/5] Detecting exact duplicate variants and boundary variants...", flush=True)
        step_time = time.time()

        with gzip.open(variants_path, 'rt') as f:
            f.readline()  # Skip header
            variants = []
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    rsid, chrom, pos, a1, a2 = parts[:5]
                    variants.append((chrom, pos, a1, a2, rsid))

        # Check if first variant is at a boundary position (ends with 00001)
        skip_first = False
        if len(variants) > 0:
            first_pos = variants[0][1]  # position is second element
            if first_pos.endswith('00001'):
                print(f"  First variant at boundary position {first_pos}, will skip it", flush=True)
                skip_first = True

        # Find exact duplicates (same chr:pos:ref:alt)
        seen = {}
        keep_indices = []
        dup_indices = []

        for idx, var in enumerate(variants):
            # Skip first variant if it's at a boundary position
            if skip_first and idx == 0:
                dup_indices.append(idx)
                continue

            var_key = var[:4]  # (chrom, pos, ref, alt)
            if var_key in seen:
                dup_indices.append(idx)
            else:
                seen[var_key] = idx
                keep_indices.append(idx)

        if dup_indices:
            if skip_first:
                print(f"  Removed 1 boundary variant + {len(dup_indices)-1} exact duplicates", flush=True)
            else:
                print(f"  Found {len(dup_indices)} exact duplicate variants, removing rows/cols...", flush=True)
            print(f"  Matrix will be reduced from {shape[0]}x{shape[0]} to {len(keep_indices)}x{len(keep_indices)}", flush=True)
        else:
            print(f"  No exact duplicates found ({time.time()-step_time:.1f}s)", flush=True)
            keep_indices = None

    # Check and fix diagonal
    diag_mask = row == col
    unique_diag = np.unique(values[diag_mask])

    if len(unique_diag) == 1 and unique_diag[0] == 0.5:
        values = values.copy()
        values[diag_mask] = 1.0

    # Reconstruct sparse matrix
    sparse_mat = coo_matrix((values, (row, col)), shape=shape)

    # Convert to dense matrix
    print(f"[3/5] Converting sparse → dense...", flush=True)
    step_time = time.time()
    dense_mat = sparse_mat.toarray().astype(np.float32)
    print(f"  Done ({time.time()-step_time:.1f}s)", flush=True)

    # Remove duplicate rows/columns if needed
    if keep_indices is not None:
        print(f"[4/5] Removing duplicate rows/columns from matrix...", flush=True)
        step_time = time.time()
        dense_mat = dense_mat[np.ix_(keep_indices, keep_indices)]
        print(f"  Done ({time.time()-step_time:.1f}s)", flush=True)
    else:
        print(f"[4/5] No duplicates to remove", flush=True)

    # Mirror triangular matrix
    print(f"[5/5] Creating symmetric matrix and writing...", flush=True)
    step_time = time.time()

    upper_count = np.sum(row < col)
    lower_count = np.sum(row > col)

    if lower_count > 0 and upper_count == 0:
        lower_indices = np.tril_indices_from(dense_mat, k=-1)
        dense_mat[lower_indices[1], lower_indices[0]] = dense_mat[lower_indices[0], lower_indices[1]]
    elif upper_count > 0 and lower_count == 0:
        upper_indices = np.triu_indices_from(dense_mat, k=1)
        dense_mat[upper_indices[1], upper_indices[0]] = dense_mat[upper_indices[0], upper_indices[1]]

    # Save as binary float32
    dense_mat.tofile(output_path)

    file_size_mb = dense_mat.nbytes / (1024 * 1024)
    print(f"  {file_size_mb:.1f} MB written ({time.time()-step_time:.1f}s)", flush=True)

    # Create .vars file if variants_path provided
    if variants_path:
        vars_output = output_path + '.vars'
        create_vars_file(variants_path, vars_output, shape[0], keep_indices)

    total_time = time.time() - start_time
    print(f"Completed in {total_time:.1f}s ({file_size_mb/total_time:.1f} MB/s)\n", flush=True)

    return dense_mat.shape

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python npz_to_binary_matrix.py input.npz variants.gz output.bin")
        print("")
        print("Arguments:")
        print("  input.npz    : Input NPZ file containing sparse LD matrix")
        print("  variants.gz  : Variant list file (rsid, chrom, pos, ref, alt)")
        print("  output.bin   : Output binary matrix file")
        sys.exit(1)

    npz_path = sys.argv[1]
    variants_path = sys.argv[2]
    output_path = sys.argv[3]

    if not os.path.exists(npz_path):
        print(f"ERROR: NPZ file not found: {npz_path}")
        sys.exit(1)

    if not os.path.exists(variants_path):
        print(f"ERROR: Variants file not found: {variants_path}")
        sys.exit(1)

    convert_npz_to_binary(npz_path, output_path, variants_path)
