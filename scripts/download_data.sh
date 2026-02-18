#!/bin/bash

# Script to download 10k Human PBMCs, Multiome v1.0, Chromium X dataset
# Source: 10x Genomics

DATA_DIR="data"
mkdir -p $DATA_DIR

echo "Downloading 10k PBMC Multiome dataset to $DATA_DIR..."

# Filtered Feature Barcode Matrix (HDF5)
echo "Downloading Filtered Feature Barcode Matrix..."
curl -o "$DATA_DIR/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5" \
    https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5

# ATAC Fragments file
echo "Downloading ATAC Fragments..."
curl -o "$DATA_DIR/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz" \
    https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz

# ATAC Fragments index
echo "Downloading ATAC Fragments Index..."
curl -o "$DATA_DIR/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi" \
    https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi

echo "Download complete."
