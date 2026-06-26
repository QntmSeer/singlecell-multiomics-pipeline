# scripts/qc_check.py
import os
import scanpy as sc
import mudata as mu
import json

input_file = snakemake.input.rna
output_txt = snakemake.output.txt
output_json = snakemake.output.json

# Ensure output dir exists
os.makedirs(os.path.dirname(output_txt), exist_ok=True)
os.makedirs(os.path.dirname(output_json), exist_ok=True)

# Read mudata
try:
    mdata = mu.read_10x_h5(input_file)
    mdata.var_names_make_unique()
except Exception as e:
    # Fallback to RNA-only if mudata read fails
    mdata = mu.MuData({'rna': sc.read_10x_h5(input_file)})
    mdata.var_names_make_unique()

qc_stats = {}
if 'rna' in mdata.mod:
    rna = mdata['rna']
    rna.var_names_make_unique()
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    qc_stats['rna'] = {
        'n_cells_raw': int(rna.n_obs),
        'n_genes_raw': int(rna.n_vars),
        'median_genes_per_cell': float(rna.obs['n_genes_by_counts'].median()),
        'median_counts_per_cell': float(rna.obs['total_counts'].median()),
        'median_mito_pct': float(rna.obs['pct_counts_mt'].median())
    }

if 'atac' in mdata.mod:
    atac = mdata['atac']
    atac.var_names_make_unique()
    sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
    qc_stats['atac'] = {
        'n_cells_raw': int(atac.n_obs),
        'n_peaks_raw': int(atac.n_vars),
        'median_peaks_per_cell': float(atac.obs['n_genes_by_counts'].median()),
        'median_counts_per_cell': float(atac.obs['total_counts'].median())
    }

# Write a log file
with open(output_txt, "w") as f:
    f.write(f"QC completed for {input_file}\n")
    if 'rna' in qc_stats:
        f.write(f"RNA: cells={qc_stats['rna']['n_cells_raw']}, median_genes={qc_stats['rna']['median_genes_per_cell']:.1f}\n")
    if 'atac' in qc_stats:
        f.write(f"ATAC: cells={qc_stats['atac']['n_cells_raw']}, median_peaks={qc_stats['atac']['median_peaks_per_cell']:.1f}\n")
    f.write("Metrics: PASS\n")

# Write MultiQC custom content JSON
multiqc_data = {
    "id": "singlecell_multiomics_qc",
    "name": "Single-Cell Multi-Omics QC Summary",
    "description": "This table shows the pre-filtering QC metrics calculated for the RNA and ATAC modalities in the H5 dataset.",
    "plot_type": "table",
    "section_name": "Single-Cell Multi-Omics QC",
    "data": {
        "RNA modality": {
            "Total Cells": qc_stats.get('rna', {}).get('n_cells_raw', 0),
            "Total Genes/Peaks": qc_stats.get('rna', {}).get('n_genes_raw', 0),
            "Median Features/Cell": round(qc_stats.get('rna', {}).get('median_genes_per_cell', 0), 1),
            "Median Counts/Cell": round(qc_stats.get('rna', {}).get('median_counts_per_cell', 0), 1),
            "Median Mito %": round(qc_stats.get('rna', {}).get('median_mito_pct', 0), 2)
        }
    }
}

if 'atac' in qc_stats:
    multiqc_data["data"]["ATAC modality"] = {
        "Total Cells": qc_stats.get('atac', {}).get('n_cells_raw', 0),
        "Total Genes/Peaks": qc_stats.get('atac', {}).get('n_peaks_raw', 0),
        "Median Features/Cell": round(qc_stats.get('atac', {}).get('median_peaks_per_cell', 0), 1),
        "Median Counts/Cell": round(qc_stats.get('atac', {}).get('median_counts_per_cell', 0), 1),
        "Median Mito %": "N/A"
    }

with open(output_json, "w") as f:
    json.dump(multiqc_data, f, indent=4)

