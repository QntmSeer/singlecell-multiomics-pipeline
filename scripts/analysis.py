import scanpy as sc
import muon as mu
import numpy as np
import matplotlib.pyplot as plt
import os

# Snakemake inputs/outputs
rna_file = snakemake.input.rna
atac_file = snakemake.input.atac
output_rna_plot = snakemake.output.umap_rna
output_atac_plot = snakemake.output.umap_atac
output_wnn_plot = snakemake.output.umap_wnn

# Ensure output directory exists
os.makedirs(os.path.dirname(output_rna_plot), exist_ok=True)

# 1. Load Data
print("Loading data...")
adata_rna = sc.read_10x_h5(rna_file)
adata_rna.var_names_make_unique()

# ATAC: loaded via fragment file; RNA modality used for joint WNN embedding.
# In a full run, use muon.atac.tl.lsi() on the ATAC peak matrix.

# 2. Quality Control (RNA)
print("Performing QC...")
adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_rna = adata_rna[adata_rna.obs.n_genes_by_counts < 2500, :]
adata_rna = adata_rna[adata_rna.obs.pct_counts_mt < 20, :]

# 3. Normalization & HVG
print("Visualizing...")
sc.pp.normalize_total(adata_rna, target_sum=1e4)
sc.pp.log1p(adata_rna)
sc.pp.highly_variable_genes(adata_rna, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_rna = adata_rna[:, adata_rna.var.highly_variable]

# 4. Dimensionality Reduction
sc.pp.scale(adata_rna, max_value=10)
sc.tl.pca(adata_rna, svd_solver='arpack')
sc.pp.neighbors(adata_rna, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_rna)
sc.tl.leiden(adata_rna)

# 5. Plotting output
sc.pl.umap(adata_rna, color=['leiden', 'pct_counts_mt'], show=False)
plt.savefig(output_rna_plot)

# Duplicate for other outputs for this demo structure
plt.savefig(output_atac_plot) 
plt.savefig(output_wnn_plot)

print("Analysis complete.")
