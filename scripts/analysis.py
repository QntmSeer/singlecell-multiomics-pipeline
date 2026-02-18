import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

matplotlib.use("Agg")

# Snakemake inputs/outputs
rna_file = snakemake.input.rna
atac_file = snakemake.input.atac
output_rna_plot = snakemake.output.umap_rna
output_atac_plot = snakemake.output.umap_atac
output_wnn_plot = snakemake.output.umap_wnn

os.makedirs(os.path.dirname(output_rna_plot), exist_ok=True)

# ── 1. Load & QC ──────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_10x_h5(rna_file)
adata.var_names_make_unique()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# Filter
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

# ── 2. Normalise & HVG ────────────────────────────────────────────────────────
print("Normalising...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Store raw (for marker gene plotting later)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# ── 3. Dimensionality reduction ───────────────────────────────────────────────
print("Running PCA + UMAP...")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=0.5)

# ── 4. Plot 1: RNA UMAP — Leiden clusters ─────────────────────────────────────
print("Saving RNA UMAP (Leiden clusters)...")
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(
    adata,
    color="leiden",
    title="scRNA-seq UMAP — Leiden Clusters\n10k PBMC Multiome (10x Genomics)",
    legend_loc="on data",
    legend_fontsize=9,
    frameon=False,
    ax=ax,
    show=False,
)
fig.tight_layout()
fig.savefig(output_rna_plot, dpi=150, bbox_inches="tight")
plt.close(fig)

# ── 5. Plot 2: Marker gene UMAP — immune cell identity ────────────────────────
print("Saving marker gene UMAP...")
# Key immune markers: CD3D=T-cell, CD14=Monocyte, MS4A1=B-cell, GNLY=NK cell
markers = ["CD3D", "CD14", "MS4A1", "GNLY"]
# Keep only markers present in the dataset
markers_present = [m for m in markers if m in adata.raw.var_names]

if markers_present:
    fig = sc.pl.umap(
        adata,
        color=markers_present,
        use_raw=True,
        ncols=2,
        title=[f"{m} (marker)" for m in markers_present],
        frameon=False,
        show=False,
        return_fig=True,
    )
    fig.suptitle("Immune Cell Marker Genes — PBMC Multiome", y=1.02, fontsize=12, fontweight="bold")
    fig.savefig(output_atac_plot, dpi=150, bbox_inches="tight")
    plt.close(fig)
else:
    # Fallback: colour by n_genes
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color="n_genes_by_counts", title="Genes per Cell", frameon=False, ax=ax, show=False)
    fig.savefig(output_atac_plot, dpi=150, bbox_inches="tight")
    plt.close(fig)

# ── 6. Plot 3: QC overlay UMAP — data quality landscape ───────────────────────
print("Saving QC overlay UMAP...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

sc.pl.umap(
    adata,
    color="pct_counts_mt",
    title="% Mitochondrial Reads",
    frameon=False,
    ax=axes[0],
    show=False,
    color_map="YlOrRd",
)
sc.pl.umap(
    adata,
    color="n_genes_by_counts",
    title="Genes Detected per Cell",
    frameon=False,
    ax=axes[1],
    show=False,
    color_map="viridis",
)

fig.suptitle("QC Metrics Overlay — 10k PBMC Multiome", fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(output_wnn_plot, dpi=150, bbox_inches="tight")
plt.close(fig)

print("✓ All 3 plots saved successfully.")
