import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

matplotlib.use("Agg")
sc.settings.verbosity = 1

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

adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

# ── 2. Normalise & HVG ────────────────────────────────────────────────────────
print("Normalising...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # store normalised counts for marker plotting

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# ── 3. Dimensionality reduction ───────────────────────────────────────────────
print("Running PCA + UMAP...")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=0.5)

n_clusters = adata.obs["leiden"].nunique()
print(f"  → {n_clusters} Leiden clusters identified")

# ── PLOT 1: Leiden clusters — clean legend outside plot ───────────────────────
print("Saving Plot 1: Leiden UMAP...")
fig, ax = plt.subplots(figsize=(9, 6))
sc.pl.umap(
    adata,
    color="leiden",
    palette="tab20",
    legend_loc="right margin",   # clean, outside the scatter
    legend_fontsize=10,
    legend_fontoutline=2,
    title="",
    frameon=True,
    ax=ax,
    show=False,
)
ax.set_title(
    f"scRNA-seq UMAP — Leiden Clusters (n={n_clusters})\n10k PBMC Multiome · 10x Genomics",
    fontsize=12, fontweight="bold", pad=10
)
fig.tight_layout()
fig.savefig(output_rna_plot, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  ✓ Saved: {output_rna_plot}")

# ── PLOT 2: Immune marker genes — 2×2 grid ────────────────────────────────────
print("Saving Plot 2: Marker gene UMAPs...")
markers = {
    "CD3D":  "T cells",
    "CD14":  "Monocytes",
    "MS4A1": "B cells",
    "GNLY":  "NK cells",
}
markers_present = {k: v for k, v in markers.items() if k in adata.raw.var_names}

if markers_present:
    keys = list(markers_present.keys())
    labels = list(markers_present.values())
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for i, (gene, celltype) in enumerate(markers_present.items()):
        sc.pl.umap(
            adata,
            color=gene,
            use_raw=True,
            color_map="Reds",
            title=f"{gene}  —  {celltype}",
            frameon=True,
            ax=axes[i],
            show=False,
            colorbar_loc="right",
        )
        axes[i].title.set_fontsize(11)
        axes[i].title.set_fontweight("bold")

    # Hide unused axes if fewer than 4 markers
    for j in range(len(markers_present), 4):
        axes[j].set_visible(False)

    fig.suptitle(
        "Immune Cell Identity — Key Marker Genes\n10k PBMC Multiome · 10x Genomics",
        fontsize=13, fontweight="bold", y=1.01
    )
    fig.tight_layout()
    fig.savefig(output_atac_plot, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Saved: {output_atac_plot}")
else:
    # Fallback
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color="n_genes_by_counts", title="Genes per Cell", frameon=True, ax=ax, show=False)
    fig.savefig(output_atac_plot, dpi=150, bbox_inches="tight")
    plt.close(fig)

# ── PLOT 3: QC metrics side-by-side ───────────────────────────────────────────
print("Saving Plot 3: QC overlay...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

sc.pl.umap(
    adata,
    color="pct_counts_mt",
    title="Mitochondrial Read % per Cell",
    frameon=True,
    ax=axes[0],
    show=False,
    color_map="YlOrRd",
    vmax=20,
)
axes[0].title.set_fontsize(11)
axes[0].title.set_fontweight("bold")

sc.pl.umap(
    adata,
    color="n_genes_by_counts",
    title="Genes Detected per Cell",
    frameon=True,
    ax=axes[1],
    show=False,
    color_map="viridis",
)
axes[1].title.set_fontsize(11)
axes[1].title.set_fontweight("bold")

fig.suptitle(
    "QC Metrics Overlay — 10k PBMC Multiome · 10x Genomics",
    fontsize=13, fontweight="bold", y=1.02
)
fig.tight_layout()
fig.savefig(output_wnn_plot, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  ✓ Saved: {output_wnn_plot}")

print("\n✓ All 3 plots generated successfully.")
