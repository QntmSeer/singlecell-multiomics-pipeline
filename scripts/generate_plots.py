"""
Standalone plot generator — run directly with:
    python scripts/generate_plots.py

No Snakemake needed. Reads from data/ and writes to results/plots/
"""
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

DATA_FILE = "data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
OUT_DIR   = "results/plots"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load ──────────────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_10x_h5(DATA_FILE)
adata.var_names_make_unique()

# ── QC ────────────────────────────────────────────────────────────────────────
print("QC filtering...")
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts > 200]
adata = adata[adata.obs.n_genes_by_counts < 6000]
adata = adata[adata.obs.pct_counts_mt < 20]
print(f"  → {adata.n_obs} cells after QC")

# ── Normalise ─────────────────────────────────────────────────────────────────
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# ── Dim reduction ─────────────────────────────────────────────────────────────
print("PCA + UMAP...")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=0.5)
n_clusters = adata.obs["leiden"].nunique()
print(f"  → {n_clusters} Leiden clusters")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Leiden clusters (clean right-margin legend)
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Plot 1: Leiden clusters...")
fig, ax = plt.subplots(figsize=(10, 7))
sc.pl.umap(
    adata,
    color="leiden",
    palette="tab20",
    legend_loc="right margin",
    legend_fontsize=11,
    legend_fontoutline=2,
    frameon=True,
    title=f"scRNA-seq UMAP — Leiden Clusters (n={n_clusters})\n10k PBMC Multiome · 10x Genomics",
    ax=ax,
    show=False,
)
fig.savefig(f"{OUT_DIR}/umap_rna.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  ✓ {OUT_DIR}/umap_rna.png")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Immune marker genes (2×2 grid, Reds colormap)
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Plot 2: Marker genes...")
markers = {"CD3D": "T cells", "CD14": "Monocytes", "MS4A1": "B cells", "GNLY": "NK cells"}
present = {g: c for g, c in markers.items() if g in adata.raw.var_names}
print(f"  → Markers found: {list(present.keys())}")

fig, axes = plt.subplots(2, 2, figsize=(13, 10))
axes = axes.flatten()

for i, (gene, celltype) in enumerate(present.items()):
    sc.pl.umap(
        adata,
        color=gene,
        use_raw=True,
        color_map="Reds",
        frameon=True,
        ax=axes[i],
        show=False,
        colorbar_loc="right",
    )
    axes[i].set_title(f"{gene}  ·  {celltype}", fontsize=12, fontweight="bold", pad=6)

for j in range(len(present), 4):
    axes[j].set_visible(False)

fig.suptitle(
    "Immune Cell Identity — Key Marker Genes\n10k PBMC Multiome · 10x Genomics",
    fontsize=14, fontweight="bold"
)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{OUT_DIR}/umap_atac.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  ✓ {OUT_DIR}/umap_atac.png")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 3 — QC metrics side-by-side (mito% + n_genes)
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Plot 3: QC overlay...")
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

sc.pl.umap(
    adata, color="pct_counts_mt",
    title="Mitochondrial Read % per Cell",
    frameon=True, ax=axes[0], show=False,
    color_map="YlOrRd", vmax=20,
)
axes[0].title.set_fontsize(12)
axes[0].title.set_fontweight("bold")

sc.pl.umap(
    adata, color="n_genes_by_counts",
    title="Genes Detected per Cell",
    frameon=True, ax=axes[1], show=False,
    color_map="viridis",
)
axes[1].title.set_fontsize(12)
axes[1].title.set_fontweight("bold")

fig.suptitle(
    "QC Metrics Overlay — 10k PBMC Multiome · 10x Genomics",
    fontsize=14, fontweight="bold"
)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{OUT_DIR}/umap_wnn.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  ✓ {OUT_DIR}/umap_wnn.png")

print("\n✓ Done. Check results/plots/")
