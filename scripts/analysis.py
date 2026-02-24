"""
Single-Cell Multi-Omics Analysis Pipeline
==========================================
Performs joint RNA + ATAC analysis using Muon/Scanpy.

Stages:
    1. Data loading (10x Multiome HDF5)
    2. QC & filtering
    3. RNA processing (normalisation, HVG, PCA, UMAP, Leiden clustering)
    4. ATAC processing (TF-IDF, LSI, UMAP, Leiden clustering)
    5. Weighted Nearest Neighbor (WNN) multi-modal integration
    6. Plotting (RNA, ATAC, WNN UMAPs & immune marker expression)
    7. Trajectory inference on Monocyte subset (Diffusion Pseudotime)
"""

import scanpy as sc
import muon as mu
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

matplotlib.use("Agg")
sc.settings.verbosity = 1
np.random.seed(42)

# ── I/O from Snakemake ────────────────────────────────────────────────────────
rna_file             = snakemake.input.rna
atac_file            = snakemake.input.atac
output_rna_plot      = snakemake.output.umap_rna
output_atac_plot     = snakemake.output.umap_atac
output_wnn_plot      = snakemake.output.umap_wnn
output_markers_plot  = snakemake.output.umap_markers
output_paga_plot     = snakemake.output.trajectory_paga
output_pseudotime_plot = snakemake.output.trajectory_pseudotime

for p in [output_rna_plot, output_atac_plot, output_wnn_plot,
          output_markers_plot, output_paga_plot, output_pseudotime_plot]:
    os.makedirs(os.path.dirname(p), exist_ok=True)

# ── Logging ───────────────────────────────────────────────────────────────────
console_log = []
def log(msg):
    print(msg)
    console_log.append(msg)

# ── 1. Load Data ──────────────────────────────────────────────────────────────
log("Loading 10x Multiome data...")
try:
    mdata = mu.read_10x_h5(rna_file)
    mdata.var_names_make_unique()
except Exception as e:
    log(f"Muon read failed ({e}). Falling back to RNA-only mode.")
    mdata = mu.MuData({'rna': sc.read_10x_h5(rna_file)})
    mdata.var_names_make_unique()

has_atac = 'atac' in mdata.mod
has_rna  = 'rna'  in mdata.mod
log(f"Modalities detected: {list(mdata.mod.keys())}")

# ── 2. QC & Filtering ─────────────────────────────────────────────────────────
if has_rna:
    log("Running QC...")
    rna = mdata['rna']
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], percent_top=None,
                                log1p=False, inplace=True)
    keep_cells = (
        (rna.obs.n_genes_by_counts > 200)  &
        (rna.obs.n_genes_by_counts < 6000) &
        (rna.obs.pct_counts_mt    < 20)
    )
    # Filter globally to keep all modalities in sync
    mdata = mdata[keep_cells, :].copy()
    log(f"  {mdata.n_obs} cells retained after QC")

# ── 3. RNA Processing ─────────────────────────────────────────────────────────
if has_rna:
    log("Processing RNA modality...")
    rna = mdata['rna']
    rna.var_names_make_unique()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    rna.raw = rna  # preserve pre-HVG expression for marker plotting
    sc.pp.highly_variable_genes(rna, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.pca(rna, svd_solver="arpack")
    sc.pp.neighbors(rna, n_neighbors=15, n_pcs=40)
    sc.tl.umap(rna)
    sc.tl.leiden(rna, resolution=0.5, key_added='leiden_rna')

# ── 4. ATAC Processing (LSI) ──────────────────────────────────────────────────
if has_atac:
    log("Processing ATAC modality (LSI)...")
    atac = mdata['atac']
    mu.atac.pp.tfidf(atac, scale_factor=1e4)
    mu.atac.tl.lsi(atac)
    sc.pp.neighbors(atac, use_rep='X_lsi', n_neighbors=15, n_pcs=40)
    sc.tl.umap(atac)
    sc.tl.leiden(atac, resolution=0.5, key_added='leiden_atac')

# ── 5. WNN Integration ────────────────────────────────────────────────────────
final_clusters = 'leiden_rna'
if has_rna and has_atac:
    log("Running WNN integration...")
    mu.pp.neighbors(mdata)
    mu.tl.umap(mdata)
    mu.tl.leiden(mdata, key_added='leiden_wnn', resolution=0.5)

    # Propagate joint clusters to individual modalities for downstream plotting
    if 'leiden_wnn' in mdata.obs:
        mdata['rna'].obs['leiden_wnn']  = mdata.obs['leiden_wnn']
        mdata['atac'].obs['leiden_wnn'] = mdata.obs['leiden_wnn']

    final_clusters = 'leiden_wnn'

# ── Helper: Save Figure ───────────────────────────────────────────────────────
def save_plot(fig, path, title):
    if fig is None:
        return
    if isinstance(fig, matplotlib.figure.Figure):
        fig.suptitle(title, y=1.02, fontsize=14, fontweight='bold')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.suptitle(title, y=1.02, fontsize=14, fontweight='bold')
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()

# ── 6. Plots ──────────────────────────────────────────────────────────────────

# 6a. RNA UMAP
log("Saving RNA UMAP...")
if has_rna:
    try:
        fig = sc.pl.umap(mdata['rna'], color=final_clusters,
                         title="RNA UMAP", return_fig=True, show=False)
        save_plot(fig, output_rna_plot, "RNA Modality")
    except Exception as e:
        log(f"RNA UMAP failed: {e}")
        save_plot(plt.figure(), output_rna_plot, "Plotting Error")
else:
    save_plot(plt.figure(), output_rna_plot, "No RNA")

# 6b. ATAC UMAP
log("Saving ATAC UMAP...")
if has_atac:
    try:
        fig = sc.pl.umap(mdata['atac'], color=final_clusters,
                         title="ATAC UMAP", return_fig=True, show=False)
        save_plot(fig, output_atac_plot, "ATAC Modality")
    except Exception as e:
        log(f"ATAC UMAP failed: {e}")
        save_plot(plt.figure(), output_atac_plot, "Plotting Error")
else:
    save_plot(plt.figure(), output_atac_plot, "No ATAC")

# 6c. WNN UMAP
log("Saving WNN UMAP...")
if has_rna and has_atac:
    try:
        mu.pl.umap(mdata, color=final_clusters, title="WNN UMAP", show=False)
        save_plot(plt.gcf(), output_wnn_plot, "WNN Integrated")
    except Exception as e:
        log(f"WNN UMAP failed: {e}")
        save_plot(plt.figure(), output_wnn_plot, "Plotting Error")
else:
    save_plot(plt.figure(), output_wnn_plot, "WNN Unavailable")

# 6d. Immune Marker Genes
log("Saving immune marker plot...")
markers = {"CD3D": "T cells", "CD14": "Monocytes", "MS4A1": "B cells", "GNLY": "NK cells"}
if has_rna:
    use_basis = 'X_umap'
    if has_rna and has_atac and 'X_umap' in mdata.obsm:
        mdata['rna'].obsm['X_umap_wnn'] = mdata.obsm['X_umap']
        use_basis = 'X_umap_wnn'

    present = {k: v for k, v in markers.items()
               if k in mdata['rna'].raw.var_names}
    if present:
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        axes = axes.flatten()
        for i, (gene, label) in enumerate(present.items()):
            sc.pl.embedding(mdata['rna'], basis=use_basis, color=gene,
                            use_raw=True, color_map='Reds', ax=axes[i],
                            show=False, frameon=False,
                            title=f"{gene} ({label})")
        for j in range(len(present), 4):
            axes[j].axis('off')
        save_plot(fig, output_markers_plot, "Immune Markers")
    else:
        save_plot(plt.figure(), output_markers_plot, "No Markers Found")
else:
    save_plot(plt.figure(), output_markers_plot, "No RNA Data")

# ── 7. Trajectory Inference: Monocyte Differentiation ─────────────────────────
# Subsets CD14+ cells and models the Classical → Non-Classical transition
# using Partition-based Graph Abstraction (PAGA) and Diffusion Pseudotime (DPT).

log("Starting trajectory inference (Monocyte subset)...")
if has_rna and 'leiden_wnn' in mdata.obs:
    try:
        rna = mdata['rna']
        if 'CD14' in rna.var_names:
            cd14_expr = rna[:, 'CD14'].X.toarray().flatten()
            is_mono   = cd14_expr > 0.5

            if np.sum(is_mono) > 50:
                mono = rna[is_mono, :].copy()
                log(f"  {mono.n_obs} monocytes selected for trajectory.")

                sc.pp.highly_variable_genes(mono, n_top_genes=2000)
                sc.pp.pca(mono)
                sc.pp.neighbors(mono, n_neighbors=20)
                sc.tl.leiden(mono, resolution=0.5, key_added='leiden_mono')

                n_clusters = len(mono.obs['leiden_mono'].unique())
                log(f"  Monocyte sub-clusters: {n_clusters}")

                # PAGA — abstract connectivity graph between sub-clusters
                if n_clusters > 1:
                    log("  Computing PAGA...")
                    try:
                        sc.tl.paga(mono, groups='leiden_mono')
                        sc.pl.paga(mono, color=['leiden_mono', 'CD14', 'FCGR3A'],
                                   show=False)
                        save_plot(plt.gcf(), output_paga_plot,
                                  "Monocyte PAGA Trajectory")
                    except Exception as e:
                        log(f"  PAGA skipped: {e}")
                        save_plot(plt.figure(), output_paga_plot, "PAGA Failed")
                else:
                    save_plot(plt.figure(), output_paga_plot,
                              "Single Cluster (No PAGA)")

                # Diffusion Pseudotime — root at cell with highest CD14 expression
                log("  Computing Diffusion Pseudotime...")
                cd14_idx = np.argmax(mono[:, 'CD14'].X.toarray().flatten())
                mono.uns['iroot'] = cd14_idx
                sc.tl.diffmap(mono)
                sc.tl.dpt(mono)

                # Plot: pseudotime + CD14 + FCGR3A side-by-side
                fig_dpt, ax = plt.subplots(1, 3, figsize=(15, 5))
                sc.pl.umap(mono, color='dpt_pseudotime', ax=ax[0],
                           show=False, title="Pseudotime")
                for i, gene in enumerate(['CD14', 'FCGR3A']):
                    if gene in mono.var_names:
                        sc.pl.umap(mono, color=gene, ax=ax[i + 1],
                                   show=False, title=gene)
                    elif mono.raw and gene in mono.raw.var_names:
                        sc.pl.umap(mono, color=gene, ax=ax[i + 1],
                                   show=False, title=gene, use_raw=True)
                    else:
                        ax[i + 1].text(0.5, 0.5, f"{gene} not found",
                                       ha='center', va='center')
                        ax[i + 1].axis('off')

                plt.tight_layout()
                save_plot(fig_dpt, output_pseudotime_plot,
                          "Monocyte Differentiation")
                log("  Trajectory analysis complete.")

            else:
                log(f"  Too few CD14+ cells (n={np.sum(is_mono)}); skipping.")
                save_plot(plt.figure(), output_paga_plot, "Insufficient Monocytes")
                save_plot(plt.figure(), output_pseudotime_plot,
                          "Insufficient Monocytes")
        else:
            log("  CD14 not detected in dataset; skipping trajectory.")
            save_plot(plt.figure(), output_paga_plot, "CD14 Missing")
            save_plot(plt.figure(), output_pseudotime_plot, "CD14 Missing")

    except Exception as e:
        import traceback
        log(f"Trajectory inference failed: {e}")
        traceback.print_exc()
        save_plot(plt.figure(), output_paga_plot, f"Error: {str(e)[:30]}")
        save_plot(plt.figure(), output_pseudotime_plot, "Trajectory Error")

else:
    log("Skipping trajectory (no RNA or WNN clusters detected).")
    save_plot(plt.figure(), output_paga_plot, "Skipped")
    save_plot(plt.figure(), output_pseudotime_plot, "Skipped")

log("✓ Analysis complete.")
