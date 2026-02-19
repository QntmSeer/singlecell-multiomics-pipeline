import scanpy as sc
import muon as mu
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

matplotlib.use("Agg")
sc.settings.verbosity = 1
np.random.seed(42)  # Fix: muon.set_seed removed in v0.1.5

# Snakemake inputs/outputs
rna_file = snakemake.input.rna
atac_file = snakemake.input.atac
output_rna_plot = snakemake.output.umap_rna
output_atac_plot = snakemake.output.umap_atac
output_wnn_plot = snakemake.output.umap_wnn
output_markers_plot = snakemake.output.umap_markers

for p in [output_rna_plot, output_atac_plot, output_wnn_plot, output_markers_plot]:
    os.makedirs(os.path.dirname(p), exist_ok=True)

# ── 1. Load Data (Muon) ──────────────────────────────────────────────────────────
console_log = []
def log(msg):
    print(msg)
    console_log.append(msg)

log("Loading h5 data with Muon...")
try:
    mdata = mu.read_10x_h5(rna_file)
    mdata.var_names_make_unique()
except Exception as e:
    log(f"Error reading with Muon: {e}. Fallback to Scanpy RNA only.")
    mdata = mu.MuData({'rna': sc.read_10x_h5(rna_file)})
    mdata.var_names_make_unique()

has_atac = 'atac' in mdata.mod
has_rna = 'rna' in mdata.mod
log(f"Modalities found: {list(mdata.mod.keys())}")

# ── 2. QC & Filtering (Global) ───────────────────────────────────────────────
if has_rna:
    log("QC & Filtering...")
    rna = mdata['rna'] # this is a view
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    
    # Define cells to keep based on RNA metrics
    keep_cells = (rna.obs.n_genes_by_counts > 200) & \
                 (rna.obs.n_genes_by_counts < 6000) & \
                 (rna.obs.pct_counts_mt < 20)
    
    # Filter mdata globally to keep modalities synced!
    mdata = mdata[keep_cells, :].copy()
    log(f"  → {mdata.n_obs} cells after QC")

# ── 3. RNA Processing ────────────────────────────────────────────────────────
if has_rna:
    log("Processing RNA...")
    # Re-access modality from the subsetted mdata
    rna = mdata['rna']
    rna.var_names_make_unique()
    
    # Norm & Log1p (in-place on mdata['rna'])
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    rna.raw = rna # Save raw for markers
    
    # HVG & PCA
    sc.pp.highly_variable_genes(rna, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.pca(rna, svd_solver="arpack")
    
    # RNA Neighbors & UMAP
    sc.pp.neighbors(rna, n_neighbors=15, n_pcs=40)
    sc.tl.umap(rna)
    sc.tl.leiden(rna, resolution=0.5, key_added='leiden_rna')

# ── 3. ATAC Processing (LSI) ─────────────────────────────────────────────────
if has_atac:
    log("Processing ATAC...")
    atac = mdata['atac']
    
    # LSI (TF-IDF + SVD)
    mu.atac.pp.tfidf(atac, scale_factor=1e4)
    mu.atac.tl.lsi(atac)
    
    # ATAC Neighbors & UMAP
    sc.pp.neighbors(atac, use_rep='X_lsi', n_neighbors=15, n_pcs=40)
    sc.tl.umap(atac)
    sc.tl.leiden(atac, resolution=0.5, key_added='leiden_atac')

# ── 4. WNN Integration ───────────────────────────────────────────────────────
final_clusters = 'leiden_rna'
if has_rna and has_atac:
    log("Running WNN integration (default neighbors)...")
    mu.pp.neighbors(mdata) # stores to .uns['neighbors']
    mu.tl.umap(mdata)      # uses .uns['neighbors']
    mu.tl.leiden(mdata, key_added='leiden_wnn', resolution=0.5) # uses .uns['neighbors']
    
    # CRITICAL FIX: Copy WNN clusters to individual modalities for plotting
    # mdata.obs contains 'leiden_wnn', but mdata['rna'].obs might not automatically have it accessible via plotting
    if 'leiden_wnn' in mdata.obs:
        mdata['rna'].obs['leiden_wnn'] = mdata.obs['leiden_wnn']
        mdata['atac'].obs['leiden_wnn'] = mdata.obs['leiden_wnn']
        
    final_clusters = 'leiden_wnn'

# Helper to save plot
def save_plot(fig, path, title):
    if fig is None: return
    if isinstance(fig, matplotlib.figure.Figure):
        fig.suptitle(title, y=1.02, fontsize=14, fontweight='bold')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.suptitle(title, y=1.02, fontsize=14, fontweight='bold')
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()

# ── PLOTS ────────────────────────────────────────────────────────────────────

# 1. RNA UMAP
log("Saving RNA UMAP...")
if has_rna:
    try:
        fig = sc.pl.umap(mdata['rna'], color=final_clusters, title="RNA UMAP", return_fig=True, show=False)
        save_plot(fig, output_rna_plot, "RNA Modality")
    except Exception as e:
        log(f"Failed to plot RNA UMAP: {e}")
        save_plot(plt.figure(), output_rna_plot, "Plotting Error")
else:
    save_plot(plt.figure(), output_rna_plot, "No RNA")

# 2. ATAC UMAP
log("Saving ATAC UMAP...")
if has_atac:
    try:
        fig = sc.pl.umap(mdata['atac'], color=final_clusters, title="ATAC UMAP", return_fig=True, show=False)
        save_plot(fig, output_atac_plot, "ATAC Modality")
    except Exception as e:
        log(f"Failed to plot ATAC UMAP: {e}")
        save_plot(plt.figure(), output_atac_plot, "Plotting Error")
else:
    save_plot(plt.figure(), output_atac_plot, "No ATAC")

# 3. WNN UMAP
log("Saving WNN UMAP...")
if has_rna and has_atac:
    try:
        # mu.pl.umap doesn't always return figure object same way as scanpy
        # It usually plots using mdata.obsm['X_umap']
        mu.pl.umap(mdata, color=final_clusters, title="WNN UMAP", show=False)
        # Grab current figure since muon might not return it
        fig = plt.gcf()
        save_plot(fig, output_wnn_plot, "WNN Integrated")
    except Exception as e:
        log(f"Failed to plot WNN UMAP: {e}")
        save_plot(plt.figure(), output_wnn_plot, "Plotting Error")
else:
    save_plot(plt.figure(), output_wnn_plot, "WNN Failed")

# 4. Marker Genes
log("Saving Markers...")
markers = {"CD3D": "T cells", "CD14": "Monocytes", "MS4A1": "B cells", "GNLY": "NK cells"}
if has_rna:
    # Use RNA UMAP coordinates for clarity if WNN not available, or copy WNN coords
    use_basis = 'X_umap'
    if has_rna and has_atac:
        # Copy WNN UMAP to RNA object explicitly
        if 'X_umap' in mdata.obsm:
            mdata['rna'].obsm['X_umap_wnn'] = mdata.obsm['X_umap']
            use_basis = 'X_umap_wnn'

    present = {k:v for k,v in markers.items() if k in mdata['rna'].raw.var_names}
    if present:
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        axes = axes.flatten()
        for i, (g, desc) in enumerate(present.items()):
            sc.pl.embedding(mdata['rna'], basis=use_basis, color=g, use_raw=True, 
                          color_map='Reds', ax=axes[i], show=False, frameon=False, title=f"{g} ({desc})")
        for j in range(len(present), 4): axes[j].axis('off')
        save_plot(fig, output_markers_plot, "Immune Markers")
    else:
        save_plot(plt.figure(), output_markers_plot, "No Markers Found")
else:
    save_plot(plt.figure(), output_markers_plot, "No RNA Data")

log("✓ Done.")
