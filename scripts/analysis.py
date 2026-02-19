import scanpy as sc
import muon as mu
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

matplotlib.use("Agg")
sc.settings.verbosity = 1
mu.set_seed(42)

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
# Try reading as MuData (RNA + ATAC)
try:
    mdata = mu.read_10x_h5(rna_file)
    mdata.var_names_make_unique()
except Exception as e:
    log(f"Error reading with Muon: {e}. Fallback to Scanpy RNA only.")
    mdata = mu.MuData({'rna': sc.read_10x_h5(rna_file)})
    mdata.var_names_make_unique()

# Check modalities
has_atac = 'atac' in mdata.mod
has_rna = 'rna' in mdata.mod

log(f"Modalities found: {list(mdata.mod.keys())}")

# ── 2. RNA Processing ────────────────────────────────────────────────────────
if has_rna:
    log("Processing RNA...")
    rna = mdata['rna']
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    
    # Filter
    sc.pp.filter_cells(rna, min_genes=200)
    sc.pp.filter_genes(rna, min_cells=3)
    rna = rna[rna.obs.pct_counts_mt < 20, :]
    
    # Norm
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    rna.raw = rna # Save raw for markers
    
    # HVG + PCA
    sc.pp.highly_variable_genes(rna, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.pca(rna, svd_solver="arpack")
    
    # RNA-only UMAP (for comparison)
    sc.pp.neighbors(rna, n_neighbors=15, n_pcs=40)
    sc.tl.umap(rna)
    sc.tl.leiden(rna, resolution=0.5, key_added='leiden_rna')

# ── 3. ATAC Processing (LSI) ─────────────────────────────────────────────────
if has_atac:
    log("Processing ATAC...")
    atac = mdata['atac']
    # Filter highly abundant peaks/cells if needed, but standard LSI is robust
    sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
    
    # LSI (TF-IDF + SVD)
    mu.atac.pp.tfidf(atac, scale_factor=1e4)
    mu.atac.tl.lsi(atac)
    
    # ATAC-only UMAP
    sc.pp.neighbors(atac, use_rep='X_lsi', n_neighbors=15, n_pcs=40)
    sc.tl.umap(atac)
    sc.tl.leiden(atac, resolution=0.5, key_added='leiden_atac')
else:
    log("WARNING: No ATAC modality found. Skipping ATAC steps.")

# ── 4. WNN Integration ───────────────────────────────────────────────────────
if has_rna and has_atac:
    log("Running WNN integration...")
    mu.pp.neighbors(mdata, key_added='wnn')
    mu.tl.umap(mdata, neighbors_key='wnn')
    mu.tl.leiden(mdata, neighbors_key='wnn', key_added='leiden_wnn', resolution=0.5)
    
    # Use WNN Leiden clusters for all plots
    final_clusters = 'leiden_wnn'
else:
    log("Skipping WNN (missing modalities). Using RNA clusters.")
    final_clusters = 'leiden_rna' if has_rna else None

# Helper to save plot
def save_plot(fig, path, title):
    if fig is None: return
    # If figure was returned
    if isinstance(fig, matplotlib.figure.Figure):
        fig.suptitle(title, y=1.02, fontsize=14, fontweight='bold')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    else:
        # If sc.pl returns Axes or list of Axes, handle current fig
        plt.suptitle(title, y=1.02, fontsize=14, fontweight='bold')
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()

# ── PLOT 1: RNA UMAP (colored by WNN clusters) ────────────────────────────────
log("Saving RNA UMAP...")
if has_rna:
    fig = sc.pl.umap(
        mdata['rna'], 
        color=final_clusters, 
        legend_loc='on data', 
        title="RNA UMAP (WNN Clusters)",
        return_fig=True,
        show=False
    )
    save_plot(fig, output_rna_plot, "RNA Modality UMAP")
else:
    # Placeholder
    fig = plt.figure()
    plt.text(0.5, 0.5, "No RNA Data", ha='center')
    save_plot(fig, output_rna_plot, "No RNA Data")

# ── PLOT 2: ATAC UMAP (colored by WNN clusters) ──────────────────────────────
log("Saving ATAC UMAP...")
if has_atac:
    fig = sc.pl.umap(
        mdata['atac'], 
        color=final_clusters, 
        legend_loc='on data', 
        title="ATAC UMAP (LSI)", 
        return_fig=True,
        show=False
    )
    save_plot(fig, output_atac_plot, "ATAC Modality UMAP")
else:
     # Fallback if no ATAC
    fig = plt.figure()
    plt.text(0.5, 0.5, "No ATAC Data or LSI Failed", ha='center')
    save_plot(fig, output_atac_plot, "No ATAC Data")

# ── PLOT 3: WNN UMAP ─────────────────────────────────────────────────────────
log("Saving WNN UMAP...")
if has_rna and has_atac:
    # Muon stores WNN UMAP in mdata.obsm['X_umap'] usually if run on mdata
    # But let's check keys
    fig = mu.pl.umap(
        mdata, 
        color=final_clusters, 
        legend_loc='right margin', 
        title="WNN Integrated UMAP", 
        return_fig=True,
        show=False
    )
    save_plot(fig, output_wnn_plot, "Weighted Nearest Neighbors (WNN) Integration")
else:
    # Fallback to QC plot or RNA
    fig = plt.figure()
    plt.text(0.5, 0.5, "WNN Requires Both RNA & ATAC", ha='center')
    save_plot(fig, output_wnn_plot, "WNN Integration Failed")

# ── PLOT 4: Marker Genes (on WNN or RNA UMAP) ────────────────────────────────
log("Saving Marker Gene Plot...")
markers = {"CD3D": "T cells", "CD14": "Monocytes", "MS4A1": "B cells", "GNLY": "NK cells"}
if has_rna:
    # Use RNA expression but overlay on WNN UMAP if available, else RNA UMAP
    # To plot genes on mdata UMAP, we need to ensure mdata has the expression in .X or .raw
    # mdata['rna'].obs.index should match mdata.obs.index
    
    basis = 'umap' if (has_rna and has_atac) else 'X_umap' # mdata.obsm key
    # Default mu.pl.umap looks for 'X_umap' in mdata.obsm
    
    # We want to plot RNA genes. 
    # easiest is to plot on mdata['rna'] but use WNN coordinates if possible?
    # Or just plot on RNA UMAP. The user liked the previous plot.
    # Let's keep it simple: Plot markers on RNA UMAP for clarity of expression.
    # OR better: Plot on WNN UMAP to show how they separate there.
    # If WNN exists, copy WNN coordinates to RNA object for plotting
    if has_rna and has_atac:
        mdata['rna'].obsm['X_umap_wnn'] = mdata.obsm['X_umap']
        use_basis = 'X_umap_wnn'
        title_suffix = "(WNN Coordinates)"
    else:
        use_basis = 'X_umap'
        title_suffix = "(RNA Coordinates)"

    # Check matches
    present = {k:v for k,v in markers.items() if k in mdata['rna'].raw.var_names}
    if present:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        for i, (g, desc) in enumerate(present.items()):
            sc.pl.embedding(
                mdata['rna'], 
                basis=use_basis, 
                color=g, 
                use_raw=True, 
                color_map='Reds', 
                ax=axes[i], 
                show=False, 
                frameon=False,
                title=f"{g} · {desc}"
            )
        for j in range(len(present), 4): axes[j].axis('off')
        save_plot(fig, output_markers_plot, f"Immune Markers {title_suffix}")
    else:
        fig = plt.figure()
        plt.text(0.5,0.5,"Markers not found", ha='center')
        save_plot(fig, output_markers_plot, "Marker Genes Missing")

log("✓ Done.")
