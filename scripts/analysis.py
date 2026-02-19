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
# Phase 3: Trajectory Outputs
output_paga_plot = snakemake.output.trajectory_paga
output_pseudotime_plot = snakemake.output.trajectory_pseudotime

for p in [output_rna_plot, output_atac_plot, output_wnn_plot, output_markers_plot, output_paga_plot, output_pseudotime_plot]:
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

# ── PLOTS (Global) ──────────────────────────────────────────────────────────

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

# ── PHASE 3: TRAJECTORY INFERENCE (Monocytes) ────────────────────────────────

log("Starting Phase 3: Trajectory Inference...")
if has_rna and 'leiden_wnn' in mdata.obs:
    try:
        # 1. Identify Monocyte Clusters (High CD14)
        # We need to access RNA expression data for CD14
        # mdata['rna'].X might be sparse
        rna = mdata['rna']
        if 'CD14' in rna.var_names:
            # Score CD14 expression per cluster
            # Simple approach: Subset to CD14+ cells (expr > 0.5) then see which clusters they belong to
            # Better: Calculate mean expression per cluster
            
            # Create a small dataframe of cluster vs CD14
            # (Keep it simple and robust for this portfolio piece)
            
            # Strategy: Just subset broadly on CD14 expression first to capture all Monocytes
            # Then re-cluster to find the trajectory
            
            log("  Subsetting Monocytes (CD14+)...")
            
            # Get cells with CD14 expression > 0.5 (log1p scale)
            cd14_expr = rna[:, 'CD14'].X.toarray().flatten()
            is_mono = cd14_expr > 0.5
            
            if np.sum(is_mono) > 50: # Ensure we have enough cells
                mono = rna[is_mono, :].copy()
                log(f"  Found {mono.n_obs} potential Monocytes.")
                
                # Re-process Monocytes
                sc.pp.highly_variable_genes(mono, n_top_genes=2000)
                sc.pp.pca(mono)
                sc.pp.neighbors(mono, n_neighbors=20)
                sc.tl.leiden(mono, resolution=0.5, key_added='leiden_mono')
                
                # Check clusters
                n_clusters = len(mono.obs['leiden_mono'].unique())
                log(f"  Monocyte Clusters found: {n_clusters}")

                # PAGA (needs > 1 cluster connectivity)
                if n_clusters > 1:
                    log("  Running PAGA...")
                    try:
                        sc.tl.paga(mono, groups='leiden_mono')
                        # Plot PAGA
                        fig_paga = plt.figure()
                        sc.pl.paga(mono, color=['leiden_mono', 'CD14', 'FCGR3A'], show=False) 
                        save_plot(plt.gcf(), output_paga_plot, "Monocyte PAGA Trajectory")
                    except Exception as e:
                        log(f"  PAGA Failed (skipping graph): {e}")
                        save_plot(plt.figure(), output_paga_plot, "PAGA Failed")
                else:
                    log("  Skipping PAGA (Only 1 cluster found).")
                    save_plot(plt.figure(), output_paga_plot, "Single Cluster (No PAGA)")
                
                # Diffusion Pseudotime (DPT)
                # Set root to max CD14 cell (Classical Monocyte)
                log("  Calculating Pseudotime (DPT)...")
                if 'iroot' not in mono.uns:
                     # Find max CD14 cell index
                     # .X is sparse CSR matrix, verify shape
                     cd14_idx = np.argmax(mono[:, 'CD14'].X.toarray().flatten())
                     mono.uns['iroot'] = cd14_idx
                     log(f"  Root cell index: {cd14_idx}")

                sc.tl.diffmap(mono) # Required for DPT
                sc.tl.dpt(mono)
                
                # Plot Pseudotime & Markers
                fig_dpt, ax = plt.subplots(1, 3, figsize=(15, 5))
                sc.pl.umap(mono, color='dpt_pseudotime', ax=ax[0], show=False, title="Pseudotime")
                
                # Gene Trends (Plot individually to avoid axis errors)
                for i, gene in enumerate(['CD14', 'FCGR3A']):
                    if gene in mono.var_names:
                        sc.pl.umap(mono, color=gene, ax=ax[i+1], show=False, title=gene)
                    elif mono.raw is not None and gene in mono.raw.var_names:
                        sc.pl.umap(mono, color=gene, ax=ax[i+1], show=False, title=gene, use_raw=True)
                    else:
                        ax[i+1].text(0.5, 0.5, f"{gene} not found", ha='center', va='center')
                        ax[i+1].axis('off')
                
                plt.tight_layout()
                save_plot(fig_dpt, output_pseudotime_plot, "Monocyte Differentiation")
                log("  Trajectory analysis complete.")
                
            else:
                log(f"  Not enough CD14+ cells found for trajectory (n={np.sum(is_mono)}).")
                save_plot(plt.figure(), output_paga_plot, f"Insufficient Monocytes (n={np.sum(is_mono)})")
                save_plot(plt.figure(), output_pseudotime_plot, "Insufficient Monocytes")
        else:
            log(f"  CD14 not found in dataset. Top var_names: {list(rna.var_names[:10])}")
            save_plot(plt.figure(), output_paga_plot, "CD14 Missing")
            save_plot(plt.figure(), output_pseudotime_plot, "CD14 Missing")
            
    except Exception as e:
        import traceback
        log(f"Trajectory Inference Failed: {e}")
        traceback.print_exc()
        # Create placeholders so Snakemake doesn't crash on subsequent runs if this fails
        save_plot(plt.figure(), output_paga_plot, f"Error: {str(e)[:20]}")
        save_plot(plt.figure(), output_pseudotime_plot, "Trajectory Error")

else:
    log("Skipping Trajectory (No RNA or Clusters).")
    save_plot(plt.figure(), output_paga_plot, "Skipped")
    save_plot(plt.figure(), output_pseudotime_plot, "Skipped")

log("✓ Done.")
