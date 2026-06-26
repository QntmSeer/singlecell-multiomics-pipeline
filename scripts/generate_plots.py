"""
Standalone plot generator — run directly with:
    python scripts/generate_plots.py

No Snakemake needed. Reads from data/ and writes to results/plots/
"""
import scanpy as sc
import muon as mu
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import scipy.sparse

matplotlib.use("Agg")
sc.settings.verbosity = 1
np.random.seed(42)

DATA_FILE = "data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
OUT_DIR   = "results/plots"
os.makedirs(OUT_DIR, exist_ok=True)

if __name__ == "__main__":
    # ── Load ──────────────────────────────────────────────────────────────────────
    print("Loading data...")
    try:
        mdata = mu.read_10x_h5(DATA_FILE)
        mdata.var_names_make_unique()
    except Exception as e:
        print(f"Muon read failed ({e}). Falling back to RNA-only mode.")
        mdata = mu.MuData({'rna': sc.read_10x_h5(DATA_FILE)})
        mdata.var_names_make_unique()

    has_atac = 'atac' in mdata.mod
    has_rna  = 'rna'  in mdata.mod
    print(f"Modalities detected: {list(mdata.mod.keys())}")

    # ── QC ────────────────────────────────────────────────────────────────────────
    if has_rna:
        print("QC filtering...")
        rna = mdata['rna']
        rna.var["mt"] = rna.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        keep_cells = (
            (rna.obs.n_genes_by_counts > 200)  &
            (rna.obs.n_genes_by_counts < 6000) &
            (rna.obs.pct_counts_mt    < 20)
        )
        mdata = mdata[keep_cells, :].copy()
        print(f"  -> {mdata.n_obs} cells after QC")

    # ── RNA Normalization & PCA ───────────────────────────────────────────────────
    if has_rna:
        print("Processing RNA modality...")
        rna = mdata['rna']
        rna.var_names_make_unique()
        sc.pp.normalize_total(rna, target_sum=1e4)
        sc.pp.log1p(rna)
        rna.raw = rna
        sc.pp.highly_variable_genes(rna, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pp.pca(rna, svd_solver="arpack")
        sc.pp.neighbors(rna, n_neighbors=15, n_pcs=40)
        sc.tl.umap(rna)
        sc.tl.leiden(rna, resolution=0.5, key_added='leiden_rna')

    # ── ATAC Normalization & LSI ──────────────────────────────────────────────────
    if has_atac:
        print("Processing ATAC modality (LSI)...")
        atac = mdata['atac']
        mu.atac.pp.tfidf(atac, scale_factor=1e4)
        mu.atac.tl.lsi(atac)
        if 'X_lsi' in atac.obsm:
            atac.obsm['X_lsi_clean'] = atac.obsm['X_lsi'][:, 1:]
            sc.pp.neighbors(atac, use_rep='X_lsi_clean', n_neighbors=15, n_pcs=39)
        else:
            sc.pp.neighbors(atac, use_rep='X_lsi', n_neighbors=15, n_pcs=40)
        sc.tl.umap(atac)
        sc.tl.leiden(atac, resolution=0.5, key_added='leiden_atac')

    # ── WNN Integration ───────────────────────────────────────────────────────────
    final_clusters = 'leiden_rna'
    if has_rna and has_atac:
        print("Running WNN integration...")
        try:
            mu.pp.neighbors(mdata)
            mu.tl.umap(mdata)
            mu.tl.leiden(mdata, key_added='leiden_wnn', resolution=0.5)
            if 'leiden_wnn' in mdata.obs:
                mdata['rna'].obs['leiden_wnn']  = mdata.obs['leiden_wnn']
                mdata['atac'].obs['leiden_wnn'] = mdata.obs['leiden_wnn']
            final_clusters = 'leiden_wnn'
        except Exception as e:
            print(f"WNN failed: {e}")

    # ── MultiVI Integration ───────────────────────────────────────────────────────
    has_scvi = False
    try:
        import scvi
        has_scvi = True
    except ImportError:
        print("scvi-tools is not installed. Skipping MultiVI.")

    if has_scvi and has_rna and has_atac:
        print("Running MultiVI integration...")
        try:
            scvi.model.MULTIVI.setup_mudata(mdata, rna_layer=None, atac_layer=None, modalities={"rna_layer": "rna", "atac_layer": "atac"})
            import torch
            torch.set_num_threads(min(8, os.cpu_count() or 4))
            model = scvi.model.MULTIVI(mdata)
            model.train(max_epochs=2, batch_size=256)
            mdata.obsm["X_multivi"] = model.get_latent_representation()
            sc.pp.neighbors(mdata, use_rep="X_multivi", key_added="multivi")
            sc.tl.umap(mdata, neighbors_key="multivi")
            sc.tl.leiden(mdata, resolution=0.5, key_added="leiden_multivi", neighbors_key="multivi")
            if 'leiden_multivi' in mdata.obs:
                mdata['rna'].obs['leiden_multivi'] = mdata.obs['leiden_multivi']
                mdata['atac'].obs['leiden_multivi'] = mdata.obs['leiden_multivi']
                final_clusters = 'leiden_multivi'
        except Exception as e:
            print(f"MultiVI failed: {e}")

    # ── Helper: Save Figure ───────────────────────────────────────────────────────
    def save_plot(fig, filename, title):
        path = os.path.join(OUT_DIR, filename)
        if isinstance(fig, matplotlib.figure.Figure):
            fig.suptitle(title, y=1.02, fontsize=14, fontweight='bold')
            fig.savefig(path, dpi=150, bbox_inches='tight')
            plt.close(fig)
        else:
            plt.suptitle(title, y=1.02, fontsize=14, fontweight='bold')
            plt.savefig(path, dpi=150, bbox_inches='tight')
            plt.close()
        print(f"  [OK] {path}")

    # ── Plotting ──────────────────────────────────────────────────────────────────
    print("Saving plots...")

    # RNA UMAP
    if has_rna:
        fig = sc.pl.umap(mdata['rna'], color=final_clusters, title="RNA UMAP", return_fig=True, show=False)
        save_plot(fig, "umap_rna.png", "RNA Modality")

    # ATAC UMAP
    if has_atac:
        fig = sc.pl.umap(mdata['atac'], color=final_clusters, title="ATAC UMAP", return_fig=True, show=False)
        save_plot(fig, "umap_atac.png", "ATAC Modality")

    # WNN UMAP
    if has_rna and has_atac and 'X_umap' in mdata.obsm:
        mu.pl.umap(mdata, color='leiden_wnn', title="WNN UMAP", show=False)
        save_plot(plt.gcf(), "umap_wnn.png", "WNN Integrated")

    # MultiVI UMAP
    if 'X_multivi' in mdata.obsm:
        sc.pl.umap(mdata, color='leiden_multivi', title="MultiVI UMAP", show=False)
        save_plot(plt.gcf(), "umap_multivi.png", "MultiVI Integrated")

    # Markers
    markers = {"CD3D": "T cells", "CD14": "Monocytes", "MS4A1": "B cells", "GNLY": "NK cells"}
    if has_rna:
        use_basis = 'X_umap'
        if has_rna and has_atac and 'X_umap' in mdata.obsm:
            mdata['rna'].obsm['X_umap_wnn'] = mdata.obsm['X_umap']
            use_basis = 'X_umap_wnn'
        present = {k: v for k, v in markers.items() if k in rna.var_names}
        if present:
            fig, axes = plt.subplots(2, 2, figsize=(10, 8))
            axes = axes.flatten()
            for i, (gene, label) in enumerate(present.items()):
                sc.pl.embedding(mdata['rna'], basis=use_basis, color=gene, use_raw=True, color_map='Reds', ax=axes[i], show=False, frameon=False, title=f"{gene} ({label})")
            for j in range(len(present), 4):
                axes[j].axis('off')
            save_plot(fig, "umap_markers.png", "Immune Markers")

    # Trajectory & CellRank
    if has_rna and final_clusters in mdata.obs:
        try:
            def get_gene_expr(adata, gene_name):
                if gene_name not in adata.var_names:
                    return np.zeros(adata.n_obs)
                expr = adata[:, gene_name].X
                if scipy.sparse.issparse(expr):
                    return expr.toarray().flatten()
                return np.asarray(expr).flatten()

            cd14_expr = get_gene_expr(rna, 'CD14')
            fcgr3a_expr = get_gene_expr(rna, 'FCGR3A')
            is_mono = (cd14_expr > 0.5) | (fcgr3a_expr > 0.5)

            if np.sum(is_mono) > 50:
                mono = rna[is_mono, :].copy()
                sc.pp.highly_variable_genes(mono, n_top_genes=2000)
                sc.pp.pca(mono)
                sc.pp.neighbors(mono, n_neighbors=20)
                sc.tl.umap(mono)
                sc.tl.leiden(mono, resolution=0.5, key_added='leiden_mono')

                # PAGA
                sc.tl.paga(mono, groups='leiden_mono')
                sc.pl.paga(mono, color=['leiden_mono', 'CD14', 'FCGR3A'], show=False)
                save_plot(plt.gcf(), "trajectory_paga.png", "Monocyte PAGA Trajectory")

                # Pseudotime
                cd14_idx = np.argmax(get_gene_expr(mono, 'CD14'))
                mono.uns['iroot'] = cd14_idx
                sc.tl.diffmap(mono)
                sc.tl.dpt(mono)

                # CellRank
                has_cellrank = False
                try:
                    import cellrank as cr
                    has_cellrank = True
                except ImportError:
                    pass
                if has_cellrank:
                    pk = cr.kernels.PseudotimeKernel(mono, time_key="dpt_pseudotime")
                    pk.compute_transition_matrix()
                    estimator = cr.estimators.GPCCA(pk)
                    estimator.compute_macrostates(n_states=2, cluster_key="leiden_mono")
                    estimator.predict_terminal_states()
                    estimator.compute_fate_probabilities()
                    fig_cr, axes_cr = plt.subplots(1, 2, figsize=(12, 5))
                    if len(estimator.fate_probabilities.names) >= 3:
                        cr.pl.circular_projection(mono, keys="leiden_mono", ax=axes_cr[0], show=False)
                        axes_cr[0].set_title("Transition Connectivity Graph")
                    else:
                        estimator.plot_macrostates(which="all", ax=axes_cr[0], show=False)
                        axes_cr[0].set_title("Predicted Macrostates")
                    # Custom drawing for fate probabilities on axes_cr[1] (CellRank 2 same_plot ignores ax)
                    def plot_fate_on_ax(adata, data, ax, basis="umap"):
                        coords = adata.obsm[f"X_{basis}"]
                        s = (120_000 / adata.n_obs + 20) / 2
                        ax.scatter(coords[:, 0], coords[:, 1], c="lightgrey", s=s, alpha=0.2, marker=".", edgecolors="none")
                        vals = np.clip(data.X, 0, None)
                        names = list(data.names)
                        lin_colors = list(data.colors)
                        n_lineages = len(names)
                        sorted_idx = np.argsort(vals, axis=1)[:, ::-1][:, :2]
                        for id0 in range(n_lineages):
                            for id1 in range(id0 + 1, n_lineages):
                                cell_mask = np.array([(id0 in row and id1 in row) for row in sorted_idx], dtype=bool)
                                if np.sum(cell_mask) < 2:
                                    continue
                                c_vals = vals[cell_mask, id1] - vals[cell_mask, id0]
                                c_coords = coords[cell_mask]
                                from matplotlib.colors import to_rgba, LinearSegmentedColormap
                                rgba_a = np.array(to_rgba(lin_colors[id0]))
                                rgba_b = np.array(to_rgba(lin_colors[id1]))
                                cmap_colors = [
                                    (*rgba_a[:3], 1.0),
                                    (1.0, 1.0, 1.0, 0.0),
                                    (*rgba_b[:3], 1.0),
                                ]
                                cmap = LinearSegmentedColormap.from_list(f"_cr_{id0}_{id1}", cmap_colors, N=256)
                                abs_max = np.max(np.abs(c_vals)) if np.max(np.abs(c_vals)) > 0 else 1.0
                                c_normed = (c_vals / abs_max + 1) / 2
                                order = np.argsort(np.abs(c_vals))
                                ax.scatter(c_coords[order, 0], c_coords[order, 1], c=c_normed[order], cmap=cmap, vmin=0, vmax=1, s=s, marker=".", edgecolors="none")
                        handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c, label=n, markersize=8) for n, c in zip(names, lin_colors)]
                        ax.legend(handles=handles, frameon=False, loc="center left", bbox_to_anchor=(1.04, 0.5))
                        ax.set_xticks([])
                        ax.set_yticks([])
                        ax.set_xlabel("UMAP1")
                        ax.set_ylabel("UMAP2")

                    plot_fate_on_ax(mono, estimator.fate_probabilities, axes_cr[1])
                    axes_cr[1].set_title("Absorption Probabilities")
                    plt.tight_layout()
                    save_plot(fig_cr, "cellrank_trajectory.png", "CellRank 2 Unified Fate Mapping")

                # Timeline plot
                fig_dpt, ax = plt.subplots(1, 3, figsize=(15, 5))
                sc.pl.umap(mono, color='dpt_pseudotime', ax=ax[0], show=False, title="Pseudotime")
                for i, gene in enumerate(['CD14', 'FCGR3A']):
                    if gene in mono.var_names:
                        sc.pl.umap(mono, color=gene, ax=ax[i+1], show=False, title=gene)
                plt.tight_layout()
                save_plot(fig_dpt, "trajectory_pseudotime.png", "Monocyte Differentiation Timeline")
        except Exception as e:
            print(f"Trajectory plotting failed: {e}")

    print("\n[OK] Standalone plotting complete.")
