"""
MultiVelo Analysis Template
===========================
Reference implementation showing how to integrate RNA velocity and chromatin
accessibility to model cell fate dynamics using the MultiVelo model.

This script expects:
  1. A loom file (e.g. from velocyto) containing spliced and unspliced counts.
  2. The peaks count matrix or promoter accessibility scores (from scATAC-seq).
"""
import scvelo as scv
import multivelo as mv
import scanpy as sc
import numpy as np

def run_multivelo_pipeline(rna_loom_path, atac_h5_path, output_plot_path):
    print("1. Loading raw RNA velocity data (spliced/unspliced)...")
    adata_rna = scv.read(rna_loom_path)
    adata_rna.var_names_make_unique()
    
    print("2. Loading ATAC peaks matrix...")
    adata_atac = sc.read_10x_h5(atac_h5_path)
    adata_atac.var_names_make_unique()
    
    print("3. Aligning RNA and ATAC modalities...")
    common_barcodes = np.intersect1d(adata_rna.obs_names, adata_atac.obs_names)
    adata_rna = adata_rna[common_barcodes].copy()
    adata_atac = adata_atac[common_barcodes].copy()
    
    print("4. Preprocessing RNA...")
    scv.pp.filter_and_normalize(adata_rna, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=30)
    
    print("5. Preprocessing ATAC peaks...")
    # MultiVelo aggregates ATAC peaks near promoter regions to model chromatin opening kinetics
    adata_atac = mv.aggregate_peaks_to_genes(adata_atac, promoter_only=True)
    
    print("6. Initializing and running MultiVelo...")
    # Recover dynamics of transcription rates and chromatin opening/closing times
    mv.recover_dynamics(adata_rna, adata_atac)
    
    # Calculate cell-specific velocity vectors
    mv.compute_velocity(adata_rna, adata_atac)
    
    # Project velocity vectors onto joint UMAP
    scv.tl.velocity_graph(adata_rna)
    
    print("7. Visualizing chromatin-coupled velocity stream...")
    scv.pl.velocity_embedding_stream(
        adata_rna, 
        basis='umap', 
        title='MultiVelo: Chromatin-Coupled Splicing Dynamics',
        save=output_plot_path
    )
    print(f"✓ MultiVelo trajectory plot saved to: {output_plot_path}")

if __name__ == "__main__":
    print("MultiVelo Template script loaded. Run this script in a python session with raw data.")
