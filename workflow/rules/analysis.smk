rule analysis:
    input:
        rna = config["data"]["rna_matrix"]
    output:
        umap_rna  = "results/plots/umap_rna.png",
        umap_atac = "results/plots/umap_atac.png",
        umap_wnn  = "results/plots/umap_wnn.png",
        umap_multivi = "results/plots/umap_multivi.png",
        cellrank_trajectory = "results/plots/cellrank_trajectory.png",
        umap_markers = "results/plots/umap_markers.png",
        trajectory_paga = "results/plots/trajectory_paga.png",
        trajectory_pseudotime = "results/plots/trajectory_pseudotime.png"
    script:
        "../../scripts/analysis.py"
