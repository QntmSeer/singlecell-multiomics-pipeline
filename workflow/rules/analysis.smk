rule analysis:
    input:
        rna = config["data"]["rna_matrix"],
        atac = config["data"]["atac_fragments"]
    output:
        umap_rna = "results/plots/umap_rna.png",
        umap_atac = "results/plots/umap_atac.png",
        umap_wnn = "results/plots/umap_wnn.png"
    script:
        "../../scripts/analysis.py"
