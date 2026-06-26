rule download_data:
    output:
        config["data"]["rna_matrix"]
    shell:
        "curl -L -o {output} https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"

rule qc_check:
    input:
        rna = config["data"]["rna_matrix"]
    output:
        txt = "results/qc/qc_metrics.txt",
        json = "results/qc/qc_metrics_mqc.json"
    script:
        "../../scripts/qc_check.py"

rule multiqc:
    input:
        txt = "results/qc/qc_metrics.txt",
        json = "results/qc/qc_metrics_mqc.json"
    output:
        "results/qc/multiqc_report.html"
    shell:
        ".venv\\Scripts\\python -m multiqc results/qc/ -o results/qc/ -n multiqc_report.html --force"

