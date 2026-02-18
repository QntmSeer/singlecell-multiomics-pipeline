rule qc_check:
    input:
        rna = config["data"]["rna_matrix"]
    output:
        "results/qc/qc_metrics.txt"
    script:
        "../../scripts/qc_check.py"

rule multiqc:
    input:
        "results/qc/qc_metrics.txt"
    output:
        "results/qc/multiqc_report.html"
    shell:
        "touch {output}"
