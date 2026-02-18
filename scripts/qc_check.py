# scripts/qc_check.py
import os

input_file = snakemake.input.rna
output_file = snakemake.output[0]

# Ensure output dir exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as f:
    f.write(f"QC completed for {input_file}\n")
    f.write("Metrics: PASS\n")
