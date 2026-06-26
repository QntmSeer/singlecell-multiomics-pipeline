# Single-Cell Multi-Omics Pipeline for Inflammation Analysis

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Snakemake](https://img.shields.io/badge/Snakemake-Workflow-blue.svg)
![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Version](https://img.shields.io/badge/version-v2.0.0--stable--local-green)

## 📽️ Project Overview

This repository contains a reproducible, containerized bioinformatics pipeline for the analysis of **single-cell multi-omics data (scRNA-seq + scATAC-seq)**. The workflow is built using **Snakemake** and utilizes **Scanpy**, **Muon**, and **scvi-tools** for the joint integration of transcriptomic and chromatin accessibility modalities.

### 🔬 Biological Context

**Research Question:** *How does chromatin accessibility landscape correlate with gene expression heterogeneity in human peripheral blood mononuclear cells (PBMCs) under inflammatory conditions?*

This pipeline targets the **10x Genomics Multiome (RNA + ATAC)** dataset to:
1. Identify distinct immune cell subsets (T-cells, B-cells, Monocytes, NK cells).
2. Integrate transcriptomics and chromatin accessibility using **WNN (Weighted Nearest Neighbor)** and **MultiVI (Deep Generative Variational Autoencoders)**.
3. Perform trajectory inference and probabilistic fate mapping on the monocyte differentiation continuum using **CellRank 2** (GPCCA).
4. Highlight regulatory elements potentially driving inflammatory responses in specific cell populations.

## 🛠️ Pipeline Architecture

```
Raw H5 Data Download (Automated)
        │
        ▼
   ┌─────────┐
   │   QC    │  <- Calculate cell QC metrics, export MultiQC JSON
   └────┬────┘
        ├──────────────────────────┐
        ▼                          ▼
 ┌───────────────┐          ┌───────────────┐
 │ RNA Modality  │          │ ATAC Modality │
 │ Normalize/HVG │          │ TF-IDF / LSI  │
 └──────┬────────┘          └──────┬────────┘
        │                          │
        ▼                          ▼
 ┌──────────────────────────────────────────┐
 │           Joint Integration              │
 ├──────────────────────────────────────────┤
 │  • WNN Integration (Heuristic Neighbors) │
 │  • MultiVI Deep Generative Model (VAE)   │
 └──────────────────┬───────────────────────┘
                    │
                    ▼
          ┌───────────────────┐
          │ Trajectory & Fate │
          ├───────────────────┤
          │ • Diffusion Time  │
          │ • PAGA Graph      │
          │ • CellRank 2      │
          └─────────┬─────────┘
                    │
                    ▼
          ┌───────────────────┐
          │     Plotting      │  <- RNA, ATAC, WNN, MultiVI,
          └───────────────────┘    CellRank Trajectories, MultiQC
```

### Advanced Publication Integrations
- **MultiVI (Nature Methods, 2022):** A deep generative model (VAE) that jointly embeds RNA + ATAC modalities to capture non-linear interactions, correct batch effects, and handle data sparsity.
- **CellRank 2 (Nature Methods, 2024):** A unified fate mapping framework using a Markov chain-based GPCCA estimator to calculate absorption/fate probabilities towards specific endpoints.
- **MultiVelo Template (Nature Biotechnology, 2023):** Included template for chromatin-coupled velocity modeling.
- **SCENIC+ Template (Nature Methods, 2023):** Included template for enhancer-driven Gene Regulatory Network (eGRN) mapping.

## 📂 Repository Structure

```
singlecell-multiomics-pipeline/
├── config/             # Configuration files
│   └── config.yaml
├── workflow/           # Snakemake workflow definition
│   ├── Snakefile
│   └── rules/
│       ├── qc.smk      # Automated download and MultiQC
│       └── analysis.smk # Modality processing and trajectory mapping
├── envs/               # Conda environment definitions
│   └── sc-omics.yaml
├── scripts/            # Analysis and plotting scripts
│   ├── analysis.py     # Main pipeline processing and trajectory
│   ├── qc_check.py     # Performs QC and exports MultiQC stats
│   ├── generate_plots.py # Standalone plotting script
│   ├── multivelo_analysis.py # Template for chromatin velocity
│   └── scenicplus_analysis.py # Template for eGRN inference
├── results/            # Output plots and tables (not committed)
├── data/               # Input data directory (not committed)
├── Dockerfile
└── README.md
```

## 🛰️ Usage

### Tested Environment
- Python 3.10 - 3.14
- `scanpy>=1.10`
- `muon` / `mudata`
- `scvi-tools>=1.1`
- `cellrank>=2.0`
- `pytorch` / `jax` / `jaxlib`
- `multiqc`

All dependencies are defined in `envs/sc-omics.yaml`.

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/QntmSeer/singlecell-multiomics-pipeline.git
   cd singlecell-multiomics-pipeline
   ```

2. Create the environment:
   ```bash
   conda env create -f envs/sc-omics.yaml
   conda activate sc-omics
   ```

### Running the Workflow

The Snakemake workflow is fully automated. The HDF5 dataset will download automatically on the first run:

```bash
snakemake --cores 8
```

To re-run specifically the analysis and plotting:
```bash
snakemake --cores 8 --force analysis
```

To run the standalone plotting script directly without Snakemake:
```bash
python scripts/generate_plots.py
```

## 📊 Expected Outputs

| Output | Location | Description |
|--------|----------|-------------|
| QC Metrics Table | `results/qc/qc_metrics_mqc.json` | Custom JSON table for MultiQC reporting |
| MultiQC Report | `results/qc/multiqc_report.html` | Aggregated HTML quality control report |
| RNA UMAP | `results/plots/umap_rna.png` | UMAP based on RNA expression (WNN clusters) |
| ATAC UMAP | `results/plots/umap_atac.png` | UMAP based on Chromatin Accessibility (LSI) |
| WNN UMAP | `results/plots/umap_wnn.png` | Weighted Nearest Neighbors joint integration |
| MultiVI UMAP | `results/plots/umap_multivi.png` | Deep generative joint latent representation (MultiVI) |
| Marker Genes | `results/plots/umap_markers.png` | Immune marker expression profiles (`CD3D`, `CD14`, `MS4A1`, `GNLY`) |
| PAGA Graph | `results/plots/trajectory_paga.png` | Cluster-level connectivity graph for monocytes |
| CellRank 2 Trajectory | `results/plots/cellrank_trajectory.png` | GPCCA macrostates and fate absorption probabilities |
| Pseudotime Timeline | `results/plots/trajectory_pseudotime.png` | Continuous developmental timeline (Classical -> Non-Classical) |

## 🐳 Docker

```bash
docker build -t sc-omics-pipeline .
docker run -v $(pwd)/data:/app/data -v $(pwd)/results:/app/results sc-omics-pipeline snakemake --cores 8
```

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## 🔗 Data Source

[10k Human PBMCs, Multiome v1.0, Chromium X - 10x Genomics](https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-x-1-standard-2-0-0)  
Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
