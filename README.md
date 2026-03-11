# SCLC Network Analysis
## Multi-platform transcriptomic network analysis revealing parallel oncogenic programs and prognostic biomarkers in small cell lung cancer

### Author
Aryan Pyakurel — BSc Bioinformatics, Kathmandu University, Nepal

### Description
This repository contains all code and analysis pipelines for a multi-platform bioinformatics study identifying parallel oncogenic programs and biomarkers in Small Cell Lung Cancer (SCLC) through consensus transcriptomic differential gene expression (DEG) network analysis, proteomics validation, and survival analysis.

### Key Findings
- Identified four biologically distinct parallel and relevant oncogenic programs in SCLC through MCL clustering of consensus DEG-derived PPI network
- Cross-platform proteomics validation of network-derived hub genes in 112 treatment-naive SCLC patients
- CASP1, IL18 and TLR2 identified as novel prognostic biomarkers (CASP1: HR=0.372, adjusted p=0.0044; IL18: HR=0.382, adjusted p=0.0044; TLR2: HR=0.490, adjusted p=0.0427) — surviving Benjamini-Hochberg multiple testing correction
- Demonstrated proliferation topological bias in cancer PPI networks and proposed MCL clustering as methodological solution

### Datasets
| Dataset | Platform | Samples | Source |
|---|---|---|---|
| GSE30219 | Microarray | 21 SCLC / 14 Normal | GEO |
| GSE149507 | Microarray | 14 paired | GEO |
| GSE60052 | RNA-seq | 79 SCLC / 7 Normal | GEO |
| GSE43346 | Microarray | 23 SCLC / 1 Normal | GEO |
| GSE69091 | Mouse Notch model | 3 Notch ICD / 3 SCLC | GEO |
| Liu et al. 2024 | Proteomics | 112 SCLC | Cell |

### Repository Structure
```
SCLC-network-analysis/
├── data/              # GEO accession numbers and data documentation
├── notebooks/         # Exploratory Jupyter notebooks
├── results/           # Key output files
├── figures/           # All generated figures
└── workflow/          # Snakemake pipeline (coming soon)
```

### Requirements
See requirements.txt for full dependency list.

### Usage
Snakemake pipeline coming soon. For now refer to notebooks/ for step by step analysis.

### Citation
Manuscript under review. Citation will be added upon publication.

### Acknowledgments
Code written in Python. Syntax assistance provided by Claude (Anthropic). All analytical decisions and biological interpretations are solely those of the authors.

### License
MIT License — see LICENSE file for details.
