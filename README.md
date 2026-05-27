# SCLC Network Analysis
## Prognostic biomarkers and parallel oncogenic programs in SCLC concealed by network topology bias uncovered through cluster-resolved analysis
### Reproducible computational pipeline for consensus DEG analysis, PPI network clustering, proteomics validation, and survival analysis in SCLC.

### Authors
Aryan Pyakurel — Lead author, BSc. Bioinformatics, Kathmandu University, Nepal
Saksham Shakya — Co-author, BSc. Bioinformatics, Kathmandu University, Nepal
Dikshant Regmi — Co-author, BSc. Bioinformatics, Kathmandu University, Nepal
Amrita Acharya — Supervisor, Kathmandu University, Nepal

### Description
This repository contains all code and analysis pipelines for a multi-platform bioinformatics study identifying parallel oncogenic programs and biomarkers in Small Cell Lung Cancer (SCLC) through consensus transcriptomic differential gene expression (DEG) network analysis, proteomics validation, and survival analysis.

### Key Findings
- Identified four biologically distinct parallel oncogenic programs in SCLC through MCL clustering of consensus DEG-derived PPI network, with downregulation dominance observed across three of four clusters and the broader consensus DEG landscape
- Cross-platform proteomics validation of network-derived hub genes in 112 treatment-naive SCLC patients
- CASP1, IL18 and TLR2 identified as novel prognostic biomarkers: CASP1 (HR=0.372, 95% CI: 0.207–0.667, adjusted p=0.0044), IL18 (HR=0.382, 95% CI: 0.213–0.685, adjusted p=0.0044), TLR2 (HR=0.490, 95% CI: 0.278–0.863, adjusted p=0.0427) — surviving Benjamini-Hochberg multiple testing correction
- Demonstrated proliferation topological bias in cancer PPI networks and proposed MCL clustering as methodological solution

### Datasets
| Dataset | Platform | Samples | Source |
|---|---|---|---|
| GSE30219 | Microarray | 21 SCLC / 14 Normal | GEO |
| GSE149507 | Microarray | 14 paired | GEO |
| GSE60052 | RNA-seq | 79 SCLC / 7 Normal | GEO |
| GSE43346 | Microarray | 23 SCLC / 1 Normal | GEO |
| GSE69091 | Mouse Notch model | 3 Notch ICD / 3 Control | GEO |
| Liu et al. 2024 | Proteomics | 112 SCLC | Cell |

### Repository Structure
```
SCLC-network-analysis/
├── data/                    # Raw data documentation (raw files not tracked — see DECISIONS.md)
│   ├── raw/                 # GEO SOFT.gz files (not tracked, must be obtained manually)
│   └── processed/           # Processed expression matrices and metadata
├── notebooks/               # Jupyter notebooks (source — authoritative analytical code)
├── results/                 # Pipeline outputs
│   ├── DEG/                 # Differential expression results
│   ├── PPI/                 # Network files (STRING/Cytoscape outputs)
│   └── enrichment/          # DAVID functional enrichment results
├── figures/                 # All generated figures
├── logs/                    # Snakemake execution logs (not tracked)
├── checkpoints/             # Snakemake stage checkpoints (not tracked)
├── Snakefile                # Snakemake pipeline definition
├── requirements.txt         # Python dependencies
├── DECISIONS.md             # Analytical decision log
└── README.md                # This file
```

### Requirements
See requirements.txt for full dependency list.

### Usage
#### Prerequisites
- Python 3.11 (see `requirements.txt` for dependencies)
- R with limma and BioConductor installed (required for GSE69091 DEG analysis via rpy2)
- Snakemake 9.x

#### Setup
1. Clone the repository
2. Install Python dependencies: `pip install -r requirements.txt`
3. Obtain required raw data files and place in the correct directories (see `DECISIONS.md` Section 10 for full list and sources)
4. Complete manual checkpoints (STRING, Cytoscape, DAVID) as documented in `DECISIONS.md` Section 10

#### Running the Pipeline
```bash
# From the repository root
snakemake --cores 1 --dry-run   # Verify pipeline before running
snakemake --cores 1              # Run full pipeline
```

The pipeline will pause at manual checkpoints — follow instructions in `DECISIONS.md` for each checkpoint before continuing.

#### Notebooks
All source notebooks are in `notebooks/`. They can be run interactively in Jupyter or executed via Snakemake. See `DECISIONS.md` for analytical decisions and justifications.

### Citation
Preprint available on bioRxiv (link to be added upon submission). Citation will be updated upon publication.

### Acknowledgments
We thank the Kathmandu University Biotechnology Department. We thank our supervisor Ms. Amrita Acharya for her continued supervision and feedback on the manuscript. We thank Richard Burchmore, University of Glasgow, for suggesting proteomics-level validation. We also thank our peer Nirja Krishna Pandit for his contribution to the project’s preliminary version. We also acknowledge Claude (Anthropic) for providing assistance in coding syntax and the Snakemake pipeline construction. 
Note: All analytical decisions and biological interpretations are solely those of the authors.

### License
MIT License — see LICENSE file for details.
