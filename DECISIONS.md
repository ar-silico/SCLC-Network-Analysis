# DECISIONS.md
# Analytical Decision Log — SCLC Network Analysis
# Every major methodological choice is documented here with justification.
# This file is the authoritative reference for manuscript methods and
# reviewer responses.
# Last updated: draft 2 — incorporates reviewer-anticipation review

---

## 1. DATASET SELECTION

### Included Datasets
| Dataset | Type | Design | Samples | Role |
|---|---|---|---|---|
| GSE30219 | Microarray (GPL570) | Unpaired | 21 SCLC, 14 Normal | Primary DEG |
| GSE60052 | RNA-seq | Unpaired | 79 SCLC, 7 Normal | Primary DEG |
| GSE149507 | Microarray (GPL) | Paired | 14 pairs (28 after QC) | Primary DEG |
| GSE43346 | Microarray (GPL570) | Unpaired | 23 SCLC, 1 Normal | Directional validation only |
| GSE69091 | Microarray (mouse) | Unpaired | 3 Control, 3 Notch-ICD | Notch functional validation |

### Excluded Datasets
- **GSE6044**: Excluded — insufficient gene coverage (8,085 genes vs 18,362–21,655 in included sets). Low-density array would introduce severe downward bias in consensus gene space and disproportionately reduce the number of detectable consensus DEGs.
- **GSE40275**: Excluded — incorrect probe summarization and malformed gene indexing identified during QC. Results not reproducible or interpretable.
- Justification notebooks retained in repository as `GSE6044_not_usuable_low_gene_screen.ipynb` and `GSE40275_quality_issue_not_usuable.ipynb`.

---

## 2. PREPROCESSING

### Sample Curation (GSE30219)
GSE30219 contains multiple lung cancer subtypes (ADK, Squamous, Basaloid, Carcinoid, LCNE). Only SCLC ('Small Cell Lung Carcinoma' in title) and Normal ('Non tumoral lung' in title) samples were retained. All other subtypes excluded as non-SCLC, non-normal comparators — including them would confound the SCLC vs Normal comparison.

### Sample Curation (GSE149507)
Three rounds of filtering applied:
1. Select SCLC (`_ca`) and Normal (`_n`) samples by title string matching.
2. Remove neoadjuvant-treated samples (neoadjuvant = YES) — treatment alters the transcriptome and would confound SCLC vs Normal comparison.
3. Remove samples with unknown neoadjuvant status (NA) — unverifiable treatment status carries the same confounding risk as confirmed treatment. Specifically GSM4504115 (`SCLC_08_ca`) and GSM4504116 (`SCLC_08_n`) removed.
4. Remove outlier pair GSM4504121 (`SCLC_11_ca`) and GSM4504122 (`SCLC_11_n`) — SCLC sample clusters with normal group in PCA despite neoadjuvant = NO. Consistent with either very early stage disease or experimental contamination of tumor sample with normal tissue. Cannot be verified; removed to avoid confounding paired analysis.
- In all cases, both members of a pair were removed to preserve paired test integrity.
- Pair 12 absent from dataset — not a removal artifact, pair 12 was never in the original GEO submission.

### Pairing Verification (GSE149507)
Pairs verified via GEO sample titles — samples interleaved as `_XX_ca` (SCLC) and `_XX_n` (Normal) in sequential numerical order. Splitting by group label preserves pair correspondence by index position.

### Gene Annotation (GSE149507)
Platform (GPL) provides only Entrez IDs without gene symbols. mygene API used to map Entrez IDs to HGNC gene symbols for consistency across datasets.

### log2(x+1) Transformation (GSE43346)
GSE43346 provided as raw counts — not pre-normalized. log2(x+1) transformation applied for consistency with pre-normalized log2 datasets. The +1 pseudocount prevents log(0) for zero-count genes.

### GSE43346 Normal Sample
Single normal sample confirmed as normal lung tissue. Sample classification by 'Lung' title string is unambiguous given only one normal sample exists in this dataset.

### GSE60052 Low-Expression Filter
RNA-seq is prone to technical zeros for lowly expressed genes. Genes expressed (>0) in fewer than 10% of samples were removed. Rationale: genes expressed in ≥10% of samples are unlikely to be universally absent — zeros in <10% of samples are more likely technical artifacts than true absence of expression.

### GSE60052 Microarray Union Filter
RNA-seq matrix filtered to the union of GSE30219 and GSE149507 gene spaces. Rationale: consensus DEGs are bounded by the smallest dataset — harmonizing the RNA-seq set to the microarray gene universe prevents RNA-seq-exclusive genes from inflating the consensus. This ensures all three datasets contribute comparably to the consensus.

### GSE69091 Platform Annotation
This microarray uses `gene_assignment` field with `//` delimiter — gene symbol is the second field. The `///` multi-gene filter applied to other datasets was intentionally skipped here. In this platform, `///` separates transcript isoforms of the same gene, not distinct gene symbols — confirmed by inspection of raw annotation.

### Duplicate Gene Collapse
All datasets: duplicate probe-to-gene mappings collapsed by mean expression. Mean chosen over max to avoid upward bias from outlier probes.

---

## 3. DIFFERENTIAL EXPRESSION ANALYSIS

### Statistical Test Selection
- **GSE30219, GSE60052**: Mann-Whitney U test (two-sided, unpaired). Non-parametric test chosen because normality of gene expression distributions cannot be assumed, particularly with unequal and relatively small group sizes.
- **GSE149507**: Wilcoxon signed-rank test (paired). Paired test used to account for inter-patient variability and increase statistical power. Appropriate because matched SCLC/Normal samples come from the same patients.
- **GSE69091**: limma with empirical Bayes shrinkage (via rpy2). 3v3 sample size too small for reliable non-parametric testing. limma's eBayes shrinkage is specifically designed for small-sample microarray data and borrows information across genes to stabilize variance estimates. eBayes proportion parameter = 0.01 (default value, stated explicitly for transparency).

### Multiple Testing Correction
Benjamini-Hochberg (BH) FDR correction applied to all datasets. BH chosen over Bonferroni because Bonferroni is overly conservative for large-scale genomic data with correlated tests. BH < 0.05 is the standard threshold for genomic analyses.

### DEG Threshold
Strict threshold: padj < 0.05 AND |log2FC| ≥ 1. The FC filter removes genes that are statistically significant but biologically trivial (small magnitude changes). |log2FC| ≥ 1 represents a twofold change, a standard cutoff in microarray literature.

### Stats-Only Version
A stats-only version (padj < 0.05, no FC filter) was also saved for each dataset. This was produced for exploratory purposes and reviewed but the strict version was used for consensus computation.

---

## 4. CONSENSUS DEG COMPUTATION

### Consensus Definition
Consensus DEGs defined as genes meeting the strict threshold (padj < 0.05, |log2FC| ≥ 1) in the same direction across all three primary datasets simultaneously. This is the most conservative approach — only genes with robust, cross-platform, directionally consistent signals are retained.

### Gene Universe
The three-way overlap of all three expression matrices = 15,504 genes. This represents the testable gene space common to all three platforms and was used as the universe for Fisher's exact test.

### GSE43346 Directional Validation
GSE43346 (23 SCLC vs 1 Normal) cannot yield reliable p-values due to single normal sample. Used solely for directional concordance validation. Fold change computed as mean of SCLC samples minus single normal sample (already log2 scale). Result: 99.0% concordance for downregulated genes, 98.4% for upregulated genes — confirms consensus directionality in an independent dataset.

---

## 5. PPI NETWORK AND CLUSTERING

### STRING Confidence Threshold
Confidence threshold = 0.700. This is STRING's recommended "high confidence" threshold, balancing network completeness with interaction reliability.

### Initial Hub Network Analysis
A full consensus network was first built and hub genes identified via CytoHubba (MCC + Betweenness). This revealed upregulation dominance among hub genes — discordant with the downregulation dominance observed in the raw DEG data. This discordance motivated the MCL clustering approach and is itself a key finding: conventional network hub analysis is topologically biased toward proliferation genes, numerically masking biologically significant suppressed programs. The initial hub network figure is retained in the paper as methodological justification showing the analytical decision was data-driven and to illustrate this topological bias explicitly.

### MCL Clustering
MCL clustering run at inflation = 1.5. Inflation = 1.2 produced overly homogeneous clusters; inflation = 3.0 produced 250+ fragments. Inflation = 1.5 gave biologically interpretable, stable clusters and was selected as optimal.

### Cluster Selection (4 of 99)
Top 4 clusters analyzed. Clusters 5+ show the same suppression-dominant pattern as clusters 2, 3, 4 (cluster 5: 96.2% downregulated — confirmed in `direction_check_for_cluster_4_and_5.ipynb`). Study limited to 4 clusters for scope manageability. Cluster 4 specifically retained because ECM/metastasis biology is clinically critical to SCLC. Cluster 5 data retained as supplementary justification for the suppression dominance finding.

### Hub Gene Selection (Clusters 1–3)
Hub genes identified via CytoHubba using MCC and Betweenness centrality scoring. Final hub gene list = union of top 10 by MCC and top 10 by Betweenness per cluster, deduplicated.

### Hub Gene Selection (Cluster 4)
Cluster 4 contains only 37 nodes — below the minimum network density required for reliable CytoHubba scoring. The CytoHubba step was therefore foregone entirely for this cluster. Hub genes selected as top 5 by STRING node degree. Degree centrality used as a direct proxy for topological importance. Top 5 chosen to maintain proportionality with the effective hub gene yield from clusters 1–3 given the cluster's smaller size.

### Unresolved ENSP IDs
One ENSP ID per cluster was not resolved by mygene API query:
- Cluster 1: `9606.ENSP00000481380` → CCNA2
- Cluster 2: `9606.ENSP00000498321` → MYD88
- Cluster 3: `9606.ENSP00000471239` → ZFP36
- All resolved via cross-reference with STRING annotation TSV files. Resolution is coded, documented, and reproducible.

---

## 6. NOTCH FUNCTIONAL VALIDATION

### Mouse-to-Human Ortholog Mapping
GSE69091 uses mouse genes. mygene API used to map mouse gene symbols to human orthologs (scopes='symbol', species='human'). Genes without human orthologs excluded. Full mapping statistics available in the executed notebook. Ortholog mapping rate not reported in the main manuscript — the Notch validation is a supporting functional observation, not a primary finding, and the paper is already sufficiently large in scope. The core contribution (parallel oncogenic programs, downregulation dominance, topological bias of network analysis, biomarker identification) is established independently.

### Rescue Overlap Significance
Fisher's exact test used to confirm Notch rescue and suppression overlaps are non-random:
- Consensus DOWN ∩ Notch UP: OR = 3.874, p ≈ 0 (n=48 genes rescued)
- Consensus UP ∩ Notch DOWN: OR = 2.789, p = 0.017 (n=7 genes suppressed)
- Universe = 15,504 genes (consensus gene pool)

### Rescued Hub Genes
CD44, CCND1, KDR identified as hub genes rescued by Notch activation in initial analysis using whole-network hub list. Results to be updated after dry run — reanalysis now uses full 50 hub genes across all four clusters.

---

## 7. FUNCTIONAL ENRICHMENT

### DAVID Enrichment Categories
Four categories retained: GOTERM_BP_DIRECT, GOTERM_CC_DIRECT, GOTERM_MF_DIRECT, KEGG_PATHWAY. GO Direct terms provide direct annotation without inferred parent terms, giving cleaner and more interpretable results than indirect GO terms. KEGG captures pathway-level context.

### BH Correction for Enrichment
Benjamini-Hochberg correction applied to DAVID enrichment results. BH < 0.05 used as significance threshold — standard for functional enrichment analysis.

### Enrichment Figure vs Saved CSV
Enrichment bar charts show top 15 terms sorted by raw p-value for visual clarity. Saved filtered CSVs apply BH < 0.05 filter — these are the values reported in the manuscript.

---

## 8. PROTEOMICS VALIDATION

### Dataset
Liu et al. 2024 (Cell) CPTAC proteogenomics dataset, Supplementary Table S3A. 112 treatment-naive SCLC patients. Tumor vs NAT (non-adjacent tissue) log2 fold changes used for validation.

### Two-Step Validation Logic
Proteomics validation applied as a two-step filter:
1. **Presence**: Is the protein detected in the CPTAC mass spectrometry dataset at all? Not all transcripts are translated at detectable levels — mass spectrometry capture is experimentally bounded. Hub genes absent from the proteomics dataset are classified as "not detected", not discordant. Absence does not contradict transcriptomic findings.
2. **Directionality**: For detected proteins, does the protein-level fold change (tumor vs NAT) show directional concordance with the transcriptomic dysregulation? Concordance confirms the transcriptional signal propagates to the protein level.

### Hub Gene Attrition (50 → 27)
Of 50 hub genes taken to proteomics validation, 27 were detected and directionally concordant in S3A. The remaining 23 were absent from the mass spectrometry dataset — not detected, not discordant. This distinction is important: absence reflects the experimental sensitivity of the proteomics assay, not a failure of the transcriptomic signal.

### S4 Non-Applicability
Liu et al. Supplementary Table S4 was assessed and documented as non-applicable by design. S4A uses a survival-first concordance-filtered pipeline orthogonal to our network topology approach. S4C/S4E are clinical covariate tables for Liu et al.'s own markers only. Not a limitation — different analytical frameworks answer different questions.

---

## 9. SURVIVAL ANALYSIS

### Patient Cohort
104 patients after QC (missing survival time removed), 56 events.

### Median Split Justification
Hartigan's dip test applied to protein expression distributions for all 27 proteomics-validated hub genes. All genes unimodal (all p > 0.05) — no bimodal distribution warranting alternative splitting strategy. Median split is appropriate and standard for continuous biomarker survival analysis.

### Gene Selection for Survival
Only the 27 S3A-proteomics-validated hub genes were taken to survival analysis. Proteomics confirmation required before survival testing — ensures we test only genes with protein-level evidence in the same patient cohort, and reduces the multiple testing burden from 50 to 27.

### Initial Screen Threshold (p < 0.10)
An initial KM screen was applied across all 27 genes using p < 0.10 as a liberal threshold to identify significant and borderline candidates for BH correction. The p < 0.10 threshold was chosen rather than p < 0.05 because SCLC is a rare disease — available cohorts are inherently underpowered and borderline findings in these datasets may reach significance in larger future cohorts. Furthermore, directional concordance of borderline genes across platforms independently supports the parallel oncogenic programs hypothesis regardless of survival significance threshold. BH correction was then applied to the 11 genes passing this screen.

### Multiple Testing Correction
BH correction applied to the 11 genes with p < 0.10 from initial KM screen. The staged filtering approach is explicitly justified: 50 hub genes → 27 proteomics-validated → 11 at p < 0.10 → BH correction. Genes with p ≥ 0.10 excluded from correction — standard practice to avoid inflating the correction burden with clearly non-significant results. The correction denominator of 11 is appropriate given the prior biological filtering steps.

### Proportional Hazards Assumption
Schoenfeld residuals test applied to all 11 genes. All p > 0.05 — PH assumption met for all genes. Cox PH model is appropriate.

### Survival Results
- BH-significant: CASP1 (padj=0.0044), IL18 (padj=0.0044), TLR2 (padj=0.0427)
- Borderline (p<0.10, did not survive BH): KIF11, CENPF, TLR3, CEBPB, ITGB5, CD44, CDK1, DLGAP5
- All three BH-significant genes are members of Cluster 2 (immune-suppressed group) — confirming that immune suppression carries a stronger independent survival signal than proliferation in SCLC, consistent with the core thesis.

---

## 10. SOFTWARE AND REPRODUCIBILITY

### LLM Assistance
Claude (Anthropic) used for code syntax assistance and Snakemake pipeline construction. To be disclosed in manuscript acknowledgments.

### Proteomics Raw Data
Raw proteomics files (CPTAC_S1_clinical.xlsb, CPTAC_S3_tumor_NAT.xlsx) not pushed to GitHub — files belong to Liu et al. 2024 and are available from the original publication.

### Manual Steps
The following steps require manual execution and cannot be automated in the Snakemake pipeline:
1. **STRING (Checkpoint 1a)**: Upload consensus DEGs, set confidence=0.700, run full network analysis, download hub network edges/nodes/degrees and CytoHubba top 10 MCC and Betweenness CSVs.
2. **STRING (Checkpoint 1b)**: Run MCL clustering (inflation=1.5) on same full network, download cluster TSVs/edges/nodes/annotations/gene lists for clusters 1–5 and consensus annotations.
3. **Cytoscape + CytoHubba (Checkpoint 2)**: Load cluster{N}_edges.tsv for clusters 1–3, run MCC scoring (top 10) and Betweenness scoring (top 10), export CSVs. Cluster 4 handled automatically by pipeline.
4. **DAVID (Checkpoint 3)**: Upload per-cluster gene lists for clusters 1–4, run functional annotation, download enrichment CSVs.

## 11. POST-CLEANUP FLAGS

### Notch Rescue Reanalysis
`gene_rescued_by_notch.ipynb` currently uses the whole-network hub gene list 
from the initial full network analysis. Should be rerun using the full 50 hub 
genes across all four clusters from `hub_gene_directional_verification.ipynb`. 
Update Snakefile inputs and DECISIONS.md rescue results accordingly once rerun 
is complete.

### ENSP Dynamic Lookup for Full Network
`PPI_plot.ipynb` depends on `hub_network_annotations.tsv` which must be 
downloaded manually from STRING (full hub network, confidence=0.700). 
File to be saved as `results/PPI/hub_network_annotations.tsv` before dry run.
