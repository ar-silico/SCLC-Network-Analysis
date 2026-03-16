# Excluded Datasets

## GSE6044
Excluded due to insufficient gene coverage — dataset contained 8,085 
detected genes, compared to 18,362–21,655 genes in the three included 
datasets (GSE30219: 21,655; GSE149507: 20,172; GSE60052: 18,362 after 
harmonization to microarray gene space). Gene coverage of approximately 
40% relative to included datasets was deemed insufficient for reliable 
consensus DEG derivation, as low coverage risks systematic 
underrepresentation of lower-abundance transcripts that are pontentially
critical to SCLC biology.

## GSE40275
Excluded due to data quality issues identified during preprocessing — 
specifically incorrect probe summarization and malformed gene indexing 
that could not be reliably corrected without introducing analytical bias.

## Exclusion Criteria
Datasets were assessed for:
- Sufficient gene detection rate
- Sample quality metrics evaluated through boxplot distribution 
  analysis and Principal Component Analysis (PCA)
- Appropriate experimental design (SCLC vs Normal)
- Probe summarization integrity and gene index validity
- Platform compatibility with downstream analysis


Only datasets meeting all criteria were retained for consensus 
DEG derivation to ensure analytical reliability.