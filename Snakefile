# =============================================================================
# Project: Multi-platform SCLC biomarker identification via consensus DEG
#          → PPI network → functional enrichment → proteomics validation
#          → survival analysis
#
# Execution order:
#   Stage 1: Preprocessing (download + QC + normalization)
#   Stage 2: DEG
#   Stage 3: Consensus DEGs + Notch validation + GSE43346 directional check
#   --- MANUAL CHECKPOINT 1a: STRING full network + CytoHubba ---
#   Stage 4a: PPI_plot (full network hub visualization)
#   --- MANUAL CHECKPOINT 1b: STRING MCL clustering ---
#   --- MANUAL CHECKPOINT 2: Cytoscape + CytoHubba per cluster ---
#   Stage 4b: PPI cluster notebooks (1-4)
#   Stage 5: Hub gene verification + Notch rescue + Overlap + Direction check
#   --- MANUAL CHECKPOINT 3: DAVID functional enrichment ---
#   Stage 6: Enrichment (full table + per cluster visualization)
#   Stage 7: Proteomics validation + Survival + PH validation
#
# NOTE: GSE40275 and GSE6044 excluded from pipeline.
# Justification documented in data/EXCLUSIONS.md
#
# To run: snakemake --cores 1
# Dry run: snakemake --cores 1 --dry-run
# =============================================================================


# Final targets — tells Snakemake what the complete pipeline produces
rule all:
    input:
        # Stage 1
        "checkpoints/stage1_preprocessing.done",
        # Stage 2
        "checkpoints/stage2_deg.done",
        # Stage 3
        "checkpoints/stage3_consensus.done",
        # Stage 4
        "figures/hub_gene_network.png",
        "figures/cluster1_hub_network.png",
        "figures/cluster2_hub_network.png",
        "figures/cluster3_hub_network.png",
        "figures/cluster3_hub_network_direction.png",
        "figures/cluster4_network.png",
        "checkpoints/stage4_ppi.done",
        # Stage 5
        "results/DEG/hub_gene_directional_verification.csv",
        "results/DEG/notch_rescued_hub_genes.csv",
        "figures/gene_overlap_venn.png",
        "results/PPI/cluster4_5_direction_summary.csv",
        # Stage 6
        "results/enrichment/cluster1_enrichment_filtered.csv",
        "results/enrichment/cluster2_enrichment_filtered.csv",
        "results/enrichment/cluster3_enrichment_filtered.csv",
        "results/enrichment/cluster4_enrichment_filtered.csv",
        "figures/cluster1_enrichment.png",
        "figures/cluster2_enrichment.png",
        "figures/cluster3_enrichment.png",
        "figures/cluster4_enrichment.png",
        # Stage 7
        "results/proteomics/km_survival_hr_results.csv",
        "results/proteomics/bh_correction_results.csv",
        "results/proteomics/schoenfeld_ph_test_results.csv",
        "figures/km_plots_significant.png"


# =============================================================================
# STAGE 1: PREPROCESSING
# Each notebook: downloads raw data via GEOparse, QC, normalization,
# gene collapse, NA removal, multi-gene probe removal, expression matrix export
# =============================================================================

# --- MAIN DEG DATASETS ---
# These three datasets feed into consensus DEG computation.
# GSE60052 depends on GSE30219 and GSE149507 being processed first.

rule preprocess_GSE30219:
    # NOTE: GSE30219 contains multiple lung cancer subtypes.
    # Samples curated to SCLC and Normal only via title string matching
    # ('Small Cell Lung Carcinoma', 'Non tumoral lung').
    # All other subtypes (ADK, Squamous, Basaloid, Carcinoid, LCNE) excluded
    # as non-SCLC, non-normal comparators. See DECISIONS.md.
    output:
        matrix   = "data/processed/GSE30219_expression_matrix.csv",
        metadata = "data/metadata/GSE30219_metadata.csv",
        boxplot  = "figures/GSE30219_boxplot.png",
        pca      = "figures/GSE30219_PCA.png"
    log:
        "logs/preprocess_GSE30219.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/GSE30219.ipynb \
            --output notebooks/GSE30219_executed.ipynb \
            > {log} 2>&1
        """

rule preprocess_GSE149507:
    # NOTE: Paired microarray dataset (14 pairs, 28 samples after QC).
    # Three-round filtering applied:
    #   1. Select SCLC (_ca) and Normal (_n) samples by title string
    #   2. Remove neoadjuvant-treated samples (neoadjuvant = YES)
    #   3. Remove samples with unknown neoadjuvant status (NA) —
    #      unverifiable treatment status carries same confounding risk
    #      as confirmed treatment (GSM4504115, GSM4504116 removed)
    #   Also removed: GSM4504121/GSM4504122 — SCLC clusters with normal
    #      in PCA, neoadjuvant=NO but unverifiable cause (early stage or
    #      experimental contamination). Whole pair removed to preserve
    #      paired test integrity.
    # Platform (GPL) lacks gene symbols — mygene API used to map
    # Entrez IDs to HGNC gene symbols for consistency. See DECISIONS.md.
    # Final outputs are _final suffixed files.
    output:
        matrix   = "data/processed/GSE149507_filtered_expression_matrix_final.csv",
        metadata = "data/processed/GSE149507_filtered_metadata_final.csv",
        boxplot  = "figures/GSE149507_filtered_boxplot_final.png",
        pca      = "figures/GSE149507_filtered_PCA_final.png"
    log:
        "logs/preprocess_GSE149507.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/GSE149507.ipynb \
            --output notebooks/GSE149507_executed.ipynb \
            > {log} 2>&1
        """

rule preprocess_GSE60052:
    # NOTE: RNA-seq dataset (79 SCLC, 7 Normal).
    # Loads expression from supplementary TSV file (not SOFT file) —
    # pre-normalized log2 values provided by authors.
    # Two-step filtering:
    #   1. Low-expression filter: keep genes expressed (>0) in >=10% of
    #      samples — RNA-seq zeros in <10% of samples likely technical
    #      artifacts rather than true absence of expression.
    #   2. Microarray union filter: RNA-seq matrix filtered to union of
    #      GSE30219 and GSE149507 gene spaces — consensus is bounded by
    #      the smallest dataset; harmonizing prevents RNA-seq-exclusive
    #      genes from inflating the consensus. See DECISIONS.md.
    # DEPENDS ON GSE30219 and GSE149507 preprocessing completing first.
    input:
        gse30219_matrix  = "data/processed/GSE30219_expression_matrix.csv",
        gse149507_matrix = "data/processed/GSE149507_filtered_expression_matrix_final.csv"
    output:
        matrix   = "data/processed/GSE60052_expression_matrix.csv",
        metadata = "data/metadata/GSE60052_metadata.csv",
        boxplot  = "figures/GSE60052_boxplot_coding.png",
        pca      = "figures/GSE60052_PCA_coding.png"
    log:
        "logs/preprocess_GSE60052.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/GSE60052.ipynb \
            --output notebooks/GSE60052_executed.ipynb \
            > {log} 2>&1
        """

# --- INDEPENDENT DIRECTIONAL VALIDATION ---
# GSE43346: 23 SCLC vs 1 Normal. Insufficient for statistical testing.
# Used solely for directional concordance validation of consensus DEGs.
# Fold change computed against mean of cancer samples vs single normal.
# Not included in consensus DEG computation. See DECISIONS.md.

rule preprocess_GSE43346:
    # NOTE: Data provided as raw counts — log2(x+1) applied for
    # consistency with pre-normalized log2 datasets. Single normal
    # sample confirmed as normal lung tissue. See DECISIONS.md.
    output:
        matrix   = "data/processed/GSE43346_expression_matrix.csv",
        metadata = "data/metadata/GSE43346_metadata.csv",
        boxplot  = "figures/GSE43346_boxplot.png",
        pca      = "figures/GSE43346_PCA.png"
    log:
        "logs/preprocess_GSE43346.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/GSE43346.ipynb \
            --output notebooks/GSE43346_executed.ipynb \
            > {log} 2>&1
        """

# --- NOTCH FUNCTIONAL VALIDATION (mouse model) ---
# GSE69091: Mouse SCLC cell line, Notch1-ICD vs Empty vector (3v3).
# Used to validate functional relevance of downregulated consensus DEGs.
# Notch activation rescues suppressed genes — supports biological
# interpretation of transcriptional suppression as oncogenic driver.
# DEG computed via limma (rpy2) due to small sample size (3v3).
# Not included in consensus DEG computation. See DECISIONS.md.

rule preprocess_GSE69091:
    # NOTE: Different platform annotation format — gene symbol extracted
    # from gene_assignment field (// delimiter, second field).
    # /// filter intentionally skipped: in this platform /// separates
    # transcript isoforms of the same gene, not distinct gene symbols.
    # Confirmed by inspection of raw annotation. See DECISIONS.md.
    output:
        matrix   = "data/processed/GSE69091_expression_matrix.csv",
        metadata = "data/metadata/GSE69091_metadata.csv",
        boxplot  = "figures/GSE69091_boxplot.png",
        pca      = "figures/GSE69091_PCA.png"
    log:
        "logs/preprocess_GSE69091.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/GSE69091.ipynb \
            --output notebooks/GSE69091_executed.ipynb \
            > {log} 2>&1
        """

# Checkpoint: all preprocessing complete
rule stage1_complete:
    input:
        "data/processed/GSE30219_expression_matrix.csv",
        "data/processed/GSE149507_filtered_expression_matrix_final.csv",
        "data/processed/GSE60052_expression_matrix.csv",
        "data/processed/GSE43346_expression_matrix.csv",
        "data/processed/GSE69091_expression_matrix.csv"
    output:
        touch("checkpoints/stage1_preprocessing.done")


# =============================================================================
# STAGE 2: DEG
# Mann-Whitney U (unpaired) for GSE30219 and GSE60052.
# Wilcoxon signed-rank (paired) for GSE149507.
# limma via rpy2 for GSE69091 (3v3, too small for non-parametric).
# All use BH correction. Threshold: padj < 0.05, |log2FC| >= 1.
# Stats-only version (padj < 0.05, no FC filter) saved separately
# for use in consensus DEG computation.
# GSE43346 directional validation handled in Consensus_DEGs notebook.
# =============================================================================

rule deg_GSE30219:
    # NOTE: Mann-Whitney U test (two-sided), unpaired.
    # Non-parametric — normality of expression distributions not assumed,
    # especially with unequal group sizes (21 SCLC vs 14 Normal).
    # BH correction applied. Threshold: padj < 0.05, |log2FC| >= 1.
    # Stats-only version saved for consensus DEG computation.
    # See DECISIONS.md.
    input:
        matrix   = "data/processed/GSE30219_expression_matrix.csv",
        metadata = "data/metadata/GSE30219_metadata.csv"
    output:
        full    = "results/DEG/GSE30219_DEG_results.csv",
        stats   = "results/DEG/GSE30219_DEG_results_stats_only.csv",
        strict  = "results/DEG/GSE30219_DEG_results_strict.csv",
        volcano = "figures/GSE30219_volcano.png"
    log:
        "logs/deg_GSE30219.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/DEG_GSE30219.ipynb \
            --output notebooks/DEG_GSE30219_executed.ipynb \
            > {log} 2>&1
        """

rule deg_GSE60052:
    # NOTE: Mann-Whitney U test (two-sided), unpaired.
    # Same methodology as GSE30219. BH correction applied.
    # Threshold: padj < 0.05, |log2FC| >= 1.
    # Stats-only version saved for consensus. See DECISIONS.md.
    input:
        matrix   = "data/processed/GSE60052_expression_matrix.csv",
        metadata = "data/metadata/GSE60052_metadata.csv"
    output:
        full    = "results/DEG/GSE60052_DEG_results.csv",
        stats   = "results/DEG/GSE60052_DEG_results_stats_only.csv",
        strict  = "results/DEG/GSE60052_DEG_results_strict.csv",
        volcano = "figures/GSE60052_volcano.png"
    log:
        "logs/deg_GSE60052.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/DEG_GSE60052.ipynb \
            --output notebooks/DEG_GSE60052_executed.ipynb \
            > {log} 2>&1
        """

rule deg_GSE149507:
    # NOTE: Wilcoxon signed-rank test (paired).
    # GSE149507 is a paired dataset — matched SCLC/Normal from same patients.
    # Paired test used to account for inter-patient variability and increase
    # statistical power. Sample order in metadata preserves pairing —
    # verified via GEO title indexing (_XX_ca / _XX_n interleaved order).
    # BH correction applied. Threshold: padj < 0.05, |log2FC| >= 1.
    # See DECISIONS.md.
    input:
        matrix   = "data/processed/GSE149507_filtered_expression_matrix_final.csv",
        metadata = "data/processed/GSE149507_filtered_metadata_final.csv"
    output:
        full    = "results/DEG/GSE149507_DEG_results.csv",
        stats   = "results/DEG/GSE149507_DEG_results_stats_only.csv",
        strict  = "results/DEG/GSE149507_DEG_results_strict.csv",
        volcano = "figures/GSE149507_volcano.png"
    log:
        "logs/deg_GSE149507.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/DEG_GSE149507.ipynb \
            --output notebooks/DEG_GSE149507_executed.ipynb \
            > {log} 2>&1
        """

rule deg_GSE69091:
    # NOTE: limma used instead of Mann-Whitney — 3v3 sample size
    # too small for non-parametric testing. limma empirical Bayes
    # shrinkage specifically designed for small sample microarray data.
    # eBayes proportion=0.01 is the default value, stated explicitly
    # for transparency. Group order verified against metadata:
    # GSM1692708-710 = Control, GSM1692711-713 = Notch_ICD.
    # See DECISIONS.md.
    input:
        matrix = "data/processed/GSE69091_expression_matrix.csv"
    output:
        full    = "results/DEG/GSE69091_limma_results.csv",
        strict  = "results/DEG/GSE69091_DEG_results_strict.csv",
        volcano = "figures/GSE69091_volcano.png"
    log:
        "logs/deg_GSE69091.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/R_DEG_GSE69091.ipynb \
            --output notebooks/R_DEG_GSE69091_executed.ipynb \
            > {log} 2>&1
        """

# Checkpoint: all DEG complete
rule stage2_complete:
    input:
        "results/DEG/GSE30219_DEG_results.csv",
        "results/DEG/GSE60052_DEG_results.csv",
        "results/DEG/GSE149507_DEG_results.csv",
        "results/DEG/GSE69091_limma_results.csv"
    output:
        touch("checkpoints/stage2_deg.done")


# =============================================================================
# STAGE 3: CONSENSUS DEGs
# Consensus = intersection of strict DEGs (padj < 0.05, |log2FC| >= 1)
# in same direction across all three primary datasets.
# Also handles:
#   - Notch rescue validation (GSE69091 mouse-to-human ortholog mapping)
#   - Fisher's exact test for rescue/suppression overlap significance
#   - GSE43346 directional concordance validation
# Gene universe = 15,504 genes present across all three expression matrices.
# =============================================================================

rule consensus_degs:
    # NOTE: Consensus defined as genes significant (padj < 0.05,
    # |log2FC| >= 1) in same direction across ALL THREE datasets.
    # Most conservative approach — ensures cross-platform robustness.
    # Strict DEGs used (not stats_only) — requires both statistical
    # significance AND biological relevance in all three platforms.
    # Mouse-to-human ortholog mapping via mygene for Notch validation.
    # Fisher's exact test confirms rescue/suppression overlaps are
    # non-random (OR=3.874 p~0 for rescue; OR=2.789 p=0.017 suppression).
    # Gene universe = 15,504 (consensus gene pool across all platforms).
    # WARNING: Fisher's exact test uses hardcoded numbers — must be
    # replaced with dynamic variables during notebook cleanup.
    # WARNING: concordance_pct in validation summary is hardcoded —
    # must be computed dynamically during cleanup. See DECISIONS.md.
    input:
        deg_30219  = "results/DEG/GSE30219_DEG_results_strict.csv",
        deg_149507 = "results/DEG/GSE149507_DEG_results_strict.csv",
        deg_60052  = "results/DEG/GSE60052_DEG_results_strict.csv",
        deg_69091  = "results/DEG/GSE69091_DEG_results_strict.csv",
        gse43346   = "data/processed/GSE43346_expression_matrix.csv",
        meta_43346 = "data/metadata/GSE43346_metadata.csv"
    output:
        consensus  = "results/DEG/consensus_DEGs.csv",
        rescued    = "results/DEG/notch_rescued_consensus_down.csv",
        suppressed = "results/DEG/notch_suppressed_consensus_up.csv",
        validation = "results/DEG/GSE43346_directional_validation.csv",
        venn       = "figures/consensus_DEG_venn.png"
    log:
        "logs/consensus_degs.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/Consensus_DEGs.ipynb \
            --output notebooks/Consensus_DEGs_executed.ipynb \
            > {log} 2>&1
        """

# Checkpoint: consensus DEGs complete
rule stage3_complete:
    input:
        "results/DEG/consensus_DEGs.csv",
        "results/DEG/notch_rescued_consensus_down.csv",
        "results/DEG/GSE43346_directional_validation.csv"
    output:
        touch("checkpoints/stage3_consensus.done")


# =============================================================================
# MANUAL CHECKPOINT 1a: STRING full network + CytoHubba
# -----------------------------------------------------------------------------
# 1. Upload results/enrichment/consensus_genes_for_STRING.txt to STRING
# 2. Set confidence threshold = 0.700
# 3. Download full network: edges TSV, nodes TSV, degrees TSV
# 4. Run CytoHubba on full network (MCC + Betweenness)
# 5. Export top 10 MCC and top 10 Betweenness as CSV
# 6. Place files in results/PPI/:
#    hub_network_edges.tsv, hub_network_degrees.tsv
#    MCC_top_10.csv, Betweenness_top_10.csv
# After completing the above, run Stage 4a (ppi_plot rule).
# =============================================================================


# =============================================================================
# STAGE 4a: PPI PLOT (full network hub visualization)
# Produced before MCL clustering — discordance motivated cluster approach.
# =============================================================================

rule ppi_plot:
    # NOTE: Full consensus network hub gene visualization.
    # Produced BEFORE MCL clustering — discordance between
    # downregulation-dominant DEGs and upregulation-dominant
    # hub network motivated MCL cluster-based approach.
    # This figure retained as methodological justification showing
    # analytical decision was data-driven. See DECISIONS.md.
    input:
        annotations = "results/PPI/hub_network_annotations.tsv",
        mcc         = "results/PPI/MCC_top_10.csv",
        betweenness = "results/PPI/Betweenness_top_10.csv",
        edges       = "results/PPI/hub_network_edges.tsv",
        degrees     = "results/PPI/hub_network_degrees.tsv"
    output:
        network = "figures/hub_gene_network.png"
    log:
        "logs/ppi_plot.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/PPI_plot.ipynb \
            --output notebooks/PPI_plot_executed.ipynb \
            > {log} 2>&1
        """


# =============================================================================
# MANUAL CHECKPOINT 1b: STRING MCL clustering
# -----------------------------------------------------------------------------
# On the same full network from Checkpoint 1a:
# 1. Run MCL clustering (inflation = 1.5)
# 2. Download for each cluster (1-4): edges TSV, nodes TSV,
#    annotations TSV, gene list TXT
# 3. Place files in results/PPI/:
#    cluster{N}_edges.tsv, cluster{N}_annotations.tsv
#    cluster{N}_degrees.tsv, cluster{N}_genes.txt
#    cluster5_annotations.tsv (for direction check)
# Also download: consensus_annotations.tsv (full network annotations)
# After completing the above, proceed to Checkpoint 2.
# =============================================================================


# =============================================================================
# MANUAL CHECKPOINT 2: Cytoscape + CytoHubba (clusters 1, 2, 3 only)
# -----------------------------------------------------------------------------
# For clusters 1, 2, 3 (cluster 4 handled automatically by ppi_cluster4):
# 1. Load cluster{N}_edges.tsv into Cytoscape
# 2. Run CytoHubba with MCC scoring → export top 10 as CSV
# 3. Run CytoHubba with Betweenness scoring → export top 10 as CSV
# 4. Place files in results/PPI/:
#    cluster{N}_MCC_top10.csv, cluster{N}_betweenness_top10.csv
# After completing the above, run Stage 4b rules.
# =============================================================================


# =============================================================================
# STAGE 4b: PPI CLUSTER NOTEBOOKS
# =============================================================================

rule ppi_cluster1:
    # NOTE: Hub genes identified via MCC + Betweenness from CytoHubba.
    # STRING annotation cross-reference. Currently hardcoded as
    # Hub genes identified via annotation TSV lookup
    # See DECISIONS.md.
    input:
        annotations = "results/PPI/cluster1_annotations.tsv",
        mcc         = "results/PPI/cluster1_MCC_top10.csv",
        betweenness = "results/PPI/cluster1_betweenness_top10.csv",
        edges       = "results/PPI/cluster1_edges.tsv",
        degrees     = "results/PPI/cluster1_degrees.tsv",
        consensus   = "results/DEG/consensus_DEGs.csv"
    output:
        network = "figures/cluster1_hub_network.png"
    log:
        "logs/ppi_cluster1.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/PPI_cluster_1.ipynb \
            --output notebooks/PPI_cluster_1_executed.ipynb \
            > {log} 2>&1
        """

rule ppi_cluster2:
    # NOTE: Hub genes identified via MCC + Betweenness from CytoHubba.
    # STRING annotation cross-reference. Currently hardcoded as
    # See DECISIONS.md.
    input:
        annotations = "results/PPI/cluster2_annotations.tsv",
        mcc         = "results/PPI/cluster2_MCC_top10.csv",
        betweenness = "results/PPI/cluster2_betweenness_top10.csv",
        edges       = "results/PPI/cluster2_edges.tsv",
        degrees     = "results/PPI/cluster2_degrees.tsv",
        consensus   = "results/DEG/consensus_DEGs.csv"
    output:
        network = "figures/cluster2_hub_network.png"
    log:
        "logs/ppi_cluster2.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/PPI_cluster_2.ipynb \
            --output notebooks/PPI_cluster_2_executed.ipynb \
            > {log} 2>&1
        """

rule ppi_cluster3:
    # NOTE: Hub genes identified via MCC + Betweenness from CytoHubba.
    # STRING annotation cross-reference (same approach as clusters 1,2).
    # Cross-reference logic for all three clusters documented here.
    # Produces two figures: method-colored and direction-colored network.
    # Direction-colored figure shows EZH2/HDAC2 as upregulated epigenetic
    # repressors within otherwise downregulated cluster. See DECISIONS.md.
    input:
        annotations      = "results/PPI/cluster3_annotations.tsv",
        annotations_c1   = "results/PPI/cluster1_annotations.tsv",
        annotations_c2   = "results/PPI/cluster2_annotations.tsv",
        annotations_cons = "results/PPI/consensus_annotations.tsv",
        mcc              = "results/PPI/cluster3_MCC_top10.csv",
        betweenness      = "results/PPI/cluster3_betweenness_top10.csv",
        edges            = "results/PPI/cluster3_edges.tsv",
        degrees          = "results/PPI/cluster3_degrees.tsv",
        consensus        = "results/DEG/consensus_DEGs.csv"
    output:
        network           = "figures/cluster3_hub_network.png",
        network_direction = "figures/cluster3_hub_network_direction.png"
    log:
        "logs/ppi_cluster3.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/PPI_cluster_3.ipynb \
            --output notebooks/PPI_cluster_3_executed.ipynb \
            > {log} 2>&1
        """

rule ppi_cluster4:
    # NOTE: Cluster 4 (37 nodes) — too small for CytoHubba.
    # Hub genes selected as top 5 by STRING node degree.
    # Degree centrality used as proxy for topological importance.
    # This notebook also handles enrichment visualization for cluster 4
    # and BH-filtered enrichment CSV output.
    # Enrichment figure shows top 15 by raw p-value; saved CSV uses
    # BH < 0.05 filter. See DECISIONS.md.
    # DEPENDS ON MANUAL CHECKPOINT 3 (DAVID enrichment CSV for cluster 4).
    input:
        annotations = "results/PPI/cluster4_annotations.tsv",
        edges       = "results/PPI/cluster4_edges.tsv",
        degrees     = "results/PPI/cluster4_degrees.tsv",
        consensus   = "results/DEG/consensus_DEGs.csv",
        enrichment  = "results/enrichment/cluster4_enrichment.csv"
    output:
        gene_list  = "results/enrichment/cluster4_genes.txt",
        network    = "figures/cluster4_network.png",
        enrichment = "figures/cluster4_enrichment.png",
        filtered   = "results/enrichment/cluster4_enrichment_filtered.csv"
    log:
        "logs/ppi_cluster4.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/cluster4_enrichment.ipynb \
            --output notebooks/cluster4_enrichment_executed.ipynb \
            > {log} 2>&1
        """

# Checkpoint: all PPI notebooks complete
rule stage4_complete:
    input:
        "figures/hub_gene_network.png",
        "figures/cluster1_hub_network.png",
        "figures/cluster2_hub_network.png",
        "figures/cluster3_hub_network.png",
        "figures/cluster4_network.png"
    output:
        touch("checkpoints/stage4_ppi.done")


# =============================================================================
# STAGE 5: HUB GENE VERIFICATION + NOTCH RESCUE + OVERLAP + DIRECTION CHECK
# =============================================================================

rule hub_gene_directional_verification:
    # NOTE: Verifies all hub genes are present in consensus DEGs
    # and confirms directional purity per cluster.
    # Hub gene lists explicitly defined — union of MCC + Betweenness
    # top 10 for clusters 1-3, top 5 by degree for cluster 4.
    # All 50 hub genes confirmed in consensus.
    # Per-cluster hub gene list files saved to results/DEG/hub_gene_lists/
    # use by all downstream notebooks. See DECISIONS.md.
    input:
        consensus = "results/DEG/consensus_DEGs.csv"
    output:
        hub_lists = directory("results/DEG/hub_gene_lists"),
        verification = "results/DEG/hub_gene_directional_verification.csv"
    log:
        "logs/hub_gene_directional_verification.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/hub_gene_directional_verification.ipynb \
            --output notebooks/hub_gene_directional_verification_executed.ipynb \
            > {log} 2>&1
        """

rule gene_rescued_by_notch:
    # NOTE: Cross-references whole-network hub genes against
    # Notch-rescued and Notch-suppressed consensus DEGs.
    # Result: CD44, CCND1, KDR rescued by Notch activation.
    # WARNING: Gene lists hardcoded from whole network —
    # to be replaced with hub gene list files in cleanup.
    # WARNING: No output file currently saved in notebook —
    # save step must be added during cleanup. See DECISIONS.md.
    input:
        rescued    = "results/DEG/notch_rescued_consensus_down.csv",
        suppressed = "results/DEG/notch_suppressed_consensus_up.csv",
        hub_lists  = "results/DEG/hub_gene_lists"
    output:
        rescued_hubs = "results/DEG/notch_rescued_hub_genes.csv",
        suppressed_hubs = "results/DEG/notch_suppressed_hub_genes.csv"
    log:
        "logs/gene_rescued_by_notch.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/gene_rescued_by_notch.ipynb \
            --output notebooks/gene_rescued_by_notch_executed.ipynb \
            > {log} 2>&1
        """

rule overlap_assessment:
    # NOTE: Documents gene space overlap across three primary datasets.
    # Three-way overlap = 15,504 genes — this is the consensus gene
    # universe used in Fisher's exact test in Consensus_DEGs.
    # GSE60052 is the limiting dataset after microarray union filtering.
    # See DECISIONS.md.
    input:
        gse30219  = "data/processed/GSE30219_expression_matrix.csv",
        gse149507 = "data/processed/GSE149507_filtered_expression_matrix_final.csv",
        gse60052  = "data/processed/GSE60052_expression_matrix.csv"
    output:
        venn = "figures/gene_overlap_venn.png"
    log:
        "logs/overlap_assessment.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/Overlap_assessment.ipynb \
            --output notebooks/Overlap_assessment_executed.ipynb \
            > {log} 2>&1
        """

rule direction_check_cluster4_5:
    # NOTE: Exploratory check confirming suppression dominance in
    # clusters 4 and 5 (94.6% and 96.2% downregulated respectively).
    # Justifies stopping analysis at cluster 4 — cluster 5 shows same
    # suppression pattern but excluded for scope manageability.
    # Cluster 4 retained because ECM/metastasis is clinically critical.
    # WARNING: No output file currently saved — save step must be
    # added during cleanup. See DECISIONS.md.
    input:
        cluster4  = "results/PPI/cluster4_annotations.tsv",
        cluster5  = "results/PPI/cluster5_annotations.tsv",
        consensus = "results/DEG/consensus_DEGs.csv"
    output:
        summary = "results/PPI/cluster4_5_direction_summary.csv"
    log:
        "logs/direction_check_cluster4_5.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute "notebooks/direction_check_for_cluster_4_and_5.ipynb" \
            --output "notebooks/direction_check_for_cluster_4_and_5_executed.ipynb" \
            > {log} 2>&1
        """


# =============================================================================
# MANUAL CHECKPOINT 3: DAVID functional enrichment
# -----------------------------------------------------------------------------
# For each cluster (1-4):
# 1. Upload results/enrichment/cluster{N}_genes.txt to DAVID
# 2. Run functional annotation
# 3. Download enrichment results as CSV
# 4. Place files in results/enrichment/:
#    cluster{N}_enrichment.csv
# Note: cluster4_genes.txt is produced by ppi_cluster4 rule above.
# After completing the above, continue pipeline from Stage 6.
# =============================================================================


# =============================================================================
# STAGE 6: ENRICHMENT
# full_enrichment_table handles clusters 1-3 (BH filtering).
# Cluster 4 enrichment handled inside ppi_cluster4 rule (Stage 4b).
# Individual cluster notebooks handle visualization per cluster.
# =============================================================================

rule full_enrichment_table:
    # NOTE: Applies category filter and BH < 0.05 to DAVID enrichment
    # results for clusters 1-3. Cluster 4 handled in ppi_cluster4 rule.
    # Categories kept: GOTERM_BP_DIRECT, GOTERM_CC_DIRECT,
    # GOTERM_MF_DIRECT, KEGG_PATHWAY.
    # Justification: GO Direct terms provide direct annotation without
    # inferred parent terms; KEGG captures pathway-level context.
    # DEPENDS ON MANUAL CHECKPOINT 3 (DAVID CSVs for clusters 1-3).
    # TODO (cleanup): extend loop to include cluster 4 for consistency.
    # See DECISIONS.md.
    input:
        c1 = "results/enrichment/cluster1_enrichment.csv",
        c2 = "results/enrichment/cluster2_enrichment.csv",
        c3 = "results/enrichment/cluster3_enrichment.csv"
    output:
        c1_filtered = "results/enrichment/cluster1_enrichment_filtered.csv",
        c2_filtered = "results/enrichment/cluster2_enrichment_filtered.csv",
        c3_filtered = "results/enrichment/cluster3_enrichment_filtered.csv"
    log:
        "logs/full_enrichment_table.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute "notebooks/full_enrichment_table.ipynb" \
            --output "notebooks/full_enrichment_table_executed.ipynb" \
            > {log} 2>&1
        """

rule enrichment_cluster1:
    # NOTE: Visualization reads from raw DAVID CSV directly —
    # shows top 15 by raw p-value for figure.
    # BH-filtered CSV produced by full_enrichment_table rule.
    # WARNING: Variable named cluster2_enrich in notebook —
    # copy-paste error, fix during cleanup.
    input:
        enrichment = "results/enrichment/cluster1_enrichment.csv",
        filtered   = "results/enrichment/cluster1_enrichment_filtered.csv"
    output:
        figure = "figures/cluster1_enrichment.png"
    log:
        "logs/enrichment_cluster1.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/cluster1_enrichment.ipynb \
            --output notebooks/cluster1_enrichment_executed.ipynb \
            > {log} 2>&1
        """

rule enrichment_cluster2:
    input:
        enrichment = "results/enrichment/cluster2_enrichment.csv",
        filtered   = "results/enrichment/cluster2_enrichment_filtered.csv"
    output:
        figure = "figures/cluster2_enrichment.png"
    log:
        "logs/enrichment_cluster2.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/cluster2_enrichment.ipynb \
            --output notebooks/cluster2_enrichment_executed.ipynb \
            > {log} 2>&1
        """

rule enrichment_cluster3:
    input:
        enrichment = "results/enrichment/cluster3_enrichment.csv",
        filtered   = "results/enrichment/cluster3_enrichment_filtered.csv"
    output:
        figure = "figures/cluster3_enrichment.png"
    log:
        "logs/enrichment_cluster3.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/cluster3_enrichment.ipynb \
            --output notebooks/cluster3_enrichment_executed.ipynb \
            > {log} 2>&1
        """


# =============================================================================
# STAGE 7: PROTEOMICS VALIDATION + SURVIVAL + PH VALIDATION
# =============================================================================

rule protein_validation_survival:
    # NOTE: Loads CPTAC proteomics (S3A) for hub gene validation.
    # Only S3A-validated genes taken to survival analysis (27 genes).
    # Median split justified by Hartigan's dip test (all 27 genes
    # unimodal, p>0.05 — no bimodal distribution warranting alternative).
    # KM plots include censored tick marks and at-risk/censored tables.
    # S4 investigated and documented as non-applicable by design —
    # Liu et al. S4 uses survival-first concordance pipeline, orthogonal
    # to our network topology approach.
    # files in cleanup. See DECISIONS.md.
    input:
        hub_lists = "results/DEG/hub_gene_lists",
        s3        = "results/proteomics/CPTAC_S3_tumor_NAT.xlsx",
        s4        = "results/proteomics/CPTAC_S4_survival.xlsx",
        s1        = "results/proteomics/CPTAC_S1_clinical.xlsb",
        consensus = "results/DEG/consensus_DEGs.csv"
    output:
        km_all = "figures/km_plots_all_hub_genes.png",
        km_sig = "figures/km_plots_significant.png",
        expr   = "figures/expression_distributions.png",
        hr     = "results/proteomics/km_survival_hr_results.csv"
    log:
        "logs/protein_validation_survival.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute notebooks/protein_validation_survival_plot.ipynb \
            --output notebooks/protein_validation_survival_plot_executed.ipynb \
            > {log} 2>&1
        """

rule proportional_hazard_bh:
    # NOTE: BH correction applied to 11 significant/borderline genes
    # (p<0.10 from initial KM screen) — not all 27 validated genes.
    # Schoenfeld residuals confirm PH assumption met for all 11 genes
    # (all p>0.05). Three genes survive BH correction: CASP1, IL18, TLR2.
    # WARNING: Duplicates data loading from survival notebook —
    # consider merging or saving intermediate merged_df in cleanup.
    # See DECISIONS.md.
    input:
        hub_lists = "results/DEG/hub_gene_lists",
        s1 = "results/proteomics/CPTAC_S1_clinical.xlsb",
        hr = "results/proteomics/km_survival_hr_results.csv"
    output:
        bh         = "results/proteomics/bh_correction_results.csv",
        schoenfeld = "results/proteomics/schoenfeld_ph_test_results.csv"
    log:
        "logs/proportional_hazard_bh.log"
    shell:
        """
        jupyter nbconvert --to notebook --execute "notebooks/proportional_hazard_assumption_and_multiple_testing_correction.ipynb" \
            --output "notebooks/proportional_hazard_assumption_and_multiple_testing_correction_executed.ipynb" \
            > {log} 2>&1
        """
