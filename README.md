# RNAseq_Isl1_Knockout

RNA-seq differential expression and pathway analysis for **Isl1 knockout vs wild-type mouse allantois** (GEO: GSE247182).

## Objective
Identify transcriptional and pathway-level effects associated with Isl1 loss using a reproducible R-based analysis workflow.

## Workflow
1. **QC**: FastQC + MultiQC
2. **Alignment**: STAR (GRCm39)
3. **Quantification**: featureCounts / Rsubread
4. **Differential expression**: DESeq2
5. **Annotation**: biomaRt
6. **Enrichment**: fgsea (MSigDB Hallmark)
7. **Visualization**: PCA, MA, volcano, heatmaps

## Repository Contents
- `RNAseq_lsl1_Knockout.R` — end-to-end analysis script
- `QC.zip` — QC reports
- `Figures.zip` — visual outputs
- `Tables.zip` — DEG/enrichment outputs
- `Project report.pdf` — detailed academic report

## Key Findings (summary)
- Significant differential expression profile between KO and WT
- Pathway shifts including EMT/hypoxia-associated signatures
- Distinct genotype clustering on PCA

## Tools
R, DESeq2, fgsea, biomaRt, ggplot2, pheatmap, EnhancedVolcano, STAR, FastQC, MultiQC

## Reproducibility Notes
- Analysis logic is consolidated in a single script for traceability.
- Figures/tables are bundled for reviewer-friendly inspection.
- Intended as an educational + portfolio demonstration of RNA-seq analysis competency.
