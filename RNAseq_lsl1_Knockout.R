# This R script is adapted from the Functional Genomics and Transcriptomics lectures 
# taught by Dr. Simon Tomlinson. Modifications were made to tailor the analysis to the current dataset, 
# including adjustment of statistical thresholds (e.g., log2 fold change and FDR cutoffs) 
# and customization of plots (PCA, MA, heatmaps, volcano plots, etc.) to ensure that all output figures 
# meet publication-quality standards. The core analytical workflow remains consistent with the original course materials.

# ----------------------------

# Data Preprocessing Summary:
# ----------------------------
# FASTQ files for GSE247182 were pre-downloaded and placed in the 'fastq' directory.
# FastQC had already been run, and results were available in the 'fastqc' folder.

# Trimming:
# No trimming was required. FastQC reports showed high-quality reads with no significant adapter contamination.
# To proceed with the alignment pipeline (which expected input from a 'fastq_tr' directory),
# a symbolic link was created to point 'fastq_tr' to the original 'fastq' directory:
# Command used (on server):
# ln -s fastq fastq_tr

# Alignment was then performed using STAR aligner on reads from 'fastq_tr'.

# Post-alignment QC:
# MultiQC was run after alignment to aggregate QC metrics from both FASTQ and BAM alignment outputs:
# multiqc -n GSE247182_multiqc_report ./fastqc ./alignment

# ----------------------------

# Load necessary R packages for RNA-seq analysis
library(Rsubread)        # For alignment and read counting (optional, depending on input)
library(edgeR)           # For expression filtering and normalization (CPM-based)
library(limma)           # For linear modeling and additional visualization
library(gplots)          # Base-level heatmap visualizations
library(DESeq2)          # Core tool for RNA-seq differential expression analysis
library(affy)            # For classical plotting tools like MA plots
library(QuasR)           # Optional: mapping and quantification tool
library(pheatmap)        # For visually appealing heatmaps
library(EnhancedVolcano) # Publication-ready volcano plots
library(stringr)         # String manipulation utilities
library(ggplot2)         # For PCA and customizable plots
library(biomaRt)         # For Ensembl annotation retrieval
library(ggrepel)         # Prevents label overlap in ggplot2
library(RColorBrewer)    # Color palettes for heatmaps and plots

# Set a seed to ensure reproducibility across random steps
set.seed(42)

# Load the featureCounts output (alignment & counting results)
load("mytable_feaures")  # This .RData contains read counts per gene and annotation

# After loading, verify the object and extract count matrix
ls()                     # List objects in workspace (should include 'mytable_feaures')

# Assume mytable_feaures is a list-like object with a 'counts' element containing the count matrix
counts_matrix <- (mytable_feaures)$counts   # raw counts matrix (genes x samples)
dim(counts_matrix)       # Check dimensions: number of genes and samples
head(rownames(counts_matrix))  # Check gene identifiers (likely Ensembl IDs)
head(colnames(counts_matrix))  # Check sample identifiers (filenames or IDs)

# Build a sample annotation data frame (adf) from the count matrix column names
sample_names <- colnames(counts_matrix)
adf <- data.frame(FileName = sample_names)         # start with filenames
# Assume filenames contain a common suffix "Aligned.sortedByCoord.out.bam" from STAR aligner
pattern <- "Aligned.sortedByCoord.out.bam"
loc <- str_locate(sample_names, pattern)           # find the position of the pattern in each filename&#8203;:contentReference[oaicite:10]{index=10}
adf$SampleID <- substr(sample_names, 1, loc[,"start"]-1)
adf

# Assign short sample names in order (assuming the counts columns are in the same order as GSM IDs)
adf$ShortName <- c("WT_1","WT_2","WT_3","Isl1_KO_1","Isl1_KO_2","Isl1_KO_3") # short labels
# Add Genotype (experimental condition) based on the short name
adf$Genotype <- factor(ifelse(grepl("^WT", adf$ShortName), "WT", "Isl1_KO"))

# Assign colors for each genotype for consistent plotting (e.g., WT = first color, MUT = second)
adf$Color <- rainbow(length(levels(adf$Genotype)))[adf$Genotype]  # color per factor level

# Finalize the annotation data frame
rownames(adf) <- adf$ShortName           # set row names to short names
# Ensure count matrix columns match these rownames
colnames(counts_matrix) <- adf$ShortName  # replace column names with short names for clarity

# Check alignment between counts and annotation
stopifnot( all(colnames(counts_matrix) == rownames(adf)) )
print(adf)

# Create DESeq2 dataset from count matrix and sample annotation
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData   = adf,
                              design    = ~ Genotype)  # design formula (one-factor design)

# Inspect the DESeqDataSet object
dim(dds)                   # number of genes and samples in the DESeq2 object
head(rownames(dds), 3)     # first 3 gene IDs
colnames(colData(dds))     # metadata columns (should include "Genotype")
table(dds$Genotype)        # sample count per genotype level

dds$Genotype <- relevel(dds$Genotype, ref="WT")  # make WT the reference level
levels(dds$Genotype)  # check levels (WT should be first)

# Filter out genes with very low counts across all samples
before <- nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1, ]   # keep genes with total count > 1
after <- nrow(dds)
cat("Genes before filtering:", before, "; after filtering:", after, "\n")

# Estimate size factors for normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Apply variance-stabilizing transformations for visualization
vsd <- vst(dds, blind=TRUE)   # variance stabilizing transformation (fast)
rld <- rlog(dds, blind=TRUE)  # regularized log transformation&
# Boxplot of transformed counts (visual QC)
png("boxplots_vst_rlog.png", width = 10, height = 5, units = "in", res = 300)

# Set side-by-side layout
par(mfrow = c(1, 2), mar = c(6, 5, 4, 2))  # increase bottom margin for rotated labels

# Custom colors
colors_vst <- c(rep("gray70", 3), rep("skyblue3", 3))
colors_rlog <- c(rep("gray70", 3), rep("skyblue3", 3))

# VST boxplot
boxplot(assay(vsd),
        main = "VST-transformed Counts",
        ylab = "VST counts",
        col = colors_vst,
        las = 2,
        outline = FALSE,
        cex.axis = 0.9,
        cex.lab = 1.2,
        cex.main = 1.4)

# rlog boxplot
boxplot(assay(rld),
        main = "rlog-transformed Counts",
        ylab = "rlog counts",
        col = colors_rlog,
        las = 2,
        outline = FALSE,
        cex.axis = 0.9,
        cex.lab = 1.2,
        cex.main = 1.4)

# Reset layout
par(mfrow = c(1, 1))

dev.off()
# Examine the assay of rld (rlog values) for a few genes and samples
assay(rld)[1:3, 1:3]

# Run the differential expression analysis
dds <- DESeq(dds)  # this performs normalization (already done), dispersion estimation, fitting, and testing&#8203;:contentReference[oaicite:33]{index=33}

# Obtain results for the Genotype effect (MUT vs WT)
res <- results(dds, contrast = c("Genotype", "Isl1_KO", "WT"))
# 'contrast' specifies comparison: here MUT vs WT, so log2FC = log2(MUT/WT)&#8203;:contentReference[oaicite:34]{index=34}

# Summary of results
summary(res)

# Remove genes with NA padj (not testable, e.g., low counts)
res <- res[!is.na(res$padj), ]

# How many significant genes at FDR < 0.05?
sig <- sum(res$padj < 0.05)
cat("Significant genes (FDR<0.05):", sig, "\n")

# Order results by p-value (or by another criterion if needed)
res_ordered <- res[order(res$padj), ]
head(res_ordered, 10)   # top 10 most significant genes

# PCA on rlog-transformed data
plotPCA(rld, intgroup="Genotype")

p <- plotPCA(rld, intgroup="Genotype", returnData=TRUE)
percentVar <- round(100 * attr(p, "percentVar"))
# `p` contains PC coordinates and metadata; we use ggplot for labeling

ggplot(p, aes(PC1, PC2, color = Genotype, label = name)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 4, fontface = "bold") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("Principal Component Analysis (PCA) of RNA-seq Samples") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  scale_color_manual(values = c("WT" = "#E76F51", "Isl1_KO" = "#2A9D8F"))

ggsave("PCA_plot.png", width = 8, height = 6, dpi = 300)

png("MA_plot_Isl1_KO_vs_WT.png", width = 8, height = 6, units = "in", res = 300)

# MA plot of the results
plotMA(res,
       ylim = c(-6, 6),# significance threshold
       main = "MA Plot: Isl1_KO vs WT",        # improved title
       cex = 1,                                 # point size
       colSig = "#E76F51",                      # color for significant points
       colLine = "black",                       # color of baseline
       lwd = 1.5)                               # thickness of baseline

dev.off()

table("padj<0.05" = res$padj < 0.05)
table("padj<0.1" = res$padj < 0.1)

# Select top 30 DEGs by smallest padj
top30_genes <- rownames(res_ordered)[1:30]
# Retrieve rlog expression values for these genes
mat <- assay(rld)[top30_genes, ]
# Scale expression by gene (z-score per row) for better visualization of relative changes
mat_scaled <- t(scale(t(mat)))  # center and scale each gene's expression across samples

# Ensure adf exists and has the proper structure
annotation_col <- data.frame(Genotype = adf$Genotype)
rownames(annotation_col) <- rownames(adf)  # This must match column names of mat_scaled


# Improved color palette for better contrast
color_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

# Better column annotation formatting
annotation_colors <- list(Genotype = c("WT" = "#E76F51", "Isl1_KO" = "#2A9D8F"))

# Generate the heatmap and capture the output
heat <- pheatmap(mat_scaled,
                 annotation_col = annotation_col,
                 annotation_colors = annotation_colors,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 color = color_palette,
                 fontsize_row = 8,
                 fontsize_col = 10,
                 main = "",  
                 cellwidth = 25,
                 cellheight = 10,
                 treeheight_row = 40,
                 treeheight_col = 40,
                 show_colnames = TRUE,
                 show_rownames = TRUE,
                 border_color = NA)

# Add centered title manually
grid.text("Top 30 Differentially Expressed Genes (rlog Z-scores)",
          x = 0.5, y = 0.96, gp = gpar(fontsize = 16, fontface = "bold"))


# Capture plot as a grob
heat_grob <- grid::grid.grab()

# Save with ggsave
ggsave("Top30_DEGs_heatmap_ggsave.png", plot = heat_grob, width = 8, height = 10, dpi = 300)

# Use biomaRt to retrieve gene annotations from Ensembl
ensembl_host <- "https://www.ensembl.org"  # use main Ensembl (could switch to archive if needed)
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host=ensembl_host)
gene_ids <- rownames(dds)  # Ensembl gene IDs present in our dataset (post-filtering)

annot <- getBM(
  filters = "ensembl_gene_id",
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene_id",  
    "description",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  ),
  values = rownames(dds),  # or rownames(res)
  mart = mart
)

# Remove duplicate Ensembl IDs (keep the first occurrence)
annot_unique <- annot[!duplicated(annot$ensembl_gene_id), ]

# Now create the data frame safely
annot_df <- data.frame(annot_unique, row.names = annot_unique$ensembl_gene_id)

# Confirm it looks good
head(annot_df)

# Ensure all our result genes are in the annotation (they should be, if genome build matches)
all(rownames(res) %in% rownames(annot_df))

# Merge annotation into results
res_annot <- as.data.frame(res)
res_annot$Gene <- annot_df[rownames(res), "external_gene_name"]
res_annot$Description <- annot_df[rownames(res), "description"]

# Create a combined ID_name for unique labeling if needed
res_annot$ID_Name <- paste(rownames(res_annot), res_annot$Gene, sep="_")
head(res_annot, 3)

write.csv(res_annot, "DESeq2_results_annotated.csv")

# Select top 10 DEGs by adjusted p-value
top10_df <- res_annot[order(res_annot$padj), ][1:10, ]

# Select and rename the relevant columns
top10_table <- top10_df[, c("Gene", "ID_Name", "log2FoldChange", "padj", "Description")]
colnames(top10_table) <- c("Gene Symbol", "Ensembl ID", "log₂ Fold Change", "Adjusted p-value", "Description")

# View the table
print(top10_table)

# Optionally, write to CSV for including in your report
write.csv(top10_table, "Top10_DEGs_table.csv", row.names = FALSE)

# Identify strictly significant DEGs
strict_genes <- res_annot %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(Gene)

# Plot with strict thresholds
volcano_strict <- EnhancedVolcano(res_annot,
                                  lab = res_annot$Gene,
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  selectLab = strict_genes,
                                  xlab = bquote(~Log[2]~ "Fold Change (Isl1_KO vs WT)"),
                                  ylab = bquote(~-log[10]~ "Adjusted P-value"),
                                  pCutoff = 0.05,
                                  FCcutoff = 1.0,
                                  pointSize = 2.0,
                                  labSize = 4.5,
                                  colAlpha = 0.7,
                                  boxedLabels = TRUE,
                                  drawConnectors = TRUE,
                                  widthConnectors = 0.4,
                                  max.overlaps = 25,
                                  title = "Strict Volcano Plot",
                                  subtitle = "padj < 0.05 and |log₂FC| > 1",
                                  caption = "Red = Statistically & biologically significant genes")

# Save Plot 1 (optional)
ggsave("volcano_strict.png", plot = volcano_strict, width = 8, height = 6, dpi = 300)

# Create the plot object
top10_genes <- res_annot %>%
  filter(pvalue < 0.05 & abs(log2FoldChange) > 0.38) %>%
  arrange(pvalue) %>%
  slice(1:10) %>%
  pull(Gene)

# Create the volcano plot for top 10 significant genes (relaxed thresholds)
volcano_top10 <- EnhancedVolcano(res_annot,
                                 lab = res_annot$Gene,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 selectLab = top10_genes,
                                 xlab = bquote(~Log[2]~ "Fold Change (Isl1_KO vs WT)"),
                                 ylab = bquote(~-log[10]~ "Adjusted P-value"),
                                 pCutoff = 0.05,
                                 FCcutoff = 0.38,  # ~1.3 fold change
                                 pointSize = 2.0,
                                 labSize = 4.5,
                                 colAlpha = 0.7,
                                 boxedLabels = TRUE,
                                 drawConnectors = TRUE,
                                 widthConnectors = 0.4,
                                 max.overlaps = 25,
                                 title = "Volcano Plot of Top 10 Significant Genes",
                                 subtitle = "padj < 0.05 and |log₂FC| > 0.38 (1.3x FC)",
                                 caption = "Only top 10 DEGs based on p-value are labeled")

# Save the volcano plot as PNG
ggsave("volcano_top10_relaxed.png", plot = volcano_top10, width = 8, height = 6, dpi = 300)


# Download mouse hallmark gene sets if not already present
if (!file.exists("mouse_H_v5p2.rdata")) {
  download.file("https://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata", 
                destfile = "mouse_H_v5p2.rdata")
}

load("mouse_H_v5p2.rdata")  # loads an object (e.g., Mm.H) with hallmark gene sets
ls()[grepl("Mm.H", ls())]

# Prepare ranked vector for GSEA
res_df <- as.data.frame(res)  
res_df$entrez <- annot_df[rownames(res_df), "entrezgene_id"]  # map Ensembl to Entrez using our annot_df
res_df <- res_df[!is.na(res_df$entrez), ]       # drop if any without Entrez ID
# Remove duplicate Entrez IDs (if multiple Ensembl map to same Entrez, keep the one with lowest p)
res_df <- res_df[order(res_df$pvalue), ]
res_df <- res_df[!duplicated(res_df$entrez), ]

# Create named vector of stat values, names = Entrez IDs
ranked_stats <- res_df$stat
names(ranked_stats) <- res_df$entrez
head(ranked_stats)

library(fgsea)
fgseaRes <- fgsea(pathways = Mm.H, 
                  stats    = ranked_stats,
                  minSize  = 15,    # require at least 15 genes in a gene set
                  maxSize  = 500)   # and at most 500 genes (hallmarks typically ~100-200 genes)

# Examine FGSEA output
fgseaRes <- fgseaRes[order(fgseaRes$pval), ]
head(fgseaRes[, c("pathway","NES","padj")], 10)

sigPathways <- fgseaRes$pathway[fgseaRes$padj < 0.05]
sigPathways

# Create the enrichment plot object using fgsea's plotEnrichment
emt_plot <- plotEnrichment(Mm.H$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, ranked_stats) +
  labs(
    title = "GSEA: Epithelial-Mesenchymal Transition (EMT)",
    subtitle = "Significantly enriched in MUT vs WT",
    x = "Rank of Genes",
    y = "Enrichment Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = c("green")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

# Save the polished plot
ggsave("GSEA_EMT_Enrichment_MUT_vs_WT.png", plot = emt_plot, width = 8, height = 6, dpi = 300)



