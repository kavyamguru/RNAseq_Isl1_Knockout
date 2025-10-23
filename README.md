# 🧬 RNAseq_Isl1_Knockout  
**Functional Genomic Technologies (PGBI11040)**  

---

## 🔍 Overview  
This repository contains the full workflow and analysis for the **Isl1 knockout vs wild-type allantois RNA-seq project** (GEO: GSE247182), completed as part of the *Functional Genomic Technologies* module.  
The project explores how **Isl1** regulates transcriptional programs during early mesoderm development using an R-based RNA-seq analysis pipeline.

---

## ⚙️ Workflow Summary  
1. **Quality Control** – FastQC and MultiQC  
2. **Alignment** – STAR (GRCm39 genome build)  
3. **Quantification** – featureCounts (Rsubread)  
4. **Normalization & DEG analysis** – DESeq2  
5. **Annotation** – biomaRt  
6. **Pathway Enrichment** – fgsea (MSigDB Hallmark sets)  
7. **Visualization** – ggplot2, pheatmap, EnhancedVolcano  

---

## 📊 Key Results  
- **9 significant DEGs (FDR < 0.05)** including *Prl8a2*, *Perp*, *Lefty2*  
- Enriched pathways: **Epithelial-Mesenchymal Transition (EMT)**, **Hypoxia**, **Cholesterol Homeostasis**  
- PCA revealed strong genotype-dependent clustering  
- Demonstrated a complete and reproducible RNA-seq analysis workflow  

---

## 📁 Repository Structure  

| File / Folder | Description |
|----------------|-------------|
| `RNAseq_Isl1_Knockout.R` | Main R script implementing DESeq2, fgsea, and visualization |
| `Figures.zip` | PCA, MA, Volcano, Heatmap, and GSEA result plots |
| `Tables.zip` | DESeq2 and fgsea summary tables |
| `QC.zip` | MultiQC summary report (open in browser) |
| `Project report.pdf` | Final MSc project report detailing analysis and results |

---

## 📈 Example Figures  
*(To view all figures, download `Figures.zip`)*  
- **PCA Plot:** Clear separation of knockout vs wild-type samples  
- **Volcano Plot:** Differential expression visualization  
- **Heatmap:** Top 50 DEGs by variance  

---

## 🧠 Learning Outcome  
This project demonstrates:  
- Proficiency in RNA-seq data processing and statistical analysis in R  
- Familiarity with genomic alignment and annotation tools  
- Integration of pathway enrichment and biological interpretation  
- Reproducible research and documentation practices  

---

## 🧩 Tools & Packages  
`R`, `DESeq2`, `Rsubread`, `biomaRt`, `fgsea`, `ggplot2`, `pheatmap`, `EnhancedVolcano`, `STAR`, `FastQC`, `MultiQC`

---

