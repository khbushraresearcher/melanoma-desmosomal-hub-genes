# melanoma-desmosomal-hub-genes
Code and processed data for the study of prognostic desmosomal hub genes in melanoma
# Integrated Transcriptomic and Single-Cell Analysis of Desmosomal Hub Genes in Melanoma

## Overview

This repository contains the computational scripts, intermediate analysis objects, result tables, and figures associated with the manuscript:

**“WGCNA and Single-Cell Analysis Identify FLG, DSG3 and DSG1 as Prognostic Desmosomal Hub Genes in Melanoma and Their Modulation by Dinaciclib”**

The study integrates transcriptomic, network, clinical, prognostic, immune-infiltration, and single-cell analyses to characterize desmosome-associated genes in melanoma.

The computational workflow includes differential expression analysis, batch correction, weighted gene co-expression network analysis (WGCNA), meta-analysis, protein–protein interaction analysis, TCGA-SKCM prognostic modelling, GSEA, TIMER-based immune infiltration analysis, and single-cell RNA-seq analysis.

---

## Study Objectives

The main objectives of the study were:

1. To identify differentially expressed genes associated with melanoma.
2. To identify melanoma-associated co-expression modules and hub genes using WGCNA.
3. To evaluate the reproducibility of candidate gene associations across independent melanoma transcriptomic datasets.
4. To investigate the prognostic relevance of desmosomal genes in the TCGA-SKCM cohort.
5. To develop and evaluate a multigene prognostic RiskScore using LASSO-Cox regression.
6. To investigate the association of candidate genes with tumor immune-cell infiltration using TIMER.
7. To characterize the cellular distribution of candidate genes using single-cell RNA sequencing.
8. To investigate the effects of Dinaciclib on melanoma cell viability, cell-cycle distribution, and expression of selected desmosomal genes.

---

## Datasets

### Bulk transcriptomic datasets

The following Gene Expression Omnibus (GEO) datasets were used for bulk transcriptomic analyses:

* **GSE3189**
* **GSE7553**
* **GSE46517**

These datasets were analyzed individually and, where applicable, after batch correction and integration.

### Single-cell RNA-seq dataset

* **GSE72056**

This dataset was used to investigate the cellular distribution and expression patterns of FLG, DSG1, and DSG3 in melanoma and associated cell populations.

### TCGA dataset

Prognostic analyses were performed using:

* **TCGA-SKCM**

TCGA-SKCM expression and clinical information were used for survival analysis, LASSO-Cox modelling, multivariable Cox regression, calibration analysis, and related prognostic analyses.

---

## Analysis Workflow

The computational workflow comprised the following major steps:

```text
Public melanoma transcriptomic datasets
                │
                ▼
      Data preprocessing
                │
                ▼
 Differential expression analysis
                │
                ▼
       ComBat batch correction
                │
                ▼
              WGCNA
                │
                ▼
     DExMA meta-analysis
                │
                ▼
        PPI / hub analysis
                │
                ▼
       Candidate desmosomal genes
                │
        ┌───────┼────────┐
        ▼       ▼        ▼
     TCGA    TIMER    scRNA-seq
   Prognosis  immune   cellular
              infiltration distribution
        │
        ▼
    LASSO-Cox
        │
        ▼
    RiskScore
        │
        ▼
 Multivariable Cox
        │
        ▼
 Calibration / validation

Additional pathway analysis:
GSEA → GO Biological Process + Hallmark pathways
```

---

## Differential Expression Analysis

Differential expression analysis was performed on the melanoma transcriptomic datasets using the `limma` framework.

Genes were considered differentially expressed based on:

* Adjusted P value < 0.05
* |log2 fold change| > 1

The individual datasets were analyzed separately before integration of selected datasets.

---

## Batch Correction

Batch effects among the integrated transcriptomic datasets were addressed using **ComBat**.

The batch-corrected expression matrix was subsequently used for downstream analyses where appropriate.

---

## Weighted Gene Co-expression Network Analysis

Weighted gene co-expression network analysis (WGCNA) was used to identify groups of co-expressed genes and melanoma-associated modules.

The analysis included:

* selection of the WGCNA soft-thresholding power;
* construction of the gene co-expression network;
* identification of gene modules;
* correlation of modules with melanoma-related phenotypes;
* identification of highly connected genes within relevant modules.

The WGCNA analyses were performed separately for the relevant GEO datasets and the integrated dataset.

---

## DExMA Meta-analysis

A random-effects meta-analysis was performed using **DExMA** to evaluate the consistency of gene-expression changes across independent melanoma datasets.

The analysis was used to assess the reproducibility of the expression patterns of candidate desmosomal genes across datasets.

---

## Protein–Protein Interaction and Hub Gene Analysis

Protein–protein interaction analysis was performed to further characterize the candidate genes identified through the transcriptomic and network analyses.

The desmosomal candidate panel included:

* **FLG**
* **DSG1**
* **DSG3**
* **DSC1**
* **DSP**

Network-based analysis was used to identify genes with prominent connectivity within the candidate interaction network.

---

## TCGA-SKCM Prognostic Analysis

The prognostic relevance of the candidate desmosomal genes was evaluated using the TCGA-SKCM cohort.

Survival analyses included:

* individual-gene Cox regression;
* LASSO-Cox regression;
* multivariable Cox regression;
* RiskScore-based survival analysis;
* time-dependent evaluation;
* calibration analysis;
* bootstrap-based assessment.

The multivariable prognostic model incorporated:

* RiskScore
* Age
* Gender
* Stage

Calibration was evaluated at **1-, 3-, and 5-year** time points using bootstrap resampling. The calibration analysis used 200 bootstrap repetitions. 

---

## LASSO-Cox RiskScore

The candidate desmosomal genes evaluated for prognostic modelling were:

```text
FLG
DSG3
DSG1
DSP
DSC1
```

LASSO-Cox regression was used to derive the prognostic RiskScore.

The final model coefficients and associated model information are provided in the corresponding result files and RData objects.

---

## Immune Infiltration Analysis

Tumor immune-cell infiltration was additionally investigated using the **TIMER online platform**.

TIMER-based analysis was used to examine the relationship between candidate desmosomal gene expression and tumor-infiltrating immune-cell populations.

These analyses provided complementary information regarding the potential relationship between the identified desmosomal gene program and the melanoma tumor microenvironment.

---

## Gene Set Enrichment Analysis

Gene set enrichment analysis (GSEA) was performed to investigate biological pathways associated with melanoma transcriptomic changes.

GSEA included:

* **Gene Ontology (GO) Biological Process**
* **MSigDB Hallmark gene sets**

The analysis was performed using ranked gene lists derived from differential-expression statistics.

The GSEA workflow included the relevant individual melanoma datasets and the integrated analysis where applicable.

---

## Single-Cell RNA-seq Analysis

Single-cell RNA-seq analysis was performed using **GSE72056**.

The workflow included:

* normalization;
* identification of variable features;
* scaling;
* principal component analysis;
* clustering;
* UMAP visualization;
* cell-type annotation;
* expression analysis of candidate genes;
* melanoma-versus-keratinocyte comparison;
* DSG3-positive versus DSG3-negative melanoma-cell analysis;
* pathway analysis;
* pseudotime analysis.

Particular attention was given to **FLG, DSG1, and DSG3** to assess their cellular distribution and expression heterogeneity.

---

## Prognostic Model Calibration

Calibration analysis was performed for the multivariable Cox model at:

* 1 year
* 3 years
* 5 years

The calibration model included RiskScore, Age, Gender, and Stage_group. 

The generated calibration datasets and figures are saved in the `Results/Prognostic/` and `Figures/` directories. 

---

## Repository Structure

The repository is organized into the following main directories:

```text
Melanoma_Project2/
│
├── Data/
│   └── Input and supporting datasets
│
├── Scripts/
│   └── R scripts used for the computational analyses
│
├── Results/
│   ├── Prognostic/
│   ├── GSEA/
│   └── Other analysis-specific results
│
├── Figures/
│   └── Generated figures and plots
│
├── RData/
│   └── Saved R analysis objects
│
└── README.md
```

The individual R scripts are organized according to the analytical workflow used in the study.

---

## Software and R Packages

The computational analyses were performed using R and established CRAN/Bioconductor packages.

Major packages used across the analyses include:

* `limma`
* `affy`
* `WGCNA`
* `DExMA`
* `clusterProfiler`
* `org.Hs.eg.db`
* `msigdbr`
* `glmnet`
* `survival`
* `rms`
* `Seurat`
* `slingshot`
* `enrichplot`
* `ggplot2`
* `dplyr`

The specific packages required for each analysis are loaded within the corresponding R script.

---

## Reproducibility

The repository provides the computational scripts and supporting analysis objects used for the analyses reported in the manuscript.

The scripts cover:

* transcriptomic preprocessing;
* differential expression analysis;
* batch correction;
* WGCNA;
* meta-analysis;
* PPI and hub-gene analysis;
* TCGA-SKCM prognostic modelling;
* LASSO-Cox RiskScore development;
* GSEA;
* single-cell analysis;
* pseudotime analysis;
* generation of analysis outputs and figures.

Intermediate RData objects and result tables are provided where applicable to facilitate reproduction of the computational workflow.

---

## Results

The repository contains analysis outputs generated during the study, including:

* differential-expression results;
* WGCNA module results;
* DExMA meta-analysis results;
* PPI/hub-gene results;
* TCGA-SKCM prognostic results;
* LASSO model information;
* RiskScore coefficients;
* calibration results;
* GSEA results;
* single-cell analysis results;
* generated figures.

---

## Data Availability

The transcriptomic datasets used in this study are publicly available through the **NCBI Gene Expression Omnibus (GEO)** under the accession numbers listed above.

The TCGA-SKCM data used for prognostic analysis are publicly available through the **TCGA** data resources.

The immune-infiltration analyses were performed using the **TIMER online platform**.

---

## Citation

If you use the scripts or analyses from this repository, please cite the associated manuscript:

**Khan B, Imtiyaz K, Athar A, Nasar MA, Abbasi B, Aslam A, Rizvi MA.**

*WGCNA and Single-Cell Analysis Identify FLG, DSG3 and DSG1 as Prognostic Desmosomal Hub Genes in Melanoma and Their Modulation by Dinaciclib.*

---

## Contact

**Bushra Khan**
Department of Biosciences
Jamia Millia Islamia
New Delhi, India

