#Step 1: Create the project structure
############################################################
## Create Melanoma Project Folder
############################################################

setwd("~/Desktop")   # Change if your Desktop path is different

dir.create("Melanoma_Project", showWarnings = FALSE)

setwd("Melanoma_Project")

dir.create("Data", showWarnings = FALSE)
dir.create("Scripts", showWarnings = FALSE)
dir.create("Results", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)
dir.create("RData", showWarnings = FALSE)


#Step 2: Move your datasets into the Data folder

file.rename("~/Desktop/GSE3189_expression_data.csv",
            "~/Desktop/Melanoma_Project/Data/GSE3189_expression_data.csv")

file.rename("~/Desktop/GSE46517_expression_data.csv",
            "~/Desktop/Melanoma_Project/Data/GSE46517_expression_data.csv")

file.rename("~/Desktop/GSE7553_expression_data.csv",
            "~/Desktop/Melanoma_Project/Data/GSE7553_expression_data.csv")

file.rename("~/Desktop/TCGA-SKCM.star_fpkm.tsv",
            "~/Desktop/Melanoma_Project/Data/TCGA-SKCM.star_fpkm.tsv")

file.rename("~/Desktop/TCGA-SKCM.clinical 3.tsv",
            "~/Desktop/Melanoma_Project/Data/TCGA-SKCM.clinical 3.tsv")

file.rename("~/Desktop/GSE72056_melanoma_single_cell_revised_v2.txt",
            "~/Desktop/Melanoma_Project/Data/GSE72056_melanoma_single_cell_revised_v2.txt")


#Step 3: Set the working directory

setwd("~/Desktop/Melanoma_Project")


#Step 4: Read files

install.packages("data.table")   # Run only once if not already installed
library(data.table)

GSE3189 <- fread("Data/GSE3189_expression_data.csv",
                 data.table = FALSE)

GSE46517 <- fread("Data/GSE46517_expression_data.csv",
                  data.table = FALSE)

GSE7553 <- fread("Data/GSE7553_expression_data.csv",
                 data.table = FALSE)

TCGA.exp <- fread("Data/TCGA-SKCM.star_fpkm.tsv",
                  data.table = FALSE)

TCGA.clin <- fread("Data/TCGA-SKCM.clinical 3.tsv",
                   data.table = FALSE)

GSE72056 <- fread("Data/GSE72056_melanoma_single_cell_revised_v2.txt",
                  sep = "\t",
                  data.table = FALSE)


save.image("RData/Step01_Loaded_Data.RData")

load("RData/Step01_Loaded_Data.RData")



#STEP 2: Preprocess GSE3189
#2.1 Load annotation packages

############################################################
## Load Annotation Packages
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "AnnotationDbi",
  "hgu133a.db",
  "limma"
))

library(AnnotationDbi)
library(hgu133a.db)
library(limma)


#2.2 Check dimensions

dim(GSE3189)

head(GSE3189[,1:5])


#2.3 Extract probe IDs

############################################################
## Extract Probe IDs
############################################################

probe3189 <- GSE3189[,1]

head(probe3189)
length(probe3189)


#2.4 Map probes to gene symbols

############################################################
## Probe Annotation
############################################################

gene3189 <- mapIds(
  hgu133a.db,
  keys = probe3189,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

head(gene3189)


#2.5 Add gene symbols

GSE3189$GeneSymbol <- gene3189

head(GSE3189)


#2.6 Remove probes without symbols

############################################################
## Remove NA genes
############################################################

GSE3189 <- GSE3189[
  !is.na(GSE3189$GeneSymbol),
]

dim(GSE3189)


#2.7 Create expression matrix

############################################################
## Expression Matrix
############################################################

expr3189 <- as.matrix(
  
  GSE3189[
    ,
    !(names(GSE3189) %in%
        c("V1","GeneSymbol"))
  ]
  
)

rownames(expr3189) <- GSE3189$GeneSymbol

mode(expr3189) <- "numeric"

dim(expr3189)


#2.8 Collapse duplicated genes

############################################################
## Average duplicate probes
############################################################

expr3189 <- avereps(
  expr3189,
  ID = rownames(expr3189)
)

dim(expr3189)


#2.9 Check whether log2 transformation is needed

summary(as.numeric(expr3189))

quantile(
  as.numeric(expr3189),
  probs=c(0,0.25,0.5,0.75,0.99,1)
)

max(expr3189)


expr3189 <- log2(expr3189 + 1)


#2.10 Verify after transformation

summary(as.numeric(expr3189))

max(expr3189)

boxplot(expr3189,
        outline=FALSE,
        las=2,
        main="GSE3189")


#2.11 Save the processed expression matrix

save(
  expr3189,
  file="RData/expr3189_processed.RData"
)

write.csv(
  expr3189,
  "Results/expr3189_processed.csv"
)

max(expr3189)

summary(as.numeric(expr3189))

save(expr3189,
     file = "RData/expr3189_processed.RData")


png("Figures/GSE3189_Boxplot.png",
    width = 1800,
    height = 1200,
    res = 300)

boxplot(expr3189,
        outline = FALSE,
        las = 2,
        main = "GSE3189")

dev.off()

graphics.off()

png("Figures/GSE3189_Boxplot.png",
    width = 1800,
    height = 1200,
    res = 300)

boxplot(expr3189,
        outline = FALSE,
        las = 2,
        main = "GSE3189")

dev.off()



dim(GSE46517)

head(GSE46517[,1])

dim(GSE7553)

head(GSE7553[,1])



############################################################
## GSE46517 PREPROCESSING
############################################################

# Extract probe IDs

probe46517 <- GSE46517[,1]


# Annotate probes

gene46517 <- mapIds(
  hgu133a.db,
  keys = probe46517,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)


# Add gene symbols

GSE46517$GeneSymbol <- gene46517


# Remove probes without gene symbols

GSE46517 <- GSE46517[!is.na(GSE46517$GeneSymbol), ]


# Create expression matrix

expr46517 <- as.matrix(
  GSE46517[, !(names(GSE46517) %in% c("V1","GeneSymbol"))]
)

rownames(expr46517) <- GSE46517$GeneSymbol

mode(expr46517) <- "numeric"


# Collapse duplicate probes

expr46517 <- avereps(
  expr46517,
  ID = rownames(expr46517)
)


# Check whether log2 transformation is required

summary(as.numeric(expr46517))

quantile(as.numeric(expr46517),
         probs=c(0,0.25,0.5,0.75,0.99,1))

max(expr46517)


# Log2 transform only if needed

if(max(expr46517) > 100){
  expr46517 <- log2(expr46517 + 1)
}


# Verify

summary(as.numeric(expr46517))

max(expr46517)


# Save

save(expr46517,
     file="RData/expr46517_processed.RData")

write.csv(expr46517,
          "Results/expr46517_processed.csv")


#Preprocess GSE7553

############################################################
## GSE7553 PREPROCESSING
############################################################

# Extract probe IDs

probe7553 <- GSE7553[,1]


# Annotate probes

gene7553 <- mapIds(
  hgu133a.db,
  keys = probe7553,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)


# Add gene symbols

GSE7553$GeneSymbol <- gene7553


# Remove probes without gene symbols

GSE7553 <- GSE7553[!is.na(GSE7553$GeneSymbol), ]


# Create expression matrix

expr7553 <- as.matrix(
  GSE7553[, !(names(GSE7553) %in% c("V1","GeneSymbol"))]
)

rownames(expr7553) <- GSE7553$GeneSymbol

mode(expr7553) <- "numeric"


# Collapse duplicate probes

expr7553 <- avereps(
  expr7553,
  ID = rownames(expr7553)
)


# Check whether log2 transformation is required

summary(as.numeric(expr7553))

quantile(as.numeric(expr7553),
         probs=c(0,0.25,0.5,0.75,0.99,1))

max(expr7553)


# Log2 transform only if needed

if(max(expr7553) > 100){
  expr7553 <- log2(expr7553 + 1)
}


# Verify

summary(as.numeric(expr7553))

max(expr7553)


# Save

save(expr7553,
     file="RData/expr7553_processed.RData")

write.csv(expr7553,
          "Results/expr7553_processed.csv")


# After running both scripts, please paste the outputs of:

max(expr46517)

summary(as.numeric(expr46517))

max(expr7553)

summary(as.numeric(expr7553))

save.image("RData/Step03_AllDatasets_Preprocessed.RData")


colnames(expr3189)

colnames(expr46517)

colnames(expr7553)


head(colnames(expr3189), 20)

head(colnames(expr46517), 20)

head(colnames(expr7553), 20)



#GSE3189 (Normal vs Melanoma)

############################################################
## GSE3189
############################################################

group3189 <- ifelse(
  grepl("^Normal", colnames(expr3189)),
  "Normal",
  ifelse(
    grepl("^Melanoma", colnames(expr3189)),
    "Melanoma",
    "Remove"
  )
)

table(group3189)

keep3189 <- group3189 != "Remove"

expr3189.DEG <- expr3189[, keep3189]

group3189.DEG <- factor(
  group3189[keep3189],
  levels=c("Normal","Melanoma")
)

table(group3189.DEG)


#GSE46517 (Primary + Metastatic vs Normal)

############################################################
## GSE46517
############################################################

group46517 <- ifelse(
  grepl("^Primary Melanoma", colnames(expr46517)),
  "Melanoma",
  ifelse(
    grepl("^Metastatic Melanoma", colnames(expr46517)),
    "Melanoma",
    ifelse(
      grepl("^Normal Skin", colnames(expr46517)),
      "Normal",
      ifelse(
        grepl(
          "^Normal Epithelial Melanocytes",
          colnames(expr46517)
        ),
        "Normal",
        "Remove"
      )
    )
  )
)

table(group46517)

keep46517 <- group46517 != "Remove"

expr46517.DEG <- expr46517[, keep46517]

group46517.DEG <- factor(
  group46517[keep46517],
  levels=c("Normal","Melanoma")
)

table(group46517.DEG)


#GSE7553 (Primary + Metastatic vs Normal)

############################################################
## GSE7553
############################################################

group7553 <- ifelse(
  grepl("^Primary Melanoma", colnames(expr7553)),
  "Melanoma",
  ifelse(
    grepl("^Metastatic Melanoma", colnames(expr7553)),
    "Melanoma",
    ifelse(
      grepl("^Normal Skin", colnames(expr7553)),
      "Normal",
      ifelse(
        grepl(
          "^normal human epidermal melanocytes",
          colnames(expr7553),
          ignore.case=TRUE
        ),
        "Normal",
        "Remove"
      )
    )
  )
)

table(group7553)

keep7553 <- group7553 != "Remove"

expr7553.DEG <- expr7553[, keep7553]

group7553.DEG <- factor(
  group7553[keep7553],
  levels=c("Normal","Melanoma")
)

table(group7553.DEG)


save(
  expr3189.DEG,
  expr46517.DEG,
  expr7553.DEG,
  group3189.DEG,
  group46517.DEG,
  group7553.DEG,
  file="RData/Filtered_DEG_Datasets.RData"
)


table(group3189.DEG)

table(group46517.DEG)

table(group7553.DEG)



#DEG for GSE3189

############################################################
## DEG Analysis - GSE3189
############################################################

design3189 <- model.matrix(
  ~0 + group3189.DEG
)

colnames(design3189) <- levels(group3189.DEG)

design3189

contrast3189 <- makeContrasts(
  Melanoma - Normal,
  levels = design3189
)

fit3189 <- lmFit(
  expr3189.DEG,
  design3189
)

fit3189 <- contrasts.fit(
  fit3189,
  contrast3189
)

fit3189 <- eBayes(
  fit3189
)

deg3189 <- topTable(
  fit3189,
  number = Inf,
  adjust.method = "BH"
)

head(deg3189)

write.csv(
  deg3189,
  "Results/GSE3189_DEG.csv"
)

save(
  deg3189,
  file = "RData/GSE3189_DEG.RData"
)


#DEG for GSE46517

############################################################
## DEG Analysis - GSE46517
############################################################

design46517 <- model.matrix(
  ~0 + group46517.DEG
)

colnames(design46517) <- levels(group46517.DEG)

contrast46517 <- makeContrasts(
  Melanoma - Normal,
  levels = design46517
)

fit46517 <- lmFit(
  expr46517.DEG,
  design46517
)

fit46517 <- contrasts.fit(
  fit46517,
  contrast46517
)

fit46517 <- eBayes(
  fit46517
)

deg46517 <- topTable(
  fit46517,
  number = Inf,
  adjust.method = "BH"
)

write.csv(
  deg46517,
  "Results/GSE46517_DEG.csv"
)

save(
  deg46517,
  file = "RData/GSE46517_DEG.RData"
)


#DEG for GSE7553

############################################################
## DEG Analysis - GSE7553
############################################################

design7553 <- model.matrix(
  ~0 + group7553.DEG
)

colnames(design7553) <- levels(group7553.DEG)

contrast7553 <- makeContrasts(
  Melanoma - Normal,
  levels = design7553
)

fit7553 <- lmFit(
  expr7553.DEG,
  design7553
)

fit7553 <- contrasts.fit(
  fit7553,
  contrast7553
)

fit7553 <- eBayes(
  fit7553
)

deg7553 <- topTable(
  fit7553,
  number = Inf,
  adjust.method = "BH"
)

write.csv(
  deg7553,
  "Results/GSE7553_DEG.csv"
)

save(
  deg7553,
  file = "RData/GSE7553_DEG.RData"
)


#Count significant DEGs

#After all three analyses finish, run:

sum(
  deg3189$adj.P.Val < 0.05 &
    abs(deg3189$logFC) >= 1
)

sum(
  deg46517$adj.P.Val < 0.05 &
    abs(deg46517$logFC) >= 1
)

sum(
  deg7553$adj.P.Val < 0.05 &
    abs(deg7553$logFC) >= 1
)


#Next: Save only significant DEGs

#Run:

DEG3189 <- subset(
  deg3189,
  adj.P.Val < 0.05 &
    abs(logFC) >= 1
)

DEG46517 <- subset(
  deg46517,
  adj.P.Val < 0.05 &
    abs(logFC) >= 1
)

DEG7553 <- subset(
  deg7553,
  adj.P.Val < 0.05 &
    abs(logFC) >= 1
)

write.csv(
  DEG3189,
  "Results/GSE3189_Significant_DEGs.csv"
)

write.csv(
  DEG46517,
  "Results/GSE46517_Significant_DEGs.csv"
)

write.csv(
  DEG7553,
  "Results/GSE7553_Significant_DEGs.csv"
)


#Get common DEGs across all three datasets

genes3189 <- rownames(DEG3189)

genes46517 <- rownames(DEG46517)

genes7553 <- rownames(DEG7553)

common_DEGs <- Reduce(
  intersect,
  list(
    genes3189,
    genes46517,
    genes7553
  )
)

length(common_DEGs)

write.csv(
  common_DEGs,
  "Results/Common_DEGs_Three_Datasets.csv",
  row.names = FALSE
  savehistory("Scripts/02_DEG_Analysis.R")
)

savehistory("Scripts/02_DEG_Analysis.R")
file.exists("Scripts/02_DEG_Analysis.R")

