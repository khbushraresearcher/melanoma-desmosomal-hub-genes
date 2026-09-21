###############################################################
setwd("~/Desktop/Melanoma_Project2")
###############################################################
## Load libraries
###############################################################
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
###############################################################
## Create output folders
###############################################################
dir.create("Results", showWarnings = FALSE)
dir.create("Results/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Tables", recursive = TRUE, showWarnings = FALSE)
dir.create("RData", showWarnings = FALSE)
###############################################################
## Load workspace from Script 02
###############################################################
load("RData/Script2_DEG_Workspace.RData")
###############################################################
## Expression matrix
###############################################################
datExpr <- t(expr7553.DEG)
###############################################################
## Quality control
###############################################################
gsg <- goodSamplesGenes(datExpr)
print(gsg)
datExpr <- datExpr[
gsg$goodSamples,
gsg$goodGenes
]
print(dim(datExpr))
###############################################################
## Sample clustering
###############################################################
sampleTree <- hclust(
dist(datExpr),
method = "average"
)
png(
"Results/Figures/GSE7553_Figure1_SampleClustering.png",
width = 2200,
height = 1800,
res = 300
)
plot(
sampleTree,
main = "Sample Clustering",
sub = "",
xlab = ""
)
dev.off()
###############################################################
## Trait data
###############################################################
traitData <- data.frame(
Melanoma = ifelse(
group7553.DEG == "Melanoma",
1,
0
)
)
rownames(traitData) <- rownames(datExpr)
print(head(traitData))
###############################################################
## Soft-threshold selection
###############################################################
powers <- 1:20
sft <- pickSoftThreshold(
datExpr,
powerVector = powers,
networkType = "signed",
verbose = 5
)
###############################################################
## Soft-threshold plots
###############################################################
png(
"Results/Figures/GSE7553_Figure2_SoftThreshold.png",
width = 3200,
height = 1600,
res = 300
)
par(mfrow = c(1,2))
plot(
sft$fitIndices[,1],
-sign(sft$fitIndices[,3]) *
sft$fitIndices[,2],
type = "n",
xlab = "Soft Threshold (Power)",
ylab = "Scale-Free Topology Fit (Signed R²)",
main = "Scale Independence"
)
text(
sft$fitIndices[,1],
-sign(sft$fitIndices[,3]) *
sft$fitIndices[,2],
labels = powers,
col = "red",
cex = 1
)
abline(
h = 0.90,
col = "blue",
lwd = 2
)
plot(
sft$fitIndices[,1],
sft$fitIndices[,5],
type = "n",
xlab = "Soft Threshold (Power)",
ylab = "Mean Connectivity",
main = "Mean Connectivity"
)
text(
sft$fitIndices[,1],
sft$fitIndices[,5],
labels = powers,
col = "red",
cex = 1
)
dev.off()
###############################################################
## Construct signed co-expression network
###############################################################
## Original analysis used power = 15
softPower <- 15
net <- blockwiseModules(
datExpr,
power = softPower,
networkType = "signed",
TOMType = "signed",
maxBlockSize = 20000,
deepSplit = 2,
minModuleSize = 30,
mergeCutHeight = 0.25,
reassignThreshold = 0,
pamRespectsDendro = FALSE,
numericLabels = TRUE,
saveTOMs = FALSE,
verbose = 3
)
###############################################################
## Module colors
###############################################################
moduleColors <- labels2colors(
net$colors
)
print(table(moduleColors))
print(
sort(
table(moduleColors),
decreasing = TRUE
)
)
###############################################################
## Check number of blocks
###############################################################
print(
length(net$dendrograms)
)
print(
sapply(
net$blockGenes,
length
)
)
###############################################################
## Module eigengenes
###############################################################
MEs <- moduleEigengenes(
datExpr,
colors = moduleColors
)$eigengenes
MEs <- orderMEs(MEs)
###############################################################
## Gene dendrograms
###############################################################
for(i in seq_along(net$dendrograms)){
png(
paste0(
"Results/Figures/GSE7553_Dendrogram_Block",
i,
".png"
),
width = 4500,
height = 2500,
res = 300
)
plotDendroAndColors(
net$dendrograms[[i]],
moduleColors[
net$blockGenes[[i]]
],
"Module Colors",
dendroLabels = FALSE,
hang = 0.03,
addGuide = TRUE,
guideHang = 0.05
)
dev.off()
}
###############################################################
## Module-Trait Relationships
###############################################################
moduleTraitCor <- cor(
MEs,
traitData,
use = "p"
)
moduleTraitPvalue <- corPvalueStudent(
moduleTraitCor,
nrow(datExpr)
)
print(moduleTraitCor)
print(moduleTraitPvalue)
###############################################################
## Module-Trait Heatmap
###############################################################
textMatrix <- paste(
signif(moduleTraitCor, 2),
"\n(",
signif(moduleTraitPvalue, 1),
")",
sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)
png(
"Results/Figures/GSE7553_ModuleTraitHeatmap.png",
width = 2200,
height = 3200,
res = 300
)
par(
mar = c(8,10,4,3)
)
labeledHeatmap(
Matrix = moduleTraitCor,
xLabels = "Melanoma",
yLabels = rownames(moduleTraitCor),
ySymbols = rownames(moduleTraitCor),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.9,
zlim = c(-1,1),
main = "Module-Trait Relationships"
)
dev.off()
###############################################################
## Gene Significance
###############################################################
GS <- as.data.frame(
cor(
datExpr,
traitData$Melanoma,
use = "p"
)
)
colnames(GS) <- "GS"
###############################################################
## Module Membership
###############################################################
MM <- as.data.frame(
cor(
datExpr,
MEs,
use = "p"
)
)
colnames(MM) <- paste0(
"MM.",
colnames(MM)
)
###############################################################
## Hub Gene Table
###############################################################
hubTable7553 <- data.frame(
Gene = colnames(datExpr),
Module = moduleColors,
GS = GS$GS,
MM,
stringsAsFactors = FALSE
)
write.csv(
hubTable7553,
"Results/Tables/HubTable_GSE7553.csv",
row.names = FALSE
)
###############################################################
## Check candidate genes
###############################################################
candidateGenes <- c(
"FLG",
"DSG3",
"DSG1"
)
candidateHubTable <- subset(
hubTable7553,
Gene %in% candidateGenes
)
print(candidateHubTable)
write.csv(
candidateHubTable,
"Results/Tables/GSE7553_Candidate_Hub_Genes.csv",
row.names = FALSE
)
###############################################################
## GS-MM Correlation
###############################################################
modules <- setdiff(
unique(moduleColors),
"grey"
)
GSMM7553 <- data.frame()
for(module in modules){
moduleGenes <- moduleColors == module
MMcolumn <- paste0(
"MM.ME",
module
)
ct <- cor.test(
MM[
moduleGenes,
MMcolumn
],
GS[
moduleGenes,
"GS"
]
)
GSMM7553 <- rbind(
GSMM7553,
data.frame(
Module = module,
Correlation = unname(
ct$estimate
),
Pvalue = ct$p.value,
Genes = sum(moduleGenes)
)
)
}
GSMM7553 <- GSMM7553[
order(
-abs(
GSMM7553$Correlation
)
),
]
print(GSMM7553)
write.csv(
GSMM7553,
"Results/Tables/GSMM_Correlation_GSE7553.csv",
row.names = FALSE
)
###############################################################
## GS vs MM plots
###############################################################
for(module in modules){
moduleGenes <- moduleColors == module
MMcolumn <- paste0(
"MM.ME",
module
)
r <- cor(
MM[
moduleGenes,
MMcolumn
],
GS[
moduleGenes,
"GS"
],
use = "p"
)
p <- cor.test(
MM[
moduleGenes,
MMcolumn
],
GS[
moduleGenes,
"GS"
]
)$p.value
png(
paste0(
"Results/Figures/GSE7553_GSvsMM_",
module,
".png"
),
width = 2200,
height = 2200,
res = 300
)
verboseScatterplot(
MM[
moduleGenes,
MMcolumn
],
GS[
moduleGenes,
"GS"
],
xlab = paste(
"Module Membership (",
module,
")",
sep = ""
),
ylab = "Gene Significance",
main = paste0(
tools::toTitleCase(module),
" Module\nr = ",
round(r, 3),
", P = ",
signif(p, 3)
),
col = module,
pch = 19,
cex = 0.8
)
dev.off()
}
###############################################################
## Blue Module Summary
###############################################################
DEG7553 <- subset(
deg7553,
adj.P.Val < 0.05 &
abs(logFC) >= 1
)
###############################################################
## Genes in blue module
###############################################################
blueGenes <- hubTable7553$Gene[
hubTable7553$Module == "blue"
]
###############################################################
## DEGs in blue module
###############################################################
blueDEGs <- DEG7553[
rownames(DEG7553) %in% blueGenes,
]
###############################################################
## Blue module summary
###############################################################
BlueSummary7553 <- data.frame(
Dataset = "GSE7553",
Module = "blue",
Trait = "Melanoma",
Correlation = round(
moduleTraitCor[
"MEblue",
"Melanoma"
],
3
),
p.value = signif(
moduleTraitPvalue[
"MEblue",
"Melanoma"
],
3
),
Genes_in_Module = length(
blueGenes
),
DEGs_in_Module = nrow(
blueDEGs
),
Upregulated = sum(
blueDEGs$logFC > 1
),
Downregulated = sum(
blueDEGs$logFC < -1
)
)
print(BlueSummary7553)
write.csv(
BlueSummary7553,
"Results/Tables/GSE7553_Blue_Module_Summary.csv",
row.names = FALSE
)
###############################################################
## Save WGCNA workspace
###############################################################
save(
net,
moduleColors,
MEs,
GS,
MM,
hubTable7553,
candidateHubTable,
GSMM7553,
moduleTraitCor,
moduleTraitPvalue,
datExpr,
traitData,
BlueSummary7553,
file = "RData/GSE7553_WGCNA.RData"
)
###############################################################
## Save complete Script 05 console history as .R
###############################################################
savehistory(
"Scripts/05_WGCNA_GSE7553.R"
)
