# Single cell RNAseq analysis
# Reva Shenwai
# 23-Feb-2026
# ---------------------
# Lecture: https://satijalab.org/seurat/articles/hashing_vignette

# libraries
library(data.table)
library(tidyverse)
library(Seurat)
library(glmGamPoi, warn.conflicts = FALSE) # BiocManager::install("glmGamPoi")  Bioconuctor package for SCTransform.  
library(harmony)

library(patchwork)
library(ggplot2)
library(dittoSeq) # # BiocManager::install("dittoSeq") 
library(future)
library(scales)
library(ggrastr) # To make raster featureplot.  CRAN 

library(harmony)  # CRAN.   install.packages("harmony") to use RunHarmony() https://github.com/immunogenomics/harmony

##################################################################################################
# settings
# setwd("~/RStudio_files/PHF19/scRNAseq_Eileen")
dir <- dirname(rstudioapi::getSourceEditorContext()$path); print(dir) # [1] "/Users/lees130/Library/CloudStorage/OneDrive-NYULangoneHealth/N19_PHF19_HashTagSeq_Hussein/05_SeuratObj_Plasmablast"
setwd(dir)
create_UMAP <- TRUE

LibraryFileDir <- "/Users/lees130/Library/CloudStorage/OneDrive-NYULangoneHealth/N19_PHF19_HashTagSeq_Hussein/04_Data_Invitro_iPC_PHF19"
MouseGeneInfoFile <- "/Users/lees130/Library/CloudStorage/OneDrive-NYULangoneHealth/N04_CodingStudy/00b_ReferenceGenomeGeneInfo_MusMusculus/Mus_musculus.gene_info_20260515.txt"
##################################################################################################

# load data
# if(create_UMAP){
######### library 1
## paste0(LibraryFileDir, "/Library1/") has the "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz" files of library 1. 
# lib1.data <- Read10X(data.dir=paste0(getwd(),"/data/Invitro_iPC-PHF19/Library1/")) 
lib1.data <- Read10X(data.dir=paste0(LibraryFileDir, "/Library1/")); class(lib1.data); # list   # this takes 3 min. 
names(lib1.data) # ] "Gene Expression"  "Antibody Capture"
class(lib1.data[["Antibody Capture"]]); dim(lib1.data[["Antibody Capture"]]) # 9 24231
rownames(lib1.data[["Antibody Capture"]])
# [1] "TotalSeq-B0301anti-mouseHashtag1" "TotalSeq-B0302anti-mouseHashtag2" "TotalSeq-B0303anti-mouseHashtag3" "TotalSeq-B0304anti-mouseHashtag4" "TotalSeq-B0305anti-mouseHashtag5"
# [6] "TotalSeq-B0306anti-mouseHashtag6" "TotalSeq-B0307anti-mouseHashtag7" "TotalSeq-B0308anti-mouseHashtag8" "TotalSeq-B0309anti-mouseHashtag9"

# rename oligos
custom_names <- c("Control1", "Control2", "Control3",     "Long4", "Long5", "Long6",       "Mutant7", "Mutant8", "Mutant9")

# rename rows of Antibody Capture matrix
rownames(lib1.data[["Antibody Capture"]]) <- custom_names

# create seurat obj
lib1 <- CreateSeuratObject(counts=lib1.data$`Gene Expression`,  project="Lib1")

# add HTO matrix (hashtag oligos) as independent assay
lib1[["HTO"]] <- CreateAssayObject(counts=lib1.data$`Antibody Capture`)
names(lib1) # "RNA" "HTO"

saveRDS(lib1, "lib1.rds")

########### ++++++++++++ library 2 +++++++++++++++  ############
## paste0(LibraryFileDir, "/Library2/") has the "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz" files of library2. 

# lib2.data <- Read10X(data.dir=paste0(getwd(),  "/data/Invitro_iPC-PHF19/Library2/"))  
lib2.data <- Read10X(data.dir=paste0(LibraryFileDir, "/Library2/")); class(lib2.data)

# rename rows of Antibody Capture matrix
rownames(lib2.data[["Antibody Capture"]]) <- custom_names

# create seurat obj
lib2 <- CreateSeuratObject(counts=lib2.data$`Gene Expression`, project="Lib2")
# add HTO matrix (hashtag oligos) as independent assay
lib2[["HTO"]] <- CreateAssayObject(counts=lib2.data$`Antibody Capture`)
names(lib2) # [1] "RNA" "HTO"

saveRDS(lib2, "lib2.rds")

################################################################################################################
############## Reva's original code:  merge into one object and then demultiplexed.   ##############   ############## 
# phf19.seurat <- merge(lib1, y=lib2, add.cell.ids=c("Lib1","Lib2"))
# Assays(phf19.seurat) # "RNA" "HTO"
# 
# head(colnames(phf19.seurat))
# table(phf19.seurat$orig.ident)
# 
# # clean
# rm(lib1,lib2,lib1.data,lib2.data)
# 
# # Switch to HTO assay
# DefaultAssay(phf19.seurat) <- "HTO"
# # Normalize HTO counts
# phf19.seurat <- NormalizeData(phf19.seurat, assay="HTO",  normalization.method="CLR")
# 
# # Demultiplex
# phf19.seurat <- HTODemux(phf19.seurat,  assay="HTO",  positive.quantile=0.99)
# 
# # Check tbl for the 9 samples
# table(phf19.seurat$hash.ID)
# # Doublet    Long5    Long6  Mutant9  Mutant8  Mutant7 Control1    Long4 Control2 Negative Control3 
# # 8404     4853     4259     3933     3705     3931     3942     4307     3679     2576     3937 
##############   ##############   ##############   ##############   ##############   ##############   ############## 

################################################################################################################
### Step1. Demultiplex first and then merge . 
################################################################################################################
### Demultiplex 
lib1 <- HTODemux(lib1, assay="HTO")
lib2 <- HTODemux(lib2, assay="HTO")

## merge 
phf19.seurat <- merge(lib1, y=lib2, add.cell.ids=c("Lib1","Lib2"))

Assays(phf19.seurat) # "RNA" "HTO"

head(colnames(phf19.seurat))  # "Lib1_AAACCCAAGGATACGC-1" "Lib1_AAACCCAAGGGACCAT-1" "Lib1_AAACCCAAGTGTACAA-1"
table(phf19.seurat$orig.ident)
# Lib1  Lib2 
# 24231 23295 

# Switch to HTO assay
DefaultAssay(phf19.seurat) <- "HTO"
# Normalize HTO counts
phf19.seurat <- NormalizeData(phf19.seurat, assay="HTO",  normalization.method="CLR")

# Check tbl for the 9 samples
table(phf19.seurat$hash.ID)
# Control1 Control2 Control3  Doublet    Long4    Long5    Long6  Mutant7  Mutant8  Mutant9 Negative 
# 4035     3975     3993     7535     2641     3045     4379     3998     3756     3972     6197 

#saveRDS(phf19.seurat, "phf19.seurat_Merged.rds")
phf19.seurat <- readRDS("phf19.seurat_Merged.rds")
################################################################################################################

################################################################################################################
### Step2. Keep signlets and visualization
################################################################################################################
# Clean data: keep only 'Singlets'
phf19.seurat <- subset(phf19.seurat,    idents="Negative", invert=TRUE)
phf19.seurat <- subset(phf19.seurat,    idents="Doublet",  invert=TRUE)

# Check tbl for the 9 samples
table(phf19.seurat$hash.ID)
# Control1 Control2 Control3    Long4    Long5    Long6  Mutant7  Mutant8  Mutant9 
# 4035     3975     3993     2641     3045     4379     3998     3756     3972 

# visualize enrichment for HTOs with ridge plots
# Group cells based on max HTO signal
Idents(phf19.seurat) <- "HTO_maxID"
# order
Idents(phf19.seurat) <- factor(Idents(phf19.seurat),
                                levels=c("Control1",  "Control2",  "Control3",  "Long4", "Long5", "Long6",   "Mutant7","Mutant8","Mutant9"))
phf19.seurat$hash.ID <- factor(phf19.seurat$hash.ID,
                               levels=c("Control1", "Control2", "Control3",  "Long4", "Long5","Long6",  "Mutant7", "Mutant8", "Mutant9"))

# write PNG plot
# jpeg(filename = paste0(getwd(),"/results/demultiplex.samples_ridge.plots.jpg"),  width = 1400,   height = 1000,    units = "px")
#pdf(filename = "Fig1A_demultiplex.samples_ridge.plots.pdf",  width = 10,   height = 10)
# plot
MyRidgePlot <- RidgePlot(phf19.seurat, assay="HTO",   features=rownames(phf19.seurat[["HTO"]])[9:1],  ncol=3)
# MyRidgePlot <- RidgePlot(phf19.seurat, assay="HTO",   features = rev(rownames(phf19.seurat[["HTO"]])),  ncol=3)
# save
# dev.off()
ggsave("Fig1A_demultiplex.samples_ridge.plots.pdf", plot = MyRidgePlot, width=10, height=10)


# Compare number of UMIs for singlets, doublets and negative cells
Idents(phf19.seurat) <- "HTO_classification.global"
VlnPlot(phf19.seurat, features="nCount_RNA", pt.size=0.1, log=TRUE)   #  It shows singlet cells. 

## Process RNA data
# Switch back to Gene Expression
DefaultAssay(phf19.seurat) <- "RNA"

# Add group info
sample_to_group <- c(
  "Control1" = "Control", "Control2" = "Control",   "Control3" = "Control",
  "Long4" = "Long",  "Long5" = "Long",  "Long6" = "Long",  
  "Mutant7" = "Mutant",    "Mutant8" = "Mutant",    "Mutant9" = "Mutant"  )

# map hash.ID to group name
group_vector <- sample_to_group[phf19.seurat$hash.ID]

# assign barcodes to groups
names(group_vector) <- colnames(phf19.seurat)

# add group to metadata
phf19.seurat <- AddMetaData(phf19.seurat, metadata=group_vector, col.name="group")

# Calculate mitochondrial %
# Check first:
 grep("^mt-|^MT-", rownames(phf19.seurat), value = TRUE)[1:10] # [1] "mt-Nd1"  "mt-Nd2"  "mt-Co1"  "mt-Co2"  "mt-Atp8" "mt-Atp6" "mt-Co3"  "mt-Nd3"  "mt-Nd4l" "mt-Nd4" 
phf19.seurat[["percent.mt"]] <- PercentageFeatureSet(phf19.seurat, pattern="^mt-")

# write PNG plot
#jpeg(filename = paste0(getwd(),"/results/biological.qc.metrics_vln.plot_pre.filter.jpg"),   width = 800,  height = 400,   units = "px")
# visualize key biological QC metrics
QCMetricsPlot<-VlnPlot(phf19.seurat,   features=c("nFeature_RNA","nCount_RNA","percent.mt"), group.by="hash.ID",   ncol=3,   pt.size=0.1)
# save
# dev.off()
ggsave("Fig1B_QCMetricsPlot_nFeature_nCount_percentmt.pdf", plot = QCMetricsPlot, width=10, height=6)

# Check correlation (RNA)
plot1 <- FeatureScatter(phf19.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "hash.ID")
plot2 <- FeatureScatter(phf19.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "hash.ID")
# write plot
# jpeg(filename = paste0(getwd(),"/results/featureplot_pre.filter.jpg"), width = 900,  height = 500,   units = "px")
(plot1 + plot2) + plot_layout(guides = "collect") & theme(legend.position = "right")
# save
# dev.off()

# Apply filters: Filter low-quality cells using nFeature_RNA, nCount_RNA, and percent.mt.
phf19.seurat <- subset(phf19.seurat,  subset = nFeature_RNA > 200 &  nFeature_RNA < 6000 &  percent.mt < 10)

# write PNG plot
# jpeg(filename = paste0(getwd(),"/results/biological.qc.metrics_vln.plot_post.filter.jpg"),  width = 800,  height = 400,   units = "px")
# visualize after filtering
VlnPlot(phf19.seurat, features=c("nFeature_RNA","nCount_RNA","percent.mt"), group.by="hash.ID",  ncol=3,   pt.size=0.1)
# dev.off()

#  Check correlation post-filtering
plot1 <- FeatureScatter(phf19.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "hash.ID")
plot2 <- FeatureScatter(phf19.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "hash.ID")

# jpeg(filename = paste0(getwd(),"/results/featureplot_post.filter.jpg"),  width = 900, height = 500, units = "px") # High resolution for crisp text
(plot1 + plot2) + plot_layout(guides = "collect") & theme(legend.position = "right")
# dev.off()

table(phf19.seurat$hash.ID) 
# Control1 Control2 Control3    Long4    Long5    Long6  Mutant7  Mutant8  Mutant9 
# 3465     3482     3521     2291     2602     3640     3405     3315     3473 
# # clean
# rm(plot1,plot2)

saveRDS(phf19.seurat, "phf19.seurat_MergedQualityControl.rds")

########################################################################################################################################
### 20260701
### Step3. remove immunoglobulin variable genes, cell-cycle genes, mitochondrial genes, ribosomal genes, or stress-related genes 
###                 from the variable features before PCA and clustering.
###     ## I did not permanently remove these genes from Seurat Object. They can still be useful for quality control, cell annotation, or biological interpretation
########################################################################################################################################
## 1) Mitochondrial genes
mito.genes <- grep("^mt-|^MT-", rownames(phf19.seurat), value=TRUE)

## 2) Ribosomal genes
ribo.genes <- grep("^Rpl|^Rps", rownames(phf19.seurat), value=TRUE)

## 3) Ribosomal genes
hla.genes <- grep("^H2-", rownames(phf19.seurat), value=TRUE); length(hla.genes); hla.genes[1:5] # "H2-K1"   "H2-Ke6"  "H2-Oa"   "H2-DMa"  "H2-DMb2"

## 4) Immunogloblulin gene variable segments 
ig.genes <- grep("^Igh|^Igk|^Igl", rownames(phf19.seurat), value = TRUE) 
## If I want to remove TCR (variable segments)
# ig.genes <- grep("^Igh|^Igk|^Igl|^Trav|^Trbv|^Trgv|^Trdv", rownames(phf19.seurat), value = TRUE)

## 5) Cell cycle genes
cc.genes <- unique(c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes ))
## Convert human genes to mouse genes. 
cc.genes <- intersect(rownames(phf19.seurat), tools::toTitleCase(tolower(cc.genes)) )
length(cc.genes); (cc.genes)[1:5] # "Hjurp" "Mcm6"  "Nuf2"  "Exo1"  "Lbr"  

## 6) stress-response genes 
stress.genes <- c("Fos","Fosb","Jun","Junb","Jund",  "Atf3","Egr1","Egr2","Egr3",  "Ddit3","Hspa1a","Hspa1b","Hsp90aa1",
                  "Hspb1","Dnajb1","Xbp1"); length(stress.genes) # 16
stress.genes <- intersect(rownames(phf19.seurat), stress.genes); length(stress.genes) # 16

#### +++ Combine all genes to exclude  ++++ ###### 
exclude.genes <- unique(c(mito.genes, ribo.genes, hla.genes,ig.genes, cc.genes, stress.genes ))
length(exclude.genes); head(exclude.genes) # 656 #[1] "mt-Nd1"  "mt-Nd2"  "mt-Co1"  "mt-Co2"  "mt-Atp8" "mt-Atp6"

########################################################################################################################################
### Step4. Normalization with SCTransform, run PCA and Harmony   ## If I use SCTransform, I don't need to run ScaleData() and NormalizeData()
## SCTransform, it already does: normalization,  variance, stabilization,  scaling (internally), regression (e.g., percent.mt)
########################################################################################################################################
## SCTransform. ### Memory issue. I should go to bigpurple
# DefaultAssay(phf19.seurat) <- "RNA"
# phf19.seurat_SCT <- SCTransform(phf19.seurat, method="glmGamPoi", vars.to.regress = "percent.mt",    verbose = FALSE)   
# phf19.seurat_SCT_3000 <- SCTransform(phf19.seurat, method="glmGamPoi", vars.to.regress="percent.mt",  
#                                verbose=FALSE, variable.features.n=3000) ### default is 2000. Let's start with 3000 because I exclude 656 genes

#saveRDS(phf19.seurat_SCT_3000, phf19.seurat_MergedQualityControl_SCT_3000g.rds)
phf19.seurat_SCT_3000g <- readRDS("phf19.seurat_MergedQualityControl_SCT_3000g.rds")

#### Before removing exclude.gene
phf19.seurat_SCT_3000g; 
# An object of class Seurat 
# 53368 features across 33794 samples within 3 assays 
# Active assay: SCT (21074 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 2 other assays present: RNA, HTO

#### Remove exclude.genes 
VariableFeatures(phf19.seurat_SCT_3000g) <- setdiff(VariableFeatures(phf19.seurat_SCT_3000g), exclude.genes)
phf19.seurat_SCT_3000g; length(VariableFeatures(phf19.seurat_SCT_3000g)) # 2627 
# An object of class Seurat 
# 53368 features across 33794 samples within 3 assays 
# Active assay: SCT (21074 features, 2627 variable features)
# 3 layers present: counts, data, scale.data
# 2 other assays present: RNA, HTO

##################    ##################    ##################    ##################    ##################    ##################  
##################  Without Harmoney, Normalization and PCA  ################
##################    ##################    ##################    ##################    ##################    ##################  

phf19.seurat_NoHarmony <- phf19.seurat
# Standard RNA workflow: normalize & find top diffexp genes
phf19.seurat_NoHarmony <- NormalizeData(phf19.seurat_NoHarmony,   normalization.method="LogNormalize",  scale.factor=10000)
phf19.seurat_NoHarmony <- FindVariableFeatures(phf19.seurat_NoHarmony, selection.method="vst",  nfeatures=2000)
# visualize top 10 most variable genes
top10 <- head(VariableFeatures(phf19.seurat_NoHarmony), 10)
# save plot
jpeg(filename = paste0(getwd(),"/results/variableFeaturePlot.top10.jpg"),    width = 900,     height = 500,   units = "px")
LabelPoints(plot = VariableFeaturePlot(phf19.seurat_NoHarmony), points = top10, repel = TRUE)
dev.off()
# scale
all.genes <- rownames(phf19.seurat_NoHarmony)
phf19.seurat_NoHarmony <- ScaleData(phf19.seurat_NoHarmony, features=all.genes)
# clean data
# rm(all.genes)

# Run PCA
phf19.seurat_NoHarmony <- RunPCA(phf19.seurat_NoHarmony,  features=VariableFeatures(object=phf19.seurat_NoHarmony))   ### This takes 10 min
saveRDS(phf19.seurat_NoHarmony, "phf19.seurat_NoHarmony_PCA.rds")

  # save plot
jpeg(filename = paste0(getwd(),"/results/topgenes_pc1_pc2.jpg"), width = 800, height = 600, units = "px")
# Visualize which genes are driving the first two PCs
VizDimLoadings(phf19.seurat_NoHarmony, dims=1:2, reduction="pca")
dev.off()

# save plot
jpeg(filename = paste0(getwd(),"/results/elbowplot.jpg"),   width = 800,    height = 350,     units = "px")
# elbow plot to decide # dimensions
ElbowPlot(phf19.seurat_NoHarmony, ndims=30)
dev.off()

# PCA visualization (technical batch: lib 1, 2)
# jpeg(filename = paste0(getwd(),"/results/pca_libraries.jpg"), width = 550, height = 400,  units = "px")
MyDimplot_NoHarmony <- DimPlot(phf19.seurat_NoHarmony, reduction = "pca", group.by = "orig.ident")
ggsave("SuppleFig1A_Dimplot_WithoutHarmony.pdf", plot = MyDimplot_NoHarmony, width=5, height=5)

# dev.off()
# biological batch effects
# save plot
jpeg(filename = paste0(getwd(),"/results/pca_samples.jpg"), 
     width = 550,
     height = 400, 
     units = "px")
DimPlot(phf19.seurat_NoHarmony, reduction = "pca", group.by = "hash.ID")
dev.off()
# groups
jpeg(filename = paste0(getwd(),"/results/pca_groups.jpg"), 
     width = 550,
     height = 400, 
     units = "px")
DimPlot(phf19.seurat_NoHarmony, reduction = "pca", group.by = "group")
dev.off()
##################    ##################    ##################    ##################    ##################    ##################  
##################    ##################    ##################    ##################    ##################    ##################  

########################################################################################################################################
### Step4. run PCA    ## If I use SCTransform, I don't need to run ScaleData()           This is quick
########################################################################################################################################
## PCA  ### <<<<<=========   
phf19.seurat_PCA <- Seurat::RunPCA(phf19.seurat_SCT_3000g, npcs=30, verbose=FALSE, features=VariableFeatures(object=phf19.seurat_SCT_3000g))  ##### <<<<===== exclude.gene
saveRDS(phf19.seurat_PCA, "phf19.seurat_PCA_AfterSCT_ExcGene.rds")

####### Dimplot before Harmony Batch adjusting 
MyDimplot_NoHarmony <- DimPlot(phf19.seurat_PCA, reduction = "pca", group.by = "orig.ident")
ggsave("SuppleFig1A_Dimplot_AfterSCT_WithoutHarmony_ExcGene.pdf", plot = MyDimplot_NoHarmony, width=5, height=5)

########################################################################################################################################
### Step5. run Harmony   This is quick
########################################################################################################################################
phf19.seurat_SCTPCAHarmony <- RunHarmony(phf19.seurat_PCA,  group.by.vars = "orig.ident",  plot_convergence=T)
# Error in check_legacy_args(...) : Argument assay.use is unhandled. Please refer to the documentation for the valid harmony options!
packageVersion("harmony") # 2.0.5 # If you're using Harmony ≥1.2 (or a recent version), you should remove assay.use.

saveRDS(phf19.seurat_SCTPCAHarmony, "phf19.seurat_SCTPCAHarmony_ExcGene.rds")

### Dimplot checking
MyDimplot_AfterHarmony <- DimPlot(phf19.seurat_SCTPCAHarmony, reduction = "pca", group.by = "orig.ident")

##################### Dimplot plot After Harmony Integration  #############################
p1 <- Seurat::DimPlot(object=phf19.seurat_SCTPCAHarmony, reduction="harmony", pt.size=.1, group.by="orig.ident") + NoLegend()
p2 <- Seurat::VlnPlot(object=phf19.seurat_SCTPCAHarmony, features="harmony_1", group.by ="orig.ident", pt.size = .1) + NoLegend()
# plot_grid(p1,p2)
ggsave(p1, height=5,width=5, dpi=300, filename=paste0("SuppleFig1A_Dimplot_AfterSCT_AfterHarmony_ExcGene.pdf"), useDingbats=FALSE)

########################################################################################################################################
### Step6.  Checking Before UMAP reduction
########################################################################################################################################
## Get metadata of seurat object and inner_join with Agegroup Annot data. 
PHF19RNAseqMetadata <- phf19.seurat_SCTPCAHarmony@meta.data %>% data.frame %>% dplyr::mutate(CaseID = orig.ident); dim(PHF19RNAseqMetadata) # 33794    16
table(PHF19RNAseqMetadata$CaseID) #
# Lib1  Lib2 
# 18024 15770 
dim(PHF19RNAseqMetadata); PHF19RNAseqMetadata[1:2,]  # 33794    16
#                         orig.ident nCount_RNA nFeature_RNA nCount_HTO nFeature_HTO HTO_maxID HTO_secondID HTO_margin HTO_classification HTO_classification.global hash.ID group percent.mt nCount_SCT nFeature_SCT CaseID
# Lib1_AAACCCAAGGGACCAT-1       Lib1      41635         5492        454            9     Long5     Control2        257              Long5                   Singlet   Long5  Long   1.981506      17797         3981   Lib1
# Lib1_AAACCCAAGTGTACAA-1       Lib1      17348         3347        440            9     Long6     Control2        226              Long6                   Singlet   Long6  Long   2.150104      17193         3327   Lib1

table(PHF19RNAseqMetadata$HTO_maxID)
# Control1 Control2 Control3    Long4    Long5    Long6  Mutant7  Mutant8  Mutant9 
# 4035     3975     3993     2641     3045     4379     3998     3756     3972
table(PHF19RNAseqMetadata$HTO_secondID)

table(PHF19RNAseqMetadata$HTO_maxID, PHF19RNAseqMetadata$HTO_secondID)
#           Control1 Control2 Control3 Long4 Long5 Long6 Mutant7 Mutant8 Mutant9
# Control1        1     2562       59   242   214   143     387     330      97
# Control2      514        2       76   530   288   418     931     894     322
# Control3      157     2430        2   259   216   125     385     325      94
# Long4         118     1651       30     3   137   124     262     254      62
# Long5         161     1856       61   194     1   111     263     305      93
# Long6         180     2783       39   254   220     0     406     379     118
# Mutant7       167     2741       31   262   208   140       0     346     103
# Mutant8       151     2538       29   240   163   121     425       0      89
# Mutant9       156     2473       30   250   184   134     395     350       0


# Control1 row → Control2 = 2562 cells
# Control3 row → Control2 = 2430 cells
# Long6 row    → Control2 = 2783 cells          ### Control2 is frequently the 2nd strongest signal across many cells

#   1) Ambient HTO contamination (most common)
#     One hashtag (here likely Control2) is “bleeding” into many droplets
#     Happens when:   excess antibody, insufficient washing, overloading

#   2) Uneven HTO signal strength
#       Control2 might simply have stronger capture efficiency
#       So it often ranks #2 even when not biologically present

#########   [[[[[[ Sanity chekc ]]]]]]    ########
#### Check classification confidence
table(phf19.seurat_SCTPCAHarmony$HTO_classification.global) # Singlet: 33794 

#### Visualize signal separation
FeatureScatter(phf19.seurat_SCTPCAHarmony, feature1="Control1", feature2="Control2")

#### Visualize signal separation
MyRidgePlot <- RidgePlot(phf19.seurat_SCTPCAHarmony, assay = "HTO", features = rownames(phf19.seurat_SCTPCAHarmony[["HTO"]])[9:1], ncol = 3  )  
ggsave("SuppleFig1B_RidgePlot_AfterSCT_AfterHarmony_ExcGene.pdf", plot = MyRidgePlot, width=10, height=8)

#### Check margin between max and second
hto <- GetAssayData(phf19.seurat_SCTPCAHarmony, assay = "HTO", layer = "data")  
max_val <- apply(hto, 2, max)
second_val <- apply(hto, 2, function(x) sort(x, decreasing=TRUE)[2])

summary(max_val - second_val) ### The result is pretty good.   This is the key separation metric used in HTODemux logic.
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000618 1.1250122 1.4964358 1.5555125 1.9312854 6.3047520 

################################################################################
## Step7.	UMAP Reduction and clustering
################################################################################
phf19.seurat_SCTPCAHarmonyUMAP <- phf19.seurat_SCTPCAHarmony %>% Seurat::RunUMAP(reduction="harmony", dims=1:30, verbose=F) %>% FindNeighbors(reduction="harmony", k.param=15, dim=1:30)  
phf19.seurat_SCTPCAHarmonyUMAP <- phf19.seurat_SCTPCAHarmonyUMAP %>% Seurat::FindClusters(resolution=0.50) %>% identity()
table(phf19.seurat_SCTPCAHarmonyUMAP@active.ident) # Res1.0, 31 clusters; Res1.5, 40 clusters
## resolution 0.2  dims1:30 
#     0     1     2     3     4     5     6     7     8     9    10 
# 15773  5553  5286  3584  1626   447   410   356   341   290   128
## resolution 0.5  dims1:30 
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 7004 5706 5437 3573 2920 1896 1653 1642 1601  585  413  360  353  329  194  128 

saveRDS(phf19.seurat_SCTPCAHarmonyUMAP, file="SeuratUMAPCluster_PHF19_Res0.5PC30KP15.rds") 
# SeuratObject_MouseLiver_snRNAseq <- readRDS(file="SeuratUMAPCluster_MouseLiver_snRNAseq_Res1.5PC30KP15.rds")

### FindMarkers() for all clusters or for each cluster vs all others??? 

## +++++++++. The most reliable “cell ID” is HTO_classification.global (and your filtered hash.ID), not HTO_maxID.
phf19.seurat_SCTPCAHarmonyUMAP$CaseID <- as.character(phf19.seurat_SCTPCAHarmonyUMAP$hash.ID) ; table(phf19.seurat_SCTPCAHarmonyUMAP$CaseID )   
# Control1 Control2 Control3    Long4    Long5    Long6  Mutant7  Mutant8  Mutant9 
# 4035     3975     3993     2641     3045     4379     3998     3756     3972

table(phf19.seurat_SCTPCAHarmonyUMAP$HTO_classification %in% phf19.seurat_SCTPCAHarmonyUMAP$CaseID) # All TRUE: 33794
identical(phf19.seurat_SCTPCAHarmonyUMAP$HTO_classification[1:10] ,  phf19.seurat_SCTPCAHarmonyUMAP$CaseID[1:10] ) # TRUE

### UMAP by cluster ID
UMAP_ByClusterNumber <- Seurat::DimPlot(phf19.seurat_SCTPCAHarmonyUMAP,reduction="umap", label=TRUE, label.size=8) + patchwork::plot_annotation(title="UMAP_Harmony and Clustering number") +
        theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))
UMAP_ByClusterNumber
ggsave(UMAP_ByClusterNumber, height=8,width=11, dpi=300, filename=paste0("Fig1A_UMAPplot_ByClusterNumber_PHF19_Res0.5.pdf"), useDingbats=FALSE)

##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  
##### $$$$$ ###### $$$$$$ ###### Supplementary Figure 1A. UMAP plot of integrated datasets by CaseID - Checking Harmony Integration ##### $$$$$ ###### $$$$$ #####
MyDimplot_SeuratMergeSCTHarmony_ByCaseID <- Seurat::DimPlot(phf19.seurat_SCTPCAHarmonyUMAP,reduction="umap",group.by="CaseID", label=FALSE) + patchwork::plot_annotation(title="UMAP_Harmony and CaseID") +
  theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))
ggsave(MyDimplot_SeuratMergeSCTHarmony_ByCaseID, height=8,width=11, dpi=300, filename=paste0("FigS1C_UMAPplot_PHF19_ByCaseID.pdf"), useDingbats=FALSE)

##### $$$$$ ###### $$$$$$ ###### Supplementary Figure 1A. UMAP plot of integrated datasets by orig.ident - Checking Harmony Integration ##### $$$$$ ###### $$$$$ #####
MyDimplot_SeuratMergeSCTHarmony_ByOrigIdent <- Seurat::DimPlot(phf19.seurat_SCTPCAHarmonyUMAP,reduction="umap",group.by="orig.ident", label=FALSE) + patchwork::plot_annotation(title="UMAP_Harmony and CaseID") +
  theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))
ggsave(MyDimplot_SeuratMergeSCTHarmony_ByOrigIdent, height=8,width=11, dpi=300, filename=paste0("FigS1D_UMAPplot_PHF19_Byorigident.pdf"), useDingbats=FALSE)

##### $$$$$ ###### $$$$$$ ###### Supplementary Figure 1A. UMAP plot of integrated datasets by group - Checking Harmony Integration ##### $$$$$ ###### $$$$$ #####
MyDimplot_SeuratMergeSCTHarmony_Bygroup <- Seurat::DimPlot(phf19.seurat_SCTPCAHarmonyUMAP,reduction="umap",group.by="group",  label=FALSE) + patchwork::plot_annotation(title="UMAP_Harmony and CaseID") +
  theme(axis.text.x=element_text(vjust=0.6, size=15,angle=0), axis.text.y=element_text(vjust=0.6, size=15,angle=0))
ggsave(MyDimplot_SeuratMergeSCTHarmony_Bygroup, height=8,width=10, dpi=300, filename=paste0("FigS1E_UMAPplot_PHF19_Bygroup.pdf"), useDingbats=FALSE)

MyDimplot_SeuratMergeSCTHarmony_BygroupSplit <- Seurat::DimPlot(phf19.seurat_SCTPCAHarmonyUMAP,reduction="umap",group.by="group", split.by="group",ncol=3, label=FALSE) + patchwork::plot_annotation(title="UMAP_Harmony and CaseID") +
  theme(axis.text.x=element_text(vjust=0.6, size=15,angle=0), axis.text.y=element_text(vjust=0.6, size=15,angle=0))
ggsave(MyDimplot_SeuratMergeSCTHarmony_BygroupSplit, height=5,width=12, dpi=300, filename=paste0("FigS1E_UMAPplot_PHF19_BygroupSplit.pdf"), useDingbats=FALSE)
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  


################################################################################
## Step8. Define cell types by marker gene expression   ####
################################################################################
## Check existence of gene in single cell data. 
phf19.seurat_SCTPCAHarmonyUMAP[["SCT"]]@data[1:5,1:3] # Xkr4  Rp1, Mrpl15, Lypla1,  Tcea1  

# grep("Hla", rownames(phf19.seurat_SCTPCAHarmonyUMAP[["SCT"]]@data), value=TRUE)
grep("^H2-", rownames(phf19.seurat_SCTPCAHarmonyUMAP[["SCT"]]@data), value=TRUE) # [1] "H2-K1"   "H2-Ke6"  "H2-Oa"   "H2-DMa"  "H2-DMb2" "H2-DMb1" "H2-Ob" 
grep("Ighg", rownames(phf19.seurat_SCTPCAHarmonyUMAP[["SCT"]]@data), value=TRUE) # [1] "Ighg2c" "Ighg2b" "Ighg1"  "Ighg3" 
grep("Kfp", rownames(phf19.seurat_SCTPCAHarmonyUMAP[["SCT"]]@data), value=TRUE) # None
                                                                     
### Main markers 
MyMainMarker <- c(Plasmacell=c("Prdm1","Derl3","Mzb1","Fcrl5", "Sdc1"),  ## Plasma cell: "Sdc1(Cd138)"  Not exist  "Fkprp11", "Ighg4","Ighgp"
                  Plasmabalst=c("Irf4","Jchain","Ighg1","Tnfrsf17","Cd27", "Mki67"),   ## plasmablast  Tnfrsf17(BCMA) "Cd38","CD27(Human)",    Not exist: , "Hla-dra", "Hla-drb1","Hla-drb5"
                  ActivatedBcell=c("Cd19", "Cd69","Cd83","Cd86","Irf8","Aicda","Cd40","Cd80","Rel","Nfkb1","Nfkbia","Cxcr4","Cxcr5", "Pax5") )
                  ## Activated B-cell)   # Prdm1(Blimp1)  ,"Xbp1","Slpi", # don't exist"B220"

## Markers provided by Hussein 
MyMainMarker <- c(Plasmacell=c("Ighg1"),  ## Plasma cell: "Sdc1(Cd138)"  Not exist  "Fkprp11", "Ighg4","Ighgp"
                  Plasmabalst=c("Cd19","Ptprc","Cxcr3","Xbp1"),   ## plasmablast  Tnfrsf17(BCMA) "Cd38","CD27(Human)",    Not exist: , "Hla-dra", "Hla-drb1","Hla-drb5"
                  ActivatedBcell=c("Bach2", "Ciita","Lgals7","Cd22","Fas") )

MyMainMarkerUnlist <- rev(unlist(MyMainMarker)) ### reverse the order of element
names(MyMainMarkerUnlist)<-NULL ##

MyDotPlot_MainType <- Seurat::DotPlot(phf19.seurat_SCTPCAHarmonyUMAP, features=MyMainMarkerUnlist)+ theme(axis.text.x=element_text(vjust=0.6, size=15,angle=0), axis.text.y=element_text(vjust=0.6, size=15,angle=0)) +
  scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + coord_flip()    # scale_color_viridis_c()    #ffe272
MyDotPlot_MainType

ggsave(MyDotPlot_MainType, height=8,width=8, dpi=300, filename=paste0("OutDotplot_PHF19_scRNAseq_ByHusseinMarker_Res0.5_31Cluster_v1.0.pdf"), useDingbats=FALSE)
## height=16,width=17   by Main markers
################################################################################
## Step9. Assign cell type
################################################################################
CellTypeAssign<- c("Plasmablast", "Plasmablast", "Plasma", "Plasmablast", "Plasmablast", "Plasmablast",
                   "Plasmablast", "ActivatedBcell", "Plasmablast", "Plasmablast", "Plasmablast",
                   "Plasmablast", "Plasmablast", "Plasmablast",  "Plasmablast", "Plasmablast")
names(CellTypeAssign) <- levels(phf19.seurat_SCTPCAHarmonyUMAP)

###################### +++++++++++++++++++ ###################### +++++++++++++++++++ ###################### +++++++++++++++++++ ###################### +++++++++++++++++++ 
###################### +++++++++++++++++++   <<<<<==== I changed Seurat object name here. +++++++++++++++++++++ ###########################
###################### +++++++++++++++++++ ###################### +++++++++++++++++++ ###################### +++++++++++++++++++ 
SeuratObject_PHF19Hashtag_CellType <- Seurat::RenameIdents(phf19.seurat_SCTPCAHarmonyUMAP, CellTypeAssign)
table(SeuratObject_PHF19Hashtag_CellType@active.ident); sum(table(SeuratObject_PHF19Hashtag_CellType@active.ident)) # 33794
# Plasmablast         Plasma ActivatedBcell 
# 26715           5437           1642 

## ==== ##### ====== ##### Add the new cell type definition to Metdata  ## ==== ##### ====== ##### 
SeuratObject_PHF19Hashtag_CellType <- Seurat::AddMetaData(object=SeuratObject_PHF19Hashtag_CellType, metadata=c(SeuratObject_PHF19Hashtag_CellType@active.ident), col.name=c("CellTypeByMarker"))
saveRDS(SeuratObject_PHF19Hashtag_CellType, file="SeuratObject_PHF19Hashtag_CellTypeAssign.rds")    #### <<<<====== use this one for scGSEA input 

##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  
##### $$$$$ ###### $$$$$$ ###### Figure 2B. UMAP plot of integrated datasets by cell type definition ##### $$$$$ ###### $$$$$ #####
# Define your 7 custom colors
my_colors <- c("ActivatedBcell"="#4DAF4A",  "Plasmablast"="#00CCCC",   "Plasma"="#E41A1C")

MyDimplot_CellType <- Seurat::DimPlot(SeuratObject_PHF19Hashtag_CellType, reduction="umap", group.by="CellTypeByMarker", 
                                      cols=my_colors, pt.size=.1, label=TRUE, label.size=8) # + NoLegend(); # label=TRUE,
ggsave(MyDimplot_CellType, height=5.8,width=7.2, dpi=300, filename=paste0("Fig2B_OutUMAP_PHF19_MainCellType_v1.1.pdf"), useDingbats=FALSE)

##### $$$$$ ###### $$$$$$ ###### Figure 2C. UMAP Feature plot by marker gene expression.   ##### $$$$$ ###### $$$$$ #####
# MyFeaturePlot <- Seurat::FeaturePlot(SeuratObject_PHF19Hashtag_CellType,raster=TRUE, cols=c("lightgrey","blue","darkblue"), order=TRUE, pt.size=0.1,
#                                      features=c("Irf8","Pax5", "Tnfrsf17","Mki67", "Prdm1","Derl3" ))  # "Nfkb1": Bcell
MyFeaturePlot <- Seurat::FeaturePlot(SeuratObject_PHF19Hashtag_CellType,raster=TRUE, cols=c("lightgrey","darkgrey","blue"), order=TRUE, pt.size=0.05,
                                     features=c("Irf8","Pax5", "Tnfrsf17","Mki67", "Prdm1","Derl3" ))  # "Nfkb1": Bcell
MyFeaturePlot
ggsave(MyFeaturePlot, height=12.1,width=10, dpi=300, filename=paste0("Fig2C_OutFeatureUMAP_ByMainMarkergeneExpression_v1.1.pdf"), useDingbats=FALSE)

##### $$$$$ ###### $$$$$$ ###### Figure 2CC. UMAP Feature plot by ADDITIONAL marker gene expression by Gareth 20260723  ##### $$$$$ ###### $$$$$ #####

MyFeaturePlotV2 <- Seurat::FeaturePlot(SeuratObject_PHF19Hashtag_CellType,raster=TRUE, cols=c("lightgrey","darkgrey","blue"), order=TRUE, pt.size=0.05,
                                     features=c("Irf8","Pax5","Cd38", ## Activated B cell markers
                                                "Irf4","Tnfrsf17","Ptprc","Xbp1","Jchain","Ighg1",  ## Plasmablast markers
                                                "Prdm1","Mzb1","Scd1"))  # Plasma cell markers
MyFeaturePlotV2
ggsave(MyFeaturePlotV2, height=9,width=14, dpi=300, filename=paste0("Fig2CC_OutFeatureUMAP_ByAdditionalMarker_v1.1.pdf"), useDingbats=FALSE)








############ ================= ############ ================= ############ ================= ############ ================= ############ ================= 
### Violin plot by marker genes. 
SeuratObject_PHF19Hashtag_CellType$CellTypeByMarker <- factor(SeuratObject_PHF19Hashtag_CellType$CellTypeByMarker,
                                                      levels = c("ActivatedBcell", "Plasmablast", "Plasma"))
Idents(SeuratObject_PHF19Hashtag_CellType) <- "CellTypeByMarker"

MyViolinplot_Main <- Seurat::VlnPlot(object=SeuratObject_PHF19Hashtag_CellType, features=c("Irf8","Nfkb1", "Tnfrsf17","Mki67", "Prdm1","Derl3" ), 
                                     stack=TRUE,flip=TRUE ) 
MyViolinplot_Main
ggsave(MyViolinplot_Main, height=6,width=6, dpi=300, filename=paste0("OutViolin_PHF19_MainCellTypeMarker_Res0.2PC30KP15.pdf"), useDingbats=FALSE)

################################################################################
## Step10. Cell type fraction per samples. Make dotplot for macrophage fraction. 
## number of cells per cluster: https://github.com/satijalab/seurat/issues/2825
################################################################################
table(SeuratObject_PHF19Hashtag_CellType$CellTypeByMarker)
# ActivatedBcell    Plasmablast         Plasma 
#       26715           5437           1642 

############ ++++++   A. Count the number of cells in each CaseID 
ContingencyTable_ByCaseID <- table( SeuratObject_PHF19Hashtag_CellType$CaseID, SeuratObject_PHF19Hashtag_CellType$CellTypeByMarker) %>%
                            as.data.frame.matrix; ContingencyTable_ByCaseID[1:2,]
#           Plasmablast Plasma ActivatedBcell
# Control1        3063    586            386
# Control2        2955    653            367

# Calculate row-wise percentages
# margin = 1 calculates the proportion per row
ContingencyTable_ByCaseIDPercent <- as.data.frame(prop.table(as.matrix(ContingencyTable_ByCaseID), margin = 1) * 100)

# View the result
PercentByCaseID_Round2 <- (round(ContingencyTable_ByCaseIDPercent, 2))
fwrite(PercentByCaseID_Round2, file="Table_CellFraction_ByCaseID.txt", col.names=TRUE, row.names=TRUE, sep="\t",quote=FALSE)

############ ++++++   B. Count the number of cells in each cluster per each CaseID 
ContingencyTable_ByCaseIDByClusterNumb <- table( SeuratObject_PHF19Hashtag_CellType$CaseID, SeuratObject_PHF19Hashtag_CellType$seurat_clusters) %>%
  as.data.frame.matrix; ContingencyTable_ByCaseIDByClusterNumb[1:2,]
#           0   1   2   3   4   5   6   7   8  9 10 11 12 13 14 15
# Control1 934 687 586 227 328 203 171 386 251 76 41 42 30 26 28 19

# Calculate row-wise percentages
# margin = 1 calculates the proportion per row
ContingencyTable_ByCaseIDByClusterNumbPercent <- as.data.frame(prop.table(as.matrix(ContingencyTable_ByCaseIDByClusterNumb), margin = 2) * 100)

# View the result
PercentByCaseIDByClusterNumb_Round2 <- (round(ContingencyTable_ByCaseIDByClusterNumbPercent, 2))
fwrite(PercentByCaseIDByClusterNumb_Round2, file="Table_CellFraction_ByCaseIDByClusterNumb.txt", col.names=TRUE, row.names=TRUE, sep="\t",quote=FALSE)

############ ++++++   C. Count the number of cells in each cluster per each Variant type 
ContingencyTable_ByGroupByClusterNumb <- table( SeuratObject_PHF19Hashtag_CellType$group, SeuratObject_PHF19Hashtag_CellType$seurat_clusters) %>%
              as.data.frame.matrix; ContingencyTable_ByGroupByClusterNumb[1:2,]
#           0    1    2    3    4   5   6    7   8   9  10  11  12  13 14 15
# Control 2569 2045 1880  859  936 617 529 1161 667 172 136 130 101  82 75 44

# Calculate row-wise percentages
# margin = 1 calculates the proportion per row
ContingencyTable_ByGroupByClusterNumbPercent <- as.data.frame(prop.table(as.matrix(ContingencyTable_ByGroupByClusterNumb), margin = 2) * 100)

# View the result
ContingencyTable_ByGroupByClusterNumbPercent_Round2 <- (round(ContingencyTable_ByGroupByClusterNumbPercent, 2))
fwrite(ContingencyTable_ByGroupByClusterNumbPercent_Round2, file="Table_CellFraction_ByGroupByClusterNumb.txt", col.names=TRUE, row.names=TRUE, sep="\t",quote=FALSE)

############ ++++++   D.  Count the number of cells in each Group 
ContingencyTable_ByGroup <- table( SeuratObject_PHF19Hashtag_CellType$group, SeuratObject_PHF19Hashtag_CellType$CellTypeByMarker) %>%
              as.data.frame.matrix; ContingencyTable_ByGroup
#         Plasmablast Plasma ActivatedBcell
# Control        8962   1880           1161
# Long           8135   1616            314
# Mutant         9618   1941            167

# Calculate row-wise percentages
# margin = 1 calculates the proportion per row
ContingencyTable_ByGroupPercent <- as.data.frame(prop.table(as.matrix(ContingencyTable_ByGroup), margin = 1) * 100)

# View the result
PercentByGroup_Round2 <- (round(ContingencyTable_ByGroupPercent, 2))
fwrite(PercentByGroup_Round2, file="Table_CellFraction_ByGroup.txt", col.names=TRUE, row.names=TRUE, sep="\t",quote=FALSE)

####### ================ Fraction barplot  =========== #############
SeuratObject_PHF19Hashtag_CellType@meta.data[1,]
#                         orig.ident nCount_RNA nFeature_RNA nCount_HTO nFeature_HTO HTO_maxID HTO_secondID HTO_margin HTO_classification HTO_classification.global hash.ID group percent.mt
# Lib1_AAACCCAAGGGACCAT-1       Lib1      41635         5492        454            9     Long5     Control2        257              Long5                   Singlet   Long5  Long   1.981506
#                         nCount_SCT nFeature_SCT SCT_snn_res.0.5 seurat_clusters CaseID CellTypeByMarker
# Lib1_AAACCCAAGGGACCAT-1      17833         3995               3               3  Long5      Plasmablast

### A.    By CaseID and By CellType
my_colors <- c("ActivatedBcell"="#4DAF4A",  "Plasmablast"="#00CCCC",   "Plasma"="#E41A1C")

CellTypeFractionBarplot <- dittoBarPlot(object = SeuratObject_PHF19Hashtag_CellType, var="CellTypeByMarker",
                                        group.by="CaseID", color.panel = my_colors)   #  group.by="hash.ID"
CellTypeFractionBarplot
## https://github.com/satijalab/seurat/issues/962
ggsave(CellTypeFractionBarplot, height=5, width=4, dpi=300, filename=paste0("Fig2E_OutFractionBarplot_CelltypeBySample.pdf"),  useDingbats=FALSE)


### B.   By CaseID and By CellType in each cluster 
my_colors_Group <- c("Control"="#FF6666",  "Long"="#00CC00", "Mutant"="#0080FF")
# 2.order factor
SeuratObject_PHF19Hashtag_CellType$seurat_clusters <- factor(SeuratObject_PHF19Hashtag_CellType$seurat_clusters,levels=0:15)
table(SeuratObject_PHF19Hashtag_CellType$seurat_clusters)
Idents(SeuratObject_PHF19Hashtag_CellType) <- SeuratObject_PHF19Hashtag_CellType$seurat_clusters

# 3. Set identities to this ordered factor
numeric_order <- c(1, 2, 9, 10, 11, 12, 13, 14, 15, 16, 3, 4, 5, 6, 7, 8)
CellTypeFractionBarplot_EachClstByGroup <- dittoBarPlot(object = SeuratObject_PHF19Hashtag_CellType, var="group",
                                        group.by="seurat_clusters", color.panel = my_colors_Group, x.reorder= numeric_order)   #  group.by="hash.ID"
CellTypeFractionBarplot_EachClstByGroup
ggsave(CellTypeFractionBarplot_EachClstByGroup, height=4, width=6, dpi=300, filename=paste0("Fig2EE_OutFractionBarplot_ByGroupEachCluster.pdf"),  useDingbats=FALSE)

### D.   By GroupID and By CellType
CellTypeFractionBarplot_ByGroup <- dittoBarPlot(object = SeuratObject_PHF19Hashtag_CellType, var="CellTypeByMarker",
                                        group.by="group", color.panel = my_colors)   #  group.by="hash.ID"
CellTypeFractionBarplot_ByGroup
## https://github.com/satijalab/seurat/issues/962
ggsave(CellTypeFractionBarplot_ByGroup, height=5, width=4, dpi=300, filename=paste0("Fig2E_OutFractionBarplot_CelltypeByGroup.pdf"),  useDingbats=FALSE)

################# ======================== ################# ========================################# ========================################# ========================
################# ======================== ################# ========================################# ========================################# ========================

################# ================== ################# ================== ################# ================== 
### +++++++++ If I truly want to remove them from the expression matrix +++++++++++ =============
# genes.keep <- setdiff(rownames(phf19.seurat),  exclude.genes)  #### exclude.genes was made in line 263. 
# SeuratObject_PHF19Hashtag_CellType <- subset( SeuratObject_PHF19Hashtag_CellType, features = genes.keep)
################# ================== ################# ================== ################# ================== 

#################################################################################################
### =========== Step11. Finding diferentially expressed features (cluster biomarkers)
#################################################################################################
table(Idents(SeuratObject_PHF19Hashtag_CellType))
# ActivatedBcell    Plasmablast         Plasma 
# 1642          26715           5437 

SeuratObject_PHF19Hashtag_CellType_PrepSCT <- PrepSCTFindMarkers(SeuratObject_PHF19Hashtag_CellType)  # This takes ~3 min. 

## 1) # find all markers distinguishing  ActivatedBcell vs PlasmaBalst 
Markers_ActivatedBcellvsPlasmablast <- FindMarkers(SeuratObject_PHF19Hashtag_CellType_PrepSCT, ident.1="ActivatedBcell", ident.2 = c("Plasmablast"))  # This takes ~25 min
head(Markers_ActivatedBcellvsPlasmablast, n = 5)
#       p_val avg_log2FC pct.1 pct.2 p_val_adj
# Ly86      0   5.091519 0.967 0.109         0
# Kynu      0   7.284751 0.867 0.010         0
# Pax5      0   4.912969 0.936 0.100         0
saveRDS(Markers_ActivatedBcellvsPlasmablast, "Markers_ActivatedBcellvsPlasmablast.rds")

## 2) # find all markers distinguishing  ActivatedBcell vs PlasmaCell
Markers_ActivatedBcellvsPlasmaCell <- FindMarkers(SeuratObject_PHF19Hashtag_CellType_PrepSCT, ident.1="ActivatedBcell", ident.2 = c("Plasma"))
head(Markers_ActivatedBcellvsPlasmaCell, n = 5)
saveRDS(Markers_ActivatedBcellvsPlasmaCell, "Markers_ActivatedBcellvsPlasmaCell.rds")

## 3) # find all markers distinguishing  ActivatedBcell vs PlasmaCell
Markers_PlasmablastvsPlasmaCell <- FindMarkers(SeuratObject_PHF19Hashtag_CellType_PrepSCT, ident.1="Plasmablast", ident.2 = c("Plasma"))
head(Markers_PlasmablastvsPlasmaCell, n = 5)
saveRDS(Markers_PlasmablastvsPlasmaCell, "Markers_PlasmablastvsPlasmaCell.rds")

#################################################################################################
### =========== Step12. Make CellPhoneDB and GSVA input files  
#################################################################################################
## Mouse Gene Info
MouseGeneInfo <- fread(MouseGeneInfoFile, header=TRUE, stringsAsFactors=FALSE); MouseGeneInfo[1:2,]; dim(MouseGeneInfo) # 112316     16
MouseGeneInfo_Proc <- MouseGeneInfo %>% dplyr::select(c(Symbol,type_of_gene )) %>% dplyr::filter(type_of_gene=="protein-coding"); dim(MouseGeneInfo_Proc); MouseGeneInfo_Proc[1:2,] #  26365     2
#       Symbol   type_of_gene
#       <char>         <char>
#   1:   Pzp2 protein-coding

### 1) Get SCT count data 
CountData_MouseBcell_CellType <- as.data.frame(SeuratObject_PHF19Hashtag_CellType[["SCT"]]@counts); dim(CountData_MouseBcell_CellType)  #  21074 33794 cells

### Remove genes and cells of low expression 
CountData_MouseBcell_CellType_NoLowCountGene <- CountData_MouseBcell_CellType[rowSums(CountData_MouseBcell_CellType)>10, ]; dim(CountData_MouseBcell_CellType_NoLowCountGene); # 10: 14275 33794
## this doesn't filter out any cells. 
# CountData_MouseBcell_CellType_NoLowCountCell <- CountData_MouseBcell_CellType_NoLowCountGene[, colSums(CountData_MouseBcell_CellType_NoLowCountGene)>5000 ]; dim(CountData_MouseBcell_CellType_NoLowCountCell); # 10/1000: 15488 22990

## remove genes of non-proteincoding
CountData_MouseBcell_CellType_ProtCoding <- CountData_MouseBcell_CellType_NoLowCountGene %>% dplyr::filter(rownames(CountData_MouseBcell_CellType_NoLowCountGene) %in% MouseGeneInfo_Proc$Symbol)
dim(CountData_MouseBcell_CellType_ProtCoding) # 11760 33794
saveRDS(CountData_MouseBcell_CellType_ProtCoding, "CountDataSCT_MouseBcell_CellType_ProtCoding.rds")

### 2) Get SCT normalized expression data 
NormalizedExp_MouseBcell_CellType <- as.data.frame(SeuratObject_PHF19Hashtag_CellType[["SCT"]]@data); dim(NormalizedExp_MouseBcell_CellType)  #  21074 33794 cells  # This takes 3 minutes
## View(NormalizedExp_MouseBcell_CellType[1:50,1:20])

## remove genes of non-proteincoding
NormalizedExp_MouseBcell_CellType_ProtCoding <- NormalizedExp_MouseBcell_CellType %>% dplyr::filter(rownames(NormalizedExp_MouseBcell_CellType) %in% MouseGeneInfo_Proc$Symbol)
dim(NormalizedExp_MouseBcell_CellType_ProtCoding) # 15356 33794
saveRDS(NormalizedExp_MouseBcell_CellType_ProtCoding, "ExpNormalizedDataSCT_MouseBcell_CellType_ProtCoding.rds")


### Metadata with cell type
MetadataCelltype <- SeuratObject_PHF19Hashtag_CellType@meta.data %>% data.frame %>% dplyr::mutate(CaseID = orig.ident); dim(MetadataCelltype) # 33794    19
MetadataCelltype[1:2,]
table(MetadataCelltype$HTO_classification)
# Control1 Control2 Control3    Long4    Long5    Long6  Mutant7  Mutant8  Mutant9 
# 4035     3975     3993     2641     3045     4379     3998     3756     3972 
table(MetadataCelltype$seurat_clusters)
#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 7004 5706 5437 3573 2920 1896 1653 1642 1601  585  413  360  353  329  194  128 

MetadataCelltype_Proc <- MetadataCelltype %>% dplyr::select(orig.ident, hash.ID, group, seurat_clusters,CellTypeByMarker )
table(colnames(CountData_MouseBcell_CellType_ProtCoding) %in% rownames(MetadataCelltype_Proc) ) # All TRUE: 33794 
MetadataCelltype_Proc[1:2,]
#               orig.ident hash.ID group seurat_clusters CellTypeByMarker
# Lib1_AAACCCAAGGGACCAT-1       Lib1   Long5  Long               3      Plasmablast
# Lib1_AAACCCAAGTGTACAA-1       Lib1   Long6  Long               0      Plasmablast

saveRDS(MetadataCelltype_Proc, "MetadataCelltype_Proc.rds")

#################################################################################################
### =========== Step14. From SCT normalized expression or Count data, make Pseudobulk data. 
#################################################################################################

# # This automatically calculates rowSums split by both Cell Type and SampleID
# pseudobulk_list <- AggregateExpression(
#   seurat_object,
#   group.by = c("cell_type", "sampleID"),
#   return.seurat = FALSE,
#   slot = "counts" # Essential: DESeq2/EdgeR require unnormalized raw counts
# )
# 
# # Your pseudobulk matrix for RNA will be here:
# raw_counts_matrix <- pseudobulk_list$RNA


################## A From SCT normalized expression
## 1) Long_PlamsaCell 
MetadataCelltype_LongPlasmaCell <- MetadataCelltype_Proc %>% dplyr::filter(grepl("Long", hash.ID)) %>% dplyr::filter(CellTypeByMarker=="Plasma")
dim(MetadataCelltype_LongPlasmaCell) # 1616  3
table(MetadataCelltype_LongPlasmaCell$hash.ID); table(MetadataCelltype_LongPlasmaCell$CellTypeByMarker)
# Control1 Control2 Control3    Long4    Long5    Long6  Mutant7  Mutant8  Mutant9 
# 0        0        0      421      530      665        0        0        0 
# 
# Plasmablast         Plasma ActivatedBcell 
#           0           1616              0 

NormalizedExp_MouseBcell_CellType_LongPlasmaCell <- NormalizedExp_MouseBcell_CellType_ProtCoding %>% 
                  dplyr::select_if(colnames(NormalizedExp_MouseBcell_CellType_ProtCoding) %in% rownames(MetadataCelltype_LongPlasmaCell) )
dim(NormalizedExp_MouseBcell_CellType_LongPlasmaCell) # 15356  1616

## Checking 
table(colnames(NormalizedExp_MouseBcell_CellType_LongPlasmaCell) %in% rownames(MetadataCelltype_LongPlasmaCell)) # 1616 

saveRDS(NormalizedExp_MouseBcell_CellType_LongPlasmaCell, "NormalizedExp_MouseBcell_CellType_LongPlasmaCell.rds")

## 2) Mutant_PlasmaCell
MetadataCelltype_MutantPlasmaCell <- MetadataCelltype_Proc %>% dplyr::filter(grepl("Mutant", hash.ID)) %>% dplyr::filter(CellTypeByMarker=="Plasma")
table(MetadataCelltype_MutantPlasmaCell$hash.ID); table(MetadataCelltype_MutantPlasmaCell$CellTypeByMarker)
# Control1 Control2 Control3    Long4    Long5    Long6  Mutant7  Mutant8  Mutant9 
#     0        0        0        0        0        0      620      679      642 
# 
# Plasmablast         Plasma ActivatedBcell 
#     0           1941              0 

NormalizedExp_MouseBcell_CellType_MutantPlasmaCell <- NormalizedExp_MouseBcell_CellType_ProtCoding %>% 
  dplyr::select_if(colnames(NormalizedExp_MouseBcell_CellType_ProtCoding) %in% rownames(MetadataCelltype_MutantPlasmaCell) )
dim(NormalizedExp_MouseBcell_CellType_MutantPlasmaCell) # 15356g  1941c

## Checking 
table(colnames(NormalizedExp_MouseBcell_CellType_MutantPlasmaCell) %in% rownames(MetadataCelltype_MutantPlasmaCell)) # 1930 

saveRDS(NormalizedExp_MouseBcell_CellType_MutantPlasmaCell, "NormalizedExp_MouseBcell_CellType_MutantPlasmaCell.rds")

## 3) Long_Plamsablast 
MetadataCelltype_LongPlasmablast <- MetadataCelltype_Proc %>% dplyr::filter(grepl("Long", hash.ID)) %>% dplyr::filter(CellTypeByMarker=="Plasmablast")
dim(MetadataCelltype_LongPlasmablast) # 8135  3
table(MetadataCelltype_LongPlasmablast$hash.ID); table(MetadataCelltype_LongPlasmablast$CellTypeByMarker)
# Control1 Control2 Control3    Long4    Long5    Long6  Mutant7  Mutant8  Mutant9 
# 0        0        0     2097     2411     3627        0        0        0 
# 
# ActivatedBcell    Plasmablast         Plasma 
#         0           8135              0 

NormalizedExp_MouseBcell_CellType_LongPlasmablast <- NormalizedExp_MouseBcell_CellType_ProtCoding %>% 
            dplyr::select_if(colnames(NormalizedExp_MouseBcell_CellType_ProtCoding) %in% rownames(MetadataCelltype_LongPlasmablast) )
dim(NormalizedExp_MouseBcell_CellType_LongPlasmablast) # 15356  8135

## Checking 
table(colnames(NormalizedExp_MouseBcell_CellType_LongPlasmablast) %in% rownames(MetadataCelltype_LongPlasmablast)) # 8135

saveRDS(NormalizedExp_MouseBcell_CellType_LongPlasmablast, "NormalizedExp_MouseBcell_CellType_LongPlasmablast.rds")

## 4) Mutant_Plasmablast
MetadataCelltype_MutantPlasmablast <- MetadataCelltype_Proc %>% dplyr::filter(grepl("Mutant", hash.ID)) %>% dplyr::filter(CellTypeByMarker=="Plasmablast")
table(MetadataCelltype_MutantPlasmablast$hash.ID); table(MetadataCelltype_MutantPlasmablast$CellTypeByMarker)
# Control1 Control2 Control3    Long4    Long5    Long6  Mutant7  Mutant8  Mutant9 
#     0        0        0        0        0        0     3329     3015     3274 
# 
# Plasmablast         Plasma ActivatedBcell 
#     9618              0              0 

NormalizedExp_MouseBcell_CellType_MutantPlasmablast <- NormalizedExp_MouseBcell_CellType_ProtCoding %>% 
  dplyr::select_if(colnames(NormalizedExp_MouseBcell_CellType_ProtCoding) %in% rownames(MetadataCelltype_MutantPlasmablast) )
dim(NormalizedExp_MouseBcell_CellType_MutantPlasmablast) # 15356g  9618c

## Checking 
table(colnames(NormalizedExp_MouseBcell_CellType_MutantPlasmablast) %in% rownames(MetadataCelltype_MutantPlasmablast)) # 9618 
saveRDS(NormalizedExp_MouseBcell_CellType_MutantPlasmablast, "NormalizedExp_MouseBcell_CellType_MutantPlasmablast.rds")



################## B From SCT Count data 
Func_RowSumCountPerSmp <- function(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, MetadataCellType=MetadataCelltype_Proc, 
                                   SampleID=c("Long4","Long5","Long6"), SmpGroup="Long", CellType="Plasma", ClusterNumb="All") {
  
    ### Get subset of Metadata by SmpGroup, CellType, and Cluster #
    if(ClusterNumb[1]=="All") {
          MetadataCelltype_SmpGroupCellType <- MetadataCelltype_Proc %>% dplyr::filter(grepl(SmpGroup, group)) %>% dplyr::filter(CellTypeByMarker==CellType)
          dim(MetadataCelltype_SmpGroupCellType) #  1616    5
    } else if ( SmpGroup == "ControlLongMutant") {
          MetadataCelltype_SmpGroupCellType <- MetadataCelltype_Proc %>% dplyr::filter(CellTypeByMarker %in% CellType)
          MetadataCelltype_SmpGroupCellType <- MetadataCelltype_SmpGroupCellType[MetadataCelltype_SmpGroupCellType$seurat_clusters %in% ClusterNumb, ] 
          dim(MetadataCelltype_SmpGroupCellType) #  2524    5
    } else if ( SmpGroup != "ControlLongMutant") {
          MetadataCelltype_SmpGroupCellType <- MetadataCelltype_Proc %>% dplyr::filter(grepl(SmpGroup, group)) %>% dplyr::filter(CellTypeByMarker %in% CellType)
          MetadataCelltype_SmpGroupCellType <- MetadataCelltype_SmpGroupCellType[MetadataCelltype_SmpGroupCellType$seurat_clusters %in% ClusterNumb, ] 
          dim(MetadataCelltype_SmpGroupCellType) #  2524    5
    }
  
    CountData_MouseBcell_CellType_LongPlasmaCell <- CountData_Protcoding %>% 
                    dplyr::select_if(colnames(CountData_Protcoding) %in% rownames(MetadataCelltype_SmpGroupCellType) )
    dim(CountData_MouseBcell_CellType_LongPlasmaCell) # 11772  1635c # Cluster#0: 11772  7004
    
    ## Checking 
    table(colnames(CountData_MouseBcell_CellType_LongPlasmaCell) %in% rownames(MetadataCellType)) # All TRUE 1635 
    # saveRDS(CountData_MouseBcell_CellType_LongPlasmaCell, paste0("CountData_MouseBcell_CellType_",  SmpGroup, CellType, ".rds")) 
    
    CountData_LongAllPlasmaCell <- data.frame(); LoopNumb<-0;
    for(EachSampleID in SampleID) { 
          # EachSampleID <- SampleID[1]; print(paste0("each smp ID:", EachSampleID))
          LoopNumb <- LoopNumb+1; 
          CountData_MouseBcell_CellType_Long4PlasmaCell <- CountData_MouseBcell_CellType_LongPlasmaCell %>% 
                    dplyr::select_if(colnames(CountData_MouseBcell_CellType_LongPlasmaCell) %in% 
                               rownames(MetadataCellType[MetadataCellType$hash.ID==EachSampleID,] ))
          dim(CountData_MouseBcell_CellType_Long4PlasmaCell) # 11760  485c
          
          ### rowSums of gene expression count for all cells in a specific sample
          CountData_Long4PlasmaCell <- data.frame(GeneSymb= rownames(CountData_MouseBcell_CellType_Long4PlasmaCell), 
                                                  ColumnName= rowSums(CountData_MouseBcell_CellType_Long4PlasmaCell)); 
          rownames(CountData_Long4PlasmaCell) <- NULL;
          
          if(ClusterNumb==2) {
                MyCellType<-"Plasmacell" 
          } else if (ClusterNumb==7) {
                MyCellType<-"ActivatedBCell"
          } else {
                MyCellType<-"Plasmablast"
          }
          
          colnames(CountData_Long4PlasmaCell) <- c("GeneSymb",  paste0(EachSampleID, "_", MyCellType))
          dim(CountData_Long4PlasmaCell); CountData_Long4PlasmaCell[1:3,]
          
          if (ClusterNumb[1]!="All") {  ## change column name from "Control1_Plasmablast" to "Control1_Plasmablast_Clst349"
                if(ClusterNumb < 10) {
                      MyClusterNumb <- paste0(0, ClusterNumb)
                } else {
                      MyClusterNumb <- ClusterNumb
                }
            
                colnames(CountData_Long4PlasmaCell)[2] <- paste0(colnames(CountData_Long4PlasmaCell)[2], "_Clst", paste(MyClusterNumb, collapse=""))
          }
          
          ## bind_cols for CountData_LongAllPlasmaCell
          if(LoopNumb==1) { 
            CountData_LongAllPlasmaCell <-  CountData_Long4PlasmaCell
          } else {
            CountData_LongAllPlasmaCell <-  dplyr::inner_join(CountData_LongAllPlasmaCell, CountData_Long4PlasmaCell)
          }
    }
    return(CountData_LongAllPlasmaCell)
}

## =========  1) Long_PlamsaCell 
CountDataPseudo_LongAllPlasmacell <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                            MetadataCellType=MetadataCelltype_Proc, 
                                                            SampleID=c("Long4","Long5","Long6"), SmpGroup="Long", CellType="Plasma") 
dim(CountDataPseudo_LongAllPlasmacell); CountDataPseudo_LongAllPlasmacell[1:2,] # 11772     4
#     GeneSymb Long4_Plasmablast Long5_Plasmablast Long6_Plasmablast
# 1     Xkr4               13                6               16
# 2   Mrpl15              380              417              580
saveRDS(CountDataPseudo_LongAllPlasmacell, "CountDataPseudo_LongAllPlasmacell.rds")


## ========= 2) Mutant_PlasmaCell
CountDataPseudo_MutantAllPlasmacell <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                              MetadataCellType=MetadataCelltype_Proc, 
                                                              SampleID=c("Mutant7","Mutant8","Mutant9"), SmpGroup="Mutant", CellType="Plasma") 
dim(CountDataPseudo_MutantAllPlasmacell); CountDataPseudo_MutantAllPlasmacell[1:2,] # 11772     4
#   GeneSymb Mutant7_PlasmaCell Mutant8_PlasmaCell Mutant9_PlasmaCell
# 1     Xkr4                  6                  9                 12
# 2   Mrpl15                452                549                468
saveRDS(CountDataPseudo_MutantAllPlasmacell, "CountDataPseudo_MutantAllPlasmacell.rds")


## =========  3) Long_Plamsablast 
CountDataPseudo_LongAllPlasmablast <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                             MetadataCellType=MetadataCelltype_Proc, 
                                                             SampleID=c("Long4","Long5","Long6"), SmpGroup="Long", CellType="Plasmablast") 
dim(CountDataPseudo_LongAllPlasmablast); CountDataPseudo_LongAllPlasmablast[1:2,] # 11772     4
#   GeneSymb Long4_Plasmablast Long5_Plasmablast Long6_Plasmablast
# 1     Xkr4                35                41                60
# 2   Mrpl15              2109              2214              3430
saveRDS(CountDataPseudo_LongAllPlasmablast, "CountDataPseudo_LongAllPlasmablast.rds")

## ========= 4) Mutant_Plasmablast
CountDataPseudo_MutantAllPlasmablast <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                               MetadataCellType=MetadataCelltype_Proc, 
                                                               SampleID=c("Mutant7","Mutant8","Mutant9"), SmpGroup="Mutant", CellType="Plasmablast") 
dim(CountDataPseudo_MutantAllPlasmablast); CountDataPseudo_MutantAllPlasmablast[1:2,] # 11772     4
#   GeneSymb Mutant7_Plasmablast Mutant8_Plasmablast Mutant9_Plasmablast
# 1     Xkr4                  61                  48                  42
# 2   Mrpl15                2784                2610                2846
saveRDS(CountDataPseudo_MutantAllPlasmablast, "CountDataPseudo_MutantAllPlasmablast.rds")


############# ================== ############# ================== ############# ================== ############# ==================
##### Pseudo for VariantType_Plasmablast_ClusterNumb
############# ================== ############# ================== ############# ================== ############# ================== 
## ========= 5) Control_Plasmablast for specific clusters #3,#4,#9
CountDataPseudo_ControlPlasmablast_Clst349 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                     MetadataCellType=MetadataCelltype_Proc, 
                                                     SampleID=c("Control1","Control2","Control3"), SmpGroup="Control", CellType="Plasmablast",  ClusterNumb=c(3,4,9)) 
dim(CountDataPseudo_ControlPlasmablast_Clst349); CountDataPseudo_ControlPlasmablast_Clst349[1:2,] # 11772     4
#     GeneSymb Control1_Plasmablast_Clst349 Control2_Plasmablast_Clst349 Control3_Plasmablast_Clst349
# 1     Xkr4                           13                           10                           10
# 2   Mrpl15                          648                          737                          717
saveRDS(CountDataPseudo_ControlPlasmablast_Clst349, "CountDataPseudo_ControlPlasmablast_Clst349.rds")


## ========= 6) Long_Plasmablast for specific clusters #3,#4,#9
CountDataPseudo_LongPlasmablast_Clst349 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                 MetadataCellType=MetadataCelltype_Proc, 
                                                                 SampleID=c("Long4","Long5","Long6"), SmpGroup="Long", CellType="Plasmablast",  ClusterNumb=c(3,4,9)) 
dim(CountDataPseudo_LongPlasmablast_Clst349); CountDataPseudo_LongPlasmablast_Clst349[1:2,] # 11772     4
#   GeneSymb Long4_Plasmablast_Clst349 Long5_Plasmablast_Clst349 Long6_Plasmablast_Clst349
# 1     Xkr4                        10                        18                        19
# 2   Mrpl15                       782                       805                      1345
saveRDS(CountDataPseudo_LongPlasmablast_Clst349, "CountDataPseudo_LongPlasmablast_Clst349.rds")

## ========= 7) Mutant_Plasmablast for specific clusters #3,#4,#9
CountDataPseudo_MutantPlasmablast_Clst349 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                   MetadataCellType=MetadataCelltype_Proc, 
                                                                   SampleID=c("Mutant7","Mutant8","Mutant9"), SmpGroup="Mutant", CellType="Plasmablast",  ClusterNumb=c(3,4,9)) 
dim(CountDataPseudo_MutantPlasmablast_Clst349); CountDataPseudo_MutantPlasmablast_Clst349[1:2,] # 11772     4
#   GeneSymb Mutant4_Plasmablast_Clst349 Mutant5_Plasmablast_Clst349 Mutant6_Plasmablast_Clst349
# 1     Xkr4                          18                          17                          14
# 2   Mrpl15                         951                         925                         981
saveRDS(CountDataPseudo_MutantPlasmablast_Clst349, "CountDataPseudo_MutantPlasmablast_Clst349.rds")


## ========= 8) Control_Plasmablast for specific clusters #0,#1,  #3,#4,#5,#6,   #8,#9,#10,#11,  #13,#14,#15    (Excluding #2,#7,#12)
CountDataPseudo_ControlPlasmablast_NoClst2Clst7Clst12 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                                 MetadataCellType=MetadataCelltype_Proc, 
                                                                                 SampleID=c("Control1","Control2","Control3"), SmpGroup="Control", CellType="Plasmablast", ClusterNumb=c(0,1,  3,4,5,6,  8,9,10,11,  13,14,15)) 
dim(CountDataPseudo_ControlPlasmablast_NoClst2Clst7Clst12); CountDataPseudo_ControlPlasmablast_NoClst2Clst7Clst12[1:2,] # 11772     4
#     GeneSymb Control1_Plasmablast_Clst013456891011131415 Control2_Plasmablast_Clst013456891011131415 Control3_Plasmablast_Clst013456891011131415
# 1     Xkr4                                          48                                          39                                          41
# 2   Mrpl15                                        2414                                        2352                                        2317
saveRDS(CountDataPseudo_ControlPlasmablast_NoClst2Clst7Clst12, "CountDataPseudo_ControlPlasmablast_NoClst2Clst7Clst12.rds")

## ========= 9) Long_Plasmablast for specific clusters  #0,#1,  #3,#4,#5,#6,   #8,#9,#10,#11,  #13,#14,#15    (Excluding #2,#7,#12)
CountDataPseudo_LongPlasmablast_NoClst2Clst7Clst12 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                              MetadataCellType=MetadataCelltype_Proc, 
                                                                              SampleID=c("Long4","Long5","Long6"), SmpGroup="Long", CellType="Plasmablast", ClusterNumb=c(0,1,  3,4,5,6,  8,9,10,11,  13,14,15)) 
dim(CountDataPseudo_LongPlasmablast_NoClst2Clst7Clst12); CountDataPseudo_LongPlasmablast_NoClst2Clst7Clst12[1:2,] # 11772     4
#   GeneSymb Long4_Plasmablast_Clst013456891011131415 Long5_Plasmablast_Clst013456891011131415 Long6_Plasmablast_Clst013456891011131415
# 1     Xkr4                                       35                                       41                                       59
# 2   Mrpl15                                     2102                                     2203                                     3415
saveRDS(CountDataPseudo_LongPlasmablast_NoClst2Clst7Clst12, "CountDataPseudo_LongPlasmablast_NoClst2Clst7Clst12.rds")

## ========= 10) Mutant_Plasmablast for specific clusters  #0,#1,  #3,#4,#5,#6,   #8,#9,#10,#11,  #13,#14,#15    (Excluding #2,#7,#12)
CountDataPseudo_MutantPlasmablast_NoClst2Clst7Clst12 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                                MetadataCellType=MetadataCelltype_Proc, 
                                                                                SampleID=c("Mutant7","Mutant8","Mutant9"), SmpGroup="Mutant", CellType="Plasmablast", ClusterNumb=c(0,1,  3,4,5,6,  8,9,10,11,  13,14,15)) 
dim(CountDataPseudo_MutantPlasmablast_NoClst2Clst7Clst12); CountDataPseudo_MutantPlasmablast_NoClst2Clst7Clst12[1:2,] # 11772     4
#   GeneSymb Mutant7_Plasmablast_Clst013456891011131415 Mutant8_Plasmablast_Clst013456891011131415 Mutant9_Plasmablast_Clst013456891011131415
# 1     Xkr4                                         58                                         48                                         41
# 2   Mrpl15                                       2740                                       2559                                       2768
saveRDS(CountDataPseudo_MutantPlasmablast_NoClst2Clst7Clst12, "CountDataPseudo_MutantPlasmablast_NoClst2Clst7Clst12.rds")



## ========= 11) Control_Plasmablast for Cluster #3
CountDataPseudo_ControlPlasmablast_Clst3 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                   MetadataCellType=MetadataCelltype_Proc, 
                                                                   SampleID=c("Control1","Control2","Control3"), SmpGroup="Control", CellType="Plasmablast",  ClusterNumb=c(3)) 
dim(CountDataPseudo_ControlPlasmablast_Clst3); CountDataPseudo_ControlPlasmablast_Clst3[1:2,] # 11772     4
#     GeneSymb Control1_Plasmablast_Clst3 Control2_Plasmablast_Clst3 Control3_Plasmablast_Clst3
# 1     Xkr4                          9                          9                          4
# 2   Mrpl15                        259                        359                        359
saveRDS(CountDataPseudo_ControlPlasmablast_Clst3, "CountDataPseudo_ControlPlasmablast_Clst3.rds")

## ========= 12) Long_Plasmablast for Cluster #3
CountDataPseudo_LongPlasmablast_Clst3 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                MetadataCellType=MetadataCelltype_Proc, 
                                                                SampleID=c("Long4","Long5","Long6"), SmpGroup="Long", CellType="Plasmablast",  ClusterNumb=c(3)) 
dim(CountDataPseudo_LongPlasmablast_Clst3); CountDataPseudo_LongPlasmablast_Clst3[1:2,] # 11772     4
#   GeneSymb Long4_Plasmablast_Clst3 Long5_Plasmablast_Clst3 Long6_Plasmablast_Clst3
# 1     Xkr4                       8                      17                      15
# 2   Mrpl15                     522                     508                     732
saveRDS(CountDataPseudo_LongPlasmablast_Clst3, "CountDataPseudo_LongPlasmablast_Clst3.rds")

## ========= 13) Mutant_Plasmablast for Cluster #3
CountDataPseudo_MutantPlasmablast_Clst3 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                  MetadataCellType=MetadataCelltype_Proc, 
                                                                  SampleID=c("Mutant7","Mutant8","Mutant9"), SmpGroup="Mutant", CellType="Plasmablast",  ClusterNumb=c(3)) 
dim(CountDataPseudo_MutantPlasmablast_Clst3); CountDataPseudo_MutantPlasmablast_Clst3[1:2,] # 11772     4
#   GeneSymb Mutant4_Plasmablast_Clst3 Mutant5_Plasmablast_Clst3 Mutant6_Plasmablast_Clst3
# 1     Xkr4                         9                        13                        12
# 2   Mrpl15                       455                       582                       600
saveRDS(CountDataPseudo_MutantPlasmablast_Clst3, "CountDataPseudo_MutantPlasmablast_Clst3.rds")


## ========= 14) Control_Plasmablast for Cluster #12
CountDataPseudo_ControlPlasmablast_Clst12 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                    MetadataCellType=MetadataCelltype_Proc, 
                                                                    SampleID=c("Control1","Control2","Control3"), SmpGroup="Control", CellType="Plasmablast",  ClusterNumb=c(12)) 
dim(CountDataPseudo_ControlPlasmablast_Clst12); CountDataPseudo_ControlPlasmablast_Clst12[1:2,] # 11772     4
#     GeneSymb Control1_Plasmablast_Clst12 Control2_Plasmablast_Clst12 Control3_Plasmablast_Clst12
# 1     Xkr4                           0                           1                           1
# 2   Mrpl15                          16                          15                          23
saveRDS(CountDataPseudo_ControlPlasmablast_Clst12, "CountDataPseudo_ControlPlasmablast_Clst12.rds")

## ========= 15) Long_Plasmablast for Cluster #12
CountDataPseudo_LongPlasmablast_Clst12 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                 MetadataCellType=MetadataCelltype_Proc, 
                                                                 SampleID=c("Long4","Long5","Long6"), SmpGroup="Long", CellType="Plasmablast",  ClusterNumb=c(12)) 
dim(CountDataPseudo_LongPlasmablast_Clst12); CountDataPseudo_LongPlasmablast_Clst12[1:2,] # 11772     4
#   GeneSymb Long4_Plasmablast_Clst12 Long5_Plasmablast_Clst12 Long6_Plasmablast_Clst12
# 1     Xkr4                        0                        0                        1
# 2   Mrpl15                        7                       11                       15
saveRDS(CountDataPseudo_LongPlasmablast_Clst12, "CountDataPseudo_LongPlasmablast_Clst12.rds")

## ========= 16) Mutant_Plasmablast for Cluster #12
CountDataPseudo_MutantPlasmablast_Clst12 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountData_MouseBcell_CellType_ProtCoding, 
                                                                   MetadataCellType=MetadataCelltype_Proc, 
                                                                   SampleID=c("Mutant7","Mutant8","Mutant9"), SmpGroup="Mutant", CellType="Plasmablast",  ClusterNumb=c(12)) 
dim(CountDataPseudo_MutantPlasmablast_Clst12); CountDataPseudo_MutantPlasmablast_Clst12[1:2,] # 11772     4
#   GeneSymb Mutant4_Plasmablast_Clst12 Mutant5_Plasmablast_Clst12 Mutant6_Plasmablast_Clst12
# 1     Xkr4                          3                          0                          1
# 2   Mrpl15                         44                         51                         78
saveRDS(CountDataPseudo_MutantPlasmablast_Clst12, "CountDataPseudo_MutantPlasmablast_Clst12.rds")


############# ================== ############# ================== ############# ================== ############# ==================
##### Pseudo for Each cluster with Control, Long, and Mutant. 
############# ================== ############# ================== ############# ================== ############# ================== 

AllClusterNumber <- c(0,1, 2, 3,4,5,6, 7, 8,9,10,11,12,13,14,15)
Func_PseudobulkPerCluster <- function(CountNormalizedData, AllClusterNumber, MetadataCelltype_Proc, DataType) { 
      for(EachClusterNumber in AllClusterNumber) {
              # EachClusterNumber <- AllClusterNumber[1]
              ## ========= 20) Control, Long, and Mutant Cluster #0
              CountDataPseudo_ControlLongMutantPlasmablast_Clst0 <- Func_RowSumCountPerSmp(CountData_Protcoding= CountNormalizedData, 
                                                              MetadataCellType=MetadataCelltype_Proc, 
                                                              SampleID=c("Control1","Control2","Control3", "Long4","Long5","Long6",
                                                                         "Mutant7","Mutant8","Mutant9"), SmpGroup="ControlLongMutant", CellType=c("ActivatedBcell", "Plasmablast", "Plasma"), 
                                                                                                                              ClusterNumb=EachClusterNumber)
              
              dim(CountDataPseudo_ControlLongMutantPlasmablast_Clst0); CountDataPseudo_ControlLongMutantPlasmablast_Clst0[1:2,] # 11772     4
              #   GeneSymb Mutant4_Plasmablast_Clst12 Mutant5_Plasmablast_Clst12 Mutant6_Plasmablast_Clst12
              # 1     Xkr4                          3                          0                          1
              # 2   Mrpl15                         44                         51                         78
              
              if(EachClusterNumber < 10) {
                    MyClusterNumb <- paste0(0, EachClusterNumber)
              } else {
                    MyClusterNumb <- EachClusterNumber
              }
              saveRDS(CountDataPseudo_ControlLongMutantPlasmablast_Clst0, paste0(DataType, "DataPseudo_ControlLongMutantPlasmablast_Clst", MyClusterNumb,".rds") )
      }
}

### 1) CountData 
Func_PseudobulkPerCluster(CountNormalizedData=CountData_MouseBcell_CellType_ProtCoding, AllClusterNumber=AllClusterNumber, MetadataCelltype_Proc=MetadataCelltype_Proc, DataType="Count") 

### 2) Normalized Exp Datta
Func_PseudobulkPerCluster(CountNormalizedData=NormalizedExp_MouseBcell_CellType_ProtCoding, AllClusterNumber=AllClusterNumber, MetadataCelltype_Proc=MetadataCelltype_Proc, DataType="NormalizedExp") 





  