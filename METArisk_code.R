# R package
library(Seurat)
library(stringr)
library(ggplot2)
library(grid)
library(gtools)
library(clustree)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(cols4all)
library(viridis)
library(tidyverse)
library(magrittr)
library(reshape2)
library(ggsci)
library(data.table)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(clusterProfiler)
library(GSEABase)
library(survival)
library(survminer)
library(limma)
library(My.stepwise)
library(tinyarray)
library(oncoPredict)
library(readxl)
library(NMF)
library(estimate)

# Data import----
# Deposited data was converted to Seurat objects (one example shown for brevity):
# Metastatic:ME01,ME02...
# Primary:PC01,PC02...

sample1.data <- Read10X(data.dir = "/path/to/data/folder/")
sample1 <- CreateSeuratObject(counts = sample1.data, project = "sample1")

sceList <- list(sample1,...,sample33)

samples <- c("sample1",...,"sample33")

sce <- merge(x = sceList[[1]], y = sceList[-1], 
                    add.cell.ids = samples,
                    merge.data = TRUE)
sce$Group = ifelse(str_detect(sce$orig.ident,"^ME"),"Metastatic","Primary")

# Data quality control----
# Datasets underwent quality control as follows:
mito_genes <- rownames(sce)[grep("^MT-", rownames(sce))]
sce <- PercentageFeatureSet(sce, "^MT-", col.name = "pMT")
nFeature_lower <- 200
nFeature_upper <- 8000
selected_c <- WhichCells(sce, expression = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper) 
selected_f <- rownames(sce)[Matrix::rowSums(sce@assays$RNA@counts > 0 ) > 3]
sce <- subset(sce, features = selected_f, cells = selected_c) 
pMT_upper <- 15
selected_mito <- WhichCells(sce, expression = pMT < pMT_upper)
sce <- subset(sce, cells = selected_mito)
sce <- subset(sce, cells = selected_mito)

# Data Preparation----
scelist <- SplitObject(sce, split.by = "orig.ident")
scelist
  
for (i in 1:length(scelist)) {
    print(i)
    scelist[[i]] <- NormalizeData(scelist[[i]], verbose = FALSE)
    scelist[[i]] <- FindVariableFeatures(scelist[[i]], verbose = FALSE,selection.method = "vst", nfeatures = 2000)
}

features <- SelectIntegrationFeatures(object.list = scelist)
  
scelist <- lapply(X = scelist, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
alldata.anchors <- FindIntegrationAnchors(object.list = scelist, dims = 1:30, 
                                            reduction = "rpca")

sce <- IntegrateData(anchorset = alldata.anchors, dims = 1:30)
  

sce <- ScaleData(sce,verbose = FALSE)
sce <- RunPCA(sce, npcs = 30,verbose = FALSE)
sce <- RunUMAP(sce, dims = 1:30)
sce <- FindNeighbors(sce,dims = 1:30) 
for (i in seq(0.1,1,0.1)) {
  sce <- FindClusters(sce, resolution = i)
  print(DimPlot(sce, reduction = "umap") + 
          labs(title = paste0("resolution: ", i)))
}

clustree(sce)

# Cell Classification----
DefaultAssay(sce) <- "RNA"
all.markers <- FindAllMarkers(sce, 
                              only.pos = TRUE, 
                              min.pct = 0.1, 
                              logfc.threshold = 0.25)
significant.markers  <- all.markers[all.markers$p_val_adj < 0.05, ]

Idents(sce) <- sce$RNA_snn_res.1

annotation_curated_main <- read_excel("./annotation/annotation_main.xlsx")
new_ids_main <- annotation_curated_main$Main_cell_type
names(new_ids_main) <- levels(sce)
sce <- RenameIdents(sce, new_ids_main)
sce <- subset(sce,idents = "Other",invert = TRUE)
levels(sce) <- c("Malignant",
                 "Fibroblast",
                 "Endothelial",
                 "T",
                 "B_Plasma",
                 "Monocytic",
                 "Mast")
sce@meta.data$Main_cell_type <- Idents(sce)

# Malignant Subset & METArisk----
Mali=sce[,sce$Main_cell_type %in% c('Malignant')] 

Malilist <- SplitObject(Mali, split.by = "orig.ident")
Malilist

for (i in 1:length(Malilist)) {
  print(i)
  Malilist[[i]] <- NormalizeData(Malilist[[i]], verbose = FALSE)
  Malilist[[i]] <- FindVariableFeatures(Malilist[[i]], verbose = FALSE,selection.method = "vst", nfeatures = 2000)
}

features <- SelectIntegrationFeatures(object.list = Malilist)

Malilist <- lapply(X = Malilist, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

alldata.anchors <- FindIntegrationAnchors(object.list = Malilist, dims = 1:30, 
                                          reduction = "rpca")

Mali <- IntegrateData(anchorset = alldata.anchors, dims = 1:30)


Mali <- ScaleData(Mali,verbose = FALSE)
Mali <- RunPCA(Mali, npcs = 30,verbose = FALSE)
Mali <- RunUMAP(Mali, dims = 1:30)
Mali <- FindNeighbors(Mali,dims = 1:30) 
for (i in seq(0.1,1,0.1)) {
  Mali <- FindClusters(Mali, resolution = i)
  print(DimPlot(Mali, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}

clustree(Mali)

# Cell Classification
DefaultAssay(Mali) <- "RNA"
all.markers <- FindAllMarkers(Mali, 
                              only.pos = TRUE, 
                              min.pct = 0.1, 
                              logfc.threshold = 0.25)
significant.markers  <- all.markers[all.markers$p_val_adj < 0.05, ]

Idents(Mali) <- Mali$RNA_snn_res.0.8

annotation_curated_main <- read_excel("./annotation/annotation_Mali.xlsx")
new_ids_main <- annotation_curated_main$Cluster
names(new_ids_main) <- levels(Mali)
Mali <- RenameIdents(Mali, new_ids_main)

levels(Mali) <- c(paste0(rep("PrimaryC",5),1:5),
                  paste0(rep("MixedC",2),1:2),
                  paste0(rep("MetastaticC",2),1:2))
Mali@meta.data$Cluster <- Idents(Mali)

mycolblue <- cols4all::c4a("brewer.blues",10)[5:9]
mycolred <- cols4all::c4a("brewer.reds",10)[7:8]
mycolpurple <- c(cols4all::c4a("brewer.bu_pu",10)[7],
                cols4all::c4a("brewer.purples",10)[7])

mycol <- c(mycolblue,mycolpurple,mycolred)

# Pseudobulk
Cell <- rep("Malignant",nrow(Mali@meta.data))
Mali$Cell <- Cell
DE <- run_de(Mali,replicate_col = "orig.ident", cell_type_col = "Cell", label_col = "Group",
            de_family = "pseudobulk",
            de_method = "edgeR")
logFC_cutoff <- 0.25
DE$change <- as.factor(ifelse(DE$p_val < 0.05 & abs(DE$avg_logFC) > logFC_cutoff,
                                   ifelse(DE$avg_logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
DE$symbol <- DE$gene
s2e <- bitr(DE$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
DE <- inner_join(DE,s2e,by=c("symbol"="SYMBOL"))

# Functional Analysis
gene_up <- DE[DE$change %in% "UP","symbol"]
gene_down <- DE[DE$change %in% "DOWN","symbol"]

gene_up <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,keys = gene_up,columns = 'ENTREZID',keytype = 'SYMBOL')[,2]))
gene_down <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,keys = gene_down,columns = 'ENTREZID',keytype = 'SYMBOL')[,2]))

run_kegg <- function(gene_up,gene_down,geneList=F,pro='test'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  kk=kk.up
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.up.csv'))
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  kk=kk.down
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.down.csv'))
}

run_kegg(gene_up,gene_down,pro=paste0("./Functional_Analysis/KEGG"))

Keggall <- read.gmt("./kegghsa.gmt") # download from GSEA database
Keggmeta <- filter(Keggall,Keggall$term %in% c("hsa01100_Metabolic_pathways")|str_detect(Keggall$term,"metabolism$"))

head(KEGG)
head(Keggmeta)

gene <- unique(as.vector(DE$symbol))
gene <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_df <- data.frame(logFC= DE$avg_log2FC,
                      SYMBOL = DE$gene)
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList <- gene_df$logFC
names(geneList) <- gene_df$SYMBOL
geneList <- sort(geneList, decreasing = TRUE)
gsea <- GSEA(geneList,TERM2GENE = META_ACTIVATE)

# AUCell
geneSets <- lapply(unique(META_ACTIVATE$term), function(x){Hallmarker$gene[META_ACTIVATE$term == x]})
names(geneSets) <- unique(META_ACTIVATE$term)
names(geneSets) <- "META_ACTIVATE"

exprMatrix <- as.matrix(Mali@assays$RNA@data)
cells_rankings <- AUCell_buildRankings(exprMatrix,splitByBlocks=TRUE) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=T, assign=TRUE)
warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)

warningMsg[which(warningMsg =="")]
sapply(cells_assignment, function(x) x$aucThr$selected)

selectedThresholds <- getThresholdSelected(cells_assignment)
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"

head(assignmentTable)

assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
selectedThresholds <- getThresholdSelected(cells_assignment)
dat <- data.frame(Mali@meta.data, Mali@reductions$umap@cell.embeddings)
umap <- data.frame(Mali@reductions$umap@cell.embeddings)
umap <- as.matrix(umap)
plot(umap, cex=.3)
assignmentTable[1:4,]

for(geneSetName in names(selectedThresholds)){
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(5)
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(5)
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0 )  {
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=5)], names(aucSplit[[1]])),
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=5)], names(aucSplit[[2]])))
    plot(umap, main=geneSetName,
         sub="Pink/red cells pass the threshold",
         col=cellColor[rownames(umap)], pch=16)
  }
}
par(mfrow=c(1,1))
AUCell_plotTSNE(tSNE=umap, exprMat=exprMatrix, 
                cellsAUC=cells_AUC[c(1),], thresholds=selectedThresholds[c(1)])
geneSet <- "META_ACTIVATE"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
Mali$AUC  <- aucs
par(mfrow=c(1,1))
umap <- data.frame(Mali@meta.data, Mali@reductions$umap@cell.embeddings)

# Cibersortx
# Cibersortx follows the analysis flow of https://cibersortx.stanford.edu/
# Here we showed the downstream analyses
# TCGA-HNSC was set as the bulk-RNA data
# exp: Expression of TCGA-HNSC all samples
# exprSet: Expression of TCGA-HNSC tumor samples
# meta: etadata for patients

ciberx <- read.csv("CIBERSORTx_Results.csv",row.names = 1)
ciberx <- merge(ciberx,meta)
ciberx$METArisk <- ciberx$META_Activate/(ciberx$META_Activate + ciberx$META_Silence + ciberx$Non_Malignant)
res.cut <- surv_cutpoint(ciberx, 
                         time = "time", 
                         event = "event", 
                         variables = c("METArisk"))
summary(res.cut) 
plot(res.cut, "METArisk", palette = "npg")

ciberx$Group <- ifelse(ciberx$METArisk > res.cut$cutpoint$cutpoint,"METArisk_high","METArisk_low")
  
# PYGL selection & relevant mechanism----
# Tumor vs Normal
table(str_sub(colnames(exp),14,15))
Group <- ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
GroupTN <- factor(Group,levels = c("normal","tumor"))
design <- model.matrix(~GroupTN)
fit <- lmFit(exprSet,design)
fit <- eBayes(fit)
deg <- topTable(fit,coef=2,number = Inf)
logFC_t <- 1
adj.P.value_t <- 0.05
k1 <- (deg$adj.P.value < adj.P.value_t)&(deg$logFC < -logFC_t)
k2 <- (deg$adj.P.value < adj.P.value_t)&(deg$logFC > logFC_t)
deg <- mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
TNdeg <- deg 

# METArisk_high vs METArisk_low
meta <- merge(ciberx,meta)
Group <- meta$Group
GroupMETA <- factor(Group,levels = c("METArisk_low","METArisk_high"))
design <- model.matrix(~GroupMETA)
fit <- lmFit(exprSet,design)
fit <- eBayes(fit)
deg <- topTable(fit,coef=2,number = Inf)
logFC_t <- 1
adj.P.value_t <- 0.05
k1 <- (deg$adj.P.value < adj.P.value_t)&(deg$logFC < -logFC_t)
k2 <- (deg$adj.P.value < adj.P.value_t)&(deg$logFC > logFC_t)
deg <- mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
METAdeg <- deg


# SVM-RFE
source("./SVM-RFE-master/msvmRFE.R")
candicate <- intersect(Diff(TNdeg),Diff(METAdeg),MRG)
set.seed(12345) 
exprSetsvm <- exprSet[candicate,]
input <- cbind(meta[,"event"],t(exprSetsvm))
colnames(input)[1] <- "event"
svmRFE(input, k=10, halve.above=20)

nfold <- 10
nrows <- nrow(input)
folds <- rep(1:nfold, len=nrows)[sample(nrows)]
folds <- lapply(1:nfold, function(x) which(folds == x))

results <- lapply(folds, svmRFE.wrap, input, k=10, halve.above=20)
top.features <- WriteFeatures(results, input, save=F)

featsweep <- lapply(1:10, FeatSweep.wrap, results, input)
no.info <- min(prop.table(table(input[,1])))
errors <- sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
PlotErrors(errors, no.info=no.info)

candicate <- top.features$FeatureName[1:16]
  
# cox-forest
e <- t(exprSet[candicate,])
dat <- cbind(meta,e)
colnames(dat)
vl <- colnames(dat)[c(3:ncol(dat))]
My.stepwise.coxph(Time = "time",
                  Status = "event",
                  variable.list = vl,
                  data = dat)
model <- coxph(formula = Surv(time, event) ~  EPHX3+FDCSP+IDO1+PYGL+FAM3B, 
               data = dat)

ggforest(model,fontsize = 1,data = dat)

# OncoPredict
dir='./DataFiles/Training Data/'
GDSC2_Expr <- readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res <- readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

testExpr <- exprSet
colnames(testExpr) <- paste0('test',colnames(testExpr))
dim(testExpr)  

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )


# GSVA
geneSets <- getGmt('metabolism_associated_signatures.gmt')# Compiled from the pathcard website
GSVA_hall <- gsva(expr=as.matrix(exprSet), 
                  gset.idx.list=geneSets, 
                  mx.diff=T, 
                  kcdf="Gaussian", 
                  parallel.sz=4) 

p <- identical(meta$Samples,colnames(GSVA_hall));p
if(!p) GSVA_hall <- GSVA_hall[,match(meta$Samples,colnames(GSVA_hall))]

design <- model.matrix(~0+meta$PYGL)
colnames(design) = levels(factor(meta$PYGL))
rownames(design) = colnames(GSVA_hall)

compare <- makeContrasts(PYGL_high - PYGL_low, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)

GSHgene <- c("GCLC","GCLM","GSS","OPLAH","GGT1","GGT5","GGT7",  
             "GSTA1","GSTM1","GSTM2","GGCT","GSTO2","GGT6","GSTA3","GSTA4", 
             "GSTM3","GSTM4","GSTM5","GSTP1","GSTT2","GSTZ1","GSTK1","MGST1", 
             "MGST2","MGST3","CHAC2","CNDP2","GSTT2B","CHAC1","GSTO1","GSR",  
             "ANPEP","GPX1","GPX2","GPX3","GPX4","AKR1A1","ESD","HPGDS", 
             "IDH1","G6PD","GPX7")
target_gene <- 'PYGL'

expst <- exprSet[c(target_gene,GSHgene),]
t_exp <- t(expst)
t_exp <- as.data.frame(t_exp)

cor_list <- list()
for (i in colnames(t_exp)) {
  tar <- t_exp[,target_gene]
  cor_res <- cor(x = tar,y = t_exp[,i],method = 'pearson')
  cor_pval <- cor.test(x = tar,y = t_exp[,i])$p.value
  final_res <- data.frame(tar_genename = target_gene,
                          gene_name = i,
                          cor_results = cor_res,
                          cor_pvalue = cor_pval)
  cor_list[[i]] <- final_res
}

gene_corres <- do.call('rbind',cor_list)
head(gene_corres,4)

high_cor <- gene_corres %>% filter(abs(cor_results) >= 0.5 ,cor_pvalue < 0.05 )
high_corgene <- high_cor$gene_name
length(high_corgene)

high_corgene <- high_corgene[-match("PYGL",high_corgene)]

# NMF
data <- exprSet[GSHgene,]
res <- nmf(data, rank=2:10, method="brunet", nrun=10, seed=123456)
plot(res)
consensusmap(res,
             annRow=NA,
             annCol=NA,
             main="Consensus matrix",
             info=FALSE)
clusterNum <- 2       
res <- nmf(data, rank=clusterNum, method="brunet", nrun=10, seed=123456)
Cluster <- predict(res)
Cluster <- as.data.frame(Cluster)
Cluster$Cluster <- paste0("Cluster", Cluster$Cluster)
clusterOut <- rbind(ID=colnames(Cluster), Cluster)
consensusmap(res,
             annRow=NA,
             annCol=NA,
             main="Consensus matrix", 
             info=FALSE)

Cluster$Samples <- rownames(Cluster)
Cluster <- merge(Cluster,meta)

# ESTIMATE
write.table(exprSet,file = "expression.txt",sep = "\t",quote = F)

exp.file <- "expression.txt"
in.gct.file <- "ESTIMATE_input.gct"

outputGCT(exp.file, in.gct.file)

filterCommonGenes(input.f = "expression.txt",
                  output.f = "ESTIMATE_input.gct",
                  id = "GeneSymbol")


out.score.file = "ESTIMATE_score.gct"
estimateScore(in.gct.file, 
              out.score.file, 
              platform = "illumina")

ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[,2:ncol(ESTIMATE_score)]))
ESTIMATE_score$Samples = rownames(ESTIMATE_score)
ESTIMATE_score = ESTIMATE_score[,c(ncol(ESTIMATE_score),2:ncol(ESTIMATE_score)-1)]

ESTIMATE_score$Samples = str_replace_all(ESTIMATE_score$Samples,fixed("."),"-")

ESTIMATE_score = merge(ESTIMATE_score,Cluster)

