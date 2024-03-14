# Section 1: Load Required Packages

library('recount3')    # To process gene expression
library('limma')      # To average replicates
library('DESeq2')     # To do the mean-variance correction
library('vsn')        # To visualize mean-variance plots
library('biomaRt')    # To convert ensembl gene/protein ids to hugo (hgnc) symbols
#library('netZooR')    # To load PANDA
library('visNetwork') # For network visualization 
library('fgsea')      # For gene set enrichment analysis 
library('ggplot2')    # To build the bubble plot for gene enrichment
library('dplyr')
library('edgeR')
library('ggrepel')
library('gridExtra')
#library('tidyverse')


# Section 2: Import Pancreatic Tissue Data from TCGA and GTEx and also lists of X and Y chromosome genes

# Import TCGA Human Pancreatic Adenocarcinoma Data From Recount 3

rse_geneTumorDataOriginal <- recount3::create_rse_manual( project = "PAAD",
                                                          project_home = "data_sources/tcga",
                                                          organism = "human",
                                                          annotation = "gencode_v26",
                                                          type = "gene"
)

# Import GTEx Pancreas Data from Recount 3

rse_geneNormalDataOriginal <- recount3::create_rse_manual(
  project = "PANCREAS",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

# Import Human TF List

human_tf <- as.data.frame(read.table("http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt", header = FALSE, sep = "\n", dec = "."))
human_tf=human_tf[["V1"]]

# Import Y chromosome genes

Ygene <- as.vector(read.table("https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&status=Approved&chr=Y&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&where=(gd_pub_chrom_map%20not%20like%20%27%25patch%25%27%20and%20gd_pub_chrom_map%20not%20like%20%27%25alternate%20reference%20locus%25%27)&submit=submit", header = TRUE, sep = "\n"))

# Import X chromosomes genes

Xgene <- as.vector(read.table("https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&status=Approved&chr=X&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&where=(gd_pub_chrom_map%20not%20like%20%27%25patch%25%27%20and%20gd_pub_chrom_map%20not%20like%20%27%25alternate%20reference%20locus%25%27)&submit=submit", header = TRUE, sep = "\n"))

# Section 3: Create VSD creator function (also filter RSE Objects)

vsdcreator <- function(rseobject, genetype, includePlots, includecountMat, isTumor){
  rse_gene <- rseobject
  rowDataRse=rowData(rse_gene)
  colDataRse=colData(rse_gene)
  
  geneIds = c()
  geneNames = c()
  geneEntrezIds = c()
  
  for(i in 1:dim(rowDataRse)[1]){
    if(rowDataRse@listData$gene_type[i][[1]][1] %in% genetype & !(rowDataRse@listData$gene_name[i][[1]][1] %in% Xgene$Approved.symbol) & !(rowDataRse@listData$gene_name[i][[1]][1] %in% Ygene$Approved.symbol)){
      geneIds       = c(geneIds,i)
      geneNames     = c(geneNames,rowDataRse@listData$gene_name[i][[1]][1])
      geneEntrezIds = c(geneEntrezIds, substr(rowDataRse@listData$gene_id[i],1,15))
    }
  }
  
  if(isTumor == TRUE){
    tumorIds=c()
    patientSymbol=c()
    
    for(i in 1:length(colDataRse@listData$tcga.gdc_cases.samples.sample_type)){
      if(colDataRse@listData$tcga.gdc_cases.samples.sample_type[i] %in% c("Primary Tumor") && colDataRse@listData$tcga.gdc_cases.project.name[i] == "Pancreatic Adenocarcinoma"){
        tumorIds = c(tumorIds,i)
        patientSymbol = c(patientSymbol, colDataRse@listData$tcga.gdc_cases.samples.portions.analytes.aliquots.submitter_id[i])
      }
    }
    
    rse_gene=rse_gene[geneIds,tumorIds]
  }
  
  if(isTumor == FALSE){
    rse_gene=rse_gene[geneIds,]
  }
  
  
  countMat=SummarizedExperiment::assay(rse_gene, 1)
  
  indFilter = which(rowSums(countMat) <= 1)
  
  rse_gene=rse_gene[-indFilter,]
  
  rowDataRse=rowData(rse_gene)
  colDataRse=colData(rse_gene)
  
  countMat=SummarizedExperiment::assay(rse_gene, 1)
  meanSdPlotBeforeCorrection <- meanSdPlot(countMat)
  meanSdPlotBeforeCorrectionLogTrans <- meanSdPlot(log(countMat+1), ranks = FALSE)
  
  if (identical(genetype,c("miRNA"))){
    vsd <- varianceStabilizingTransformation(countMat, blind=FALSE)
  } else{
    vsd <- vst(countMat, blind=FALSE)
  }
  
  meanSdPlot <- meanSdPlot(vsd)
  
  
  if(isTumor == TRUE){
    colnames(vsd)=patientSymbol
  }
  rownames(vsd)=geneNames[-indFilter]
  
  vsd=avearrays(t(vsd))
  vsd=as.data.frame(t(vsd))
  
  if(isTumor == TRUE){
    colnames(countMat)=patientSymbol
  }
  rownames(countMat)=geneNames[-indFilter]
  countMat=avearrays(t(countMat))
  countMat=as.data.frame(t(countMat))
  
  if (includePlots == FALSE & includecountMat == FALSE){
    return(list(rse = rse_gene, vsdoutput = vsd))
  }
  if (includePlots == TRUE & includecountMat == FALSE){
    return(list(rse = rse_gene,vsdoutput = vsd, meanSdPlotBeforeCorrectionoutput = meanSdPlotBeforeCorrection, meanSdPlotBeforeCorrectionLogTransoutput = meanSdPlotBeforeCorrectionLogTrans, meanSdPlotoutput = meanSdPlot))
  }
  
  if (includePlots == TRUE & includecountMat == TRUE){
    return(list(rse = rse_gene, vsdoutput = vsd, meanSdPlotBeforeCorrectionoutput = meanSdPlotBeforeCorrection, meanSdPlotBeforeCorrectionLogTransoutput = meanSdPlotBeforeCorrectionLogTrans, meanSdPlotoutput = meanSdPlot, countMatoutput = countMat))
  }
  
  if (includePlots == FALSE & includecountMat == TRUE){
    return(list(rse = rse_gene, vsdoutput = vsd, countMatoutput = countMat))
  }
}

# Section 4: Pancreatic Cancer Data

vsdcreatorTumorOutput <- vsdcreator(rseobject = rse_geneTumorDataOriginal, genetype = c("protein_coding"), includePlots = TRUE, includecountMat = TRUE, isTumor = TRUE)

rse_geneTumorData <- vsdcreatorTumorOutput[[1]]
vsdTumorData <- vsdcreatorTumorOutput[[2]]
totalTumormeanSdPlotBeforeCorrectionoutput <- vsdcreatorTumorOutput[[3]]
totalTumormeanSdPlotBeforeCorrectionLogTrans <- vsdcreatorTumorOutput[[4]]
totalTumormeanSdPlot <- vsdcreatorTumorOutput[[5]]
totalTumorcountMat <- vsdcreatorTumorOutput[[6]]

vsdcreatorTumorOutputlincRNA <- vsdcreator(rseobject = rse_geneTumorDataOriginal, genetype = c("lincRNA"), includePlots = TRUE, includecountMat = TRUE, isTumor = TRUE)

rse_geneTumorDatalincRNA <- vsdcreatorTumorOutputlincRNA[[1]]
vsdTumorDatalincRNA <- vsdcreatorTumorOutputlincRNA[[2]]
totalTumormeanSdPlotBeforeCorrectionoutputlincRNA <- vsdcreatorTumorOutputlincRNA[[3]]
totalTumormeanSdPlotBeforeCorrectionLogTranslincRNA <- vsdcreatorTumorOutputlincRNA[[4]]
totalTumormeanSdPlotlincRNA <- vsdcreatorTumorOutputlincRNA[[5]]
totalTumorcountMatlincRNA <- vsdcreatorTumorOutputlincRNA[[6]]

vsdcreatorTumorOutputmiRNA <- vsdcreator(rseobject = rse_geneTumorDataOriginal, genetype = c("miRNA"), includePlots = TRUE, includecountMat = TRUE, isTumor = TRUE)

rse_geneTumorDatamiRNA <- vsdcreatorTumorOutputmiRNA[[1]]
vsdTumorDatamiRNA <- vsdcreatorTumorOutputmiRNA[[2]]
totalTumormeanSdPlotBeforeCorrectionoutputmiRNA <- vsdcreatorTumorOutputmiRNA[[3]]
totalTumormeanSdPlotBeforeCorrectionLogTransmiRNA <- vsdcreatorTumorOutputmiRNA[[4]]
totalTumormeanSdPlotmiRNA <- vsdcreatorTumorOutputmiRNA[[5]]
totalTumorcountMatmiRNA <- vsdcreatorTumorOutputmiRNA[[6]]

# Section 5: Normal Pancreas Data

vsdcreatorNormalOutput <- vsdcreator(rseobject = rse_geneNormalDataOriginal, genetype = c("protein_coding"), includePlots = TRUE, includecountMat = TRUE, isTumor = FALSE)

rse_geneNormalData <- vsdcreatorNormalOutput[[1]]
vsdNormalData <- vsdcreatorNormalOutput[[2]]
totalNormalmeanSdPlotBeforeCorrectionoutput <- vsdcreatorNormalOutput[[3]]
totalNormalmeanSdPlotBeforeCorrectionLogTrans <- vsdcreatorNormalOutput[[4]]
totalNormalmeanSdPlot <- vsdcreatorNormalOutput[[5]]
totalNormalcountMat <- vsdcreatorNormalOutput[[6]]

vsdcreatorNormalOutputlincRNA <- vsdcreator(rseobject = rse_geneNormalDataOriginal, genetype = c("lincRNA"), includePlots = TRUE, includecountMat = TRUE, isTumor = FALSE)

rse_geneNormalDatalincRNA <- vsdcreatorNormalOutputlincRNA[[1]]
vsdNormalDatalincRNA <- vsdcreatorNormalOutputlincRNA[[2]]
totalNormalmeanSdPlotBeforeCorrectionoutputlincRNA <- vsdcreatorNormalOutputlincRNA[[3]]
totalNormalmeanSdPlotBeforeCorrectionLogTranslincRNA <- vsdcreatorNormalOutputlincRNA[[4]]
totalNormalmeanSdPlotlincRNA <- vsdcreatorNormalOutputlincRNA[[5]]
totalNormalcountMatlincRNA <- vsdcreatorNormalOutputlincRNA[[6]]

vsdcreatorNormalOutputmiRNA <- vsdcreator(rseobject = rse_geneNormalDataOriginal, genetype = c("miRNA"), includePlots = TRUE, includecountMat = TRUE, isTumor = FALSE)

rse_geneNormalDatamiRNA <- vsdcreatorNormalOutputmiRNA[[1]]
vsdNormalDatamiRNA <- vsdcreatorNormalOutputmiRNA[[2]]
totalNormalmeanSdPlotBeforeCorrectionoutputmiRNA <- vsdcreatorNormalOutputmiRNA[[3]]
totalNormalmeanSdPlotBeforeCorrectionLogTransmiRNA <- vsdcreatorNormalOutputmiRNA[[4]]
totalNormalmeanSdPlotmiRNA <- vsdcreatorNormalOutputmiRNA[[5]]
totalNormalcountMatmiRNA <- vsdcreatorNormalOutputmiRNA[[6]]


vsdcreatorTumorplm <- vsdcreator(rseobject = rse_geneTumorDataOriginal, genetype = c("protein_coding", "lincRNA", "miRNA"), includePlots = TRUE, includecountMat = TRUE, isTumor = TRUE)
rse_geneTumorplm <- vsdcreatorTumorplm[[1]]
vsdTumorplm <- vsdcreatorTumorplm[[2]]

vsdcreatorNormalplm <- vsdcreator(rseobject = rse_geneNormalDataOriginal, genetype = c("protein_coding", "lincRNA", "miRNA"), includePlots = TRUE, includecountMat = TRUE, isTumor = FALSE)
rse_geneNormalplm <- vsdcreatorNormalplm[[1]]
vsdNormalplm <- vsdcreatorNormalplm[[2]]

# Section 6: Tumor Sex Subset Function

tumorsexsubsetter <- function(rseobject, sex){
  if (sex == "M"){
    tcgaidentifier = "male"
    colnamesidentifier = "; gender=m"
  }
  if (sex == "F"){
    tcgaidentifier = "female"
    colnamesidentifier = "; gender=f"
  }
  
  rse_gene <- rseobject
  colDataRse <- colData(rseobject)
  tumorIds = c()
  patientSymbol = c()
  
  for(i in 1:length(colDataRse@listData$tcga.gdc_cases.demographic.gender)){
    if((colDataRse@listData$tcga.gdc_cases.demographic.gender[i] %in% c(tcgaidentifier)) && (colDataRse@listData$tcga.gdc_cases.project.name[i] == "Pancreatic Adenocarcinoma")){
      tumorIds = c(tumorIds,i)
      patientSymbol = c(patientSymbol, colDataRse@listData$tcga.gdc_cases.samples.portions.analytes.aliquots.submitter_id[i])
    }
  }
  
  rse_gene=rse_gene[,tumorIds]
  colnames(rse_gene) = patientSymbol
  colnames(rse_gene) = paste0(colnames(rse_gene), colnamesidentifier)
  
  return(list(rse = rse_gene, Ids = tumorIds))
  
}

tumorsexsubsetter2 <- function(rseobject, sex){
  if (sex == "M"){
    tcgaidentifier = "male"
  }
  if (sex == "F"){
    tcgaidentifier = "female"
  }
  
  rse_gene <- rseobject
  colDataRse <- colData(rseobject)
  tumorIds = c()
  patientSymbol = c()
  
  for(i in 1:length(colDataRse@listData$tcga.gdc_cases.demographic.gender)){
    if((colDataRse@listData$tcga.gdc_cases.demographic.gender[i] %in% c(tcgaidentifier)) && (colDataRse@listData$tcga.gdc_cases.project.name[i] == "Pancreatic Adenocarcinoma")){
      tumorIds = c(tumorIds,i)
      patientSymbol = c(patientSymbol, colDataRse@listData$tcga.gdc_cases.samples.portions.analytes.aliquots.submitter_id[i])
    }
  }
  
  rse_gene=rse_gene[,tumorIds]
  colnames(rse_gene) = patientSymbol
  
  return(list(rse = rse_gene, Ids = tumorIds))
  
}

# VSD Sex subsetter

sexsubsettervsd <- function(rse_gene, vsd){
  list <- rse_gene@colData@rownames
  for (i in 1:length(list)){
    list[i] <- substr(list[i], 1, nchar(list[i])-10)
  }
  vsdoutput <- vsd[,list]
  
  return(vsdoutput)
  
}

sexsubsettervsd2 <- function(rse_gene, vsd){
  list <- rse_gene@colData@rownames
  vsdoutput <- vsd[,list]
  
  return(vsdoutput)
  
}

# Section 7: Normal Sex Subset Function

normalsexsubsetter <- function(rseobject, sex){
  if (sex == "M"){
    gtexidentifier = 1
    colnamesidentifier = "; gender=m"
  }
  if (sex == "F"){
    gtexidentifier = 2
    colnamesidentifier = "; gender=f"
  }
  
  rse_gene <- rseobject
  colDataRse <- colData(rseobject)
  normalIds = c()
  patientSymbol = c()
  
  for(i in 1:length(colDataRse@listData$gtex.sex)){
    if((colDataRse@listData$gtex.sex[i] == gtexidentifier)){
      normalIds = c(normalIds,i)     
      patientSymbol = c(patientSymbol, colDataRse@listData$subjid[i])
    }
  }
  rse_gene=rse_gene[,normalIds]
  colnames(rse_gene) = paste0(colnames(rse_gene), colnamesidentifier)
  
  return(list(rse = rse_gene, Ids = normalIds))
  
}

normalsexsubsetter2 <- function(rseobject, sex){
  if (sex == "M"){
    gtexidentifier = 1
  }
  if (sex == "F"){
    gtexidentifier = 2
  }
  
  rse_gene <- rseobject
  colDataRse <- colData(rseobject)
  normalIds = c()
  patientSymbol = c()
  
  for(i in 1:length(colDataRse@listData$gtex.sex)){
    if((colDataRse@listData$gtex.sex[i] == gtexidentifier)){
      normalIds = c(normalIds,i)     
      patientSymbol = c(patientSymbol, colDataRse@listData$subjid[i])
    }
  }
  rse_gene=rse_gene[,normalIds]
  
  return(list(rse = rse_gene, Ids = normalIds))
  
}

# Section 8: topTable and VolcanoPlot generator

topTableAndVolcanoPlotGenerator <- function(efit, logFCbound, adjpvaluebound){
  toptable <- topTable(efit, sort.by = "logFC", n = Inf, coef = 1)
  toptable <- toptable[order(toptable$P.Value),]
  
  toptable$diffexpressed <- "NO"
  toptable$diffexpressed[toptable$logFC > logFCbound & toptable$P.Value < adjpvaluebound] <- "UP"
  toptable$diffexpressed[toptable$logFC < -logFCbound & toptable$P.Value < adjpvaluebound] <- "DOWN"
  toptable$delabel <- NA
  toptable$delabel[toptable$diffexpressed != "NO"] <- rownames(toptable)[toptable$diffexpressed != "NO"]
  volcanoplot <- ggplot(data=toptable, aes(x=logFC, y= -log10(P.Value), col=diffexpressed, label=delabel)) + geom_point() + theme_minimal() + geom_text() + scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-logFCbound, logFCbound), col="red") +
    geom_hline(yintercept=-log10(adjpvaluebound), col="red")
  
  toptable1 <- toptable
  toptable <- filter(toptable, abs(logFC) > logFCbound & P.Value < adjpvaluebound)
  geneListComparison <- as.vector(rownames(toptable))
  geneListComparisonUP <- as.vector(rownames(toptable)[toptable$diffexpressed == "UP"])
  geneListComparisonDOWN <- as.vector(rownames(toptable)[toptable$diffexpressed == "DOWN"])
  
  
  return(list(tt=toptable, vp=volcanoplot, geneList=geneListComparison, geneListUP=geneListComparisonUP, geneListDOWN=geneListComparisonDOWN, tt1=toptable1))
}

# Section 9

rse_geneMaleTumor <- tumorsexsubsetter(rse_geneTumorData, "M")[[1]]
vsdMaleTumor <- sexsubsettervsd(rse_geneMaleTumor, vsdTumorData)

rse_geneMaleNormal <- normalsexsubsetter(rse_geneNormalData, "M")[[1]]
vsdMaleNormal <- sexsubsettervsd(rse_geneMaleNormal, vsdNormalData)

vsdMaleTumorCopy <- vsdMaleTumor
colnames(vsdMaleTumorCopy) <- paste0("tumor", colnames(vsdMaleTumorCopy))
vsdMaleNormalCopy <- vsdMaleNormal
colnames(vsdMaleNormalCopy) <- paste0("norma", colnames(vsdMaleNormalCopy))

vsdDataMale <- merge(x=cbind(gene=rownames(vsdMaleTumorCopy),vsdMaleTumorCopy),y=cbind(gene=rownames(vsdMaleNormalCopy),vsdMaleNormalCopy),by="gene",all.x=TRUE, all.y=TRUE)
vsdDataMale <- na.omit(vsdDataMale)
rownames(vsdDataMale) <- vsdDataMale[,1]
vsdDataMale <- vsdDataMale[-c(1)]

namesMale = colnames(vsdDataMale)
groupMale = substr(namesMale, 1, 5)
vsdDataMale$samples$group = groupMale #Assign samples to appropriate group
designMale = model.matrix(~0 + groupMale)
colnames(designMale) = gsub("group", "", colnames(designMale))
contr.matrixMale = makeContrasts(TumorMaleandNormalMale = Maletumor-Malenorma, levels = colnames(designMale))

# Fit a linear model using weighted least square for each gene

vfitDataMale = lmFit(vsdDataMale[,1:ncol(vsdDataMale)-1], designMale)

vfitDataMaleContr = contrasts.fit(vfitDataMale, contrasts = contr.matrixMale)

#  Perform empirical Bayes smoothing of standard errors
head(coef(vfitDataMale))
efitDataMale = eBayes(vfitDataMaleContr)

ttvpMale <- topTableAndVolcanoPlotGenerator(efitDataMale, 5, 0.05)

topTableMale <- ttvpMale[[1]]
volcanPlotMale <- ttvpMale[[2]]
geneComparisonListMaleAll <- ttvpMale[[3]]
geneComparisonListMaleUP <- ttvpMale[[4]]
geneComparisonListMaleDOWN <- ttvpMale[[5]]
topTableMale1 <- ttvpMale[[6]]

# Gene Regulatory Network

Normaloutput <- vsdcreator(rseobject = rse_geneNormalDataOriginal, genetype = c("protein_coding","lincRNA","miRNA"), includePlots = FALSE, includecountMat = FALSE, isTumor = FALSE)
Tumoroutput <- vsdcreator(rseobject = rse_geneTumorDataOriginal, genetype = c("protein_coding","lincRNA","miRNA"), includePlots = FALSE, includecountMat = FALSE, isTumor = TRUE)
NormalMalerse <- normalsexsubsetter2(Normaloutput[[1]], "M")[[1]]
NormalMalevsd <- sexsubsettervsd2(NormalMalerse, Normaloutput[[2]])
TumorMalerse <- tumorsexsubsetter2(Tumoroutput[[1]], "M")[[1]]
TumorMalevsd <- sexsubsettervsd2(TumorMalerse, Tumoroutput[[2]])

NormalMaleMatrix <- as.data.frame(t(NormalMalevsd))
TumorMaleMatrix <- as.data.frame(t(TumorMalevsd))

ListNMTMComp <- intersect(colnames(NormalMaleMatrix), colnames(TumorMaleMatrix))

CorNormalMaleMatrix <- as.data.frame(cor(NormalMaleMatrix[,ListNMTMComp]))
CorTumorMaleMatrix <- as.data.frame(cor(TumorMaleMatrix[,ListNMTMComp]))

CorNormalMaleMatrix <- CorNormalMaleMatrix %>% select(as.vector(intersect(ListNMTMComp, rownames(vsdMaleNormal))))
CorNormalMaleMatrix <- as.data.frame(subset(CorNormalMaleMatrix, !(rownames(CorNormalMaleMatrix) %in% colnames(CorNormalMaleMatrix))))

CorTumorMaleMatrix <- CorTumorMaleMatrix %>% select(as.vector(intersect(ListNMTMComp, rownames(vsdMaleTumor))))
CorTumorMaleMatrix <- as.data.frame(subset(CorTumorMaleMatrix, !(rownames(CorTumorMaleMatrix) %in% colnames(CorTumorMaleMatrix))))

MaleComparisonMatrix <- CorTumorMaleMatrix - CorNormalMaleMatrix
MaleComparisonMatrixrownames = rownames(MaleComparisonMatrix)

num_edges <- 100

edges = matrix(0L, num_edges, 3)
colnames(edges) = c("from","to","value")
edges = as.data.frame(edges)

sort_mat = order(as.matrix(abs(MaleComparisonMatrix)), decreasing = TRUE)
edges$value  = as.matrix(MaleComparisonMatrix)[sort_mat[1:num_edges]]

pcgeneIdsTop = (sort_mat[1:num_edges] %/% dim(MaleComparisonMatrix)[1]) + 1
rnaIdsTop = sort_mat[1:num_edges] %% dim(MaleComparisonMatrix)[1]
nrnas = dim(MaleComparisonMatrix)[1]
rnaIdsTop[rnaIdsTop == 0] = nrnas

edges$to = colnames(MaleComparisonMatrix)[pcgeneIdsTop]
edges$from = rownames(MaleComparisonMatrix)[rnaIdsTop]
edges$arrows = "to"
edges$color = ifelse(edges$value > 0, "green", "red")
edges$value = abs(edges$value)

nodes = data.frame(id = unique(as.vector(as.matrix(edges[,c(1,2)]))), label=unique(as.vector(as.matrix(edges[,c(1,2)]))))
#nodes$group = ifelse(nodes$id %in% edges$from, "miRNA/lincRNA", "Protein-coding Gene")
idmirna <- c()
idlincrna <- c()
idtfs <- c()
mirnavector <- as.vector(union(rownames(vsdNormalDatamiRNA),rownames(vsdTumorDatamiRNA)))
lincrnavector <- as.vector(union(rownames(vsdNormalDatalincRNA),rownames(vsdTumorDatalincRNA)))

for(i in 1:dim(nodes)[1]){
  if(nodes[i,1] %in% human_tf){
    idtfs = c(idtfs, i)
  }
  if(nodes[i,1] %in% mirnavector){
    idmirna = c(idmirna, i)
  }
  if(nodes[i,1] %in% lincrnavector){
    idlincrna = c(idlincrna, i)
  }
}
nodes$group = 'Protein-coding Gene'
nodes$group[idmirna] = 'miRNA'
nodes$group[idlincrna] = 'lincrna'
nodes$group[idtfs] = 'tfs'

net <- visNetwork(nodes, edges, width = "100%")
net <- visGroups(net, groupname = "tfs", shape = "triangle",
                 color = list(background = "green", border="green"))
net <- visGroups(net, groupname = "lincrna", shape = "triangle",
                 color = list(background = "pink", border="pink"))
net <- visGroups(net, groupname = "miRNA", shape = "triangle",
                 color = list(background = "orange", border="orange"))
net <- visGroups(net, groupname = "Protein-coding Gene", shape = "dot",       
                 color = list(background = "darkblue", border="black"))
visLegend(net, main="Legend", position="right", ncol=1)