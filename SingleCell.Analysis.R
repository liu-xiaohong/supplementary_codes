
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)

logFCfilter=1              
adjPvalFilter=0.05          
inputFile=" "       
setwd("")       


rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

pbmc=CreateSeuratObject(counts = data,project = "seurat", min.cells=3, min.features=50, names.delim = "_")
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")(file="01.featureViolin.pdf", width=10, height=6)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    #??file="01.featureCor.pdf",width=10,height=6)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#???pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf",width=10,height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()





pbmc=ScaleData(pbmc)         CA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA????
?P02.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#???????ɷ02.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#???ɷַ??02.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#?õ?ÿ??P
ckStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:15)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:15)
dev.off()





######4
pbmc <- indNeighbors(object = pbmc, dims = 1:pcSelect)       #?????ڽӾindClusters(object = pbmc, resolution = 0.5)         #??ϸ????unTSNE(object = pbmc, dims = 1:pcSelect)             #TSNE??"03.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    #TSNE???ӻ
write.table(pbmc$seurat_clusters,file="03.tsneCluster.txt",quote=F,sep="\t",col.names=F)

##????ÿ??ers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#????marke03.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()

#????marke03.markerViolin.pdf",width=10,height=6)
VlnPlot(object = pbmc, features = row.names(sig.markers)[1:2])
dev.off()

#??Ҫչʾ?03.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = pbmc, features = showGenes, cols = c("green", "red"))
dev.off()

#????marke03.markerBubble.pdf",width=12,height=6)
cluster10Marker=showGenes
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()




######
mc@assays$RNA@counts
clusters<-pbmc@meta.data$seurat_clusters
ann=pbmc@meta.data$orig.ident
#ref=get(load("ref_Human_all.RData"))
ref=celldex::HumanPrimaryCellAtlasData()
singler=SingleR(test=counts, ref =ref,
                labels=ref$label.main, clusters = clusters)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
singler2=SingleR(test=counts, ref =ref, 
                labels=ref$label.main)
cellAnn=as.data.frame(singler2)
cellAnn=cbind(id=row.names(cellAnn), cellAnn)
cellAnn=cellAnn[,c("id", "labels")]
write.table(cellAnn, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#clusterע?ͺ??Ŀ??ӻ?
newLabels=singler$labels
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, newLabels)
pdf(file="04.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    #TSNE???

##cluster�rs=FindAllMarkers(object = pbmc,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.cellMarkers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.cellMarkers,file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)



#########trix=as.matrix(pbmc@assays$RNA@data)
monocle.sample=pbmc@meta.data
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.clusterAnn=clusterAnn
monocle.markers=sig.markers

#??Seurat?(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#????ϸ???=as.character(monocle.clusterAnn[,2])
names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

#ϸ???켣?imateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, as.vector(sig.markers$gene))
#plot_ordeduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)
#??????֦?"05.trajectory.State.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "State")
dev.off()
#????ʱ???"05.trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()
#????ϸ???"05.trajectory.cellType.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
#?????????"05.trajectory.cluster.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
dev.off()

#ϸ???켣
set(pData(cds),select='State')
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
geneList=list()
for(i in levels(factor(groups$State))){
	pbmc.markers=FindMarkers(pbmc, ident.1 = i, group.by = 'group')
	sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
	sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
	write.table(sig.markers,file=paste0("05.monocleDiff.", i, ".txt"),sep="\t",row.names=F,quote=F)
	geneList[[i]]=row.names(sig.markers)
}
#???潻???s=Reduce(union,geneList)
write.table(file="05.monocleDiff.union.txt",unionGenes,sep="\t",quote=F,col.names=F,row.names=F)


######Vid