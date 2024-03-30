
library(limma)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(CellChat)

expFile=   
annFile=
geneFile=
setwd()   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

meta=read.table(annFile, header=T, sep="\t", check.names=F, row.names=1)

cellchat <- createCellChat(object = data, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use="labels")
groupSize <- as.numeric(table(cellchat@idents))      


CellChatDB <- CellChatDB.human      
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

pdf(file="COMM01.DatabaseCategory.pdf", width=7, height=5)
showDatabaseCategory(CellChatDB)
dev.off()
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)      
cellchat <- identifyOverExpressedInteractions(cellchat)      
cellchat <- projectData(cellchat, PPI.human)  


cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net=subsetCommunication(cellchat)
write.table(file="COMM02.Comm.network.xls", df.net, sep="\t", row.names=F, quote=F)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
file="COMM03.cellNetworkCount.pdf", width=7, height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
#???(file="COMM04.cellNetworkWeight.pdf", width=7, height=6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength")
dev.off()

#??�="COMM05.singleCell.pdf", width=8, height=6)
weight_mat <- cellchat@net$weight
par(mfrow = c(2,3), mgp=c(0,0,0), xpd=TRUE)
for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=0.8,title.name=cel)
}
dev.off()

#???????=paste0("COMM06.bubble.pdf"), width=8, height=5.5)
netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x = 45)
dev.off()

#??ȡ???ead.csv(geneFile, header=T, sep=",", check.names=F)
hubGenes=as.vector(geneRT[,1])
def.hub=df.net[((df.net$ligand %in% hubGenes) | (df.net$receptor %in% hubGenes)),]
write.table(file="COMM07.Comm.hubNetwork.xls", def.hub, sep="\t", row.names=F, quote=F)

#ͨ·ˮ�@netP$pathways     #չʾ??s.show="SPP1"       #ѡ???�e=paste0("COMM08.", pathways.show , ".circle.pdf"), width=8, height=6)
circle=netVisual_aggregate(cellchat, signaling=pathways.show, layout="circle", vertex.size = groupSize)
print(circle)
dev.off()
#ʹ?ò??e=paste0("COMM09.", pathways.show , ".hierarchy.pdf"), width=12, height=6)
hierarchy=netVisual_aggregate(cellchat, signaling=pathways.show, layout="hierarchy",  vertex.receiver=seq(1,4), vertex.size = groupSize)
print(hierarchy)
dev.off()
#ϸ??ͨ�e=paste0("COMM10.", pathways.show , ".heatmap.pdf"), width=8, height=6)
heatmap=netVisual_heatmap(cellchat, signaling=pathways.show, color.heatmap = "Reds", measure= 'weight')	
print(heatmap)
dev.off()
#ϸ?????e=paste0("COMM11.", pathways.show , ".netAnalysis.pdf"), width=6, height=5)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis=netAnalysis_signalingRole_network(cellchat, signaling =pathways.show, width = 8, height = 5, font.size = 12)
print(netAnalysis)
dev.off()

#?鿴??�=paste0("COMM12.", pathways.show , ".contribution.pdf"), width=8, height=6)
contribution=netAnalysis_contribution(cellchat, signaling= pathways.show)
print(contribution)
dev.off()
#?鿴ͨ·?л????????ı???ˮƽ
pdf(file=paste0("COMM13.", pathways.show , ".geneExp.pdf"), width=8, height=6)
geneExp=plotGeneExpression(cellchat, signaling=pathways.show)
print(geneExp)
dev.off()

#??????<- extractEnrichedLR(cellchat, signaling=pathways.show, geneLR.return=FALSE)
pdf(file=paste0("COMM14.", pathways.show , ".pairLR.pdf"), width=9, height=8)
pairCircos=netVisual_individual(cellchat, signaling=pathways.show, pairLR.use=pairLR[1] , layout="circle" )
print(pairCircos)
dev.off()
#??ͨ·?n 1:nrow(pairLR)){
	pdf(file=paste0("COMM15.", pairLR[i,], ".pairLR.pdf"), width=8, height=6)
	pairChord=netVisual_individual(cellchat, signaling=pathways.show, pairLR.use=pairLR[i,], layout="chord" )
	print(pairChord)
	dev.off()
}


######?