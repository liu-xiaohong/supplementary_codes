
library(pheatmap)
setwd("")      
rt=read.table(" ",sep="\t",header=T,row.names=1,check.names=F)   

Type=read.table("cluster.txt",sep="\t",check.names=F,header=T)
Type[,2]=paste0("Cluster",Type[,2])
Type=Type[order(Type[,2]),]
rt=rt[,as.vector(Type[,1])]

cluster=as.data.frame(Type[,2])
row.names(cluster)=Type[,1]
colnames(cluster)="Cluster"

pdf("heatmap.pdf",height=5,width=9)
pheatmap(rt, annotation=cluster, 
         color = colorRampPalette(c("#BA55D3", "white", "#FF8C00"))(100),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3)
dev.off()


ClusterC1="ERS_Lipid_H"
ClusterC2="ERS_Lipid_L"
ClusterC3="ERS_Lipid_H"

#????????????
a=c()
a[Type[,2]=="ClusterC1"]=ClusterC1
a[Type[,2]=="ClusterC2"]=ClusterC2
a[Type[,2]=="ClusterC3"]=ClusterC3
clusterOut=cbind(Type,a)
write.table(clusterOut,file="cluster.txt",sep="\t",quote=F,col.names=F,row.names=F)
