
library(limma)
library(ConsensusClusterPlus)

workDir=""     
setwd(workDir)

rt=read.table("H_ssgseaOut_gse.txt", header=T, sep="\t", check.names=F, row.names=1)
data=rt[,(3:(ncol(rt)))]
data=t(data)

maxK=9     
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="pam",         
              distance="euclidean",     
              seed=123456,      
              plot="png")     


clusterNum=2      ֵ
Cluster=results[[clusterNum]][["consensusClass"]]
outTab=cbind(rt, Cluster)
outTab[,"Cluster"]=paste0("C", outTab[,"Cluster"])
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="cluster.txt", sep="\t", quote=F, col.names=F)

