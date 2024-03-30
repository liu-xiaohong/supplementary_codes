
library(limma)
library(reshape2)
library(ggplot2)

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

eRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), rownames(data))
data=t(data[sameGene,])

#??È=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-1),drop=F]

#???ab=data.frame()
for(checkpiont in colnames(data)){
	for(gene in colnames(risk)){
		x=as.numeric(data[,checkpiont])
		y=as.numeric(risk[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, checkpiont=checkpiont, cor, text, pvalue))
	}
}

#???Tab$Gene=factor(outTab$Gene, levels=colnames(risk))
outTab$cor=as.numeric(outTab$cor)
pdf(file="checkpointCor.pdf", width=8, height=7)
ggplot(outTab, aes(Gene, checkpiont)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +    #È¥?ô±³¾?
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),   #x??????
	 .text.y = element_text(size = 10, face = "bold")) +       #y??????
	l=paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #????Í¼??
	screte(position = "bottom")      #X????????Ê
######Video 