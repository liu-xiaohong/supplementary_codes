
library(limma)
library(GSEABase)
library(GSVA)
library(reshape2)
library(ggplot2)

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)

data=t(gsvaResult)

data=avereps(data)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-1),drop=F]

outTab=data.frame()
for(Geneset in colnames(data)){
	for(gene in colnames(risk)){
		x=as.numeric(data[,Geneset])
		y=as.numeric(risk[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Geneset=Geneset, cor, text, pvalue))
	}
}

outTab$Gene=factor(outTab$Gene, levels=colnames(risk))
outTab$cor=as.numeric(outTab$cor)
pdf(file="GSVAcor2.pdf", width=10, height=7)
ggplot(outTab, aes(Gene, Geneset)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#40E0D0", mid = "white", high = "#FF6347") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +    #È¥?ô±³¾?
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),   #x?????axis.text.y = element_text(size = 10, face = "bold")) +       #y??????
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #????Í¼x_discrete(position = "bottom")      #X?????()


######Vi