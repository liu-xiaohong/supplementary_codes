
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

tmb=read.table("", header=T, sep="\t", check.names=F, row.names=1)
	
risk=read.table("", header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)
data$TMB[data$TMB>quantile(data$TMB,0.99)]=quantile(data$TMB,0.99)

data$Risk=ifelse(data$Risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

boxplot=ggboxplot(data, x="Risk", y="TMB", color="Risk",
			      xlab="",
			      ylab="Tumor tmbation burden",
			      legend.title="",
			      palette = c("blue", "red"),
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	

pdf(file="riskTMB.pdf", width=5, height=4.5)
print(boxplot)
dev.off()


xlab="riskScore"
ylab="TMB"
x=as.numeric(data[,xlab])
y=as.numeric(data[,ylab])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
		xlab("Risk score") + ylab("Tumor tmbation burden")+
		geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))

pdf(file="cor.pdf", width=5.2, height=5)
print(p2)
dev.off()


