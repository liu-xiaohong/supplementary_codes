
library(limma)
library(parallel)
set.seed(12345)

expFile=" "     
setwd(" ")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]


data=t(data)

data=avereps(data)
data=t(data)

GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 

calcPhenotype(trainingExprData = GDSC2_Expr,    
              trainingPtype = GDSC2_Res,        
              testExprData = data,              
              batchCorrect = 'eb',              
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,     
              minNumSamples = 10,               
              printOutput = TRUE,               
              removeLowVaringGenesFrom = 'rawData')


