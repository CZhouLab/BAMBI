library(edgeR)
args = commandArgs(trailingOnly=TRUE)
Readcount_input=args[1]
Condtion_input=args[2]
output=args[3]

readCount<-read.table(file=Readcount_input, header = T, row.names = 1, stringsAsFactors = F,check.names = F)
conditions<-read.table(file=Condtion_input, header = F)
conditions<-factor(t(conditions))

y <- DGEList(counts=readCount,group=conditions)
##Remove rows conssitently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
##Perform TMM normalization and transfer to CPM (Counts Per Million)
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)

pvalues <- sapply(1:nrow(count_norm),function(i){
     data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
     p=wilcox.test(gene~conditions, data)$p.value
     return(p)
   })
fdr=p.adjust(pvalues,method = "fdr")

conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

outRst<-data.frame(log2FoldChange=foldChanges, pvalue=pvalues, padj=fdr)
rownames(outRst)=rownames(count_norm)
## outRst=na.omit(outRst)
## fdrThres=0.05
write.table(outRst, file=output,sep="\t", quote=F,row.names = T,col.names = T)




