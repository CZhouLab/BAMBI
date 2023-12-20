.libPaths("/home/peng.zhou-umw/R_packeges") 
library(limma)
args = commandArgs(trailingOnly=TRUE)
matrix_input_file=args[1]
label_input_file = args[2]
output_file = args[3]

## matrix_input_file="/nl/umw_chan_zhou/Billy/Biomarker_Detection_Paper/Gentles2015_limma_NewPipeline/output/summary_step1.txt.microarray"
## label_input_file ="/nl/umw_chan_zhou/Billy/Biomarker_Detection_Paper/Gentles2015_limma_NewPipeline/output/summary_step1.txt.microarray.condition.csv"

data.matrix <- read.table(matrix_input_file, header=TRUE, sep = "\t", row.names =1)

data.label <- read.table(label_input_file, header=TRUE, sep = "\t")

data.label <- data.label[,ncol(data.label)]

targets <- data.frame(Cy3 = rep("null",length = length(data.label)),Cy5 = data.label)
rownames(targets) <- rownames(df)
design <- modelMatrix(targets,ref = "null")
colnames(design)  <- c("group1","group2")

fit  <- lmFit(data.matrix,design)
contrast.matrix <- makeContrasts(group2-group1, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

DE   <- topTable(fit2,adjust = "BH",coef = 1, number = nrow(data.matrix))

output_DE <- data.frame(baseMean = rep("NA",length = nrow(DE)),log2FoldChange = DE["logFC"], lfcSE = rep("NA",length = nrow(DE)), stat = rep("NA",length = nrow(DE)), pvalue = DE["P.Value"], padj = DE["adj.P.Val"])

colnames(output_DE) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

write.table(output_DE,output_file,sep="\t",quote=F)


