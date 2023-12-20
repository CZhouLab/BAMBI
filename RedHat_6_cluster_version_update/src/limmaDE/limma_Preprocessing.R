.libPaths("/home/peng.zhou-umw/R_packeges") 
library(limma)
library(impute)
args = commandArgs(trailingOnly=TRUE)
microarray_input_file=args[1]
logtransform = as.integer(args[2])
microarray_output_file=args[3]
microarray_output_label_file=args[4]


df <- read.table(microarray_input_file, header=TRUE, sep = "\t", row.names =1)
dim(df)

data.matrix = t(df)

if (logtransform == 1){
tem.numbers <- as.numeric(unlist(data.matrix))
Med <- median(tem.numbers,na.rm = T) 
if(Med > 16) data.matrix <- log2(data.matrix) ## judge the data -- log2 transformed or not
rm(Med)}

if(exists(".Random.seed")) rm(.Random.seed)
set.seed(7)
saved.state <- .Random.seed
na.length <- length(which(is.na(data.matrix)==T))
if(na.length > 0) data.matrix <- impute.knn(data.matrix)$data
rm(na.length)
data.matrix <- normalizeBetweenArrays(data.matrix)
output.matrix <- cbind(Gene = rownames(data.matrix), data.matrix)
write.table(output.matrix,microarray_output_file,sep="\t",quote=F, row.names=F)

# data.label = df[,ncol(df)]
# write.table(data.label,microarray_output_label_file,sep="\t",quote=F, row.names=F)

