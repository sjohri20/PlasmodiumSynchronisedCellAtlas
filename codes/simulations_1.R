library(Seurat)

options(stringsAsFactors=F)
setwd("C:/Users/HP/Desktop/sem 8/project/IG/plasmodium/synchronised plasmodium cell atlas")
seurat_obj <- readRDS("combined.rds")

obj.t = subset(x = seurat_obj,subset = condition=="treatment")
obj.c = subset(x = seurat_obj,subset = condition=="control")

original.c = table(obj.c@meta.data[,7])
original.t = table(obj.t@meta.data[,7])

sampled.c = data.frame(matrix(ncol=9, nrow=5))
sampled.t = data.frame(matrix(ncol=9, nrow=5))

for (i in c(1:100000))
{
  df = obj.c@meta.data[sample(nrow(obj.c@meta.data),3375),]
  sampled.c[i,]=table(df[,7])
  df1 = obj.t@meta.data[sample(nrow(obj.t@meta.data),4700),]
  sampled.t[i,]=table(df1[,7])
}

prob.t = data.frame(matrix(ncol=5, nrow = 1))
for (i in c(1:9))
{
  prob.t[i,1] = max(sampled.t[,i])
  prob.t[i,2] = min(sampled.t[,i])
  prob.t[i,3] = mean(sampled.t[,i])
  prob.t[i,4] = median(sampled.t[,i])
  prob.t[i,5] = sd(sampled.t[,i])
}
colnames(prob.t) = c("Max","Min","Mean","Median", "Stdev")

hist(sampled.t[,9],main = "Size distribution of cluster unique to treated cells", 
     xlab = "No of cells in cluster", ylab = "Freqeuncy of occurence", ylim = c(0,30000), xlim = c(30,65),
     labels = F, col="gray")




