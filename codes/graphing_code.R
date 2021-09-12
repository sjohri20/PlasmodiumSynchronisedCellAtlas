setwd("C:/Users/HP/Desktop/sem 8/project/IG/plasmodium/simulation2_re")

data_clust = data.frame(matrix(nrow = 1, ncol=15))
data_tot = data.frame(matrix(nrow = 1, ncol=15))
obj = readRDS("C:/Users/HP/Desktop/sem 8/project/IG/plasmodium/synchronised plasmodium cell atlas/combined.rds")
barcodes = rownames(obj@meta.data[which(obj@meta.data[,7]==10),])

data_df = data.frame(matrix(nrow=50, ncol = 2))
for (i in c(1:100))
{
  set.seed(i)
  x = readRDS(paste(paste("C:/Users/HP/Desktop/sem 8/project/IG/plasmodium/simulation2_re/run_",i,sep=""),".rds",sep = ""))
  l1 = table(na.omit(x[c(barcodes),7]))
  l2 = table(na.omit(x[,7]))
  l3=data.frame(matrix(nrow = 1, ncol = length(l1)))
  for (k in c(1:length(l1)))
  {
    l3[1,k] = l1[k]/l2[k]
  }
  ind = which(l3[1,]==max(l3[1,]))
  indi=1
  if (length(ind)>1)
  {
    indi = ind[1]
  }
  else
  {
    indi = ind
  }
  data_df[i,1] = l1[indi]/sum(l1)
  data_df[i,2] = l1[indi]/l2[indi]
}

ggplot(data_df, aes(x=data_df[,1], y=data_df[,2]))+
  geom_point(pch=19, cex=2, alpha=0.3, position = "jitter")+
  # geom_text_repel(aes(label=paste("run_",rownames(data_df),sep="")))+
  labs(x="Largest percentage of Cell in Cluster", y="Percentage of treated cells in a cluster")

