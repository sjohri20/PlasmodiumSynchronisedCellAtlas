####### Monocle data prep #################

data_df_control = as.data.frame(reducedDims(cds)[["UMAP"]][which(tt=="_1"),])
data_df_treatment = as.data.frame(reducedDims(cds)[["UMAP"]][which(tt=="_2"),])

data_df_control["clust"] = cds@colData$seurat_clusters[which(tt=="_1")]
data_df_treatment["clust"] = cds@colData$seurat_clusters[which(tt=="_2")]

##################


data_df = data.frame("Total RNA of gene in Control Cells"= sum(df[input$var1_1,intersect(which(seurat_obj@meta.data[,6]==input$var1_2),which(tt=="_1"))]),
                     "Average RNA of gene per Control Cell" = sum(df[input$var1_1,intersect(which(seurat_obj@meta.data[,6]==input$var1_2),which(tt=="_1"))])/df_frequency_controls[as.integer(input$var1_2) + 1, 2])

df = as.data.frame(seurat_obj@assays$RNA@counts)
gene_1 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==0),])
gene_set_1 = as.data.frame(df[,gene_1])
gene_2 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==1),])
gene_set_2 = as.data.frame(df[,gene_2])
gene_3 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==2),])
gene_set_3 = as.data.frame(df[,gene_3])
gene_4 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==3),])
gene_set_4 = as.data.frame(df[,gene_4])
gene_5 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==4),])
gene_set_5 = as.data.frame(df[,gene_5])
gene_6 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==5),])
gene_set_6 = as.data.frame(df[,gene_6])
gene_7 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==6),])
gene_set_7 = as.data.frame(df[,gene_7])
gene_8 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==7),])
gene_set_8 = as.data.frame(df[,gene_8])
gene_9 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==8),])
gene_set_9 = as.data.frame(df[,gene_9])
gene_10 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==9),])
gene_set_10 = as.data.frame(df[,gene_10])
gene_11 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==10),])
gene_set_11 = as.data.frame(df[,gene_11])
gene_12 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[,5]==11),])
gene_set_12 = as.data.frame(df[,gene_12])

final_data = as.data.frame(rowSums(gene_set_1))
final_data["Cluster1"] = rowSums(gene_set_2)
final_data["Cluster2"] = rowSums(gene_set_3)
final_data["Cluster3"] = rowSums(gene_set_4)
final_data["Cluster4"] = rowSums(gene_set_5)
final_data["Cluster5"] = rowSums(gene_set_6)
final_data["Cluster6"] = rowSums(gene_set_7)
final_data["Cluster7"] = rowSums(gene_set_8)
final_data["Cluster8"] = rowSums(gene_set_9)
final_data["Cluster9"] = rowSums(gene_set_10)
final_data["Cluster10"] = rowSums(gene_set_11)
final_data["Cluster11"] = rowSums(gene_set_12)

colnames(final_data)[1]="Cluster0"

data_df = as.data.frame(rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==0)]))
data_df["Cluster1"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==1)])
data_df["Cluster2"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==2)])
data_df["Cluster3"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==3)])
data_df["Cluster4"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==4)])
data_df["Cluster5"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==5)])
data_df["Cluster6"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==6)])
data_df["Cluster7"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==7)])
data_df["Cluster8"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==8)])
data_df["Cluster9"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==9)])
data_df["Cluster10"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==10)])
data_df["Cluster11"] = rowSums(cds@assays@data@listData$counts[,which(monocle_data[,3]==11)])


######## Generate Gene Counts #################

df = seurat_obj@assays$RNA@counts

final_df = data.frame(matrix(nrow=5224, ncol=17))

gene_1 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_1"),7]==0),])
gene_set_1 = as.data.frame(df[,gene_1])
final_df[,1] = rowSums(gene_set_1)
tgene_1 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_2"),7]==0),])
tgene_set_1 = as.data.frame(df[,tgene_1])
final_df[,2] = rowSums(tgene_set_1)

gene_2 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_1"),7]==1),])
gene_set_2 = as.data.frame(df[,gene_2])
final_df[,3] = rowSums(gene_set_2)
tgene_2 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_2"),7]==1),])
tgene_set_2 = as.data.frame(df[,tgene_2])
final_df[,4] = rowSums(tgene_set_2)

gene_3 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_1"),7]==2),])
gene_set_3 = as.data.frame(df[,gene_3])
final_df[,5] = rowSums(gene_set_3)
tgene_3 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_2"),7]==2),])
tgene_set_3 = as.data.frame(df[,tgene_3])
final_df[,6] = rowSums(tgene_set_3)

gene_4 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_1"),7]==3),])
gene_set_4 = as.data.frame(df[,gene_4])
final_df[,7] = rowSums(gene_set_4)
tgene_4 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_2"),7]==3),])
tgene_set_4 = as.data.frame(df[,tgene_4])
final_df[,8] = rowSums(tgene_set_4)

gene_5 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_1"),7]==4),])
gene_set_5 = as.data.frame(df[,gene_5])
final_df[,9] = rowSums(gene_set_5)
tgene_5 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_2"),7]==4),])
tgene_set_5 = as.data.frame(df[,tgene_5])
final_df[,10] = rowSums(tgene_set_5)

gene_6 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_1"),7]==5),])
gene_set_6 = as.data.frame(df[,gene_6])
final_df[,11] = rowSums(gene_set_6)
tgene_6 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_2"),7]==5),])
tgene_set_6 = as.data.frame(df[,tgene_6])
final_df[,12] = rowSums(tgene_set_6)

gene_7 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_1"),7]==6),])
gene_set_7 = as.data.frame(df[,gene_7])
final_df[,13] = rowSums(gene_set_7)
tgene_7 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_2"),7]==6),])
tgene_set_7 = as.data.frame(df[,tgene_7])
final_df[,14] = rowSums(tgene_set_7)

gene_8 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_1"),7]==7),])
gene_set_8 = as.data.frame(df[,gene_8])
final_df[,15] = rowSums(gene_set_8)
tgene_8 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_2"),7]==7),])
tgene_set_8 = as.data.frame(df[,tgene_8])
final_df[,16] = rowSums(tgene_set_8)

tgene_9 = rownames(seurat_obj@meta.data[which(seurat_obj@meta.data[which(tt=="_2"),7]==8),])
tgene_set_9 = as.data.frame(df[,tgene_9])
final_df[,17] = rowSums(tgene_set_9)

rownames(final_df)=rownames(seurat_obj@assays$RNA@counts)
colnames(final_df)=c("Cluster0 : Controls","Cluster0: Treated","Cluster1 : Controls","Cluster1: Treated",
                     "Cluster2 : Controls","Cluster2: Treated","Cluster3 : Controls","Cluster3: Treated",
                     "Cluster4 : Controls","Cluster4: Treated","Cluster5 : Controls","Cluster5: Treated",
                     "Cluster6 : Controls","Cluster6: Treated","Cluster7 : Controls","Cluster7: Treated",
                     "Cluster8 : Treated")


################

df1 = combined@reductions$umap@cell.embeddings
t2=lapply(rownames(df1),function(x){substr(x,17,18)})

mdf=as.data.frame(df1[which(t2==".1"),])
rownames(mdf) = lapply(rownames(mdf), function(x){substr(x,1,16)})


mdf["bulk"] = as.data.frame(mca_df[,3])
par(mfrow=c(3,1))
m2 = as.data.frame(df1[which(t2=="_1"),])
m1 = as.data.frame(df1[which(t2=="_3"),])
p1 = ggplot(data = m2,mapping = aes(x=UMAP_1,y=UMAP_2))+geom_point(size=1,aes(color="#48d1cc"))+ggtitle("Controls")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",panel.grid = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+ xlab("UMAP 1")+ ylab("UMAP 2")
p2 = ggplot(data = mdf,mapping = aes(x=UMAP_1,y=UMAP_2))+geom_point(size=1,aes(color=bulk)) + ggtitle("Malaria Cell Atlas")+
  theme(plot.title = element_text(hjust = 0.6),panel.grid = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("UMAP 1")+ ylab("UMAP 2")
p3 = ggplot(data = m1,mapping = aes(x=UMAP_1,y=UMAP_2))+geom_point(size=1,aes(color="#48d1cc"))+ ggtitle("Treatment")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",panel.grid = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+ xlab("UMAP 1")+ ylab("UMAP 2")
cowplot::plot_grid(p1,p2,p3)

df2 = combined@meta.data["seurat_clusters"]
rn = rownames(df2)

cdf = as.data.frame(df1[which(t2=="_1"),])
rownames(cdf) = lapply(rownames(cdf), function(x){substr(x,1,16)})
tt=substr(rownames(df2),17,18)

trn = lapply(rownames(rn[which(tt=="_2")]),function(x){substr(x,1,16)})
crn = lapply(rownames(rn[which(tt=="_1")]),function(x){substr(x,1,16)})

cdf["seurat_clusters"] = df2[which(tt=="_1"),1]

tdf = as.data.frame(df1[which(t2=="_3"),])
rownames(tdf) = lapply(rownames(tdf), function(x){substr(x,1,16)})

for (i in c(1:nrow(tdf)))
{
  tdf[i,3] = df2[c(paste(rownames(tdf)[i],"_2", sep = "")),]
}
colnames(tdf)[3] = "seurat_clusters"

for (i in c(1:nrow(cdf)))
{
  cdf[i,3] = df2[c(paste(rownames(cdf)[i],"_1", sep = "")),]
}

for(i in c(1:nrow(cdf)))
{
  if(cdf[i,3]==0)
    cdf[i,4] = 0
  else if (cdf[i,3]==1)
    cdf[i,4] = 1
  else if (cdf[i,3]==4)
    cdf[i,4] = 2
  else if (cdf[i,3]==2)
    cdf[i,4] = 3
  else if (cdf[i,3]==3)
    cdf[i,4] = 4
  else if (cdf[i,3]==7)
    cdf[i,4] = 5
  else if (cdf[i,3]==5)
    cdf[i,4] = 6
  else if (cdf[i,3]==6)
    cdf[i,4] = 7
}

colnames(cdf)[4] = "s_clusters"

xdf = cdf
xdf = as.data.frame(xdf[,-3])

for (i in c(1:nrow(tdf)))
{
  if(tdf[i,3]==0)
    tdf[i,4] = 0
  else if (tdf[i,3]==1)
    tdf[i,4] = 1
  else if (tdf[i,3]==4)
    tdf[i,4] = 2
  else if (tdf[i,3]==2)
    tdf[i,4] = 3
  else if (tdf[i,3]==3)
    tdf[i,4] = 4
  else if (tdf[i,3]==7)
    tdf[i,4] = 5
  else if (tdf[i,3]==5)
    tdf[i,4] = 6
  else if (tdf[i,3]==6)
    tdf[i,4] = 7
  else if (tdf[i,3]==8)
    tdf[i,4] = 8
}

colnames(tdf)[4] = "s_clusters"

p1 <- ggplot(data = xdf,mapping = aes(x=UMAP_1,y=UMAP_2,color=as.factor(s_clusters)))+
  geom_point(size=1)+
  scale_color_manual(values = c("#B284BE","#00308F","#7CB9E8","#00CC22","#FFDD33","#FF8C19", "#FF0000","#C46210"))+
  ggtitle("Controls")+
  theme(plot.title = element_text(hjust = 0.6),panel.grid = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_fill_discrete(name="Cluster")

p2 <- ggplot(data = tdf,mapping = aes(x=UMAP_1,y=UMAP_2,color=as.factor(s_clusters)))+
  geom_point(size=1)+
  scale_color_manual(values = c("#B284BE","#00308F","#7CB9E8","#00CC22","#FFDD33","#FF8C19", "#FF0000","#C46210","#000000"))+
  ggtitle("Treatment")+
  theme(plot.title = element_text(hjust = 0.6),panel.grid = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


plot_grid(p1,p2)

p3 <- ggplot(data = mdf,mapping = aes(x=UMAP_1,y=UMAP_2))+
  geom_point(size=2, color="grey", alpha=0.3)+
  ggtitle("Controls")+
  theme(plot.title = element_text(hjust = 0.6),panel.grid = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = cdf, mapping=aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters), size=1)+
  scale_color_manual(values = c("#990099","#00308F","#3399FF","#00CC22","#FFF700","#FF8C19", "#FF0000","#C46210"))+
  xlab("UMAP 1")+ylab('UMAP 2')

p4 <- ggplot(data = mdf,mapping = aes(x=UMAP_1,y=UMAP_2))+
  geom_point(size=2, color="grey", alpha=0.3)+
  ggtitle("Treatment")+
  theme(plot.title = element_text(hjust = 0.6),panel.grid = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = tdf, mapping=aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters), size=1)+
  scale_color_manual(values = c("#990099","#00308F","#3399FF","#00CC22","#FFF700","#FF8C19", "#FF0000","#C46210","#000000"))+
  xlab("UMAP 1")+ylab("UMAP 2")+labs("Clusters")
