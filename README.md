# code-for-prostate-cancer
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(data.table)
library(patchwork)
library(copykat)
library(monocle3)
library(monocle)
library(ggpubr)
library(CytoTRACE2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(BiocManager)
library(ComplexHeatmap)
library(ggunchull)
library(jjAnno)
library(scRNAtoolVis)
library(decoupleR)
library(scMetabolism)
library(CellChat)
library(NMF)
library(ggalluvial)
library(ClusterGVis)
library(tidyverse)
library(circlize)
library(cols4all)
library(harmony)
library(tibble)
sce.all <- CreateSeuratObject(counts = GSE143791tumor, min.cells = 5, min.features = 300)
sce.all <- Pcadt


#group----
table(sce.all$orig.ident)

sce.all@meta.data$Group<-ifelse(sce.all@meta.data$orig.ident %in% c("Pca01N","Pca04N","Pca06N"),"Normal",
                     ifelse(sce.all@meta.data$orig.ident %in% c("Pca01T1", "Pca02T", "Pca03TA", "Pca03TB", "Pca04T1", "Pca06T"),"PCa",
                            ifelse(sce.all@meta.data$orig.ident %in% c("P1", "P2", "P3", "P4", "P5", "P6"),"CRPC", "mCRPC")))


sce.all@meta.data$group<-ifelse(sce.all@meta.data$orig.ident %in% c("Pca01N","Pca04N","Pca06N"),"Normal",
                                ifelse(sce.all@meta.data$orig.ident %in% c("Pca01T1","Pca01T2","Pca02T", "Pca03TA", "Pca03TB", "Pca04T1", "Pca06T"),"PCa","CRPC"))


table(sce.all@meta.data$group)





#expressedmatrix----
expdt <- as.matrix(sce.all[["RNA"]]@counts)
colnames(expdt) <- gsub("-", "_", colnames(expdt))

expdt[1:4,1:4]

#copycat----
sc.cnv <- copykat(rawmat = expdt, ngene.chr = 5, win.size = 25, KS.cut = 0.2,
                  distance = "euclidean", n.cores = 60)

malignant <- read.delim("_copykat_prediction.txt", row.names=1)

pred.test <- data.frame(sc.cnv$prediction)
cna.test <- data.frame(sc.cnv$CNAmat)

sce.all@meta.data$copykat.pred <- pred.test$copykat.pred

tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]

tumor.mat <- cna.test[, which(colnames(cna.test) %in% tumor.cells)]

DimPlot(ca, group.by = "copykat.pred", reduction = "umap")

saveRDS(sc.cnv, "sc.cnvresults.RDS")


#降维----
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
sce.all <- subset(sce.all, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)



#normalization----
sce.all <- NormalizeData(sce.all,normalization.method = "LogNormalize", scale.factor = 10000)

#features selection----
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 3000)

top10 <- head(VariableFeatures(sce.all), 10)

plot3 <- VariableFeaturePlot(sce.all)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = T)

#均一化----

all.genes <- row.names(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)

#PCA----

sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))

print(sce.all[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(sce.all, dims = 1:2, reduction = "pca")

DimPlot(sce.all, reduction = "pca")

DimHeatmap(sce.all, dims = 1:15, cells = 500, balanced = T)

ElbowPlot(sce.all)

#不进行批次校正----
sce.all<-RunUMAP(sce.all,reduction="pca",dims=1:10,group.by.vars = "orig.ident")
DimPlot(sce.all,reduction="umap")

#进行批次校正----
sce.all <- RunHarmony(sce.all,reduction="pca",group.by.vars = "orig.ident")

sce.all<-RunUMAP(sce.all,reduction="harmony",dims=1:10)

DimPlot(sce.all,reduction="umap")

sce.all <- FindNeighbors(sce.all, dims=1:10)

seq <- seq(0.1,1, by=0.1)
for(res in seq){
  sce.all <- FindClusters(sce.all, resolution = res)
}
p11 <- clustree(sce.all, prefix = "RNA_snn_res.")+coord_flip()
p12 <- DimPlot(sce.all, group.by = "RNA_snn_res.0.3", label = T)
p11+p12

sce.all <- FindClusters(sce.all, resolution = 0.3)
head(Idents(sce.all))

sce.all <- RunUMAP(sce.all, dims = 1:10)

DimPlot(sce.all, reduction = "umap", label = T)

saveRDS(sce.all, file = "Pcadt.RDS")

sce.all <- CRPCdt

#美化----
pc12<-Embeddings(object=sce.all,reduction='umap')%>%data.frame()
lower<-floor(min(min(pc12$UMAP_1),min(pc12$UMAP_2)))-2
linelen<-abs(0.3*lower)+lower
mid<-abs(0.3*lower)/2+lower
axes<-data.frame(x=c(lower,lower,lower,linelen),y=c(lower,linelen,lower,lower),
                   group=c(1,1,2,2),
                   label=rep(c('UMAP_2','UMAP_1'),each=2))

label<-data.frame(lab=c('UMAP_2','UMAP_1'),angle=c(90,0),
                    x=c(lower-4,mid),y=c(mid,lower-1.5))
library(ggsci)
library(randomcoloR)
color_map <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
               "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
               "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
               "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
               "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
               "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")

palette <- distinctColorPalette(25)
DimPlot(sce.all,reduction='umap',label=T)+
  NoAxes()+NoLegend()+
  theme(aspect.ratio=1)+scale_color_discrete(type = palette)+
  geom_line(data=axes,
            aes(x=x,y=y,group=group),
            arrow=arrow(length=unit(0.1,"inches"),
                          ends="last",type="closed"))+
  geom_text(data=label,
            aes(x=x,y=y,angle=angle,label=lab),fontface='italic')


sce.all.markers <- FindAllMarkers(sce.all,only.pos = T)

#sce.all.markers.a <- FindMarkers(sce.all, ident.1 = "Club cell (Pca)", ident.2 = "Club cell (Nor)", min.pct = 0.25)
#sce.all.markers.b <- FindMarkers(sce.all, ident.1 = "Club cell (crpc)", ident.2 = "Club cell (Pca)", min.pct = 0.25)
#sce.all.markers.c <- FindMarkers(sce.all, ident.1 = "Club cell (crpc)", ident.2 = "Club cell (Nor)", min.pct = 0.25)

#write.table(sce.all.markers.a[sce.all.markers.a$p_val_adj<=0.05,], file = "./pcavsnor.markers.csv", sep = ",",col.names = T, row.names = T, quote = F)
#write.table(sce.all.markers.b[sce.all.markers.b$p_val_adj<=0.05,], file = "./crpcvspca.markers.csv", sep = ",",col.names = T, row.names = T, quote = F)
#write.table(sce.all.markers.c[sce.all.markers.c$p_val_adj<=0.05,], file = "./crpcvsnor.markers.csv", sep = ",",col.names = T, row.names = T, quote = F)


all.markers <- sce.all.markers |> group_by(cluster) 
write.table(all.markers[all.markers$p_val_adj<=0.05,], file = "./2216030.31.markers.csv", sep = ",",col.names = T, row.names = F, quote = F)

VlnPlot(ca, features = c("AR", "KRT8", "PSCA", "KLK3", "EPCAM", "KRT18", "KLK2", "ACPP"))

FeaturePlot(ca, features = c("AR"), cols = c("#eaf2f9","#e85a40"))

markers <- c("CD14", "CD68", #Myeloidcells
             "EPCAM", "KRT18", #Epithelial cells
             "ACTA2",  "CSF1R", "COL1A1","PDPN", #Fibroblast
             "VWF", "PECAM1", #Endothelialcells
             "SOX9", "PTPRC","CD3E","CD24",
             "CD3D"#STEMNESS,
)

DotPlot(sce.all, features = markers) + coord_flip()

new.cluster.ids <- c("Tumor", #0
                     "T cell", #1
                     "Tumor", #2
                     "Macrophage", #3
                     "Endothelial cell", #4
                     "Tumor", #5
                     "Epithelia", #6
                     "Fibroblast", #7
                     "Macrophage", #8
                     "Epithelia", #9
                     "B cell", #10
                     "B cell", #11
                     "Fibroblast",#12
                     "Mast cell", #13
                     "Fibroblast")

names(new.cluster.ids) <- levels(sce.all)

attributes(new.cluster.ids)

sce.all <- RenameIdents(sce.all, new.cluster.ids)

DimPlot(sce.all, reduction = "umap", label = F, pt.size = 0.5) 


sce.all$celltype <- Idents(sce.all)


#小提琴图美化----

markers <- c("MMP7", "OLFM4", "PIGR", "SCGB1A1", "SCGB3A1")
library(paletteer)

colormap<-c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))

p1<-VlnPlot(ca,features=markers,
              group.by="celltype",
              flip=T,stack=T,cols=colormap
)

p1+NoLegend()

p1+theme_bw()+
  theme(
    axis.text.x.bottom=element_text(angle=45,hjust=1,vjust=1),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.position="none",
    axis.text.x=element_text(color='black',size=11),
    axis.text.y=element_blank(),
    axis.title.x=element_text(color='black',size=15),
    axis.ticks.x=element_line(color='black'),
    axis.ticks.y=element_blank(),
  )


#heatmap----

top5=subset_MIA_IAC.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
g=unique(top5$gene)

p<-DotPlot(subset_MIA_IAC,features=g)
View(p)

dotplot_data<-p$data
class(dotplot_data)
heatmap_data<-dotplot_data%>%
  dplyr::select(features.plot,id,avg.exp.scaled)%>%
  pivot_wider(names_from=features.plot,values_from=avg.exp.scaled)

heatmap_data=column_to_rownames(heatmap_data,var="id")

#1----
pheatmap(heatmap_data,
         cluster_rows=F,
         cluster_cols=F,color=colorRampPalette(c("#c3ddf0","#eaf2f9","#e85a40"))(100))

annotation<-data.frame(CellType=row.names(heatmap_data))
row.names(annotation)<-row.names(heatmap_data)


#2----
pheatmap(heatmap_data,
         cellwidth=18,
         cellheight=24,
         cluster_rows=F,
         cluster_cols=F,
         color=colorRampPalette(c("#F6F5F1FF","#F9F9F9FF","#C8102EFF"))(100),
         annotation_row=annotation)


#3----
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
mean_gene_exp<-AverageExpression(subset.f.t.c,
                                   features=top5$gene,
                                   group.by='celltype',
                                   slot='data')%>%
  data.frame()%>%
  as.matrix()

colnames(mean_gene_exp)<- c("Epithelia", "Fibroblast", "NK cell", "T cell", "Macrophage", "Club cell", "Endothelial cell", "Mast cell", "Neuroendocrine cell", "B cell")

htdf<-t(scale(t(mean_gene_exp),scale=T,center=T))
anno_col<-brewer.pal(10,"Paired")
names(anno_col) <- c("Epithelia", "Fibroblast", "NK cell", "T cell", "Macrophage", "Club cell", "Endothelial cell", "Mast cell", "Neuroendocrine cell", "B cell")
col_fun=colorRamp2(c(-2,0,2),c("#0099CC","white","#CC0033"))

column_ha=HeatmapAnnotation(cluster=colnames(htdf),
                              col=list(cluster=anno_col))
Heatmap(htdf,
        name="Z-score",
        cluster_columns=F,cluster_rows=F,
        row_title="Top 5 Marker genes",
        column_title="Cell type",
        row_names_gp=gpar(fontface='italic',fontsize=10),
        row_names_side='left',
        border=T,
        rect_gp=gpar(col="white",lwd=1),
        column_names_side='top',
        column_names_rot=45,
        top_annotation=column_ha,
        #column_split=paste('clsuter',0:8,sep=''),
        col=col_fun)


#4----
AverageHeatmap(object=subset_T,
               markerGene=top5$gene)




#percentage----

PropPlot<-function(object,groupBy){
  #(1)获取绘图数据
  plot_data=object@meta.data%>%
    dplyr::select(orig.ident,{{groupBy}})%>%
    dplyr::rename(group=as.name(groupBy))
  
  #(2)绘图
  figure=ggbarstats(data=plot_data,
                      x=group,y=orig.ident,
                      package='ggsci',
                      palette='category20c_d3',
                      results.subtitle=FALSE,
                      bf.message=FALSE,
                      proportion.test=FALSE,
                      label.args=list(size=2,
                                        fill='white',
                                        alpha=0.85,
                                        family='Arial',
                                        fontface='bold'),
                      perc.k=2,
                      title='',
                      xlab='',
                      legend.title='Seurat Cluster',
                      ggtheme=ggpubr::theme_pubclean())+
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(color='black',lineend='round'),
          legend.position='right',
          axis.text.x=element_text(size=15,color='black',family='Arial'),
          axis.text.y=element_text(size=15,color='black',family='Arial'),
          legend.text=element_text(family='Arial',size=10,color='black'),
          legend.title=element_text(family='Arial',size=13,color='black'))
  
  #(3)去除柱子下面的样本量标识：
  gginnards::delete_layers(x=figure,match_type='GeomText')
}
sce.all$celltype <- Idents(sce.all)
library(ggstatsplot)
source('PropPlot.R')
table(sce.all$seurat_clusters,sce.all$celltype)
PropPlot(object=sce.all,groupBy='celltype')+
  PropPlot(object=sce.all,groupBy='seurat_clusters')


head(sce.all@meta.data)
table(sce.all$celltype)
mynames<-table(sce.all$celltype)%>%names()
myratio<-table(sce.all$celltype)%>%as.numeric()
pielabel<-paste0(mynames,"(",round(myratio/sum(myratio)*100,2),"%)")

cols<-c('#E64A35','#4DBBD4','#01A187','#6BD66B','#3C5588',
         '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
         '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)
pie(myratio,labels=pielabel,
    radius=1.0,clockwise=F,
    main="Cell Percentage",col=cols)

library(plotrix)
pie3D(myratio,labels=pielabel,explode=0.1,
      main="Cell Proption",
      height=0.3,
      labelcex=1)






#tumor cells----
sce.all <- GSE137829cluster

ca <- sce.all[,sce.all@meta.data$celltype %in% c("Tumor")]

ca <- NormalizeData(ca,normalization.method = "LogNormalize", scale.factor = 10000)

ca <- FindVariableFeatures(ca, selection.method = "vst", nfeatures = 2000)

all.genes.ca <- row.names(ca)
ca <- ScaleData(ca, features = all.genes.ca)

ca <- RunPCA(ca, features = VariableFeatures(object = ca))

ca <-  FindNeighbors(ca, dims = 1:10)

seq <- seq(0.1,1, by=0.1)
for(res in seq){
  ca <- FindClusters(ca, resolution = res)
}
p11 <- clustree(ca, prefix = "RNA_snn_res.")+coord_flip()
p12 <- DimPlot(ca, group.by = "RNA_snn_res.0.3", label = T)
p11+p12

ca <- FindClusters(ca, resolution = 0.1)

ca <- RunUMAP(ca, dims = 1:10)
DimPlot(ca, reduction = "umap")

ca.markers <- FindAllMarkers(ca, min.pct = 0.25)


ca.markers <- ca.markers |> group_by(cluster) 
write.table(ca.markers[ca.markers$p_val_adj<=0.05,], file = "./ca221603.markers.csv", sep = ",",col.names = T, row.names = F, quote = F)
ca.1 <- FindMarkers(sce.all, ident.1 = "C10", ident.2 = "Epithelial cell",in.pct = 0.25)
write.table(ca.1, file = "./ca.1.csv", sep = ",",col.names = T, row.names = T, quote = F)

VlnPlot(ca, features = c("MMP7", "OLFM4", "PIGR", "SCGB1A1", "SCGB3A1"))





new.cluster.ids <- c("C1", "C2", "C3", "C4", "C5")

names(new.cluster.ids) <- levels(ca)

attributes(new.cluster.ids)

ca <- RenameIdents(ca, new.cluster.ids)

DimPlot(ca, reduction = "umap", label = T, pt.size = 0.5) 


ca$celltype <- Idents(ca)


FeaturePlot(ca, features = c("HLA-DRA",
                             "HLA-DRB1",
                             "CD74",
                             "HLA-DRB5",
                             "CFB"
))





#美化-----
library(cols4all)
library(tidydr)
umapdata <- as.data.frame(ca@reductions$umap@cell.embeddings)

Cancer.cells.type <- Idents(ca)

umapdata <- cbind(umapdata, Cancer.cells.type)
head(umapdata)

mytheme <- theme_void() + theme(plot.margin = margin(5.5, 10, 5.5, 5.5))

pca <- ggplot(data = umapdata, aes(x = UMAP_1, y = UMAP_2)) + geom_point(aes(color=Cancer.cells.type), size=0.4, alpha=0.8)

pca2 <- pca+stat_ellipse(aes(color=Cancer.cells.type), level = 0.95, linetype=2, show.legend = F) + mytheme


pca3 <- pca2+theme_dr(xlength=0.2, ylength=0.2, arrow= grid::arrow(length=unit(0.1, "inches"),
                                                                   ends="last", type = "closed"))+
                                              theme(panel.grid = element_blank())

mycol <- c( "#F6CF71","#66C55C", "#F89C74", "#DCB0F2", "#9EB9F3", "#C9DB74","#7b34c1" ,"#0cf29a")


pca4 <- pca3 + scale_color_manual(values = mycol) + scale_fill_manual(values = mycol) 


ca.markers <- FindAllMarkers(ca, only.pos = T)

ca.markers.0 <- FindMarkers(ca, ident.1 = "C1", min.pct = 0.25)
ca.markers.1 <- FindMarkers(ca, ident.1 = "C2", min.pct = 0.25)
ca.markers.2 <- FindMarkers(ca, ident.1 = "C3", min.pct = 0.25)
ca.markers.3 <- FindMarkers(ca, ident.1 = "C4", min.pct = 0.25)
ca.markers.4 <- FindMarkers(ca, ident.1 = "C5", min.pct = 0.25)
ca.markers.5 <- FindMarkers(ca, ident.1 = "C6", min.pct = 0.25)

DotPlot(ca, features = unique(genes.check),)
VlnPlot(ca, features = c("RAC1", "CDC42", "TFAM"))

#gokegg1----

top200 <- ca.markers %>% group_by(cluster) %>% top_n(n=200, wt = avg_log2FC)
top200_gene <- bitr(unique(top200$gene), "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")
write.table(top200_gene, "geneid.txt", sep="\t")
markers <- full_join(top200, top200_gene, by=c("gene"="SYMBOL"))

keggpolot <- compareCluster(ENTREZID~cluster, data = markers, fun = "enrichKEGG")

dotplot(keggpolot, showCategory=40, label_format=40) 

go.results <- enrichGO(markers$ENTREZID, keyType = "ENTREZID", ont = "ALL", OrgDb = "org.Hs.eg.db", readable = T)
barplot(go.results, showCategory=10, label_format=40)


mygene <- c("SPINK1", "SCNN1B", 
            "HIST1H4C", "UBE2C", "KIAA0101",
            "SLPI",  "S100A8",
            "C15orf48", "GDF15",  
            "IGFBP3", "HILPDA",
            "TYROBP", "CCL17", "HSPA1A")
jjVolcano(diffData = ca.markers, myMarkers = mygene, log2FC.cutoff = 0.3, pSize = 1, tile.col = corrplot::COL2("RdBu", 15)[4:12],
          topGeneN = 3, base_size = 14, fontface="italic")


#GOKEGG2----

ca.markers.5 <- ca.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)

gvisdata <- prepareDataFromscRNA(object = ca, diffData = ca.markers.5, showAverage = T)

enrichdf <- enrichCluster(object = gvisdata, OrgDb = org.Hs.eg.db, type = "KEGG",
                          organism = "hsa", pvalueCutoff = 0.05, topn = 8, seed = 123123)

head(enrichdf,3)

markgenes <- unique(ca.markers.5$gene)[sample(1:length(unique(ca.markers.5$gene)), 40, replace = T)]

visCluster(object = gvisdata, plot.type = "line")

visCluster(object = gvisdata, plot.type = "heatmap", column_names_rot=45, markGenes = markgenes, cluster.order = c(1:9))

pdf("keggheat.pdf", height = 9, width = 14, onefile = F)
visCluster(object = gvisdata, plot.type = "both", column_names_rot=45,
           show_row_dend = F, markGenes = markgenes, markGenes.side = "left",
           annoTerm.data = enrichdf, line.side = "left", cluster.order = c(1:6),
           go.col = rep(jjAnno::useMyCol("stallion", n=5), each=5),
           add.bar = T)
dev.off()

#gokegg3----

ca.markers.1.1 <- ca.markers %>% filter(cluster=="2")

giddf <- bitr(unique(ca.markers.1.1$gene), "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")

markers.ca <- full_join(ca.markers.1.1, giddf, by=c("gene"="SYMBOL"))

GO <- enrichGO(markers.ca$ENTREZID,  OrgDb = "org.Hs.eg.db", keyType ="ENTREZID",
               ont = "ALL", pvalueCutoff = 0.01, pAdjustMethod = "BH",
               qvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500, readable = T)

barplot(GO, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scales = "free") + 
  scale_y_discrete(labels= function(x)stringr::str_wrap(x, width=80))

KEGG <- enrichKEGG(markers.ca$ENTREZID,  keyType ="kegg", organism = "hsa",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH",
                     qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200, use_internal_data = F)

pdf("C1keggheat.pdf", height = 9, width = 10, onefile = F)
barplot(KEGG, showCategory = 20, title = "KEGG pathway")
dev.off()

#gokegg4----

library(singleseqgset)

h.human<-msigdbr(species="Homo sapiens",category="H")
h.names<-unique(h.human$gs_name)

h.sets<-vector("list",length=length(h.names))
names(h.sets)<-h.names

for(i in names(h.sets)){
  h.sets[[i]]<-pull(h.human[h.human$gs_name==i,"gene_symbol"])
}

logfc.data<-logFC(cluster.ids=ca@meta.data$celltype,expr.mat=ca@assays$RNA@data)
names(logfc.data)

gse.res<-wmw_gsea(expr.mat=ca@assays$RNA@data,
                    cluster.cells=logfc.data[[1]],
                    log.fc.cluster=logfc.data[[2]],
                    gene.sets=h.sets)

names(gse.res)

res.stats<-gse.res[["GSEA_statistics"]]
res.pvals<-gse.res[["GSEA_p_values"]]

res.pvals<-apply(res.pvals,2,p.adjust,method="fdr")

res.stats[order(res.stats[,1],decreasing=TRUE)[1:20],]

res.pvals[order(res.stats[,1],decreasing=TRUE)[1:20],]

heatmap3(res.stats,Colv=NA,cexRow=1,cexCol=1,scale="row")
dev.off()
#heatmap美化----
mycol<- colorRamp2(c(-1.5, 0, 1.5), c("#0da9ce", "white", "#e74a32"))

cellwidth= 1
cellheight= 0.35
cn= dim(res.stats)[2]
rn= dim(res.stats)[1]
w= cellwidth*cn
h= cellheight*rn

Heatmap(as.matrix(res.stats),
        width= unit(w, "cm"),
        height= unit(h, "cm"),
        name= 'Gene Ratio', 
        col= mycol,
        cluster_columns= F,
        cluster_rows= F,
        column_names_side= c('top'), #列名置于顶部
        column_names_rot= 60,
        row_names_gp= gpar(fontsize = 7, fontface = 'italic'), #行名大小，改为斜体
        rect_gp= gpar(col = "white", lwd = 1.5), #格子描边颜色和粗细
        heatmap_legend_param= list(legend_height = unit(5, "cm"),
                                    grid_width= unit(0.4, "cm"),
                                    labels_gp= gpar(fontsize = 10))
                      )


#gsva5----

library(GSVA)
library(GSEABase)

table(Idents(ca))
expgsva=AverageExpression(ca)[[1]]
expgsva=as.matrix(expgsva)
expgsva=expgsva[rowSums(expgsva)>0,]
expgsva[1:4,1:4]

h_df=read.gmt("~/h.all.v2024.1.Hs.symbols.gmt")[,c(2,1)]
h_list=unstack(h_df)

GSVA <- gsva(expgsva, h_list)
GSVA[1:4,1:4]

pheatmap(GSVA,scale="row",angle_col="45",
         color=colorRampPalette(c("#0da9ce","white","#C8102EFF"))(50),
                                cluster_rows=F,
                                cluster_cols=F,cellwidth = 15, cellheight = 10, 
         heatmap_legend_param = list(title="GSVA score"),fontsize_row = 7)






#monocle----
cacluster <- ca[,ca@meta.data$seurat_clusters %in% c("0", "1", "2", "3")]
cacluster <- ca[,ca@meta.data$celltype %in% c("C1", "C2", "C3", "C4", "C5")]

saveRDS(cacluster, "cacluster.RDS")
expmatrix <- as(as.matrix(cacluster@assays$RNA@counts), "sparseMatrix")
pdata <- cacluster@meta.data
fdata <- data.frame(gene.short.name=row.names(cacluster), row.names = rownames(cacluster))
pd <- new("AnnotatedDataFrame", data=pdata)
fd <- new("AnnotatedDataFrame", data=fdata)

cds <- newCellDataSet(cellData = expmatrix, phenoData = pd, featureData = fd, 
                      expressionFamily = negbinomial.size())
pData(cds)$celltype=cacluster$celltype

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = .1) #过滤

ca.markers.df <- FindAllMarkers(ca, only.pos = T, min.pct = 0.20, logfc.threshold = 0.20)
ca.markers.df <- ca.markers.df[ca.markers.df$p_val_adj<0.05,]

cds <- setOrderingFilter(cds, ca.markers.df$gene) #轨迹定义基因


cds <- reduceDimension(cds, max_components = 2, method="DDRTree")


cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave("tree_Pseudotime.pdf", width = 6, height = 6)


plot_cell_trajectory(cds, color_by = "celltype")
ggsave("tree_Pseudotimeclusters.pdf", width = 6, height = 6)

plot_cell_trajectory(cds, color_by = "celltype")
ggsave("tree_PseudotimeState.pdf", width = 6, height = 6)


plot_cell_trajectory(cds, color_by = "seurat_clusters")+facet_wrap("~seurat_clusters", nrow = 1)

plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "celltype")


Time_diff <- differentialGeneTest(cds[ca.markers.df$gene,], cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[, c(0,1,2,3,4)]

Time.genes <- as.character(Time_diff$gene.short.name[1:30])

plot_pseudotime_heatmap(cds[Time.genes,], show_rownames = T, return_heatmap = T)


pData(cds)[, "CRYAB"]=cacluster@assays$RNA@scale.data["CRYAB",]

plot_cell_trajectory(cds, color_by = 'CRYAB') + scale_color_gsea()

colnames(pData(cds))
pData(cds)$CRYAB=log2(exprs(cds)['CRYAB',]+1)
plot_cell_trajectory(cds,color_by="CRYAB")+scale_color_gsea()

gs <- c("CRYAB", "HIF1A")

plot_genes_in_pseudotime(cds[gs,],color_by="celltype",ncol=1)


#cytotrace2----

expmatrixcy <- as.matrix(ca@assays$RNA@counts)

phenotype <- data.frame(phenotype=ca@active.ident)

names(phenotype) <- colnames(ca)
resultcyto <- cytotrace2(input = ca,  species = "human", is_seurat = T, slot_type = "counts")

plots <- plotData(cytotrace2_result = resultcyto, annotation = phenotype, is_seurat = T)

p1 <- plots$CytoTRACE2_UMAP

p2 <- plots$CytoTRACE2_Potency_UMAP

p3 <- plots$CytoTRACE2_Relative_UMAP

p4 <- plots$Phenotype_UMAP

p5 <- plots$CytoTRACE2_Boxplot_byPheno

plot_grid(p1,p2,p4, p5, nrow = 2)

#cellcycle----

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
caclustercellcycle <- CellCycleScoring(cacluster, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
                                       
                                       
VlnPlot(caclustercellcycle, features = c("S.Score", "G2M.Score"), group.by = "old.ident", pt.size = 0)                                       
                                       
                                       
                                       
#TFdecouple----

tfnetpro <- get_progeny(organism = "human", top = 500)
tfmat <- as.matrix(ca@assays$RNA@data)
tfacct <- run_mlm(mat = tfmat, network = tfnetpro, 
                  .source = "source", .target = "target", .mor = "weight", minsize = 5)

ca[["pathwaysmlm"]] <- tfacct %>% pivot_wider(id_cols = "source", names_from = "condition",
                                                      values_from = "score") %>% column_to_rownames("source") %>% Seurat::CreateAssayObject(.)

DefaultAssay(object = ca) <- "pathwaysmlm"

ca <- ScaleData(ca)

ca@assays$pathwaysmlm@data <- ca@assays$pathwaysmlm@scale.data

P1 <- DimPlot(ca, reduction = "umap", label = T, pt.size = 0.5)

p2 <- (FeaturePlot(ca, features = c("HIF1A")) & scale_color_gradient2(low = "blue", mid = "grey", high = "red")) + ggtitle("HIF1A activity")

DefaultAssay(object = ca) <- "RNA"

P3 <- FeaturePlot(ca, features = c("HIF1A")) + ggtitle("HIF1A expression")

p2+P3

tfaccthot <- t(as.matrix(ca@assays$pathwaysmlm@data)) %>% as.data.frame() %>% 
  mutate(cluster=Idents(ca)) %>% 
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>% 
  dplyr::group_by(cluster, source) %>% dplyr::summarize(mean=mean(score))



top.acts.mat <- tfaccthot %>% pivot_wider(id_cols = "cluster", names_from = "source", values_from = "mean") %>% 
                column_to_rownames("cluster") %>% as.matrix()

palette_length=100

my_color <- colorRampPalette(c("#FDDFA4FF", "white", "#C2697FFF"))(palette_length)
my.breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

pheatmap(top.acts.mat, border_color = NA, color = my_color, breaks = my.breaks)



tfnetcolltri <- get_collectri(organism = "human", split_complexes = F)
tfmat <- as.matrix(ca@assays$RNA@data)

tfactcoll <- run_ulm(mat = tfmat, network = tfnetcolltri, 
                  .source = "source", .target = "target", .mor = "mor", minsize = 5)


ca[["tfsulm"]] <- tfactcoll %>% pivot_wider(id_cols = "source", names_from = "condition",
                                                   values_from = "score") %>% column_to_rownames("source") %>% Seurat::CreateAssayObject(.)

DefaultAssay(object = ca) <- "tfsulm"

ca <- ScaleData(ca)

ca@assays$tfsulm@data <- ca@assays$tfsulm@scale.data

p1 <- DimPlot(ca, reduction = "umap", label = T, pt.size = 0.5)

p2 <- (FeaturePlot(ca, features = c("AR")) & scale_color_gradient2(low = "#339933", mid = "white", high = "red")) + ggtitle("AR activity")

DefaultAssay(object = ca) <- "RNA"

p3 <- (FeaturePlot(ca, features = c("AR")) & scale_color_gradient2(low = "#339933", mid = "grey", high = "red"))+ ggtitle("AR expression")


n.tfs <- 50

tfactulm <- t(as.matrix(ca@assays$tfsulm@data)) %>% as.data.frame() %>%
  mutate(cluster=Idents(ca)) %>% 
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>% 
  dplyr::group_by(cluster, source) %>% dplyr::summarise(mean=mean(score))


tfalls <- tfactulm %>% dplyr::group_by(source) %>% dplyr::summarise(std=sd(mean)) %>% arrange(-abs(std)) %>% head(n.tfs) %>% pull(source)

top.acts.mat.ulm <- tfactulm %>% filter(source %in% tfalls) %>%
  pivot_wider(id_cols = "cluster", names_from = "source", values_from = "mean") %>% 
  column_to_rownames("cluster") %>% as.matrix()

palette_length=100

my_color <- colorRampPalette(c("#FDDFA4FF", "white", "#C2697FFF"))(palette_length)
my.breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

pheatmap(top.acts.mat.ulm, border_color = NA, color = my_color, breaks = my.breaks)





#metabolism----

scmt.se <- sc.metabolism.Seurat(obj = ca, method = "AUCell", imputation=F, ncores=8, metabolism.type = "KEGG")

scmt.se.matrix <- scmt.se@assays$METABOLISM$score

input.pathway <- rownames(scmt.se@assays[["METABOLISM"]][["score"]])[1:30]

DotPlot.metabolism(obj = scmt.se, pathway = input.pathway, phenotype = "celltype", norm = "y")


input.pathway <- rownames(scmt.se@assays[["METABOLISM"]][["score"]])[1:10]

countexp.Seurat <- scmt.se
BoxPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "celltype", ncol = 1)

DimPlot.metabolism(obj = countexp.Seurat, pathway = "Synthesis and degradation of ketone bodies", 
                   dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1)


#通路活性分析----
library(fgsea)
all_cell_types <- as.vector(ca$celltype)
cell_types <- unique(all_cell_types)
#data2 metabolism pathway gene set
pathway_file <- "./KEGG_metabolism_nc.gmt"
pathways <- gmtPathways(pathway_file)

ca$newcelltype <- paste0(ca$orig.ident,"_", ca$celltype)
all_cell_types <- as.vector(ca$newcelltype)
cell_types <- unique(all_cell_types)

metabolism_activaty <- Pathway_act_Score(ca,
                                         pathways=pathways,
                                         assay = "RNA",
                                         filterGene=T,
                                         Mean_cut = 0.001,
                                         percent_cut =0.1,
                                         all_cell_types = all_cell_types,
                                         cell_types = cell_types)







#亚群放回----


Idents(sce.all, cells=colnames(ca)) <- Idents(ca)
DimPlot(sce.all, pt.size = 1, label = F)

#cellchat----

cellchatdata <- subset(sce.all, idents=c("C2","T cell","NK cell","B cell", "Monocyte", "Myeloid cell", "Plasma cell", "Monocyte progenitor cell", "Dendritic cell"))
cellchatdata <- subset.f.t.c
data.input <- cellchatdata@assays$RNA@data
meta.datainout <- cellchatdata@meta.data
meta.datainout$celltype <- Idents(cellchatdata)

meta.datainout=meta.datainout[!is.na(meta.datainout$celltype),]
data.input=data.input[, row.names(meta.datainout)]

cellchat <- createCellChat(object = data.input, meta = meta.datainout, group.by = "celltype")
CellChatDB = CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use=CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(object = cellchat, raw.use = T)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
par(mfrow=c(1,2), xpd=T)
groupsize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupsize, 
                 weight.scale = T, label.edge = F, title.name = "Number of Interaction")

netVisual_circle(cellchat@net$weight, vertex.weight = groupsize, 
                 weight.scale = T, label.edge = F, title.name = "Interaction Weight")

mat <- cellchat@net$weight
par(mfrow=c(1,2), xpd=T)

for (i in 1:nrow(mat)){
  mat2 <- matrix(0, nrow=nrow(mat), ncol=ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2, vertex.weight = groupsize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


cellchat@netP$pathways

pathways.show <- "FN1"

table(cellchat@idents)

cellchat@meta$celltype %>% head()

vertex.receiver=c(1,2)

pathways.show <- "FN1"
# "MK", "CD99", MIF, APP, GALECTIN, VISFATIN, 
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver,layout = "hierarchy")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = pathways.show)

netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8:10), remove.isolate = F)

netVisual_chord_gene(cellchat, sources.use = 7, targets.use = c(1:14), lab.cex = 0.3, legend.pos.y = 30, legend.pos.x = 10)


plotGeneExpression(cellchat, signaling = "FN1")

cellchatnet <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

netAnalysis_signalingRole_network(cellchatnet, signaling = pathways.show, width = 8, height = 3, font.size = 10)

gg1 <- netAnalysis_signalingRole_scatter(cellchatnet)
gg2 <- netAnalysis_signalingRole_scatter(cellchatnet, signaling = "FN1")
gg1+gg2


hh1 <- netAnalysis_signalingRole_heatmap(cellchatnet, pattern = "outgoing", height = 14)
hh2 <- netAnalysis_signalingRole_heatmap(cellchatnet, pattern = "incoming", height = 14)
hh1+hh2

netAnalysis_signalingRole_heatmap(cellchatnet, signaling = "FN1")


#nichenetr----
library(nichenetr)

lr_network<-readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix<-readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks<-readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

lr_network<-lr_network %>% distinct(from,to)

receiver="T cell"
expressed_genes_receiver<-get_expressed_genes(receiver, sce.all, pct=0.05)

all_receptors<-unique(lr_network$to)
expressed_receptors<-intersect(all_receptors,expressed_genes_receiver)
potential_ligands<-lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes<-c("C2", "C3", "C4")

list_expressed_genes_sender<-sender_celltypes%>%unique()%>%lapply(get_expressed_genes,sce.all,0.05)
expressed_genes_sender<-list_expressed_genes_sender%>%unlist()%>%unique()
potential_ligands_focused<-intersect(potential_ligands,expressed_genes_sender)

length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)

condition_oi<-"Tumor"
condition_reference<-"Involved"

seurat_obj_receiver<-subset(sce.all,idents=receiver)
DE_table_receiver<-FindMarkers(object=seurat_obj_receiver,
                                  ident.1=condition_oi,ident.2=condition_reference,
                                  group.by="group",
                                  min.pct=0.05)%>%rownames_to_column("gene")

geneset_oi<-DE_table_receiver %>% filter(p_val_adj<=0.05 & abs(avg_log2FC)>=0.25) %>% pull(gene)
geneset_oi<-geneset_oi %>%.[.%in% rownames(ligand_target_matrix)]

background_expressed_genes<-expressed_genes_receiver%>%.[.%in% rownames(ligand_target_matrix)]


ligand_activities<-predict_ligand_activities(geneset=geneset_oi,
                                               background_expressed_genes=background_expressed_genes,
                                               ligand_target_matrix=ligand_target_matrix,
                                               potential_ligands=potential_ligands)
ligand_activities<-ligand_activities%>%arrange(-aupr_corrected)%>%mutate(rank=rank(-aupr_corrected))
ligand_activities


p_hist_lig_activity<-ggplot(ligand_activities,aes(x=aupr_corrected))+
  geom_histogram(color="black",fill="darkorange")+
  geom_vline(aes(xintercept=min(ligand_activities%>%top_n(30,aupr_corrected)%>%pull(aupr_corrected))),
             color="red",linetype="dashed",size=1)+
  labs(x="ligand activity (PCC)",y="#ligands")+
  theme_classic()

p_hist_lig_activity

best_upstream_ligands<-ligand_activities%>%top_n(100,aupr_corrected)%>%arrange(-aupr_corrected)%>%pull(test_ligand)

vis_ligand_aupr<-ligand_activities%>%dplyr::filter(test_ligand%in%best_upstream_ligands)%>%
  column_to_rownames("test_ligand")%>%dplyr::select(aupr_corrected)%>%arrange(aupr_corrected)%>%as.matrix(ncol=1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands","Ligand activity",
                     legend_title="AUPR",color="darkorange")+
    theme(axis.text.x.top=element_blank()))

active_ligand_target_links_df<-best_upstream_ligands%>%
  lapply(get_weighted_ligand_target_links,
         geneset=geneset_oi,
         ligand_target_matrix=ligand_target_matrix,
         n=100)%>%
  bind_rows()%>%drop_na()
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links<-prepare_ligand_target_visualization(
  ligand_target_df=active_ligand_target_links_df,
  ligand_target_matrix=ligand_target_matrix,
  cutoff=0.33)

nrow(active_ligand_target_links)
head(active_ligand_target_links)

order_ligands<-intersect(best_upstream_ligands,colnames(active_ligand_target_links))%>%rev()
order_targets<-active_ligand_target_links_df$target%>%unique()%>%intersect(rownames(active_ligand_target_links))

vis_ligand_target<-t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target,"Prioritized ligands","Predicted target genes",
                    color="purple",legend_title="Regulatory potential")+
  scale_fill_gradient2(low="whitesmoke",high="purple")

p_dotplot<-DotPlot(subset(sce.all,celltype %in% sender_celltypes),
                     features=rev(best_upstream_ligands),cols="RdYlBu")+
  coord_flip()+
  scale_y_discrete(position="right")

p_dotplot



#CommPath----
library(CommPath)

expr=as.matrix(subset.f.t.c@assays$RNA@data)
label=Idents(subset.f.t.c)


comm.obj<-createCommPath(expr.mat=expr,
                          cell.info=label,
                           species="hsapiens")

comm.obj<-findLRmarker(object=comm.obj,method="wilcox.test")
comm.obj<-findLRpairs(object=comm.obj,
                        logFC.thre=0,
                        p.thre=0.01)
circosPlot(object=comm.obj,filter=FALSE)
circosPlot(object=comm.obj,plot="intensity",filter=FALSE)
comm.obj<-findLRpath(object=comm.obj,category="kegg")

comm.obj<-scorePath(object=comm.obj,method="gsva",min.size=10,parallel.sz=40)
acti.path.dat<-diffAllPath(object=comm.obj,only.posi=TRUE,only.sig=TRUE)
head(acti.path.dat)
write.csv(acti.path.dat,"crpc02.acti.path.dat.csv")
pathHeatmap(object=comm.obj,col =c("#0da9ce", "white", "#e74a32"),
            acti.path.dat=acti.path.dat,pathway.label.size = 10,
            top.n.pathway=10,
            cell.aver=TRUE)


comm.obj<-filterLR(object=comm.obj,acti.path.dat=acti.path.dat)
acti.path.filtered.dat<-filterPath(object=comm.obj,acti.path.dat=acti.path.dat)
circosPlot(object=comm.obj,filter=TRUE)
circosPlot(object=comm.obj,plot="intensity",filter=TRUE)
select.ident='Exhausted T cell'
circosPlot(object=comm.obj,select.ident=select.ident,filter=TRUE)

select.ident="Exhausted T cell"
select.receptor="CTLA4"
ident.up.dat <- findLigand(object=comm.obj,
                           select.ident=select.ident,
                           select.receptor=select.receptor)
head(ident.up.dat)

select.ident="Exhausted T cell"
select.ligand="CD44"
ident.down.dat<-findReceptor(object=comm.obj,
                               select.ident=select.ident,
                               select.ligand=select.ligand)

dotPlot.LR(object=comm.obj,receptor.ident=select.ident)
dotPlot.LR(object=comm.obj,ligand.ident=select.ident)

comm.obj<-pathNet(object=comm.obj,acti.path.filtered.dat=acti.path.filtered.dat)
kup=comm.obj@pathway.net$upstream
kup$type=rep("upstream",length(rownames(kup)))
kdown=comm.obj@pathway.net$downstream
kdown$type=rep("downstream",length(rownames(kdown)))
write.csv(kup,"03.upstream.path.csv")
write.csv(kdown,"03.downstream.path.csv")

set.seed(1234)
pathNetPlot(object=comm.obj,select.ident=select.ident,plot="upstream",
            layout='layout.davidson.harel',
            vert.size.LR=3,vert.size.path.adj=10,
            LR.label='R',vertex.label.cex.LR=0.25,vertex.label.cex.path=0.3)

set.seed(1234)
pathNetPlot(object=comm.obj,select.ident=select.ident,plot="upstream",
            top.n.path=5,
            layout='layout.davidson.harel',
            vert.size.LR=3,vert.size.path.adj=10,
            LR.label='R',vertex.label.cex.LR=0.25,vertex.label.cex.path=0.3)


pathNetPlot(object=comm.obj,select.ident=select.ident,plot="upstream",#限定upstream
            select.path=c("Cholesterol metabolism","Herpessimplex virus 1 infection","Th1 and Th2 cell differentiation","Th17 cell differentiation","Neutrophil extracellular trap formation"),
            layout='layout.davidson.harel',
            vert.size.LR=3,vert.size.path.adj=10,
            LR.label='R',vertex.label.cex.LR=0.25,vertex.label.cex.path=0.3)

set.seed(1234)
pathNetPlot(object=comm.obj,select.ident=select.ident,plot="downstream", #downstream
            layout='layout.davidson.harel',
            vert.size.LR=3,vert.size.path.adj=10,
            LR.label='L',vertex.label.cex.LR=0.25,vertex.label.cex.path=0.3)

dotPlot.pathway(object=comm.obj,pathway=pathway,acti.path.filtered.dat=acti.path.filtered.dat,receptor.ident=select.ident,top.n.inter=10)


select.ident='Regulatory T cell'
pathPlot(object=comm.obj,
         select.ident=select.ident,top.n.receptor=30,
         acti.path.dat=acti.path.filtered.dat)

pathPlot(object=comm.obj,
         select.ident=select.ident,
         up.ident=c("T cell","NK cell","B cell","Plasma cell","Dendritic cell","Myeloid cell"),
         select.path=c("Huntington disease"),
         top.n.receptor=5,
         acti.path.dat=acti.path.filtered.dat)


pathPlot(object=comm.obj,
         select.ident=select.ident,
         down.ident=c("T cell","NK cell","B cell","Plasma cell","Dendritic cell","Myeloid cell"),
         select.path=c("PI3K-Akt signaling pathway","MAPK signaling pathway"),
         top.n.receptor=20,
         acti.path.dat=acti.path.filtered.dat)

select.ident='Regulatory T cell'
pathChainPlot(object=comm.obj,
              select.ident=select.ident,
              up.ident=c("iCAF","C3","B cell", "Cytotoxic T cell",
                         "Exhausted T cell", "myCAF", "Mast cell", "Macrophage"),
              down.ident=c("iCAF","C3","B cell", "Cytotoxic T cell",
                           "Exhausted T cell", "Macrophage", "Mast cell", "myCAF"),
              #select.path=c("MAPK signaling pathway","HIF-1 signaling pathway","AMPK signaling pathway"),
              top.n.path=30,
              acti.path.dat=acti.path.filtered.dat)





#Tcell-----
subset_T <- sce.all[,sce.all@meta.data$celltype %in% c("T cell")]
subset_T[["percent.mt"]] <- PercentageFeatureSet(subset_T, pattern = "^MT-")
VlnPlot(subset_T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(subset_T, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(subset_T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
subset_T <- subset(subset_T, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 15)


#normalization----
subset_T <- NormalizeData(subset_T,normalization.method = "LogNormalize", scale.factor = 10000)

#features selection----
subset_T <- FindVariableFeatures(subset_T, selection.method = "vst", nfeatures = 3000)

#均一化----

all.genes.T <- row.names(subset_T)
subset_T <- ScaleData(subset_T, features = all.genes.T)

#PCA----

subset_T <- RunPCA(subset_T, features = VariableFeatures(object = subset_T))

print(subset_T[["pca"]], dims = 1:5, nfeatures = 5)

subset_T <- RunHarmony(subset_T,reduction="pca",group.by.vars = "orig.ident")

subset_T<-RunUMAP(subset_T,reduction="pca",dims=1:10)

DimPlot(subset_T,reduction="umap")

subset_T <- FindNeighbors(subset_T, dims=1:10)

seq <- seq(0.1,1, by=0.1)
for(res in seq){
  subset_T <- FindClusters(subset_T, resolution = res)
}
p11 <- clustree(subset_T, prefix = "RNA_snn_res.")+coord_flip()
p12 <- DimPlot(subset_T, group.by = "RNA_snn_res.0.3", label = T)
p11+p12

subset_T <- FindClusters(subset_T, resolution = 0.5)
head(Idents(subset_T))

subset_T <- RunUMAP(subset_T, dims = 1:10)
DimPlot(subset_T, reduction = "umap", label = T)

subset_T.markers <- FindAllMarkers(subset_T, only.pos = T)

VlnPlot(subset_T, features = c("AR", "KRT8", "PSCA", "KLK3", "EPCAM", "KRT18", "KLK2", "ACPP"))

subset_T.all.markers <- subset_T.markers |> group_by(cluster) 
write.table(subset_T.all.markers[subset_T.all.markers$p_val_adj<=0.05,], file = "./Tttt.markers.csv", sep = ",",col.names = T, row.names = F, quote = F)

new.cluster.ids.t <- c("Naive T cell", #0
                     "Regulatory T cell", #1
                     "Cytotoxic T cell", #2
                     "Cytotoxic T cell", #3
                     "Regulatory T cell", #4
                     "Exhausted T cell", #5
                     "Cytotoxic T cell", #6
                     "Cytotoxic T cell", #7
                     "Cytotoxic T cell", #8
                     "Regulatory T cell")


names(new.cluster.ids.t) <- levels(subset_T)

attributes(new.cluster.ids.t)

subset_T <- RenameIdents(subset_T, new.cluster.ids.t)

DimPlot(subset_T, reduction = "umap", label = T, pt.size = 0.5) 

Idents(sce.all, cells=colnames(subset_T)) <- Idents(subset_T)
DimPlot(sce.all, pt.size = 1, label = F)



subset_T$celltype <- Idents(subset_T)




#fibroblst----

subset_f <- sce.all[,sce.all@meta.data$celltype %in% c("Fibroblast")]
subset_f[["percent.mt"]] <- PercentageFeatureSet(subset_f, pattern = "^MT-")
VlnPlot(subset_f, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(subset_f, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(subset_f, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
subset_f <- subset(subset_f, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)


#normalization----
subset_f <- NormalizeData(subset_f,normalization.method = "LogNormalize", scale.factor = 10000)

#features selection----
subset_f <- FindVariableFeatures(subset_f, selection.method = "vst", nfeatures = 3000)

#均一化----

all.genes.f <- row.names(subset_f)
subset_f <- ScaleData(subset_f, features = all.genes.f)

#PCA----

subset_f <- RunPCA(subset_f, features = VariableFeatures(object = subset_f))

print(subset_f[["pca"]], dims = 1:5, nfeatures = 5)

subset_f <- RunHarmony(subset_f,reduction="pca",group.by.vars = "orig.ident")

subset_f<-RunUMAP(subset_f,reduction="harmony",dims=1:10)

DimPlot(subset_f,reduction="umap")

subset_f <- FindNeighbors(subset_f, dims=1:10)

seq <- seq(0.1,1, by=0.1)
for(res in seq){
  subset_f <- FindClusters(subset_f, resolution = res)
}
p11 <- clustree(subset_f, prefix = "RNA_snn_res.")+coord_flip()
p12 <- DimPlot(subset_f, group.by = "RNA_snn_res.0.3", label = T)
p11+p12

subset_f <- FindClusters(subset_f, resolution = 0.2)
head(Idents(subset_f))

subset_f <- RunUMAP(subset_f, dims = 1:10)
DimPlot(subset_f, reduction = "umap", label = T)

subset_f.markers <- FindAllMarkers(subset_f, only.pos = T)

VlnPlot(subset_f, features = c("AR", "KRT8", "PSCA", "KLK3", "EPCAM", "KRT18", "KLK2", "ACPP"))

subset_f.all.markers <- subset_f.markers |> group_by(cluster) 
write.table(subset_f.all.markers[subset_f.all.markers$p_val_adj<=0.05,], file = "./f.markers.csv", sep = ",",col.names = T, row.names = F, quote = F)

new.cluster.ids <- c("myCAF", #0
                     "iCAF", #1
                     "iCAF", #2
                     "SMC", #3
                     "myCAF", #4
                     "iCAF", #5
                     "apCAF", #6
                     "apCAF", #7
                     "pCAF")


names(new.cluster.ids) <- levels(subset_f)

attributes(new.cluster.ids)

subset_f <- RenameIdents(subset_f, new.cluster.ids)

DimPlot(subset_f, reduction = "umap", label = T, pt.size = 0.5) 

subset_f$celltype <- Idents(subset_f)


Idents(sce.all, cells=colnames(subset_f)) <- Idents(subset_f)

DimPlot(subset.all, pt.size = 1, label = F)

sce.all$celltype <- Idents(sce.all)


subset.all <- sce.all[,sce.all@meta.data$celltype %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8",
                                                           "myCAF", "iCAF", "apCAF", "SMC", "pCAF", 
                                                           "Cytotoxic T cell","Naive T cell","Epithelia",                                                            "Regulatory T cell", 
                                                           "Exhausted T cell", 
                                                           "Macrophage", "Mast cell","Endothelial cell","B cell")]


subset.f.t.c <- sce.all[,sce.all@meta.data$celltype %in% c("C3","myCAF", "iCAF", "apCAF",
                                                           "Cytotoxic T cell", 
                                                           "Regulatory T cell", 
                                                           "Exhausted T cell", 
                                                           "Macrophage", "Mast cell","B cell")]
table(subset.f.t.c$celltype)



Idents(subset.f.t.c)






