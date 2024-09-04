

path<-getwd()

library(Seurat)
library(ggpubr)
library(dittoSeq)
library(dplyr)
library(dbplyr)

#Diet object (to reduce size)####
embryo_d <- DietSeurat(embryo, scale.data = FALSE, dimreducs = re)

saveRDS(embryo_d, file = "Embryo_IG.rds")


######Selected resolution (number of clusters 0.35) ####
Idents(embryo_d)<-'RNA_snn_res.0.35'
Idents(embryo_d)<-'labels' 

levels(ob@active.ident)
clus_col <- c('C0-Mesenchymal cells2'= '#ba6917', 'C01-Osteoblast and fibroblasts'= '#f76dde',
              'C02-Mesenchymal cells'='#f45328','C03-Endothelial'='#985fe0',
              'C04-Chondrocytes'='#ccced1', 'C05-Epithelial'='#04c5fd', 'C06-Fetal brain (CNS)'='red2',
              'C07-SmoothMuscle'='orange', 'C08-Blood progenitors'='#579316',
              'C09-Myocytes'='#a7baa8', 'C10-Glial cells'='#0cfa66',
              'C11-Hepatocytes'='#0ba4fb', 'C12-Erythrocytes'='#fc0387', 
              'C13-NeuralCrest'='#606980', 'C14-Lymphatic progenitors'='#f16915',
              'C15-Cardiomyocytes'='#aff9fc', 'C16-PNS'='#c8fa6a')
only_cols <- c('#ba6917', '#f76dde', '#f45328','#985fe0',  '#ccced1','#04c5fd','red2',
                  'orange', '#579316','#a1f115','#0cfa66',
               '#0ba4fb','#fb0b2c', '#606980','#f16915',
                  '#1ab228', '#c8fa6a')

clus_col_new <- c('C15-Mesenchymal cells2'= '#ba6917', 'C16-Osteoblast and fibroblasts'= '#f76dde',
                  'C14-Mesenchymal cells'='#f45328','C06-Endothelial'='#985fe0',
                  'C13-Chondrocytes'='#ccced1', 'C04-Epithelial'='#04c5fd',
                  'C01-Fetal brain (CNS)'='red2',
                  'C12-SmoothMuscle'='orange', 'C09-Blood progenitors'='#579316',
                  'C11-Myocytes'='#a7baa8', 'C03-Glial cells'='#0cfa66',
                  'C10-Hepatocytes'='#0ba4fb', 'C08-Erythrocytes'='#fc0387', 
                  'C0-NeuralCrest'='#606980', 'C05-Lymphatic progenitors'='#f16915',
                  'C07-Cardiomyocytes'='#aff9fc', 'C02-PNS'='#c8fa6a')

Bestholtz_palette <- c("#DEDAD6","#FEE392","#FEC44E","#FE9929","#ED6F11","#CC4C17","#993411","#65260C")

clus_id <- c('C0-Mesenchymal cells2','C01-Osteoblast and fibroblasts','C02-Mesenchymal cells','C03-Endothelial',
             'C04-Chondrocytes','C05-Epithelial','C06-Fetal brain (CNS)','C07-SmoothMuscle',
             'C08-Blood progenitors','C09-Myocytes','C10-Glial cells','C11-Hepatocytes', 'C12-Erythrocytes',
             'C13-NeuralCrest', 'C14-Lymphatic progenitors',
             'C15-Cardiomyocytes', 'C16-PNS')

new_clus_id <- c('C15-Mesenchymal cells2','C16-Osteoblast and fibroblasts','C14-Mesenchymal cells','C06-Endothelial',
             'C13-Chondrocytes','C04-Epithelial','C01-Fetal brain (CNS)','C12-SmoothMuscle',
             'C09-Blood progenitors','C11-Myocytes','C03-Glial cells','C10-Hepatocytes', 'C08-Erythrocytes',
             'C0-NeuralCrest', 'C05-Lymphatic progenitors',
             'C07-Cardiomyocytes', 'C02-PNS')



#Generate new name for labeled clusters######
Idents(ob)<-'RNA_snn_res.0.35'
names(clus_id) <- levels(ob)
ob <- RenameIdents(ob, clus_id)
clus_id <- ob@active.ident
ob$new_labels <- as.character(clus_id)

saveRDS(ob, file = "Labelled_IreneFeb2023.rds")

#Generate new name for labeled clusters to reorder ######
Idents(ob)<-'RNA_snn_res.0.35'
names(new_clus_id) <- levels(ob)
ob <- RenameIdents(ob, new_clus_id)
new_clus_id <- ob@active.ident
ob$new_labels_rb <- as.character(new_clus_id)

saveRDS(ob, file = "Labelled_IreneMarch2023.rds")



Idents(ob) <- 'new_labels_rb'

levels(ob) <- c('C0-NeuralCrest','C01-Fetal brain (CNS)',
                'C02-PNS','C03-Glial cells','C04-Epithelial',
                'C05-Lymphatic progenitors','C06-Endothelial',
                'C07-Cardiomyocytes','C08-Erythrocytes',
                'C09-Blood progenitors','C10-Hepatocytes',
                'C11-Myocytes','C12-SmoothMuscle',
                'C13-Chondrocytes','C14-Mesenchymal cells',
                'C15-Mesenchymal cells2','C16-Osteoblast and fibroblasts')



###Plot dimensions clusters #####
Idents(ob) <- 'new_labels'


DimPlot(ob, reduction = 'umap_.19.dims',
        group.by = 'new_labels',
        cols=clus_col,
        label = T, label.box = T, repel = T, pt.size =0.7) + ggtitle('scRNAseq Clusters (Rbpj floxed)')


DimPlot(ob, reduction = 'umap_.19.dims',
        group.by = 'new_labels', split.by = 'Condition',
        cols=clus_col,pt.size =1.1) + ggtitle('scRNAseq Clusters (Rbpj floxed)')

DimPlot(ob, reduction = 'umap_.19.dims',
        group.by = 'Condition', cols=c('#80e13d','#f83d14'),pt.size =0.9) + ggtitle('scRNAseq Clusters (Rbpj floxed)')



FeaturePlot(ob, features = 'Myh11', pt.size = 1, reduction = 'umap_.19.dims')+
  scale_colour_gradientn(colours = Bestholtz_palette)




dittoBarPlot(ob, var='new_labels', group.by = 'Condition',
             color.panel = clus_col, x.reorder = c(2,1)) + 
  ggtitle('Cluster relative proportion by cell-type')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5, vjust = 0, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

dittoBarPlot(ob, var='Condition', group.by = 'new_labels', color.panel = c('#f83d14', '#80e13d')) + 
  ggtitle('Contribution per cluster')+
  theme(axis.text.x =  element_text(angle=90, hjust = 1,  vjust = 0, face = "italic", size = 10), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

temp <- c('C07-SmoothMuscle'='orange')
DimPlot(ob, reduction = 'umap_.19.dims',
        group.by = 'new_labels',
        cols=temp,
        label = T, label.box = T, repel = T, pt.size =0.7) + ggtitle('scRNAseq Clusters (Rbpj floxed)')



DimPlot(ob, reduction = 'umap_.19.dims',
        group.by = 'new_labels') + ggtitle('scRNAseq Clusters (Rbpj floxed)')



#Positive markers per cluster####
Idents(ob) <- 'new_labels'

markers_res0.35 <- FindAllMarkers(ob, logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf)

Idents(ob) <- 'new_labels_rb'



#I ordered after this step
markers <- FindAllMarkers(ob, logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf)

write.csv(markers, file='UpRegulated_ClusterLabelledRB.csv')


Idents(ob)<- 'Condition'

Ctrl <- subset(ob, idents="RbpjWt")
Idents(Ctrl) <- 'new_labels'
markers_ctrl <- FindAllMarkers(Ctrl, logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf)


write.csv(markers_ctrl, file='UpRegulated_ClusterLabelled_OnlyCtrl.csv')

write.csv(markers_res0.35, file='UpRegulated_ClusterLabelled.csv')


# Select top genes per cluster ####
markers$order <- 1:nrow(markers)
markers$cluster <- gsub("CO-Mesenchymal cells2", "C0-Mesenchymal cells2", markers$cluster)
Top5 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = -order)
s <- unique(Top5$gene)
Top5 <- Top5[match(s, Top5$gene),]
Top5 <- Top5 %>%  arrange(cluster) %>% group_by(cluster) %>%top_n(n = 5, wt = -order)
Top5genes <- Top5$gene

markers$order <- 1:nrow(markers)
Top5 <- markers %>%group_by(cluster) %>% top_n(n = 15, wt = -order)
s <- unique(Top5$gene)
Top5 <- Top5[match(s, Top5$gene),]
Top5 <- Top5 %>%  arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes <-  Top5$gene

markers_ctrl$order <- 1:nrow(markers_ctrl)
Top5 <- markers_ctrl %>%group_by(cluster) %>% top_n(n = 10, wt = -order)
s <- unique(Top5$gene)
Top5 <- Top5[match(s, Top5$gene),]
Top5 <- Top5 %>%  arrange(cluster) %>% group_by(cluster) %>%top_n(n = 5, wt = -order)
Top5genes <- Top5$gene


####Gene markers per cell type in ALL populations ####
myo <- c('Neb', 'Myh3', 'Tpm2', 'Acta2') #Myocytes
card <- c('Myl2', 'Myocd', 'Hcn4', 'Ctnna3', 'Ryr2', 'Tbx20') #Cardiac progenitors
hep <- c('Afp','Alb', 'Apoa2', 'Afp29', 'Pik3c2g', 'Hoga1', 'TFEB') #Hepatocytes
cont <- c('Il1rapl2', 'Meox2', 'Tgfb2', 'Adamts9', 'Postn', 'Ror1')#connective tissue progenitors
condroste <- c('Runx2', 'Twist2', 'Prrx1')#Condrocytes and Osteoblasts         
jaw_tooth <- c('Sox9', 'Foxp2', 'Col2a1','Col9a1', 'Col11a1', 'Pax9') #Jaw adn tooth progenitors
neuronsall <- c('Prmt8', 'Gadd45g', 'Cdkn1c',
                'Btg2','Pax2', 'Slc6a5', 'Ntng1',
                'Car10', 'Nrn1', 'Slc17a6', 'Grem2',
                'Syt13', 'Shox2', 'Ptprr', 'Pcbp3')
megaka <- c('Pf4', 'Itgb3', 'Itga2b', 'Ppbp', 'Cd226')#Megakaryocytes
megaka2 <- c('Pf4', 'Ppbp', 'Ctla2a', 'Gp1bb', 'Gp9', 'Plek', 'Clec1b', 'Itga2b', 'Itgb3', 'Fermt3')
melanocytes <- c('Tyr','TRPM1' , 'Pmel')
erythro <- c('Snca', 'Hbb-bs', 'Abcb4', 'Slc4a1', 'Kel', 'HBB-y',' Hba-x', 'Hba-a1', 'Hbb-bh1')
intermesoderm <- c('Wt1', 'Mylk', 'Ednra') #Intermediate mesoderm
condrocytes <- c('Itga11', 'Atp1a2', 'Lamc3', 'Epha7') #Chondrocyte progenitor
osteoblast <- c(  'Ibsp',  'Ifitm5',  'Smpd3',  'Sgms2',  'Sp7',  'Panx3',  'Alpl',  'Satb2',  'Fabp3',  'Hck',  'Col22a1',  'Dlx3',  'Cfh',  'Slc8a3')
osteoblastembryo <- c('Col1a1', 'Camk1d', 'Rbm8a')
epithelial <- c('Epcam', 'Trp63', 'Grhl2')
chondosteo <- c('Runx2', 'Twist2', 'Prrx1')

notochordcells <- c('Shh', 'Slit1', 'Ntn1')
neutrophils <- c('Ngp', 'S100A8')
whitebloodcells <- c('Apoe', 'Lyz2', 'Selenop', 'Ptprc', 'Ly86', 'CTSS')
lymphatic <- c('Pdpn', 'Flt4', 'Lyve1')

colinergicneurons <- c('Slit2', 'Slit3', 'Chat')
astrocytes <- c('Nfia','Glast', 'Fabp7','Sox9')
glial <- c('Gfap', 'Nestin', 'Nf1a', 'Pax6', 'Pth2r')
microglia <- c('Cd11b', 'Iba1', 'C1qa', 'Cd68', 'Cd45', 'Cx3cr1')
oligodendrocityes <- c('Olig1', 'Olig2','Sox10', 'Cnp', 'O4', 'Apc', 'Galc','Mbp')
astrocytes <- c('Gfap','S100b', 'Bfabp','Aldh1l1')
periferalglia <- c('Sox9', 'Sox10', 'Ap2a', 'Pax3', 'Sox2', 'Foxd3', 'Bfabp', 'S100b')
radialglia <- c('Sox2', 'Pax6', 'Hes1', 'Hes5','Nes', 'Vim', 'Pth2r', 'Fabp7', 'Pax3', 'Fzd10' )
neuprog <- c('Prmt8', 'Gadd45g', 'Cdkn1c', 'Btg2') #Neural progenitors
radglia <- c('Pth2r', 'Fabp7', 'Pax3', 'Fzd10', 'Hes5') #Radial glia
lens <- c('Gja8', 'Cryba1', 'Cryaa')
motoneurons <- c('Isl1', 'Isl2', 'Olig2')
sensoryneurons <- c('Syt13', 'Shox2', 'Ptprr', 'Pcbp3')
inhneurons <- c('Pax2', 'Slc6a5') #Inhibitory neurons
neuraltube <- c('Foxb1', 'Scube2', 'Prtg')
limbmesenchyme <- c('Msx1', 'Fgf10', 'Wnt5a', 'Lmx1b')
granuleneurons <- c('Neurod2', 'Tiam2')




VlnPlot(ob, group.by = 'new_labels', features=c)

FeaturePlot(ob, reduction = 'umap_.19.dims', features='Ptprc')+
  scale_colour_gradientn(colours = Bestholtz_palette)

DimPlot(ob, reduction = 'umap_.19.dims', group.by = 'new_labels')

t <- subset(Top10, Top10$cluster=='0')

Idents(ob) <- 'new_labels'
Idents(ob) <- factor(x = Idents(ob), levels = sort(levels(ob))) # to order levels

df <- data.frame(
  x = rep(c(3,8,13,18,23, 28, 33,38,43,48,53, 58,63,68,73,78,83)),
  y = rep(c(0),1),
  z = factor(rep(1:17)),
  color = only_cols)

DotPlot(Ctrl, group.by = 'new_labels', features=Top5genes, col.min = 0, dot.scale = 5, scale = T)+ 
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 7),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=ob$new_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = only_cols)
 

library(dittoSeq)



ob$Condition <- factor(ob$Condition, levels= c("RbpjWt", "RbpjLOF"))






DimPlot(ob, reduction = 'umap_.19.dims', group.by = 'new_labels',cols=clus_col)+
  ggtitle('scRNAseq Clusters (Rbpj floxed)')
DimPlot(ob, reduction = 'umap_.19.dims', group.by = 'Condition', cols=c('#80e13d','#f83d14'))+
  ggtitle('scRNAseq Clusters (Rbpj floxed)')
 


#### To subset clusters and get higher resolution #####

Idents(ob) <- 'new_labels'
Ec <- subset(ob, idents="C03-Endothelial")
saveRDS(Ec, file = "EndothelialCells.rds")

Idents(ob) <- 'new_labels'
Cardio <- subset(ob, idents="C15-Cardiomyocytes")
saveRDS(Cardio, file = "CardiomyocytesCells.rds")

Idents(ob) <- 'new_labels'
BloodP <- subset(ob, idents="C08-Blood progenitors")
saveRDS(BloodP, file = "BloodProgenitorCells.rds")

Idents(ob) <- 'new_labels'
smc <- subset(ob, idents="C07-SmoothMuscle")
saveRDS(smc, file = "SmoothMuscleCells.rds")


Idents(ob) <- 'new_labels'
epi <- subset(ob, idents="C05-Epithelial")
saveRDS(Ec, file = "EpithelialCells.rds")

#With new label names
Idents(ob) <- 'new_labels_rb'
osfi <- subset(ob, idents="C16-Osteoblast and fibroblasts")
saveRDS(osfi, file = "C16_Osteoblast_fibroblasts.rds")

Idents(ob) <- 'new_labels_rb'
mesenc <- subset(ob, idents="C15-Mesenchymal cells2")
saveRDS(mesenc, file = "C15-Mesenchymal cells2.rds")

Idents(ob) <- 'new_labels_rb'
mesenc1 <- subset(ob, idents="C14-Mesenchymal cells")
saveRDS(mesenc1, file = "C14-Mesenchymal cells.rds")

Idents(ob) <- 'new_labels_rb'
chond <- subset(ob, idents="C13-Chondrocytes")
saveRDS(chond, file = "C13-Chondrocytes.rds")

Idents(ob) <- 'new_labels_rb'
myoc <- subset(ob, idents="C11-Myocytes")
saveRDS(myoc, file = "C11-Myocytes.rds")

Idents(ob) <- 'new_labels_rb'
glial <- subset(ob, idents="C03-Glial cells"  )
saveRDS(glial, file = "C03_Glial_cells.rds")

Idents(ob) <- 'new_labels_rb'
pns <- subset(ob, idents= "C02-PNS")
saveRDS(pns, file = "C02_pns_cells.rds")


Idents(ob) <- 'new_labels_rb'
neurl <- subset(ob, idents= "C0-NeuralCrest")
saveRDS(neurl, file = "C0-NeuralCrest.rds")


### Different reduction testing #####
dim <- c(10, 13, 16, 19, 22, 25, 30)
n <- as.data.frame(list(c(1:7)))

myplots <- vector('list', nrow(n))

for (i in 1:length(dim)){
  FindNeighbors(Ec, dims = 1:dim[i])
  Ec <- FindClusters(Ec, resolution = 0.3)
  Ec <- RunUMAP(Ec, dims = 1:dim[i])
  myplots[[i]] <- local({
    i <- i
    p1 <- DimPlot(Ec, label = TRUE, repel = T, label.box = TRUE) + ggtitle(paste(dim[i],"PCAs"))+ theme(plot.title = element_text(hjust = 0.5))+ NoLegend()
  })
  
}

p25 <- patchwork::wrap_plots(myplots, ncol = 2)

##########To do reduction and cluster analysis #############

# Epithelial cells
ElbowPlot(epi, ndims = 45)

epi <- FindNeighbors(epi, dims = 1:18)
epi <- FindClusters(epi, resolution = 0.25)
epi <- RunUMAP(epi, dims = 1:18, reduction.name = "epi_umap18")


Idents(epi) <- 'RNA_snn_res.0.25'
epi_label <- c('C00','C01','C02','C03')
epi_cols <- c('C00'='#f9d520','C01'='#6bc115','C02'='#f93420','C03'='#a912f9')

DimPlot(epi, reduction="epi_umap18", pt.size = 2, cols = epi_cols)

names(epi_label) <- levels(epi)
epi <- RenameIdents(epi, epi_label)
epi_label <- epi@active.ident
epi$epithelial_labels <- as.character(epi_label)

# Endothelial cells
ElbowPlot(smc, ndims = 30)

Ec <- FindNeighbors(Ec, dims = 1:18)
Ec <- FindClusters(Ec, resolution = 0.35)
Ec <- RunUMAP(Ec, dims = 1:18, reduction.name = "ec_umap18")


Idents(Ec) <- 'RNA_snn_res.0.3'
end_label <- c('C00','C01','C02','C03',
             'C04','C05')

DimPlot(Ec, reduction="ec_umap20") + ggtitle('Blood progenitors subcluster')

names(end_label) <- levels(Ec)
Ec <- RenameIdents(Ec, end_label)
end_label <- Ec@active.ident
Ec$endothelial_labels <- as.character(end_label)


# Blood progenitors

ElbowPlot(BloodP, ndims = 30)

DefaultAssay(BloodP) <- "integrated"
BloodP <- FindNeighbors(BloodP, dims = 1:18)
BloodP <- FindClusters(BloodP, resolution = 0.2)
BloodP <- RunUMAP(BloodP, dims = 1:18, reduction.name = "blood_res0.2")

DimPlot(BloodP, reduction="blood_res0.2", group.by = 'blood_labels') + ggtitle('Blood progenitors subcluster')


Idents(BloodP) <- 'RNA_snn_res.0.2'
blood_label <- c('C00-Mieloid_Progenitors','C01-Macrophages','C02-Microglia','C03-Monocytes')
names(blood_label) <- levels(BloodP)
BloodP <- RenameIdents(BloodP, blood_label)
blood_label <- BloodP@active.ident
BloodP$blood_labels <- as.character(blood_label)


# Cardiomyocytes
ElbowPlot(Cardio, ndims = 30)
DefaultAssay(Cardio) <- "integrated"
Cardio <- FindNeighbors(Cardio, dims = 1:10)
Cardio <- FindClusters(Cardio, resolution = 0.15)
Cardio <- RunUMAP(Cardio, dims = 1:10, reduction.name = "card_res0.1")

DimPlot(Cardio, reduction="card_res0.1", group.by = 'Cardio_labels', pt.size = 2)+ ggtitle('Cardiomyocytes subclustering')


Idents(BloodP) <- 'RNA_snn_res.0.15'
Cardio_label <- c('C00','C01','C02')
Cardio_cols <- c('C00'='red','C01'='green','C02'='blue') 
names(Cardio_label) <- levels(Cardio)
Cardio <- RenameIdents(Cardio, Cardio_label)
Cardio_label <- Cardio@active.ident
Cardio$Cardio_labels <- as.character(Cardio_label)

# Brain cluster 6 
ElbowPlot(fbrain, ndims = 30)
fbrain <- FindNeighbors(fbrain, dims = 1:18)
fbrain <- FindClusters(fbrain, resolution = 0.25)
fbrain <- RunUMAP(fbrain, dims = 1:18, reduction.name = "fbrain_res0.25")

p<- DimPlot(fbrain, reduction="fbrain_res0.25", pt.size = 2, group.by = 'fbrain_labels', cols = fbrain_cols)


Idents(fbrain) <- 'RNA_snn_res.0.25'
fbrain_label <- c('C00','C01','C02', 'C03', 'C04')
fbrain_cols <- c('C00'='#f1490e','C01'='#16923d','C02'='#e8ce28', 'C03'='#0e7cf1', 'C04'='#410c7c') 
names(fbrain_label) <- levels(fbrain)
fbrain <- RenameIdents(fbrain, fbrain_label)
fbrain_label <- fbrain@active.ident
fbrain$fbrain_labels <- as.character(fbrain_label)

saveRDS(fbrain, file='BrainC06Cluster.rds')

# SMC
ElbowPlot(smc, ndims =20)
smc <- FindNeighbors(smc, dims = 1:15)
smc <- FindClusters(smc, resolution = 0.3)
smc <- RunUMAP(smc, dims = 1:15, reduction.name = "smc_15")

DimPlot(smc, reduction="smc_15", group.by = 'smc_labels', pt.size = 2, cols = smc_cols)+  ggtitle('Subcluster of SMCs')


Idents(smc) <- 'RNA_snn_res.0.3'
smc_label <- c('C00','C01','C02','C03', 'C04')
smc_cols <- c('C00'='#9ef57b','C01'='#f69521','C02'='#1899e7', 'C03'='#f5e710', 'C04'='#d67bf5') 
names(smc_label) <- levels(smc)
smc <- RenameIdents(smc, smc_label)
smc_label <- smc@active.ident
smc$smc_labels <- as.character(smc_label)
saveRDS(smc, file='SmoothMuscleCells_v2.rds')

# C16_Osteoblas and fibroblas
ElbowPlot(osfi, ndims = 45)

osfi <- FindNeighbors(osfi, dims = 1:25)
osfi <- FindClusters(osfi, resolution = 0.18)
osfi <- RunUMAP(osfi, dims = 1:25, reduction.name = "osfi_umap0.25")


DimPlot(osfi, reduction="osfi_umap0.35", pt.size = 2, group.by = 'osfi_labels', cols=osfi_cols)

Idents(osfi) <- 'RNA_snn_res.0.25'
osfi_label <- c('C00','C01','C02','C03')
osfi_cols <- c('C00'='#f9d520','C01'='#26c556','C02'='#f93420','C03'='#a912f9')
names(osfi_label) <- levels(osfi)
osfi <- RenameIdents(osfi, osfi_label)
osfi_label <- osfi@active.ident
osfi$osfi_labels <- as.character(osfi_label)
saveRDS(osfi, file='C16_Osteoblast_fibroblasts_v2.rds')

# C15-Mesenchymal cells2
ElbowPlot(mesenc, ndims = 30)
mesenc <- FindNeighbors(mesenc, dims = 1:18)
mesenc <- FindClusters(mesenc, resolution = 0.3)
mesenc <- RunUMAP(mesenc, dims = 1:18, reduction.name = "mesenc_res0.3")

DimPlot(mesenc, reduction="mesenc_res0.3", group.by = 'mesenc_labels', cols=mesenc_cols)

Idents(mesenc) <- 'RNA_snn_res.0.3'
mesenc_label <- c('C00','C01','C02','C03')
mesenc_cols <- c('C00'='#f37c13','C01'='#26c556','C02'='#1361f3','C03'='#f92aaa')
names(mesenc_label) <- levels(mesenc)
mesenc <- RenameIdents(mesenc, mesenc_label)
mesenc_label <- mesenc@active.ident
mesenc$mesenc_labels <- as.character(mesenc_label)
saveRDS(mesenc, file='C15-Mesenchymal cells2_v2.rds')

# C14-Mesenchymal cells
ElbowPlot(mesenc1, ndims = 30)
mesenc1 <- FindNeighbors(mesenc1, dims = 1:18)
mesenc1 <- FindClusters(mesenc1, resolution = 0.2)
mesenc1 <- RunUMAP(mesenc1, dims = 1:18, reduction.name = "mesenc1_res0.2")

DimPlot(mesenc1, reduction="mesenc1_res0.2", group.by='mesenc1_labels', cols=mesenc1_cols)

Idents(mesenc1) <- 'RNA_snn_res.0.2'
mesenc1_label <- c('C00','C01','C02','C03', 'C04', 'C05')
mesenc1_cols <- c('C00'='#e8150f','C01'='#26c556','C02'= '#0ff6d3','C03'='#eff60f', 'C04'='#c94bfc', 'C05'= '#f6890f')
names(mesenc1_label) <- levels(mesenc1)
mesenc1 <- RenameIdents(mesenc1, mesenc1_label)
mesenc1_label <- mesenc1@active.ident
mesenc1$mesenc1_labels <- as.character(mesenc1_label)
saveRDS(mesenc1, file='C15-mesenc1hymal cells2_v2.rds')

# C13-Chondrocytes
ElbowPlot(chond, ndims = 30)
chond <- FindNeighbors(chond, dims = 1:18)
chond <- FindClusters(chond, resolution = 0.2)
chond <- RunUMAP(chond, dims = 1:18, reduction.name = "chond_res0.2")

DimPlot(chond, reduction="chond_res0.2", group.by = 'chond_labels', cols=chond_cols)

Idents(chond) <- 'RNA_snn_res.0.2'
chond_label <- c('C00','C01','C02','C03')
chond_cols <- c('C00'='#e8150f','C01'='#26c556','C02'= '#f3db0e','C03'='#9312e8')
names(chond_label) <- levels(chond)
chond <- RenameIdents(chond, chond_label)
chond_label <- chond@active.ident
chond$chond_labels <- as.character(chond_label)
saveRDS(chond, file='C13_chondrocytes_v2.rds')

# C11-myocytes
ElbowPlot(myoc, ndims = 30)
myoc <- FindNeighbors(myoc, dims = 1:18)
myoc <- FindClusters(myoc, resolution = 0.2)
myoc <- RunUMAP(myoc, dims = 1:18, reduction.name = "myoc_res0.2")

DimPlot(myoc, reduction="myoc_res0.2")
DimPlot(myoc, reduction="myoc_res0.2", group.by = 'myoc_labels', cols=myoc_cols)

Idents(myoc) <- 'RNA_snn_res.0.2'
myoc_label <- c('C00','C01','C02','C03')
myoc_cols <- c('C00'='#f7a5fa','C01'='#26c556','C02'= '#f97c06','C03'='#9312e8')
names(myoc_label) <- levels(myoc)
myoc <- RenameIdents(myoc, myoc_label)
myoc_label <- myoc@active.ident
myoc$myoc_labels <- as.character(myoc_label)
saveRDS(myoc, file='C11_myocytes_v2.rds')

# C05 Glial cells
ElbowPlot(glial, ndims = 30)
glial <- FindNeighbors(glial, dims = 1:18)
glial <- FindClusters(glial, resolution = 0.2)
glial <- RunUMAP(glial, dims = 1:18, reduction.name = "glial_res0.2")

DimPlot(glial, reduction="glial_res0.2")
DimPlot(glial, reduction="glial_res0.2", group.by = 'glial_labels', cols=glial_cols)

Idents(glial) <- 'RNA_snn_res.0.2'
glial_label <- c('C00','C01','C02','C03')
glial_cols <- c('C00'='#f7a5fa','C01'='#1206f9','C02'= '#f90678','C03'='#30c507')
names(glial_label) <- levels(glial)
glial <- RenameIdents(glial, glial_label)
glial_label <- glial@active.ident
glial$glial_labels <- as.character(glial_label)
saveRDS(glial, file='C05_GlialCells_v2.rds')

# C00 neural crest 
ElbowPlot(neurl, ndims = 40)
neurl <- FindNeighbors(neurl, dims = 1:20)
neurl <- FindClusters(neurl, resolution = 0.25)
neurl <- RunUMAP(neurl, dims = 1:20, reduction.name = "neurl_res0.25")

DimPlot(neurl, reduction="neurl_res0.25")
DimPlot(neurl, reduction="neurl_res0.25", group.by = 'neurl_labels', cols=neurl_cols)

Idents(neurl) <- 'RNA_snn_res.0.2'
neurl_label <- c('C00','C01','C02','C03')
neurl_cols <- c('C00'='#27a2d4','C01'='#0925c7','C02'= '#fb901f')
names(neurl_label) <- levels(neurl)
neurl <- RenameIdents(neurl, neurl_label)
neurl_label <- neurl@active.ident
neurl$neurl_labels <- as.character(neurl_label)
saveRDS(neurl, file='C00_NeuralCrest_v2.rds')




####Subset control cells and find clusters and markers ##############

# Endothelial cells
Idents(Ec) <- 'Condition'
Ec_ctrl <- subset(Ec, idents="RbpjWt")


ec_mark_res0.3 <- FindAllMarkers(Ec,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
Idents(Ec_ctrl) <- 'endothelial_labels3'
ec_mark_Ctrl <- FindAllMarkers(Ec_ctrl,logfc.threshold = 0.5, only.pos = T, verbose = T, min.diff.pct = -Inf )


Idents(Ec) <- 'endothelial_labels3'
ec_mark <- FindAllMarkers(Ec,logfc.threshold = 0.5, only.pos = T, verbose = T, min.diff.pct = -Inf )

write.csv(ec_mark, file='DEGs_ECs_PerCluster.csv')

# Epithelial cells
Idents(epi) <- 'Condition'
epi_down <- subset(epi, downsample=364)

epi_ctrl <- subset(epi, idents="RbpjWt")
epi_markers_rbpj <- FindAllMarkers(epi,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(epi_markers_rbpj, file='Markers_DEGS_Rbpj.csv')

Idents(epi) <- 'epithelial_labels'
epi_markers <- FindAllMarkers(epi,logfc.threshold = 0.5, only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(epi_markers, file='Markers_DEGS_Clusters_Epithelial.csv')


epi_C0_Rbpj <- FindMarkers(epi,ident.1 ="RbpjWt_C00", ident.2="RbpjLOF_C00")
epi_C1_Rbpj <- FindMarkers(epi,ident.1 ="RbpjWt_C01", ident.2="RbpjLOF_C01") #fold change indicates more expressed in ident.2
epi_C2_Rbpj <- FindMarkers(epi,ident.1 ="RbpjWt_C02", ident.2="RbpjLOF_C02")
epi_C3_Rbpj <- FindMarkers(epi,ident.1 ="RbpjWt_C03", ident.2="RbpjLOF_C03")


#Blood progenitors
Idents(BloodP) <- 'Condition'
BloodP_ctrl <- subset(BloodP, idents="RbpjWt")

Idents(BloodP) <- 'blood_labels'
levels(BloodP) <- c( "C00-Mieloid_Progenitors", "C01-Macrophages"  ,"C02-Microglia", "C03-Monocytes"   )
BloodP_mark <- FindAllMarkers(BloodP,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(BloodP_mark, file='DEGs_BloodPsAll.csv')

Idents(BloodP_ctrl) <- 'blood_labels'
Idents(BloodP_ctrl) <- factor(x = Idents(BloodP_ctrl), levels = sort(levels(BloodP_ctrl))) # to order levels
BloodP_mark_Ctrl <- FindAllMarkers(BloodP_ctrl,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )

write.csv(BloodP_mark_Ctrl, file='DEGs_BloodPs_Wtbased.csv')
Idents(BloodP) <- 'Condition'
blood_markers <- FindMarkers(BloodP, ident.1 = 'RbpjLOF')
write.csv(blood_markers, file='DEGs_Blood_Rbpj.csv')

Idents(BloodP) <- 'Condition'
blood_downs <- subset(BloodP, downsample=228)


Idents(BloodP) <- 'Condition_clustering'

blood_C0_Rbpj <- FindMarkers(BloodP,ident.1 ="RbpjWt_C00-Mieloid_Progenitors", ident.2="RbpjLOF_C00-Mieloid_Progenitors")
blood_C1_Rbpj <- FindMarkers(BloodP,ident.1 ="RbpjWt_C01-Macrophages", ident.2="RbpjLOF_C01-Macrophages") #fold change indicates more expressed in ident.2
blood_C2_Rbpj <- FindMarkers(BloodP,ident.1 ="RbpjWt_C02-Microglia", ident.2="RbpjLOF_C02-Microglia")
blood_C3_Rbpj <- FindMarkers(BloodP,ident.1 ="RbpjWt_C03-Monocytes", ident.2="RbpjLOF_C03-Monocytes")


#Cardiomyocytes progenitors
Idents(Cardio) <- 'Condition'
Cardio_ctrl <- subset(Cardio, idents="RbpjWt")

Idents(Cardio_ctrl) <- 'Cardio_labels'
Idents(Cardio_ctrl) <- factor(x = Idents(Cardio_ctrl), levels = sort(levels(Cardio_ctrl))) # to order levels
Cardio_mark_Ctrl <- FindAllMarkers(Cardio_ctrl,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )

write.csv(Cardio_mark_Ctrl, file='DEGs_Cardios_Wtbased.csv')

# Brain central nervous system cells
Idents(ob) <- 'new_labels'
fbrain <- subset(ob, idents = "C06-Fetal brain (CNS)")
Idents(fbrain) <- 'Condition'
brain_markers <- FindMarkers(fbrain, ident.1 = 'RbpjLOF')
write.csv(brain_markers, file='DEGs_Brain_Rbpj.csv')

Idents(brain) <- 'Condition_Clusterting' #to see DEGs Rbpj/WT within clusters

brain_C0_Rbpj <- FindMarkers(brain,ident.1 ="RbpjLOF_C00", ident.2="RbpjWt_C00")
brain_C1_Rbpj <- FindMarkers(brain,ident.1 ="RbpjWt_C01", ident.2="RbpjLOF_C01") #fold change indicates more expressed in ident.2
brain_C2_Rbpj <- FindMarkers(brain,ident.1 ="RbpjWt_C02", ident.2="RbpjLOF_C02")
brain_C3_Rbpj <- FindMarkers(brain,ident.1 ="RbpjWt_C03", ident.2="RbpjLOF_C03")
brain_C4_Rbpj <- FindMarkers(brain,ident.1 ="RbpjWt_C04", ident.2="RbpjLOF_C04")


write.csv(brain_C4_Rbpj, file='DEGs_Brain_C4.csv')


# SMC smooth muscle cells
Idents(smc) <- 'Condition_Clusterting2' #to see DEGs Rbpj/WT within clusters

smc_C0_Rbpj <- FindMarkers(smc,ident.1 ="RbpjLOF_C00", ident.2="RbpjWt_C00")
smc_C1_Rbpj <- FindMarkers(smc,ident.1 ="RbpjWt_C01", ident.2="RbpjLOF_C01") #fold change indicates more expressed in ident.2
smc_C2_Rbpj <- FindMarkers(smc,ident.1 ="RbpjWt_C02", ident.2="RbpjLOF_C02")
smc_C3_Rbpj <- FindMarkers(smc,ident.1 ="RbpjWt_C03", ident.2="RbpjLOF_C03")
smc_C4_Rbpj <- FindMarkers(smc,ident.1 ="RbpjWt_C04", ident.2="RbpjLOF_C04")


# C16_Osteoblast_fibroblasts
Idents(osfi) <- 'Condition'
osfi_DEGs_Rbpj <- FindAllMarkers(osfi,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(osfi_DEGs_Rbpj, file='osfi_DEGs_Rbpj.csv')

Idents(osfi) <- 'osfi_labels'
osfi_DEGs_clusters <- FindAllMarkers(osfi,logfc.threshold = 0.5, only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(osfi_DEGs_clusters, file='osfi_DEGs_Cluster.csv')



# C15-Mesenchymal cells2
Idents(mesenc) <- 'Condition'
mesenc_DEGs_Rbpj <- FindAllMarkers(mesenc,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(mesenc_DEGs_Rbpj, file='mesenc_DEGs_Rbpj.csv')

Idents(mesenc) <- 'mesenc_labels'
mesenc_DEGs_clusters <- FindAllMarkers(mesenc,logfc.threshold = 0.5, only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(mesenc_DEGs_clusters, file='mesenc_DEGs_Cluster.csv')

# C14-Mesenchymal cells1
Idents(mesenc1) <- 'Condition'
mesenc1_DEGs_Rbpj <- FindAllMarkers(mesenc1,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(mesenc1_DEGs_Rbpj, file='mesenc1_DEGs_Rbpj.csv')

Idents(mesenc1) <- 'mesenc1_labels'
levels(mesenc1) <- c('C00','C01','C02','C03', 'C04', 'C05')
mesenc1_DEGs_clusters <- FindAllMarkers(mesenc1,logfc.threshold = 0.5, only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(mesenc1_DEGs_clusters, file='mesenc1_DEGs_Cluster.csv')


# C13-Chondrocytes
Idents(chond) <- 'Condition'
chond_DEGs_Rbpj <- FindAllMarkers(chond,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(chond_DEGs_Rbpj, file='chond_DEGs_Rbpj.csv')

Idents(chond) <- 'chond_labels'
levels(chond) <- c('C00','C01','C02','C03')
chond_DEGs_clusters <- FindAllMarkers(chond,logfc.threshold = 0.5, only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(chond_DEGs_clusters, file='chond_DEGs_Cluster.csv')

# C11_Myocytes
Idents(myoc) <- 'Condition'
myoc_DEGs_Rbpj <- FindAllMarkers(myoc,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(myoc_DEGs_Rbpj, file='myoc_DEGs_Rbpj.csv')

Idents(myoc) <- 'myoc_labels'
levels(myoc) <- c('C00','C01','C02','C03')
myoc_DEGs_clusters <- FindAllMarkers(myoc,logfc.threshold = 0.5, only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(myoc_DEGs_clusters, file='myoc_DEGs_Cluster.csv')

# C05-Glial cells
Idents(glial) <- 'Condition'
glial_DEGs_Rbpj <- FindAllMarkers(glial,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(glial_DEGs_Rbpj, file='glial_DEGs_Rbpj.csv')

Idents(glial) <- 'glial_labels'
levels(glial) <- c('C00','C01','C02','C03')
glial_DEGs_clusters <- FindAllMarkers(glial,logfc.threshold = 0.5, only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(glial_DEGs_clusters, file='glial_DEGs_Cluster.csv')


# C00-Neural crest cells
Idents(neurl) <- 'Condition'
neurl_DEGs_Rbpj <- FindAllMarkers(neurl,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(neurl_DEGs_Rbpj, file='neurl_DEGs_Rbpj.csv')

Idents(neurl) <- 'neurl_labels'
levels(neurl) <- c('C00','C01','C02')
neurl_DEGs_clusters <- FindAllMarkers(neurl,logfc.threshold = 0.5, only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(neurl_DEGs_clusters, file='neurl_DEGs_Cluster.csv')




# Find markers per cluster and per genotype/cell-type ####
Idents(Ec) <- 'Condition'
ec_rbpj_mark <- FindAllMarkers(Ec,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
ec_markers <- FindMarkers(Ec, ident.1 = 'RbpjLOF', logfc.threshold = 0.5)
ec_markers_t <- FindMarkers(Ec, ident.1 = 'RbpjWt', test.use = 't')
write.csv(ec_markers, file='DEGs_ECs_Rbpj.csv')
write.csv(ec_rbpj_mark, file='DEGs_ECs_Rbpj.csv')




Idents(brain) <- 'fbrain_labels'
brain_markers_cluster <- FindAllMarkers(brain,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
write.csv(brain_markers_cluster, file='DEGs_Brain_Percluster.csv')

# Find markers per cluster and per genotype/cell-type ####
Idents(smc) <- 'Condition'
smc_mark_rbpj <- FindAllMarkers(smc,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )
Idents(smc) <- 'smc_labels'
levels(smc) <- c('C00', 'C01', 'C02', 'C03', 'C04')
smc_mark_clust <- FindAllMarkers(smc,logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf )


write.csv(smc_mark_clust, file='DEGs_smc_Percluster.csv')
write.csv(smc_mark_rbpj, file='DEGs_smc_Rbpj.csv')

#######Select 5 top genes excluding duplicates.#####

#Endothelial (top n)
ec_mark_Ctrl$order <- 1:nrow(ec_mark_Ctrl)
Top_ec <- ec_mark_Ctrl %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_ec$gene)
Top_ec <- Top_ec[match(s, Top_ec$gene),]
Top10_ec <- Top_ec%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 10, wt = -order)
Top10genes_ec <- Top10_ec$gene



ec_mark$order <- 1:nrow(ec_mark)
Top_ec <- ec_mark %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 10, wt = -order)
s<-unique(Top_ec$gene)
Top_ec <- Top_ec[match(s, Top_ec$gene),]
Top5_ec <- Top_ec%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_ec <- Top5_ec$gene



#Brain (top n)
brain_markers_cluster$order <- 1:nrow(brain_markers_cluster)
Top_brain <- brain_markers_cluster %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_brain$gene)
Top_brain <- Top_brain[match(s, Top_brain$gene),]
Top10_brain <- Top_brain%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 10, wt = -order)
Top10genes_brain <- Top10_brain$gene

brain_markers_cluster$order <- 1:nrow(brain_markers_cluster)
Top_brain <- brain_markers_cluster %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 15, wt = -order)
s<-unique(Top_brain$gene)
Top_brain <- Top_brain[match(s, Top_brain$gene),]
Top5_brain <- Top_brain%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_brain <- Top5_brain$gene



mark <- brain_C4_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
head(new_mark %>% arrange(desc(avg_log2FC)))
tail(new_mark %>% arrange(desc(avg_log2FC)))

DEGs_C_Rbpj_more <- c('Atp5k')
DEGs_C_Rbpj_less <- c('Ncam1','Mrps23', 'Bcl11b')

#Epithelial (top n)
epi_markers$order <- 1:nrow(epi_markers)
Top_epi <- epi_markers %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_epi$gene)
Top_epi <- Top_epi[match(s, Top_epi$gene),]
Top10_epi <- Top_epi%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 10, wt = -order)
Top10genes_epi <- Top10_epi$gene

epi_markers$order <- 1:nrow(epi_markers)
Top_epi <- epi_markers %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 15, wt = -order)
s<-unique(Top_epi$gene)
Top_epi <- Top_epi[match(s, Top_epi$gene),]
Top5_epi <- Top_epi%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_epi <- Top5_epi$gene



mark <- epi_markesr_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
head(new_mark %>% arrange(desc(avg_log2FC)))
Top5_Rbpj_more <- c('Fabp7','Cxcl14', 'Col1a1', 'Vim', 'Krt8', 'H1f5')


mark <- epi_C4_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
head(new_mark %>% arrange(desc(avg_log2FC)))
tail(new_mark %>% arrange(desc(avg_log2FC)))

DEGs_C_Rbpj_more <- c('H1f5')
DEGs_C_Rbpj_less <- c()

#Blood progenitors
BloodP_mark_Ctrl$order <- 1:nrow(BloodP_mark_Ctrl)
Top_BloodP <- BloodP_mark_Ctrl %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 10, wt = -order)
s<-unique(Top_BloodP$gene)
Top_BloodP <- Top_BloodP[match(s, Top_BloodP$gene),]
Top5_BloodP <- Top_BloodP%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_BloodP <- Top5_BloodP$gene

BloodP_mark$order <- 1:nrow(BloodP_mark)
Top_BloodP <- BloodP_mark %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_BloodP$gene)
Top_BloodP <- Top_BloodP[match(s, Top_BloodP$gene),]
Top10_BloodP <- Top_BloodP%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 10, wt = -order)
Top10genes_BloodP <- Top10_BloodP$gene



mark <- blood_C3_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
head(new_mark %>% arrange(desc(avg_log2FC)))
tail(new_mark %>% arrange(desc(avg_log2FC)))

DEGs_C_Rbpj_more <- c('Col1a2','Ptn', 'Vim')
DEGs_C_Rbpj_less <- c('Ifitm1', 'Cx3cr1')


#Cardio progenitors
Cardio_mark_Ctrl$order <- 1:nrow(Cardio_mark_Ctrl)
Top_Cardio <- Cardio_mark_Ctrl %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 10, wt = -order)
s<-unique(Top_Cardio$gene)
Top_Cardio <- Top_Cardio[match(s, Top_Cardio$gene),]
Top5_Cardio <- Top_Cardio%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_Cardio <- Top5_Cardio$gene

#SMC progenitors
smc_mark_clust$order <- 1:nrow(smc_mark_clust)
Top_smc <- smc_mark_clust %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 30, wt = -order)
s<-unique(Top_smc$gene)
Top_smc <- Top_smc[match(s, Top_smc$gene),]
Top5_smc <- Top_smc%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_smc <- Top5_smc$gene

mark <- smc_C4_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
head(new_mark %>% arrange(desc(avg_log2FC)))
tail(new_mark %>% arrange(desc(avg_log2FC)))

DEGs_C_Rbpj_more <- c('Actg2', 'Fabp7', 'Acta2')
DEGs_C_Rbpj_less <- c()


#C16_Osteoblast_fibroblasts (top n)

osfi_DEGs_clusters$order <- 1:nrow(osfi_DEGs_clusters)
Top_osfi <- osfi_DEGs_clusters %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_osfi$gene)
Top_osfi <- Top_osfi[match(s, Top_osfi$gene),]
Top5_osfi <- Top_osfi%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_osfi <- Top5_osfi$gene

mark <- osfi_DEGs_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
head(new_mark %>% arrange(desc(avg_log2FC)))
tail(new_mark %>% arrange(desc(avg_log2FC)))

#C15-Mesenchymal cells2 (top n)

mesenc_DEGs_clusters$order <- 1:nrow(mesenc_DEGs_clusters)
Top_mesenc <- mesenc_DEGs_clusters %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_mesenc$gene)
Top_mesenc <- Top_mesenc[match(s, Top_mesenc$gene),]
Top5_mesenc <- Top_mesenc%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_mesenc <- Top5_mesenc$gene

mark <- mesenc_DEGs_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
head(new_mark %>% arrange(desc(avg_log2FC)))
tail(new_mark %>% arrange(desc(avg_log2FC)))

#C14-Mesenchymal cells1 (top n)

mesenc1_DEGs_clusters$order <- 1:nrow(mesenc1_DEGs_clusters)
Top_mesenc1 <- mesenc1_DEGs_clusters %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_mesenc1$gene)
Top_mesenc1 <- Top_mesenc1[match(s, Top_mesenc1$gene),]
Top5_mesenc1 <- Top_mesenc1%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_mesenc1 <- Top5_mesenc1$gene

mark <- mesenc1_DEGs_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
t<-rownames(head(new_mark %>% arrange(desc(avg_log2FC))))
tail(new_mark %>% arrange(desc(avg_log2FC)))

#C13-Chondrocytes (top n)

chond_DEGs_clusters$order <- 1:nrow(chond_DEGs_clusters)
Top_chond <- chond_DEGs_clusters %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_chond$gene)
Top_chond <- Top_chond[match(s, Top_chond$gene),]
Top5_chond <- Top_chond%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_chond <- Top5_chond$gene

mark <- chond_DEGs_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
t<-rownames(head(new_mark %>% arrange(desc(avg_log2FC))))
tail(new_mark %>% arrange(desc(avg_log2FC)))


#C11-Myocytes (top n)
myoc_DEGs_clusters$order <- 1:nrow(myoc_DEGs_clusters)
Top_myoc <- myoc_DEGs_clusters %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_myoc$gene)
Top_myoc <- Top_myoc[match(s, Top_myoc$gene),]
Top5_myoc <- Top_myoc%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_myoc <- Top5_myoc$gene

mark <- myoc_DEGs_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
t<-rownames(head(new_mark %>% arrange(desc(avg_log2FC))))
tail(new_mark %>% arrange(desc(avg_log2FC)))

#C05-Glial cells (top n)
glial_DEGs_clusters$order <- 1:nrow(glial_DEGs_clusters)
Top_glial <- glial_DEGs_clusters %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_glial$gene)
Top_glial <- Top_glial[match(s, Top_glial$gene),]
Top5_glial <- Top_glial%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_glial <- Top5_glial$gene

mark <- glial_DEGs_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
t<-rownames(head(new_mark %>% arrange(desc(avg_log2FC))))
tail(new_mark %>% arrange(desc(avg_log2FC)))

#C00 Neural crest (top n)
neurl_DEGs_clusters$order <- 1:nrow(neurl_DEGs_clusters)
Top_neurl <- neurl_DEGs_clusters %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 20, wt = -order)
s<-unique(Top_neurl$gene)
Top_neurl <- Top_neurl[match(s, Top_neurl$gene),]
Top5_neurl <- Top_neurl%>% arrange(cluster) %>% group_by(cluster) %>% top_n(n = 5, wt = -order)
Top5genes_neurl <- Top5_neurl$gene

mark <- neurl_DEGs_Rbpj
new_mark <- subset(mark, p_val_adj<0.05)
t<-rownames(head(new_mark %>% arrange(desc(avg_log2FC))))
tail(new_mark %>% arrange(desc(avg_log2FC)))



### labels and colours of subclusters ####


#Endothelial
end_cols <- c('C00'='#befb0e','C01'='#c778f1',
              'C02-endoEMT'='#ffa600','C03-BBB_ECs'='#479608','C04-Liver_ECs'='#3175de')

end_cols3 <- c('C00-mesenchymal-proliferative'='#cb4cf3', 'C01-endoEMT'='#ffa600',
               'C02-BBB_ECs'='#479608','C03-Liver_ECs'='#3175de')
t1 <- c('#befb0e','#c778f1','#ffa600','#479608','#3175de')

end_label <- c('C00','C01','C02','C03',
               'C04','C05')
#BloodProgenitors
blood_cols <-  c('C00-Mieloid_Progenitors'='#ffa600','C01-Macrophages'='#e718b5','C02-Microglia'='#18e799','C03-Monocytes'='#3d18e7')
t2 <- c('#ffa600','#e718b5','#18e799','#3d18e7')



### Rename clusters ######


Idents(Ec)<-'endothelial_labels'
new_names <- c('C00-mesenchymal-proliferative', 'C00-mesenchymal-proliferative','C01-endoEMT',
               'C02-BBB_ECs','C03-Liver_ECs')

names(new_names) <- levels(Ec)
Ec <- RenameIdents(Ec, new_names)
new_names <- Ec@active.ident
Ec$endothelial_labels3 <- as.character(new_names)
saveRDS(Ec, file='EndothlelialCellsv3.rds')

#Downsample smc to have same amount of control and mutant cells
Idents(smc) <- 'Condition'
smc_down <- subset(smc, downsample=297)

#Downsample fibroblas to have same amount of control and mutant cells
Idents(osfi) <- 'Condition'
osfi_down <- subset(osfi, downsample=1161 )

#Downsample EC to have same amount of control and mutant cells
Idents(Ec) <- 'Condition'
Ec_same <- subset(Ec, downsample=350)

#Downsample Mesenchymal to have same amount of control and mutant cells
Idents(mesenc) <- 'Condition'
mesenc_down<- subset(mesenc, downsample= 1406)
#Downsample Mesenchymal 2to have same amount of control and mutant cells
Idents(mesenc1) <- 'Condition'
mesenc1_down <- subset(mesenc1, downsample= 1474)

#Downsample chondrocytes to have same amount of control and mutant cells
Idents(chond) <- 'Condition'
chond_down <- subset(chond, downsample= 739 )

#Downsample myocytes to have same amount of control and mutant cells
Idents(myoc) <- 'Condition'
myoc_down <- subset(myoc, downsample= 115 )

#Downsample glial cells to have same amount of control and mutant cells
Idents(glial) <- 'Condition'
glial_down <- subset(glial, downsample= 233)

#Downsample neural crest to have same amount of control and mutant cells
Idents(neurl) <- 'Condition'
neurl_down <- subset(neurl, downsample= 158)


DimPlot(Ec, reduction='ec_umap20', group.by = 'endothelial_labels3',
        pt.size = 2, split.by = 'Condition', cols = end_cols3)



dittoBarPlot(Ec, var='endothelial_labels2', group.by = 'Condition', color.panel = end_cols, x.reorder = c(2,1))+
  ggtitle('Endothelial subclusters (contribution by genotype)')+
  theme(axis.text.x = element_text(angle=0, size=15, hjust = 0.5))

dittoBarPlot(Ec, var='Condition', group.by = 'endothelial_labels2', color.panel =c('#f83d14','#80e13d') )+
  ggtitle('Endothelial subclusters (contribution by genotype)')+
  theme(axis.text.x = element_text(angle=0, size=10, hjust = 0.5))


dittoBarPlot(Cardio, var='Condition', group.by = 'Cardio_labels', color.panel =c('#f83d14','#80e13d') )+
  ggtitle('Cardiomyocytes subclusters (contribution by genotype)')+
  theme(axis.text.x = element_text(angle=0, size=10, hjust = 0.5))

####### Plots for ECs order levels ##########

Idents(Ec_ctrl) <- 'endothelial_labels2'
Idents(Ec_ctrl) <- factor(x = Idents(Ec_ctrl), levels = sort(levels(Ec_ctrl))) # to order levels

Idents(BloodP_ctrl) <- 'blood_labels'
Idents(BloodP_ctrl) <- factor(x = Idents(BloodP_ctrl), levels = sort(levels(BloodP_ctrl))) # to order levels
##__________ Dot plot ####
df <- data.frame(
  x = rep(c(3,8,13,18,23)),
  y = rep(c(0),1),
  z = factor(rep(1:5)),
  color = t1)


DotPlot(Ec_ctrl, features = Top5genes_ec, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label = Top5genes_ec)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = t1)


#Blood cells
df <- data.frame(
  x = rep(c(3,8,13,18)),
  y = rep(c(0),1),
  z = factor(rep(1:4)),
  color = t2)


DotPlot(BloodP_ctrl, features = Top5genes_BloodP, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label = Top5genes_BloodP)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = t2)






# Other graphs ####

dittoBarPlot(Ec, var='Condition', group.by = 'endothelial_labels2',
             color.panel = c('#f83d14','#80e13d')) + ggtitle('KO/WT distribution by endothelial cluster')

dittoBarPlot(Ec, var='endothelial_labels2',group.by = 'Condition',
             color.panel = end_cols, x.reorder = c(2,1))+ 
            ggtitle('Cluster proportion by cell type')+
              theme(axis.title.y.left = element_text(size=15),
                    axis.title.x.bottom = element_blank(),
                    title = element_text(size = 17),
                    axis.text.x = element_text(angle=0, size=20, hjust=0.5))

VlnPlot(Ec, feature='Oit3',group.by = 'endothelial_labels2',
        split.by = 'Condition', cols = c('#80e13d','#f83d14'))

notch <- c('Hes1', 'Hes5', 'Hey1', 'Hey3','HeyL' )
lig <- c('Dll4', 'Dll1', 'Jag1')
liv <- c('Oit3', 'Gata4','Stab2' )
muc <- c('Fabp4', 'Actb', 'Aqp1')

VlnPlot(Ec_ctrl, feature='Esm1', group.by = 'endothelial_labels2', cols = end_cols)+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90, hjust=1,vjust=0, size=13))

FeaturePlot(Ec_ctrl, features='Dlk1', reduction='ec_umap25')


DimPlot(Ec, reduction='ec_umap25', group.by = 'endothelial_labels2' , pt.size = 1.1, cols = end_cols)

DimPlot(BloodP, reduction='blood_res0.2', group.by = 'blood_labels', pt.size = 1.5, cols = blood_cols)+ ggtitle('Sub-cluster of blood progenitor cells')


dittoBarPlot(BloodP, var='blood_labels', group.by = 'Condition', color.panel = blood_cols, x.reorder = c(2,1))+
  ggtitle('Sub-cluster of blood progenitor cells')+
  theme(axis.title.y.left = element_text(size=15),
        axis.title.x.bottom = element_blank(),
        title = element_text(size = 17),
        axis.text.x = element_text(angle=0, size=20, hjust=0.5))

dittoBarPlot(BloodP, var='Condition', group.by = 'blood_labels', color.panel = c('#f83d14', '#80e13d'))+
  ggtitle('Sub-cluster of blood progenitor cells')+
    theme(axis.title.y.left = element_text(size=15),
        axis.title.x.bottom = element_blank(),
        title = element_text(size = 25),
        axis.text.x = element_text(angle=0, size=10, hjust=0.5))



FeaturePlot(Ec, features = 'Gata4', pt.size = 1, reduction = 'ec_umap25')+
  scale_colour_gradientn(colours = Bestholtz_palette)

FeaturePlot(BloodP, features = 'Itgam', pt.size = 1, reduction = 'blood_res0.2')+
  scale_colour_gradientn(colours = Bestholtz_palette)

VlnPlot(Ec, feature='Bmx', group.by = 'endothelial_labels2', split.by = 'Condition', cols = c( '#80e13d','#f83d14'))+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90, hjust=1,vjust=0, size=13))
  ggtitle('CD14 (Monocyte marker)')
  
  
VlnPlot(ob, feature=notch, group.by = 'new_labels', cols = c( '#80e13d','#f83d14'), split.by = 'Condition')
VlnPlot(Ec, feature=notch, group.by = 'endothelial_labels', cols = c( '#80e13d','#f83d14'), split.by = 'Condition')

#Get percentage of cells expressing a gene of interest ####
  
  test_genes <- c("Gja4", "Gja5", "Bmx", 'Efnb2')
  Idents(Ec) <- 'endothelial_labels2'
  a <- DotPlot(Ec, features = test_genes)
  a$data
  a<- as.data.frame(a$data)
  
  
ggplot(a, aes(x=features.plot, y=pct.exp, fill=id))+ geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values= end_cols)+
  ggtitle('Absolute percentage of cells with the expression of the gene per cluster in cells')+
  theme(axis.title.x = element_blank(), panel.background = element_blank())+
  ylab('Percentage of cells gene+ (Within each cluster)')+ labs(fill='Cell type')

Idents(Ec) <- 'Condition'
a <- DotPlot(Ec, features = test_genes)
a$data
a<- as.data.frame(a$data)

ggplot(a, aes(x=features.plot, y=pct.exp, fill=id))+ geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values= c('#80e13d','#f83d14'))+
  ggtitle('Absolute percentage of cells with the expression of the gene per cells type')+
  theme(axis.title.x = element_blank(), panel.background = element_blank())+
  ylab('Percentage of cells gene+ (separated by cell type)')+ labs(fill='Cell type')

ggplot(a, aes(x=features.plot, y=pct.exp))+ geom_bar(stat='identity') +
  scale_fill_manual(values= c('#80e13d','#f83d14'))+
  ggtitle('Absolute percentage of cells with the expression of the gene')+
  theme(axis.title.x = element_blank(), panel.background = element_blank())+
  ylab('Percentage of cells gene+ (separated by cell type)')+ labs(fill='Cell type')
  
  



# Get the expression modules for gene expression to identify cluster ID ####

arterial_genes <- list('Gja4', 'Bmx', 'Gja5', 'Cxcr4')
arterEC <- list('Gja4', 'Sat1', 'Hey1', 'Serpinf1', 'Igfbp3', 'Vegfc')
muscleEC <- list('Actb', 'Fabp4', 'Aqp1')
brainEC <- list('Slc2a1', 'Vwa1', 'Pglyrp1', 'Foxq1', 'Pltp', 'lsr', 'Foxf2', 'Mfsd2a')
kidneyEc <- list('Rsad2', 'Pbx1', 'Chchd10', 'Egln3')
endoEMT <- list('Dlk1', 'Meg3', 'Ptn', 'Col3a1', 'Ndufa4l2', 'Col1a2')
endCM <- list('Pdlim3', 'Gsta4', 'Smtnl2', 'Nkain4', 'Fn1', 'Igf2r')
vein <- list('Lbp', 'Spint2', 'Entpd1', 'Cthrc1', 'Bmp4')
liver <- list('Ehd3', 'Oit3', 'Gpr182', 'Tspan7', 'Bmp2', 'Stab2', 'Mt1', 'Mrc1', 'Maf', 'Lyve1', 'Akr1b8', 'Npl', 'Il1a')
lung <- list('Foxf1', 'Scn7a', 'Hapln1','Hoxa5', 'Car14', 'Phgdh')





Ec <- AddModuleScore(Ec, features=list(arterial_genes), name = "Art_g")
Ec <- AddModuleScore(Ec, features=list(muscleEC), name = "MusEcg")
Ec <- AddModuleScore(Ec, features=list(brainEC), name = "BrainEcg")
Ec <- AddModuleScore(Ec, features=list(arterEC), name = "ArtEcg")
Ec <- AddModuleScore(Ec, features=list(kidneyEc), name = "KidnEcg")
Ec <- AddModuleScore(Ec, features=list(endoEMT), name = "endoEMTg")
Ec <- AddModuleScore(Ec, features=list(endCM), name = "endCMg")
Ec <- AddModuleScore(Ec, features=list(vein), name = "veing")
Ec <- AddModuleScore(Ec, features=list(liver), name = "liverg")
Ec <- AddModuleScore(Ec, features=list(lung), name = "lungg")


# For blood analysis


macrophage_markers <- list('Itgam', 'Adgre1', 'Cd80','Cd86','Cd163', 'Mrc1')
BloodP <- AddModuleScore(BloodP, features=list(macrophage_markers), name = "Macroph")



# For epithelial analysis


Mesenchymal_markers <- list('Cdh2', 'Cdh1', 'Snail1', 'Snail2', 'Vim', 'Fn1')
EMT_markers <- list('Ctnnb1', 'Cdh11' , 'Gsk3b', 'Mmp2', 'Sumo2', 'Twist1', 'Zeb1')
epithelial_markers <- list('Col4a1', 'Dsg3', 'Lama2', 'Muc1', 'Sdc1')
hallmarks_EMT <- list(S)


epi <- AddModuleScore(epi, features=list(EMT_markers), name = "EMT")
epi <- AddModuleScore(epi, features=list(Mesenchymal_markers), name = "Mesenchymal")
epi <- AddModuleScore(epi, features=list(epithelial_markers), name = "Epithe")
epi <- AddModuleScore(epi, features=list(hallmarks_EMT), name = "Hallmark_EMT")

# For smc analysis


heart_aorta <- list('Ramp1', 'Synpo2', 'Pde3a', 'Ppp1r14a', 'Mrvi1', 'Rrad', 'Rasl12', 'Slmap', 'Mustn1', 'Filip1l')


smc <- AddModuleScore(smc, features=list(heart_aorta), name = "HA")

# Set variable condition + cluster for dotplot analysis ####




ncells = ncol(ob)

for (i in 1:ncells){
  b = as.character(ob@meta.data$Condition[i])
  a = as.character(ob@meta.data$new_labels_rb[i])
  
  ob@meta.data$Condition_Clustering[i] = paste0(b,"_",a)
  
}

Idents(ob) <- 'Condition_Clustering'

levels(ob) <- c(  "RbpjWt_C0-Mesenchymal cells2" ,"RbpjLOF_C0-Mesenchymal cells2" ,
                  "RbpjWt_C01-Osteoblast and fibroblasts" ,"RbpjLOF_C01-Osteoblast and fibroblasts",
                  "RbpjWt_C02-Mesenchymal cells"  ,"RbpjLOF_C02-Mesenchymal cells" ,
                  "RbpjWt_C03-Endothelial" , "RbpjLOF_C03-Endothelial",
                  "RbpjWt_C04-Chondrocytes"  ,   "RbpjLOF_C04-Chondrocytes"  ,
                  "RbpjWt_C05-Epithelial",   "RbpjLOF_C05-Epithelial",
                  "RbpjWt_C06-Fetal brain (CNS)"     ,  "RbpjLOF_C06-Fetal brain (CNS)" ,
                  "RbpjWt_C07-SmoothMuscle", "RbpjLOF_C07-SmoothMuscle"  ,
                  "RbpjWt_C08-Blood progenitors" , "RbpjLOF_C08-Blood progenitors"   ,  
                  "RbpjWt_C09-Myocytes" , "RbpjLOF_C09-Myocytes"       , 
                  "RbpjWt_C10-Glial cells"           , "RbpjLOF_C10-Glial cells"       ,
                  "RbpjWt_C11-Hepatocytes" ,"RbpjLOF_C11-Hepatocytes"   ,
                  "RbpjWt_C12-Erythrocytes"   ,"RbpjLOF_C12-Erythrocytes"       ,
                  "RbpjWt_C13-NeuralCrest"  ,"RbpjLOF_C13-NeuralCrest"  ,
                  "RbpjWt_C14-Lymphatic progenitors"    ,  "RbpjLOF_C14-Lymphatic progenitors" , 
                  "RbpjWt_C15-Cardiomyocytes"    ,"RbpjLOF_C15-Cardiomyocytes"       ,      
                  "RbpjWt_C16-PNS"   ,"RbpjLOF_C16-PNS" )


levels(ob) <- c(  "RbpjWt_C0-NeuralCrest"                 , "RbpjLOF_C0-NeuralCrest"                ,
                  "RbpjWt_C01-Fetal brain (CNS)"         , "RbpjLOF_C01-Fetal brain (CNS)"         ,
                  "RbpjWt_C02-PNS"                        , "RbpjLOF_C02-PNS"                        ,
                  "RbpjWt_C03-Glial cells"                 , "RbpjLOF_C03-Glial cells"                ,
                  "RbpjWt_C04-Epithelial"  , "RbpjLOF_C04-Epithelial"                ,
                  "RbpjWt_C05-Lymphatic progenitors"       , "RbpjLOF_C05-Lymphatic progenitors"      ,
                  "RbpjWt_C06-Endothelial"                 , "RbpjLOF_C06-Endothelial"               ,
                  "RbpjWt_C07-Cardiomyocytes"             , "RbpjLOF_C07-Cardiomyocytes"             ,
                  "RbpjWt_C08-Erythrocytes"               , "RbpjLOF_C08-Erythrocytes"              ,
                  "RbpjWt_C09-Blood progenitors"           , "RbpjLOF_C09-Blood progenitors"          ,
                  "RbpjWt_C10-Hepatocytes"                , "RbpjLOF_C10-Hepatocytes"               ,
                  "RbpjWt_C11-Myocytes"                   , "RbpjLOF_C11-Myocytes"                   ,
                  "RbpjWt_C12-SmoothMuscle"                , "RbpjLOF_C12-SmoothMuscle"              ,
                  "RbpjWt_C13-Chondrocytes"               , "RbpjLOF_C13-Chondrocytes"              ,
                  "RbpjWt_C14-Mesenchymal cells"          , "RbpjLOF_C14-Mesenchymal cells"         ,
                  "RbpjWt_C15-Mesenchymal cells2"          , "RbpjLOF_C15-Mesenchymal cells2"         ,
                  "RbpjWt_C16-Osteoblast and fibroblasts" , "RbpjLOF_C16-Osteoblast and fibroblasts" )






ncells = ncol(Ec)

for (i in 1:ncells){
  b = as.character(Ec@meta.data$Condition[i])
  a = as.character(Ec@meta.data$endothelial_labels3[i])
  
  Ec@meta.data$Condition_Clustering2[i] = paste0(b,"_",a)
  
}


Idents(Ec) <- 'Condition_Clustering2'
 levels(Ec)<-c("RbpjWt_C00-mesenchymal-proliferative","RbpjLOF_C00-mesenchymal-proliferative" ,
 "RbpjWt_C01-endoEMT" ,"RbpjLOF_C01-endoEMT" ,
 "RbpjWt_C02-BBB_ECs" , "RbpjLOF_C02-BBB_ECs",
 "RbpjWt_C03-Liver_ECs" , "RbpjLOF_C03-Liver_ECs" )


levels(Ec) <- c("RbpjWt_C00" ,"RbpjLOF_C00" ,
                "RbpjWt_C01" , "RbpjLOF_C01",
                "RbpjWt_C02-endoEMT", "RbpjLOF_C02-endoEMT" ,
                "RbpjWt_C03-BBB_ECs" ,"RbpjLOF_C03-BBB_ECs"  ,
                "RbpjWt_C04-Liver_ECs", "RbpjLOF_C04-Liver_ECs")

#____________________________________________________________________________________________________________

ncells = ncol(smc)

for (i in 1:ncells){
  b = as.character(smc@meta.data$Condition[i])
  a = as.character(smc@meta.data$smc_labels[i])
  
  smc@meta.data$Condition_Clustering2[i] = paste0(b,"_",a)
  
}

Idents(smc)<- 'Condition_Clustering2'

levels(smc)<-c("RbpjWt_C00","RbpjLOF_C00",
               "RbpjWt_C01","RbpjLOF_C01",
               "RbpjWt_C02","RbpjLOF_C02",
               "RbpjWt_C03","RbpjLOF_C03",
               "RbpjWt_C04","RbpjLOF_C04")

saveRDS(smc, file='SmoothMuscleCells_v3.rds')

Idents(smc)<- 'smc_labels'
levels(smc) <- c("C04" ,"C03", "C02","C01", "C00" )
#____________________________________________________________________________________________________________


ncells = ncol(brain)

for (i in 1:ncells){
  b = as.character(brain@meta.data$Condition[i])
  a = as.character(brain@meta.data$fbrain_labels[i])
  
  brain@meta.data$Condition_Clustering[i] = paste0(b,"_",a)
  
}

Idents(brain) <- 'Condition_Clustering'

levels(brain) <- c("RbpjWt_C00","RbpjLOF_C00","RbpjWt_C01","RbpjLOF_C01",
                   "RbpjWt_C02","RbpjLOF_C02","RbpjWt_C03","RbpjLOF_C03",
                   "RbpjWt_C04" ,"RbpjLOF_C04")

levels(brain) <- c("C00", "C01", "C02", "C03","C04")

saveRDS(brain, file='BrainC06Cluster_v2.rds')
#____________________________________________________________________________________________________________


Idents(BloodP) <- 'Condition_Clustering'

levels(BloodP) <- c("RbpjWt_C00-Mieloid_Progenitors" , "RbpjLOF_C00-Mieloid_Progenitors",
"RbpjWt_C01-Macrophages"  , "RbpjLOF_C01-Macrophages",
"RbpjWt_C02-Microglia" ,  "RbpjLOF_C02-Microglia" ,
"RbpjWt_C03-Monocytes"    ,"RbpjLOF_C03-Monocytes"  )

Idents(BloodP) <- 'blood_labels'
levels(BloodP) <-  c("C00-Mieloid_Progenitors","C01-Macrophages","C02-Microglia" ,"C03-Monocytes" )
levels(BloodP) <-  c("C03-Monocytes","C02-Microglia" ,"C01-Macrophages","C00-Mieloid_Progenitors")
#____________________________________________________________________________________________________________
ncells = ncol(epi)

for (i in 1:ncells){
  b = as.character(epi@meta.data$Condition[i])
  a = as.character(epi@meta.data$epithelial_labels[i])
  
  epi@meta.data$Condition_Clustering[i] = paste0(b,"_",a)
  
}
saveRDS(epi, file="Epithelial_v2.rds")
#order levels
Idents(epi) <- 'Condition_Clustering'
levels(epi) <- c("RbpjWt_C00" , "RbpjLOF_C00" ,"RbpjWt_C01","RbpjLOF_C01", "RbpjWt_C02","RbpjLOF_C02" , "RbpjWt_C03","RbpjLOF_C03"   )
Idents(epi) <- 'epithelial_labels'
levels(epi) <- c( "C00" ,"C01","C02", "C03")
levels(epi) <- c( "C03","C02","C01", "C00" )


#______________________________________________________________________________________________________________


ncells = ncol(osfi)

for (i in 1:ncells){
  b = as.character(osfi@meta.data$Condition[i])
  a = as.character(osfi@meta.data$osfi_labels[i])
  
  osfi@meta.data$Condition_Clustering[i] = paste0(b,"_",a)
  
}

Idents(osfi) <- 'Condition_Clustering'
levels(osfi) <- c("RbpjWt_C00","RbpjLOF_C00",
                  "RbpjWt_C01","RbpjLOF_C01",
                  "RbpjWt_C02" ,"RbpjLOF_C02",
                  "RbpjWt_C03" ,"RbpjLOF_C03")

#____________C15-Mesenchymal cells2_____________

ncells = ncol(mesenc)

for (i in 1:ncells){
  b = as.character(mesenc@meta.data$Condition[i])
  a = as.character(mesenc@meta.data$mesenc_labels[i])
  
  mesenc@meta.data$Condition_Clustering[i] = paste0(b,"_",a)
  
}

Idents(mesenc) <- 'Condition_Clustering'
levels(mesenc) <- c("RbpjWt_C00","RbpjLOF_C00",
                  "RbpjWt_C01","RbpjLOF_C01",
                  "RbpjWt_C02" ,"RbpjLOF_C02",
                  "RbpjWt_C03" ,"RbpjLOF_C03")

#____________C14-Mesenchymal cells_____________

ncells = ncol(mesenc1)

for (i in 1:ncells){
  b = as.character(mesenc1@meta.data$Condition[i])
  a = as.character(mesenc1@meta.data$mesenc1_labels[i])
  
  mesenc1@meta.data$Condition_Clustering[i] = paste0(b,"_",a)
  
}

Idents(mesenc1) <- 'Condition_Clustering'
levels(mesenc1) <- c("RbpjWt_C00","RbpjLOF_C00",
                    "RbpjWt_C01","RbpjLOF_C01",
                    "RbpjWt_C02" ,"RbpjLOF_C02",
                    "RbpjWt_C03" ,"RbpjLOF_C03",
                    "RbpjWt_C04" ,"RbpjLOF_C04",
                    "RbpjWt_C05" ,"RbpjLOF_C05")

#____________C13_Chondrocytes_____________

ncells = ncol(chond)

for (i in 1:ncells){
  b = as.character(chond@meta.data$Condition[i])
  a = as.character(chond@meta.data$chond_labels[i])
  
  chond@meta.data$Condition_Clustering[i] = paste0(b,"_",a)
  
}

Idents(chond) <- 'Condition_Clustering'
levels(chond) <- c("RbpjWt_C00","RbpjLOF_C00",
                     "RbpjWt_C01","RbpjLOF_C01",
                     "RbpjWt_C02" ,"RbpjLOF_C02",
                     "RbpjWt_C03" ,"RbpjLOF_C03")

#____________C11_Myocytes_____________

ncells = ncol(myoc)

for (i in 1:ncells){
  b = as.character(myoc@meta.data$Condition[i])
  a = as.character(myoc@meta.data$myoc_labels[i])
  
  myoc@meta.data$Condition_Clustering[i] = paste0(b,"_",a)
  
}

Idents(myoc) <- 'Condition_Clustering'
levels(myoc) <- c("RbpjWt_C00","RbpjLOF_C00",
                   "RbpjWt_C01","RbpjLOF_C01",
                   "RbpjWt_C02" ,"RbpjLOF_C02",
                   "RbpjWt_C03" ,"RbpjLOF_C03")

#____________C05-Glial cells _____________

ncells = ncol(glial)

for (i in 1:ncells){
  b = as.character(glial@meta.data$Condition[i])
  a = as.character(glial@meta.data$glial_labels[i])
  
  glial@meta.data$Condition_Clustering[i] = paste0(b,"_",a)
  
}

Idents(glial) <- 'Condition_Clustering'
levels(glial) <- c("RbpjWt_C00","RbpjLOF_C00",
                  "RbpjWt_C01","RbpjLOF_C01",
                  "RbpjWt_C02" ,"RbpjLOF_C02",
                  "RbpjWt_C03" ,"RbpjLOF_C03")


#____________C00 neural crest _____________

ncells = ncol(neurl)

for (i in 1:ncells){
  b = as.character(neurl@meta.data$Condition[i])
  a = as.character(neurl@meta.data$neurl_labels[i])
  
  neurl@meta.data$Condition_Clustering[i] = paste0(b,"_",a)
  
}

Idents(neurl) <- 'Condition_Clustering'
levels(neurl) <- c("RbpjWt_C00","RbpjLOF_C00",
                   "RbpjWt_C01","RbpjLOF_C01",
                   "RbpjWt_C02" ,"RbpjLOF_C02")

