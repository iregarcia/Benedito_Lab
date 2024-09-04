
#colours and variables for plots####

clus_col <- c('C0-Mesenchymal cells2'= '#ba6917', 'C01-Osteoblast and fibroblasts'= '#f76dde',
              'C02-Mesenchymal cells'='#f45328','C03-Endothelial'='#985fe0',
              'C04-Chondrocytes'='#ccced1', 'C05-Epithelial'='#04c5fd', 'C06-Fetal brain (CNS)'='red2',
              'C07-SmoothMuscle'='orange', 'C08-Blood progenitors'='#579316',
              'C09-Myocytes'='#a7baa8', 'C10-Glial cells'='#0cfa66',
              'C11-Hepatocytes'='#0ba4fb', 'C12-Erythrocytes'='#fc0387', 
              'C13-NeuralCrest'='#606980', 'C14-Lymphatic progenitors'='#f16915',
              'C15-Cardiomyocytes'='#aff9fc', 'C16-PNS'='#c8fa6a')


clus_col_new <- c('C15-Mesenchymal cells2'= '#ba6917', 'C16-Osteoblast and fibroblasts'= '#f76dde',
                  'C14-Mesenchymal cells'='#f0d8f8','C06-Endothelial'='#985fe0',
                  'C13-Chondrocytes'='#fbf405', 'C04-Epithelial'='#04c5fd',
                  'C01-Fetal brain (CNS)'='red2',
                  'C12-SmoothMuscle'='orange', 'C09-Blood progenitors'='#579316',
                  'C11-Myocytes'='#a7baa8', 'C03-Glial cells'='#0cfa66',
                  'C10-Hepatocytes'='#0ba4fb', 'C08-Erythrocytes'='#fc0387', 
                  'C0-NeuralCrest'='#606980', 'C05-Lymphatic progenitors'='#f16915',
                  'C07-Cardiomyocytes'='#aff9fc', 'C02-PNS'='#c8fa6a')

clus_rb_cols <- c('#606980',  'red2', '#c8fa6a', '#0cfa66', '#04c5fd', '#f16915', '#985fe0',
                  '#aff9fc', '#fc0387', '#579316',   '#0ba4fb', '#a7baa8','orange','#fbf405',
                  '#f0d8f8', '#ba6917', '#f76dde')

Bestholtz_palette <- c("#DEDAD6","#FEE392","#FEC44E","#FE9929","#ED6F11","#CC4C17","#993411","#65260C")


notch <- c('Hes1', 'Hes5', 'Hey1','Heyl','Hey2', 'Rbpj')

end_cols <- c('C00'='#befb0e','C01'='#c778f1','C02-endoEMT'='#ffa600','C03-BBB_ECs'='#479608','C04-Liver_ECs'='#3175de')

only_cols <- c( '#ba6917',  '#f76dde', '#f45328','#985fe0', '#ccced1','#04c5fd',
               'red2', 'orange', '#579316', '#a7baa8', '#0cfa66','#0ba4fb',
               '#fc0387','#606980','#f16915', '#aff9fc', '#c8fa6a')
blood_cols <-  c('C00-Mieloid_Progenitors'='#ffa600','C01-Macrophages'='#e718b5','C02-Microglia'='#18e799','C03-Monocytes'='#3d18e7')

fbrain_cols <- c('C00'='#f1490e','C01'='#16923d','C02'='#e8ce28', 'C03'='#0e7cf1', 'C04'='#410c7c') 

notchligands <- c('Notch1', 'Notch2', 'Notch3', 'Notch4', 'Dll1', 'Dll2', 'Dll3', 'Dll4', 'Jag1', 'Jag2')

end_cols3 <- c('C00-mesenchymal-proliferative'='#cb4cf3', 'C01-endoEMT'='#ffa600',
               'C02-BBB_ECs'='#479608','C03-Liver_ECs'='#3175de')
epi_cols <- c('C00'='#f9d520','C01'='#6bc115','C02'='#f93420','C03'='#a912f9')
cols3 <- c('#cb4cf3', '#ffa600','#479608','#3175de')
arterial_genes <- c('Gja4', 'Bmx', 'Gja5', 'Cxcr4')
smc_cols <- c('C00'='#9ef57b','C01'='#f69521','C02'='#1899e7', 'C03'='#f5e710', 'C04'='#d67bf5') 

#Dim plots ####


DimPlot(ob, reduction = 'umap_.19.dims',
        group.by = 'new_labels_rb',
        cols=clus_col_new,
        label = T, label.box = T, repel = T, pt.size =0.7) + ggtitle('scRNAseq Clusters (Rbpj floxed)')

DimPlot(ob, reduction = 'umap_.19.dims',
        group.by = 'new_labels',
        cols=clus_col,
        label = T, label.box = T, repel = T, pt.size =0.9) + ggtitle('scRNAseq Clusters (Rbpj floxed)')

p2 <- DimPlot(ob, reduction = 'umap_.19.dims',
        group.by = 'new_labels_rb', split.by = 'Condition',
        cols=clus_col_new,pt.size =0.7) + ggtitle('scRNAseq Clusters (Rbpj floxed)')

p2<-DimPlot(ob, reduction = 'umap_.19.dims', shuffle = TRUE,
        group.by = 'Condition', cols=c('#80e13d','#f83d14'),pt.size =0.7) + ggtitle('scRNAseq Clusters (Rbpj floxed)')

DimPlot(ob, reduction = 'umap_.19.dims',
        group.by = 'Condition', cols=c('#80e13d','#f83d14'),pt.size =0.9) +
  ggtitle('scRNAseq Clusters (Rbpj floxed)')


p<- DimPlot(ob, reduction = 'umap_.19.dims', split.by = 'Condition',
        group.by = 'Condition', cols=c('#80e13d','#f83d14'),pt.size =0.7) + ggtitle('scRNAseq Clusters (Rbpj floxed)')

p<-DimPlot(Ec, reduction = 'ec_umap20',group.by = 'endothelial_labels3',
        cols=end_cols3,pt.size =1) + ggtitle('Endothelial subclusters')

p<-DimPlot(Ec_down, reduction = 'ec_umap20',group.by = 'endothelial_labels3', split.by = 'Condition',
        cols=end_cols3,pt.size =0.7) + ggtitle('Endothelial subclusters')

DimPlot(Ec, reduction = 'ec_umap20',group.by = 'Condition',shuffle = T,
        cols=c('#80e13d','#f83d14'),pt.size =1.5) + ggtitle('Endothelial subclusters')

p<-DimPlot(blood_downs, reduction = 'blood_res0.2',group.by = 'blood_labels',split.by = 'Condition',
        cols=blood_cols,pt.size =0.7) + ggtitle('Blood progenitor subclusters')

p10<-DimPlot(BloodP, reduction = 'blood_res0.2',group.by = 'blood_labels',
        cols=blood_cols,pt.size =1.5) + ggtitle('Blood progenitor subclusters')

p10<-DimPlot(BloodP, reduction = 'blood_res0.2',group.by = 'Condition', shuffle = T,
        cols=c('#80e13d','#f83d14'),pt.size =1.5) + ggtitle('Blood progenitor subclusters per cell type')


p<- DimPlot(fbrain, reduction="fbrain_res0.25", pt.size = 3, group.by = 'fbrain_labels', cols = fbrain_cols)+
  ggtitle('Neuron subclusters')

p<- DimPlot(fbrain, reduction="fbrain_res0.25", pt.size = 3, group.by = 'fbrain_labels', cols = fbrain_cols, split.by = 'Condition')+
  ggtitle('Neuron subclusters')



p1<-DimPlot(brain, reduction = "fbrain_res0.25",group.by = 'Condition',shuffle = T,
        cols=c('#80e13d','#f83d14'),pt.size =2) + ggtitle('Brain CNS - Cell distribution')

p2<-DimPlot(brain, reduction="fbrain_res0.25", pt.size = 1.5, group.by = 'fbrain_labels', cols = fbrain_cols, split.by = 'Condition')+
  ggtitle('Brain CNS subclusters')


p<-DimPlot(epi_down, reduction="epi_umap18", pt.size = 0.7, group.by = 'epithelial_labels', cols =epi_cols, split.by = 'Condition')+
  ggtitle('Epithelial subclusters')

p<-DimPlot(epi, reduction="epi_umap18", pt.size = 1.5, group.by = 'Condition', cols =c('#80e13d','#f83d14'), shuffle=T)+
  ggtitle('Epithelial subclusters')
  
p<-DimPlot(epi, reduction="epi_umap18", pt.size = 1.5, group.by = 'epithelial_labels', cols =epi_cols)+
  ggtitle('Epithelial subclusters')


t<-DimPlot(smc_down, pt.size = 0.7, reduction="smc_0.15" , group.by = 'smc_labels', cols =smc_cols, split.by = 'Condition')+
  ggtitle('SMC clustering')

p<-DimPlot(smc_down, pt.size = 0.7, reduction="smc_15" , group.by = 'Condition',cols =c('#80e13d','#f83d14'), shuffle=T)+
  ggtitle('Epithelial subclusters')

p<-DimPlot(mesenc_down, pt.size = 0.9, reduction="mesenc_res0.3" , group.by = 'mesenc_labels',cols =mesenc_cols, split.by = 'Condition')

p<-DimPlot(mesenc_down, pt.size = 0.9, reduction="mesenc_res0.3" , group.by = 'mesenc_labels',cols =mesenc_cols, split.by = 'Condition')

p<-DimPlot(mesenc1_down, pt.size = 0.9, reduction="mesenc1_res0.2" , group.by = 'mesenc1_labels',cols =mesenc1_cols, split.by = 'Condition')

p<-DimPlot(chond_down, pt.size = 0.9, reduction="chond_res0.2" , group.by = 'chond_labels',cols =chond_cols, split.by = 'Condition')

p<-DimPlot(myoc_down, pt.size = 0.9, reduction="myoc_res0.2" , group.by = 'myoc_labels',cols =myoc_cols, split.by = 'Condition')


p<-DimPlot(glial_down, pt.size = 0.9, reduction="glial_res0.2" , group.by = 'glial_labels',cols =glial_cols, split.by = 'Condition')


p<-DimPlot(neurl_down, pt.size = 0.9, reduction="neurl_res0.25" , group.by = 'neurl_labels',cols =neurl_cols, split.by = 'Condition')



  
#Dot plots ####

Idents(ob) <- 'new_labels'
Idents(ob) <- factor(x = Idents(ob), levels = sort(levels(ob))) # to order levels
#______________________________________________________________________________
df <- data.frame(
  x = rep(c(3,8,13,18,23, 28, 33,38,43,48,53, 58,63,68,73,78,83)),
  y = rep(c(0),1),
  z = factor(rep(1:17)),
  color = clus_rb_cols)

DotPlot(ob, group.by = 'new_labels_rb', features=Top5genes, col.min = 0, dot.scale = 4, scale = T)+ 
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size =8), legend.text  = element_text(size = 10),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=ob$new_labels_rb)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = clus_rb_cols)


#______________________________________________________________________________
s<-seq(from=1.5, to =34.5, by=2)
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:17)),
  color = clus_rb_cols)

p10<-DotPlot(ob,  features=c('Cxcr4', 'Gja4', 'Gja5'), col.min = 0, dot.scale = 4, scale = T)+
  coord_flip()+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5, face = "italic", size = 5),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), axis.ticks.y = ) +
    scale_x_discrete(label=ob$new_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = clus_rb_cols)
#______________________________________________________________________________
s<-seq(from=1.5, to =8.5, by=2)
cols3 <- c("#cb4cf3", "#ffa600", "#479608", "#3175de" )
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = cols3)

p15<-DotPlot(Ec, features=c('Gja4', 'Bmx', 'Gja5', 'Cxcr4', 'Mki67', 'Cdkn1a'), col.min = 0, dot.scale = 8, scale = T)+
  coord_flip()+ggtitle('General Arterial embryonic markers')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels2)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = cols3)

DotPlot(Ec, features=c('Runx1', 'Runx2','Myc', 'Jag1'), col.min = 0, dot.scale = 8, scale = T)+
  coord_flip()+ggtitle('Cell cycle markers')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels2)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = cols3)

#_________________________________________________________________________________________

s<-seq(from=1.5, to =8.5, by=2)
cols3 <- c("#cb4cf3", "#ffa600", "#479608", "#3175de" )
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = cols3)

DotPlot(Ec, features=c('Gja4', 'Sat1', 'Hey1', 'Serpinf1', 'Igfbp3', 'Vegfc'), col.min = 0, dot.scale = 8, scale = T)+
  coord_flip()+ggtitle('Arterial embryonic markers')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels3)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = cols3)


p<-DotPlot(Ec, features=c('Esm1', 'Apln', 'Kcne3', 'Cdkn1a', 'Rbpj', 'Hes1'), col.min = 0, dot.scale = 6, scale = T)+
  coord_flip()+ggtitle('Notch ligands and receptors')+
  theme(strip.text.x = element_blank(), title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = cols3)

p13 <- DotPlot(Ec, features=notch, col.min = 0, dot.scale = 6, scale = T)+
  coord_flip()+ggtitle('Notch ligands and receptors')+
  theme(strip.text.x = element_blank(), title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = cols3)



#______________________________________________________________________________
s<-seq(from=5.5, to =35.5, by=10)
onencols <- c("#cb4cf3", "#ffa600", "#479608", "#3175de" )
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:4)),
  color = onencols)

t <- DotPlot(Ec_ctrl, features=Top10genes_ec, col.min = 0, dot.scale = 7, scale = T)+
  ggtitle('Markers per cluster (using Ctrl cells)')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)
#______________________________________________________________________________
s<-seq(from=3, to =18, by=5)
onencols <- c("#cb4cf3", "#ffa600", "#479608", "#3175de" )
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:4)),
  color = onencols)

t <- DotPlot(Ec, features=Top5genes_ec, col.min = 0, dot.scale = 7, scale = T)+
  ggtitle('Top 5 Markers per cluster (Using all cells)')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)
#______________________________________________________________________________
s<-seq(from=3, to =23, by=5)
onencols <- c("#9ef57b", "#f69521" ,"#1899e7" ,"#f5e710" ,"#d67bf5" )
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:5)),
  color = onencols)

p<-DotPlot(smc, features=Top5genes_smc, col.min = 0, dot.scale = 7, scale = T)+
  ggtitle('Top5 Markers per cluster')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)



s<-seq(from=1.5, to =10.5, by=2)
cols3 <- c("#9ef57b", "#f69521" ,"#1899e7" ,"#f5e710" ,"#d67bf5")
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:5)),
  color = cols3)

p<-DotPlot(smc, features=DEGs_C_Rbpj_more, col.min = 0, dot.scale = 8, scale = T)+
  coord_flip()+ggtitle('Selected markers')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels3)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = cols3)



#______________________________________________________________________________
s<-seq(from=3, to =22.5, by=5)
onencols <- c('#f9d520','#6bc115','#f93420','#a912f9')
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:4)),
  color = onencols)

p<-DotPlot(epi, features=Top5genes_epi, col.min = 0, dot.scale = 6, scale = T)+
  ggtitle('Top5 Markers per cluster in epithelial cells')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=epi$epithelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)


s<-seq(from=1.5, to =8.5, by=2)
cols3 <- c('#f9d520','#6bc115','#f93420','#a912f9')
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = cols3)

p<-DotPlot(epi, features=c('Cdkn1a', 'Mki67'), col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('General Notch ligands')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=epi$epithelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = cols3)
#_____________________________________________________________________________________

s<-seq(from=1, to =4.5, by=1)
cols3 <- c('#f9d520','#6bc115','#f93420','#a912f9')
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = cols3)

p<-DotPlot(epi, features=c('EMT1', 'Mesenchymal1', 'Epithe1', 'Hallmark_EMT1'), col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('Markers for epithelia to mesenchymal transition')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=epi$epithelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = cols3)

#_____________________________________________________________________________________
s<-seq(from=1.5, to =8.5, by=2)
cols3 <- c('#f9d520','#6bc115','#f93420','#a912f9')
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = cols3)

p<-DotPlot(epi, features=Top5_Rbpj_more, col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('General Notch ligands')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=epi$epithelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = cols3)


#______________________________________________________________________________
s<-seq(from=5.5, to =45.5, by=10)
onencols <- c("#befb0e"   ,  "#c778f1"   ,  "#ffa600"  ,   "#479608"   ,  "#3175de" )
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:5)),
  color = onencols)

p<-DotPlot(Ec, features=Top10genes_ec, col.min = 0, dot.scale = 6, scale = T)+
  ggtitle('Markers per cluster (using Ctrl cells)')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)

#______________________________________________________________________________
s<-seq(from=3, to =20, by=5)
onencols <- c( "#ffa600"    ,"#e718b5", "#18e799","#3d18e7" )
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:4)),
  color = onencols)

levels(BloodP)<- c("C00-Mieloid_Progenitors","C01-Macrophages"  , "C02-Microglia"  , "C03-Monocytes"                  )

t<-DotPlot(BloodP, features=Top5genes_BloodP, col.min = 0, dot.scale = 5, scale = T)+
  ggtitle('Top5 Markers per cluster')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size =11),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)


#______________________________________________________________________________
s<-seq(from=1.5, to =8.5, by=2)
onencols <- c( "#ffa600"    ,"#e718b5", "#18e799","#3d18e7" )
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = onencols)

bcells <-c('Pik3r1', 'Pik3cd', 'Btk')
macromono <- c('Itgam', 'Adgre1', 'Cd86', 'Mrc1')
cellcycle <- c('Mki67', 'Cdkn1a', 'Myc', 'Mycn', 'Odc1')



p<-DotPlot(BloodP, features=c(DEGs_C_Rbpj_more, DEGs_C_Rbpj_less), col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('Selected marker genes')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle =90, hjust = 1, vjust =0,  face = "italic", size =5),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)


#______________________________________________________________________________
s<-seq(from=3, to =25, by=5)
onencols <- c( "#f1490e" ,"#16923d" ,"#e8ce28" ,"#0e7cf1" ,"#410c7c" )


df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:5)),
  color = onencols)



t<-DotPlot(brain, features=Top5genes_brain, col.min = 0, dot.scale = 5, scale = T)+
  ggtitle('Top5 Markers per cluster')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size =11),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)
#______________________________________________________________________________
s<-seq(from=1.5, to =10.5, by=2)
onencols <- c( "#f1490e" ,"#16923d" ,"#e8ce28" ,"#0e7cf1" ,"#410c7c" )
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:5)),
  color = onencols)


p<-DotPlot(brain, features=c(DEGs_C_Rbpj_more,DEGs_C_Rbpj_less), col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('MostDEGs identified per cluster (WT vs RbpjLOF))')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle =90, hjust = 1, vjust =0,  face = "italic", size =5),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=Ec$endothelial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)

#______________________________________________________________________________
s<-seq(from=3, to =18, by=5)
onencols <- c('#f9d520','#26c556','#f93420','#a912f9' )
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:4)),
  color = onencols)

levels(osfi) <- c("C03","C02","C01","C00")
p<-DotPlot(osfi, features=Top5genes_osfi, col.min = 0, dot.scale = 7, scale = T)+
  ggtitle('Top 5 Markers per cluster (Using all cells)')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=osfi$osfi_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)


s<-seq(from=1.5, to =8.5, by=2)
onencols <- c('#f9d520','#26c556','#f93420','#a912f9' )
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = onencols)

Idents(osfi) <- 'Condition_Clustering'
p<-DotPlot(osfi, features=c('Cdkn1a', 'Mki67','Hes1','Rbpj'), col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('Selected marker genes')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle =90, hjust = 1, vjust =0,  face = "italic", size =5),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=osfi$osfi_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)

#C15-Mesenchymal cells2 ______________________________________________________________________________
s<-seq(from=3, to =18, by=5)
onencols <- c('#f37c13','#26c556','#1361f3','#f92aaa')
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:4)),
  color = onencols)

Idents(mesenc) <- 'mesenc_labels'
levels(mesenc) <- c("C03","C02","C01","C00")
p<-DotPlot(mesenc, features=Top5genes_mesenc, col.min = 0, dot.scale = 7, scale = T)+
  ggtitle('Top 5 Markers per cluster (Using all cells)')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=mesenc$mesenc_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)


s<-seq(from=1.5, to =8.5, by=2)
onencols <- c('#f37c13','#26c556','#1361f3','#f92aaa')
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = onencols)

selected_genes <- c('Cdkn1a', 'Mki67','Hes1','Rbpj')
DEGS_Rbpj <- c("Ulk4" ,  "Rbpj"  , "Hes1" ,"Col1a1")
p<-DotPlot(osfi, features=DEGS_Rbpj, col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('Selected marker genes')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle =90, hjust = 1, vjust =0,  face = "italic", size =5),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=osfi$osfi_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)

#C14-Mesenchymal cells1 ______________________________________________________________________________

s<-seq(from=3, to =28, by=5)
onencols <- mesenc1_cols <- c('#e8150f','#26c556', '#0ff6d3','#eff60f','#c94bfc', '#f6890f')
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:6)),
  color = onencols)

Idents(mesenc1) <- 'mesenc1_labels'
levels(mesenc1) <- c("C05","C04", "C03","C02","C01","C00")
p<-DotPlot(mesenc1, features=Top5genes_mesenc1, col.min = 0, dot.scale = 7, scale = T)+
  ggtitle('Top 5 Markers per cluster (Using all cells)')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=mesenc1$mesenc1_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)


s<-seq(from=1.5, to =13, by=2)
onencols <-  c('#e8150f','#26c556', '#0ff6d3','#eff60f','#c94bfc', '#f6890f')
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:6)),
  color = onencols)

selected_genes <- c('Cdkn1a', 'Mki67','Hes1','Rbpj')
DEGS_Rbpj <- c("Fabp7" , "H2ac24" ,"Ulk4"  , "Notch2", "Col1a1" )

p<-DotPlot(mesenc1, features=selected_genes, col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('Selected marker genes')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle =90, hjust = 1, vjust =0,  face = "italic", size =5),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=mesenc1$mesenc1_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)

#_____________C13_Chondrocytes________________________________________________________________________
s<-seq(from=3, to =18, by=5)
onencols <-c('#e8150f','#26c556','#f3db0e','#9312e8')
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:4)),
  color = onencols)

Idents(chond)<-'chond_labels'
levels(chond) <- c("C03","C02","C01","C00")
p<-DotPlot(chond, features=Top5genes_chond, col.min = 0, dot.scale = 7, scale = T)+
  ggtitle('Top 5 Markers per cluster (Using all cells)')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=chond$chond_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)


s<-seq(from=1.5, to =8.5, by=2)
onencols <- c('#e8150f','#26c556','#f3db0e','#9312e8')
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = onencols)

Idents(chond) <- 'Condition_Clustering'
cond_DEGS <- c("H1f5"  , "Vtn"  ,  "Ulk4"   ,"Dcn"    ,"Col1a1")
p<-DotPlot(chond, features=cond_DEGS, col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('Selected marker genes')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle =90, hjust = 1, vjust =0,  face = "italic", size =5),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=chond$chond_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)
#_____________C11-Myocytes________________________________________________________________________
s<-seq(from=3, to =18, by=5)
onencols <-c('#f7a5fa','#26c556', '#f97c06','#9312e8')
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:4)),
  color = onencols)

Idents(myoc)<-'myoc_labels'
levels(myoc) <- c("C03","C02","C01","C00")
p<-DotPlot(myoc, features=Top5genes_myoc, col.min = 0, dot.scale = 7, scale = T)+
  ggtitle('Top 5 Markers per cluster (Using all cells)')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=myoc$myoc_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)


s<-seq(from=1.5, to =8.5, by=2)
onencols <- c('#f7a5fa','#26c556', '#f97c06','#9312e8')
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = onencols)

Idents(myoc) <- 'Condition_Clustering'
levels(myoc) <- c("RbpjWt_C00","RbpjLOF_C00", "RbpjWt_C01","RbpjLOF_C01", "RbpjWt_C02" ,"RbpjLOF_C02","RbpjWt_C03" ,"RbpjLOF_C03")
cond_DEGS <- 
p<-DotPlot(myoc, features=selected_genes, col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('Selected marker genes')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle =90, hjust = 1, vjust =0,  face = "italic", size =5),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=myoc$myoc_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)

#_____________C05 Glial cells ________________________________________________________________________
s<-seq(from=3, to =18, by=5)
onencols <- c('#f7a5fa','#1206f9','#f90678','#30c507')
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:4)),
  color = onencols)

Idents(glial)<-'glial_labels'
levels(glial) <- c("C03","C02","C01","C00")
p<-DotPlot(glial, features=Top5genes_glial, col.min = 0, dot.scale = 7, scale = T)+
  ggtitle('Top 5 Markers per cluster (Using all cells)')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=glial$glial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)


s<-seq(from=1.5, to =8.5, by=2)
onencols <- c('#f7a5fa','#1206f9','#f90678','#30c507')
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:4)),
  color = onencols)

Idents(glial) <- 'Condition_Clustering'
levels(glial) <- c("RbpjWt_C00","RbpjLOF_C00", "RbpjWt_C01","RbpjLOF_C01", "RbpjWt_C02" ,"RbpjLOF_C02","RbpjWt_C03" ,"RbpjLOF_C03")
cond_DEGS <- 
  p2<-DotPlot(glial, features=t, col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('Selected marker genes')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle =90, hjust = 1, vjust =0,  face = "italic", size =5),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=glial$glial_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)

#_____________C00 Neural crest ________________________________________________________________________
s<-seq(from=3, to =15, by=5)
onencols <-c('#27a2d4','#0925c7', '#fb901f')
df <- data.frame(
  x = s,
  y = rep(c(0)),
  z = factor(rep(1:3)),
  color = onencols)

Idents(neurl)<-'neurl_labels'
levels(neurl) <- c("C02","C01","C00")
p<-DotPlot(neurl, features=Top5genes_neurl, col.min = 0, dot.scale = 7, scale = T)+
  ggtitle('Top 5 Markers per cluster (Using all cells)')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=neurl$neurl_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)


s<-seq(from=1.5, to =5.5, by=2)
onencols <- c('#27a2d4','#0925c7', '#fb901f')
df <- data.frame(
  x = rep(c(0)),
  y = s,
  z = factor(rep(1:3)),
  color = onencols)

Idents(neurl) <- 'Condition_Clustering'
levels(neurl) <- c("RbpjWt_C00","RbpjLOF_C00", "RbpjWt_C01","RbpjLOF_C01", "RbpjWt_C02" ,"RbpjLOF_C02")
cond_DEGS <- 
  p2<-DotPlot(neurl, features=selected_genes, col.min = 0, dot.scale = 7, scale = T)+
  coord_flip()+ggtitle('Selected marker genes')+
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle =90, hjust = 1, vjust =0,  face = "italic", size =5),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size =8), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label=neurl$neurl_labels)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = onencols)

# Feature plots ####


FeaturePlot(ob, features = "Crem", pt.size = 2, split.by='Condition')

FeaturePlot(smc, features = 'HA1', reduction="smc_15", pt.size = 1)+
  scale_colour_gradientn(colours = Bestholtz_palette)

p<-FeaturePlot(Ec_down, features = c('Kcne3', 'Esm1', 'Apln'), reduction="ec_umap20", pt.size = 1, split.by = 'Condition', cols =Bestholtz_palette )

# Violin plots ####

VlnPlot(Ec, features = c('Gja5', 'Gja4', 'Efnb2'), group.by = 'endothelial_labels3', cols = end_cols3)
t <-VlnPlot(Ec, features = c('Cdkn1a', 'Rbpj','Hes1' ),
        group.by = 'endothelial_labels3',
        split.by =  'Condition',cols = c('#80e13d', '#f83d14')) 

p<-VlnPlot(brain, features = c('Cdkn1a', 'Rbpj','Hes1' ),
        group.by = 'fbrain_labels',
        split.by =  'Condition',cols = c('#80e13d', '#f83d14')) 
p<-VlnPlot(brain, features = c('Mki67','Gja4', 'Gja5' ),
           group.by = 'fbrain_labels',
           split.by =  'Condition',cols = c('#80e13d', '#f83d14')) 


p17 <- VlnPlot(BloodP, features = 'Macroph1',  cols =  blood_cols)+
  ggtitle('Macrophage markers module')+
  ylab('Gene cluster enrichment')+ 
  theme(axis.title.x = element_blank(), legend.position = 'none')

p17<-VlnPlot(BloodP, features = 'Gja4',  cols =  blood_cols)+
  ggtitle('Monocyte markers module')+
  ylab('Gene cluster enrichment')+ 
  theme(axis.title.x = element_blank(), legend.position = 'none')

VlnPlot(ob, features = c('Gja4', 'Bmx'),  group.by = 'Condition')+
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=0, hjust=0.5))


p15<-VlnPlot(ob, features = 'Rbpj',  group.by = 'new_labels', split.by = 'Condition',  cols =  c('#80e13d', '#f83d14'))+
  theme(axis.title.x = element_blank(),legend.position = 'none',
        axis.text.x = element_text(angle=90, vjust=0.5, size=10))

p16<-VlnPlot(ob, features = 'Hes1',  group.by = 'Condition',  cols =  c('#80e13d', '#f83d14'))+
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=0, hjust=0.5, size=15))

VlnPlot(Ec, features = 'Tbx20', group.by = 'Condition', cols =  c('#80e13d', '#f83d14'))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0, hjust = 0.5, size=15), legend.position='none')+
  ylab('Expression level')

VlnPlot(Ec, features = 'Slc22a8', group.by = 'endothelial_labels3',
        cols = end_cols3)+
  theme(axis.title.x = element_blank(),legend.position = 'none',
        axis.text.x = element_text(angle=0, hjust = 0.5, size=7))+
  ylab('Expression level')


VlnPlot(epi, features = 'Rbpj', group.by = 'Condition',
        cols = c('#80e13d', '#f83d14'))+
  theme(axis.title.x = element_blank(),legend.position = 'none',
        axis.text.x = element_text(angle=0, hjust = 0.5, size=15))+
  ylab('Expression level')



p2<-VlnPlot(epi, features = 'Bcc7', group.by = 'Condition',
            cols = c('#80e13d', '#f83d14'))+
  theme(axis.title.x = element_blank(),legend.position = 'none',
        axis.text.x = element_text(angle=0, hjust = 0.5, size=15))+
  ylab('Expression level')

p1<-VlnPlot(epi, features = c('Rbpj', 'Hes1', 'Mki67' ), group.by = 'epithelial_labels', split.by = 'Condition',
        cols = c('#80e13d', '#f83d14'))

p2<-VlnPlot(epi, features = c('Cdkn1a', 'Gja4', 'Gja5'), group.by = 'epithelial_labels', split.by = 'Condition',
        cols = c('#80e13d', '#f83d14'))

VlnPlot(smc, features = c('Tbx18','S1pr1'), group.by = 'smc_labels', cols = smc_cols)



p1<-VlnPlot(BloodP, features = c('Rbpj', 'Hes1', 'Mki67' ), group.by = 'blood_labels', split.by = 'Condition',
            cols = c('#80e13d', '#f83d14'))

p2<-VlnPlot(BloodP, features = c('Cdkn1a', 'Gja4', 'Gja5'), group.by = 'blood_labels', split.by = 'Condition',
            cols = c('#80e13d', '#f83d14'))



p1<-VlnPlot(smc, features = c('Rbpj', 'Hes1', 'Mki67' ), group.by = 'smc_labels', split.by = 'Condition',
            cols = c('#80e13d', '#f83d14'))

p2<-VlnPlot(smc, features = c('Cdkn1a', 'Gja4', 'Gja5'), group.by = 'smc_labels', split.by = 'Condition',
            cols = c('#80e13d', '#f83d14'))


p2<-VlnPlot(Ec, features = c('Gja4', 'Gja5'), group.by = 'Condition',
            cols = c('#80e13d', '#f83d14'))

p<-VlnPlot(brain, features = c('Gja1', 'Gjc1'), group.by = 'fbrain_labels', split.by = 'Condition',cols = c('#80e13d', '#f83d14'))

# Bar plots ####


dittoBarPlot(Ec, var='endothelial_labels3', group.by = 'Condition', color.panel = end_cols3, x.reorder = c(2,1)) + 
  ggtitle('Contribution per ECs cluster')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p<- dittoBarPlot(ob, var='new_labels', group.by = 'Condition', color.panel = clus_col, x.reorder = c(2,1)) + 
  ggtitle('Contribution per ECs cluster')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p2<-dittoBarPlot(ob, var='new_labels_rb', group.by = 'Condition', color.panel = clus_col_new, x.reorder = c(2,1)) + 
  ggtitle('Contribution per ECs cluster')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p10<-dittoBarPlot(ob, var='Condition', group.by = 'new_labels_rb', color.panel = c( '#f83d14','#80e13d')) + 
  ggtitle('Contribution per ECs cluster')+
  theme(axis.text.x =  element_text(angle=90, hjust = 0.5,  vjust = 0, face = "italic", size = 10), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p10<- dittoBarPlot(Ec_same, var='Condition', group.by = 'endothelial_labels3', color.panel = c( '#f83d14','#80e13d')) + 
  ggtitle('Contribution per ECs cluster')+
  theme(axis.text.x =  element_text(angle=90, hjust = 1,  vjust = 0.5, face = "italic", size = 10), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p11<-dittoBarPlot(blood_downs, var='Condition', group.by = 'blood_labels', color.panel = c( '#f83d14','#80e13d')) + 
  ggtitle('Contribution per Blood progenitor cluster (Equal Ctrl and Mut number of cells)')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0.5, face = "italic", size = 13), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p21<-dittoBarPlot(BloodP, var='blood_labels', group.by = 'Condition', color.panel = blood_cols, x.reorder = c(2,1)) + 
  ggtitle('Contribution per Blood progenitor cluster')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0.5, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p<-dittoBarPlot(brain_down, var='Condition', group.by = 'fbrain_labels', color.panel = c( '#f83d14','#80e13d')) + 
  ggtitle('Contribution per Neuron clusters')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0.5, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p3<-dittoBarPlot(brain, var='fbrain_labels', group.by = 'Condition', color.panel = fbrain_cols) + 
  ggtitle('Contribution per Neuron clusters')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0.5, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p10<- dittoBarPlot(epi_down, var='Condition', group.by = 'epithelial_labels', color.panel = c( '#f83d14','#80e13d')) + 
  ggtitle('Contribution per Epithelial clusters (equal nº WT and RbpjLOF)')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0.5, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p10<-dittoBarPlot(epi, var='epithelial_labels', group.by = 'Condition', color.panel = epi_cols, x.reorder = c(2,1)) + 
  ggtitle('Contribution per genotype of Epithelial clusters')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0.5, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p<- dittoBarPlot(smc, var='smc_labels', group.by = 'Condition', color.panel = smc_cols, x.reorder = c(2,1)) + 
  ggtitle('Contribution per genotype of SMC clusters')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0.5, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))


p<-dittoBarPlot(smc_down, var='Condition', group.by = 'smc_labels', color.panel = c( '#f83d14','#80e13d')) + 
  ggtitle('Contribution per genotype of SMC clusters (same nº Ctrl and KO cells)')
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0.5, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))
  
  
p<- dittoBarPlot(osfi_down, var='osfi_labels', group.by = 'Condition', color.panel = osfi_cols, x.reorder = c(2,1)) + 
    ggtitle('Contribution per osteoclast/fibroblast cluster')+
    theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0, face = "italic", size = 15), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=15))

p<-dittoBarPlot(mesenc_down, var='mesenc_labels', group.by = 'Condition', color.panel = mesenc_cols, x.reorder = c(2,1)) + 
  ggtitle('Contribution per meshencymal cluster')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p<-dittoBarPlot(mesenc1_down, var='mesenc1_labels', group.by = 'Condition', color.panel = mesenc1_cols, x.reorder = c(2,1)) + 
  ggtitle('Contribution per meshencymal cluster')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p<-dittoBarPlot(chond_down, var='chond_labels', group.by = 'Condition', color.panel = chond_cols, x.reorder = c(2,1)) + 
  ggtitle('Contribution per meshencymal cluster')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))


p<-dittoBarPlot(myoc_down, var='myoc_labels', group.by = 'Condition', color.panel = myoc_cols, x.reorder = c(2,1)) + 
  ggtitle('Contribution per meshencymal cluster')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p<-dittoBarPlot(glial_down, var='glial_labels', group.by = 'Condition', color.panel = glial_cols, x.reorder = c(2,1)) + 
  ggtitle('Contribution per meshencymal cluster')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))

p<-dittoBarPlot(neurl_down, var='neurl_labels', group.by = 'Condition', color.panel = neurl_cols, x.reorder = c(2,1)) + 
  ggtitle('Contribution per meshencymal cluster')+
  theme(axis.text.x =  element_text(angle=0, hjust = 0.5,  vjust = 0, face = "italic", size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15))


# Bar plots for GSEA NEW ####
p15<- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES, fill=)) +
  scale_fill_gradient2(low='red3', mid='yellow2', high='green3')+
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
  
  ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=NES)) +
    scale_fill_gradient2(low='red3', mid='yellow2', high='green3')+ scale_x_continuous(limits=c(-2,2))+
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways NES from GSEA") + 
    theme_minimal() 
  


p16<-ggplot(fgsplot, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  scale_fill_gradient2(low='red3', mid='yellow2', high='green3')+
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA - C08 Blood progenitors") + 
  theme_minimal()

p17<- ggplot(fgsplot, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA - C08 Blood progenitors") + 
  theme_minimal()


p1<-ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  scale_fill_gradient2(low="#10c62c", mid='yellow2', high='#e82212')+
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

p1<-ggplot(fgsplot, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  scale_fill_gradient2(low="#10c62c", mid='yellow2', high='#e82212', )+
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA (Epithelial)") + 
  theme_minimal()

# Percentage of cells expressing a gene ####

#from dotplot
a <- p5$data
a<- as.data.frame(a)
p6<-  ggplot(a, aes(x=features.plot, y=pct.exp, fill=id))+ geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values= clus_col)+
  ggtitle('% of Gene+ cells in every cluster')+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(size=15), axis.title.y = element_text(size=15))+
  ylab('Percentage Gene+ cells (separated by cluster)')+ labs(fill='Cluster ID')


a <- p7$data
a<- as.data.frame(a)
p8<-  ggplot(a, aes(x=features.plot, y=pct.exp, fill=id))+ geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values= end_cols)+
  ggtitle('% of Gene+ cells in every EC cluster')+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(size=15), axis.title.y = element_text(size=15))+
  ylab('Percentage Gene+ cells (Endothelial cluster)')+ labs(fill='Cluster ID')

a <- p15$data
a<- as.data.frame(a)

p<-ggplot(a2, aes(x=features.plot, y=pct.exp, fill=id))+ geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values= c('#80e13d','#f83d14','#80e13d','#f83d14','#80e13d','#f83d14','#80e13d','#f83d14'))+
  ggtitle('% of positive cells expressing the gene')+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(size=15), axis.title.y = element_text(size=15))+
  ylab('Percentage Gene+ cells')+ labs(fill='Cluster ID')


p10<-DotPlot(BloodP, features=c('Gja4', 'Gja5','Gjc1', 'Gja1', 'Gjb2', 'Gjb6'), col.min = 0, dot.scale = 5, scale = T)+
  ggtitle('Markers per cluster (using Ctrl cells)')

a <- p10$data

p<-ggplot(a, aes(x=features.plot, y=pct.exp, fill=id))+ geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values= blood_cols)+
  ggtitle('Percentage Gene+ cells')+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(size=15), axis.title.y = element_text(size=15))+
  ylab('Percentage Gene+ cells')+ labs(fill='Cluster ID')



p10<-DotPlot(epi, features=c('EMT1', 'Mesenchymal1'), col.min = 0, dot.scale = 5, scale = T)+
  ggtitle('Markers per cluster (using Ctrl cells)')

a <- p10$data

p<-ggplot(a, aes(x=features.plot, y=avg.exp, fill=id))+ geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values= c('#80e13d','#f83d14','#80e13d','#f83d14','#80e13d','#f83d14','#80e13d','#f83d14'))+
  ggtitle('Average gene expression')+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(size=10), axis.title.y = element_text(size=10),
        title = element_text(size=15))+
  ylab('Average gene expression')+ labs(fill='Cluster ID')

p1<-ggplot(a, aes(x=features.plot, y=pct.exp, fill=id))+ geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values= c('#80e13d','#f83d14','#80e13d','#f83d14','#80e13d','#f83d14','#80e13d','#f83d14'))+
  ggtitle('Percentage positive cells')+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(size=10), axis.title.y = element_text(size=10),
        title = element_text(size=15))+
  ylab('Gene+ cells %')+ labs(fill='Cluster ID')

# To save as PDF in folder ####


setwd('S:/LAB_RB/LAB/Irene/iFlpMosaics Paper Summer 2022/scRNAseq_Embryos/scRNA_Embryos_2023/CLUSTERS/C13_Chondrocyte')
cairo_pdf(width=8, height =7, filename = 'GSEA_chondrocyte.pdf')
p
dev.off()


