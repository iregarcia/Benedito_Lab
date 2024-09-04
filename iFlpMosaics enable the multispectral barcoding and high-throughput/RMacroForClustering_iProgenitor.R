path<-choose.dir(default = "", caption = "Select folder with csv file for heatmap analysis")

path2<- paste(path, "//ReadyFor_Seurat.csv", sep = "")


table2 <- read.csv(path2)


ta2<-as.data.frame(table2)

library(stringr)


table$Celltype <- str_c(table$Reporter1.positivity, '_', table$Reporter2.positivity)  


table[table=='No_Yes'] <- 'YFP'
table[table=='Yes_No'] <- 'Tom'
table[table=='No_No'] <- 'Neg'
table[table=='Yes_Yes'] <- 'Double'

install.packages("fpc")
install.packages("dbscan")
library(dbscan)
library(fpc)
library(tidyr)
library(ggthemes)
library(dbplyr)

pos_tb <- subset(table, table$Celltype!='Neg')


n <- unique(pos_tb$Image)


im_tb11 <- subset(pos_tb, pos_tb$Image==n[1])
im_tb11 <- drop_na(im_tb11)

df <- im_tb11[,1:2]


    
      dbscan::kNNdistplot(df, k = 10)
      abline(h = 0.4, lty = 2)
      
      res1 <- fpc::dbscan(df, eps =8, MinPts = 2) #Generates the unsupervised clusters
      
      res1 <- as.factor(res1$cluster)

im_tb11$Clone <- res1


ggplot(im_tb11, aes(x=Relative.X, y=Relative.Y, colour=Clone)) + 
  geom_point(size=2)+ ggtitle('clones identified by distribution')



im_tb1$Cluster <- str_c(im_tb1$Clone, '_', 'Im1')
im_tb2$Cluster <- str_c(im_tb2$Clone, '_', 'Im2')
im_tb3$Cluster <- str_c(im_tb3$Clone, '_', 'Im3')
im_tb4$Cluster <- str_c(im_tb4$Clone, '_', 'Im4')
im_tb5$Cluster <- str_c(im_tb5$Clone, '_', 'Im5')
im_tb6$Cluster <- str_c(im_tb6$Clone, '_', 'Im6')
im_tb7$Cluster <- str_c(im_tb7$Clone, '_', 'Im7')
im_tb8$Cluster <- str_c(im_tb8$Clone, '_', 'Im8')
im_tb9$Cluster <- str_c(im_tb9$Clone, '_', 'Im9')
im_tb10$Cluster <- str_c(im_tb10$Clone, '_', 'Im10')
im_tb11$Cluster <- str_c(im_tb11$Clone, '_', 'Im11')

im_tb12$Cluster <- str_c(im_tb11$Clone, '_', 'Im12')
im_tb13$Cluster <- str_c(im_tb13$Clone, '_', 'Im13')
im_tb14$Cluster <- str_c(im_tb14$Clone, '_', 'Im14')
im_tb15$Cluster <- str_c(im_tb15$Clone, '_', 'Im15')
im_tb16$Cluster <- str_c(im_tb16$Clone, '_', 'Im16')
im_tb17$Cluster <- str_c(im_tb17$Clone, '_', 'Im17')





total_tb <- rbind(im_tb1, im_tb2, im_tb3, im_tb4, im_tb5, im_tb6, im_tb7, im_tb8, im_tb9, im_tb10, im_tb11 )

iPro_all <- rbind(im_tb1, im_tb2, im_tb3, im_tb4, im_tb5, im_tb6, im_tb7, im_tb8, im_tb9, im_tb10, im_tb11,im_tb12, im_tb13, im_tb14, im_tb15, im_tb16, im_tb17)

write.csv(x=total_tb, file='AllClonesAGOMLitter.csv' )

write.csv(x=iPro_all, file='AllClonesAGOM_AGMG_AGLF.csv' )

#Filter clone 0 in every image (does not constitute a real clone)

tot_nozero <- total_tb %>% filter(!grepl(pattern = '0_', x=total_tb$Cluster))


iPro_NOC0 <- iPro_all %>% filter(!grepl(pattern = '0_', x=iPro_all$Cluster))

write.csv(x=tot_nozero, file='AllClonesNoZERO_AGOMLitter.csv')

write.csv(x=iPro_NOC0, file='AllClonesNoZeroAGOM_AGMG_AGLF.csv' )

#Count cells per cluster
percluster <- tot_nozero %>% count(tot_nozero$Cluster)
perclusterall <- iPro_NOC0 %>% count(iPro_NOC0$Cluster)
#Select clusters with less than 30 cells
percluster30 <- subset(percluster, percluster$n<30)
percluster30all <- subset(perclusterall, perclusterall$n<30)


percluster30_3_all <- subset(percluster30all, percluster30all$n>3)


#Subset clusters smaller than 30
iPro_less30 <- tot_nozero[tot_nozero$Cluster %in% percluster30$`tot_nozero$Cluster`, ]

iPro_less30all <- iPro_NOC0[iPro_NOC0$Cluster %in% percluster30all$`iPro_NOC0$Cluster`, ]
iPro_less30More3 <- iPro_NOC0[iPro_NOC0$Cluster %in% percluster30_3_all$`iPro_NOC0$Cluster`, ]


write.csv(iPro_less30More3, file='iProgenitor_ClonesLess30More3Cells.csv')
write.csv(x=iPro_less30, file='ClonesLess30Cells_AGOMLitter.csv')

write.csv(x=iPro_less30all, file='iProgenitor_ClonesLess30Cells_AGOM_AGMG_AGLF.csv')

im <- unique(iPro_less30all$Image)



#plot clones smaller than 30 cells

ggplot(iPro_less30all, aes(x=Relative.X, y=Relative.Y, colour=Cluster)) + 
  geom_point(size=1)+
  ggtitle('clones identified by distribution')+
  facet_wrap(~Image)





#Barplots organize by size





ggplot(iPro_NOC0, aes(x=forcats::fct_infreq(Cluster), fill=Celltype))+ 
  geom_histogram(stat ='count')+
  scale_fill_manual(values=c('orange', 'red3', 'green3')) + 
  xlab('Clone ID')+
  ylab('Nº cells per clone')+
  ggtitle('Cells per clone and cell type',subtitle = 'FlpoERT2 x iProgenitor x iFlpMosaic')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1), panel.grid = element_blank())
  

ggplot(iPro_less30all, aes(x=forcats::fct_infreq(Clone), fill=Celltype))+ 
  geom_histogram(stat ='count')+
  scale_fill_manual(values=c('orange', 'red3', 'green3')) + 
  xlab('Clone ID')+
  ylab('Nº cells per clone')+
  ggtitle('Cells per clone and cell type',subtitle = 'FlpoERT2 x iProgenitor x iFlpMosaic')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size=5),panel.grid = element_blank())+
  facet_wrap(~Image)


ggplot(iPro_less30, aes(x=Cluster,  fill=Celltype) ) +geom_bar(position='fill')+
  ggtitle('Contribution of cell type per clone')+
  xlab('Clone ID')+
  ylab('Relative %')+
  scale_fill_manual(values=c('orange', 'red3', 'green3'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1),panel.grid = element_blank())

ggplot(tot_nozero, aes(x=Cluster,  fill=Region) ) +geom_bar(position='fill')+
  ggtitle('Contribution of regions per clone')+
  xlab('Clone ID')+
  ylab('Relative %')+
  scale_fill_gdocs()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1),panel.grid = element_blank())







# to get images separately:




x=15

a <- subset(iPro_NOC0, iPro_NOC0$Image==im[x])

a2 <- subset(a, a$Clone!=1)

#To generate cluster contour
hull_data <-  a %>%
  group_by(Clone) %>%
  slice(chull(Relative.X, Relative.Y))

ggplot(a, aes(x=Relative.X, y=Relative.Y))+ 
   ggtitle(im[x])+
  theme_classic()+
  geom_polygon(data = hull_data, alpha = 0.7, aes(fill=Clone))+
  scale_fill_brewer(palette = 'Set3')+
  geom_point(size=1.5, aes(colour=Celltype))+
  scale_colour_manual(values=c('orange2','red2','green3'))

ggplot(a, aes(x=Relative.X, y=Relative.Y, colour=Clone))+ geom_point() 


str_c('ClonesTomYFP_',im[x],'.pdf')

###Plot all the clusters of all images in one graph####

hull_data2 <-  iPro_less30 %>%
  group_by(Cluster) %>%
  slice(chull(Relative.X, Relative.Y))

ggplot(iPro_less30, aes(x=Relative.X, y=Relative.Y))+ 
  ggtitle('All images analysed (clones/clusters)')+
  theme_classic()+
  geom_polygon(data = hull_data2, alpha = 0.7, aes(fill=Cluster))+
  geom_point(size=0.5)+
  facet_wrap(~Image)


ggplot(iPro_less30, aes(x=Relative.X, y=Relative.Y, colour=Cluster))+ geom_point()

#####Dendogram #####
library(ggraph)
library(igraph)
library(tidyverse)


dist <- dist( a[,1:2], diag=TRUE)
hc <- hclust(dist)


plot(hc)

a$ID <- rownames(a)

ggplot(subset(a, a$Clone=='4'), aes(x=Relative.X, y=Relative.Y,label=ID, colour=Celltype))+ 
  geom_point()+
  geom_text(hjust=0, vjust=0)




######numer of cells per clone and image####

y=15

p2 <- subset(iPro_NOC0, iPro_NOC0$Image==im[y] & iPro_NOC0$Clone!=3)
p <- subset(iPro_NOC0, iPro_NOC0$Image==im[y])


freqtable <- data.frame(table(p$Cluster))

ggplot(p, aes(x=Relative.X, y=Relative.Y, colour=Clone))+ geom_point()           
           
ggplot(p, aes(x=forcats::fct_infreq(Cluster), fill=Celltype))+ 
  geom_bar(stat='count')+
  scale_fill_manual(values=c('orange2','red3', 'green3')) + 
  xlab('Clones')+
  ylab('% each  cell type per clone')+
  ggtitle('Cells per clone and cell type',subtitle = im[y])+
  geom_text(data = freqtable, aes(x = Var1, y = Freq, label = Freq, size=50, vjust=-0.5),inherit.aes = FALSE)

str_c('Totalnumberclones_',im[y],'.pdf')



ggplot(subset(pos_tb, pos_tb$Sample_Type!='AGMG'), aes(x=Sample_Type, fill=Celltype))+ geom_bar(position='fill')+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3, hjust=1, size=30),panel.grid = element_blank())+
  scale_fill_manual(values=c('orange', 'red3', 'green3'))

ggplot(pos_tb, aes(x=Cluster, fill=Celltype))+ geom_bar(position='fill')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size=10),panel.grid = element_blank())+
  scale_fill_manual(values=c('orange', 'red3', 'green3'))

