

# ANALYSIS OF RAW OUTPUT OBTAINED FROM FIJI USING CELL COUNTER TO LABEL 
# RECOMBINED CELLS ACCORDING TO THEIR COLOR CODE (see color codes in line 37)

# install.packages("plotly")
# install.packages('roperators')
#install.packages("rlang")
require(roperators)
library(ggplot2)
library(plotly)
library(gapminder)
library(ggplotify)
library(ggpubr)
library(plotrix)
library(png)

#opening and reading the data

path <- choose.dir()
setwd(path)
file <- file.choose()
#data <- read.csv(file)


# Read with custom delimiter
data <- read.csv(file, sep = ";")
data$X.micron. <- as.numeric(as.character(data$X.micron.))
data$Y.micron. <- as.numeric(as.character(data$Y.micron.))

# get file name without extension

name <- tools::file_path_sans_ext(file)

color <- unlist(list("red", "green","red", "red", "green", "green","green","red","red", "red", "green", "green","green", "green", "red","red", "red", "green","green","red","red","green", "green","red","red", "red","nan","red", "nan","red", "red","red", "red","red","green", "red", "green","red", "red", "red", "red","red","red", "red","red","red", "red", "red", "red","red","red","green","green","red")) 

esp <- unlist(list("Mutant-Tomato+DAPI+","WT-YFP+DAPI+","Mutant-Tomato+H2B-GFP+","Mutant-Tomato+H2B-Cer+","WT-YFP+H2B-GFP+","WT-YFP+H2B-Cer+","WT-YFP+H2B-Cherry+","Mutant-Tomato+YFP+","Mutant-Tomato+Tfp1 HA+","Mutant-Tomato+DAPI-"
         ,"WT-YFP+DAPI-","WT-H2B-Cherry+","WT-H2B-GFP+","WT-H2B-Cer+","Mutant-Tomato+YFP+H2B-Cherry+","Mutant-Tomato+YFP+H2B-GFP+","Mutant-Tomato+YFP+H2B-Cer+","WT-YFP+H2B-Cherry+H2B-Cer+","WT-YFP+H2B-Cherry+H2B-GFP+","Mutant-Tomato+Tfp1 HA+H2B-GFP+"
         ,"Mutant-Tomato+Tfp1 HA+H2B-Cer+","WT-H2B-Cherry+H2B-GFP+","WT-H2B-Cherry+H2B-Cer+","Mutant-Tomato+YFP+DAPI-","Mutant-Tomato+Tfp1 HA+DAPI-","Mutant-Tomato+H2B-Cherry+","Nan","Mutant-Tomato+YFP+H2B-Cherry+H2B-GFP+","Nan","Mutant-Tomato+Tfp1 HA+YFP+H2B-Cherry+","Mutant-Tomato+Tfp1 HA+YFP+DAPI-","Mutant-Tomato+Tfp1 HA+H2B-Cherry+","Mutant-Tomato+H2B-Cherry+H2B-Cer+","Mutant-Tomato+H2B-GFP+H2B-Cer+","WT-YFP+H2B-GFP+H2B-Cer+","Mutant-Tomato+H2B-Cherry+H2B-GFP+H2B-Cer+","WT-YFP+H2B-Cherry+H2B-GFP+H2B-Cer+","Mutant-Tomato+YFP+H2B-Cherry+H2B-Cer+","Mutant-Tomato+YFP+H2B-GFP+H2B-Cer+","Mutant-Tomato+YFP+H2B-Cherry+H2B-GFP+H2B-Cer+"
         ,"Mutant-Tomato+Tfp1 HA+H2B-Cherry+H2B-GFP+","Mutant-Tomato+Tfp1 HA+H2B-Cherry+H2B-Cer+","Mutant-Tomato+Tfp1 HA+H2B-GFP+H2B-Cer","Mutant-Tomato+Tfp1 HA+H2B-Cherry+H2B-GFP+H2B-Cer+","Mutant-Tomato+Tfp1 HA+YFP+","Mutant-Tomato+Tfp1 HA+YFP+H2B-GFP+","Mutant-Tomato+Tfp1 HA+YFP+H2B-Cer+","Mutant-Tomato+Tfp1 HA+YFP+H2B-Cherry+H2B-GFP+","Mutant-Tomato+Tfp1 HA+YFP+H2B-Cherry+H2B-Cer+","Mutant-Tomato+Tfp1 HA+YFP+H2B-GFP+H2B-Cer+"
         ,"Mutant-Tomato+Tfp1 HA+YFP+H2B-Cherry+H2B-GFP+H2B-Cer+","WT-H2B-Cherry+H2B-GFP+","WT-H2B-Cherry+H2B-GFP+H2B-Cer+","Mutant-Tomato+H2B-Cherry+H2B-GFP"))  

n <- unlist(1:54)  # number possible celltypes

del <- unique(data$Type)

species = data.frame(unlist(n), unlist(esp),unlist(color))
colnames(species) <- c("Num","Species_Name", "Colour")

species <- species[species$Num %in% del, ]

# subset data according to cell_fluo_markers

for (i in species[ , 1]) {

  assign(paste0('data_type',i), subset(data, data$Type==i))

}

#-------------------------------------------------------------------
# FUNCTION TO DETERMINE CLONALITY
#-------------------------------------------------------------------
Clone_ID <- function(data, type, max_dist){
  # function determines whether cells with same color code belong
  # to the same clone or not and labels them accordingly (eg. 
  # cells labelled 1.1, 1.2 and 1.3 belong to 3 different clones
  # but have color code type 1)
  
  # VARIABLES
  # -------------------------------------------------------
  # data -> data.frame,
  #         parameters for cells with a given color code
  # type -> integer number,
  #         defines the color code
  # max_dist -> integer number,
  #             maximal distance (in microns) allowed
  #                 between clonal cells
  
  # RETURNS
  # -------------------------------------------------------
  # vector -> vector,
  #           with the clone IDs for each cell in dataframe
  #               input
  # -------------------------------------------------------
  
  if (nrow(data)>0){
    vector <- c(rep(0, nrow(data)))
    ID = type
    for (i in 1:nrow(data)){
      SizeClone = 1
      if (vector[i] == 0){
        ID = ID + 0.001
        vector[i] = ID
        count = 1
        for (j in 1:nrow(data)){
          if (i != j){    #new addition
            x1 <- data$X.micron.[i]
            y1 <- data$Y.micron.[i]
            x2 <- data$X.micron.[j]
            y2 <- data$Y.micron.[j]
            dist <- sqrt((x1-x2)^2+(y1-y2)^2)
            if (dist < max_dist){
              vector[j] = ID
              count = count + 1
            }
          }
        }
        while (SizeClone < count){
          initialSizeClone = count
          for (k in 1:nrow(data)){
            for (l in 1:nrow(data)){
              if ((vector[l] ==0) & (vector[k] == ID)){
                x1 <- data$X.micron.[k]
                y1 <- data$Y.micron.[k]
                x2 <- data$X.micron.[l]
                y2 <- data$Y.micron.[l]
                dist <- sqrt((x1-x2)^2+(y1-y2)^2)
                if (dist < max_dist){
                  vector[l] = ID
                  count = count + 1
                }
              }
            }
          }
          SizeClone = initialSizeClone
          
        }
      }
    }
  }
  return(vector)
}
#-------------------------------------------------------------------
#-------------------------------------------------------------------


## Get clonal info on all celltypes
##----------------------------------------------------------------------

a <- Clone_ID(data_type1,1,20)
b <- Clone_ID(data_type2,2,20)
for (i in species[ , 1]) {
  
  assign(paste0('Type',i,'_clone_ID'), Clone_ID(eval(parse(text = paste0('data_type',i))), i, 15)) # Last number is the max distance for clones
  
}


# join CloneIDs in a single vector

CloneID <-c()
for (i in species[ , 1]) {
  if (nrow(eval(parse(text = paste0('data_type',i)))) > 0){
    
    CloneID <- c(CloneID, eval(parse(text = paste0('Type',i,'_clone_ID'))))
    
  }
}

# Create new table with relevant data
Comp_data <-c()
for (i in species[ , 1]) {
  if (nrow(eval(parse(text = paste0('data_type',i)))) > 0){
    
    Comp_data <- rbind(Comp_data, eval(parse(text = paste0('data_type',i))))
    
  }
}

Comp_data <- cbind(Comp_data, CloneID)
Comp_data <- Comp_data[c('Type', 'X.micron.', 'Y.micron.', 'CloneID')]
colnames(Comp_data) <- c("Cell_Type", "X_microns", "Y_microns", "CloneID")
Comp_data <- as.data.frame(Comp_data)


#------------------------------------------------------------------
# FUNCTION TO DETERMINE DISTANCE TO CLOSEST CELL WITHIN CLONE
#-------------------------------------------------------------------

Dist_closest <- function(data){
  # function determines the distance of each cell to its 
  # nearest neighbour of the same clone, and store the 
  # coordinates of this closest clonal cell
  
  # VARIABLES
  # -------------------------------------------------------
  # Comp_data -> data.frame,
  #              with parameters CloneID, FluoMarkers, and 
  #                 X,Y coordinates for all cells
  
  # RETURNS
  # -------------------------------------------------------
  # dist_frame -> data.frame,
  #               with distance (in microns) to nearest clonal
  #                 cell and its coordinates
  # -------------------------------------------------------
  Dist_Closest_Cell <- c(rep(0,length(CloneID)))
  X_Closest_Cell <- c(rep(0,length(CloneID)))
  Y_Closest_Cell <- c(rep(0,length(CloneID)))
  
  for (i in 1:length(CloneID)){
    min_dist = 200
    for (j in 1:length(CloneID)){
      if ((i != j) & (CloneID[i] == CloneID[j])){
        x1 <- data$X_microns[i]
        y1 <- data$Y_microns[i]
        x2 <- data$X_microns[j]
        y2 <- data$Y_microns[j]
        dist <- sqrt((x1-x2)^2+(y1-y2)^2)
        if (dist < min_dist){
          min_dist = dist
          Dist_Closest_Cell[i] = min_dist
          X_Closest_Cell[i] = data$X_microns[j]
          Y_Closest_Cell[i] = data$Y_microns[j]
        }
      }
    }
  }  
  dist_frame <- cbind(data, Dist_Closest_Cell, X_Closest_Cell, Y_Closest_Cell)
  return (dist_frame)
}
#-------------------------------------------------------------------
#-------------------------------------------------------------------


## Complete clonal data
##----------------------------------------------------------------------
Comp_data_new <- as.data.frame(Dist_closest(Comp_data))



# subset data according to Cell_Type

for (i in species[ , 1]) {
  
  assign(paste0('cell_type',i), subset(Comp_data_new, Comp_data_new$Cell_Type==i))
  
}


#----------------------------------------------------------------------
# FUNCTION TO DETERMINE PARAMETERS FOR CENTROID
#----------------------------------------------------------------------

Centroid_dispersion <- function(cell_type){
  # VARIABLES
  #--------------------------------------------------------------------
  # cell_type -> data.frame,
  #               with clonal info from 1 type of color combination
  
  if (nrow(cell_type) > 0){
    # order the dataset by CloneID so that these are grouped
    cell_type <- cell_type[order(cell_type$CloneID), ]
    ID_size <- as.data.frame(table(cell_type$CloneID))
    colnames(ID_size) <- c("Clone_ID", "CloneSize")
    # calculate centroid and dispersion coordinates
    start = 1
    x_mean <- c()
    y_mean <- c()
    for (i in 1:nrow(ID_size)){
      max_perp_dist = 0
      # store x coordinates of polygon
      x_poly <-c()
      # store y coordinates of polygon
      y_poly <-c()
      end = start + ID_size$CloneSize[i] - 1
      if (ID_size$CloneSize[i] > 1){
        for (j in start:end){
          x_poly <-c(x_poly, cell_type$X_microns[j])
          y_poly <-c(y_poly, cell_type$Y_microns[j])
        }
        start = end + 1
        # coordinates of centroid (average of x-coord and y-coord)
        x_mean <- c(x_mean, mean(x_poly))
        y_mean <- c(y_mean, mean(y_poly))
        
        
      }else if(ID_size$CloneSize[i] == 1){
        x_mean <- c(x_mean, cell_type$X_microns[start])
        y_mean <- c(y_mean, cell_type$Y_microns[start])
        start = start + 1
        
      }
    }
    
    
    CentroidDisp <- as.data.frame(cbind(x_mean, y_mean))
    return(cbind(ID_size, CentroidDisp))
  }
}

#-------------------------------------------------------------------
#-------------------------------------------------------------------

for (i in species[ , 1]) {
  
  assign(paste0('Cen_disp_',i), Centroid_dispersion(eval(parse(text = paste0('cell_type',i)))))
  
}


pdf("Clone size frequencies of all cell types.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

counts <- c()

for (i in species[ , 1]) {
  
  if (nrow(eval(parse(text = paste0('Cen_disp_',i))))>0) {
    barplot(table(eval(parse(text = paste0('Cen_disp_',i,'$CloneSize')))), col= species$Colour[species$Num==i], main = species$Species_Name[species$Num==i])
  }
  
  assign(paste0('counts_',i), as.data.frame(table(eval(parse(text = paste0('Cen_disp_',i,'$CloneSize'))))))
  assign(paste0('counts_p_',i), as.data.frame(prop.table(table(eval(parse(text = paste0('Cen_disp_',i,'$CloneSize')))))*100))
  assign(paste0('counts_',i), cbind(rep(species$Species_Name[species$Num==i], nrow(eval(parse(text = paste0('counts_',i))))),  eval(parse(text = paste0('counts_',i))), eval(parse(text = paste0('counts_p_',i)))))

  counts <- rbind(counts,eval(parse(text = paste0('counts_',i))))  
  
}

colnames(counts) <- c('CloneType','CloneSize','AbsFrequency', 'CloneSize','RelatFrequency')


# Closing the graphical device
dev.off()

mutant <- species[species$Colour=="red",]
WT <- species[species$Colour=="green",]


mut_cen_disp <- lapply(paste0('Cen_disp_',mutant[,1]),get)
WT_cen_disp <- lapply(paste0('Cen_disp_',WT[,1]),get)

# Histogram of frequencies of clone size ALL Mutant cells

mut_clonesize <- lapply(mut_cen_disp, "[", , "CloneSize")

if (nrow(rbind(mut_cen_disp))>0){
  
  pdf("Total mutant clones_Frequency.pdf",         # File name
      width = 11.69, height = 8.27,   # Width and height in inches
      paper = "A4r")          # Paper size
  barplot(table(unlist(c(mut_clonesize))), col='red',
          main='Mutant clones')

  barplot(prop.table(table(unlist(c(mut_clonesize))))*100,
          col='red',
          main='Mutant clones_percentage')
  # Closing the graphical device
  dev.off()
}

WT_clonesize <- lapply(WT_cen_disp, "[", , "CloneSize")

# Histogram of frequencies of clone size WT cells
if (nrow(rbind(WT_cen_disp))>0){
  
  pdf("Total WT clones_Frequency.pdf",         # File name
      width = 11.69, height = 8.27,   # Width and height in inches
      paper = "A4r")          # Paper size
  barplot(table(unlist(c(WT_clonesize))), col='green',
          main='WT clones')

  barplot(prop.table(table(unlist(c(WT_clonesize))))*100,
          col='green', main='WT clones_percentage')
  # Closing the graphical device
  dev.off()
}

#----------------------------------------------------------------------
# FUNCTION TO DRAW CIRCLE AROUND CLONES
#----------------------------------------------------------------------

plot_Centroid_Disp <- function(data, color){
  if (nrow(data) > 0){
    for (i in 1:nrow(data)){
      draw.circle(data$x_mean[i],data$y_mean[i],data$CloneSize[i]*4,
                  border=color,
                  col=color, lwd=1)
      
    }
  }
}

#----------------------------------------------------------------------
#----------------------------------------------------------------------

cen_disp_tot <- c()
cell_type_tot <- c()

for (j in species[,1]) {
  
  cen_disp_tot <- rbind(cen_disp_tot,eval(parse(text = paste0('Cen_disp_',j))))
  cell_type_tot <- rbind(cell_type_tot,eval(parse(text = paste0('cell_type',j))))
  
  }

write.csv(cen_disp_tot,file="Full_Clonal_Data.csv")
write.csv(cell_type_tot,file="Full_Cell_Data.csv")
write.csv(counts, file="CloneSizeCounts.csv")

#----------------------------------------------------------------------
#----------------------------------------------------------------------

# MULTICOLOR CELLS INFO

tot_nucleus_mul <- aggregate(. ~ CloneType, data=counts, FUN=sum)
colnames(tot_nucleus_mul)[1] <- c("Species_Name")
tot_mul <- merge(tot_nucleus_mul, species, how="inner", on="Species_Name")
tot_mul <- tot_mul[,-3:-5]

species_mul <- species
uni <- select.list(species_mul[,2], preselect = NULL, multiple = TRUE, title = "Select the monocolor species", graphics = TRUE)
species_mul <- species_mul[! species_mul$Species_Name %in% uni,]
bi <- select.list(species_mul[,2], preselect = NULL, multiple = TRUE, title = "Select the bicolor species", graphics = TRUE)
species_mul <- species_mul[! species_mul$Species_Name %in% bi,]
mul <- select.list(species_mul[,2], preselect = NULL, multiple = TRUE, title = "Select the multicolor species", graphics = TRUE)
 species_mul <- species_mul[! species_mul$Species_Name %in% mul,]

bn <- askYesNo("Are there any species you do not want to add for normalization?", default = TRUE, 
          prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
if (bn == TRUE){
  
  normal <- select.list(species_mul[,2], preselect = NULL, multiple = TRUE, title = "Are there any species you do not want to add for normalization?", graphics = TRUE)
  normal <- tot_mul[! tot_mul$Species_Name %in% mul,]
  
} else {
  
  normal <- tot_mul
  
}



uni <- tot_mul[tot_mul$Species_Name %in% uni,]
bi <- tot_mul[tot_mul$Species_Name %in% bi,]
mul <- tot_mul[tot_mul$Species_Name %in% mul,]


uni_Mut <- uni[uni$Colour=="red",]
uni_WT <- uni[uni$Colour=="green",]
bi_Mut <- bi[bi$Colour=="red",]
bi_WT <- bi[bi$Colour=="green",]
mul_Mut <- mul[mul$Colour=="red",]
mul_WT <- mul[mul$Colour=="green",]
normal_Mut <- normal[normal$Colour=="red",]
normal_WT <- normal[normal$Colour=="green",]

uni_Mut <- sum(uni_Mut$CloneSize)
uni_WT <- sum(uni_WT$CloneSize)
bi_Mut <- sum(bi_Mut$CloneSize)
bi_WT <- sum(bi_WT$CloneSize)
mul_Mut <- sum(mul_Mut$CloneSize)
mul_WT <- sum(mul_WT$CloneSize)
normal_Mut <- sum(normal_Mut$CloneSize)
normal_WT <- sum(normal_WT$CloneSize)

normal_uni_Mut <- uni_Mut/normal_Mut
normal_uni_WT <- uni_WT/normal_WT
normal_bi_Mut <- bi_Mut/normal_Mut
normal_bi_WT <- bi_WT/normal_WT
normal_mul_Mut <- mul_Mut/normal_Mut
normal_mul_WT <- mul_WT/normal_WT

nucrat_Mut <- (bi_Mut+mul_Mut)
nucrat_WT <- (bi_WT+mul_WT)

normal_nucrat_Mut <- (normal_bi_Mut+normal_mul_Mut)
normal_nucrat_WT <- (normal_bi_WT+normal_mul_WT)


multinuc <- data.frame(Number = c(uni_Mut,bi_Mut,mul_Mut,uni_WT,bi_WT,mul_WT),
                       Type = rep(c("Mutant","WT"), each = 3), Group = c("A - Monocolor","B - Bicolor","C - Multicolor"))

plot1 <- ggplot(multinuc, aes(x = Type, y = Number, fill = Group))+
  geom_bar(stat = "identity", position = "dodge")+
  ggtitle("Absolute")

nucrat <- data.frame(Number = c( uni_Mut, nucrat_Mut, uni_WT, nucrat_WT), Type = rep(c("Mutant","WT"), each = 2),
                     Group = c("A - Monocolor", "B - Multicolor"))

plot2 <- ggplot(nucrat, aes(x = Type, y = Number, fill = Group ))+
       geom_bar(stat = "identity", position = "dodge")+
  ggtitle("Absolute")

normal_multinuc <- data.frame(Number = c(normal_uni_Mut,normal_bi_Mut,normal_mul_Mut,normal_uni_WT,normal_bi_WT,normal_mul_WT),
                       Type = rep(c("Mutant","WT"), each = 3), Group = c("A - Monocolor","B - Bicolor","C - Multicolor"))

plot3 <- ggplot(normal_multinuc, aes(x = Type, y = Number, fill = Group))+
  geom_bar(stat = "identity", position = "dodge")+
  ggtitle("Normalizated")

normal_nucrat <- data.frame(Number = c( normal_uni_Mut, normal_nucrat_Mut, normal_uni_WT, normal_nucrat_WT), Type = rep(c("Mutant","WT"), each = 2),
                            Group = c("A - Monocolor", "B - Multicolor"))

plot4 <- ggplot(normal_nucrat, aes(x = Type, y = Number, fill = Group ))+
  geom_bar(stat = "identity", position = "dodge")+
  ggtitle("Normalizated")

pdf("Multicolor Classification.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4")          # Paper size)

print(plot1)
print(plot3)
print(plot2)
print(plot4)

# Closing the graphical device
dev.off()

multicolor_cell_info <- data.frame(list(cbind(tot_mul,normal_Mut,normal_WT, uni_Mut, uni_WT, bi_Mut, bi_WT, mul_Mut, mul_WT)))
multicolor_cell_info <- multicolor_cell_info[,-3]
colnames(multicolor_cell_info) <- c("Species_Name", "CloneSize", "Total Mutant Nucleus (rep)", "Total WT Nucleus (rep)", "Monocolor Mutant Nucleus (rep)", "Monocolor WT Nucleus (rep)", "Bicolor Mutant Nucleus (rep)", "Bicolor WT Nucleus (rep)", "Multicolor Mutant Nucleus (rep)", "Multicolor WT Nucleus (rep)")

write.csv(multicolor_cell_info,file="MULTICOLOR CELLS INFO.csv")

#----------------------------------------------------------------------
# DRAW CELLS AND CLONE DELIMITING CIRCLE
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# FUNCTION TO REPRESENT INDIVIDUAL CELLS WITH TINY CIRCLES
#----------------------------------------------------------------------

plot_Cell <- function(data, color){
  if (nrow(data) > 0){
    for (i in 1:nrow(data)){
      draw.circle(data$X_microns[i],data$Y_microns[i],1,
                  border=color,
                  col=color, lwd=1)
      
    }
  }
}

#----------------------------------------------------------------------
# PDF SCATTER PLOT
#----------------------------------------------------------------------

tot_nucleus <- aggregate(. ~ CloneType, data=counts, FUN=sum)
tot_nucleus <- tot_nucleus[,-3:-4]

tot_clones1 <- data.frame(table(counts$CloneType))
tot_clones1[] <- lapply(tot_clones1, as.character)
colnames(tot_clones1) <- c("CloneType","Freq")

tot <- merge(tot_nucleus, tot_clones1, how="inner", on="CloneType")

tot <-  cbind(tot, paste0(tot$CloneType,"   Clone Numbers: ",tot$Freq, "   Cell Numbers: ", tot$CloneSize))
colnames(tot)[4] <- c("Selection")
del3 <- select.list(tot[,4], preselect = NULL, multiple = TRUE, title = "Select the species you don't want to represent", graphics = TRUE)

del4 <- tot[tot$Selection %in% del3,]
colnames(del4)[1] <- c("Species_Name")
del4 <- del4[,-2:-4]

species<- species[! species$Species_Name %in% del4,]

mutant_plot <- species[species$Colour=="red",]
WT_plot <- species[species$Colour=="green",]

for (j in mutant_plot[,1]) {

  first_plot <- c(eval(parse(text = paste0('Cen_disp_',i))))
  break

}

hot <- rainbow(nrow(species),alpha = 0.5)

pdf("Cells with Cluster labelling_IMAGE.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4")          # Paper size

plot(first_plot$x_mean, first_plot$y_mean, pch = '.',
     main = 'Clones centroids (Mutant: multicolor and WT: grey)', xlab='x',
     ylab='y', col= "black", xlim=c(0,1500), ylim=c(0,1500))
legend("topright", legend = c(mutant_plot$Species_Name), col = hot, pch = 16, xpd = TRUE, cex = 0.5)

grid(nx = NULL, ny = NULL,
     lty = 3,      # Grid line type
     col = "#e6e6e6", # Grid line color
     lwd = 1)      # Grid line width

c <- 0

for (i in mutant_plot[ , 1]){

  c %+=% 1
  plot_Centroid_Disp(eval(parse(text = paste0('Cen_disp_',i))), hot[c])

}

for (i in WT_plot[ , 1]){

  plot_Centroid_Disp(eval(parse(text = paste0('Cen_disp_',i))), rgb(0.5,0.5,0.5, 0.5))

}

for (i in mutant_plot[ , 1]){

  plot_Cell(eval(parse(text = paste0('cell_type',i))), rgb(55/255,136/255,5/255))

}

for (i in WT_plot[ , 1]){

  plot_Cell(eval(parse(text = paste0('cell_type',i))), rgb(0,0,0))

}

# Closing the graphical device
dev.off()

#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# INTERACTIVE SCATTER PLOT
#----------------------------------------------------------------------

species_plot_clones <- c()
species_plot_cells <- c()

for (j in species[,1]) {
  
  assign(paste0("species_plot_clones",j), data.frame(eval(parse(text = paste0('Cen_disp_',j)))))
  assign(paste0("species_plot_clones",j), cbind(rep(species$Species_Name[species$Num==j], nrow(eval(parse(text = paste0('species_plot_clones',j))))),  eval(parse(text = paste0('species_plot_clones',j)))))
  assign(paste0("species_plot_clones",j), cbind(rep(species$Colour[species$Num==j], nrow(eval(parse(text = paste0('species_plot_clones',j))))),  eval(parse(text = paste0('species_plot_clones',j)))))
  species_plot_clones <- rbind(species_plot_clones,eval(parse(text = paste0('species_plot_clones',j))))
  
  assign(paste0("species_plot_cells",j), data.frame(eval(parse(text = paste0('cell_type',j)))))
  assign(paste0("species_plot_cells",j), cbind(rep(species$Species_Name[species$Num==j], nrow(eval(parse(text = paste0('species_plot_cells',j))))),  eval(parse(text = paste0('species_plot_cells',j)))))
  assign(paste0("species_plot_cells",j), cbind(rep(species$Colour[species$Num==j], nrow(eval(parse(text = paste0('species_plot_cells',j))))),  eval(parse(text = paste0('species_plot_cells',j)))))
  species_plot_cells <- rbind(species_plot_cells,eval(parse(text = paste0('species_plot_cells',j))))
  
}

colnames(species_plot_clones)[2] <- 'CloneType'
colnames(species_plot_cells)[2] <- 'CloneType' 
colnames(species_plot_clones)[1] <- 'Colour'
colnames(species_plot_cells)[1] <- 'Colour' 

species_plot_clones_mutant <- subset(species_plot_clones,species_plot_clones$Colour=="red")
species_plot_cells_mutant <- subset(species_plot_cells,species_plot_cells$Colour=="red")
species_plot_clones_WT <- subset(species_plot_clones,species_plot_clones$Colour=="green")
species_plot_cells_WT <- subset(species_plot_cells,species_plot_cells$Colour=="green")

paleta <- rainbow(nrow(species_plot_clones_mutant),alpha = 0.5)
# imagen <- file.choose()
# imagen <- readPNG(imagen)

plot_ly() %>%
    
    add_trace(data = species_plot_clones_mutant,
              x = species_plot_clones_mutant$x_mean,
              y = species_plot_clones_mutant$y_mean,
              name = "Mutant Clones",
              type = "scatter", mode = "markers",
              marker = list(size = species_plot_clones_mutant$CloneSize*10, 
                            colors = ~species_plot_clones_mutant),
              text = ~paste('Species: ', species_plot_clones_mutant$CloneType)
              ) %>%
    
    add_trace(data = species_plot_clones_WT,
              x = species_plot_clones_WT$x_mean,
              y = species_plot_clones_WT$y_mean,
              name = "WT Clones",
              type = "scatter", mode = "markers",
              marker = list(size = species_plot_clones_WT$CloneSize*10, 
                            color = "rgba(0.5, 0.5, 0.5, .5)",
                            line = list(
                              color = "Lime",
                              width = 1
                            )),
              text = ~paste('Species: ', species_plot_clones_WT$CloneType)
              ) %>%
  
    add_trace(data = species_plot_cells_mutant,
              x = species_plot_cells_mutant$X_microns,
              y = species_plot_cells_mutant$Y_microns,
              name = "Mutant Cells",
              type = "scatter", mode = "markers",
              marker = list(size = 4,
                            color = "Red")
              ) %>%
  
    
    add_trace(data = species_plot_cells_WT,
              x = species_plot_cells_WT$X_microns,
              y = species_plot_cells_WT$Y_microns,
              name = "WT Cells",
              type = "scatter", mode = "markers",
              marker = list(size = 4,
                            color = "Green"))
  


#----------------------------------------------------------------------
#----------------------------------------------------------------------


mutant_plot <- species[species$Colour=="red",]
WT_plot <- species[species$Colour=="green",]


mut_cen_disp_plot <- lapply(paste0('Cen_disp_',mutant_plot[,1]),get)
WT_cen_disp_plot <- lapply(paste0('Cen_disp_',WT_plot[,1]),get)

# Histogram of frequencies of clone size Plotted Mutant cells

mut_clonesize_plot <- lapply(mut_cen_disp_plot, "[", , "CloneSize")

if (nrow(rbind(mut_cen_disp_plot))>0){
  
  pdf("Plotted mutant clones_Frequency.pdf",         # File name
      width = 11.69, height = 8.27,   # Width and height in inches
      paper = "A4r")          # Paper size
  barplot(table(unlist(c(mut_clonesize_plot))), col='red',
          main='Mutant clones')

  barplot(prop.table(table(unlist(c(mut_clonesize_plot))))*100,
          col='red',
          main='Mutant clones_percentage')
  # Closing the graphical device
  dev.off()
}

WT_clonesize_plot <- lapply(WT_cen_disp_plot, "[", , "CloneSize")

# Histogram of frequencies of clone size Plotted WT cells
if (nrow(rbind(WT_cen_disp_plot))>0){
  
  pdf("Plotted WT clones_Frequency.pdf",         # File name
      width = 11.69, height = 8.27,   # Width and height in inches
      paper = "A4r")          # Paper size
  barplot(table(unlist(c(WT_clonesize_plot))), col='green',
          main='WT clones')

  barplot(prop.table(table(unlist(c(WT_clonesize_plot))))*100,
          col='green', main='WT clones_percentage')
  # Closing the graphical device
  dev.off()
}


















