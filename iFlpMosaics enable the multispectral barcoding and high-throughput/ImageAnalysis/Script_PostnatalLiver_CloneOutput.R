# ANALYSIS OF RAW OUTPUT OBTAINED FROM FIJI USING CELL COUNTER TO LABEL 
# RECOMBINED CELLS ACCORDING TO THEIR COLOR CODE (color codes in lines 19-23)
# Although script was initially though for 9 different colour combinations,
# only used 5 colours for this specific experiment.

# opening and reading the data
path <- choose.dir()
setwd(path)
file <- file.choose()
data <- read.csv(file, sep = ";")

# get file name without extension
name <- tools::file_path_sans_ext(file)


# subset data according to cell_fluo_markers
data_type1 <- subset(data, data$Type=="1")
data_type2 <- subset(data, data$Type=="2")
data_type3 <- subset(data, data$Type=="3")  # MbTomato+,H2B-GFP+
data_type4 <- subset(data, data$Type=="4")  # MbTomato+,H2B-Cer+
data_type5 <- subset(data, data$Type=="5")  # MbYFP+,H2B-GFP+
data_type6 <- subset(data, data$Type=="6")  # MbYFP+,H2B-Cer+
data_type7 <- subset(data, data$Type=="7")  # MbYFP+,H2B-Cherry+
data_type8 <- subset(data, data$Type=="8")
data_type9 <- subset(data, data$Type=="9")


library(ggplot2)
library(ggpubr)


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
  #         defines the celltype code
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
              if ((vector[l] == 0) & (vector[k] == ID)){
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
Type1_clone_ID <- Clone_ID(data_type1, 1, 40)
Type2_clone_ID <- Clone_ID(data_type2, 2, 40)
Type3_clone_ID <- Clone_ID(data_type3, 3, 40)
Type4_clone_ID <- Clone_ID(data_type4, 4, 40)
Type5_clone_ID <- Clone_ID(data_type5, 5, 40)
Type6_clone_ID <- Clone_ID(data_type6, 6, 40)
Type7_clone_ID <- Clone_ID(data_type7, 7, 40)
Type8_clone_ID <- Clone_ID(data_type8, 8, 40)
Type9_clone_ID <- Clone_ID(data_type9, 9, 40)

# join CloneIDs in a single vector
CloneID <-c()
if (nrow(data_type1) > 0){
  CloneID <- c(CloneID, Type1_clone_ID)
}
if (nrow(data_type2) > 0){
  CloneID <- c(CloneID, Type2_clone_ID)
}
if (nrow(data_type3) > 0){
  CloneID <- c(CloneID, Type3_clone_ID)
}
if (nrow(data_type4) > 0){
  CloneID <- c(CloneID, Type4_clone_ID)
}
if (nrow(data_type5) > 0){
  CloneID <- c(CloneID, Type5_clone_ID)
}
if (nrow(data_type6) > 0){
  CloneID <- c(CloneID, Type6_clone_ID)
}
if (nrow(data_type7) > 0){
  CloneID <- c(CloneID, Type7_clone_ID)
}
if (nrow(data_type8) > 0){
  CloneID <- c(CloneID, Type8_clone_ID)
}
if (nrow(data_type9) > 0){
  CloneID <- c(CloneID, Type9_clone_ID)
}

# Create new table with relevant data 
Comp_data <- cbind(data$Type, data$X.micron., data$Y.micron., CloneID)
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
cell_type1 <- subset(Comp_data_new, Comp_data_new$Cell_Type=="1")
cell_type2 <- subset(Comp_data_new, Comp_data_new$Cell_Type=="2")
cell_type3 <- subset(Comp_data_new, Comp_data_new$Cell_Type=="3")
cell_type4 <- subset(Comp_data_new, Comp_data_new$Cell_Type=="4")
cell_type5 <- subset(Comp_data_new, Comp_data_new$Cell_Type=="5")
cell_type6 <- subset(Comp_data_new, Comp_data_new$Cell_Type=="6")
cell_type7 <- subset(Comp_data_new, Comp_data_new$Cell_Type=="7")
cell_type8 <- subset(Comp_data_new, Comp_data_new$Cell_Type=="8")
cell_type9 <- subset(Comp_data_new, Comp_data_new$Cell_Type=="9")

library(plotrix)




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

Cen_disp_1 <- Centroid_dispersion(cell_type1)
Cen_disp_2 <- Centroid_dispersion(cell_type2)
Cen_disp_3 <- Centroid_dispersion(cell_type3)
Cen_disp_4 <- Centroid_dispersion(cell_type4)
Cen_disp_5 <- Centroid_dispersion(cell_type5)
Cen_disp_6 <- Centroid_dispersion(cell_type6)
Cen_disp_7 <- Centroid_dispersion(cell_type7)
Cen_disp_8 <- Centroid_dispersion(cell_type8)
Cen_disp_9 <- Centroid_dispersion(cell_type9)

# Histogram of frequencies of clone size
if (nrow(Cen_disp_1)>0){
    barplot(table(Cen_disp_1$CloneSize))
}

# Histogram of frequencies of clone size
if (nrow(Cen_disp_2)>0){
  barplot(table(Cen_disp_2$CloneSize))
}

# Histogram of frequencies of clone size
pdf("Clone size frequencies_Mutant_nGFP.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
if (nrow(Cen_disp_3)>0){
  barplot(table(Cen_disp_3$CloneSize), col='red', main = 'KO-nGFP', ylim=c(0,30))
}
# Closing the graphical device
dev.off()

counts <- as.data.frame(table(Cen_disp_3$CloneSize))
counts_p <- as.data.frame(prop.table(table(Cen_disp_3$CloneSize))*100)
counts <- cbind(rep('KO_nGFP', nrow(counts)), counts, counts_p)
colnames(counts) <- c('CloneType','CloneSize','AbsFrequency', 'CloneSize',
                      'RelatFrequency')

# Histogram of frequencies of clone size
pdf("Clone size frequencies_Mutant_nCer.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
if (nrow(Cen_disp_4)>0){
  barplot(table(Cen_disp_4$CloneSize), col='red', main = 'KO-nCer', ylim=c(0,30))
}
# Closing the graphical device
dev.off()

counts1 <- as.data.frame(table(Cen_disp_4$CloneSize))
counts_p <- as.data.frame(prop.table(table(Cen_disp_4$CloneSize))*100)
counts1 <- cbind(rep('KO_nCer', nrow(counts1)), counts1, counts_p)
colnames(counts1) <- c('CloneType','CloneSize','AbsFrequency', 'CloneSize',
                       'RelatFrequency')
counts <- rbind(counts, counts1)

# Histogram of frequencies of clone size
pdf("Clone size frequencies_Control_nGFP.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
if (nrow(Cen_disp_5)>0){
  barplot(table(Cen_disp_5$CloneSize), col='green', main = 'WT-nGFP', ylim=c(0,30))
}
# Closing the graphical device
dev.off()

counts1 <- as.data.frame(table(Cen_disp_5$CloneSize))
counts_p <- as.data.frame(prop.table(table(Cen_disp_5$CloneSize))*100)
counts1 <- cbind(rep('Control_nGFP', nrow(counts1)), counts1, counts_p)
colnames(counts1) <- c('CloneType','CloneSize','AbsFrequency', 'CloneSize',
                       'RelatFrequency')
counts <- rbind(counts, counts1)

# Histogram of frequencies of clone size
pdf("Clone size frequencies_Control_nCer.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
if (nrow(Cen_disp_6)>0){
  barplot(table(Cen_disp_6$CloneSize), col='green', main = 'WT-nCer', ylim=c(0,30))
}
# Closing the graphical device
dev.off()

counts1 <- as.data.frame(table(Cen_disp_6$CloneSize))
counts_p <- as.data.frame(prop.table(table(Cen_disp_6$CloneSize))*100)
counts1 <- cbind(rep('Control_nCer', nrow(counts1)), counts1, counts_p)
colnames(counts1) <- c('CloneType','CloneSize','AbsFrequency', 'CloneSize',
                       'RelatFrequency')
counts <- rbind(counts, counts1)


# Histogram of frequencies of clone size
pdf("Clone size frequencies_Control_nCherry.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
if (nrow(Cen_disp_7)>0){
  barplot(table(Cen_disp_7$CloneSize), col='green', main = 'WT-nCherry', ylim=c(0,30))
}
# Closing the graphical device
dev.off()

counts1 <- as.data.frame(table(Cen_disp_7$CloneSize))
counts_p <- as.data.frame(prop.table(table(Cen_disp_7$CloneSize))*100)
counts1 <- cbind(rep('Control_nCherry', nrow(counts1)), counts1, counts_p)
colnames(counts1) <- c('CloneType','CloneSize','AbsFrequency', 'CloneSize',
                       'RelatFrequency')
counts <- rbind(counts, counts1)

# Histogram of frequencies of clone size
if (nrow(Cen_disp_8)>0){
  barplot(table(Cen_disp_8$CloneSize))
}

# Histogram of frequencies of clone size
if (nrow(Cen_disp_9)>0){
  barplot(table(Cen_disp_9$CloneSize))
}

# Histogram of frequencies of clone size ALL Mutant cells
if (nrow(rbind(Cen_disp_3, Cen_disp_4))>0){
  
  pdf("Total mutant clones_AbsFrequency.pdf",         # File name
      width = 11.69, height = 8.27,   # Width and height in inches
      paper = "A4r")          # Paper size
  barplot(table(c(Cen_disp_3$CloneSize, Cen_disp_4$CloneSize)), col='red',
          main='Mutant clones')
  # Closing the graphical device
  dev.off()
  
  pdf("Total mutant clones_RelativeFrequency.pdf",         # File name
      width = 11.69, height = 8.27,   # Width and height in inches
      paper = "A4r")          # Paper size
  barplot(prop.table(table(c(Cen_disp_3$CloneSize, Cen_disp_4$CloneSize)))*100,
          col='red',
          main='Mutant clones_percentage')
  # Closing the graphical device
  dev.off()
}

# Histogram of frequencies of clone size WT cells
if (nrow(rbind(Cen_disp_5, Cen_disp_6, Cen_disp_7))>0){
  
  pdf("Total WT clones_AbsFrequency.pdf",         # File name
      width = 11.69, height = 8.27,   # Width and height in inches
      paper = "A4r")          # Paper size
  barplot(table(c(Cen_disp_5$CloneSize, Cen_disp_6$CloneSize, Cen_disp_7$CloneSize)), col='green',
          main='WT clones')
  # Closing the graphical device
  dev.off()
  
  pdf("Total WT clones_RelativeFrequency.pdf",         # File name
      width = 11.69, height = 8.27,   # Width and height in inches
      paper = "A4r")          # Paper size
  barplot(prop.table(table(c(Cen_disp_5$CloneSize, Cen_disp_6$CloneSize, Cen_disp_7$CloneSize)))*100,
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


Full_Clonal_Data <- rbind(Cen_disp_1, Cen_disp_2, Cen_disp_3, Cen_disp_4, Cen_disp_5, Cen_disp_6, Cen_disp_7, Cen_disp_8, Cen_disp_9)
write.csv(Full_Clonal_Data,file="Full_Clonal_Data.csv")

Cell_clone_ID <- rbind(cell_type1, cell_type2, cell_type3, cell_type4, cell_type5, cell_type6, cell_type7, cell_type8, cell_type9)
Cell_clone_ID <- cbind(Cell_clone_ID, data)
write.csv(Cell_clone_ID,file="Full_Cell_Data.csv")

write.csv(counts, file="CloneSizeCount.csv")




#----------------------------------------------------------------------
# DRAW CELLS AND CLONE DELIMITING CIRCLE
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# FUNCTION TO DRAW TINY CIRCLES FOR EACH CELL
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
#----------------------------------------------------------------------


pdf("Cells with Cluster labelling.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4")          # Paper size

plot(Cen_disp_1$x_mean, Cen_disp_1$y_mean, pch = '.', 
     main = 'Clones centroids (KO: red and WT: green)', xlab='x',
     ylab='y', col= "black", xlim=c(0,4000), ylim=c(0,4000))
grid(nx = NULL, ny = NULL,
     lty = 3,      # Grid line type
     col = "#e6e6e6", # Grid line color
     lwd = 1)      # Grid line width


plot_Centroid_Disp(Cen_disp_3, rgb(1,0,0, 0.5))
plot_Centroid_Disp(Cen_disp_4, rgb(1,0,0, 0.5))
plot_Centroid_Disp(Cen_disp_5, rgb(0,1,0, 0.5))
plot_Centroid_Disp(Cen_disp_6, rgb(0,1,0, 0.5))
plot_Centroid_Disp(Cen_disp_7, rgb(0,1,0, 0.5))

plot_Cell(cell_type3, rgb(55/255,136/255,5/255))
plot_Cell(cell_type4, rgb(0,68/255,1))
plot_Cell(cell_type5, rgb(55/255,136/255,5/255))
plot_Cell(cell_type6, rgb(0,68/255,1))
plot_Cell(cell_type7, rgb(0.1,0.1,0.1))


# Closing the graphical device
dev.off()




