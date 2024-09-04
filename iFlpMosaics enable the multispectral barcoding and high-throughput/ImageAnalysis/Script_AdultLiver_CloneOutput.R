# ANALYSIS OF RAW OUTPUT OBTAINED FROM FIJI USING CELL COUNTER TO LABEL 
# RECOMBINED CELLS ACCORDING TO THEIR COLOR CODE (see color codes in line 39)

# ANALYSIS OF EACH INDIVIDUAL IMAGE

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
data <- read.csv(file)
data$X.micron. <- as.numeric(as.character(data$X.micron.))
data$Y.micron. <- as.numeric(as.character(data$Y.micron.))

# get file name without extension

name <- tools::file_path_sans_ext(file)

color <- unlist(list("red", "red", "red", "red", "red", "red",
                          "red", "red", "red","green", "green","green","green",
                          "red","red", "red", "red")) 

esp <- unlist(list("Mutant-Tomato+nGFP+", "Mutant-Tomato+nCherry+",
                   "Mutant-Tomato+nCer+", "Mutant-Tomato+YFP+nGFP+",
                   "Mutant-Tomato+YFP+nCherry+", "Mutant-Tomato+YFP+nCer+",
                   "Mutant-Tomato+TFP+nGFP+", "Mutant-Tomato+TFP+nCherry+",
                   "Mutant-Tomato+TFP+nCer+", "WT-YFP+nGFP+", "WT-YFP+nCherry+",
                   "WT-YFP+nCer+", "WT-YFP+nGFP+nCherry",
                   "Mutant-Tomato+nGFP+nCherry", "Mutant-Tomato+TFP+YFP+nGFP+",
                   "Mutant-Tomato+YFP+nCherry+nCer+",
                   "Mutant-Tomato+YFP+nGFP+nCherry+"))

polyploidy <- unlist(list(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4))

nuc_plot_colours <- c("#11bc11", "#dd0b6b", "#0000ff", "#11bc11", "#dd0b6b",
                               "#0000ff", "#11bc11", "#dd0b6b", "#0000ff",
                               "#11bc11", "#dd0b6b", "#0000ff", "#f7e00b", 
                               "#dd0b6b", "#11bc11", "#9108c4", "#f7e00b")
                               
memb_plot_colours <- c("#f4a8b3", "#f4a8b3", "#f4a8b3", "#f4b664", "#f4b664",
                               "#f4b664", "#bd72f2", "#bd72f2", "#bd72f2",
                               "#b4ef8e", "#b4ef8e", "#b4ef8e", "#b4ef8e",
                               "#f4a8b3", "#bd72f2", "#f4b664", "#f4b664")
                   
                  
n <- unlist(1:17)  # number possible celltypes

del <- unique(data$Type)

species = data.frame(unlist(n), unlist(esp),unlist(color), unlist(polyploidy))
colnames(species) <- c("Num","Species_Name", "Colour", "Polyploidy")

species <- species[species$Num %in% del, ]

#add column to data with polyploid feature
polyploid <- c()
for (i in 1:nrow(data)){
  for (j in 1:nrow(species)){
    if (data$Type[i] == species$Num[j]){
      polyploid <- c(polyploid, species$Polyploidy[j])
    }
  } 
}

data <- cbind(data, polyploid)

# subset data according to cell_fluo_markers

for (i in species[ , 1]) {

  assign(paste0('data_type',i), subset(data, data$Type==i))

}

#-------------------------------------------------------------------
# FUNCTION TO DETERMINE IF NUCLEI BELONG TO THE SAME CELL (BINUCLEATED)
#-------------------------------------------------------------------
Nuc_ID <- function(data, type, max_dist){
  # function determines whether nuclei with same color code belong
  # to the same cell or not and labels them accordingly (eg. 
  # cells labelled 1.1, 1.2 and 1.3 belong to 3 different cells
  # but have color code type 1)
  
  # VARIABLES
  # -------------------------------------------------------
  # data -> data.frame,
  #         parameters for cells with a given color code
  # type -> integer number,
  #         defines the color code
  # max_dist -> integer number,
  #             maximal distance (in microns) allowed
  #                 between nuclei of the same cell
  
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

#a <- Clone_ID(data_type1,1,40)
#b <- Clone_ID(data_type2,2,20)
for (i in species[ , 1]) {
  
  assign(paste0('Type',i,'_Nuc_ID'), Nuc_ID(eval(parse(text = paste0('data_type',i))), i, 12.5)) # Last number is the max distance for clones
  
}


# join NucIDs in a single vector

NucID <-c()
for (i in species[ , 1]) {
  if (nrow(eval(parse(text = paste0('data_type',i)))) > 0){
    
    NucID <- c(NucID, eval(parse(text = paste0('Type',i,'_Nuc_ID'))))
    
  }
}

# Create new table with relevant data
Comp_data <-c()
for (i in species[ , 1]) {
  if (nrow(eval(parse(text = paste0('data_type',i)))) > 0){
    
    Comp_data <- rbind(Comp_data, eval(parse(text = paste0('data_type',i))))
    
  }
}

Comp_data <- cbind(Comp_data, NucID)
Comp_data <- Comp_data[c('Type', 'X.micron.', 'Y.micron.', 'polyploid',
                         'NucID')]
colnames(Comp_data) <- c("Cell_Type", "X_microns", "Y_microns", 'polyploidy',
                         "NucID")
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
  #               with distance (in microns) to cell with same
  #                 color code
  # -------------------------------------------------------
  Dist_Closest_Cell <- c(rep(0,length(data$NucID)))
  X_Closest_Cell <- c(rep(0,length(data$NucID)))
  Y_Closest_Cell <- c(rep(0,length(data$NucID)))
  
  for (i in 1:length(data$NucID)){
    min_dist = 200
    for (j in 1:length(data$NucID)){
      if (i != j){
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
colnames(Comp_data_new)<- c("Cell_Type", "X_microns", "Y_microns", 'polyploidy',
                            "NucID", "Dist_Closest_Nuc", "X_Closest_Nuc",
                            "Y_Closest_Nuc")

pdf("Distribution of nuclei distance.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

#histogram of distances to closest cell
plot(density(Comp_data_new$Dist_Closest_Nuc), main="Distribution of distances between nuclei with same color code", col = "red")
abline(v=12.5, col="blue")
abline(v=25, col="green")

# Closing the graphical device
dev.off()

Bi <- subset(Comp_data_new, Comp_data_new$Dist_Closest_Nuc < 12.5)
Mono_2_clone <- subset(Comp_data_new, Comp_data_new$Dist_Closest_Nuc > 12.5 & Comp_data_new$Dist_Closest_Nuc < 25)
Mono_Single_cell <- subset(Comp_data_new, Comp_data_new$Dist_Closest_Nuc > 25)

No_binuc_cells <- nrow(Bi)/2  #(number of cells with 2 nuc -> numb_nuc/2)
No_clone_2_mononuc_cells <- nrow(Mono_2_clone)/2  #(number of clones with 2 mononuc cell)
No_single_mononuc_cells <- nrow(Mono_Single_cell)

Per_No_binuc_cells <- No_binuc_cells/(No_binuc_cells+No_clone_2_mononuc_cells+No_single_mononuc_cells)*100
Per_No_clone_2_mononuc_cells <- No_clone_2_mononuc_cells/(No_binuc_cells+No_clone_2_mononuc_cells+No_single_mononuc_cells)*100
Per_No_single_mononuc_cells <- No_single_mononuc_cells/(No_binuc_cells+No_clone_2_mononuc_cells+No_single_mononuc_cells)*100

cell_polyploidy <- cbind(Per_No_binuc_cells, Per_No_clone_2_mononuc_cells, Per_No_single_mononuc_cells)
colnames(cell_polyploidy) <- c("Binuc", "2 Mononuc clone", "Mononuc")

pdf("Polyploidy frequencies.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r") 

barplot(cell_polyploidy, col=c("blue", "red", "green"), ylim=c(0,60))

# Closing the graphical device
dev.off()

# subset data according to Cell_Type

for (i in species[ , 1]) {
  
  assign(paste0('cell_type',i), subset(Comp_data_new, Comp_data_new$Cell_Type==i))
  
}


#----------------------------------------------------------------------
# FUNCTION TO DETERMINE CENTROID OF CELL
#----------------------------------------------------------------------

Centroid_dispersion <- function(cell_type){
  # VARIABLES
  #--------------------------------------------------------------------
  # cell_type -> data.frame,
  #               with clonal info from 1 type of color combination
  
  if (nrow(cell_type) > 0){
    # order the dataset by CloneID so that these are grouped
    cell_type <- cell_type[order(cell_type$NucID), ]
    ID_size <- as.data.frame(table(cell_type$NucID))
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

# place nucleus number and cell centroid data in cell_type_i dataframes 
for (i in species[ , 1]) {
  No_Nuc <- c()
  Cell_center_X <- c()
  Cell_center_Y <- c()
  for (j in 1:nrow(eval(parse(text = paste0('cell_type',i))))){
    value <- eval(parse(text = paste0('cell_type',i)))$NucID[j]
    x <- eval(parse(text = paste0('Cen_disp_',i)))
    ind <- match(value, x$Clone_ID)
    No_Nuc <- c(No_Nuc, x$CloneSize[ind])
    Cell_center_X <- c(Cell_center_X, x$x_mean[ind])
    Cell_center_Y <- c(Cell_center_Y, x$y_mean[ind])
  }
  new_data <- eval(parse(text = paste0('cell_type',i)))
  new_data <- cbind(new_data, No_Nuc, Cell_center_X, Cell_center_Y)
  assign(paste0('cell_type',i), new_data)
}



#-------------------------------------------------------------------
# FUNCTION TO DETERMINE CELLULAR CLONALITY
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
          if (i != j) {    #guarantee that comparing different cells
            if (data$NucID[i] == data$NucID[j]){
              vector[j] = 'NA'
            } else {
              x1 <- data$Cell_center_X[i]
              y1 <- data$Cell_center_Y[i]
              x2 <- data$Cell_center_X[j]
              y2 <- data$Cell_center_Y[j]
              dist <- sqrt((x1-x2)^2+(y1-y2)^2)
              if (dist < max_dist){
                vector[j] = ID
                count = count + 1
              }
            }
          }
        }
        while (SizeClone < count){
          initialSizeClone = count
          for (k in 1:nrow(data)){
            for (l in 1:nrow(data)){
              if ((vector[l] ==0) & (vector[k] == ID)){
                x1 <- data$Cell_center_X[k]
                y1 <- data$Cell_center_Y[k]
                x2 <- data$Cell_center_X[l]
                y2 <- data$Cell_center_Y[l]
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

for (i in species[ , 1]) {
  
  assign(paste0('Type',i,'_clone_ID'), Clone_ID(eval(parse(text = paste0('cell_type',i))), i, 25)) # Last number is the max distance for clones
  
}


# join CloneIDs into cell_type data
for (i in species[ , 1]) {
  assign('CloneID', eval(parse(text = paste0('Type',i,'_clone_ID'))))
  new_data <- eval(parse(text = paste0('cell_type',i)))
  new_data <- cbind(new_data, CloneID)
  assign(paste0('cell_type',i), new_data)
}


#----------------------------------------------------------------------
# FUNCTION TO DETERMINE CENTROID OF CLONE
#----------------------------------------------------------------------

Clone_Centroid <- function(cell_type){
  # VARIABLES
  #--------------------------------------------------------------------
  # cell_type -> data.frame,
  #               with clonal info from 1 type of color combination
  
  if (nrow(cell_type) > 0){
    # order the dataset by CloneID so that these are grouped
    cell_type <- cell_type[order(cell_type$CloneID), ]
    ID_size <- as.data.frame(table(cell_type$CloneID))
    colnames(ID_size) <- c("Clone_ID", "CloneSize")
    # calculate centroid coordinates
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
          x_poly <-c(x_poly, cell_type$Cell_center_X[j])
          y_poly <-c(y_poly, cell_type$Cell_center_Y[j])
        }
        start = end + 1
        # coordinates of centroid (average of x-coord and y-coord)
        x_mean <- c(x_mean, mean(x_poly))
        y_mean <- c(y_mean, mean(y_poly))
        
        
      }else if(ID_size$CloneSize[i] == 1){
        x_mean <- c(x_mean, cell_type$Cell_center_X[start])
        y_mean <- c(y_mean, cell_type$Cell_center_Y[start])
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
  
  assign(paste0('CloneCent_',i), Clone_Centroid(eval(parse(text = paste0('cell_type',i)))))
  
}

# join Clone centroids into cell_type data
for (i in species[ , 1]) {
  Clone_Cent_X <- c()
  Clone_Cent_Y <- c()
  CloneSize <- c()
  x <- eval(parse(text = paste0('cell_type',i)))
  y <- eval(parse(text = paste0('CloneCent_',i))) 
  for (j in 1:nrow(x)){
    if (x$CloneID[j] != 'NA'){
      ind <- match(x$CloneID[j], y$Clone_ID)
      Clone_Cent_X <- c(Clone_Cent_X, y$x_mean[ind])
      Clone_Cent_Y <- c(Clone_Cent_Y, y$y_mean[ind])
      CloneSize <- c(CloneSize, y$CloneSize[ind])
    } else {
      Clone_Cent_X <- c(Clone_Cent_X, 'NA')
      Clone_Cent_Y <- c(Clone_Cent_Y, 'NA')
      CloneSize <- c(CloneSize, 'NA')
    }
      
  }
  x <- cbind(x, Clone_Cent_X, Clone_Cent_Y, CloneSize)
  assign(paste0('cell_type',i), x)
}



# collecting final numbers
for (i in species[ , 1]) {
  # data to work with
  x <- eval(parse(text = paste0('cell_type',i)))
  x$CloneID <- as.numeric(x$CloneID)
  x$Clone_Cent_X <- as.numeric(x$Clone_Cent_X)
  x$Clone_Cent_Y <- as.numeric(x$Clone_Cent_Y)
  x$No_Nuc <- as.numeric(x$No_Nuc)
  x$CloneSize <- as.numeric(x$CloneSize)
  
  Clonal_dat <- data.frame()

  # subset for total single cells
  Single <- subset(x, x$CloneSize == 1)
  Single_mono_2n <- nrow(subset(Single, Single$No_Nuc == 1 & Single$polyploidy == 2))
  Single_mono_4n <- nrow(subset(Single, Single$No_Nuc == 1 & Single$polyploidy == 4))
  Single_bi_2n_x_2n <- nrow(subset(Single, Single$No_Nuc == 2 & Single$polyploidy == 2))
  Single_bi_4n_x_4n <- nrow(subset(Single, Single$No_Nuc == 2 & Single$polyploidy == 4))
  
  Clonal_dat <- cbind(rep('Single', 4),
                      rbind(Single_mono_2n, Single_mono_4n, Single_bi_2n_x_2n,
                            Single_bi_4n_x_4n))

  # subset for total clonal cells with 2n nuclei
  Clonal_2n <- subset(x, x$CloneSize > 1 & x$polyploidy == 2)
  max_clone <- 0
  
  if (nrow(Clonal_2n)> 0){
    max_clone <- max(Clonal_2n$CloneSize)
    
    for (m in 2:max_clone){
      clone_type <- c()
      cent_coord <- c()
      Clonal_mono_2n = 0
      Clonal_bi_2n = 0
      Clonal_mixed_2n = 0
      
      Clonal_pac <- c()
      
      clone_type <- subset(Clonal_2n, Clonal_2n$CloneSize == m)
      if (nrow(clone_type) > 0){
        cent_coord <- unique(clone_type$Clone_Cent_X)
        for (n in 1:length(cent_coord)){
          sub <- subset(Clonal_2n, Clonal_2n$Clone_Cent_X == cent_coord[n])
          if (sum(sub$No_Nuc) == nrow(sub)){
            Clonal_mono_2n = Clonal_mono_2n + 1
          }else if (sum(sub$No_Nuc) == nrow(sub)*2){
            Clonal_bi_2n = Clonal_bi_2n + 1
          }else {
            Clonal_mixed_2n = Clonal_mixed_2n + 1
          }
        }
        
        Clonal_pac <- cbind(rep(paste0('Clone_with_', m, '_cells'), 3),
                            rbind(Clonal_mono_2n, Clonal_bi_2n, Clonal_mixed_2n))
        Clonal_dat <- rbind(Clonal_dat, Clonal_pac)
      }
      
      
    }
  }
  
  
  # subset for total clonal cells with 4n nuclei
  Clonal_4n <- subset(x, x$CloneSize > 1 & x$polyploidy == 4)
  max_clone <- 0
  
  if (nrow(Clonal_4n)> 0){
    max_clone <- max(Clonal_4n$CloneSize)
    
    for (m in 2:max_clone){
      clone_type <- c()
      cent_coord <- c()
      Clonal_mono_4n = 0
      Clonal_bi_4n = 0
      Clonal_mixed_4n = 0
      
      Clonal_pac <- c()
      
      clone_type <- subset(Clonal_4n, Clonal_4n$CloneSize == m)
      if (nrow(clone_type) > 0){
        cent_coord <- unique(clone_type$Clone_Cent_X)
        for (n in 1:length(cent_coord)){
          sub <- subset(Clonal_4n, Clonal_4n$Clone_Cent_X == cent_coord[n])
          if (sum(sub$No_Nuc) == nrow(sub)){
            Clonal_mono_4n = Clonal_mono_4n + 1
          }else if (sum(sub$No_Nuc) == nrow(sub)*2){
            Clonal_bi_4n = Clonal_bi_4n + 1
          }else {
            Clonal_mixed_4n = Clonal_mixed_4n + 1
          }
        }
        
        Clonal_pac <- cbind(rep(paste0('Clone_with_', m, '_cells'), 3),
                            rbind(Clonal_mono_4n, Clonal_bi_4n, Clonal_mixed_4n))
        Clonal_dat <- rbind(Clonal_dat, Clonal_pac)
      }
      
      
    }
  }
  
  Descrip <- row.names(Clonal_dat)
  Clonal_dat <- cbind(Descrip, Clonal_dat)
  assign(paste0('clonal_data_cell_type', i), as.data.frame(Clonal_dat))
}  
  

Colour_com <- rep(species$Species_Name[i], 4)
Geno <- rep(species$Colour[i], 4)

# Save full clonal data
FullClonalData <- c()
for (i in species[ , 1]) {
  # data to work with
  x <- eval(parse(text = paste0('clonal_data_cell_type',i)))
  z <- match(i, species$Num)
  Colour_com <- rep(species$Species_Name[z], nrow(x))
  Genot <- rep(species$Colour[z], nrow(x))
  x <- cbind(Colour_com, Genot, x)
  FullClonalData <- rbind(FullClonalData, x)
}
write.csv(FullClonalData,file="FullClonalData.csv")


#summary of clonal data
mutant <- subset(FullClonalData, FullClonalData$Genot == 'red')
comb <- paste0(mutant$Descrip, '+', mutant$V2)
mutant <- cbind(mutant, comb)
types <- unique(mutant$comb)
abs_no <- c()
for (i in 1:length(types)){
  a <- sum(subset(as.numeric(mutant$V3), mutant$comb == types[i]))
  abs_no <- c(abs_no, a)
}
mutant_abs <- cbind(types, abs_no)
mutant_abs <- cbind(rep('Mut', nrow(mutant_abs)), mutant_abs)




wt <- subset(FullClonalData, FullClonalData$Genot == 'green')
comb_w <- paste0(wt$Descrip, '+', wt$V2)
wt <- cbind(wt, comb_w)
types <- unique(wt$comb)
abs_no <- c()
for (i in 1:length(types)){
  a <- sum(subset(as.numeric(wt$V3), wt$comb == types[i]))
  abs_no <- c(abs_no, a)
}
wt_abs <- cbind(types, abs_no)
wt_abs <- cbind(rep('Wt', nrow(wt_abs)), wt_abs)


summary_data <- as.data.frame(rbind(wt_abs, mutant_abs))
write.csv(summary_data,file="SummaryClonalData.csv")


# Save full clonal data
FullCellData <- c()
for (i in species[ , 1]) {
  # data to work with
  x <- eval(parse(text = paste0('cell_type',i)))
  FullCellData <- rbind(FullCellData, x)
}
write.csv(FullCellData,file="FullCellData.csv")




# PLOTTING THE DATA


#----------------------------------------------------------------------
# FUNCTION TO DRAW CIRCLE TO REPRESENT CELLS
#----------------------------------------------------------------------
plot_Cells <- function(data, color){
  if (nrow(data) > 0){
    for (i in 1:nrow(data)){
      draw.circle(data$Cell_center_X[i],data$Cell_center_Y[i], 4,
                  border=color, col=color, lwd=1)
    }
  }
}
#----------------------------------------------------------------------
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# FUNCTION TO REPRESENT INDIVIDUAL NUCLEI WITH TINY CIRCLES
#----------------------------------------------------------------------
plot_Nuclei <- function(data, color){
  if (nrow(data) > 0){
    for (i in 1:nrow(data)){
      draw.circle(data$X_microns[i],data$Y_microns[i], radius=data$polyploidy/2,
                  border=color,
                  col=color, lwd=1)
      
    }
  }
}
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# FUNCTION TO DRAW CIRCLE TO REPRESENT CLONES
#----------------------------------------------------------------------
plot_Clones <- function(data){
  if (nrow(data) > 0){
    for (i in 1:nrow(data)){
      if (data$CloneSize[i] > 1)
        draw.circle(as.numeric(data$Clone_Cent_X)[i],
                    as.numeric(data$Clone_Cent_Y)[i],
                    radius = as.numeric(data$CloneSize)[i]*7, border="#b2beb5",
                    lty = 1, lwd=0.5)
    }
  }
}
#----------------------------------------------------------------------
#----------------------------------------------------------------------


pdf("Cells with Cluster labelling_IMAGE.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4")          # Paper size

plot(cell_type1$Cell_center_X, cell_type1$Cell_center_Y, pch = '.',
     main = 'Nuclei_cells_clones centroids', xlab ="x",
     ylab='y', col= "white", xlim=c(0,1500), ylim=c(0,1500))
grid(nx = NULL, ny = NULL,
     lty = 3,      # Grid line type
     col = "#e6e6e6", # Grid line color
     lwd = 1)      # Grid line width

# plot cells
for (j in species[ , 1]) {
  plot_Cells(eval(parse(text = paste0('cell_type',j))), memb_plot_colours[j])
}

# plot nuclei
for (j in species[ , 1]) {
  plot_Nuclei(eval(parse(text = paste0('cell_type',j))), nuc_plot_colours[j])
}

# plot clones
for (j in species[ , 1]) {
  plot_Clones(eval(parse(text = paste0('cell_type',j))))
}

# Closing the graphical device
dev.off()










