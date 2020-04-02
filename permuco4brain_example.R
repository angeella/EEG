#Example clusterMass from permuco4brain using eeguana

library(igraph)
library(permuco)
library(permuco4brain)
library(Matrix)
library(edf)
library(abind)
library(xlsx)
# Set working directory
download.file("https://zenodo.org/record/1169140/files/ERP_by_subject_by_condition_data.zip","eeg.zip")
unzip("eeg.zip")


design <- expand.grid(subject = c(111, 113, 115, 116, 117, 118, 120, 122, 123, 124, 126, 127, 128, 
                                  130, 131, 132, 134, 135, 137, 138, 139, 140, 141, 142, 143, 144, 
                                  145, 146, 147), stimuli = c("AP", "SED"), action = c("Av_App", "Av_Ev"))


edf_filname <- list.files("raw_data")

signal <- list()
for(i in 1:nrow(design)){
  # selecting the experimental contion
  which_id <- grepl(design$subject[i], edf_filname) 
  which_stimuli <- grepl(design$stimuli[i], edf_filname) 
  which_action <- grepl(design$action[i], edf_filname) 
  
  # selecting the neutral VS Not neutral
  which_neutral <- grepl("Rond", edf_filname) 
  which_nonneutral <- grepl("Task", edf_filname)
  
  # dowloading the data
  data_nonneutral <- read.edf(paste("raw_data/", edf_filname[which_id&which_stimuli&which_action&which_nonneutral], sep = ""))
  data_neutral <- read.edf(paste("raw_data/", edf_filname[which_id&which_action&which_neutral], sep = ""))
  
  # Getting the signal
  data_nonneutral <- data_nonneutral$signal[sort(names(data_nonneutral$signal))]
  data_nonneutral <- t(sapply(data_nonneutral, function(x)x$data))
  data_neutral <- data_neutral$signal[sort(names(data_neutral$signal))]
  data_neutral <- t(sapply(data_neutral, function(x)x$data))
  
  # Storing the signal relative to neutral
  signal[[i]] <- data_nonneutral - data_neutral}



# creating the 3D array we the appropriate dimension
signal <- abind(signal, along = 3)
signal <- aperm(signal, c(3, 2, 1))

# Select usefull time windows (800ms begining at 200ms at 512hz)

signal <- signal[,102:512,]


data_sr <- read.csv("https://zenodo.org/record/1169140/files/data_self_report_R_subset_zen.csv",sep=";")
design <- dplyr::left_join(design, data_sr, by = "subject")


# Reshape the design

design$stimuli <- plyr::revalue(design$stimuli, c("AP" = "pa","SED" = "sed"))
design$action <- plyr::revalue(design$action, c("Av_App" = "appr","Av_Ev" = "avoid"))
design$mvpac <- as.numeric(scale(design$MVPA, scale = F))

download.file("https://www.biosemi.com/download/Cap_coords_all.xls","Cap_coords_all.xls",mode="wb")

# import electrode coordinate: the third sheet contains the position of 64 electrodes cap.
coord <-  read.xlsx(file="Cap_coords_all.xls", sheetIndex = 3, header =T,startRow=34)

# Clean coordinate data
coord <- coord[1:64,c(1,5:7)]
colnames(coord) <- c("electrode","x","y","z")

coord$electrode <- plyr::revalue(coord$electrode, c("T7 (T3)" = "T7", 
                                                    "Iz (inion)" = "Iz", "T8 (T4)" = "T8", "Afz" = "AFz"))
coord$electrode <-  as.character(coord$electrode)

distance_matrix <- dist(coord[, 2:4])

adjacency_matrix <- as.matrix(distance_matrix) < 35
diag(adjacency_matrix) <- FALSE

dimnames(adjacency_matrix) = list(coord[,1], coord[,1])

graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
graph <- delete_vertices(graph, V(graph)[!get.vertex.attribute(graph, "name")%in%(coord[,1])])

graph <- set_vertex_attr(graph,"x", value = coord[match(vertex_attr(graph,"name"),coord[,1]),2])
graph <-set_vertex_attr(graph,"y", value = coord[match(vertex_attr(graph,"name"),coord[,1]),3])
graph <-set_vertex_attr(graph,"z", value = coord[match(vertex_attr(graph,"name"),coord[,1]),4])
plot(graph)

np <- 5000
aggr_FUN <- sum
ncores <- 5
contr <- contr.sum
formula <- signal ~ mvpac*stimuli*action + Error(subject/(stimuli*action))
pmat <- Pmat(np = np, n = nrow(design))

model <- permuco4brain::brainperm(formula = formula,
                                    data = design,
                                    graph = graph,
                                    np = np,
                                    method = NULL,
                                    type = "signflip",
                                    test = "fisher",
                                    aggr_FUN = NULL,
                                    threshold = NULL,
                                    multcomp = "clustermass",
                                    effect = 2,
                                    return_distribution = TRUE)

save(model, file = "model_permuco4brain.RData")