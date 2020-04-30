library(edf)
library(abind)

# Set working directory
setwd("C:/Users/Angela Andreella/Documents/Thesis_Doc/Permutation_based_ARI/EEG/EEG/")
download.file("https://zenodo.org/record/1169140/files/ERP_by_subject_by_condition_data.zip","eeg.zip")
unzip("eeg.zip")
head(list.files("raw_data"),n=10)
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
library(igraph)
library(permuco)
library(Matrix)
library(xlsx)
download.file("https://www.biosemi.com/download/Cap_coords_all.xls","Cap_coords_all.xls",mode="wb")

# import electrode coordinate: the third sheet contains the position of 64 electrodes cap.
coord <-  read.xlsx(file="Cap_coords_all.xls", sheetIndex = 3, header =T,startRow=34)

# Clean coordinate data
coord <- coord[1:64,c(1,5:7)]
colnames(coord) <- c("electrode","x","y","z")

coord$electrode <- plyr::revalue(coord$electrode, c("T7 (T3)" = "T7", 
                                                    "Iz (inion)" = "Iz", "T8 (T4)" = "T8", "Afz" = "AFz"))
coord$electrode <-  as.character(coord$electrode)
# Creating adjacency matrix
distance_matrix <- dist(coord[, 2:4])

adjacency_matrix <- as.matrix(distance_matrix) < 35
diag(adjacency_matrix) <- FALSE

dimnames(adjacency_matrix) = list(coord[,1], coord[,1])
graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
graph <- delete_vertices(graph, V(graph)[!get.vertex.attribute(graph, "name")%in%(coord[,1])])

graph <- set_vertex_attr(graph,"x", value = coord[match(vertex_attr(graph,"name"),coord[,1]),2])
graph <-set_vertex_attr(graph,"y", value = coord[match(vertex_attr(graph,"name"),coord[,1]),3])
graph <-set_vertex_attr(graph,"z", value = coord[match(vertex_attr(graph,"name"),coord[,1]),4])
B <- 4000
aggr_FUN <- sum
ncores <- 5
contr <- contr.sum
formula <- signal~ mvpac*stimuli*action + Error(subject/(stimuli*action))
pmat <- Pmat(np = np, n = nrow(design))

model <- brainperm(formula = formula,
                   data = design,
                   graph = graph,
                   np = B,
                   method = NULL,
                   type = "signflip",
                   test = "fisher",
                   aggr_FUN = NULL,
                   threshold = NULL,
                   multcomp = "clustermass",
                   effect = NULL,
                   return_distribution = TRUE)

eff = "`stimuli:action`"
family = "Simes"
delta = 0
ct = c(0,1)
alpha = 0.1

Test <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$uncorrected$distribution")))
dim(Test) <- c(dim(Test)[1], dim(Test)[2]*dim(Test)[3])

test_obs <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$data$statistic")))
TT <- rbind(test_obs,Test)
pvalues <- switch(alternative, 
                  "two.sided" =  matrixStats::colRanks(-abs(TT)) / nrow(TT),
                  "greater" = matrixStats::colRanks(-TT) / nrow(TT),
                  "less" = matrixStats::colRanks(TT) / nrow(TT))


pvalues <-t(pvalues)
#pvalues <- pvalues[,which(model$multiple_comparison[[1]]$clustermass$data$cluster_id!=0)]
pvalues_ord <- rowSortC(pvalues)
lambda <- lambdaOpt(pvalues = pvalues_ord, family = family, ct = ct, alpha = alpha, delta = delta) 
cvOpt = cv(pvalues = pvalues_ord, family = family, alpha = alpha, lambda= lambda, delta = delta)

#clstr_id <- model$multiple_comparison[[1]]$clustermass$data$cluster_id[which(model$multiple_comparison[[1]]$clustermass$data$cluster_id!=0)]
clstr_id <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$data$cluster_id")))

#clusters <- c(1:model$multiple_comparison[[1]]$clustermass$cluster$no)[which(model$multiple_comparison[[1]]$clustermass$cluster$pvalue<=0.1)]
#clusters <- which(model$multiple_comparison[[1]]$clustermass$cluster$pvalue<=0.1)
clusters <- c(1:eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$cluster$no"))))
#hom <-hommel(pvalues[1,])

out=lapply(clusters,function(i){
  ix= which(clstr_id == i)
  
  #cluster_ids=which(ix,arr.ind = TRUE)
  #cluster_ids=cbind(cluster_ids,Stat=StatFun(ix))
  #Error if I put pvalues[,mask] instead of pvalues in SingleStepCT
  #perm <- SingleStepCT(pvalues = pvalues,ct =ct, ix =as.vector(which(ix[mask])), alpha = alpha, shift = shift, family = 'Simes', lambda = lambda)
  #perm <- discoveriesPerm(praw = praw, ix = ix[mask], cvh = cvh)
  summary_cluster_eeg(clusters = i,model = model, 
                      cv = cvOpt,ix=ix,pvalues = pvalues, eff = eff)
  #summary_hommel_eeg(hommel = hom,ix=ix, alpha = 0.1, clusters = i)
  
  
})

out_d <- t(as.data.frame(out))
rownames(out_d) = NULL
out_d

