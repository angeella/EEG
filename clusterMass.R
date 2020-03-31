###################Perform clustermass#######################

#load data
load("data_eeg_emotion.RData")
#20 Subjects,
#32 Channels
#Stimuli: pictures. Conditions:
#  1 (f): fear (face)
#  2 (h): happiness (face)
#  3 (d): disgust (face)
#  4 (n): neutral (face)
#  5 (o): object

#We analyze the contrast neutral and happy faces

#For the clusterMass we need:
# 1. design
# 2. signal a 3 dimentional array. The the row are the observations, colomn the time and the third dimention are the nodes of the graph.
# 3. formula
# graph

#1. design
dati1 <- dati
dati1$signals <- dati$signals[,1:27]

design <- expand.grid(subject = unique(dati1$events$subj), stimuli = c(2,4))

data1 <- dplyr::left_join(dati1$timings, dati1$events, by = "epoch")

#2. signal
signal <- list()
for(i in 1:nrow(design)){
  # selecting the experimental condition
  which_id <- design$subject[i] 
  
  # Getting the signal
  data_nonneutral <- dati1$signal[data1$subj == which_id & data1$event_type == 2,]
  data_neutral <- dati1$signal[data1$subj == which_id & data1$event_type == 4,]

  # Storing the signal relative to neutral
  signal[[i]] <- data_nonneutral - data_neutral
  
  }

# creating the 3D array we the appropriate dimension
signal <- abind(signal, along = 3)
signal <- aperm(signal, c(3, 1, 2))

# Select usefull time windows (800ms begining at 200ms at 512hz)

coord <- data.frame(dati$chan_info$electrode,
                    dati$chan_info$cart_x,
                    dati$chan_info$cart_y,
                    dati$chan_info$cart_z)
colnames(coord) <- c("electrode","x","y","z")
coord <- coord[1:27,]


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


np <- 100 #number of permutations
aggr_FUN <- sum
ncores <- 10
contr <- contr.sum
formula <- ~ stimuli + Error(subject/(stimuli))
pmat <- Pmat(np = np, n = nrow(design))

#compute clusters from permuco4brain
# @param threshold The threshold used to compute the clusters.
# @param aggr_FUN The function that aggregate the cluster into a scalar (cluster mass).
# @param graph A igraph object representing the adjacency of electrod in a scalp.
model <- clustergraph_rnd(formula = formula, data = design, signal = signal, graph = graph, 
                          aggr_FUN = aggr_FUN, method = "Rd_kheradPajouh_renaud", contr = contr, 
                          return_distribution = T, threshold = NULL, ncores = ncores, P = pmat)
image(model,effect = 2)
print(model, effect = "stimuli")

