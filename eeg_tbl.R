#################Create eeg_lst object##################################
path <- "C:/Users/Angela Andreella/Downloads/"
load(paste0(path,"data_eeg_emotion.RData"))
#signal tbl
#events tbl
#segments_tbl

####signal tbl
###.id = dati$timings$epoch
## .sample = dati$timings$sample with attribute sampling_rate (dati$srate)
## .channels_dbl

.sample <- sample_int(dati$timings$sample,sampling_rate = dati$srate)

.id <- dati$timings$epoch

nam <- dati$chan_info$electrode
for(i in 1:length(dati$chan_info$electrode)){
  assign(nam[i],channel_dbl(
    values = as.numeric(unlist(dati$signals[,i])),
    .x = dati$chan_info$cart_x[i],
    .y = dati$chan_info$cart_y[i],
    .z = dati$chan_info$cart_z[i],
    reference = "",
    number = dati$chan_info$electrode[i],
    unit = "s",
    radius = dati$chan_info$radius[i],
    theta = dati$chan_info$theta[i],
    phi = dati$chan_info$phi[i]
  ))
}


.signal <- lst(.id,.sample,Fp1,Fp2,F3,F4,F7,F8,FC1,FC2,C3,C4,T7,T8,
               CP1, CP2 ,P7  ,  P8  ,  P3  ,  P4,   PO7 ,  PO8 , O1 , O2  ,  Fz  ,  FCz , 
                Cz, CPz,Pz,RM,EOGvo,EOGvu,EOGhl,EOGhr) 
class(.signal) <- c("signal_tbl", "data.table", "data.frame")
data.table::setkey(.signal, .id, .sample)

#.segments, .id, .recording, segment
####.events tbl
###.id Integers indicating to which group the row of the signal matrix belongs.
## .initial = dati$events$event_onset
## .final =  as_sample_int(dati$events$event_onset + dati$events$event_time, sampling_rate =  dati$srate, unit = "s") - 1L,
## .channel = dati$chan_info$electrode
## descriptions_dt = data.table::data.table()


.events <- new_events_tbl (
  .id = dati$events$epoch,
  .initial = dati$events$event_onset %>%
    as_sample_int(sampling_rate = dati$srate, unit = "s"),
  .final = as_sample_int(dati$events$event_onset + dati$events$event_time, 
                         sampling_rate =  dati$srate, unit = "s") - 1L,
  .channel = NA_character_,
  descriptions_dt = data.table::data.table(.recording = dati$epochs$recording, 
                                           .type = dati$events$event_type,
                                           .subj = dati$events$subj))

segments_tbl <- dplyr::tibble(.id =  dati$epochs$epoch, .recording = dati$epochs$recording)
#segments_tbl <- validate_segments(segments_tbl)


data <-   eeg_lst(
  signal_tbl = .signal,
  events_tbl = .events,
  segments_tbl = segments_tbl)

is_eeg_lst(data)

plot(data) 

summary(data)
summary(faces)
channels_tbl(data)
signal_tbl(data)
events_tbl(data)

data_ica <- data %>%
  eeg_ica(method = adapt_fast_ICA, ignore = NULL)
data_ica %>% plot_components()

faces_ica <- data_ica %>%
  eeg_artif_step(Fp1, Fp2, F3, threshold = 30, window = 200, unit = "ms", freq=c(1,10)) %>%
  eeg_artif_peak(Fp1, threshold=100, freq=c(1,10))

s_peaks <- filter(events_tbl(faces_ica),  str_starts(.description, "peak")) %>% pull(.initial)
eeguana:::plot_ica.eeg_ica_lst(faces_ica, samples = seq(s_peaks[1]-2000,s_peaks[1]+2000))



conditions <- c(rep(1,500), rep(2,500), rep(3,500), rep(4,500), rep(5,500))

data %>%
  select(O1, O2) %>%
  ggplot(aes(x = .time, y = .value))+
  geom_line(alpha = .1, aes(group = data$.signal$.id, color = conditions)) +
  stat_summary(
    fun.y = "mean", geom = "line", alpha = 1, size = 1.5,
    aes(color = 10)
  )

design <- expand.grid(subject = c(111, 113, 115, 116, 117, 118, 120, 122, 123, 124, 126, 127, 128, 
                                  130, 131, 132, 134, 135, 137, 138, 139, 140, 141, 142, 143, 144, 
                                  145, 146, 147), stimuli = c("AP", "SED"), action = c("Av_App", "Av_Ev"))

distance_matrix <- dist(cbind(dati$chan_info$cart_x,dati$chan_info$cart_y,dati$chan_info$cart_z))
adjacency_matrix <- as.matrix(distance_matrix) < 35
diag(adjacency_matrix) <- FALSE
dimnames(adjacency_matrix) = list(dati$chan_info$electrode, dati$chan_info$electrode)

graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
graph <- delete_vertices(graph, V(graph)[!get.vertex.attribute(graph, "name")%in%(dati$chan_info$electrode)])

graph <- set_vertex_attr(graph,"x", value = coord[match(vertex_attr(graph,"name"),dati$chan_info$electrode),2])
graph <-set_vertex_attr(graph,"y", value = coord[match(vertex_attr(graph,"name"),dati$chan_info$electrode),3])
graph <-set_vertex_attr(graph,"z", value = coord[match(vertex_attr(graph,"name"),dati$chan_info$electrode),4])
plot(graph)

#I need to compute the 3d array: The first dimension is the permutations, the second the samples, the third is the channels.

#devo creare 
#design: 
design <- expand.grid(subject = dati$events$subj, stimuli = c(2, 4))

#formula: ~ stimuli + Error(subject/(stimuli))
#signal:signal a 3 dimentional array. The the row are the observations, 
#colomn the time and the third dimention are the nodes of the graph.

signal <- as.array(.signal)

brainperm #compute the clustermass test based on multiple signal.

np <- 100 #number of permutations
aggr_FUN <- sum
ncores <- 5
contr <- contr.sum
formula <- ~ stimuli + Error(subject/(stimuli))
pmat <- Pmat(np = np, n = nrow(design))

#compute clusters from permuco4brain
# @param threshold The threshold used to compute the clusters.
# @param aggr_FUN The function that aggregate the cluster into a scalar (cluster mass).
# @param graph A igraph object representing the adjacency of electrod in a scalp.
model <- clustergraph_rnd(formula = formula, data = design, signal = signal, graph = graph, 
                          aggr_FUN = aggr_FUN, method = "Rd_kheradPajouh_renaud", contr = contr, 
                          return_distribution = F, threshold = NULL, ncores = ncores, P = pmat)
