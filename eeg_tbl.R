#################Create eeg_lst object##################################

source("library.R")
source("utils.R")
source("utils1.R")
source("new_events_tbl.R")
source("validate_eeg_lst.R")
source("signal_tbl.R")

#load data
load("data_eeg_emotion.RData")

#Builds an eeg_lst object composed of two data.table`objects and one tibble.
# All three are linked by a unique identifier `.id`. Amplitude values and timestamps appear in the `signal`
#table. Triggers, blinks, artifact rejection markings, and other
#events logged by the EEG recording software appear in the events table.
#Segment information and recording IDs appear in the segments tibble.

#signal tbl
#events tbl
#segments_tbl

####signal tbl

#' The `signal` table is organised into columns representing timestamps
#' (`.sample`) and individual electrodes. Each `.sample` corresponds to
#' 1 sample in the original recording, i.e. if the sampling rate of the EEG
#' recording is 500 Hz, then each `.sample` corresponds to 2 milliseconds.
#' These timestamps correspond to `.initial` in the `events` table, which
#' displays only the timestamps where logged events began.

###.id = dati$timings$epoch
## .sample = dati$timings$sample with attribute sampling_rate (dati$srate)
## .channels_dbl

.sample <- sample_int(dati$timings$sample,sampling_rate = dati$srate)

.id <- dati$timings$epoch

ch_n <- sapply(c(1:32), function(x) paste0("Ch",x))

nam <- dati$chan_info$electrode
for(i in 1:length(dati$chan_info$electrode)){
  assign(nam[i],channel_dbl(
    values = as.numeric(unlist(dati$signals[,i])),
    .x = dati$chan_info$cart_x[i],
    .y = dati$chan_info$cart_y[i],
    .z = dati$chan_info$cart_z[i],
    reference = "",
    resolution = 1,
    number = ch_n[i],
    unit = "microvolt",
    radius = dati$chan_info$radius[i],
    theta = dati$chan_info$theta[i],
    phi = dati$chan_info$phi[i]
  ))
}


.signal <- data.table::data.table(.id,.sample,Fp1,Fp2,F3,F4,F7,F8,FC1,FC2,C3,C4,T7,T8,
               CP1, CP2 ,P7  ,  P8  ,  P3  ,  P4,   PO7 ,  PO8 , O1 , O2  ,  Fz  ,  FCz , 
                Cz, CPz,Pz,RM,EOGvo,EOGvu,EOGhl,EOGhr)
#class(.signal) <- c("signal_tbl", "data.table", "data.frame")
#data.table::setkey(.signal, .id, .sample)
.signal[, .id := .id][, .sample := .sample]
data.table::setnames(.signal, make_names(colnames(.signal)))
data.table::setcolorder(.signal, c(".id", ".sample"))
data.table::setattr(.signal, "class", c("signal_tbl", class(.signal)))
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
  #.final = as_sample_int(dati$events$event_onset, 
  #                       sampling_rate =  dati$srate, unit = "s"),
  .channel = NA_character_,
  descriptions_dt = data.table::data.table(
 # .recording = dati$epochs$recording, 
                                           .type = rep("Stimulus",length(dati$events$event_type)),
                                           .description = dati$events$event_type))
 #                                          .subj = dati$events$subj))

data.table::setkey(.events, .id)
#segments_tbl <- dplyr::tibble(.id =  dati$epochs$epoch, .recording = dati$epochs$recording)
segments_tbl <- dplyr::tibble(.id = dati$epochs$epoch, 
                              .recording = dati$epochs$recording,
                              conditions = dati$events$event_type)
                              #segment = rep(1,length(dati$epochs$recording)))

segments_tbl <- validate_segments(segments_tbl)

data <-   eeg_lst(
  signal_tbl = .signal,
  events_tbl = .events,
  segments_tbl = segments_tbl)

is_eeg_lst(data)
is_events_tbl(data$.events)
is_signal_tbl(data$.signal)

plot(data) 

summary(data)


channels_tbl(data)
signal_tbl(data)
events_tbl(data)

#other plot
plot_components(data)

data_filter <- data %>%
  filter(conditions %in% c(2, 4)) %>%
  mutate(
    condition =
      if_else(conditions == 2, "happy", "neutral")
  )

data_filter %>%
  select(O1, O2, P7, P8) %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(alpha = .1, aes(group = .id, color = condition)) +
  stat_summary(
    fun.y = "mean", geom = "line", alpha = 1, size = 1.5,
    aes(color = condition)
  ) +
  facet_wrap(~.key) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = .17, linetype = "dotted") +
  theme(legend.position = "bottom")

data_filter %>%
  filter(between(as_time(.sample, unit = "milliseconds"), 100, 200)) %>%
  group_by(condition) %>%
  summarize_at(channel_names(.), mean, na.rm = TRUE) %>%
  plot_topo() +
  annotate_head() +
  geom_contour() +
  geom_text(colour = "black") +
  facet_grid(~condition)




#We consider the happy faces and netural stimuli
data_segs <- data %>%
  eeguana::eeg_segment(.description %in% c(2,4), lim = c(0,0.2))


segments_tbl(data_segs)

data1 <- merge(data$.signal,data$.events,by = ".id", all.x = TRUE) %>% filter(.description %in% c(2,4))

data_segs$.signal <- data.table::data.table(data1[1:34])
a<-as.factor(data_segs$.signal$.id)
levels(a) <- c(1:40)
data_segs$.signal$.id <- as.integer(a)
data_segs$.signal[, .id := .id][, .sample := .sample]
data.table::setnames(data_segs$.signal, make_names(colnames(data_segs$.signal)))
data.table::setcolorder(data_segs$.signal, c(".id", ".sample"))
data.table::setattr(data_segs$.signal, "class", c("signal_tbl", class(data_segs$.signal)))

data.table::setkey(data_segs$.signal, .id, .sample)



data_segs %>%
  select(Fz, FCz, Cz, CPz) %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(alpha = .1, aes(group = .id, color = condition)) +
  stat_summary(
    fun.y = "mean", geom = "line", alpha = 1, size = 1.5,
    aes(color = condition)
  ) +
  facet_wrap(~.key) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = .17, linetype = "dotted") +
  theme(legend.position = "bottom")



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

