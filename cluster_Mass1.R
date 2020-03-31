#cluster mass

signal <- 
  data_segs_some%>%
  signal_tbl()%>%
  group_by(.id)%>%
  nest()%>%
  mutate(data = map(data,~as.matrix(.x[-1])))%>%
  pull(data)%>%
  invoke(abind::abind,.,along = 3)%>%
  aperm(c(3,1,2))


np <- 100
stimuli <- data_segs_some$.segments$condition
subject <- data_segs_some$.segments$.subj
pmat <- Pmat(np = np, n = nrow(design))

formula <- signal ~ stimuli + Error(subject/(stimuli))
model <- permuco4brain::brainperm(formula = formula, data = design, graph = graph)
