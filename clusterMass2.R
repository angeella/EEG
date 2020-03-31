#permuco4brain-with-eeguana
library(tidyverse)
library(eeguana)
library(permuco)
library(permuco4brain)
library(igraph)
library(tidyr)
library(purrr)

load("data.RData")

channel_names(data)

data <- data %>% filter(condition %in% c(1,4))  %>%
  mutate(
    condition =
      if_else(condition == 1, "object", "neutral")
  ) 

signal <- 
  data%>%
  signal_tbl()%>%
  group_by(.id)%>%
  nest()%>%
  mutate(data = map(data,~as.matrix(.x[-1])))%>%
  pull(data)%>%
  invoke(abind::abind,.,along = 3)%>%
  aperm(c(3,1,2))

dim(signal)
#First dim is the observations, second time, third node
#We take a subset

#signal_sub <- signal[,100:400,]

data$.segments$condition <- data$.events$.description

design <- 
  segments_tbl(data)%>%
  select(.subj, condition)

graph <- position_to_graph(channels_tbl(data), name = .channel, delta = 53,
                           x = .x, y = .y, z = .z)

igraph::rglplot(graph)
plot(graph)

formula <- signal ~ condition +Error(.subj/(condition))
np = 200
pmat <- Pmat(np = np, n = nrow(design))
model <- permuco4brain::brainperm(  formula = formula,
                                    data = design,
                                    graph = graph,
                                    np = np,
                                    method = NULL,
                                    type = "permutation",
                                    test = "fisher",
                                    aggr_FUN = NULL,
                                    threshold = NULL,
                                    multcomp = "clustermass",
                                    effect = NULL,
                                    return_distribution = TRUE)

dim(model$model.matrix)
dim(model$multiple_comparison[[1]]$uncorrected$distribution)



