#permuco4brain-with-eeguana
library(tidyverse)
library(eeguana)
library(permuco)
library(permuco4brain)
library(igraph)
library(tidyr)
library(purrr)
library(ARIpermutation)

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

signal <- signal[,100:400,]

data$.segments$condition <- data$.events$.description

design <- 
  segments_tbl(data)%>%
  select(.subj, condition)

graph <- position_to_graph(channels_tbl(data), name = .channel, delta = 53,
                           x = .x, y = .y, z = .z)

igraph::rglplot(graph)
plot(graph)

formula <- signal ~ condition +Error(.subj/(condition))
np = 5000
#pmat <- Pmat(np = np, n = nrow(design))
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

#We visualize the results using a statistical map. 
#Statistics below the threshold are shown in white, 
#non-significant cluster in grey and significant cluster in red-yellow:
image(model)
#The statistics of the cluster are show using the print method. 
#It displays, for each cluster, the number of test, 
#its cluster-mass and the p-value based on the permutation null distribution:
print(model, effect = "condition")

#We need to use the matrix of Pvalues distribution having dimension np times k, where
#k is the number of tests (301*27)
dim(model$multiple_comparison[[1]]$uncorrected$distribution)

Pvalues <- model$multiple_comparison[[1]]$uncorrected$distribution
dim(Pvalues) <- c(dim(Pvalues)[1], dim(Pvalues)[2]*dim(Pvalues)[3])
dim(Pvalues)

lambda <- lambdaOpt(pvalues = Pvalues,family = "Finner", alpha = 0.05, delta = 200)

out <- SingleStepCT(pvalues = Pvalues, ix = c(1,8000), family = "Finner", delta = 10)
out$discoveries

str(model$multiple_comparison[[1]]$clustermass$cluster)
str(model$multiple_comparison[[1]]$clustermass$data)

