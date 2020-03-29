#cluster mass
np <- 1
stimuli <- dati$events$event_type
subject <- dati$events$subj
formula <- signal ~ stimuli + Error(subject/(stimuli))
pmat <- Pmat(np = np, n = nrow(design))

out <- permuco4brain::brainperm(formula, data = design, graph, np = 5000,method = NULL, type = "permutation", test = "fisher", aggr_FUN = NULL,
                                threshold = NULL, multcomp = "clustermass", effect = NULL, P = pmat)
