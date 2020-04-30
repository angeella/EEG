data=data_seg 
alpha=0.1
family="Simes"
ct=c(0,1)
alternative = "two.sided"
timeS = c(100:500)
dist = 53
formula = signal ~ condition + Error(.subj/(condition))
variable = c(".subj", "condition")
B = 5000
eff = "condition"
signal <- 
  data%>%
  signal_tbl()%>%
  group_by(.id)%>%
  nest()%>%
  mutate(data = map(data,~as.matrix(.x[-1])))%>%
  pull(data)%>%
  invoke(abind,.,along = 3)%>%
  aperm(c(3,1,2))

if(!is.null(timeS)){signal <- signal[,timeS,]}

design <- 
  segments_tbl(data)%>%
  select(variable)

graph <- position_to_graph(channels_tbl(data), name = .channel, delta = dist,
                           x = .x, y = .y, z = .z)

#formula <- signal ~ condition + Error(.subj/(condition))

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

lambdaV <- c(0:3254)
OUT <- list()
for(i in c(1:length(lambdaV))){
  
  lambda <- lambdaOpt(pvalues = pvalues_ord, family = family, ct = ct, alpha = alpha, delta = lambdaV[i]) 
  cvOpt = cv(pvalues = pvalues_ord, family = family, alpha = alpha, lambda= lambda, delta = lambdaV[i])
  
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
  OUT[[i]] <- as.vector(out_d[which(model$multiple_comparison$condition$clustermass$cluster$pvalue<0.05),7])
  
  
}