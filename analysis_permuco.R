library(ARIpermutation)
library(ARIeeg)
load("C:/Users/Angela Andreella/Documents/Thesis_Doc/Permutation_based_ARI/EEG/model_permuco.RData")

eff = "`stimuli:action`"
alternative = "two.sided"

Test <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$uncorrected$distribution")))
dim(Test) <- c(dim(Test)[1], dim(Test)[2]*dim(Test)[3])

#test_obs <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$data$statistic")))
#TT <- rbind(test_obs,Test)
pvalues <- switch(alternative, 
                  "two.sided" =  matrixStats::colRanks(-abs(Test)) / nrow(Test),
                  "greater" = matrixStats::colRanks(-Test) / nrow(Test),
                  "less" = matrixStats::colRanks(Test) / nrow(Test))


pvalues <-t(pvalues)
#pvalues <- pvalues[,which(model$multiple_comparison[[1]]$clustermass$data$cluster_id!=0)]
pvalues_ord <- rowSortC(pvalues)

family = "beta"
delta = NULL
ct = c(0,1)
alpha = 0.1
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

deltaV <- seq(0,1500,1)
family = "Simes"
OUT <- list()
for(i in 1:length(deltaV)){
  
  lambda <- lambdaOpt(pvalues = pvalues_ord, family = family, ct = ct, alpha = alpha, delta = deltaV[i]) 
  cvOpt = cv(pvalues = pvalues_ord, family = family, alpha = alpha, lambda= lambda, delta = deltaV[i])
  
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
  OUT[[i]] <- as.vector(out_d[5,5])
  
}


