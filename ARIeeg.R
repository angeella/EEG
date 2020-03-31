#' @title ARI Permutation-based for EEG data
#' @description Performs ARI using permutation local test for EEG data
#' @usage ARIpermCT(model, alpha, summary_stat, silent, family, delta, ct)
#' @param model clusterMass model from permuco4brain package
#' @param alpha alpha level
#' @param summary_stat Choose among "max", "center-of-mass"
#' @param silent FALSE by default.
#' @param family which family for the confidence envelope? simes, finner, beta or higher.criticism. default is simes
#' @param delta do you want to consider at least delta size set?
#' @param B number of permutation, default 1000
#' @param ct set of thresholds
#' @author Angela Andreella
#' @return Returns a list with the following objects: discoveries number of discoveries in the set selected, cluster id, maximum test statistic and relative coordinates
#' @export
#' @importFrom plyr laply
#' 
ARIeeg <- function(model, alpha, summary_stat, silent, family, delta, ct){
  
  pvalues <- model$multiple_comparison[[1]]$uncorrected$distribution
  pv_obs <- model$multiple_comparison[[1]]$clustermass$data$pvalue
  dim(pvalues) <- c(dim(pvalues)[1], dim(pvalues)[2]*dim(pvalues)[3])
  pvalues <-rbind(pv_obs,pvalues)
  
  pvalues_ord <- rowSortC(pvalues)
  lambda <- lambdaOpt(pvalues = pvalues_ord, family = family, ct = ct, alpha = alpha, delta = delta) 
  cvOpt = cv(pvalues = pvalues_ord, family = family, alpha = alpha, lambda= lambda, delta = delta)
  
  clstr_id <- model$multiple_comparison[[1]]$clustermass$data$cluster_id
  clusters <- c(1:model$multiple_comparison[[1]]$clustermass$cluster$no)
  
  out=lapply(clusters,function(i){
    ix= which(clstr_id == i)
    
    #cluster_ids=which(ix,arr.ind = TRUE)
    #cluster_ids=cbind(cluster_ids,Stat=StatFun(ix))
    #Error if I put pvalues[,mask] instead of pvalues in SingleStepCT
    #perm <- SingleStepCT(pvalues = pvalues,ct =ct, ix =as.vector(which(ix[mask])), alpha = alpha, shift = shift, family = 'Simes', lambda = lambda)
    #perm <- discoveriesPerm(praw = praw, ix = ix[mask], cvh = cvh)
    c(summary_perm_eeg(cv = cvOpt,ix=ix,pvalues = pvalues),
             summary_cluster_eeg(i,model)
    )
  })
  
  out_d <- t(as.data.frame(out))
  rownames(out_d) = NULL
  out_d
}

summary_perm_eeg <- function(cv,ix,pvalues){
  #idix <- which(ix)
  p <- pvalues[1,ix]
  Total = length(p)
  False_Null= dI(ix = ix,cv = cv,praw = pvalues[1,])
  True_Null=Total - False_Null
  Active_Proportion= False_Null / Total
  out = c(Size=Total,FalseNull=False_Null,TrueNull=True_Null,ActiveProp=Active_Proportion)
  out
}

summary_cluster_eeg <- function(clusters,model){

  info <- model$multiple_comparison[[1]]$clustermass$cluster
  clustermass <- model$multiple_comparison[[1]]$clustermass$cluster$clustermass[clusters]
  pvalue <- model$multiple_comparison[[1]]$clustermass$cluster$pvalue[clusters]
  csize <- model$multiple_comparison[[1]]$clustermass$cluster$csize[clusters]
  #membership <- toString(names(model$multiple_comparison[[1]]$clustermass$cluster$membership[model$multiple_comparison[[1]]$clustermass$cluster$membership ==clusters]))
  out=c(clustermass,pvalue,csize)
  names(out)=c("clustermass", "pvalue", "csize")
  #out=c(clustermass,pvalue,csize, membership)
  #names(out)=c("clustermass", "pvalue", "csize", "membership")
  out
}