brainperm <- function(formula, data, graph, np = 5000,method = NULL, type = "permutation", test = "fisher", aggr_FUN = NULL,
                      threshold = NULL, multcomp = "clustermass", effect = NULL,...){
  
  
  Terms <- terms(formula, special = "Error", data = data)
  indError <- attr(Terms, "specials")$Error
  dotargs = list(...)
  
  
  if (is.null(dotargs$return_distribution)) {
    dotargs$return_distribution = F
  }
  
  if (is.null(dotargs$new_method)) {
    dotargs$new_method = F
  }
  
  
  if(is.null(dotargs$ncores)){dotargs$ncores = detectCores()-1}
  
  if (is.null(dotargs$coding_sum)) {
    switch(test, t = {
      dotargs$coding_sum = F
    }, fisher = {
      dotargs$coding_sum = T
    })
  }
  
  multcomp <- match.arg(multcomp, c("clustermass", "troendle"), several.ok = F)
  
  type <- match.arg(type, c("permutation", "signflip"), several.ok = F)
  
  if (is.null(indError)) {
    result <- brainperm_fix(formula = formula, data = data, method = method, threshold = threshold, np = np, P = P,
                            graph = graph, effect = effect, coding_sum = dotargs$coding_sum, test = test,
                            aggr_FUN = aggr_FUN, multcomp = multcomp, ncores = dotargs$ncores, type = type,
                            return_distribution = dotargs$return_distribution,new_method = dotargs$new_method,
                            rnd_rotation = dotargs$rnd_rotation)
  }
  else if (!is.null(indError)) {
    if (test != "fisher") {
      warning("Random effects model only accept fisher statistics. Test statistic is set to fisher.")
      test = "fisher"
    }
    result <- brainperm_rnd(formula = formula, data = data, method = NULL, threshold = NULL, np = np, P = P,
                            graph = graph, effect = NULL, coding_sum = coding_sum, test = test,
                            aggr_FUN = aggr_FUN, multcomp = multcomp, ncores = ncores, type = type,
                            return_distribution = return_distribution,new_method = new_method)
  }
  return(result)
}