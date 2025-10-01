compute_information_criteria2 <- function(FOIfit,...){
  
  estimated_parameters <- FOIfit$model$estimated_parameters
  chains <- rstan::extract(FOIfit$fit)
  FOIs <- chains$lambda
  
  sensitivity <- FOIfit$model$se 
  specificity <- FOIfit$model$sp
  
  S <- nrow(FOIs)
  N <- FOIfit$data$N
  A <- FOIfit$data$A 
  Ncategory <- FOIfit$data$Ncategory
  NAgeGroups <- FOIfit$data$NAgeGroups
  
  LogLikelihoods <- matrix(0, nrow = S, ncol = N) # for each iteration and each individual
  Y <- FOIfit$data$Y
  age <- FOIfit$data$age
  category <- FOIfit$data$categoryindex
  
  for (s in seq(1,S)){
    lk  = chains$Like[s,]
    for (i in seq(1,N)){
      if( Y[i] == FALSE){ # if the individual is seronegative
        L = log(1-lk[i])   
      }  else{
        L=log(lk[i])
      } 
      LogLikelihoods[s,i] <- L
    }
  } 
  
  # log-likelihood on the mean lambdas 
  # posterior mean
  
  LogLikelihoodMean <- 0
  P <- (colMeans(chains$P))
  for (i in seq(1,N) ){
    
    age <- FOIfit$data$age[i]
    age_group <- FOIfit$data$age_group[i]
    cat <- category[i]
    p <- P[age,age_group, cat]
    
    if(Y[i] == TRUE){
      LogLikelihoodMean <- LogLikelihoodMean + log(sensitivity-p*(sensitivity+specificity-1) )
    }else{
      LogLikelihoodMean <- LogLikelihoodMean + log(1-sensitivity+p*(sensitivity+specificity-1) )
    }
  }
  
  LP <- rowSums(LogLikelihoods)
  # Compute the AIC
  # Assumes the maximum likelihood is reached 
  AIC <- -2*max(LP) +2*estimated_parameters
  
  # Compute the DIC
  Dbar = -2*mean(LP)
  Dthetabar = -2*LogLikelihoodMean
  pD = Dbar-Dthetabar
  
  DIC = pD+Dbar 
  # Compute the WAIC
  #variance along the column. Each individual has its own variance measured over all sampled parameters 
  
  V = ColVar(exp(LogLikelihoods))
  pwaic = sum(V)
  lpd = sum(log( colSums(exp(LogLikelihoods))/S))
  WAIC <- -2*(lpd-pwaic)
  
  ## Compute the PSIS-LOO using the package loo
  PSIS_LOO = loo(FOIfit$fit, moment_match = TRUE)
  
  information_criteria <- list(AIC = AIC,
                               DIC = DIC,
                               pD = pD,
                               Dbar=Dbar,
                               WAIC = WAIC,
                               pwaic = pwaic,
                               lpd = lpd,
                               k = estimated_parameters,
                               MLE= max(LP),
                               PSIS_LOO=PSIS_LOO)
  
  class(information_criteria) <- "information_criteria"
  return(information_criteria)
  
}

#' @export
print.information_criteria <- function(x,...){
  
  cat(sprintf('AIC:  %f, MLE: %f, k:  %f\n' , x$AIC, x$MLE, x$k))
  cat(sprintf('DIC:  %f, Dbar:  %f, pD:  %f\n' , x$DIC, x$Dbar, x$pD))
  cat(sprintf('WAIC:  %f, pwaic:  %f, lpd: %f \n' , x$WAIC, x$pwaic, x$lpd))
  cat(sprintf('looic:  %f, elpd_loo:  %f, p_loo: %f \n' , x$PSIS_LOO$looic, x$PSIS_LOO$elpd_loo, x$PSIS_LOO$p_loo))
  
}


RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

ColVar <- function(x) {
  colSums((x - colMeans(x))^2)/(dim(x)[1] - 1)
}
