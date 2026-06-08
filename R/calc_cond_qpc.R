#' Calculate Conditional Qpc
#'
#' This function calculates conditional Qpc given given using the eigen decomposition of the kinship matrix and trait values.
#' @param focal_trait vector of trait values for trait of interest. Not normalized yet.
#' @param cor_trait vector of trait values for the correlated trait. Not normalized yet
#' @param eigen_k named list from eigen decomposition of kinship matrix with 'vectors' and 'values'
#' @param test_pcs the vector of PCs you want to test for selection 
#' @param var_pcs the vector of the PCs used to estimate Va
#' @export
#' @examples
#' calc_cond_qpc()

calc_cond_qpc = function(focal_trait, cor_trait, eigen_k, test_pcs, var_pcs){
  focal_centered = (focal_trait[-length(focal_trait)] - mean(focal_trait)) %*% eigen_k$vectors/sqrt(eigen_k$values)
  cor_centered = (cor_trait[-length(cor_trait)] - mean(cor_trait)) %*% eigen_k$vectors/sqrt(eigen_k$values)
  cov_a = sum(focal_centered[var_pcs] * cor_centered[var_pcs])/length(var_pcs)
  cor_va = (sum(cor_centered[var_pcs]^2))/length(var_pcs)
  focal_va = sum(focal_centered[var_pcs]^2)/length(var_pcs)
  mu_cond = mean(focal_trait[-length(focal_trait)]) + (cov_a/cor_va) * (cor_trait[-length(cor_trait)] - mean(cor_trait))
  cond_va = focal_va - (cov_a^2)/cor_va
  myQ = ((focal_trait[-length(focal_trait)] - mu_cond) %*% eigen_k$vectors[, test_pcs])/sqrt(eigen_k$values[test_pcs] * cond_va)
  p_values = 1-pnorm(abs(myQ), mean=0, sd=1)
  return(list(myQ = myQ, mu_cond = mu_cond, cond_va = cond_va, pvals = p_values))
}