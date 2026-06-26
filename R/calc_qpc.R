#' Calculate Qpc
#'
#' This function calculates Qpc given using the eigen decomposition of the kinship matrix and trait values.
#' @param trait_values vector of traits. Not normalized yet.
#' @param eigen_k Named list from eigen decomposition of kinship matrix with 'vectors' and 'values'
#' @param test_pcs the range of PCs you want to test for selection 
#' @param var_pcs the range of PCs used to estimate Va
#' @export
#' @examples
#' calc_qpc()

calc_qpc = function(trait_values, eigen_k, test_pcs, var_pcs) {
  trait_values = trait_values[1:dim(eigen_k$vectors)[1]] - mean(trait_values)
  myCmM = (trait_values %*% eigen_k$vectors[, test_pcs])/sqrt(eigen_k$values[test_pcs])
  myCmL = (trait_values %*% eigen_k$vectors[, var_pcs])/sqrt(eigen_k$values[var_pcs])
  myQm = apply(myCmM, MARGIN = 2, function(n) {
    var0(n)/var0(myCmL)
  })
  p_values = sapply(1:length(test_pcs), function(x) {
    pf(myQm[x], 1, length(var_pcs), lower.tail = F)
  })
  retdf = list(cm = myCmM, cml = myCmL, qm = myQm, pvals = p_values)
  return(retdf)
}
