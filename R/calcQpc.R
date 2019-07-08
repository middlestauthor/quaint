#' Calculate Qpc
#'
#' This function calculates Qpc given data about the relatedness matrix, and a set of trait values
#' @param myZ vector of traits. Not normalized yet.
#' @param myU matrix of eigenvectors of the kinship matrix (each column is an eigenvector)
#' @param myLambdas vector of eigenvalues of the kinship matrix 
#' @param myM the range of PCs you want to test for selection 
#' @param myL the range of PCs used to estimate Va
#' @export
#' @examples
#' calcQpc()

calcQpc <- function(myZ, myU, myLambdas, myL, myM){
  myZ = myZ[1:dim(myU)[1]] - mean(myZ) #mean center phenotypes
  myCmM = (myZ %*% myU[,myM])/sqrt(myLambdas[myM]) #project + standardize by the eigenvalues for testing for selection
  myCmL = (myZ %*% myU[,myL])/sqrt(myLambdas[myL]) #project + standardize by the eigenvalues for estimating Va
  myQm = sapply(myM, function(n){var0(myCmM[n])/var0(myCmL) })  #test for selection
  myPs = sapply(1:pcm, function(x){pf(myQm[x], 1, length(myL), lower.tail=F)}) #get a pvalue
  retdf = list(cm = myCm, qm = myQm, pvals = myPs)
  return(retdf)
  }



 
