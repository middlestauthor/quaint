#' Calculate conditional Qpc
#'
#' This function calculates Qpc given data about the relatedness matrix, and a set of trait values
#' @param myZ two columned data frame containing vectors of two traits. Not normalized yet. The first trait is the one we're testing for selection on conditional on the second trait.
#' @param myU matrix of eigenvectors of the kinship matrix (each column is an eigenvector)
#' @param myLambdas vector of eigenvalues of the kinship matrix 
#' @param myM the range of PCs you want to test for selection 
#' @param myL the range of PCs used to estimate Va
#' @export
#' @examples
#' calcQpc()

condQpc = function (myZ, myU, myLambdas, myM, myL) 
{
  myX1centered = (myZ[-nrow(myZ), 1] - mean(myZ[, 1])) %*% 
    myU/sqrt(myLambdas)
  myX2centered = (myZ[-nrow(myZ), 2] - mean(myZ[, 2])) %*% 
    myU/sqrt(myLambdas)
  Ca12 = sum(myX1centered[myL] * myX2centered[myL])/length(myL)
  Va2 = (sum(myX2centered[myL]^2))/length(myL)
  Va1 = sum(myX1centered[myL]^2)/length(myL)
  mu1cond = mean(myZ[-nrow(myZ), 1]) + (Ca12/Va2) * (myZ[-nrow(myZ), 
                                                         2] - mean(myZ[, 2]))
  va1cond = Va1 - (Ca12^2)/Va2
  myQ = ((myZ[-nrow(myZ), 1] - mu1cond) %*% myU[, myM])/sqrt(myLambdas[myM] * 
                                                               va1cond)
  myPs = 1-pnorm(abs(myQ), mean=0, sd=1)
  return(list(myQ = myQ, mu1cond = mu1cond, va1cond = va1cond, pvals = myPs))
}