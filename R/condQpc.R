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

condQpc <- function(myZ,myU, myLambdas, myM, myL){
  
  # get Xms for each PC (for Z1, the focal trait and Z2, the correlated trait) THESE ARE MEAN CENTERED
  myX1centered = (myZ[-nrow(myZ),1]-mean(myZ[,1]))%*%myU/sqrt(myLambdas)
  myX2centered = (myZ[-nrow(myZ),2]-mean(myZ[,2]))%*%myU/sqrt(myLambdas)
  
  #get mu' for each PC
  Ca12 = sum(myX1centered[myL]*myX2centered[myL])/length(myL) #is this the right way to do this??
  Va2 = (sum(myX2centered[myL]^2))/length(myL)
  Va1 = sum(myX1centered[myL]^2)/length(myL)
  mu1cond = mean(myZ[-nrow(myZ),1]) + (Ca12/Va2)*(myZ[-nrow(myZ),2] - mean(myZ[,2])) #one value for each individual
  va1cond = Va1 - (Ca12^2)/Va2 
  
  #now test for selection
  myQ = ((myZ[-nrow(myZ),1]-mu1cond)%*%myU[,myM])/sqrt(myLambdas[myM]*va1cond) #get a vector of the projections that we'll test
  #under neutrality, my Q ~ N(0,1)
  return(myQ)
}

 
