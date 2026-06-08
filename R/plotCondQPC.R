#' Plot a Qpc result
#'
#' This function takes outputs of calc qpc and makes a plot
#' @param myCondQpc is the output of the CondQpc function
#' @param myEig is the output of the eigenvalue function applied to the kinship matrix. It includes both vectors and values.
#' @param traitValues is vector of traits. Not normalized yet.
#' @param populationLabels is a list of population labels in the same order as the other inputs
#' @param testPCs is the numbers of the PCs used to test for selection
#' @export


plotCondQpc = function(myCondQpc, myEig, traitValues, populationLabels, testPCs){
  myCI = 1.96*sqrt(myCondQpc$va1cond*myEig$values)
  myLM = lm(myCondQpc$mu1cond~myEig$vectors[,testPCs]) ##conditional expectation
  tib = as_tibble(myEig$vectors[,testPCs])
  namesTib = tibble(pc = paste("PC", testPCs, sep = ""),
                    ci = paste("ci", testPCs, sep = ""),
                    ciNeg = paste("ciNeg", testPCs, sep = ""))
  tib[,namesTib$ci] = matrix(rep(myLM$coefficients[(1:length(testPCs))+1] + myCI[(1:length(testPCs))],
                                 each = nrow(tib)), ncol = length(testPCs), byrow = F)
  tib[,namesTib$ciNeg] = matrix(rep(myLM$coefficients[(1:length(testPCs))+1] - myCI[(1:length(testPCs))],
                                    each = nrow(tib)), ncol = length(testPCs), byrow = F)
  names(tib) = c(namesTib$pc, namesTib$ci, namesTib$ciNeg)
  tib = tib %>% mutate(trait = traitValues[1:nrow(tib)],
                       pop = as.factor(populationLabels[1:nrow(tib)]),
                       mean = mean(trait))
  plots = lapply(1:length(testPCs), function(x){
    ggplot(tib, aes(x = .data[[namesTib$pc[x]]], y = trait)) +
      geom_point(aes(color = pop)) +
      geom_smooth(method = "lm", se = F, formula = y~x) +
      geom_abline(aes(slope = .data[[namesTib$ci[x]]], intercept = mean), linetype = 2, color = "red") + 
      geom_abline(aes(slope = .data[[namesTib$ciNeg[x]]], intercept = mean), linetype = 2, color = "red") + 
      theme_classic() +
      labs(x = namesTib$pc[x],
           y = "Trait value", 
           title = paste("Conditional QPC - ", namesTib$pc[x], sep = ""))
  })
  return(plots)
}