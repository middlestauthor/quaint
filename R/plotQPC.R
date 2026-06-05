#' Plot a Qpc result
#'
#' This function takes outputs of calc qpc and makes a plot
#' @param myQpc is the output of the CalcQpc function
#' @param myEig is the output of the eigenvalue function applied to the kinship matrix. It includes both vectors and values.
#' @param traitValues is vector of traits. Not normalized yet.
#' @param populationLabels is a list of population labels in the same order as the other inputs
#' @param testPCs is the numbers of the PCs used to test for selection
#' @export


plotQPC = function (myQpc, myEig, traitValues, populationLabels, testPCs) {
  myVaest = var0(myQpc$cml)
  myCI = 1.96 * sqrt(myVaest * myEig$values)
  tib = as_tibble(myEig$vectors[, testPCs])
  namesTib = tibble(pc = paste("PC", testPCs, sep = ""), 
                    ci = paste("ci", testPCs, sep = ""), 
                    ciNeg = paste("ciNeg", testPCs, 
                                  sep = ""))
  tib[, namesTib$ci] = matrix(rep(myCI[1:length(testPCs)], each = nrow(tib)), 
                              ncol = length(testPCs), byrow = F)
  tib[, namesTib$ciNeg] = matrix(-1 * rep(myCI[1:length(testPCs)], each = nrow(tib)), 
                                 ncol = length(testPCs), byrow = F)
  names(tib) = c(namesTib$pc, namesTib$ci, namesTib$ciNeg)
  tib = tib %>% mutate(trait = traitValues[1:nrow(tib)], pop = populationLabels[1:nrow(tib)], 
                       mean = mean(trait))
  plot = lapply(1:length(testPCs), function(x) {
    ggplot(tib, aes(x = .data[[namesTib$pc[x]]], y = trait)) +
      geom_point(aes(color = as.factor(pop))) + 
      geom_smooth(method = "lm", se = F, formula = y ~ x) + 
      geom_abline(aes(slope = .data[[namesTib$ci[x]]], intercept = mean), linetype = 2, color = "red") + 
      geom_abline(aes(slope = .data[[namesTib$ciNeg[x]]], intercept = mean), linetype = 2, color = "red") + 
      theme_classic() +
      labs(x = namesTib$pc[x],
           y = "Trait value", 
           title = paste("QPC - ", namesTib$pc[x], sep = ""))
  })
  return(plot)
}

