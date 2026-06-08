#' Plot Qpc results
#'
#' This function calculates Qpc given data about the relatedness matrix, and a set of trait values
#' @param qpc the output of the calc_qpc function. 
#' @param trait_values vector of traits for each individual in the same order as input for calc_qpc
#' @param eigen_k named list from eigen decomposition of kinship matrix with list items 'vectors' and 'values'
#' @param test_pcs vector of tested PCs you want to plot
#' @param pop_labels optional vector of population labels for each individual in the same order as input for calc_qpc
#' @export
#' @examples
#' plot_qpc()

plot_qpc = function(qpc, trait_values, eigen_k, test_pcs, pop_labels = NULL) {
  myVaest = var0(qpc$cml)
  conf_int = (1.96 * sqrt(myVaest * eigen_k$values))[test_pcs]
  tib = as_tibble(eigen_k$vectors[, test_pcs])
  name_tib = tibble(pc = paste("PC", test_pcs, sep = ""), 
                    ci = paste("ci", test_pcs, sep = ""), 
                    ciNeg = paste("ciNeg", test_pcs, sep = ""))
  tib[, name_tib$ci] = matrix(rep(conf_int, each = nrow(tib)), 
                              ncol = length(test_pcs), byrow = F)
  tib[, name_tib$ciNeg] = matrix(-1 * rep(conf_int, each = nrow(tib)), 
                                 ncol = length(test_pcs), byrow = F)
  names(tib) = c(name_tib$pc, name_tib$ci, name_tib$ciNeg)
  tib = tib %>% mutate(trait = trait_values[1:nrow(tib)],
                       pop = pop_labels[1:nrow(tib)],
                       mean = mean(trait))
  if(is.null(pop_labels)){
    plot = lapply(test_pcs, function(x) {
      ggplot(tib, aes(x = .data[[paste("PC", x, sep = "")]], y = trait)) +
        geom_point() +
        geom_smooth(method = "lm", se = F, formula = y ~ x) +
        geom_abline(aes(slope = .data[[paste("ci", x, sep = "")]], intercept = mean), linetype = 2, color = "red") +
        geom_abline(aes(slope = .data[[paste("ciNeg", x, sep = "")]], intercept = mean), linetype = 2, color = "red") +
        theme_classic() +
        labs(x = paste("PC", x, sep = ""),
             y = "Trait value", 
             title = paste("QPC - ", paste("PC", x, sep = ""), sep = ""))
    })
  }
  else{
    plot = lapply(test_pcs, function(x) {
      ggplot(tib, aes(x = .data[[paste("PC", x, sep = "")]], y = trait)) +
        geom_point(aes(color = as.factor(pop))) +
        geom_smooth(method = "lm", se = F, formula = y ~ x) +
        geom_abline(aes(slope = .data[[paste("ci", x, sep = "")]], intercept = mean), linetype = 2, color = "red") +
        geom_abline(aes(slope = .data[[paste("ciNeg", x, sep = "")]], intercept = mean), linetype = 2, color = "red") +
        theme_classic() +
        labs(x = paste("PC", x, sep = ""),
             y = "Trait value", 
             title = paste("QPC - ", paste("PC", x, sep = ""), sep = ""),
             color = "Population")
    })
  }
  return(plot)
}