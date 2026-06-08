#' Plot Conditional Qpc results
#'
#' This function calculates Qpc given data about the relatedness matrix, and a set of trait values
#' @param cond_qpc the output of the calc_cond_qpc function
#' @param trait_values vector of focal trait values for each individual in the same order as eigen_k
#' @param eigen_k Named list from eigen decomposition of kinship matrix with list items 'vectors' and 'values'
#' @param test_pcs vector of tested PCs you want to plot
#' @param pop_labels optional vector of population labels for each individual in the same order as input for calc_qpc
#' @export
#' @examples
#' plot_cond_qpc()

plot_cond_qpc = function(cond_qpc, eigen_k, trait_values, test_pcs, pop_labels = NULL){
  conf_int = (1.96 * sqrt(cond_qpc$cond_va * eigen_k$values))[test_pcs]
  var_pcsM = lm(cond_qpc$mu_cond~eigen_k$vectors[,test_pcs]) ##conditional expectation
  my_coefs = var_pcsM$coefficients[-1]
  names(my_coefs) = paste("PC", test_pcs, sep = "")
  tib = as_tibble(eigen_k$vectors[,test_pcs])
  name_tib = tibble(pc = paste("PC", test_pcs, sep = ""),
                    ci = paste("ci", test_pcs, sep = ""),
                    ciNeg = paste("ciNeg", test_pcs, sep = ""))
  tib[,name_tib$ci] = matrix(rep(my_coefs[paste("PC", test_pcs, sep = "")] + conf_int,
                                 each = nrow(tib)), ncol = length(test_pcs), byrow = F)
  tib[,name_tib$ciNeg] = matrix(rep(my_coefs[paste("PC", test_pcs, sep = "")] - conf_int,
                                    each = nrow(tib)), ncol = length(test_pcs), byrow = F)
  names(tib) = c(name_tib$pc, name_tib$ci, name_tib$ciNeg)
  tib = tib %>% mutate(trait = trait_values[1:nrow(tib)],
                       pop = pop_labels[1:nrow(tib)],
                       mean = mean(trait))
  if(is.null(pop_labels)){
    plots = lapply(test_pcs, function(x){
      ggplot(tib, aes(x = .data[[paste("PC", x, sep = "")]], y = trait)) +
        geom_point() +
        geom_smooth(method = "lm", se = F, formula = y~x) +
        geom_abline(aes(slope = .data[[paste("ci", x, sep = "")]], intercept = mean), linetype = 2, color = "red") + 
        geom_abline(aes(slope = .data[[paste("ciNeg", x, sep = "")]], intercept = mean), linetype = 2, color = "red") + 
        theme_classic() +
        labs(x = paste("PC", x, sep = ""),
             y = "Trait value", 
             title = paste("Conditional QPC - ", paste("PC", x, sep = ""), sep = ""))
    })
  }
  else{
    plots = lapply(test_pcs, function(x){
      ggplot(tib, aes(x = .data[[paste("PC", x, sep = "")]], y = trait)) +
        geom_point(aes(color = as.factor(pop))) +
        geom_smooth(method = "lm", se = F, formula = y~x) +
        geom_abline(aes(slope = .data[[paste("ci", x, sep = "")]], intercept = mean), linetype = 2, color = "red") + 
        geom_abline(aes(slope = .data[[paste("ciNeg", x, sep = "")]], intercept = mean), linetype = 2, color = "red") + 
        theme_classic() +
        labs(x = paste("PC", x, sep = ""),
             y = "Trait value", 
             title = paste("Conditional QPC - ", paste("PC", x, sep = ""), sep = ""))
    })
  }
  return(plots)
}