# UKA test methods

#' Function for performing two group UKA test using Goeman's global test package
#'
#' The function takes a kinase-substrate relationship table and performs the global test for each kinase, using the substrates as features and the group labels as response.
#' The weights for the substrates can be specified in the ksr_table
#' The function returns a data frame with the kinase name, p-value and number of substrates used in the test
#' @param X A matrix of substrate measurements, with substrates as columns and samples as rows
#' @param grp A vector of group labels for the samples, must have exactly two levels
#' @param ksr_table A data frame containing the kinase-substrate relationships and weights.
#' The data frame must contain columns for kinase name, substrate ID and weight, which can be specified using the kinase_id, substrate_id and weight_id parameters
#' @param kinase_id The name of the column in ksr_table that contains the kinase
#' @param substarte_id The name of the column in ksr_table that contains the substrate ID, which should match the column names in X
#' @param weight_id The name of the column in ksr_table that contains the weights
#' @param scale Whether to standardize the substrate measurements before performing the test
#' @param directional Whether to perform a directional test (1) or a non-directional
#' @import dplyr globaltest
#' @export
gtUKA = function(X, grp,  ksr_table, kinase_id = "Kinase_Name", substarte_id = "ID", weight_id = "wn", scale = TRUE, directional = 1){
  ksr_table = ksr_table %>%
    select(kinase_id = all_of(kinase_id), substrate_id = all_of(substarte_id), weight_id = all_of(weight_id), everything())

  grp = as.factor(grp)
  if (length(levels(grp)) != 2){
    stop("grp must have exactly two levels")
  }

  result = ksr_table %>%
    group_by(kinase_id) %>%
    do(gtFun(., X, grp, scale, directional)) %>%
    ungroup()
}


#'Utility function for perforiming the global test
#'@import globaltest dplyr
gtFun = function(ksr_table, X, grp, scale, directional){
  Xk = X[, colnames(X) %in% ksr_table$substrate_id, drop = FALSE]
  if (ncol(Xk) < 1){
    kinase_id = ksr_table$kinase_id[1]
    p = NaN
    n = 0
  } else {
    dfw = data.frame(xID = colnames(Xk)) %>%
      left_join(ksr_table %>% select(substrate_id, weight_id), by = c("xID" = "substrate_id"))

    aGt = gt(grp, Xk, weights = dfw$weight_id, standardize = scale, directional = directional)
    p = p.value(aGt)
    n = ncol(Xk)
  }
  result = data.frame(kinase_id = ksr_table$kinase_id[1], p_value = p, n_substrates = n)
}

