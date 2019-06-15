#' Get Useful Predictors
#'
#'Extract the variables used in subgroup definitions (partykit objects)
#'
#' @param x partykit model object
#'
#' @return Variable names that define the subgroups in the partykit model object

### Extract Predictors Used in Party/MOB Tree ###
getUsefulPredictors <- function(x) {
  varid <- nodeapply(x, ids = nodeids(x),
                     FUN = function(n) split_node(n)$varid)
  varid <- unique(unlist(varid))
  names(data_party(x))[varid]
}
