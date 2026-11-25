#' @keywords internal
"_PACKAGE"

#' @importFrom rlang .data
NULL

# tmerge() evaluates these inside its data arguments
utils::globalVariables(c("ind_id", "tstop", "event", "time_of_vaccination", "tdc"))
