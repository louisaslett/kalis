#' @import dplyr
#' @import checkmate
#' @importFrom glue glue glue_collapse
#' @importFrom rlang duplicate
#' @importFrom digest digest
#' @importFrom stringr str_split_fixed
#'
#' @useDynLib kalis, .registration = TRUE, .fixes = "CCall_"
NULL


## ADD FUNCTION TO GET DIMENSIONS OF CACHE AFTER LOADING
## RETURN EG list(num.haps = 2034, hap.length = 10000)

## ADD UNIT TESTS TO MAKE SURE ALL DIFFERENT TYPES MATCH

