#' List Available Feature Selection and Classification Approaches
#' 
#' Prints a list of keywords to use with \code{\link{crossValidate}}
#' 
#' 
#' @aliases available
#' @param what Default: \code{"classifier"}. Either \code{"classifier"}, \code{"selectionMethod"}
#' or \code{"multiViewMethod"}.
#' @return A data frame of keywords and their descriptions. For classifiers and feature selection methods, the column
#' \code{Outcomes} says whether each one accepts classes, survival outcomes or both.

#' @author Dario Strbenac
#' @examples
#' 
#'  available()
#'  
#' @export

available <- function(what = c("classifier", "selectionMethod", "multiViewMethod"))
{
  what <- match.arg(what)
  switch(what, 
    selectionMethod = .ClassifyRenvir[["selectKeywords"]],
    classifier = .ClassifyRenvir[["classifyKeywords"]],
    multiViewMethod = .ClassifyRenvir[["multiViewKeywords"]])
}