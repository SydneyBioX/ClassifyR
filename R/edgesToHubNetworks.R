#' Convert a Two-column Matrix or Data Frame into a Hub Node List
#' 
#' Interactions between pairs of features (typically a protein-protein
#' interaction, commonly abbreviated as PPI, database) are restructured into a
#' named list. The name of the each element of the list is a feature and the
#' element contains all features which have an interaction with it.
#' 
#' 
#' @param edges A two-column \code{matrix} or \code{data.frame} for which each
#' row specifies a known interaction betwen two interactors. If feature X
#' appears in the first column and feature Y appears in the second, there is no
#' need for feature Y to appear in the first column and feature X in the
#' second.
#' @param minCardinality An integer specifying the minimum number of features
#' to be associated with a hub feature for it to be present in the result.
#' @return An object of type \code{\link{FeatureSetCollection}}.
#' @author Dario Strbenac
#' @references VAN: an R package for identifying biologically perturbed
#' networks via differential variability analysis, Vivek Jayaswal, Sarah-Jane
#' Schramm, Graham J Mann, Marc R Wilkins and Yee Hwa Yang, 2010, \emph{BMC
#' Research Notes}, Volume 6 Article 430,
#' \url{https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-6-430}.
#' @examples
#' 
#'   interactor <- c("MITF", "MITF", "MITF", "MITF", "MITF", "MITF",
#'                   "KRAS", "KRAS", "KRAS", "KRAS", "KRAS", "KRAS",
#'                   "PD-1")
#'   otherInteractor <- c("HINT1", "LEF1", "PSMD14", "PIAS3", "UBE2I", "PATZ1",
#'                        "ARAF", "CALM1", "CALM2", "CALM3", "RAF1", "HNRNPC",
#'                        "PD-L1")
#'   edges <- data.frame(interactor, otherInteractor, stringsAsFactors = FALSE)
#'   
#'   edgesToHubNetworks(edges, minCardinality = 4)
#' @export
edgesToHubNetworks <- function(edges, minCardinality = 5)
{
  if(is.matrix(edges))
  {
    allFeatures <- unique(as.vector(edges))
  } else { # data.frame
    if(!all(sapply(edges, class) == "character"))
      stop("Both columns must be character type.")
    allFeatures <- unique(unlist(edges))
  }

  featuresByFirstEdge <- split(edges[, 2], factor(edges[, 1], levels = allFeatures))
  featuresBySecondEdge <- split(edges[, 1], factor(edges[, 2], levels = allFeatures))
  featuresByHub <- mapply(c, featuresByFirstEdge, featuresBySecondEdge, SIMPLIFY = FALSE)
  featuresByHub <- lapply(featuresByHub, unique)
  featuresByHub <- featuresByHub[sapply(featuresByHub, length) >= minCardinality]
  if(length(featuresByHub) == 0)
    stop("No features have at least ", minCardinality, " interactors.")
  FeatureSetCollection(featuresByHub)
}
