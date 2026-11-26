#' @importFrom stats aggregate binomial chisq.test density dist dnorm fisher.test glm hclust mad median model.frame
#' model.matrix na.omit na.pass pnorm prcomp predict quantile quasibinomial sd setNames splinefun var weighted.mean
#' @importFrom utils combn head tail
#' @importFrom S4Vectors first second mcols<-
NULL

.ClassifyRenvir <- new.env(parent = emptyenv())

# Used internally during parameter selection based on best performance.
.ClassifyRenvir[["performanceInfoTable"]] <- matrix(c("Error", "lower",
                                                      "Accuracy", "higher",
                                                      "Balanced Error", "lower",
                                                      "Balanced Accuracy", "higher",
                                                      "Micro Precision", "higher",
                                                      "Micro Recall", "higher",
                                                      "Micro F1", "higher",
                                                      "Macro Precision", "higher",
                                                      "Macro Recall", "higher",
                                                      "Macro F1", "higher",
                                                      "Matthews Correlation Coefficient", "higher",
                                                      "AUC", "higher",
                                                      "C-index", "higher"),
                                                    ncol = 2, byrow = TRUE, dimnames = list(NULL, c("type", "better"))
) |> as.data.frame()

.ClassifyRenvir[["performanceTypes"]] <- .ClassifyRenvir[["performanceInfoTable"]][, "type"]

# The classifiers and feature selection methods, their keywords and display names are in the registry (registry.R).

.ClassifyRenvir[["multiViewKeywords"]] <- matrix(
  c("none", "Keep assays separate.",
    "merge", "Concatenate all selected feaures into single table before modelling.",
    "prevalidation", "Reduce each assay into a vector and concatenate to clinical data before modelling.",
    "PCA", "Reduce each assay into a lower dimensional representation and concatenate the principal components to the clinical data before modelling."
    ),
  ncol = 2, byrow = TRUE, dimnames = list(NULL, c("multiViewMethod Keyword", "Description"))
) |> as.data.frame()

.ClassifyRenvir[["prepareDataFormals"]] <- c("useFeatures", "maxMissingProp", "maxSimilarity", "topNvariance")