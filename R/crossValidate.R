#' Cross-validation to evaluate classification performance.
#' 
#' This function has been designed to facilitate the comparison of classification
#' methods using cross-validation, particularly when there are multiple assays per biological unit.
#' A selection of typical comparisons are implemented. The \code{train} function
#' is a convenience method for training on one data set and likewise \code{predict} for predicting on an
#' independent validation data set.
#'
#' @param measurements Either a \code{\link{DataFrame}}, \code{\link{data.frame}}, \code{\link{matrix}}, \code{\link{MultiAssayExperiment}} 
#' or a list of the basic tabular objects containing the data.
#' @param x Same as \code{measurements} but only training samples.
#' @param outcome A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurements} or a character vector of length 1 containing the
#' column name in \code{measurements} if it is a \code{\link{DataFrame}}. Or a \code{\link{Surv}} object or a character vector of
#' length 2 or 3 specifying the time and event columns in \code{measurements} for survival outcome. If \code{measurements} is a
#' \code{\link{MultiAssayExperiment}}, the column name(s) in \code{colData(measurements)} representing the outcome.  If column names
#' of survival information, time must be in first column and event status in the second.
#' @param outcomeTrain For the \code{train} function, either a factor vector of classes, a \code{\link{Surv}} object, or
#' a character string, or vector of such strings, containing column name(s) of column(s)
#' containing either classes or time and event information about survival. If column names
#' of survival information, time must be in first column and event status in the second.
#' @param extraParams A list of parameters that will be used to overwrite default settings of transformation, selection, or model-building functions or
#' parameters which will be passed into the data cleaning function or cross-validation mode used for parameter tuning. Each name of a list element is a list and must be one of \code{"prepare"},
#' \code{"select"}, \code{"train"}, \code{"predict"}, \code{tuneCross}. By default, no parameter tuning is done. To use the a default parameter range for tuning (see the article titled Parameter Tuning Presets for crossValidate and Their Customisation on
#' the website), specify a list element of \code{"select"} or \code{"train"} lists named \code{"tuneParams"} with value \code{"auto"}. To specify your own range of values, specify a \code{list} with names being the parameters in the functions
#' described in the same article on the website. For the valid element names in the \code{"prepare"} list, see \code{?prepareData} for its parameter names. The list \code{"tuneCross"} can have elements named \code{"tuneMode"} and \code{"performanceType"}. Valid values for \code{"tuneMode"} are \code{"Resubstitution"} or \code{"Nested CV"}. For \code{"performanceType"}, it is any of the metrics which can be specified to \code{\link{calcPerformance}}.
#' @param nFeatures The number of features to be used for classification. If a named vector with the same names of multiple assays, 
#' a different number of features will be used for each assay. Set to \code{"all"} if all features should be used. To tune it, specify a vector or \code{list} of named vectors to \code{"tuneParams"} list of
#' \code{"select"} element list of \code{extraParams} list.
#' @param selectionMethod Default: \code{"auto"}. A character vector of feature selection methods to compare. If a named character vector with names corresponding to different assays, 
#' and performing multiview classification, the respective selection methods will be used on each assay. If \code{"auto"}, t-test (two categories) / F-test (three or more categories) ranking
#' and top \code{nFeatures} optimisation is done. Otherwise, the ranking method is per-feature Cox proportional hazards p-value. \code{"none"} is also a valid value, meaning that no
#' feature selection prior to model building will be performed (but implicit selection might still happen with the classifier).
#' @param classifier Default: \code{"auto"}. A character vector of classification methods to compare. If a named character vector with names corresponding to different assays, 
#' and performing multiview classification, the respective classification methods will be used on each assay. If \code{"auto"}, then a random forest is used for a classification
#' task or Cox proportional hazards model for a survival task.
#' @param multiViewMethod Default: \code{"none"}. A character vector specifying the multiview method or data integration approach to use. See \code{available("multiViewMethod") for possibilities.}
#' @param assayCombinations A character vector or list of character vectors proposing the assays or, in the case of a list, combination of assays to use
#' with each element being a vector of assays to combine. Special value \code{"all"} means all possible subsets of assays.
#' @param nFolds A numeric specifying the number of folds to use for cross-validation.
#' @param nRepeats A numeric specifying the the number of repeats or permutations to use for cross-validation.
#' @param nCores A numeric specifying the number of cores used if the user wants to use parallelisation. 
#' @param characteristicsLabel A character specifying an additional label for the cross-validation run.
#' @param ... For \code{train} and \code{predict} functions, parameters not used by the non-DataFrame signature functions but passed into the DataFrame signature function.
#' @param object A trained model to predict with.
#' @param newData The data to use to make predictions with.
#' @param verbose Default: 0. A number between 0 and 3 for the amount of
#' progress messages to give.  A higher number will produce more messages as
#' more lower-level functions print messages.
#'
#' @details
#' \code{classifier} can be any a keyword for any of the implemented approaches as shown by \code{available()}.
#' \code{selectionMethod} can be a keyword for any of the implemented approaches as shown by \code{available("selectionMethod")}.
#' \code{multiViewMethod} can be a keyword for any of the implemented approaches as shown by \code{available("multiViewMethod")}.
#'
#' @return An object of class \code{\link{ClassifyResult}}
#' @export
#' @aliases crossValidate crossValidate,matrix-method crossValidate,DataFrame-method
#' crossValidate,MultiAssayExperiment-method, crossValidate,data.frame-method
#' @rdname crossValidate
#'
#' @examples
#' 
#' data(asthma)
#' 
#' # Compare randomForest and SVM classifiers.
#' result <- crossValidate(measurements, classes, classifier = c("randomForest", "SVM"))
#' performancePlot(result)
#' 
#' 
#' # Compare performance of different assays. 
#' # First make a toy example assay with multiple data types. We'll randomly assign different features to be clinical, gene or protein.
#' # set.seed(51773)
#' # measurements <- DataFrame(measurements, check.names = FALSE)
#' # mcols(measurements)$assay <- c(rep("clinical", 20), sample(c("gene", "protein"), ncol(measurements) - 20, replace = TRUE))
#' # mcols(measurements)$feature <- colnames(measurements)
#' 
#' # We'll use different nFeatures for each assay. We'll also use repeated cross-validation with 5 repeats for speed in the example.
#' # set.seed(51773)
#' #result <- crossValidate(measurements, classes, nFeatures = c(clinical = 5, gene = 20, protein = 30), classifier = "randomForest", nRepeats = 5)
#' # performancePlot(result)
#' 
#' # Merge different assays. But we will only do this for two combinations. If assayCombinations is not specified it would attempt all combinations.
#' # set.seed(51773)
#' # resultMerge <- crossValidate(measurements, classes, assayCombinations = list(c("clinical", "protein"), c("clinical", "gene")), multiViewMethod = "merge", nRepeats = 5)
#' # performancePlot(resultMerge)
#' 
#' 
#' # performancePlot(c(result, resultMerge))
#' 
#' @importFrom survival Surv
#' @usage NULL
setGeneric("crossValidate", function(measurements, outcome, ...)
    standardGeneric("crossValidate"))

#' @rdname crossValidate
#' @export
setMethod("crossValidate", "DataFrame",
          function(measurements,
                   outcome,
                   nFeatures = 20,
                   selectionMethod = "auto",
                   classifier = "auto",
                   multiViewMethod = "none",
                   assayCombinations = "all",
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL, extraParams = NULL, verbose = 0)

          {
              # Check that data is in the right format, if not already done for MultiAssayExperiment input.
              if(!"assay" %in% colnames(S4Vectors::mcols(measurements))) # Assay is put there by prepareData for MultiAssayExperiment, skip if present. 
              {
                prepParams <- list(measurements, outcome)
                if("prepare" %in% names(extraParams))
                  prepParams <- c(prepParams, extraParams[["prepare"]])
                measurementsAndOutcome <- do.call(prepareData, prepParams)
                measurements <- measurementsAndOutcome[["measurements"]]
                outcome <- measurementsAndOutcome[["outcome"]]
              }
              
              # Ensure performance type is one of the ones that can be calculated by the package.
              isTuneCross <- !is.null(extraParams[["tuneCross"]])
              if(isTuneCross && !extraParams[["tuneCross"]][["performanceType"]] %in% c("auto", .ClassifyRenvir[["performanceTypes"]]))
                stop(paste("performanceType for parameter tuning must be one of", paste(c("auto", .ClassifyRenvir[["performanceTypes"]]), collapse = ", "), "but is", extraParams[["tuneCross"]][["performanceType"]]))
              
              isCategorical <- is.character(outcome) && (length(outcome) == 1 || length(outcome) == nrow(measurements)) || is.factor(outcome)
              if(isTuneCross && extraParams[["tuneCross"]][["performanceType"]] == "auto")
                if(isCategorical) extraParams[["tuneCross"]][["performanceType"]] <- "Balanced Accuracy" else extraParams[["tuneCross"]][["performanceType"]] <- "C-index"
              
              if(length(selectionMethod) == 1 && selectionMethod == "auto")
                if(isCategorical) selectionMethod <- "t-test" else selectionMethod <- "CoxPH"
              if(length(classifier) == 1 && classifier == "auto")
                if(isCategorical) classifier <- "randomForest" else classifier <- "CoxPH"
              
              
              # Which data-types or data-views are present?
              assayIDs <- unique(S4Vectors::mcols(measurements)$assay)
              if(is.null(assayIDs)) assayIDs <- 1

              # Check that other variables are in the right format and fix
              nFeaturesUse <- extraParams$select$tuneParams$nFeatures
              if(is.null(nFeaturesUse)) nFeaturesUse <- nFeatures
              nFeaturesUse <- cleanNFeatures(nFeatures = nFeaturesUse,
                                          measurements = measurements)
              selectionMethod <- cleanSelectionMethod(selectionMethod = selectionMethod,
                                                      measurements = measurements)
              classifier <- cleanClassifier(classifier = classifier,
                                            measurements = measurements, nFeatures = nFeaturesUse)
              
              ##!!!!! Do something with data combinations

              ################################
              #### No multiview
              ################################

              if(multiViewMethod == "none"){

                  # The below loops over assay and classifier and allows us to answer
                  # the following questions:
                  #
                  # 1) One assay using one classifier
                  # 2) One assay using multi classifiers
                  # 3) Multi assays individually
                  
                  # We should probably transition this to use grid instead.
                  resClassifier <-
                      sapply(assayIDs, function(assayIndex) {
                          # Loop over assays
                          sapply(classifier[[assayIndex]], function(classifierForAssay) {
                              # Loop over classifiers
                              sapply(selectionMethod[[assayIndex]], function(selectionForAssay) {
                                  # Loop over selectors
                                  measurementsUse <- measurements
                                  if(verbose > 0)
                                    message(Sys.time(), ": Running selection ", selectionForAssay, ", classifier ", classifierForAssay, '.')

                                  if(assayIndex != 1) measurementsUse <- measurements[, S4Vectors::mcols(measurements)[, "assay"] == assayIndex, drop = FALSE]
                                  CV(
                                      measurements = measurementsUse, outcome = outcome,
                                      assayIDs = assayIndex,
                                      nFeatures = nFeaturesUse[assayIndex],
                                      selectionMethod = selectionForAssay,
                                      classifier = classifierForAssay,
                                      multiViewMethod = multiViewMethod,
                                      nFolds = nFolds,
                                      nRepeats = nRepeats,
                                      nCores = nCores,
                                      characteristicsLabel = characteristicsLabel,
                                      extraParams = extraParams, verbose = verbose
                                  )
                              },
                              simplify = FALSE)
                          },
                          simplify = FALSE)
                      },
                      simplify = FALSE)
                  result <- unlist(unlist(resClassifier))
              }

              ################################
              #### Yes multiview
              ################################

              ### Merging or binding to combine data
              if(multiViewMethod == "merge"){


                  # The below loops over different combinations of assays and merges them together.
                  # This allows someone to answer which combinations of the assays might be most useful.

                  if(!is.list(assayCombinations) && assayCombinations[1] == "all") assayCombinations <- do.call("c", sapply(seq_along(assayIDs), function(nChoose) combn(assayIDs, nChoose, simplify = FALSE)))

                  result <- sapply(assayCombinations, function(assayIndex){
                      CV(measurements = measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assayIndex, drop = FALSE],
                         outcome = outcome, assayIDs = assayIndex,
                         nFeatures = nFeaturesUse[assayIndex],
                         selectionMethod = selectionMethod[assayIndex],
                         classifier = classifier[assayIndex],
                         multiViewMethod = ifelse(length(assayIndex) == 1, "none", multiViewMethod),
                         nFolds = nFolds,
                         nRepeats = nRepeats,
                         nCores = nCores,
                         characteristicsLabel = characteristicsLabel,
                         extraParams = extraParams, verbose = verbose)
                  }, simplify = FALSE)

              }


              ### Prevalidation to combine data
              if(multiViewMethod == "prevalidation"){


                  # The below loops over different combinations of assays and combines them together using prevalidation.
                  # This allows someone to answer which combinations of the assays might be most useful.


                  if(!is.list(assayCombinations) && assayCombinations[1] == "all")
                  {
                      assayCombinations <- do.call("c", sapply(seq_along(assayIDs), function(nChoose) combn(assayIDs, nChoose, simplify = FALSE)))
                      assayCombinations <- assayCombinations[sapply(assayCombinations, function(combination) "clinical" %in% combination, simplify = TRUE)]
                      if(length(assayCombinations) == 0) stop("No assayCombinations with \"clinical\" data")
                  }

                  result <- sapply(assayCombinations, function(assayIndex){
                      CV(measurements = measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assayIndex, drop = FALSE],
                         outcome = outcome, assayIDs = assayIndex,
                         nFeatures = nFeaturesUse[assayIndex],
                         selectionMethod = selectionMethod[assayIndex],
                         classifier = classifier[assayIndex],
                         multiViewMethod = ifelse(length(assayIndex) == 1, "none", multiViewMethod),
                         nFolds = nFolds,
                         nRepeats = nRepeats,
                         nCores = nCores,
                         characteristicsLabel = characteristicsLabel,
                         extraParams = extraParams, verbose = verbose)
                  }, simplify = FALSE)

              }



              ### Principal Components Analysis to combine data
              if(multiViewMethod == "PCA"){


                  # The below loops over different combinations of assays and combines them together using prevalidation.
                  # This allows someone to answer which combinations of the assays might be most useful.


                  if(!is.list(assayCombinations) && assayCombinations[1] == "all"){
                      assayCombinations <- do.call("c", sapply(seq_along(assayIDs),function(nChoose) combn(assayIDs, nChoose, simplify = FALSE)))
                      assayCombinations <- assayCombinations[sapply(assayCombinations, function(combination) "clinical" %in% combination, simplify = TRUE)]
                      if(length(assayCombinations) == 0) stop("No assayCombinations with \"clinical\" data")
                  }


                  result <- sapply(assayCombinations, function(assayIndex){
                      CV(measurements = measurements[, S4Vectors::mcols(measurements)$assay %in% assayIndex, drop = FALSE],
                         outcome = outcome, assayIDs = assayIndex,
                         nFeatures = nFeaturesUse[assayIndex],
                         selectionMethod = selectionMethod[assayIndex],
                         classifier = classifier[assayIndex],
                         multiViewMethod = ifelse(length(assayIndex) == 1, "none", multiViewMethod),
                         nFolds = nFolds,
                         nRepeats = nRepeats,
                         nCores = nCores,
                         characteristicsLabel = characteristicsLabel,
                         extraParams = extraParams, verbose = verbose)
                  }, simplify = FALSE)

              }
              if(length(result) == 1) result <- result[[1]]
              result

          })


#' @rdname crossValidate
#' @export
# One or more omics data sets, possibly with clinical data.
setMethod("crossValidate", "MultiAssayExperimentOrList",
          function(measurements,
                   outcome,
                   nFeatures = 20,
                   selectionMethod = "auto",
                   classifier = "auto",
                   multiViewMethod = "none",
                   assayCombinations = "all",
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL, extraParams = NULL, verbose = 0)
          {
              # Check that data is in the right format, if not already done for MultiAssayExperiment input.
              prepParams <- list(measurements, outcome)
              if("prepare" %in% names(extraParams))
                prepParams <- c(prepParams, extraParams[["prepare"]])
              measurementsAndOutcome <- do.call(prepareData, prepParams)

              crossValidate(measurements = measurementsAndOutcome[["measurements"]],
                            outcome = measurementsAndOutcome[["outcome"]], 
                            nFeatures = nFeatures,
                            classifier = classifier,
                            multiViewMethod = multiViewMethod,
                            assayCombinations = assayCombinations,
                            nFolds = nFolds,
                            nRepeats = nRepeats,
                            nCores = nCores,
                            characteristicsLabel = characteristicsLabel,
                            extraParams = extraParams,
                            verbose = verbose)
          })

#' @rdname crossValidate
#' @export
setMethod("crossValidate", "data.frame", # data.frame of numeric measurements.
          function(measurements,
                   outcome, 
                   nFeatures = 20,
                   selectionMethod = "auto",
                   classifier = "auto",
                   multiViewMethod = "none",
                   assayCombinations = "all",
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL, extraParams = NULL, verbose = 0)
          {
              measurements <- S4Vectors::DataFrame(measurements, check.names = FALSE)
              crossValidate(measurements = measurements,
                            outcome = outcome,
                            nFeatures = nFeatures,
                            selectionMethod = selectionMethod,
                            classifier = classifier,
                            multiViewMethod = multiViewMethod,
                            assayCombinations = assayCombinations,
                            nFolds = nFolds,
                            nRepeats = nRepeats,
                            nCores = nCores,
                            characteristicsLabel = characteristicsLabel, extraParams = extraParams, verbose = verbose)
          })

#' @rdname crossValidate
#' @export
setMethod("crossValidate", "matrix", # Matrix of numeric measurements.
          function(measurements,
                   outcome,
                   nFeatures = 20,
                   selectionMethod = "auto",
                   classifier = "auto",
                   multiViewMethod = "none",
                   assayCombinations = "all",
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL, extraParams = NULL, verbose = 0)
          {
              measurements <- S4Vectors::DataFrame(measurements, check.names = FALSE)
              crossValidate(measurements = measurements,
                            outcome = outcome,
                            nFeatures = nFeatures,
                            selectionMethod = selectionMethod,
                            classifier = classifier,
                            multiViewMethod = multiViewMethod,
                            assayCombinations = assayCombinations,
                            nFolds = nFolds,
                            nRepeats = nRepeats,
                            nCores = nCores,
                            characteristicsLabel = characteristicsLabel, extraParams = extraParams, verbose = verbose)
          })

######################################
######################################
cleanNFeatures <- function(nFeatures, measurements){
    #### Clean up
    if(!is.null(S4Vectors::mcols(measurements)$assay))
      obsFeatures <- unlist(as.list(table(S4Vectors::mcols(measurements)[, "assay"])))
    else obsFeatures <- ncol(measurements)
    if(is.null(nFeatures) || length(nFeatures) == 1 && nFeatures == "all") nFeatures <- as.list(obsFeatures)
    if(is.null(names(nFeatures)) && length(nFeatures) == 1) nFeatures <- as.list(pmin(obsFeatures, nFeatures))
    if(is.null(names(nFeatures)) && length(nFeatures) > 1) nFeatures <- sapply(obsFeatures, function(x)pmin(x, nFeatures), simplify = FALSE)
    #if(is.null(names(nFeatures)) && length(nFeatures) > 1) stop("nFeatures needs to be a named numeric vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && !all(names(obsFeatures) %in% names(nFeatures))) stop("nFeatures needs to be a named numeric vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && all(names(obsFeatures) %in% names(nFeatures)) & is(nFeatures, "numeric")) nFeatures <- as.list(pmin(obsFeatures, nFeatures[names(obsFeatures)]))
    if(!is.null(names(obsFeatures)) && all(names(obsFeatures) %in% names(nFeatures)) & is(nFeatures, "list")) nFeatures <- mapply(pmin, nFeatures[names(obsFeatures)], obsFeatures, SIMPLIFY = FALSE)
    nFeatures
}

######################################
######################################
cleanSelectionMethod <- function(selectionMethod, measurements){
    #### Clean up
    if(!is.null(S4Vectors::mcols(measurements)$assay))
      obsFeatures <- unlist(as.list(table(S4Vectors::mcols(measurements)[, "assay"])))
    else return(list(selectionMethod))

    if(is.null(names(selectionMethod)) & length(selectionMethod) == 1 & !is.null(names(obsFeatures))) selectionMethod <- sapply(names(obsFeatures), function(x) selectionMethod, simplify = FALSE)
    if(is.null(names(selectionMethod)) & length(selectionMethod) > 1 & !is.null(names(obsFeatures))) selectionMethod <- sapply(names(obsFeatures), function(x) selectionMethod, simplify = FALSE)
    #if(is.null(names(selectionMethod)) & length(selectionMethod) > 1) stop("selectionMethod needs to be a named character vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && !all(names(obsFeatures) %in% names(selectionMethod))) stop("selectionMethod needs to be a named character vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && all(names(obsFeatures) %in% names(selectionMethod)) & is(selectionMethod, "character")) selectionMethod <- as.list(selectionMethod[names(obsFeatures)])
    selectionMethod
}

######################################
######################################
cleanClassifier <- function(classifier, measurements, nFeatures){
    #### Clean up
    if(!is.null(S4Vectors::mcols(measurements)$assay))
      obsFeatures <- unlist(as.list(table(S4Vectors::mcols(measurements)[, "assay"])))
    else return(list(classifier))

    if(is.null(names(classifier)) & length(classifier) == 1 & !is.null(names(obsFeatures))) classifier <- sapply(names(obsFeatures), function(x)classifier, simplify = FALSE)
    if(is.null(names(classifier)) & length(classifier) > 1 & !is.null(names(obsFeatures))) classifier <- sapply(names(obsFeatures), function(x)classifier, simplify = FALSE)
    #if(is.null(names(classifier)) & length(classifier) > 1) stop("classifier needs to be a named character vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && !all(names(obsFeatures) %in% names(classifier))) stop("classifier needs to be a named character vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && all(names(obsFeatures) %in% names(classifier)) & is(classifier, "character")) classifier <- as.list(classifier[names(obsFeatures)])
    
    nFeatures <- nFeatures[names(classifier)]
    checkENs <- which(classifier %in% c("ridgeGLM", "elasticNetGLM", "LASSOGLM"))
    if(length(checkENs) > 0)
    {
      replacements <- sapply(checkENs, function(checkEN) ifelse(any(nFeatures[[checkEN]] == 1), "GLM", classifier[[checkEN]]))
      classifier[checkENs] <- replacements
      if(any(replacements == "GLM"))
        warning("Penalised GLM requires two or more features as input but there is only one.
Using an ordinary GLM instead.")
    }
    classifier
}

generateCrossValParams <- function(nRepeats, nFolds, nCores, extraParams){

    if(!exists(".Random.seed")) stop("Predictive modelling should always be reproducible. Please use set.seed(<number>) yourself and run 'crossValidate' again.")
    index <- ifelse(.Random.seed[2] + 2 == length(.Random.seed), 3, 3 + .Random.seed[2]) # Right after set.seed, the second number is the length of the random integer vector.
    seed <- .Random.seed[index] # Get current random number.
    
    # Set the BPparam RNGseed.
    if(nCores == 1)
    {
        BPparam <- SerialParam(RNGseed = seed)
    } else { # Parallel processing is desired.
        
        if(Sys.info()["sysname"] == "Windows") {# Only SnowParam suits Windows.
            BPparam <- BiocParallel::SnowParam(min(nCores, BiocParallel::snowWorkers("SOCK")), RNGseed = seed)
        } else if (Sys.info()["sysname"] %in% c("MacOS", "Linux")) {
            BPparam <- BiocParallel::MulticoreParam(min(nCores, BiocParallel::multicoreWorkers()), RNGseed = seed) # Multicore is faster than SNOW, but it doesn't work on Windows.
        } else { # Something weird.
            BPparam <- BiocParallel::bpparam() # BiocParallel will figure it out.
        }
    }
    tuneMode <- "none"
    performanceType <- "N/A"
    if(!is.null(extraParams[["tuneCross"]][["performanceType"]])) performanceType <- extraParams[["tuneCross"]][["performanceType"]]
    
    if(!is.null(extraParams[["tuneCross"]])) tuneMode <- extraParams[["tuneCross"]][["tuneMode"]]
    if(!any(tuneMode %in% c("Resubstitution", "Nested CV", "none"))) stop("tuneMode must be Nested CV or Resubstitution or none.")
    CrossValParams(permutations = nRepeats, folds = nFolds, parallelParams = BPparam, tuneMode = tuneMode, performanceType = performanceType)
}

# Returns a single parameter set.
generateModellingParams <- function(assayIDs,
                                    measurements,
                                    nFeatures,
                                    selectionMethod,
                                    classifier,
                                    multiViewMethod = "none",
                                    extraParams
){
    if(multiViewMethod != "none") {
        params <- generateMultiviewParams(assayIDs,
                                          measurements,
                                          nFeatures,
                                          selectionMethod,
                                          classifier,
                                          multiViewMethod, extraParams)
        return(params)
    }

    obsFeatures <- ncol(measurements)

    if(is.list(nFeatures) && any(nFeatures[[1]] > obsFeatures)) {
      warning("nFeatures greater than the maximum number of features in data. Setting to maximum.")
          nFeatures[[1]][nFeatures[[1]] > obsFeatures] <- obsFeatures
    }
    if(is.numeric(nFeatures) && any(nFeatures > obsFeatures)) {
      warning("nFeatures greater than the maximum number of features in data. Setting to maximum.")
          nFeatures[nFeatures > obsFeatures] <- obsFeatures
    }

    classifier <- unlist(classifier)
    
    # Check classifier
    knownClassifiers <- .ClassifyRenvir[["classifyKeywords"]][, "classifier Keyword"]
    if(any(!classifier %in% knownClassifiers))
        stop(paste("classifier must exactly match these options (be careful of case):", paste(knownClassifiers, collapse = ", ")))
    
    # Always return a list for ease of processing. Unbox at end if just one.
    classifierParams <- .classifierKeywordToParams(classifier, extraParams[["train"]][["tuneParams"]])

    if(!is.null(extraParams) && "train" %in% names(extraParams))
    {
      for(paramIndex in seq_along(extraParams[["train"]]))
      {
        parameter <- extraParams[["train"]][[paramIndex]]
        parameterName <- names(extraParams[["train"]])[paramIndex]
        if(length(parameter) == 1)
        {
          if(is.null(classifierParams$trainParams@otherParams)) classifierParams$trainParams@otherParams <- extraParams[["train"]][paramIndex]
          else classifierParams$trainParams@otherParams[parameterName] <- parameter
          if(parameterName %in% names(classifierParams$trainParams@tuneParams)) classifierParams$trainParams@tuneParams[[parameterName]] <- NULL
        } else if(length(parameter) > 1) {
          if(is.null(classifierParams$trainParams@tuneParams)) classifierParams$trainParams@tuneParams <- extraParams[["train"]][paramIndex]
          else classifierParams$trainParams@tuneParams[parameterName] <- parameter # Multiple values, so tune them.
        } else { # Remove the parameter
          inOther <- match(parameterName, names(classifierParams$trainParams@otherParams))
          if(!is.na(inOther)) classifierParams$trainParams@otherParams <- classifierParams$trainParams@otherParams[-inOther]
        } 
      }
    }
    if(!is.null(extraParams) && "predict" %in% names(extraParams))
    {
      for(paramIndex in seq_along(extraParams[["predict"]]))
      {
        parameter <- extraParams[["predict"]][[paramIndex]]
        parameterName <- names(extraParams[["predict"]])[paramIndex]
        if(length(parameter) == 1)
        {
          if(is.null(classifierParams$predictParams@otherParams)) classifierParams$predictParams@otherParams <- extraParams[["predict"]][paramIndex]
          else classifierParams$predictParams@otherParams[parameterName] <- parameter
          if(parameterName %in% names(classifierParams$predictParams@tuneParams)) classifierParams$predictParams@tuneParams[[parameterName]] <- NULL
        } else if(length(parameter) > 1) {
          if(is.null(classifierParams$predictParams@tuneParams)) classifierParams$predictParams@tuneParams <- extraParams[["predict"]][paramIndex]
          else classifierParams$predictParams@tuneParams[parameterName] <- parameter # Multiple values, so tune them.
        } else { # Remove the parameter
          inOther <- match(parameterName, names(classifierParams$predictParams@otherParams))
          if(!is.na(inOther)) classifierParams$predictParams@otherParams <- classifierParams$predictParams@otherParams[-inOther]
        } 
      }
    }    

    selectionMethod <- unlist(selectionMethod)

    if(selectionMethod != "none")
    {
      if(length(nFeatures[[1]]) > 1 || length(nFeatures) > 1)
      {
        extraParams[["select"]][["tuneParams"]][["nFeatures"]] <- unname(unlist(nFeatures))
        selectParams <- SelectParams(selectionMethod, nFeatures = NULL)
      } else {
          nFeatures <- unlist(nFeatures)
          selectParams <- SelectParams(selectionMethod, nFeatures = nFeatures)}
      
      if(!is.null(extraParams) && "select" %in% names(extraParams))
      {
        others <- setdiff(names(extraParams[["select"]]), "tuneParams")
        if(length(others) > 0) selectParams@otherParams <- extraParams[["select"]][others]
        selectParams@tuneParams <- extraParams[["select"]][["tuneParams"]]
      }
    } else {selectParams <- NULL}
    params <- ModellingParams(
        balancing = "none",
        selectParams = selectParams,
        trainParams = classifierParams$trainParams,
        predictParams = classifierParams$predictParams
    )

    #if(multiViewMethod != "none") stop("I haven't implemented multiview yet.")

    #
    # if(multiViewMethod == "prevalidation"){
    #     params$trainParams <- function(measurements, outcome) prevalTrainInterface(measurements, outcome, params)
    #     params$trainParams <- function(measurements, outcome) prevalTrainInterface(measurements, outcome, params)
    # }
    #

    params

}
######################################



generateMultiviewParams <- function(assayIDs,
                                    measurements,
                                    nFeatures,
                                    selectionMethod,
                                    classifier,
                                    multiViewMethod, extraParams){

    if(multiViewMethod == "merge"){

        if(length(classifier) > 1) classifier <- classifier[[1]]

        # Split measurements up by assay.
        assayTrain <- sapply(assayIDs, function(assayID) if(assayID == 1) measurements else measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assayID, drop = FALSE], simplify = FALSE)

        # Generate params for each assay. This could be extended to have different selectionMethods for each type
        paramsAssays <- mapply(generateModellingParams,
                                 nFeatures = nFeatures[assayIDs],
                                 selectionMethod = selectionMethod[assayIDs],
                                 assayIDs = assayIDs,
                                 measurements = assayTrain[assayIDs],
                                 MoreArgs = list(
                                     classifier = classifier,
                                     multiViewMethod = "none",
                                     extraParams = extraParams),
                                 SIMPLIFY = FALSE)

        # Generate some params for merged model. Which ones?
        # Reconsider how to do this well later. 
        params <- generateModellingParams(assayIDs = assayIDs,
                                          measurements = measurements,
                                          nFeatures = max(unlist(nFeatures)),
                                          selectionMethod = selectionMethod[[1]],
                                          classifier = classifier[[1]],
                                          multiViewMethod = "none",
                                          extraParams = extraParams)

        # Update selectParams to use
        params@selectParams <- SelectParams("selectMulti",
                                            params = paramsAssays,
                                            characteristics = S4Vectors::DataFrame(characteristic = "Selection Name", value = "merge"),
                                            tuneParams = NULL)
        return(params)
    }

    if(multiViewMethod == "prevalidation"){

        # Split measurements up by assay.
        assayTrain <- sapply(assayIDs, function(assayID) measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assayID, drop = FALSE], simplify = FALSE)

        # Generate params for each assay. This could be extended to have different selectionMethods for each type
        paramsAssays <- mapply(generateModellingParams,
                                 nFeatures = nFeatures[assayIDs],
                                 selectionMethod = selectionMethod[assayIDs],
                                 assayIDs = assayIDs,
                                 measurements = assayTrain[assayIDs],
                                 classifier = classifier[assayIDs],
                                 MoreArgs = list(multiViewMethod = "none", extraParams = extraParams),
                                 SIMPLIFY = FALSE)


        params <- ModellingParams(
            balancing = "none",
            selectParams = NULL,
            trainParams = TrainParams(prevalTrainInterface, params = paramsAssays, characteristics = paramsAssays$clinical@trainParams@characteristics,
                          getFeatures = prevalFeatures),
            predictParams = PredictParams(prevalPredictInterface, characteristics = paramsAssays$clinical@predictParams@characteristics)
        )

        return(params)
    }

    if(multiViewMethod == "PCA"){

        # Split measurements up by assay.
        assayTrain <- sapply(assayIDs, function(assayID) measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assayID, drop = FALSE], simplify = FALSE)

        # Generate params for each assay. This could be extended to have different selectionMethods for each type
        paramsClinical <-  list(clinical = generateModellingParams(
                                 nFeatures = nFeatures["clinical"],
                                 selectionMethod = selectionMethod["clinical"],
                                 assayIDs = "clinical",
                                 measurements = assayTrain[["clinical"]],
                                 classifier = classifier["clinical"],
                                 multiViewMethod = "none",
                                 extraParams = extraParams))


        params <- ModellingParams(
            balancing = "none",
            selectParams = NULL,
            trainParams = TrainParams(pcaTrainInterface, params = paramsClinical, nFeatures = nFeatures, characteristics = paramsClinical$clinical@trainParams@characteristics,
                                      getFeatures = PCAfeatures),
            predictParams = PredictParams(pcaPredictInterface, characteristics = paramsClinical$clinical@predictParams@characteristics)
        )

        return(params)
    }

}

# measurements, outcome are mutually exclusive with x, outcomeTrain, measurementsTest, outcomeTest.
CV <- function(measurements, outcome, x, outcomeTrain, measurementsTest, outcomeTest,
               assayIDs,
               nFeatures,
               selectionMethod,
               classifier,
               multiViewMethod,
               nFolds,
               nRepeats,
               nCores,
               characteristicsLabel, extraParams, verbose)

{
    # Which data-types or data-views are present?
    if(is.null(characteristicsLabel)) characteristicsLabel <- "none"

    # Setup cross-validation parameters. Could be needed for independent train/test if parameter tuning
    # is specified to be done by nested cross-validation.
    crossValParams <- generateCrossValParams(nRepeats = nRepeats,
                                             nFolds = nFolds,
                                             nCores = nCores,
                                             extraParams = extraParams)
    

    # Turn text into TrainParams and TestParams objects
    modellingParams <- generateModellingParams(assayIDs = assayIDs,
                                               measurements = if(!is.null(measurements)) measurements else x,
                                               nFeatures = nFeatures,
                                               selectionMethod = selectionMethod,
                                               classifier = classifier,
                                               multiViewMethod = multiViewMethod, extraParams = extraParams)
    
    if(length(assayIDs) > 1 || length(assayIDs) == 1 && assayIDs != 1) assayText <- assayIDs else assayText <- NULL
    characteristics <- S4Vectors::DataFrame(characteristic = c(if(!is.null(assayText)) "Assay Name" else NULL, "Classifier Name", "Selection Name", "multiViewMethod", "characteristicsLabel"), value = c(if(!is.null(assayText)) paste(assayText, collapse = ", ") else NULL, paste(classifier, collapse = ", "),  paste(selectionMethod, collapse = ", "), multiViewMethod, characteristicsLabel))

    if(!is.null(measurements))
    { # Cross-validation.
      classifyResults <- runTests(measurements, outcome, crossValParams = crossValParams, modellingParams = modellingParams, characteristics = characteristics, verbose = verbose)
    } else { # Independent training and testing.
      classifyResults <- runTest(x, outcomeTrain, measurementsTest, outcomeTest, crossValParams = crossValParams, modellingParams = modellingParams, characteristics = characteristics)
      if(is.character(classifyResults)) stop(classifyResults)
      fullResult <- runTest(measurements, outcome, measurements, outcome, crossValParams = crossValParams, modellingParams = modellingParams, characteristics = characteristics, .iteration = 1)
      classifyResults@finalModel <- fullResult$models
      class(classifyResults@finalModel) <- c("trainedByClassifyR", classifyResults@finalModel)
      attr(classifyResults@finalModel, "predictFunction") <- modellingParams@trainParams@classifier
    }
    
    classifyResults
}

simplifyResults <- function(results, values = c("assay", "classifier", "selectionMethod", "multiViewMethod")){
    ch <- sapply(results, function(x) x@characteristics[x@characteristics$characteristic %in% values, "value"], simplify = TRUE)
    ch <- data.frame(t(ch))
    results[!duplicated(ch)]
}

#' @rdname crossValidate
#' @importFrom generics train
#' @method train matrix
#' @export
train.matrix <- function(x, outcomeTrain, ...)
               {
                 x <- DataFrame(x, check.names = FALSE)
                 train(x, outcomeTrain, ...)
               }

#' @rdname crossValidate
#' @method train data.frame
#' @export
train.data.frame <- function(x, outcomeTrain, ...)
                    {
                      x <- DataFrame(x, check.names = FALSE)
                      train(x, outcomeTrain, ...)
                    }

#' @rdname crossValidate
#' @param assayIDs A character vector for assays to train with. Special value \code{"all"}
#' uses all assays in the input object.
#' @method train DataFrame
#' @export
train.DataFrame <- function(x, outcomeTrain, selectionMethod = "auto", nFeatures = 20, classifier = "auto",
                            multiViewMethod = "none", assayIDs = "all", extraParams = NULL, verbose = 0, ...)
                   {
              prepParams <- list(x, outcomeTrain)
              if(!is.null(extraParams) && "prepare" %in% names(extraParams))
                prepParams <- c(prepParams, extraParams[["prepare"]])
              measurementsAndOutcome <- do.call(prepareData, prepParams)
              measurements <- measurementsAndOutcome[["measurements"]]
              outcomeTrain <- measurementsAndOutcome[["outcome"]]
              
              # Ensure performance type is one of the ones that can be calculated by the package.
              isTuneCross <- !is.null(extraParams[["tuneCross"]])
              if(isTuneCross && !extraParams[["tuneCross"]][["performanceType"]] %in% c("auto", .ClassifyRenvir[["performanceTypes"]]))
                stop(paste("performanceType for tuning must be one of", paste(c("auto", .ClassifyRenvir[["performanceTypes"]]), collapse = ", "), "but is", extraParams[["tuneCross"]][["performanceType"]]))
              
              isCategorical <- is.character(outcomeTrain) && (length(outcomeTrain) == 1 || length(outcomeTrain) == nrow(measurements)) || is.factor(outcomeTrain)
              if(isTuneCross && extraParams[["tuneCross"]][["performanceType"]] == "auto")
                if(isCategorical) extraParams[["tuneCross"]][["performanceType"]] <- "Balanced Accuracy" else extraParams[["tuneCross"]][["performanceType"]] <- "C-index"
              if(length(selectionMethod) == 1 && selectionMethod == "auto")
                if(isCategorical) selectionMethod <- "t-test" else selectionMethod <- "CoxPH"
              if(length(classifier) == 1 && classifier == "auto")
                if(isCategorical) classifier <- "randomForest" else classifier <- "CoxPH"

              classifier <- cleanClassifier(classifier = classifier, measurements = measurements)
              selectionMethod <- cleanSelectionMethod(selectionMethod = selectionMethod, measurements = measurements)
              if(assayIDs == "all") assayIDs <- unique(S4Vectors::mcols(measurements)$assay)
              if(is.null(assayIDs)) assayIDs <- 1
              names(assayIDs) <- assayIDs
              names(classifier) <- assayIDs

              if(multiViewMethod == "none"){
                  resClassifier <-
                      sapply(assayIDs, function(assayIndex) {
                          # Loop over assays
                          sapply(classifier[[assayIndex]], function(classifierForAssay) {
                              # Loop over classifiers
                                sapply(selectionMethod[[assayIndex]], function(selectionForAssay) {
                                  # Loop over selectors
                              
                                  measurementsUse <- measurements
                                  if(assayIndex != 1) measurementsUse <- measurements[, S4Vectors::mcols(measurements)[, "assay"] == assayIndex, drop = FALSE]
                                  
                                  modellingParams <- generateModellingParams(assayIDs = assayIDs, measurements = measurements, nFeatures = nFeatures,
                                                     selectionMethod = selectionMethod, classifier = classifier, multiViewMethod = "none", extraParams = extraParams)

                                  if(!is.null(modellingParams@selectParams))
                                  {
                                    topFeatures <- .doSelection(measurementsUse, outcomeTrain, CrossValParams(), modellingParams, verbose = verbose)
                                    selectedFeaturesIndices <- topFeatures[[2]] # Extract for subsetting.
                                    tuneDetailsSelect <- topFeatures[[3]]
                                    measurementsUse <- measurementsUse[, selectedFeaturesIndices]
                                  } else {
                                    tuneDetailsSelect <- NULL
                                    measurementsUse <- measurements
                                  }

                                  classifierParams <- .classifierKeywordToParams(classifierForAssay, extraParams[["train"]][["tuneParams"]])
                                  if(!is.null(extraParams) && "train" %in% names(extraParams))
                                  {
                                     for(paramIndex in seq_along(extraParams[["train"]]))
                                     {
                                        parameter <- extraParams[["train"]][[paramIndex]]
                                        parameterName <- names(extraParams[["train"]])[paramIndex]
                                        if(length(parameter) == 1)
                                        {
                                          if(is.null(classifierParams$trainParams@otherParams)) classifierParams$trainParams@otherParams <- extraParams[["train"]][paramIndex]
                                          else classifierParams$trainParams@otherParams[parameterName] <- parameter
                                        } else if (length(parameter) > 1) {
                                          if(is.null(classifierParams$trainParams@tuneParams)) classifierParams$trainParams@tuneParams <- extraParams[["train"]][paramIndex]
                                          else classifierParams$trainParams@tuneParams[parameterName] <- parameter # Multiple values, so tune them.
                                        } else { # Remove the parameter
                                          inOther <- match(parameterName, names(classifierParams$trainParams@otherParams))
                                          inTune <- match(parameterName, names(classifierParams$trainParams@tuneParams))
                                          if(!is.na(inOther)) classifierParams$trainParams@otherParams <- classifierParams$trainParams@otherParams[-inOther]
                                          if(!is.na(inTune)) classifierParams$trainParams@tuneParams <- classifierParams$trainParams@otherParams[-inTune]
                                        }
                                      }
                                    }
                                  if(!is.null(extraParams) && "predict" %in% names(extraParams))
                                  {
                                      for(paramIndex in seq_along(extraParams[["predict"]]))
                                      {
                                        parameter <- extraParams[["predict"]][[paramIndex]]
                                        parameterName <- names(extraParams[["predict"]])[paramIndex]
                                        if(length(parameter) == 1)
                                        {
                                          if(is.null(classifierParams$predictParams@otherParams)) classifierParams$predictParams@otherParams <- extraParams[["predict"]][paramIndex]
                                          else classifierParams$predictParams@otherParams[parameterName] <- parameter
                                        } else if (length(parameter) > 1) {
                                          if(is.null(classifierParams$predictParams@tuneParams)) classifierParams$predictParams@tuneParams <- extraParams[["predict"]][paramIndex]
                                          else classifierParams$predictParams@tuneParams[parameterName] <- parameter # Multiple values, so tune them.
                                        } else { # Remove the parameter
                                          inOther <- match(parameterName, names(classifierParams$predictParams@otherParams))
                                          inTune <- match(parameterName, names(classifierParams$predictParams@tuneParams))
                                          if(!is.na(inOther)) classifierParams$predictParams@otherParams <- classifierParams$predictParams@otherParams[-inOther]
                                          if(!is.na(inTune)) classifierParams$predictParams@tuneParams <- classifierParams$predictParams@otherParams[-inTune]
                                        } 
                                      }
                                    }
                                  
                                  modellingParams <- ModellingParams(balancing = "none", selectParams = NULL,
                                                                     trainParams = classifierParams$trainParams, predictParams = classifierParams$predictParams)
                                  if(!is.null(tuneDetailsSelect))
                                  {
                                    tuneDetailsSelectUse <- tuneDetailsSelect[["tuneCombinations"]][tuneDetailsSelect[["bestIndex"]], , drop = FALSE]
                                    avoidTune <- match(colnames(tuneDetailsSelectUse), names(modellingParams@trainParams@tuneParams))
                                    if(any(!is.na(avoidTune)))
                                    {
                                      modellingParams@trainParams@otherParams <- c(modellingParams@trainParams@otherParams, tuneDetailsSelectUse[!is.na(avoidTune)])
                                      modellingParams@trainParams@tuneParams <- modellingParams@trainParams@tuneParams[-na.omit(avoidTune)]
                                      if(length(modellingParams@trainParams@tuneParams) == 0) modellingParams@trainParams@tuneParams <- NULL
                                    }
                                  }
                                  
                                  trained <- .doTrain(measurementsUse, outcomeTrain, NULL, NULL, CrossValParams(), modellingParams, verbose = verbose)[["model"]]
                                  attr(trained, "predictFunction") <- classifierParams$predictParams@predictor
                                  trained
                                  ## train model
                                }, simplify = FALSE)
                          }, simplify = FALSE)
                      }, simplify = FALSE)

                  models <- unlist(unlist(resClassifier, recursive = FALSE), recursive = FALSE)
                  if(length(models) == 1) {
                      model <- models[[1]]
                      class(model) <- c("trainedByClassifyR", class(model))
                      models <- NULL
                  } else {
                      class(models) <- c("listOfModels", "trainedByClassifyR", class(models))
                  }
              }

              ################################
              #### Yes multiview
              ################################

              ### Merging or binding to combine data
              if(multiViewMethod == "merge"){
                  measurementsUse <- measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assayIDs, drop = FALSE]
                  model <- .doTrain(measurementsUse, outcomeTrain, NULL, NULL, crossValParams, modellingParams, verbose = verbose)[["model"]]
                  attr(model, "predictFunction") <- modellingParams@trainParams@classifier
                  class(model) <- c("trainedByClassifyR", class(model))
              }


              ### Prevalidation to combine data
              if(multiViewMethod == "prevalidation"){
                # Split measurements up by assay.
                 assayTrain <- sapply(assayIDs, function(assayID) measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assayID, drop = FALSE], simplify = FALSE)

               # Generate params for each assay. This could be extended to have different selectionMethods for each type
                 paramsAssays <- mapply(generateModellingParams,
                                        nFeatures = nFeatures[assayIDs],
                                        selectionMethod = selectionMethod[assayIDs],
                                        assayIDs = assayIDs,
                                        measurements = assayTrain[assayIDs],
                                        classifier = classifier[assayIDs],
                                        MoreArgs = list(multiViewMethod = "none"),
                                 SIMPLIFY = FALSE)

                 modellingParams <- ModellingParams(
                                    balancing = "none",
                                    selectParams = NULL,
                                    trainParams = TrainParams(prevalTrainInterface, params = paramsAssays, characteristics = paramsAssays$clinical@trainParams@characteristics,
                                                              getFeatures = prevalFeatures),
                                    predictParams = PredictParams(prevalPredictInterface, characteristics = paramsAssays$clinical@predictParams@characteristics))
                 model <- .doTrain(measurementsUse, outcomeTrain, NULL, NULL, crossValParams, modellingParams, verbose = verbose)[["model"]]
                 attr(model, "predictFunction") <- modellingParams@trainParams@classifier
                 class(model) <- c("trainedByClassifyR", class(model))
              }
              
              ### Principal Components Analysis to combine data
              if(multiViewMethod == "PCA"){
                measurementsUse <- measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assayIDs, drop = FALSE]
                paramsClinical <-  list(clinical = generateModellingParams(
                                        assayIDs = "clinical",
                                        measurements = measurements[, S4Vectors::mcols(measurements)[["assay"]] == "clinical", drop = FALSE],
                                        classifier = classifier["clinical"],
                                        multiViewMethod = "none"))
                
                modellingParams <- ModellingParams(balancing = "none", selectParams = NULL,
                                   trainParams = TrainParams(pcaTrainInterface, params = paramsClinical, nFeatures = nFeatures, characteristics = paramsClinical$clinical@trainParams@characteristics,
                                                             getFeatures = PCAfeatures),
                                   predictParams = PredictParams(pcaPredictInterface, characteristics = paramsClinical$clinical@predictParams@characteristics))
                model <- .doTrain(measurementsUse, outcomeTrain, NULL, NULL, crossValParams, modellingParams, verbose = verbose)[["model"]]
                attr(model, "predictFunction") <- modellingParams@trainParams@classifier
                class(model) <- c("trainedByClassifyR", class(model))
              }
              if(missing(models) || is.null(models)) return(model) else return(models)
          }

#' @rdname crossValidate
#' @method train list
#' @export
train.list <- function(x, outcomeTrain, ...)
              {
                # Check data type is valid
                if (!(all(sapply(x, function(element) is(element, "tabular")))))
                  stop("assays must be of type data.frame, DataFrame or matrix")
              
                # Check the list is named
                if (is.null(names(x)))
                  stop("Measurements must be a named list")
              
                # Check same number of samples for all datasets
                if (!length(unique(sapply(x, nrow))) == 1)
                  stop("All datasets must have the same samples")
              
                # Check the number of outcome is the same
                if (!all(sapply(x, nrow) == length(outcomeTrain)) && !is.character(outcomeTrain))
                  stop("outcome must have same number of samples as measurements")
              
              df_list <- sapply(x, S4Vectors::DataFrame)
              
              df_list <- mapply(function(meas, nam){
                  S4Vectors::mcols(meas)$assay <- nam
                  S4Vectors::mcols(meas)$feature <- colnames(meas)
                  meas
              }, df_list, names(df_list))
              
              combined_df <- do.call(cbind, df_list)
              
              # Each list of tabular data has been collapsed into a DataFrame.
              # Will be subset to relevant assayIDs inside the DataFrame method.
              
              train(combined_df, outcomeTrain, ...)
}

#' @rdname crossValidate
#' @method train MultiAssayExperiment
#' @export
train.MultiAssayExperiment <- function(x, outcome, ...)
          {
              prepArgs <- list(x, outcome)
              extraInputs <- list(...)
              prepExtras <- trainExtras <- numeric()
              if(length(extraInputs) > 0)
                prepExtras <- which(names(extraInputs) %in% .ClassifyRenvir[["prepareDataFormals"]])
              if(length(prepExtras) > 0)
                prepArgs <- append(prepArgs, extraInputs[prepExtras])
              measurementsAndOutcome <- do.call(prepareData, prepArgs)
              trainArgs <- list(measurementsAndOutcome[["measurements"]], measurementsAndOutcome[["outcome"]])
              if(length(extraInputs) > 0)
                trainExtras <- which(!names(extraInputs) %in% .ClassifyRenvir[["prepareDataFormals"]])
              if(length(trainExtras) > 0)
                trainArgs <- append(trainArgs, extraInputs[trainExtras])
              do.call(train, trainArgs)
          }

#' @rdname crossValidate
#' @param object A fitted model or a list of such models.
#' @param newData For the \code{predict} function, an object of type \code{matrix}, \code{data.frame}
#' \code{DataFrame}, \code{list} (of matrices or data frames) or \code{MultiAssayExperiment} containing
#' the data to make predictions with with either a fitted model created by \code{train} or the final model
#' stored in a \code{\link{ClassifyResult}} object.
#' @method predict trainedByClassifyR
#' @export
predict.trainedByClassifyR <- function(object, newData, outcome, ...)
{
  if(is(newData, "tabular")) # Simply tabular data.
  {
    colnames(newData) <- make.names(colnames(newData)) # Ensure that feature names are syntactically valid, like during model fitting.
  } else if(is.list(newData) && !is(object, "listOfModels")) # Don't check all those conditions that train function does.
  { # Merge the list of data tables and keep track of assay names in columns' metadata.
    newData <- mapply(function(meas, nam){
               S4Vectors::mcols(meas)$assay <- nam
               S4Vectors::mcols(meas)$feature <- colnames(meas)
               meas
               }, newData, names(newData))
    newData <- do.call(cbind, newData)
    } else if(is(newData, "MultiAssayExperiment"))
    {
              newData <- prepareData(newData, outcome)[["measurements"]]
    }
    
    predictFunctionUse <- attr(object, "predictFunction")
    class(object) <- class(object)[-1] # Now want the predict method of the specific model to be picked, so put model class first.
    if (is(object, "listOfModels")) 
         mapply(function(model, assay) predictFunctionUse(model, assay), object, newData, MoreArgs = list(...), SIMPLIFY = FALSE)
    else do.call(predictFunctionUse, list(object, newData, ...)) # Object is itself a trained model and it is assumed that a predict method is defined for it.
}
