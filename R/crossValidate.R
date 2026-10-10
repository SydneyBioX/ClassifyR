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
#' @param nFeatures The number of features to choose in the feature selection stage and use in the subsequent classifier training stage. If a named vector with the same names of multiple assays, 
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
#' Every assay, classifier, selection method and combination of assays in one call is evaluated on the same
#' training and test splits, so their performances can be compared sample by sample and split by split.
#'
#' If \code{nFeatures} has several values and no tuning mode is given in \code{extraParams}, the number of features is
#' chosen by resubstitution: the classifier is trained and evaluated on the training samples for each value. Classifiers
#' that fit their training samples perfectly (e.g. random forest) give every value the same performance, and then the
#' smallest value is chosen.
#'
#' The penalised GLM classifiers (\code{"ridgeGLM"}, \code{"elasticNetGLM"} and \code{"LASSOGLM"}) choose lambda by
#' 5-fold cross-validation of the balanced error within each training set. \code{extraParams = list(train =
#' list(lambdaTuning = "resubstitution"))} chooses it by the balanced error of the training samples instead, and
#' \code{nFoldsLambda} sets the number of folds.
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
              
              # Automatically set tuning mode if user has specified a range of nFeatures.
              if(length(nFeatures) > 1)
              {
                if(is.null(extraParams))
                {
                  message("Tune mode is \"none\" but 'nFeatures' has multiple values. Setting to resubstitution performance.")
                  extraParams <- list(tuneCross = list(tuneMode = "Resubstitution", performanceType = "auto"))
                } else if(!"tuneCross" %in% names(extraParams)) {
                  message("Tune mode is \"none\" but 'nFeatures' has multiple values. Setting to resubstitution performance.")
                  extraParams[["tuneCross"]] <- list(tuneMode = "Resubstitution", performanceType = "auto")
                }
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
              
              if(length(multiViewMethod) != 1 || !multiViewMethod %in% c("none", "merge", "prevalidation", "PCA"))
                stop("multiViewMethod must be one of \"none\", \"merge\", \"prevalidation\" or \"PCA\" (see available(\"multiViewMethod\")).")

              # Every cross-validation is prepared first, then all of their splits share one pool of workers.
              queue <- .crossValidationQueue()


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
                                      extraParams = extraParams, verbose = verbose, queue = queue
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

              # Merging, prevalidation or PCA combine the assays of each combination of assays. This allows
              # someone to answer which combinations of the assays might be most useful.
              if(multiViewMethod %in% c("merge", "prevalidation", "PCA"))
              {
                  if(!is.list(assayCombinations) && assayCombinations[1] == "all")
                  {
                      assayCombinations <- do.call("c", sapply(seq_along(assayIDs), function(nChoose) combn(assayIDs, nChoose, simplify = FALSE)))
                      if(multiViewMethod != "merge") # Prevalidation and PCA add the other assays to the clinical data.
                      {
                          assayCombinations <- assayCombinations[sapply(assayCombinations, function(combination) "clinical" %in% combination, simplify = TRUE)]
                          if(length(assayCombinations) == 0) stop("No assayCombinations with \"clinical\" data")
                      }
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
                         extraParams = extraParams, verbose = verbose, queue = queue)
                  }, simplify = FALSE)
              }

              result <- .runCrossValidationQueue(queue, result, nCores)
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
                            selectionMethod = selectionMethod,
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

generateCrossValParams <- function(nRepeats, nFolds, nCores, extraParams, seed = NULL){

    if(is.null(seed))
    {
      if(!exists(".Random.seed")) stop("Predictive modelling should always be reproducible. Please use set.seed(<number>) yourself and run 'crossValidate' again.")
      index <- ifelse(.Random.seed[2] + 2 == length(.Random.seed), 3, 3 + .Random.seed[2]) # Right after set.seed, the second number is the length of the random integer vector.
      seed <- .Random.seed[index] # Get current random number.
    }
    
    # The seed sets the random number streams of the splits. Workers come from one pool per crossValidate call
    # (.makeWorkerPool); nested cross-validations within a split run serially.
    BPparam <- BiocParallel::SerialParam(RNGseed = seed)
    tuneMode <- "none"
    performanceType <- "N/A"
    if(!is.null(extraParams[["tuneCross"]][["performanceType"]])) performanceType <- extraParams[["tuneCross"]][["performanceType"]]
    
    if(!is.null(extraParams[["tuneCross"]])) tuneMode <- extraParams[["tuneCross"]][["tuneMode"]]
    if(!any(tuneMode %in% c("Resubstitution", "Nested CV", "none"))) stop("tuneMode must be Nested CV or Resubstitution or none.")
    CrossValParams(permutations = nRepeats, folds = nFolds, parallelParams = BPparam, tuneMode = tuneMode, performanceType = performanceType)
}

# Applies the user's extraParams for one stage ("train" or "predict") to its TrainParams or PredictParams.
# A single value is a fixed setting, several values are tuned, and an empty value removes a setting.
# The element "tuneParams" holds a list of parameter ranges to tune ("auto" uses the classifier's preset ranges).
.applyStageExtras <- function(stageParams, extras)
{
  canTune <- methods::.hasSlot(stageParams, "tuneParams")
  for(parameterName in names(extras))
  {
    parameter <- extras[[parameterName]]
    if(parameterName == "tuneParams")
    {
      if(is.list(parameter))
      {
        if(!canTune) stop("Prediction parameters can't be tuned.")
        stageParams@tuneParams[names(parameter)] <- parameter
      }
    } else if(length(parameter) == 1) {
      stageParams@otherParams[parameterName] <- list(parameter)
      if(canTune) stageParams@tuneParams[[parameterName]] <- NULL
    } else if(length(parameter) > 1) {
      if(!canTune) stop("Prediction parameter '", parameterName, "' has more than one value but prediction parameters can't be tuned.")
      stageParams@tuneParams[parameterName] <- list(parameter)
    } else {
      stageParams@otherParams[[parameterName]] <- NULL
      if(canTune) stageParams@tuneParams[[parameterName]] <- NULL
    }
  }
  if(canTune && length(stageParams@tuneParams) == 0) stageParams@tuneParams <- NULL
  if(length(stageParams@otherParams) == 0) stageParams@otherParams <- NULL
  stageParams
}

# A queue of prepared cross-validations; add() returns the position of the added one.
# The first cross-validation's seed and splits are kept, so that every later one uses the same splits.
.crossValidationQueue <- function()
{
  items <- list()
  shared <- NULL
  list(add = function(item)
       {
         items[[length(items) + 1]] <<- item
         if(is.null(shared))
           shared <<- list(seed = BiocParallel::bpRNGseed(item[["crossValParams"]]@parallelParams),
                           splits = item[["splits"]], nSamples = length(item[["outcome"]]))
         length(items)
       },
       items = function() items,
       shared = function() shared)
}

# Parallel workers for nCores cores, made once per call by .poolApply after the data are prepared: forked processes
# on Linux and macOS, which share the main session's data without copying it, and socket workers on Windows, which
# are sent the data once each.
.makeWorkerPool <- function(nCores, nTasks)
{
  structure(list(workers = as.integer(max(1, min(nCores, nTasks))),
                 type = if(.Platform$OS.type == "windows") "PSOCK" else "FORK"), class = "workerPool")
}

# lapply(X, FUN) on the workers of pool. FUN is left in .ClassifyRenvir of each worker (by forking, or sent once to
# each socket worker), so only the elements of X and the results are passed between processes. The random number
# state of the main session is the same afterwards as before.
.poolApply <- function(X, FUN, pool)
{
  previousSeed <- if(exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if(is.null(previousSeed)) suppressWarnings(rm(".Random.seed", envir = globalenv())) else
            assign(".Random.seed", previousSeed, envir = globalenv()))
  if(pool[["workers"]] == 1) return(lapply(X, FUN))
  
  assign("currentTask", FUN, envir = .ClassifyRenvir)
  on.exit(rm("currentTask", envir = .ClassifyRenvir), add = TRUE)
  if(pool[["type"]] == "FORK")
  {
    workers <- parallel::makeForkCluster(pool[["workers"]])
  } else {
    workers <- parallel::makePSOCKcluster(pool[["workers"]])
    # The same package libraries as this session, then the task function. The functions sent have the global
    # environment as theirs, so that this call's data aren't sent along with them.
    setLibraries <- function(paths) { base::.libPaths(paths); NULL }
    setTask <- function(task) { assign("currentTask", task, envir = get(".ClassifyRenvir", envir = asNamespace("ClassifyR"))); NULL }
    environment(setLibraries) <- environment(setTask) <- globalenv()
    parallel::clusterCall(workers, setLibraries, .libPaths())
    parallel::clusterCall(workers, setTask, FUN)
  }
  on.exit(parallel::stopCluster(workers), add = TRUE)
  runTask <- function(element) get("currentTask", envir = .ClassifyRenvir)(element)
  environment(runTask) <- asNamespace("ClassifyR") # Sent to workers as a reference, not with this call's data.
  # Tasks are handed out in chunks as workers become free, about eight chunks per worker: enough for quick and slow
  # tasks to balance out, and few enough that workers don't wait long for the main process to send the next chunk.
  parallel::parLapplyLB(workers, X, runTask, chunk.size = max(1L, ceiling(length(X) / (8L * pool[["workers"]]))))
}

# Runs the splits of every queued cross-validation in one pool of workers and replaces each queue position in
# positions (a vector or list, possibly named) by its ClassifyResult.
.runCrossValidationQueue <- function(queue, positions, nCores)
{
  crossValidations <- queue$items()
  if(length(crossValidations) == 0) return(positions)
  nTasks <- sum(sapply(crossValidations, function(crossValidation) length(crossValidation[["splits"]][["train"]])))
  results <- .runTestsSplits(crossValidations, .makeWorkerPool(nCores, nTasks))
  classifyResults <- mapply(.assembleTests, crossValidations, results, SIMPLIFY = FALSE)
  lapply(as.list(positions), function(position) classifyResults[[position]])
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

    classifierParams$trainParams <- .applyStageExtras(classifierParams$trainParams, extraParams[["train"]])
    if(!is.null(classifierParams$predictParams))
      classifierParams$predictParams <- .applyStageExtras(classifierParams$predictParams, extraParams[["predict"]])

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

# Cross-validation of one assay or combination of assays with one classifier and selection method.
CV <- function(measurements, outcome,
               assayIDs,
               nFeatures,
               selectionMethod,
               classifier,
               multiViewMethod,
               nFolds,
               nRepeats,
               nCores,
               characteristicsLabel, extraParams, verbose, queue = NULL)

{
    # Which data-types or data-views are present?
    if(is.null(characteristicsLabel)) characteristicsLabel <- "none"

    # Cross-validations queued by one crossValidate call share the first one's seed and splits.
    shared <- if(!is.null(queue)) queue$shared() else NULL
    crossValParams <- generateCrossValParams(nRepeats = nRepeats,
                                             nFolds = nFolds,
                                             nCores = nCores,
                                             extraParams = extraParams,
                                             seed = shared[["seed"]])
    

    # Turn text into TrainParams and TestParams objects
    modellingParams <- generateModellingParams(assayIDs = assayIDs,
                                               measurements = measurements,
                                               nFeatures = nFeatures,
                                               selectionMethod = selectionMethod,
                                               classifier = classifier,
                                               multiViewMethod = multiViewMethod, extraParams = extraParams)
    if(crossValParams@tuneMode == "none" && !is.null(modellingParams@trainParams@tuneParams))
      stop("Training parameters ", paste(names(modellingParams@trainParams@tuneParams), collapse = ", "), " have several values to tune ",
           "but no tuning mode is set. Specify extraParams = list(tuneCross = list(tuneMode = ..., performanceType = ...)).")
    
    if(length(assayIDs) > 1 || length(assayIDs) == 1 && assayIDs != 1) assayText <- assayIDs else assayText <- NULL
    characteristics <- S4Vectors::DataFrame(characteristic = c(if(!is.null(assayText)) "Assay Name" else NULL, "Classifier Name", "Selection Name", "multiViewMethod", "characteristicsLabel"), value = c(if(!is.null(assayText)) paste(assayText, collapse = ", ") else NULL, paste(classifier, collapse = ", "),  paste(selectionMethod, collapse = ", "), multiViewMethod, characteristicsLabel))

    if(!is.null(queue)) # Prepare it now and run it with the others in the queue; return its position in the queue.
      return(queue$add(.prepareTests(measurements, outcome, crossValParams, modellingParams, characteristics, verbose,
                                     splits = if(!is.null(shared) && shared[["nSamples"]] == nrow(measurements)) shared[["splits"]],
                                     deferFinal = TRUE)))
    runTests(measurements, outcome, crossValParams = crossValParams, modellingParams = modellingParams, characteristics = characteristics, verbose = verbose)
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

              nFeatures <- cleanNFeatures(nFeatures = nFeatures, measurements = measurements)
              classifier <- cleanClassifier(classifier = classifier, measurements = measurements, nFeatures = nFeatures)
              selectionMethod <- cleanSelectionMethod(selectionMethod = selectionMethod, measurements = measurements)
              if(identical(assayIDs, "all")) assayIDs <- unique(S4Vectors::mcols(measurements)$assay)
              if(is.null(assayIDs)) assayIDs <- 1
              names(assayIDs) <- assayIDs

              # Parameter tuning, if requested, is done on the training samples alone.
              tuneCross <- extraParams[["tuneCross"]]
              crossValParams <- CrossValParams(parallelParams = BiocParallel::SerialParam(),
                                               tuneMode = if(is.null(tuneCross)) "none" else tuneCross[["tuneMode"]],
                                               performanceType = if(is.null(tuneCross)) "auto" else tuneCross[["performanceType"]])

              # Selects features and fits a model on all of the training samples, as runTests does for its final model.
              fitModel <- function(measurementsUse, assayIndex, multiView, selectionUse, classifierUse)
              {
                modellingParams <- generateModellingParams(assayIDs = assayIndex, measurements = measurementsUse,
                                                           nFeatures = nFeatures[assayIndex], selectionMethod = selectionUse,
                                                           classifier = classifierUse, multiViewMethod = multiView, extraParams = extraParams)
                if(is.null(modellingParams@predictParams))
                  stop("Classifier ", paste(unlist(classifierUse), collapse = ", "), " trains and predicts in one step, so it can't be trained on its own.")
                trained <- runTest(measurementsUse, outcomeTrain, measurementsUse, outcomeTrain, crossValParams = crossValParams,
                                   modellingParams = modellingParams, verbose = verbose, .iteration = 1)
                if(is.character(trained)) stop(trained)
                model <- trained[["models"]]
                attr(model, "predictFunction") <- modellingParams@predictParams@predictor
                model
              }

              if(multiViewMethod == "none")
              {
                models <- unlist(lapply(assayIDs, function(assayIndex)
                {
                  measurementsUse <- measurements
                  if(assayIndex != 1) measurementsUse <- measurements[, S4Vectors::mcols(measurements)[, "assay"] == assayIndex, drop = FALSE]
                  unlist(lapply(classifier[[assayIndex]], function(classifierForAssay)
                    lapply(selectionMethod[[assayIndex]], function(selectionForAssay)
                      fitModel(measurementsUse, assayIndex, "none", selectionForAssay, classifierForAssay))), recursive = FALSE)
                }), recursive = FALSE)
              } else { # Merge, prevalidation or PCA combine all of the chosen assays into one model.
                measurementsUse <- measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assayIDs, drop = FALSE]
                models <- list(fitModel(measurementsUse, assayIDs, multiViewMethod, selectionMethod[assayIDs], classifier[assayIDs]))
              }

              if(length(models) == 1)
              {
                model <- models[[1]]
                class(model) <- c("trainedByClassifyR", class(model))
                return(model)
              }
              class(models) <- c("listOfModels", "trainedByClassifyR", class(models))
              models
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
              
              # The samples of every table in the order of the first table's.
              if(!is.null(rownames(x[[1]])))
                x <- lapply(x, function(measurements) measurements[rownames(x[[1]]), , drop = FALSE])
              df_list <- lapply(x, S4Vectors::DataFrame, check.names = FALSE)
              
              # Features are named assay_feature, as crossValidate and predict name them.
              df_list <- mapply(function(meas, nam){
                  S4Vectors::mcols(meas)$assay <- nam
                  S4Vectors::mcols(meas)$feature <- colnames(meas)
                  colnames(meas) <- paste(nam, colnames(meas), sep = '_')
                  meas
              }, df_list, names(df_list))
              
              combined_df <- do.call(cbind, unname(df_list))
              
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
  # Name the features of the new data as train() names them.
  if(is(newData, "MultiAssayExperiment"))
  {
    newData <- prepareData(newData, outcome)[["measurements"]]
  } else if(is(newData, "tabular")) {
    newData <- S4Vectors::DataFrame(newData, check.names = FALSE)
  } else if(is.list(newData)) { # Features of several assays are named assay_feature.
    if(!is.null(rownames(newData[[1]]))) # The samples of every table in the order of the first table's.
      newData <- lapply(newData, function(measurements) measurements[rownames(newData[[1]]), , drop = FALSE])
    newData <- do.call(cbind, mapply(function(measurementsOne, assayID)
    {
      measurementsOne <- S4Vectors::DataFrame(measurementsOne, check.names = FALSE)
      colnames(measurementsOne) <- paste(assayID, colnames(measurementsOne), sep = '_')
      measurementsOne
    }, newData, names(newData), SIMPLIFY = FALSE) |> unname())
  }
  colnames(newData) <- make.names(colnames(newData), unique = TRUE)

  # Each model predicts from the features it was trained with.
  predictOne <- function(model)
  {
    predictFunction <- attr(model, "predictFunction")
    features <- attr(model, "featuresForTrain")
    class(model) <- setdiff(class(model), "trainedByClassifyR")
    if(!is.null(features)) newData <- newData[, features, drop = FALSE]
    extras <- list(...)
    if(!"verbose" %in% names(extras)) extras[["verbose"]] <- 0
    do.call(predictFunction, c(list(model, newData), extras))
  }
  if(is(object, "listOfModels")) lapply(unclass(object), predictOne) else predictOne(object)
}
