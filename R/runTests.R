#' Reproducibly Run Various Kinds of Cross-Validation
#' 
#' Enables doing classification schemes such as ordinary 10-fold, 100
#' permutations 5-fold, and leave one out cross-validation. Processing in
#' parallel is possible by leveraging the package \code{\link{BiocParallel}}.
#' 
#' 
#' @aliases runTests runTests,matrix-method runTests,DataFrame-method
#' runTests,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing all of the data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns
#' are features.
#' @param outcome Either a factor vector of classes, a \code{\link{Surv}} object, or
#' a character string, or vector of such strings, containing column name(s) of column(s)
#' containing either classes or time and event information about survival. If
#' \code{measurements} is a \code{MultiAssayExperiment}, the names of the column (class) or
#' columns (survival) in the table extracted by \code{colData(data)} that contain(s) the samples'
#' outcome to use for prediction. If column names of survival information, time must be in first
#' column and event status in the second.
#' @param crossValParams An object of class \code{\link{CrossValParams}},
#' specifying the kind of cross-validation to be done.
#' @param modellingParams An object of class \code{\link{ModellingParams}},
#' specifying the class rebalancing, transformation (if any), feature selection
#' (if any), training and prediction to be done on the data set.
#' @param characteristics A \code{\link{DataFrame}} describing the
#' characteristics of the classification used. First column must be named
#' \code{"charateristic"} and second column must be named \code{"value"}.
#' Useful for automated plot annotation by plotting functions within this
#' package.  Transformation, selection and prediction functions provided by
#' this package will cause the characteristics to be automatically determined
#' and this can be left blank.
#' @param ... Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
#' are passed into and used by the \code{DataFrame} method or passed onwards to \code{\link{prepareData}}.
#' @param verbose Default: 1. A number between 0 and 3 for the amount of
#' progress messages to give.  A higher number will produce more messages as
#' more lower-level functions print messages.
#' @return An object of class \code{\link{ClassifyResult}}.
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     data(asthma)
#'     
#'     CVparams <- CrossValParams(permutations = 5, tuneMode = "Resubstitution")
#'     tuneList <- list(nFeatures = seq(5, 25, 5))
#'     attr(tuneList, "performanceType") <- "Balanced Error"
#'     selectParams <- SelectParams("t-test", tuneParams = tuneList)
#'     modellingParams <- ModellingParams(selectParams = selectParams)
#'     runTests(measurements, classes, CVparams, modellingParams,
#'              DataFrame(characteristic = c("Assay Name", "Classifier Name"),
#'                        value = c("Asthma", "Different Means"))
#'              )
#'   #}
#'
#' @export
#' @usage NULL
setGeneric("runTests", function(measurements, ...) standardGeneric("runTests"))

#' @rdname runTests
#' @export
setMethod("runTests", c("matrix"), function(measurements, outcome, ...) # Matrix of numeric measurements.
{
  if(is.null(rownames(measurements)))
    stop("'measurements' matrix must have sample identifiers as its row names.")
  runTests(S4Vectors::DataFrame(measurements, check.names = FALSE), outcome, ...)
})

#' @rdname runTests
#' @import BiocParallel
#' @export
setMethod("runTests", "DataFrame", function(measurements, outcome, crossValParams = CrossValParams(), modellingParams = ModellingParams(),
           characteristics = S4Vectors::DataFrame(), ..., verbose = 1)
{
  crossValidation <- .prepareTests(measurements, outcome, crossValParams, modellingParams, characteristics, verbose, ...)
  results <- .runTestsSplits(list(crossValidation), crossValParams@parallelParams)[[1]]
  .assembleTests(crossValidation, results)
})

# A cross-validation is run in three steps, so that the splits of several cross-validations can share one pool of
# parallel workers (see crossValidate):
# 1. .prepareTests checks the data, makes the training and test splits and fits the final model to all samples.
# 2. .runTestsSplits runs the training and testing of every split of one or more cross-validations.
# 3. .assembleTests collects one cross-validation's results into a ClassifyResult.
# Each split uses the random number stream that bpmapply with the cross-validation's RNGseed would give it, so the
# results don't depend on how the splits are distributed among workers.

.prepareTests <- function(measurements, outcome, crossValParams, modellingParams, characteristics, verbose, ..., splits = NULL, deferFinal = FALSE)
{
  if(is.null(rownames(measurements)))
  {
    warning("'measurements' DataFrame must have sample identifiers as its row names. Generating generic ones.")
    rownames(measurements) <- paste("Sample", seq_len(nrow(measurements)))
  }
  
  if(any(vapply(as.list(measurements), anyNA, logical(1))))
    stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")            

  originalFeatures <- colnames(measurements)
  if("assay" %in% colnames(S4Vectors::mcols(measurements)))
      originalFeatures <- S4Vectors::mcols(measurements)[, c("assay", "feature")]                 
  splitDataset <- prepareData(measurements, outcome, ...)
  measurements <- splitDataset[["measurements"]]
  outcome <- splitDataset[["outcome"]]
  
  if(crossValParams@performanceType == "auto")
  {
    if(is.factor(outcome)) crossValParams@performanceType <- "Balanced Accuracy" else
      crossValParams@performanceType <- "C-index"    
  }
  if(!is.null(modellingParams@selectParams))
  {
    nFeatures <- modellingParams@selectParams@tuneParams[["nFeatures"]]
    if(is.null(nFeatures)) nFeatures <- modellingParams@selectParams@nFeatures
  }
  if(!is.null(modellingParams@selectParams) && max(nFeatures) > ncol(measurements))
  {
    warning("Attempting to evaluate more features for feature selection than in
input data. Autmomatically reducing to smaller number.")
    if(is.null(modellingParams@selectParams@nFeatures))
      modellingParams@selectParams@tuneParams[["nFeatures"]][modellingParams@selectParams@tuneParams[["nFeatures"]] > max(nFeatures)] <- max(nFeatures)
    else
      modellingParams@selectParams@nFeatures <- max(nFeatures)
  }

  # Create all partitions of training and testing sets, unless given ones shared with other cross-validations.
  if(!is.null(splits)) samplesSplitsList <- splits else
  samplesSplitsList <- samplesSplits(crossValParams@samplesSplits, crossValParams@permutations, crossValParams@folds, crossValParams@percentTest, crossValParams@leave, outcome)
  splitsTestInfoTable <- splitsTestInfo(crossValParams@samplesSplits, crossValParams@permutations, crossValParams@folds, crossValParams@percentTest, crossValParams@leave, samplesSplitsList)

  # The final model is fitted with the random number state that follows making the splits. If deferred, it is
  # fitted alongside the splits (by .runTestsSplits) from that state.
  finalState <- if(exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  if(deferFinal) fullResult <- NULL else
  fullResult <- runTest(measurements, outcome, measurements, outcome, crossValParams = crossValParams, modellingParams = modellingParams, characteristics = characteristics, .iteration = 1)

  list(measurements = measurements, outcome = outcome, originalFeatures = originalFeatures, crossValParams = crossValParams,
       modellingParams = modellingParams, characteristics = characteristics, verbose = verbose,
       splits = samplesSplitsList, splitsInfo = splitsTestInfoTable, fullResult = fullResult, finalState = finalState)
}

# Random number streams of the elements of a bplapply or bpmapply call with RNGseed equal to seed.
.splitStreams <- function(seed, nSplits)
{
  if(is.null(seed)) return(vector("list", nSplits))
  streams <- vector("list", nSplits)
  streams[[1]] <- BiocParallel:::.rng_init_stream(seed)
  for(splitIndex in seq_len(nSplits - 1))
    streams[[splitIndex + 1]] <- parallel::nextRNGSubStream(streams[[splitIndex]])
  streams
}

# Runs every split of the cross-validations, and their deferred final models (split 0), in one call of bplapply.
# Returns, for each cross-validation, its list of split results with the final model's result as attribute
# "fullResult" when it was deferred.
.runTestsSplits <- function(crossValidations, parallelParams)
{
  tasks <- do.call(rbind, lapply(seq_along(crossValidations), function(index)
  {
    crossValidation <- crossValidations[[index]]
    splitNumbers <- seq_along(crossValidation[["splits"]][["train"]])
    if(is.null(crossValidation[["fullResult"]])) splitNumbers <- c(0L, splitNumbers)
    data.frame(crossValidation = index, split = splitNumbers)
  }))
  streams <- unlist(lapply(seq_along(crossValidations), function(index)
  {
    crossValidation <- crossValidations[[index]]
    splitStreams <- .splitStreams(BiocParallel::bpRNGseed(crossValidation[["crossValParams"]]@parallelParams),
                                  length(crossValidation[["splits"]][["train"]]))
    if(is.null(crossValidation[["fullResult"]])) splitStreams <- c(list(crossValidation[["finalState"]]), splitStreams)
    splitStreams
  }), recursive = FALSE)

  results <- bplapply(seq_len(nrow(tasks)), function(taskIndex)
  {
    if(!is.null(streams[[taskIndex]])) assign(".Random.seed", streams[[taskIndex]], envir = globalenv())
    crossValidation <- crossValidations[[tasks[taskIndex, "crossValidation"]]]
    setNumber <- tasks[taskIndex, "split"]
    if(setNumber == 0) # The final model, fitted to all samples.
      return(runTest(crossValidation[["measurements"]], crossValidation[["outcome"]], crossValidation[["measurements"]], crossValidation[["outcome"]],
                     crossValParams = crossValidation[["crossValParams"]], modellingParams = crossValidation[["modellingParams"]],
                     characteristics = crossValidation[["characteristics"]], .iteration = 1))
    if(crossValidation[["verbose"]] >= 1 && setNumber %% 10 == 0)
      message(Sys.time(), ": Processing sample set ", setNumber, '.')
    
    trainingSamples <- crossValidation[["splits"]][["train"]][[setNumber]]
    testSamples <- crossValidation[["splits"]][["test"]][[setNumber]]
    # crossValParams is needed at least for nested feature tuning.
    result <- runTest(crossValidation[["measurements"]][trainingSamples, , drop = FALSE], crossValidation[["outcome"]][trainingSamples],
                      crossValidation[["measurements"]][testSamples, , drop = FALSE], crossValidation[["outcome"]][testSamples],
                      crossValidation[["crossValParams"]], crossValidation[["modellingParams"]], crossValidation[["characteristics"]],
                      crossValidation[["verbose"]], .iteration = setNumber)
    # A random forest grown only to rank features has been used by now; fold models don't keep it.
    if(is.list(result) && !is.null(attr(result[["models"]], "forImportance")))
      attr(result[["models"]], "forImportance") <- NULL
    result
  }, BPPARAM = parallelParams)
  lapply(unname(split(seq_len(nrow(tasks)), tasks[, "crossValidation"])), function(taskIndices)
  {
    isFinal <- tasks[taskIndices, "split"] == 0
    splitResults <- results[taskIndices[!isFinal]]
    if(any(isFinal)) attr(splitResults, "fullResult") <- results[[taskIndices[isFinal]]]
    splitResults
  })
}

.assembleTests <- function(crossValidation, results)
{
  measurements <- crossValidation[["measurements"]]
  outcome <- crossValidation[["outcome"]]
  crossValParams <- crossValidation[["crossValParams"]]
  modellingParams <- crossValidation[["modellingParams"]]
  characteristics <- crossValidation[["characteristics"]]
  splitsTestInfoTable <- crossValidation[["splitsInfo"]]
  fullResult <- crossValidation[["fullResult"]]
  if(is.null(fullResult)) fullResult <- attr(results, "fullResult")
  attr(results, "fullResult") <- NULL

  # Error checking and reporting.
  resultErrors <- sapply(results, function(result) is.character(result))
  if(sum(resultErrors) == length(results))
  {
      message(Sys.time(), " - Error: All cross-validations had an error.")
      if(length(unique(unlist(results))) == 1)
        stop("The common problem is: ", unlist(results)[[1]])
      return(results)
  } else if(sum(resultErrors) != 0) # Filter out cross-validations resulting in error.
  {
    warning(paste(sum(resultErrors),  "cross-validations, but not all, had an error and have been removed from the results."))
    results <- results[!resultErrors]
    iterationID <- do.call(paste, as.data.frame(splitsTestInfoTable))
    iterationIDlevels <- unique(iterationID)
    errorRows <- iterationID %in% iterationIDlevels[which(resultErrors)]
    splitsTestInfoTable <- splitsTestInfoTable[!errorRows, ]
  }
  
  validationText <- .validationText(crossValParams)
  
  modParamsList <- list(modellingParams@transformParams, modellingParams@selectParams, modellingParams@trainParams, modellingParams@predictParams)
  autoCharacteristics <- lapply(modParamsList, function(stageParams) if(!is.null(stageParams) && !is(stageParams, "PredictParams")) stageParams@characteristics)
  autoCharacteristics <- do.call(rbind, autoCharacteristics)

  # Add extra settings.
  extras <- unlist(lapply(modParamsList, function(stageParams) if(!is.null(stageParams)) stageParams@otherParams), recursive = FALSE)
  if(length(extras) > 0)
    extras <- extras[sapply(extras, is.atomic)] # Store basic variables, not complex ones.
  extrasDF <- S4Vectors::DataFrame(characteristic = names(extras), value = unname(unlist(extras)))
  characteristics <- rbind(characteristics, extrasDF)
  characteristics <- .filterCharacteristics(characteristics, autoCharacteristics)
  characteristics <- rbind(characteristics,
                             S4Vectors::DataFrame(characteristic = "Cross-validation", value = validationText))

  if(!is.data.frame(results[[1]][["predictions"]]))
  {
    if(is.numeric(results[[1]][["predictions"]])) # Survival task.
        predictsColumnName <- "risk"
    else # Classification task. A factor.
        predictsColumnName <- "class"
    predictionsTable <- S4Vectors::DataFrame(sample = unlist(lapply(results, "[[", "testSet")), splitsTestInfoTable, unlist(lapply(results, "[[", "predictions")), check.names = FALSE)
    colnames(predictionsTable)[ncol(predictionsTable)] <- predictsColumnName
  } else { # data frame
    predictionsTable <- S4Vectors::DataFrame(sample = unlist(lapply(results, "[[", "testSet")), splitsTestInfoTable, do.call(rbind, lapply(results, "[[", "predictions")), check.names = FALSE)
  }
  rownames(predictionsTable) <- NULL
  tuneList <- lapply(results, "[[", "tune")
  if(length(unlist(tuneList)) == 0)
    tuneList <- NULL
  importance <- NULL
  if(!is.null(results[[1]][["importance"]]))
    importance <- do.call(rbind, lapply(results, "[[", "importance"))
  
  if(is.character(fullResult))
  {
    warning("Unable to fit a full model: ", fullResult)
    fullResult <- list(models = NULL)
  }
  
  ClassifyResult(characteristics, rownames(measurements), crossValidation[["originalFeatures"]],
                 lapply(results, "[[", "ranked"), lapply(results, "[[", "selected"),
                 lapply(results, "[[", "models"), tuneList, predictionsTable, outcome, importance, modellingParams, fullResult$models)
}

#' @rdname runTests
#' @import MultiAssayExperiment methods
#' @export
setMethod("runTests", c("MultiAssayExperiment"),
          function(measurements, outcome, ...)
{
  prepArgs <- list(measurements, outcome)              
  extraInputs <- list(...)
  prepExtras <- numeric()
  if(length(extraInputs) > 0)
    prepExtras <- which(names(extraInputs) %in% .ClassifyRenvir[["prepareDataFormals"]])
  if(length(prepExtras) > 0)
    prepArgs <- append(prepArgs, extraInputs[prepExtras])
  measurementsAndOutcome <- do.call(prepareData, prepArgs)
  
  runTestsArgs <- list(measurementsAndOutcome[["measurements"]], measurementsAndOutcome[["outcome"]])
  if(length(extraInputs) > 0 && (length(prepExtras) == 0 || length(extraInputs[-prepExtras]) > 0))
  {
    if(length(prepExtras) == 0) runTestsArgs <- append(runTestsArgs, extraInputs) else
    runTestsArgs <- append(runTestsArgs, extraInputs[-prepExtras])
  }
  do.call(runTests, runTestsArgs)
})
