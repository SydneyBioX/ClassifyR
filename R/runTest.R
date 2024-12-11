#' Perform a Single Classification
#' 
#' For a data set of features and samples, the classification process is run.
#' It consists of data transformation, feature selection, classifier training
#' and testing.
#' 
#' This function only performs one classification and prediction. See
#' \code{\link{runTests}} for a driver function that enables a number of
#' different cross-validation schemes to be applied and uses this function to
#' perform each iteration.
#' 
#' @aliases runTest runTest,matrix-method runTest,DataFrame-method
#' runTest,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' @param outcomeTrain Either a factor vector of classes, a \code{\link{Surv}} object, or
#' a character string, or vector of such strings, containing column name(s) of column(s)
#' containing either classes or time and event information about survival. If column names
#' of survival information, time must be in first column and event status in the second.
#' @param measurementsTest Same data type as \code{measurementsTrain}, but only the test
#' samples.
#' @param outcomeTest Same data type as \code{outcomeTrain}, but for only the test
#' samples.
#' @param crossValParams An object of class \code{\link{CrossValParams}},
#' specifying the kind of cross-validation to be done, if nested
#' cross-validation is used to tune any parameters.
#' @param modellingParams An object of class \code{\link{ModellingParams}},
#' specifying the class rebalancing, transformation (if any), feature selection
#' (if any), training and prediction to be done on the data set.
#' @param outcomeColumns If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' names of the column (class) or columns (survival) in the table extracted by \code{colData(data)}
#' that contain(s) the samples' outcome to use for prediction.
#' @param ... Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
#' are passed into and used by the \code{DataFrame} method or passed onwards to \code{\link{prepareData}}.
#' @param characteristics A \code{\link{DataFrame}} describing the
#' characteristics of the classification used. First column must be named
#' \code{"charateristic"} and second column must be named \code{"value"}.
#' Useful for automated plot annotation by plotting functions within this
#' package. Transformation, selection and prediction functions provided by
#' this package will cause the characteristics to be automatically determined
#' and this can be left blank.
#' @param verbose Default: 1. A number between 0 and 3 for the amount of
#' progress messages to give.  A higher number will produce more messages as
#' more lower-level functions print messages.
#' @param .iteration Not to be set by a user. This value is used to keep track
#' of the cross-validation iteration, if called by \code{\link{runTests}}.
#' @return If called directly by the user rather than being used internally by
#' \code{\link{runTests}}, a \code{\link{ClassifyResult}} object. Otherwise a
#' list of different aspects of the result which is passed back to \code{\link{runTests}}.
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     data(asthma)
#'     CVparams <- CrossValParams(tuneMode = "Resubstitution")
#'     tuneList <- list(nFeatures = seq(5, 25, 5))
#'     attr(tuneList, "performanceType") <- "Balanced Error"
#'     selectParams <- SelectParams("limma", tuneParams = tuneList)
#'     modellingParams <- ModellingParams(selectParams = selectParams)
#'     trainIndices <- seq(1, nrow(measurements), 2)
#'     testIndices <- seq(2, nrow(measurements), 2)
#'     
#'     runTest(measurements[trainIndices, ], classes[trainIndices],
#'             measurements[testIndices, ], classes[testIndices],
#'             crossValParams = CVparams, modellingParams = modellingParams)
#'   #}
#' 
#' @importFrom S4Vectors do.call mcols
#' @usage NULL
#' @export
setGeneric("runTest", function(measurementsTrain, ...) standardGeneric("runTest"))

#' @rdname runTest
#' @export
setMethod("runTest", "matrix", # Matrix of numeric measurements.
  function(measurementsTrain, outcomeTrain, measurementsTest, outcomeTest, ...)
{
  runTest(measurementsTrain = S4Vectors::DataFrame(measurementsTrain, check.names = FALSE),
          outcomeTrain = outcomeTrain,
          measurementsTest = S4Vectors::DataFrame(measurementsTest, check.names = FALSE),
          outcomeTest = outcomeTest,
          ...)
})

#' @rdname runTest
#' @export
setMethod("runTest", "DataFrame", # Sample information data or one of the other inputs, transformed.
function(measurementsTrain, outcomeTrain, measurementsTest, outcomeTest,
         crossValParams = CrossValParams(), # crossValParams might be used for tuning optimisation.
         modellingParams = ModellingParams(), characteristics = S4Vectors::DataFrame(), ..., verbose = 1, .iteration = NULL)
{
  if(is.null(.iteration)) # Not being called by runTests but by user. So, check the user input.
  {
    if(is.null(rownames(measurementsTrain)))
      stop("'measurementsTrain' DataFrame must have sample identifiers as its row names.")
    if(any(is.na(measurementsTrain)))
      stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")                
    
    splitDatasetTrain <- prepareData(measurementsTrain, outcomeTrain, ...)
    if(crossValParams@performanceType == "auto")
    {
      if(is.factor(splitDatasetTrain[["outcome"]])) crossValParams@performanceType <- "Balanced Accuracy" else
        crossValParams@performanceType <- "C-index"    
    }
      
    # Rebalance the class sizes of the training samples by either downsampling or upsampling
    # or leave untouched if balancing is none.
    if(!is(outcomeTrain, "Surv"))
    {
      rebalancedTrain <- .rebalanceTrainingClasses(splitDatasetTrain[["measurements"]], splitDatasetTrain[["outcome"]], modellingParams@balancing)
      measurementsTrain <- rebalancedTrain[["measurementsTrain"]]
      outcomeTrain <- rebalancedTrain[["classesTrain"]]
    }
  }
  
  if("feature" %in% colnames(S4Vectors::mcols(measurementsTrain))) originalFeatures <- S4Vectors::mcols(measurementsTrain)[, na.omit(match(c("assay", "feature"), colnames(S4Vectors::mcols(measurementsTrain))))]    
  else originalFeatures <- colnames(measurementsTrain)

  if(!is.null(modellingParams@selectParams))
  {
    nFeatures <- modellingParams@selectParams@tuneParams[["nFeatures"]]
    if(is.null(nFeatures)) nFeatures <- modellingParams@selectParams@nFeatures
  }
  if(!is.null(modellingParams@selectParams) && max(nFeatures) > ncol(measurementsTrain))
  {
    warning("Attempting to evaluate more features for feature selection than in
input data. Autmomatically reducing to smaller number.")
    if(is.null(modellingParams@selectParams@nFeatures))
      modellingParams@selectParams@tuneParams[["nFeatures"]][modellingParams@selectParams@tuneParams[["nFeatures"]] > max(nFeatures)] <- max(nFeatures)
    else
      modellingParams@selectParams@nFeatures <- max(nFeatures)
  }
  
  if(!is.null(crossValParams) && !is.null(crossValParams@adaptiveResamplingDelta))
  { # Iteratively resample training samples until their class probability or risk changes, on average, less than delta.
    delta <- crossValParams@adaptiveResamplingDelta
    crossValParams@adaptiveResamplingDelta <- NULL
    scoresPrevious <- rep(1, nrow(measurementsTrain))
    repeat{
      newSamples <- sample(nrow(measurementsTrain), replace = TRUE, prob = scoresPrevious)
      measurementsTrainResampled <- measurementsTrain[newSamples, ]
      outcomeResampled <- outcomeTrain[newSamples]
      ASpredictions <- runTest(measurementsTrainResampled, outcomeResampled,
              measurementsTrain, outcomeTrain, crossValParams, modellingParams,
              .iteration = "internal")[["predictions"]]
      if(is.factor(outcomeResampled))
          scoresNew <- mapply(function(rowIndex, class) ASpredictions[rowIndex, class], 1:nrow(ASpredictions), as.character(outcomeTrain))
      else
          scoresNew <- ASpredictions[, "risk"]

      if(mean(abs(scoresNew - scoresPrevious)) < delta)
      {
        measurementsTrain <- measurementsTrainResampled
        break;
      }
      scoresPrevious <- scoresNew
    }
  }
  
  if(!is.null(modellingParams@transformParams))
  {
    if(length(modellingParams@transformParams@intermediate) != 0)
      modellingParams@transformParams <- .addIntermediates(modellingParams@transformParams)
    
    measurementsTransformedList <- tryCatch(.doTransform(measurementsTrain, measurementsTest, modellingParams@transformParams, verbose), error = function(error) error[["message"]])
    if(is.character(measurementsTransformedList[[1]])) return(measurementsTransformedList[[1]]) # An error occurred.
  }

  rankedFeaturesIndices <- NULL
  selectedFeaturesIndices <- NULL
  tuneDetailsSelect <- NULL

  if(!is.null(modellingParams@selectParams))
  {
    if(length(modellingParams@selectParams@intermediate) != 0)
      modellingParams@selectParams <- .addIntermediates(modellingParams@selectParams)
 
    #topFeatures <- .doSelection(measurementsTrain, outcomeTrain, crossValParams, modellingParams, verbose)
    topFeatures <- tryCatch(.doSelection(measurementsTrain, outcomeTrain, crossValParams, modellingParams, verbose),
                            error = function(error) error[["message"]])
    if(is.character(topFeatures)) return(topFeatures) # An error occurred.
    
    rankedFeaturesIndices <- topFeatures[[1]] # Extract for result object.
    selectedFeaturesIndices <- topFeatures[[2]] # Extract for subsetting.
    tuneDetailsSelect <- topFeatures[[3]]

    if(modellingParams@selectParams@subsetToSelections == TRUE)
    {
      measurementsTrain <- measurementsTrain[, selectedFeaturesIndices, drop = FALSE]
      measurementsTest <- measurementsTest[, selectedFeaturesIndices, drop = FALSE]
    }
  }

  # Training stage.
  if(length(modellingParams@trainParams@intermediate) > 0)
    modellingParams@trainParams <- .addIntermediates(modellingParams@trainParams)
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
  
  # Some classifiers have one function for training and testing, so that's why test data is also passed in.
  # trained <- .doTrain(measurementsTrain, outcomeTrain, measurementsTest, outcomeTest, crossValParams, modellingParams, verbose)
  trained <- tryCatch(.doTrain(measurementsTrain, outcomeTrain, measurementsTest, outcomeTest, crossValParams, modellingParams, verbose),
                      error = function(error) error[["message"]])
  if(is.character(trained)) return(trained) # An error occurred.
  
  tuneDetailsTrain <- trained[[2]] # Second element is tuning results.
  
  if(!is.null(modellingParams@trainParams@getFeatures)) # Features chosen inside classifier.
  {
    extrasList <- list()
    extras <- .methodFormals(modellingParams@trainParams@getFeatures, class(trained[[1]]))[-1]
    if(length(extras) > 0)
    {
      if("classesTrain" %in% names(extras))
      {
         whichRename <- match("classesTrain", names(extras))
         names(extras)[whichRename] <- "outcomeTrain"
         extrasList <- mget(setdiff(names(extras), "..."))
         whichRename <- match("outcomeTrain", names(extrasList))
         names(extrasList)[whichRename] <- "classesTrain"
      } else {extrasList <- mget(setdiff(names(extras), "..."))}
    }
    rankedChosenList <- do.call(modellingParams@trainParams@getFeatures, c(unname(trained[1]), extrasList))
    rankedFeaturesIndices <- rankedChosenList[[1]]
    selectedFeaturesIndices <- rankedChosenList[[2]]
  }
  
  if(!is.null(modellingParams@predictParams))
  {
    if(length(modellingParams@predictParams@intermediate) != 0)
      modellingParams@predictParams <- .addIntermediates(modellingParams@predictParams)
    
    #predictedOutcome <- .doTest(trained[["model"]], measurementsTest, modellingParams@predictParams, verbose)
    predictedOutcome <- tryCatch(.doTest(trained[["model"]], measurementsTest, modellingParams@predictParams, verbose),
                                 error = function(error) error[["message"]])
    
    if(is.character(predictedOutcome)) # An error occurred.
      return(predictedOutcome) # Return early.
    
  } else { # One function that does training and testing, so predictions were made earlier
           # by .doTrain, rather than this .doTest stage.
    predictedOutcome <- trained[[1]]
  }
  
  # Exclude one feature at a time, build model, predict test samples.
  importanceTable <- NULL
  if(is.numeric(.iteration) && modellingParams@doImportance == TRUE)
  {
    performanceMP <- modellingParams@selectParams@tuneParams[["performanceType"]]
    performanceType <- ifelse(!is.null(performanceMP), performanceMP, "Balanced Error")
    performancesWithoutEach <- sapply(selectedFeaturesIndices, function(selectedIndex)
    {
      measurementsTrainLess1 <- measurementsTrain[, -selectedIndex, drop = FALSE]
      measurementsTestLess1 <- measurementsTest[, -selectedIndex, drop = FALSE]
      modelWithoutOne <- tryCatch(.doTrain(measurementsTrainLess1, outcomeTrain, measurementsTestLess1, outcomeTest, crossValParams, modellingParams, verbose),
                                  error = function(error) error[["message"]])
      if(!is.null(modellingParams@predictParams))
      predictedOutcomeWithoutOne <- tryCatch(.doTest(modelWithoutOne[["model"]], measurementsTestLess1, modellingParams@predictParams, verbose),
                                              error = function(error) error[["message"]])
      else predictedOutcomeWithoutOne <- modelWithoutOne[["model"]]

      if(!is.null(ncol(predictedOutcomeWithoutOne)))
        predictedOutcomeWithoutOne <- predictedOutcomeWithoutOne[, na.omit(match(c("class", "risk"), colnames(predictedOutcomeWithoutOne)))]
      calcExternalPerformance(outcomeTest, predictedOutcomeWithoutOne, performanceType)
    })
    
    if(!is.null(ncol(predictedOutcome)))
        predictedOutcome <- predictedOutcome[, na.omit(match(c("class", "risk"), colnames(predictedOutcome)))]
    performanceChanges <- round(performancesWithoutEach - calcExternalPerformance(outcomeTest, predictedOutcome, performanceType), 2)

    if(is.null(S4Vectors::mcols(measurementsTrain)) || !any(c("assay", "feature") %in% colnames(S4Vectors::mcols(measurementsTrain))))
    {
      if(is.null(modellingParams@trainParams@getFeatures))
      selectedFeatures <- originalFeatures[selectedFeaturesIndices]
      else selectedFeatures <- colnames(measurementsTrain)[selectedFeaturesIndices] 
    } else {
      featureColumns <- na.omit(match(c("assay", "feature"), colnames(S4Vectors::mcols(measurementsTrain))))
      if(length(featureColumns) == 1)
      {
         if(is.null(modellingParams@trainParams@getFeatures))
            selectedFeatures <- originalFeatures[selectedFeaturesIndices]
         else selectedFeatures <- S4Vectors::mcols(measurementsTrain)[selectedFeaturesIndices, featureColumns]            
      } else {
                if(is.null(modellingParams@trainParams@getFeatures))
                  selectedFeatures <- originalFeatures[selectedFeaturesIndices, ]
                else selectedFeatures <- S4Vectors::mcols(measurementsTrain)[selectedFeaturesIndices, featureColumns]            
      }
    }
    importanceTable <- S4Vectors::DataFrame(selectedFeatures, performanceChanges)
    if(ncol(importanceTable) == 2) colnames(importanceTable)[1] <- "feature"
    colnames(importanceTable)[ncol(importanceTable)] <- paste("Change in", performanceType)
  }
  
  if(is.null(modellingParams@predictParams)) models <- NULL else models <- trained[[1]] # One function for training and testing. Typically, the models aren't returned to the user, such as Poisson LDA implemented by PoiClaClu.
  if(!is.null(tuneDetailsSelect)) tuneDetails <- tuneDetailsSelect else tuneDetails <- tuneDetailsTrain

  # Convert back into original, potentially unsafe feature identifiers unless it is a nested cross-validation.
  if(is.null(.iteration) || .iteration != "internal")
  {
    if(!is.null(rankedFeaturesIndices))
    {
      if(!is.numeric(rankedFeaturesIndices)) # Metafeatures created by prevalidation or PCA.
      {
        rankedFeatures <- rankedFeaturesIndices
      } else if(is.null(S4Vectors::mcols(measurementsTrain)) || !any(c("assay", "feature") %in% colnames(S4Vectors::mcols(measurementsTrain))))
      {
        if(is.null(modellingParams@trainParams@getFeatures))          
          rankedFeatures <- originalFeatures[rankedFeaturesIndices]
        else rankedFeatures <- colnames(measurementsTrain)[rankedFeaturesIndices]            
      } else {
        featureColumns <- na.omit(match(c("assay", "feature"), colnames(S4Vectors::mcols(measurementsTrain))))          
        if(length(featureColumns) == 1)
        {
          if(is.null(modellingParams@trainParams@getFeatures))
            rankedFeatures <- originalFeatures[rankedFeaturesIndices]
          else rankedFeatures <- S4Vectors::mcols(measurementsTrain)[rankedFeaturesIndices, featureColumns]
        } else {
          if(is.null(modellingParams@trainParams@getFeatures))
            rankedFeatures <- originalFeatures[rankedFeaturesIndices, ]
          else rankedFeatures <- S4Vectors::mcols(measurementsTrain)[rankedFeaturesIndices, featureColumns] 
        }
      }
    } else { rankedFeatures <- NULL}
    if(!is.null(selectedFeaturesIndices))
    {
      if(!is.numeric(selectedFeaturesIndices)) # Metafeatures created by prevalidation or PCA.
      {
          selectedFeatures <- selectedFeaturesIndices
      } else if(is.null(S4Vectors::mcols(measurementsTrain)) || !any(c("assay", "feature") %in% colnames(S4Vectors::mcols(measurementsTrain))))
      {
        if(is.null(modellingParams@trainParams@getFeatures))
          selectedFeatures <- originalFeatures[selectedFeaturesIndices]
        else selectedFeatures <- colnames(measurementsTrain)[selectedFeaturesIndices]
      } else {
        featureColumns <- na.omit(match(c("assay", "feature"), colnames(S4Vectors::mcols(measurementsTrain))))  
        if(length(featureColumns) == 1)
        {
          if(is.null(modellingParams@trainParams@getFeatures))
            selectedFeatures <- originalFeatures[selectedFeaturesIndices]
          else selectedFeatures <- S4Vectors::mcols(measurementsTrain)[selectedFeaturesIndices, featureColumns]            
        } else {
                if(is.null(modellingParams@trainParams@getFeatures))
                  selectedFeatures <- originalFeatures[selectedFeaturesIndices, ]
                else selectedFeatures <- S4Vectors::mcols(measurementsTrain)[selectedFeaturesIndices, featureColumns]            
        }
      }
    } else { selectedFeatures <- NULL}
  } else { # Nested use in feature selection. No feature selection in inner execution, so ignore features. 
      rankedFeatures <- selectedFeatures <- NULL
  }

  if(!is.null(.iteration)) # This function was not called by the end user.
  {
    list(ranked = rankedFeatures, selected = selectedFeatures, models = models, testSet = rownames(measurementsTest), predictions = predictedOutcome, tune = tuneDetails, importance = importanceTable)
  } else { # runTest executed by the end user. Create a ClassifyResult object.
    # Only one training, so only one tuning choice, which can be summarised in characteristics.
    modParamsList <- list(modellingParams@transformParams, modellingParams@selectParams, modellingParams@trainParams, modellingParams@predictParams)
    if(!is.null(tuneDetails)) characteristics <- rbind(characteristics, data.frame(characteristic = colnames(tuneDetails[["tuneCombinations"]]),
                                                                                   value = unlist(tuneDetails[["tuneCombinations"]][tuneDetails[["bestIndex"]], ])))
    autoCharacteristics <- do.call(rbind, lapply(modParamsList, function(stageParams) if(!is.null(stageParams) && !is(stageParams, "PredictParams")) stageParams@characteristics))
    characteristics <- .filterCharacteristics(characteristics, autoCharacteristics)
    characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Cross-validation", value = "Independent Set"))
    
    extras <- unlist(lapply(modParamsList, function(stageParams) if(!is.null(stageParams)) stageParams@otherParams), recursive = FALSE)
    if(length(extras) > 0)
        extras <- extras[sapply(extras, is.atomic)] # Store basic variables, not complex ones.
    extrasDF <- S4Vectors::DataFrame(characteristic = names(extras), value = unname(unlist(extras)))
    characteristics <- rbind(characteristics, extrasDF)
    
    allSamples <- c(rownames(measurementsTrain), rownames(measurementsTest))
    if(!is.null(ncol(outcomeTrain)))
    {
      allOutcome <- rbind(outcomeTrain, outcomeTest)
      rownames(allOutcome) <- allSamples
    } else { 
      allOutcome <- c(outcomeTrain, outcomeTest)
      names(allOutcome) <- allSamples
    }

    predictsTable <- S4Vectors::DataFrame(sample = rownames(measurementsTest), predictedOutcome, check.names = FALSE)
    if(!is(predictedOutcome, "tabular"))
    {
        if(is.factor(predictedOutcome)) colnames(predictsTable)[2] <- "class" else colnames(predictsTable)[2] <- "risk"
    }
    ClassifyResult(characteristics, allSamples, originalFeatures, list(rankedFeatures), list(selectedFeatures),
                   list(models), tuneDetails, predictsTable, allOutcome, importanceTable, modellingParams, list(models))
  }  
})

#' @rdname runTest
#' @export
setMethod("runTest", c("MultiAssayExperiment"),
          function(measurementsTrain, measurementsTest, outcomeColumns, ...)
{
  prepArgsTrain <- list(measurementsTrain, outcomeColumns)
  prepArgsTest <- list(measurementsTest, outcomeColumns)
  extraInputs <- list(...)
  prepExtras <- numeric()
  if(length(extraInputs) > 0)
    prepExtras <- which(names(extrasInputs) %in% .ClassifyRenvir[["prepareDataFormals"]])
  if(length(prepExtras) > 0)
  {      
    prepArgsTrain <- append(prepArgsTrain, extraInputs[prepExtras])
    prepArgsTest <- append(prepArgsTest, extraInputs[prepExtras])
  }
  measurementsAndOutcomeTrain <- do.call(prepareData, prepArgs)
  measurementsAndOutcomeTest <- do.call(prepareData, prepArgs)
  
  runTestArgs <- list(measurementsAndOutcomeTrain[["measurements"]], measurementsAndOutcomeTrain[["outcome"]],
                      measurementsAndOutcomeTest[["measurements"]], measurementsAndOutcomeTest[["outcome"]])
  if(length(extraInputs) > 0 && (length(prepExtras) == 0 || length(extraInputs[-prepExtras]) > 0))
  {
    if(length(prepExtras) == 0) runTestArgs <- append(runTestArgs, extraInputs) else
    runTestArgs <- append(runTestArgs, extraInputs[-prepExtras])
  }
  do.call(runTest, runTestArgs)
})
