#' Split Sample Indexes into Training and Test Partitions for Cross-validation Taking Into Account Classes.
#'
#' \code{samplesSplits} Creates two lists of lists. First has training samples, second has test samples for a range
#'  of different cross-validation schemes.
#'  
#' @aliases samplesSplits samplesSplits,CrossValParams-method samplesSplits,numeric-method
#' @param samplesSplits Default: \code{"k-Fold"}. One of \code{"k-Fold"}, \code{"Permute k-Fold"}, \code{"Permute Percentage Split"}, \code{"Leave-k-Out"}.
#' @param permutations Default: \code{100}. An integer. The number of times the samples are permuted before splitting (repetitions).
#' @param folds Default: \code{5}. An integer. The number of folds to which the samples are partitioned to. Only relevant if \code{samplesSplits} is \code{"k-Fold"} or \code{"Permute k-Fold"}.
#' @param leave Default: \code{2}. An integer. The number of samples to keep for the test set in leave-k-out cross-validation. Only relevant if \code{samplesSplits} is \code{"Leave-k-Out"}.
#' @param percentTest Default: \code{25}. A positive number between 0 and 100. The percentage of samples to keep for the test partition. Only relevant if \code{samplesSplits} is \code{"Permute Percentage Split"}.
#' @param outcome A \code{factor} vector or \code{\link{Surv}} object containing the samples to be partitioned.
#' 
#' @return For \code{samplesSplits}, two lists of the same length. First is training partitions. Second is test partitions.
#' @export
#' @rdname samplesSplitting
#'
#' @examples
#' 
#' classes <- factor(rep(c('A', 'B'), c(15, 5)))
#' splitsList <-samplesSplits(permutations = 1, outcome = classes)
#' splitsList

samplesSplits <- function(samplesSplits = c("k-Fold", "Permute k-Fold", "Permute Percentage Split", "Leave-k-Out"),
                          permutations = 100, folds = 5, percentTest = 25, leave = 2, outcome)
{
  samplesSplits <- match.arg(samplesSplits)
  if(samplesSplits %in% c("k-Fold", "Permute k-Fold"))
  {
    nPermutations <- ifelse(samplesSplits == "k-Fold", 1, permutations)
    nFolds <- folds
    samplesFolds <- lapply(1:nPermutations, function(permutation)
    {
      # Create maximally-balanced folds, so class balance is about the same in all.
      allFolds <- vector(mode = "list", length = nFolds)
      foldsIndexes <- rep(1:nFolds, length.out = length(outcome))
      
      foldsIndex = 1
      # Balance the non-censored observations across folds.
      if(is(outcome, "Surv")) outcome <- factor(outcome[, "status"])
      for(outcomeName in levels(outcome))
      {
        # Permute the indexes of samples in the class.
        whichSamples <- sample(which(outcome == outcomeName))
        whichFolds <- foldsIndexes[foldsIndex:(foldsIndex + length(whichSamples) - 1)]
        
        # Put each sample into its fold.
        for(sampleIndex in 1:length(whichSamples))
        {
          allFolds[[whichFolds[sampleIndex]]] <- c(allFolds[[whichFolds[sampleIndex]]], whichSamples[sampleIndex])
        }
        # Move the beginning index to the first new index.
        foldsIndex <- foldsIndex + length(whichSamples)
      }
      
      list(train = lapply(1:nFolds, function(index) unlist(allFolds[setdiff(1:nFolds, index)])),
           test = allFolds
           )
    })
    # Reorganise into two separate lists, no more nesting.
    list(train = unlist(lapply(samplesFolds, '[[', 1), recursive = FALSE),
         test = unlist(lapply(samplesFolds, '[[', 2), recursive = FALSE))
  } else if(samplesSplits == "Permute Percentage Split") {
    # Take the same percentage of samples from each class to be in training set.
    percent <- percentTest
    samplesTrain <- round((100 - percent) / 100 * table(outcome))
    samplesTest <- round(percent / 100 * table(outcome))
    samplesLists <- lapply(1:permutations, function(permutation)
    {
      trainSet <- unlist(mapply(function(outcomeName, number)
      {
        sample(which(outcome == outcomeName), number)
      }, levels(outcome), samplesTrain))
      testSet <- setdiff(1:length(classes), trainSet)
      list(trainSet, testSet)
    })
    # Reorganise into two lists: training, testing.
    list(train = lapply(samplesLists, "[[", 1), test = lapply(samplesLists, "[[", 2))
  } else if(samplesSplits == "Leave-k-Out") { # leave k out. 
    testSamples <- as.data.frame(utils::combn(length(outcome), leave))
    trainingSamples <- lapply(testSamples, function(sample) setdiff(1:length(outcome), sample))
    list(train = as.list(trainingSamples), test = as.list(testSamples))
  }
}

#' Create a Tabular Representation of Test Set Samples.
#' 
#' \code{splitsTestInfo} creates a table for tracking the permutation, fold number, or subset of each set
#' of test samples. Useful for column-binding to the predictions, once they are unlisted into a vector.
#' 
#' @param samplesSplits Default: \code{"k-Fold"}. One of \code{"k-Fold"}, \code{"Permute k-Fold"}, \code{"Permute Percentage Split"}, \code{"Leave-k-Out"}.
#' @param permutations Default: \code{100}. An integer. The number of times the samples are permuted before splitting (repetitions).
#' @param folds Default: \code{5}. An integer. The number of folds to which the samples are partitioned to. Only relevant if \code{samplesSplits} is \code{"k-Fold"} or \code{"Permute k-Fold"}.
#' @param leave Default: \code{2}. An integer. The number of samples to keep for the test set in leave-k-out cross-validation. Only relevant if \code{samplesSplits} is \code{"Leave-k-Out"}.
#' @param percentTest Default: \code{25}. A positive number between 0 and 100. The percentage of samples to keep for the test partition. Only relevant if \code{samplesSplits} is \code{"Permute Percentage Split"}.
#' @param splitsList The return value of the function \code{samplesSplits}.
#' 
#' @return For \code{splitsTestInfoTable}, a table with a subset of columns \code{"permutation"}, \code{"fold"} and \code{"subset"}, depending on the cross-validation scheme specified.
#' @export
#' @rdname samplesSplitting
#'
#' @examples
#' splitsTestInfo(permutations = 1, splitsList = splitsList)

splitsTestInfo <- function(samplesSplits = c("k-Fold", "Permute k-Fold", "Permute Percentage Split", "Leave-k-Out"),
                                permutations = 100, folds = 5, percentTest = 25, leave = 2, splitsList)
{
  permutationIDs <- NULL
  foldIDs <- NULL
  subsetIDs <- NULL
  samplesSplits <- match.arg(samplesSplits)
  if(samplesSplits %in% c("k-Fold", "Permute k-Fold"))
  {
    foldsSamples <- lengths(splitsList[[2]][1:folds])
    totalSamples <- sum(foldsSamples) 
    if(samplesSplits == "Permute k-Fold")
      permutationIDs <- rep(1:permutations, each = totalSamples)
    times <- ifelse(is.null(permutations), 1, permutations)
    foldIDs <- rep(rep(1:folds, foldsSamples), times = times)
  } else if(samplesSplits == "Permute Percentage Split") {
    permutationIDs <- rep(1:permutations, each = length(splitsList[[2]][[1]]))
  } else { # Leave-k-out
    totalSamples <- length(unique(unlist(splitsList[[2]])))
    subsetIDs <- rep(1:choose(totalSamples, leave), each = leave)
  } 
  summaryTable <- cbind(permutation = permutationIDs, fold = foldIDs, subset = subsetIDs)
  summaryTable
}

# Add extra variables from within runTest functions to function specified by a params object.
.addIntermediates <- function(params)
{
  intermediateName <- params@intermediate
  intermediates <- list(dynGet(intermediateName, inherits = TRUE))
  if(is.null(names(params@intermediate))) names(intermediates) <- intermediateName else names(intermediates) <- names(params@intermediate)
  params@otherParams <- c(params@otherParams, intermediates)
  params
}


# Carries out one iteration of feature selection. Basically, a ranking function is used to rank
# the features in the training set from best to worst and different top sets are used either for
# predicting on the training set (resubstitution) or nested cross-validation of the training set,
# to find the set of top features which give the best (user-specified) performance measure.
.doSelection <- function(measurementsTrain, outcomeTrain, crossValParams, modellingParams, verbose)
{
  tuneParams <- modellingParams@selectParams@tuneParams
  performanceType <- crossValParams@performanceType
  if(!is.null(tuneParams[["nFeatures"]])) topNfeatures <- tuneParams[["nFeatures"]] else topNfeatures <- modellingParams@selectParams@nFeatures
  tuneParams <- tuneParams[-match("nFeatures", names(tuneParams))] # Only used as evaluation metric.
  
  # Make selectParams NULL, since we are currently doing selection and it shouldn't call
  # itself infinitely, but save the parameters of the ranking function for calling the ranking
  # function directly using do.call below.
  featureRanking <- modellingParams@selectParams@featureRanking
  otherParams <- modellingParams@selectParams@otherParams
  doSubset <- modellingParams@selectParams@subsetToSelections 
  modellingParams@selectParams <- NULL
  betterValues <- .ClassifyRenvir[["performanceInfoTable"]][.ClassifyRenvir[["performanceInfoTable"]][, "type"] == performanceType, "better"]
  if(is.function(featureRanking)) # Not a list for ensemble selection.
  {
    paramList <- list(measurementsTrain, outcomeTrain, verbose = verbose)
    paramList <- append(paramList, otherParams) # Used directly by a feature ranking function for rankings of features.
    if(length(tuneParams) == 0) tuneParams <- list(None = "none")
    tuneCombosSelect <- expand.grid(tuneParams, stringsAsFactors = FALSE)

    # Generate feature rankings for each one of the tuning parameter combinations.
    rankings <- lapply(1:nrow(tuneCombosSelect), function(rowIndex)
    {
      tuneCombo <- tuneCombosSelect[rowIndex, , drop = FALSE]
      if(tuneCombo != "none") # Add real parameters before function call.
        paramList <- append(paramList, tuneCombo)
      if(attr(featureRanking, "name") == "randomSelection")
        paramList <- append(paramList, list(nFeatures = topNfeatures))
      do.call(featureRanking, paramList)
    })

    if(attr(featureRanking, "name") %in% c("randomSelection", "previousSelection", "Union Selection")) # Actually selection not ranking.
      return(list(NULL, rankings[[1]], NULL))

    if(crossValParams@tuneMode == "none") # No parameters to choose between.
        return(list(rankings[[1]], rankings[[1]][1:topNfeatures], NULL))

    tuneParamsTrain <- list(topN = topNfeatures)
    tuneParamsTrain <- append(tuneParamsTrain, modellingParams@trainParams@tuneParams)
    tuneCombosTrain <- expand.grid(tuneParamsTrain, stringsAsFactors = FALSE)  
    modellingParams@trainParams@tuneParams <- NULL
    
    allPerformanceTables <- lapply(rankings, function(rankingsVariety)
    {
      # Creates a matrix. Columns are top n features, rows are varieties (one row if None).
      performances <- sapply(1:nrow(tuneCombosTrain), function(rowIndex)
      {
        whichTry <- 1:tuneCombosTrain[rowIndex, "topN"]
        if(doSubset)
        {
          topFeatures <- rankingsVariety[whichTry]
          measurementsTrain <- measurementsTrain[, topFeatures, drop = FALSE] # Features in columns
        } else { # Pass along features to use.
          modellingParams@trainParams@otherParams <- c(modellingParams@trainParams@otherParams, setNames(list(rankingsVariety[whichTry]), names(modellingParams@trainParams@intermediate)))
        }
        if(ncol(tuneCombosTrain) > 1) # There are some parameters for training.
          modellingParams@trainParams@otherParams <- c(modellingParams@trainParams@otherParams, tuneCombosTrain[rowIndex, 2:ncol(tuneCombosTrain), drop = FALSE])
        modellingParams@trainParams@intermediate <- character(0)

        # Do either resubstitution classification or nested-CV classification and calculate the resulting performance metric.
        if(crossValParams@tuneMode == "Resubstitution")
        {
          # Specify measurementsTrain and outcomeTrain for testing, too.
          result <- runTest(measurementsTrain, outcomeTrain, measurementsTrain, outcomeTrain,
                            crossValParams = NULL, modellingParams = modellingParams,
                            verbose = verbose, .iteration = "internal")
          if(is.character(result)) stop(result)
          
          predictions <- result[["predictions"]]
          
          # Classifiers will use a column "class" and survival models will use a column "risk".
          if(class(predictions) == "data.frame")
           predictedOutcome <- predictions[, na.omit(match(c("class", "risk"), colnames(predictions)))]
          else
           predictedOutcome <- predictions
          calcExternalPerformance(outcomeTrain, predictedOutcome, performanceType)
        } else {
           result <- runTests(measurementsTrain, outcomeTrain, crossValParams, modellingParams, verbose = verbose)
           if(is.character(result)[[1]]) stop(result)
           result <- calcCVperformance(result, performanceType)
           median(performance(result)[[performanceType]])
         }
       })

        bestOne <- ifelse(betterValues == "lower", which.min(performances)[1], which.max(performances)[1])
        list(data.frame(tuneCombosTrain, performance = performances), bestOne)
      })

      tablesBestMetrics <- sapply(allPerformanceTables, function(tableIndexPair) tableIndexPair[[1]][tableIndexPair[[2]], "performance"])
      tunePick <- ifelse(betterValues == "lower", which.min(tablesBestMetrics)[1], which.max(tablesBestMetrics)[1])
      
      if(verbose == 3)
         message("Features selected.")

      tuneDetails <- allPerformanceTables[[tunePick]] # List of length 2.
      
      rankingUse <- rankings[[tunePick]]
      selectionIndices <- rankingUse[1:(tuneDetails[[1]][tuneDetails[[2]], "topN"])]
      
      names(tuneDetails) <- c("tuneCombinations", "bestIndex")
      colnames(tuneDetails[[1]])[ncol(tuneDetails[[1]])] <- performanceType
      list(ranked = rankingUse, selected = selectionIndices, tune = tuneDetails)
    } else if(is.list(featureRanking)) { # It is a list of functions for ensemble selection.
      featuresIndiciesLists <- mapply(function(selector, selParams)
      {
        paramList <- list(measurementsTrain, outcomeTrain, trainParams = trainParams,
                          predictParams = predictParams, verbose = verbose)
        paramList <- append(paramList, selParams)
        do.call(selector, paramList)
      }, modellingParams@selectParams@featureRanking, modellingParams@selectParams@otherParams, SIMPLIFY = FALSE)

      performances <- sapply(topNfeatures, function(topN)
      {
        topIndices <- unlist(lapply(featuresIndiciesLists, function(featuresIndicies) featuresIndicies[1:topN]))
        topIndicesCounts <- table(topIndices)
        keep <- names(topIndicesCounts)[topIndicesCounts >= modellingParams@selectParams@minPresence]
        measurementsTrain <- measurementsTrain[, as.numeric(keep), drop = FALSE] # Features in columns
        
        if(crossValParams@tuneMode == "Resubstitution")
        {
          result <- runTest(measurementsTrain, outcomeTrain,
                            measurementsTrain, outcomeTrain,
                            crossValParams = NULL, modellingParams,
                            verbose = verbose, .iteration = "internal")
          if(is.character(result)) stop(result)

          predictions <- result[["predictions"]]
          if(class(predictions) == "data.frame")
            predictedOutcome <- predictions[, "class"]
          else
            predictedOutcome <- predictions
          calcExternalPerformance(outcomeTrain, predictedOutcome, performanceType)
        } else {
          result <- runTests(measurementsSubset, outcomeTrain, crossValParams, modellingParams, verbose = verbose)
          if(is.character(result[[1]])) stop(result)
          result <- calcCVperformance(result, performanceType)
          median(performance(aResult)[[performanceType]])
        }
      })
      bestOne <- ifelse(betterValues == "lower", which.min(performances)[1], which.max(performances)[1])

      selectionIndices <- unlist(lapply(featuresLists, function(featuresList) featuresList[1:topNfeatures[bestOne]]))
      names(table(selectionIndices))[table(selectionIndices) >= modellingParams@selectParams@minPresence]

      list(NULL, selectionIndices, NULL)
    } else { # Previous selection
      selectedFeatures <- list(NULL, selectionIndices, NULL)
    }
}

# Only for transformations that need to be done within cross-validation.
.doTransform <- function(measurementsTrain, measurementsTest, transformParams, verbose)
{
  paramList <- list(measurementsTrain, measurementsTest)
  if(length(transformParams@otherParams) > 0)
    paramList <- c(paramList, transformParams@otherParams)
  paramList <- c(paramList, verbose = verbose)
  do.call(transformParams@transform, paramList)
}

# Code to create a function call to a training function. Might also do training and testing
# within the same function, so test samples are also passed in case they are needed.
.doTrain <- function(measurementsTrain, outcomeTrain, measurementsTest, outcomeTest, crossValParams, modellingParams, verbose)
{
  tuneDetails <- NULL
  if(!is.null(modellingParams@trainParams@tuneParams) && is.null(modellingParams@selectParams))
  {
    performanceType <- crossValParams@performanceType
    tuneCombos <- expand.grid(modellingParams@trainParams@tuneParams, stringsAsFactors = FALSE)
    modellingParams@trainParams@tuneParams <- NULL
    
    performances <- sapply(1:nrow(tuneCombos), function(rowIndex)
    {
      modellingParams@trainParams@otherParams <- c(modellingParams@trainParams@otherParams, as.list(tuneCombos[rowIndex, , drop = FALSE]))
      if(crossValParams@tuneMode == "Resubstitution")
      {
        result <- runTest(measurementsTrain, outcomeTrain, measurementsTrain, outcomeTrain,
                          crossValParams = NULL, modellingParams,
                          verbose = verbose, .iteration = "internal")
        if(is.character(result)) stop(result)

        predictions <- result[["predictions"]]
        if(class(predictions) == "data.frame")
          predictedOutcome <- predictions[, colnames(predictions) %in% c("class", "risk")]
        else
          predictedOutcome <- predictions
        calcExternalPerformance(outcomeTrain, predictedOutcome, performanceType)
      } else {
        result <- runTests(measurementsTrain, outcomeTrain,
                           crossValParams, modellingParams,
                           verbose = verbose, .iteration = "internal")
        if(is.character(result[[1]])) stop(result)
        result <- calcCVperformance(result, performanceType)
        median(performances(result)[[performanceType]])
      }
    })
    allPerformanceTable <- data.frame(tuneCombos, performances)
    colnames(allPerformanceTable)[ncol(allPerformanceTable)] <- performanceType
    
    betterValues <- .ClassifyRenvir[["performanceInfoTable"]][.ClassifyRenvir[["performanceInfoTable"]][, "type"] == performanceType, "better"]
    bestOne <- ifelse(betterValues == "lower", which.min(performances)[1], which.max(performances)[1])
    tuneChosen <- tuneCombos[bestOne, , drop = FALSE]
    tuneDetails <- list(tuneCombos, bestOne)
    names(tuneDetails) <- c("tuneCombinations", "bestIndex")
    modellingParams@trainParams@otherParams <- tuneChosen
  }

    if (!"previousTrained" %in% attr(modellingParams@trainParams@classifier, "name")) 
    # Don't name these first two variables. Some classifier functions might use classesTrain and others use outcomeTrain.
    paramList <- list(measurementsTrain, outcomeTrain)
  else # Don't pass the measurements and classes, because a pre-existing classifier is used.
    paramList <- list()
  if(is.null(modellingParams@predictParams)) # One function does both training and testing.
      paramList <- c(paramList, measurementsTest)
    
  if(length(modellingParams@trainParams@otherParams) > 0)
    paramList <- c(paramList, modellingParams@trainParams@otherParams)
  paramList <- c(paramList, verbose = verbose)
  trained <- do.call(modellingParams@trainParams@classifier, paramList)
  if(verbose >= 2)
    message("Training completed.")  
  
  list(model = trained, tune = tuneDetails)
}

# Creates a function call to a prediction function.
.doTest <- function(trained, measurementsTest, predictParams, verbose)
{
  if(!is.null(predictParams@predictor))
  {
    paramList <- list(trained, measurementsTest)
    if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
    paramList <- c(paramList, verbose = verbose)
    prediction <- do.call(predictParams@predictor, paramList)
  } else { prediction <- trained } # Trained is actually the predictions because only one function, not two.
    if(verbose >= 2)
      message("Prediction completed.")    
    prediction
}

# Converts the characteristics of cross-validation into a pretty character string.
.validationText <- function(crossValParams)
{
  switch(crossValParams@samplesSplits,
  `Permute k-Fold` = paste(crossValParams@permutations, "Permutations,", crossValParams@folds, "Folds"),
  `k-Fold` = paste(crossValParams@folds, "-fold cross-validation", sep = ''),
  `Leave-k-Out` = paste("Leave", crossValParams@leave, "Out"),
  `Permute Percentage Split` = paste(crossValParams@permutations, "Permutations,", crossValParams@percentTest, "% Test"),
  independent = "Independent Set")
}

# Used for ROC area under curve calculation.
# PRtable is a data frame with columns FPR, TPR and class.
# distinctClasses is a vector of all of the class names.
.calcArea <- function(PRtable, distinctClasses)
{
  do.call(rbind, lapply(distinctClasses, function(aClass)
  {
    classTable <- subset(PRtable, class == aClass)
    areaSum <- 0
    for(index in 2:nrow(classTable))
    {
      # Some samples had identical predictions but belong to different classes.
      if(classTable[index, "FPR"] != classTable[index - 1, "FPR"] && classTable[index, "TPR"] != classTable[index - 1, "TPR"])
      {
        newArea <- (classTable[index, "FPR"] - classTable[index - 1, "FPR"]) * classTable[index - 1, "TPR"] + # Rectangle part
         0.5 * (classTable[index, "FPR"] - classTable[index - 1, "FPR"]) * (classTable[index, "TPR"] - classTable[index - 1, "TPR"]) # Triangle part on top.
      } else { # Only one sample with predicted score. Line went either up or right, but not both.
        newArea <- (classTable[index, "FPR"] - classTable[index - 1, "FPR"]) * classTable[index, "TPR"]
      }
      areaSum <- areaSum + newArea
    }
    data.frame(classTable, AUC = round(areaSum, 2), check.names = FALSE)
  }))
}

# Converts features into strings to be displayed in plots.
.getFeaturesStrings <- function(importantFeatures)
{
  # Do nothing if only a simple vector of feature IDs.
  if(!is.null(ncol(importantFeatures[[1]]))) # Data set and feature ID columns.
    importantFeatures <- lapply(importantFeatures, function(features) paste(features[, 1], features[, 2]))
  else if("Pairs" %in% class(importantFeatures[[1]]))
    importantFeatures <- lapply(importantFeatures, function(features) paste(first(features), second(features), sep = '-'))
  importantFeatures
}

# Function to overwrite characteristics which are automatically derived from function names
# by user-specified values.
.filterCharacteristics <- function(characteristics, autoCharacteristics)
{
  # Overwrite automatically-chosen names with user's names.
  if(nrow(autoCharacteristics) > 0 && nrow(characteristics) > 0)
  {
    overwrite <- na.omit(match(characteristics[, "characteristic"], autoCharacteristics[, "characteristic"]))
    if(length(overwrite) > 0)
      autoCharacteristics <- autoCharacteristics[-overwrite, ]
  }
  
  # Merge characteristics tables and return.
  rbind(characteristics, autoCharacteristics)
}

# Don't just plot groups alphabetically, but do so in a meaningful order.
.addUserLevels <- function(plotData, orderingList, metric)
{
  for(orderingIndex in seq_along(orderingList))
  {
    orderingID <- names(orderingList)[orderingIndex]
    ordering <- orderingList[[orderingIndex]]
    if(length(ordering) == 1 && ordering %in% c("performanceAscending", "performanceDescending"))
    { # Order by median values of each group.
      characteristicMedians <- by(plotData[, metric], plotData[, orderingID], median)
      ordering <- names(characteristicMedians)[order(characteristicMedians, decreasing = ordering == "performanceDescending")]
    }
    plotData[, orderingID] <- factor(plotData[, orderingID], levels = ordering)
  }
  plotData
}

# Function to identify the parameters of an S4 method.
.methodFormals <- function(f, signature) {
  tryCatch({
    fdef <- getGeneric(f)
    method <- selectMethod(fdef, signature)
    genFormals <- base::formals(fdef)
    b <- body(method)
    if(is(b, "{") && is(b[[2]], "<-") && identical(b[[2]][[2]], as.name(".local"))) {
      local <- eval(b[[2]][[3]])
      if(is.function(local))
        return(formals(local))
      warning("Expected a .local assignment to be a function. Corrupted method?")
    }
    genFormals
  },
    error = function(error) {
      formals(f)
    })
}

# Find the x-axis positions where a set of density functions cross-over.
# The trivial cross-overs at the beginning and end of the data range are removed.
# Used by the mixtures of normals and naive Bayes classifiers.
.densitiesCrossover <- function(densities) # A list of densities created by splinefun.
{
  if(!all(table(unlist(lapply(densities, function(density) density[['x']]))) == length(densities)))
    stop("x positions are not the same for all of the densities.")
  
  lapply(1:length(densities), function(densityIndex) # All crossing points with other class densities.
  {
    unlist(lapply(setdiff(1:length(densities), densityIndex), function(otherIndex)
    {
      allDifferences <- densities[[densityIndex]][['y']] - densities[[otherIndex]][['y']]
      crosses <- which(diff(sign(allDifferences)) != 0)
      crosses <- sapply(crosses, function(cross) # Refine location for plateaus.
      {
        isSmall <- rle(allDifferences[(cross+1):length(allDifferences)] < 0.000001)
        if(isSmall[["values"]][1] == "TRUE")
          cross <- cross + isSmall[["lengths"]][1] / 2
        cross
      })
      if(length(crosses) > 1 && densities[[densityIndex]][['y']][crosses[1]] < 0.000001 && densities[[densityIndex]][['y']][crosses[length(crosses)]] < 0.000001)
        crosses <- crosses[-c(1, length(crosses))] # Remove crossings at ends of densities.      
      densities[[densityIndex]][['x']][crosses]
    }))
  })
}

# Samples in the training set are upsampled or downsampled so that the class imbalance is
# removed.
.rebalanceTrainingClasses <- function(measurementsTrain, classesTrain, balancing)
{
  samplesPerClassTrain <- table(classesTrain)
  downsampleTo <- min(samplesPerClassTrain)
  upsampleTo <- max(samplesPerClassTrain)
  trainBalanced <- unlist(mapply(function(classSize, className)
  {
    if(balancing == "downsample" && classSize > downsampleTo)
      sample(which(classesTrain == className), downsampleTo)
    else if(balancing == "upsample" && classSize < upsampleTo)
      sample(which(classesTrain == className), upsampleTo, replace = TRUE)
    else
      which(classesTrain == className)
  }, samplesPerClassTrain, names(samplesPerClassTrain), SIMPLIFY = FALSE))
  measurementsTrain <- measurementsTrain[trainBalanced, ]
  classesTrain <- classesTrain[trainBalanced]
  
  list(measurementsTrain = measurementsTrain, classesTrain = classesTrain)
}

.transformKeywordToFunction <- function(keyword)
{
  switch(
        keyword,
        "none" = NULL,
        "diffLoc" = subtractFromLocation
    )
}

.selectionKeywordToFunction <- function(keyword)
{
  switch(
        keyword,
        "none" = NULL,
        "t-test" = differentMeansRanking,
        "limma" = limmaRanking,
        "edgeR" = edgeRranking,
        "Bartlett" = bartlettRanking,
        "Levene" = leveneRanking,
        "DMD" = DMDranking,
        "likelihoodRatio" = likelihoodRatioRanking,
        "KS" = KolmogorovSmirnovRanking,
        "KL" = KullbackLeiblerRanking,
        "CoxPH" = coxphRanking,
        "previousSelection" = previousSelection,
        "randomSelection" = randomSelection,
        "selectMulti" = selectMulti
    )
}

.classifierKeywordToParams <- function(keyword, tuneParams)
{
    switch(
        keyword,
        "randomForest" = RFparams(tuneParams = tuneParams),
        "randomSurvivalForest" = RSFparams(tuneParams = tuneParams),
        "XGB" = XGBparams(tuneParams = tuneParams),
        "GLM" = GLMparams(),
        "ridgeGLM" = ridgeGLMparams(),
        "elasticNetGLM" = elasticNetGLMparams(),
        "LASSOGLM" = LASSOGLMparams(),
        "SVM" = SVMparams(tuneParams = tuneParams),
        "NSC" = NSCparams(),
        "DLDA" = DLDAparams(),
        "naiveBayes" = naiveBayesParams(tuneParams = tuneParams),
        "mixturesNormals" = mixModelsParams(),
        "kNN" = kNNparams(),
        "CoxPH" = coxphParams(),
        "CoxNet" = coxnetParams()
    )    
}

.dlda <- function(x, y, prior = NULL){ # Remove this once sparsediscrim is reinstated to CRAN.
  obj <- list()
  obj$labels <- y
  obj$N <- length(y)
  obj$p <- ncol(x)
  obj$groups <- levels(y)
  obj$num_groups <- nlevels(y)

  est_mean <- "mle"

  # Error Checking
  if (!is.null(prior)) {
    if (length(prior) != obj$num_groups) {
      stop("The number of 'prior' probabilities must match the number of classes in 'y'.")
    }
    if (any(prior <= 0)) {
      stop("The 'prior' probabilities must be nonnegative.")
    }
    if (sum(prior) != 1) {
      stop("The 'prior' probabilities must sum to one.")
    }
  }
  if (any(table(y) < 2)) {
    stop("There must be at least 2 observations in each class.")
  }

  # By default, we estimate the 'a priori' probabilities of class membership with
  # the MLEs (the sample proportions).
  if (is.null(prior)) {
    prior <- as.vector(table(y) / length(y))
  }

  # For each class, we calculate the MLEs (or specified alternative estimators)
  # for each parameter used in the DLDA classifier. The 'est' list contains the
  # estimators for each class.
  obj$est <- tapply(seq_along(y), y, function(i) {
    stats <- list()
    stats$n <- length(i)
    stats$xbar <- colMeans(x[i, , drop = FALSE])
    stats$var <- with(stats, (n - 1) / n * apply(x[i, , drop = FALSE], 2, var))
    stats
  })

  # Calculates the pooled variance across all classes.
  obj$var_pool <- Reduce('+', lapply(obj$est, function(x) x$n * x$var)) / obj$N

  # Add each element in 'prior' to the corresponding obj$est$prior
  for(k in seq_len(obj$num_groups)) {
    obj$est[[k]]$prior <- prior[k]
  }
  class(obj) <- "dlda"
  obj
}

#' @method predict dlda
predict.dlda <- function(object, newdata, ...) { # Remove once sparsediscrim is reinstated to CRAN.
  if (!inherits(object, "dlda"))  {
    stop("object not of class 'dlda'")
  }
  if (is.vector(newdata)) {
    newdata <- as.matrix(newdata)
  }

  scores <- apply(newdata, 1, function(obs) {
    sapply(object$est, function(class_est) {
      with(class_est, sum((obs - xbar)^2 / object$var_pool) + log(prior))
    })
  })

  if (is.vector(scores)) {
    min_scores <- which.min(scores)
  } else {
    min_scores <- apply(scores, 2, which.min)
  }

  # Posterior probabilities via Bayes Theorem
  means <- lapply(object$est, "[[", "xbar")
  covs <- replicate(n=object$num_groups, object$var_pool, simplify=FALSE)
  priors <- lapply(object$est, "[[", "prior")
  posterior <- .posterior_probs(x=newdata,
                               means=means,
                               covs=covs,
                               priors=priors)

  class <- factor(object$groups[min_scores], levels = object$groups)

  list(class = class, scores = scores, posterior = posterior)
}

.posterior_probs <- function(x, means, covs, priors) { # Remove once sparsediscrim is reinstated to CRAN.
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  x <- as.matrix(x)

  posterior <- mapply(function(xbar_k, cov_k, prior_k) {
    if (is.vector(cov_k)) {
      post_k <- apply(x, 1, function(obs) {
        .dmvnorm_diag(x=obs, mean=xbar_k, sigma=cov_k)
      })
    } else {
      post_k <- dmvnorm(x=x, mean=xbar_k, sigma=cov_k)
    }
    prior_k * post_k
  }, means, covs, priors)

  if (is.vector(posterior)) {
    posterior <- posterior / sum(posterior)
    posterior <- matrix(posterior, nrow = 1) # Ensure it's always matrix, like just below.
    colnames(posterior) <- names(priors)
  } else {
    posterior <- posterior / rowSums(posterior)
  }
  posterior
}

.dmvnorm_diag <- function(x, mean, sigma) { # Remove once sparsediscrim is reinstated to CRAN.
  exp(sum(dnorm(x, mean=mean, sd=sqrt(sigma), log=TRUE)))
}

# Function to create permutations of a vector, with the possibility to restrict values at certain positions.
# fixed parameter is a data frame with first column position and second column value.
.permutations <- function(data, fixed = NULL)
{
  items <- length(data)
  multipliedTo1 <- factorial(items)
  if(items > 1) 
    permutations <- structure(vapply(seq_along(data), function(index)
                     rbind(data[index], .permutations(data[-index])), 
                     data[rep(1L, multipliedTo1)]), dim = c(items, multipliedTo1))
  else permutations <- data
  
  if(!is.null(fixed))
  {
    for(rowIndex in seq_len(nrow(fixed)))
    {
      keepColumns <- permutations[fixed[rowIndex, 1], ] == fixed[rowIndex, 2]
      permutations <- permutations[, keepColumns]
    }
  }
  permutations
}