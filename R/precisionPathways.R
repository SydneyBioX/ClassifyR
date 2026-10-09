#' Precision Pathways for Sample Prediction Based on Prediction Confidence.
#' 
#' Precision pathways allows the evaluation of various permutations of multiomics or multiview data.
#' Samples are predicted by a particular assay if they were consistently predicted as a particular class
#' during cross-validation. Otherwise, they are passed onto subsequent assays/tiers for prediction. Balanced accuracy
#' is used to evaluate overall prediction performance and sample-specific accuracy for individual-level evaluation.
#'
#' @param measurements Either a \code{\link{MultiAssayExperiment}} or a list of the basic tabular objects containing the data.
#' @param class If a \code{\link{MultiAssayExperiment}}, a column name in \code{colData(measurements)} with the classes. If \code{measurements} is a \code{list} of tabular data, may also be
#' a vector of classes. Optional for \code{precisionPathwaysPredict}: the classes of the new samples are only needed to calculate accuracies.
#' @param useFeatures Default: \code{NULL} (i.e. use all provided features). A named list of features to use. Otherwise, the input data is a single table and this can just be a vector of feature names.
#' For any assays not in the named list, all of their features are used. \code{"clinical"} is also a valid assay name and refers to the clinical data table.
#' This allows for the avoidance of variables such spike-in RNAs, sample IDs, sample acquisition dates, etc. which are not relevant for outcome prediction.
#' @param maxMissingProp Default: 0.0. A proportion less than 1 which is the maximum
#' tolerated proportion of missingness for a feature to be retained for modelling.
#' @param topNvariance Default: NULL. An integer number of most variable features per assay to subset to.
#' Assays with less features won't be reduced in size.
#' @param fixedAssays A character vector of assay names specifying any assays which must be at the
#' beginning of the pathway, in the order given. \code{NULL} allows any assay to be first.
#' @param confidenceCutoff The minimum confidence of predictions for a sample to be predicted by a particular assay. If a sample
#' was predicted to belong to a particular class a proportion \eqn{p} times, then the confidence is \eqn{2 \times |p - 0.5|}.
#' A sample whose confidence is equal to or greater than \code{confidenceCutoff} is predicted at that tier.
#' @param minAssaySamples An integer specifying the minimum number of samples a tier may have. If a subsequent tier
#' would have less than this number of samples, the samples are incorporated into the current tier.
#' @param mode Default: \code{"stability"}. Either \code{"stability"} or \code{"combinatorial"}. If \code{"stability"}, then pathways grow by
#' passing samples with confidence below \code{confidenceCutoff} to the next assay, until the assays are exhausted or the samples are.
#' If \code{"combinatorial"}, then all possible combinations of assays respecting \code{fixedAssays}, as well as each assay individually
#' are considered.
#' @param nFeatures Default: 20. The number of features to consider during feature selection, if feature selection is done.
#' @param selectionMethod A named character vector of feature selection methods to use for the assays, one for each. The names must correspond to names of \code{measurements}.
#' @param classifier A named character vector of modelling methods to use for the assays, one for each. The names must correspond to names of \code{measurements}.
#' @param nFolds A numeric specifying the number of folds to use for cross-validation.
#' @param nRepeats A numeric specifying the the number of repeats or permutations to use for cross-validation.
#' @param nCores A numeric specifying the number of cores used if the user wants to use parallelisation.
#' @param pathways A set of pathways created by \code{precisionPathwaysTrain} which is an object of class \code{PrecisionPathways} to be used for predicting on a new data set.
#' @details \code{precisionPathwaysTrain} cross-validates a classifier for each assay. The confidence of a sample for an assay is
#' calculated from its cross-validated predictions. \code{precisionPathwaysPredict} predicts each new sample with every model fitted
#' during that cross-validation (one model per fold of every repeat), and the confidence of a new sample for an assay is calculated from
#' those predictions. Each new sample is predicted independently of the classes of the new samples; the stored models are used as they are.
#' The pathways are then formed in the same way as in training. Each pathway records the assays in the order in which samples passed through
#' them, in its element named \code{"assays"}.
#' @rdname precisionPathways
#' @aliases precisionPathwaysTrain precisionPathwaysPredict
#' @return An object of class \code{PrecisionPathways} which is basically a named list that other plotting and
#' tabulating functions can use.
#' @examples
#' # To be determined.

#' @usage NULL
setGeneric("precisionPathwaysTrain", function(measurements, class, ...)
    standardGeneric("precisionPathwaysTrain"))

#' @rdname precisionPathways
#' @export
setMethod("precisionPathwaysTrain", "MultiAssayExperimentOrList",
          function(measurements, class, useFeatures = NULL, maxMissingProp = 0.0, topNvariance = NULL,
                   fixedAssays = "clinical", confidenceCutoff = 0.8, minAssaySamples = 10, mode = c("stability", "combinatorial"),
                   nFeatures = 20, selectionMethod = setNames(c("none", rep("t-test", length(measurements))), c("clinical", names(measurements))),
                   classifier = setNames(c("elasticNetGLM", rep("randomForest", length(measurements))), c("clinical", names(measurements))),
                   nFolds = 5, nRepeats = 20, nCores = 1)
          {
            if(is.list(measurements)) # Ensure plain list has clinical data.
            {
              # One of the tables must be named "clinical".
              if (!any(names(measurements) == "clinical"))
                stop("One of the tables must be named \"clinical\".")
              clinicalColumns <- colnames(measurements[["clinical"]])
            } else {
              clinicalColumns <- colnames(MultiAssayExperiment::colData(measurements))
            }
            # A list input already has a table named "clinical", so the default vectors name it twice.
            selectionMethod <- selectionMethod[!duplicated(names(selectionMethod))]
            classifier <- classifier[!duplicated(names(classifier))]

            if(is.null(useFeatures) || !"clinical" %in% names(useFeatures))
            {
              warning("No 'useFeatures' named list element for clinical data is specified. Clinical data often has\n",
    "lots of uninformative variables. Please consider specifying useful features.")
              if(is.character(class)) clinicalColumns <- setdiff(clinicalColumns, class)
              if(is.null(useFeatures)) useFeatures <- list()
              useFeatures[["clinical"]] <- clinicalColumns
            }
            mode <- match.arg(mode)

            prepArgs <- list(measurements, class, useFeatures = useFeatures,
                             maxMissingProp = maxMissingProp, topNvariance = topNvariance)
            measurementsAndClass <- do.call(prepareData, prepArgs)

            .precisionPathwaysTrain(measurementsAndClass[["measurements"]], measurementsAndClass[["outcome"]],
                                   useFeatures = useFeatures, fixedAssays = fixedAssays, confidenceCutoff = confidenceCutoff,
                                   minAssaySamples = minAssaySamples, mode = mode, nFeatures = nFeatures,
                                   selectionMethod = selectionMethod, classifier = classifier,
                                   nFolds = nFolds, nRepeats = nRepeats, nCores = nCores)
          })

# Internal method which carries out all of the processing, obtaining reformatted data from the
# MultiAssayExperiment and list (of basic rectangular tables) S4 methods.
.precisionPathwaysTrain <- function(measurements, class, fixedAssays = "clinical",
                   useFeatures = useFeatures, confidenceCutoff = 0.8, minAssaySamples = 10, mode,
                   nFeatures = 20, selectionMethod = setNames(c("none", rep("t-test", length(setdiff(unique(mcols(measurements)$assay), "clinical")))), c("clinical", setdiff(unique(mcols(measurements)$assay), "clinical"))),
                   classifier = setNames(c("elasticNetGLM", rep("randomForest", length(setdiff(unique(mcols(measurements)$assay), "clinical")))), c("clinical", setdiff(unique(mcols(measurements)$assay), "clinical"))),
                   nFolds = 5, nRepeats = 20, nCores = 1)
          {
            # Step 1: Determine all valid permutations of assays, taking into account the
            # assays to be used and which assays, if any, must be included.
            assayIDs <- unique(S4Vectors::mcols(measurements)[["assay"]])
            if(!all(fixedAssays %in% assayIDs))
              stop("'fixedAssays' ", paste(setdiff(fixedAssays, assayIDs), collapse = ", "), " not found in the data. Assays are: ",
                   paste(assayIDs, collapse = ", "), '.')
            fixedPositions <- NULL
            if(length(fixedAssays) > 0) fixedPositions <- data.frame(seq_along(fixedAssays), fixedAssays)
            assaysPermutations <- .permutations(assayIDs, fixed = fixedPositions)
            assaysPermutations <- lapply(asplit(assaysPermutations, 2), as.vector)
            if(mode == "combinatorial")
            {
                subsetPermutations <- list()
                variableAssays <- setdiff(assayIDs, fixedAssays)
                nVariable <- 1
                while((length(fixedAssays) + nVariable) < length(fixedAssays) + length(variableAssays))
                {
                    subsetCombinations <- combn(variableAssays, nVariable)
                    subsetPermutations <- c(subsetPermutations, apply(subsetCombinations, 2, function(combination)
                    {
                        permutations <- .permutations(combination, NULL)
                        if(is.matrix(permutations)) permutations <- asplit(permutations, 2)
                        lapply(permutations, function(combinationPermutation) c(fixedAssays, combinationPermutation))
                    }))
                    nVariable <- nVariable + 1
                }
                subsetPermutations <- unlist(subsetPermutations, recursive = FALSE)
                assaysPermutations <- c(assaysPermutations, subsetPermutations) # All subsets.
                assaysPermutations <- c(assaysPermutations, assayIDs) # Each assay on its own.
            }

            # Step 2: Build a classifier for each assay using all of the samples.
            modelsList <- crossValidate(measurements, class, nFeatures, selectionMethod,
                                        classifier = classifier, nFolds = nFolds,
                                        nRepeats = nRepeats, nCores = nCores)
            if(is(modelsList, "ClassifyResult")) modelsList <- list(modelsList)
            modelsList <- lapply(modelsList, calcCVperformance, "Sample Accuracy") # Add sample accuracy, which can be subset later.
            names(modelsList) <- .pathwaysAssayNames(modelsList, assayIDs)
            modelsList <- lapply(modelsList, .pathwaysTrimModels)

            # Step 3: Loop over each pathway and each assay in order to determine which samples are used at that level
            # and which are passed onwards.
            predictionsSamplesCounts <- lapply(modelsList, function(result)
                                        {
                                          allPredictions <- predictions(result)
                                          table(allPredictions[, "sample"], allPredictions[, "class"])
                                        })
            sampleAccuracies <- lapply(modelsList, function(result) performance(result)[["Sample Accuracy"]])
            allSampleIDs <- sampleNames(modelsList[[1]])
            knownClasses <- actualOutcome(modelsList[[1]])
            precisionPathways <- lapply(assaysPermutations, .pathwayTiers, predictionsSamplesCounts, sampleAccuracies,
                                        allSampleIDs, knownClasses, confidenceCutoff, minAssaySamples)

            names(precisionPathways) <- sapply(precisionPathways, "[[", "pathway")
            attr(modelsList, "crossValParams") <- generateCrossValParams(nRepeats, nFolds, nCores, NULL)
            precisionPathways <- precisionPathways[unique(names(precisionPathways))] # In case early termination causes duplicates.
            result <- list(models = modelsList, assaysPermutations = assaysPermutations,
                           parameters = list(confidenceCutoff = confidenceCutoff, minAssaySamples = minAssaySamples, classifier = classifier),
                           useFeatures = useFeatures, pathways = precisionPathways)
            class(result) <- "PrecisionPathways"

            result
}

# The assay name of each cross-validation result, so that models are always looked up by assay and never by position.
.pathwaysAssayNames <- function(modelsList, assayIDs)
{
  assayNames <- sapply(modelsList, function(result)
  {
    characteristics <- result@characteristics
    assayName <- characteristics[characteristics[, "characteristic"] == "Assay Name", "value"]
    if(length(assayName) == 0) NA_character_ else assayName
  })
  if(anyNA(assayNames) && length(modelsList) == length(assayIDs)) # Results lacking an assay name are in data order.
    assayNames[is.na(assayNames)] <- assayIDs[is.na(assayNames)]
  assayNames
}

# Prediction of new samples uses only the fitted models. The extra forests that some classifiers
# store for calculating variable importance are not needed.
.pathwaysTrimModels <- function(result)
{
  result@models <- lapply(result@models, function(model)
  {
    attr(model, "forImportance") <- NULL
    model
  })
  result
}

# Confidence of each sample, from a table of counts of predicted classes (samples in rows, classes in columns).
# A small tolerance ensures that a confidence equal to the cut-off is not missed because of floating point error.
.pathwaysConfident <- function(predictionsSamplesCounts, confidenceCutoff)
{
  confidences <- 2 * abs(predictionsSamplesCounts[, 1] / rowSums(predictionsSamplesCounts) - 0.5)
  names(confidences)[confidences >= confidenceCutoff - sqrt(.Machine$double.eps)]
}

# Balanced accuracy, with predicted and known classes using the same factor levels in the same order.
.pathwaysBalancedAccuracy <- function(knownClasses, predictedClasses)
{
  if(is.null(knownClasses) || all(is.na(knownClasses))) return(NA_real_)
  classLevels <- union(levels(predictedClasses), levels(knownClasses))
  calcExternalPerformance(factor(as.character(knownClasses), levels = classLevels),
                          factor(as.character(predictedClasses), levels = classLevels))
}

# Pass samples along one sequence of assays. At each tier, the samples which have not yet been predicted and
# are predicted with sufficient confidence are predicted by that tier. Used for training and for prediction.
.pathwayTiers <- function(permutation, predictionsSamplesCounts, sampleAccuracies, allSampleIDs, knownClasses,
                          confidenceCutoff, minAssaySamples)
{
  assaysProcessed <- character()
  assaysTested <- character()
  samplesUsed <- character()
  individualsTableAll <- S4Vectors::DataFrame()
  tierTableAll <- S4Vectors::DataFrame()
  breakEarly = FALSE
  for(assay in permutation)
  {
    # Identify all samples which are consistently predicted.
    predictionsSamplesCountsUse <- predictionsSamplesCounts[[assay]]
    sampleIDsUse <- .pathwaysConfident(predictionsSamplesCountsUse, confidenceCutoff)
    sampleIDsUse <- setdiff(sampleIDsUse, samplesUsed)
    assaysTested <- c(assaysTested, assay) # Every sample not yet predicted has this assay done.

    # Check if too few samples left for next round. Include them in this round, if so.
    remainingIDs <- setdiff(allSampleIDs, c(samplesUsed, sampleIDsUse))
    if(length(remainingIDs) < minAssaySamples || assay == permutation[length(permutation)])
    {
      sampleIDsUse <- c(sampleIDsUse, remainingIDs)
      breakEarly = TRUE
    }
    if(length(sampleIDsUse) == 0) next # No samples have sufficient confidence at this tier.

    predictionsSamplesCountsUse <- predictionsSamplesCountsUse[sampleIDsUse, , drop = FALSE]

    # Individuals predictions and sample-wise accuracy, tier-wise error.
    maxVotes <- apply(predictionsSamplesCountsUse, 1, function(sample) which.max(sample))
    predictedClasses <- factor(colnames(predictionsSamplesCountsUse)[maxVotes],
                               levels = colnames(predictionsSamplesCountsUse))
    individualsTable <- S4Vectors::DataFrame(Tier = assay,
                                             `Sample ID` = sampleIDsUse,
                                             `Predicted` = predictedClasses,
                                             `Accuracy` = sampleAccuracies[[assay]][sampleIDsUse],
                                             check.names = FALSE)
    balancedAccuracy <- .pathwaysBalancedAccuracy(knownClasses[match(sampleIDsUse, allSampleIDs)], predictedClasses)
    tierTable <- S4Vectors::DataFrame(Tier = assay,
                                      `Balanced Accuracy` = balancedAccuracy, check.names = FALSE)

    assaysProcessed <- c(assaysProcessed, assay)
    individualsTableAll <- rbind(individualsTableAll, individualsTable)
    tierTableAll <- rbind(tierTableAll, tierTable)
    samplesUsed <- c(samplesUsed, sampleIDsUse)

    if(breakEarly == TRUE) break
  }
  pathwayString <- paste(assaysProcessed, collapse = '-')
  individualsTableAll[, "Tier"] <- factor(individualsTableAll[, "Tier"], levels = permutation)
  list(pathway = pathwayString, assays = assaysProcessed, assaysTested = assaysTested,
       individuals = individualsTableAll, tiers = tierTableAll)
}

# A nice print method to avoid flooding the screen with lots of tables
# when result is shown in console.
#' @export
print.PrecisionPathways <- function(x, ...)
{
  cat("An object of class 'PrecisionPathways'.\n")
  cat("Pathways:\n")
  cat(paste(names(x[["pathways"]]), collapse = '\n'))
  cat("\n")
  invisible(x)
}

#' @usage NULL
setGeneric("precisionPathwaysPredict", function(pathways, measurements, class, ...)
    standardGeneric("precisionPathwaysPredict"))

#' @rdname precisionPathways
#' @export
setMethod("precisionPathwaysPredict", c("PrecisionPathways", "MultiAssayExperimentOrList"),
          function(pathways, measurements, class)
          {
            if(missing(class)) class <- NULL # Classes of new samples are only used for calculating accuracy.
            if(is.list(measurements)) # Ensure plain list has clinical data.
            {
              # One of the tables must be named "clinical".
              if (!any(names(measurements) == "clinical"))
                stop("One of the tables must be named \"clinical\".")
            }

            if(is.null(class) && is(measurements, "MultiAssayExperiment"))
            { # prepareData needs an outcome column. Sample identifiers stand in for it and are then discarded.
              class <- ".pathwaysSampleID"
              MultiAssayExperiment::colData(measurements)[, class] <- rownames(MultiAssayExperiment::colData(measurements))
              measurementsAndClass <- prepareData(measurements, class, useFeatures = pathways[["useFeatures"]])
              measurementsAndClass["outcome"] <- list(NULL)
            } else {
              prepArgs <- list(measurements, class, useFeatures = pathways[["useFeatures"]])
              measurementsAndClass <- do.call(prepareData, prepArgs)
            }

            .precisionPathwaysPredict(pathways, measurementsAndClass[["measurements"]], measurementsAndClass[["outcome"]])
          })

.precisionPathwaysPredict <- function(pathways, measurements, class)
{
  # Step 1: Extract all of previously fitted models and pathways. Predict the new samples with every fitted model of each assay.
  modelsList <- pathways[["models"]]
  names(modelsList) <- .pathwaysAssayNames(modelsList, names(modelsList))
  assaysPermutations <- lapply(pathways[["pathways"]], function(pathway)
                        {
                          if(!is.null(pathway[["assays"]])) pathway[["assays"]] else as.character(pathway[["tiers"]][, "Tier"])
                        })
  confidenceCutoff <- pathways[["parameters"]][["confidenceCutoff"]]
  minAssaySamples <- pathways[["parameters"]][["minAssaySamples"]]
  classifiers <- pathways[["parameters"]][["classifier"]]
  allSampleIDs <- rownames(measurements)

  resultsNew <- mapply(function(fittedModel, assayID)
  {
    if(!is.null(fittedModel@modellingParams))
      predictParams <- fittedModel@modellingParams@predictParams
    else # Result created without its modelling parameters.
      predictParams <- PredictParams(classifiers[[assayID]])
    if(is.null(predictParams))
      stop("The classifier of assay ", assayID, " trains and predicts in one step, so its models cannot predict new samples.")
    classLevels <- levels(actualOutcome(fittedModel))
    predictedClasses <- lapply(models(fittedModel), function(model)
    {
      measurementsUse <- measurements[, attr(model, "featuresForTrain"), drop = FALSE]
      paramList <- c(list(model, measurementsUse), predictParams@otherParams, verbose = 0)
      predicted <- do.call(predictParams@predictor, paramList)
      if(!is.factor(predicted)) predicted <- predicted[, "class"]
      as.character(predicted)
    })
    predictsNew <- S4Vectors::DataFrame(sample = rep(allSampleIDs, length(predictedClasses)),
                                        class = factor(unlist(predictedClasses), levels = classLevels))
    predictionsSamplesCounts <- table(predictsNew[, "sample"], predictsNew[, "class"])
    sampleAccuracies <- setNames(rep(NA_real_, nrow(predictionsSamplesCounts)), rownames(predictionsSamplesCounts))
    if(!is.null(class))
    {
      knownClasses <- as.character(class)[match(rownames(predictionsSamplesCounts), allSampleIDs)]
      correctColumn <- match(knownClasses, colnames(predictionsSamplesCounts))
      correctCounts <- predictionsSamplesCounts[cbind(seq_len(nrow(predictionsSamplesCounts)), correctColumn)]
      correctCounts[is.na(correctCounts)] <- 0 # Known class never predicted by the models.
      sampleAccuracies[] <- correctCounts / rowSums(predictionsSamplesCounts)
    }
    list(predictionsSamplesCounts, sampleAccuracies)
  }, modelsList, names(modelsList), SIMPLIFY = FALSE)
  predictionsSamplesCounts <- lapply(resultsNew, "[[", 1)
  sampleAccuracies <- lapply(resultsNew, "[[", 2)

  # Step 2: Loop over each pathway and each assay in order to determine which samples are used at that level
  # and which are passed onwards.
  precisionPathways <- lapply(assaysPermutations, .pathwayTiers, predictionsSamplesCounts, sampleAccuracies,
                              allSampleIDs, class, confidenceCutoff, minAssaySamples)
  names(precisionPathways) <- sapply(precisionPathways, "[[", "pathway")
  precisionPathways <- precisionPathways[unique(names(precisionPathways))]
  if(is.null(class)) class <- factor(rep(NA, length(allSampleIDs)), levels = levels(actualOutcome(modelsList[[1]])))
  result <- list(models = modelsList, assaysPermutations = unname(assaysPermutations),
                 parameters = list(confidenceCutoff = confidenceCutoff, minAssaySamples = minAssaySamples, classifier = classifiers),
                 pathways = precisionPathways, testSampleInfo = DataFrame(sampleID = allSampleIDs, class = class))
  class(result) <- "PrecisionPathways"

  result
}

# Calculate accuracy and costs of each pathway.

#' Various Functions for Evaluating Precision Pathways
#'
#' These functions tabulate or plot various aspects of precision pathways, such as accuracies and costs.
#'
#' The cost of a pathway is the total cost of the assays done on all of the samples. A sample predicted at a particular tier
#' had every assay of the pathway up to and including that tier done, so it is charged the cost of each of those assays.
#' For example, a sample predicted by the second assay of the pathway clinical-RNA costs the price of clinical data plus the price of RNA.
#'
#' @param precisionPathways A pathway of class \code{PrecisionPathways}.
#' @param costs A named vector of assays with the cost of each one. An assay that is not named costs nothing.
#' @rdname precisionPathwaysEvaluations
#' @export
calcCostsAndPerformance <- function(precisionPathways, costs = NULL)
{
  if(is.null(costs))
    stop("'costs' of each assay must be specified.")

  if("testSampleInfo" %in% names(precisionPathways))
  {
      knownClasses <- precisionPathways[["testSampleInfo"]][, "class"]
      allNames <- precisionPathways[["testSampleInfo"]][, "sampleID"]
  } else{
      knownClasses <- actualOutcome(precisionPathways$models[[1]])
      allNames <- sampleNames(precisionPathways$models[[1]])
  }
  accuraciesCosts <- do.call(rbind, lapply(precisionPathways[["pathways"]], function(pathway)
  {
    predictions <- pathway[["individuals"]][, "Predicted"]
    knownClasses <- knownClasses[match(pathway[["individuals"]][, "Sample ID"], allNames)]
    balancedAccuracy <- .pathwaysBalancedAccuracy(knownClasses, predictions)

    # Each sample is charged for every assay it passed through, up to and including the tier that predicted it.
    assaysTested <- pathway[["assaysTested"]]
    if(is.null(assaysTested)) assaysTested <- as.character(pathway[["tiers"]][, "Tier"])
    assaysCosts <- unname(costs[match(assaysTested, names(costs))])
    assaysCosts[is.na(assaysCosts)] <- 0
    cumulativeCosts <- setNames(cumsum(assaysCosts), assaysTested)
    costTotal <- sum(cumulativeCosts[as.character(pathway[["individuals"]][, "Tier"])])

    data.frame(accuracy = round(balancedAccuracy, 2), cost = costTotal)
  }))

  precisionPathways$performance <- accuraciesCosts
  precisionPathways
}

# Print a summary table, including accuracy and costs.

#' @param object A set of pathways of class \code{PrecisionPathways}.
#' @param weights A numeric vector of length two named \code{"accuracy"} and \code{"cost"} specifying how to weight the predictive accuracy
#' and the cost during ranking. Must sum to 1. An unnamed vector is taken to be in the order accuracy, cost.
#' @param ... Not used but just following the S3 requirement of the generic template.
#' @rdname precisionPathwaysEvaluations
#' @export
summary.PrecisionPathways <- function(object, weights = c(accuracy = 0.5, cost = 0.5), ...)
{
  if(!"performance" %in% names(object))
    stop("Nothing to summarise. Please run function 'calcCostsAndPerformance' first.")
  if(is.null(names(weights))) names(weights) <- c("accuracy", "cost")[seq_along(weights)]
  if(length(weights) > 2 || anyNA(names(weights)) || !all(names(weights) %in% c("accuracy", "cost")) || anyDuplicated(names(weights)))
    stop("'weights' must have one element named \"accuracy\" and one named \"cost\".")
  if(abs(sum(weights) - 1) > sqrt(.Machine$double.eps))
    stop("'weights' must sum to 1.")
  weightOf <- function(name) if(name %in% names(weights)) weights[[name]] else 0

  hasCost <- "cost" %in% names(object[["performance"]])
  summaryTable <- data.frame(Pathway = rownames(object[["performance"]]),
                             `Balanced Accuracy` = object[["performance"]][, "accuracy"],
                              check.names = FALSE)
  if(hasCost) summaryTable[, "Total Cost"] <- object[["performance"]][, "cost"]

  finalScores <- rank(object[["performance"]][, "accuracy"], na.last = "keep") * weightOf("accuracy")
  if(hasCost)
    finalScores <- finalScores + rank(-object[["performance"]][, "cost"], na.last = "keep") * weightOf("cost")
  summaryTable <- cbind(summaryTable, Score = finalScores)
  summaryTable <- summaryTable[order(summaryTable[, "Score"], decreasing = TRUE), ]
  summaryTable
}

#' @export
#' @rdname precisionPathwaysEvaluations
bubblePlot <- function (precisionPathways, ...) {
   UseMethod("bubblePlot", precisionPathways)
 }

#' @param precisionPathways A pathway of class \code{PrecisionPathways}.
#' @param pathwayColours A named vector of colours with names being the names of pathways. If none is specified,
#' a default colour scheme will automatically be chosen.
#' @rdname precisionPathwaysEvaluations
#' @export
bubblePlot.PrecisionPathways <- function(precisionPathways, pathwayColours = NULL, ...)
{
  if(is.null(pathwayColours)) pathwayColours <- scales::hue_pal()(length(precisionPathways[["pathways"]]))
  performance <- precisionPathways[["performance"]]
  performance <- cbind(Sequence = rownames(performance), performance)
  ggplot2::ggplot(performance, aes(x = accuracy, y = cost, colour = Sequence, size = 4)) + ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = pathwayColours) + ggplot2::labs(x = "Balanced Accuracy", y = "Total Cost") +
    ggplot2::guides(size = "none") + ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA))
}


#' @export
#' @rdname precisionPathwaysEvaluations
flowchart <- function (precisionPathways, ...) {
   UseMethod("flowchart", precisionPathways)
 }

#' @param precisionPathways A pathway of class \code{PrecisionPathways}.
#' @param pathway A character vector of length 1 specifying which pathway to plot, e.g. "clinical-mRNA".
#' @param orientation Default: \code{"horizontal"}. Either \code{"horizontal"} or \code{"vertical"}. Specifies the layout of the flowchart.
#' @param nodeColours A named vector of colours with names being \code{"assay"}, \code{"class1"},\code{"class2"}.
#' a default colour scheme will automatically be chosen.
#' @rdname precisionPathwaysEvaluations
#' @export
flowchart.PrecisionPathways <- function(precisionPathways, pathway,
                                        orientation = c("horizontal", "vertical"),
                                        nodeColours = c(assay = "snow3", class1 = "#3F48CC", class2 = "#880015"), ...)
{
  if(!requireNamespace("data.tree", quietly = TRUE))
    stop("The package 'data.tree' could not be found. Please install it.")
  orientation <- match.arg(orientation)
  orientation <- ifelse(orientation == "horizontal", "LR", "TB")

  pathwayUse <- precisionPathways[["pathways"]][[pathway]]
  assayIDs <- pathwayUse[["tiers"]][, 1]
  possibleClasses <- levels(precisionPathways$models[[1]]@actualOutcome)
  samplesTiers <- pathwayUse$individuals
      
  pathwayTree <- data.tree::Node$new(assayIDs[1])
  currentNode <- pathwayTree
  for(assay in assayIDs)
  {
    toDo <- nrow(samplesTiers) - max(which(samplesTiers[, "Tier"] == assay))
    class1Predictions <- currentNode$AddChild(possibleClasses[1], counter = nrow(subset(samplesTiers, Predicted == possibleClasses[1] & Tier == assay)), nodeType = "Class1")
    if(assay != assayIDs[length(assayIDs)])
      uncertain <- currentNode$AddChild("Uncertain", counter = toDo, nodeType = "Uncertain")
    class2Predictions = currentNode$AddChild(possibleClasses[2], counter = nrow(subset(samplesTiers, Predicted == possibleClasses[2] & Tier == assay)), nodeType = "Class2")
    
    if(toDo == 0) {break} else {
      currentNode <- uncertain
      currentNode <- currentNode$AddChild(assayIDs[match(assay, assayIDs) + 1], nodeType = "Platform")
    }
  }
  
  .getFillColour <- function(node) {
    if(node$isRoot || node$nodeType == "Platform"){
      colour <- nodeColours[["assay"]]
    } else if(node$nodeType == "Class1"){
      colour <- nodeColours[["class1"]]
    } else if(node$nodeType == "Class2"){
      colour <- nodeColours[["class2"]]
    } else {
      colour <- "snow3"
    }
    return(colour)
  }
  
  .getEdgeLabel <- function(node)
  {
    nSamples <- nrow(samplesTiers)
    if(node$isRoot || node$nodeType == "Platform")
    {
      label <- NULL
    } else {
      value <- round((node$counter / nSamples) * 100)
      label <- paste(value,  "% (", node$counter, ")", sep = '')
    }
    return(label)
  }
          
  data.tree::SetGraphStyle(pathwayTree, rankdir = orientation)
  data.tree::SetEdgeStyle(pathwayTree, fontname = 'helvetica', label = .getEdgeLabel)
  data.tree::SetNodeStyle(pathwayTree, style = "filled", shape = .getNodeShape, fontcolor = "black", fillcolor = .getFillColour, fontname = 'helvetica')
  pathwayTree
}

.getNodeShape <- function(node){
  if(node$isRoot || node$nodeType == "Platform"){
    shape = "oval"
  } else {
    shape = "box"
  }
}


#' @export
#' @rdname precisionPathwaysEvaluations
strataPlot <- function (precisionPathways, ...) {
   UseMethod("strataPlot", precisionPathways)
 }

#' @param classColours A named vector of colours with names being \code{"class1"},\code{"class2"}, and \code{"accuracy"}.
#' a default colour scheme will automatically be chosen.
#' @rdname precisionPathwaysEvaluations
#' @export
strataPlot.PrecisionPathways <- function(precisionPathways, pathway, classColours = c(class1 = "#3F48CC", class2 = "#880015"), ...)
{
  if(missing(pathway)) stop("'pathway' is not specified. Please specify one of the ones shown by printing the trained object.")    
  pathwayUse <- precisionPathways[["pathways"]][[pathway]]
  assayIDs <- pathwayUse[["tiers"]][, 1]
  possibleClasses <- levels(precisionPathways$models[[1]]@actualOutcome)
  samplesTiers <- pathwayUse$individuals
  
  if("testSampleInfo" %in% names(precisionPathways))
  {
    knownClasses <- precisionPathways[["testSampleInfo"]][, "class"]
    allNames <- precisionPathways[["testSampleInfo"]][, "sampleID"]
    samplesTiers$trueClass <- knownClasses[match(samplesTiers[, "Sample ID"], allNames)]
  } else {
    samplesTiers$trueClass <- actualOutcome(precisionPathways$models[[1]])[match(samplesTiers[, "Sample ID"], sampleNames(precisionPathways$models[[1]]))]     
  }
  
  samplesTiers <- dplyr::arrange(as.data.frame(samplesTiers), Tier, trueClass, Accuracy)
  samplesTiers$ID = 1:nrow(samplesTiers)
  samplesTiers$Tier <- droplevels(samplesTiers$Tier)
  # The true classes are drawn as a strip of tiles above the tiers.
  samplesTiers$trueClassRow <- length(levels(samplesTiers[, "Tier"])) + 0.8
  trueClassColours <- setNames(unname(classColours[c("class1", "class2")]), levels(samplesTiers[, "Predicted"])[1:2])
  strataPlot <- ggplot2::ggplot(mapping = ggplot2::aes(x = ID, y = Tier), data = samplesTiers) +
    ggplot2::geom_tile(ggplot2::aes(fill = Accuracy)) +
    ggplot2::scale_fill_gradient2(low = "#86C57C", mid = "white", high = "#FFFF28", midpoint = 0.5,
                                  name = "Accuracy", guide = ggplot2::guide_colorbar(title.position = "top")) +
    ggplot2::labs(title = paste("Pathway:", pathway), x = "", y = "") +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(ggplot2::aes(y = trueClassRow, fill = trueClass), height = 0.6) +
    ggplot2::scale_fill_manual(values = trueClassColours, na.value = "grey70",
                               name = "True Class", guide = ggplot2::guide_legend(title.position = "top")) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          aspect.ratio = 1/4, 
          plot.title = ggplot2::element_text(face = "bold", size = 20),
          legend.title = ggplot2::element_text(face = "bold", size = 12),
          legend.text = ggplot2::element_text(size = 10),
          legend.position = "bottom",
          axis.text = ggplot2::element_text(size = 15)) +
    ggplot2::coord_cartesian(expand = FALSE) 
    
  strataPlot
}
