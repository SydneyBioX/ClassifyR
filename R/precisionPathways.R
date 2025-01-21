#' Precision Pathways for Sample Prediction Based on Prediction Confidence.
#' 
#' Precision pathways allows the evaluation of various permutations of multiomics or multiview data.
#' Samples are predicted by a particular assay if they were consistently predicted as a particular class
#' during cross-validation. Otherwise, they are passed onto subsequent assays/tiers for prediction. Balanced accuracy
#' is used to evaluate overall prediction performance and sample-specific accuracy for individual-level evaluation.
#'
#' @param measurements Either a \code{\link{MultiAssayExperiment}} or a list of the basic tabular objects containing the data.
#' @param class If a \code{\link{MultiAssayExperiment}}, a column name in \code{colData(measurements)} with the classes. If \code{measurements} is a \code{list} of tabular data, may also be
#' a vector of classes.
#' @param useFeatures Default: \code{NULL} (i.e. use all provided features). A named list of features to use. Otherwise, the input data is a single table and this can just be a vector of feature names.
#' For any assays not in the named list, all of their features are used. \code{"clinical"} is also a valid assay name and refers to the clinical data table.
#' This allows for the avoidance of variables such spike-in RNAs, sample IDs, sample acquisition dates, etc. which are not relevant for outcome prediction.
#' @param maxMissingProp Default: 0.0. A proportion less than 1 which is the maximum
#' tolerated proportion of missingness for a feature to be retained for modelling.
#' @param topNvariance Default: NULL. An integer number of most variable features per assay to subset to.
#' Assays with less features won't be reduced in size.
#' @param fixedAssays A character vector of assay names specifying any assays which must be at the
#' beginning of the pathway.
#' @param confidenceCutoff The minimum confidence of predictions for a sample to be predicted by a particular issue
#' . If a sample was predicted to belong to a particular class a proportion \eqn{p} times, then the confidence is \eqn{2 \times |p - 0.5|}.
#' @param minAssaySamples An integer specifying the minimum number of samples a tier may have. If a subsequent tier
#' would have less than this number of samples, the samples are incorporated into the current tier.
#' @param nFeatures Default: 20. The number of features to consider during feature selection, if feature selection is done.
#' @param selectionMethod A named character vector of feature selection methods to use for the assays, one for each. The names must correspond to names of \code{measurements}.
#' @param classifier A named character vector of modelling methods to use for the assays, one for each. The names must correspond to names of \code{measurements}.
#' @param nFolds A numeric specifying the number of folds to use for cross-validation.
#' @param nRepeats A numeric specifying the the number of repeats or permutations to use for cross-validation.
#' @param nCores A numeric specifying the number of cores used if the user wants to use parallelisation.
#' @param pathways A set of pathways created by \code{precisionPathwaysTrain} which is an object of class \code{PrecisionPathways} to be used for predicting on a new data set.
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
                   fixedAssays = "clinical", confidenceCutoff = 0.8, minAssaySamples = 10,
                   nFeatures = 20, selectionMethod = setNames(c("none", rep("t-test", length(measurements))), c("clinical", names(measurements))),
                   classifier = setNames(c("elasticNetGLM", rep("randomForest", length(measurements))), c("clinical", names(measurements))),
                   nFolds = 5, nRepeats = 20, nCores = 1)
          {
            if(is.list(measurements)) # Ensure plain list has clinical data.
            {
              # One of the tables must be named "clinical".
              if (!any(names(measurements) == "clinical"))
                stop("One of the tables must be named \"clinical\".")
            }
            if(is.null(useFeatures) && fixedAssays == "clinical")
            {
              warning("No 'useFeatures' named list element for clincal data is specified. Clinical data often has\n", 
    "lots of uninformative variables. Please consider specifying useful features.")        
              useFeatures <- list(clinical = colnames(MultiAssayExperiment::colData(measurements)))
            }
              
            prepArgs <- list(measurements, outcomeColumns = class, useFeatures = useFeatures,
                             maxMissingProp = maxMissingProp, topNvariance = topNvariance)
            measurementsAndClass <- do.call(prepareData, prepArgs)
              
            .precisionPathwaysTrain(measurementsAndClass[["measurements"]], measurementsAndClass[["outcome"]],
                                   useFeatures = useFeatures, fixedAssays = fixedAssays, confidenceCutoff = confidenceCutoff,
                                   minAssaySamples = minAssaySamples, nFeatures = nFeatures,
                                   selectionMethod = selectionMethod, classifier = classifier,
                                   nFolds = nFolds, nRepeats = nRepeats, nCores = nCores)
          })

# Internal method which carries out all of the processing, obtaining reformatted data from the
# MultiAssayExperiment and list (of basic rectangular tables) S4 methods.
.precisionPathwaysTrain <- function(measurements, class, fixedAssays = "clinical",
                   useFeatures = useFeatures, confidenceCutoff = 0.8, minAssaySamples = 10,
                   nFeatures = 20, selectionMethod = setNames(c(NULL, rep("t-test", length(measurements))), c("clinical", names(measurements))),
                   classifier = setNames(c("elasticNetGLM", rep("randomForest", length(measurements))), c("clinical", names(measurements))),
                   nFolds = 5, nRepeats = 20, nCores = 1)
          {
            # Step 1: Determine all valid permutations of assays, taking into account the
            # assays to be used and which assays, if any, must be included.
            assayIDs <- unique(S4Vectors::mcols(measurements)[["assay"]])
            assaysPermutations <- .permutations(assayIDs, fixed = data.frame(seq_along(fixedAssays), fixedAssays))
            permutationIDs <- apply(assaysPermutations, 2, function(permutation) paste(permutation, collapse = '-'))
            
            # Step 2: Build a classifier for each assay using all of the samples.
            modelsList <- crossValidate(measurements, class, nFeatures, selectionMethod,
                                        classifier = classifier, nFolds = nFolds,
                                        nRepeats = nRepeats, nCores = nCores)
            modelsList <- lapply(modelsList, calcCVperformance, "Sample Accuracy") # Add sample accuracy, which can be subset later.

            # Step 3: Loop over each pathway and each assay in order to determine which samples are used at that level
            # and which are passed onwards.
            precisionPathways <- lapply(as.data.frame(assaysPermutations), function(permutation)
            {
              assaysProcessed <- character()
              samplesUsed <- character()
              individualsTableAll <- S4Vectors::DataFrame()
              tierTableAll <- S4Vectors::DataFrame()
              breakEarly = FALSE
              for(assay in permutation)
              {
                # Step 3a: Identify all samples which are consistently predicted.
                modelIndex <- match(assay, assayIDs)
                allPredictions <- predictions(modelsList[[modelIndex]])
                allSampleIDs <- sampleNames(modelsList[[modelIndex]])
                predictionsSamplesCounts <- table(allPredictions[, "sample"], allPredictions[, "class"])
                confidences <- 2 * abs(predictionsSamplesCounts[, 1] / rowSums(predictionsSamplesCounts) - 0.5)
                sampleIDsUse <- names(confidences)[confidences > confidenceCutoff]
                sampleIDsUse <- setdiff(sampleIDsUse, samplesUsed)
                
                # Check if too few samples left for next round. Include them in this round, if so.
                remainingIDs <- setdiff(allSampleIDs, c(samplesUsed, sampleIDsUse))
                if(length(remainingIDs) < minAssaySamples)
                {
                  sampleIDsUse <- c(sampleIDsUse, remainingIDs)
                  breakEarly = TRUE
                }
                
                predictionsSamplesCounts <- predictionsSamplesCounts[sampleIDsUse, ]
                
                # Step 3b: Individuals predictions and sample-wise accuracy, tier-wise error.
                maxVotes <- apply(predictionsSamplesCounts, 1, function(sample) which.max(sample))
                predictedClasses <- factor(colnames(predictionsSamplesCounts)[maxVotes],
                                           levels = colnames(predictionsSamplesCounts))    
                individualsTable <- S4Vectors::DataFrame(Tier = assay,
                                                         `Sample ID` = sampleIDsUse,
                                                         `Predicted` = predictedClasses,
                                                         `Accuracy` = performance(modelsList[[modelIndex]])[["Sample Accuracy"]][sampleIDsUse],
                                                         check.names = FALSE)
                knownClasses <- actualOutcome(modelsList[[modelIndex]])[match(sampleIDsUse, allSampleIDs)]
                balancedAccuracy <- calcExternalPerformance(knownClasses, predictedClasses)
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
              list(pathway = pathwayString,
                  individuals = individualsTableAll, tiers = tierTableAll)
            })
            names(precisionPathways) <- sapply(precisionPathways, "[[", "pathway")
            precisionPathways <- precisionPathways[unique(names(precisionPathways))] # In case early termination causes duplicates.
            result <- list(models = modelsList, assaysPermutations = assaysPermutations,
                           parameters = list(confidenceCutoff = confidenceCutoff, minAssaySamples = minAssaySamples),
                           useFeatures = useFeatures, pathways = precisionPathways)
            class(result) <- "PrecisionPathways"
            
            result
}

# A nice print method to avoid flooding the screen with lots of tables
# when result is shown in console.
#' @export
print.PrecisionPathways <- function(x, ...)
{
  cat("An object of class 'PrecisionPathways'.\n")
  cat("Pathways:\n")
  cat(paste(names(x[["pathways"]]), collapse = '\n'))
}

#' @usage NULL
setGeneric("precisionPathwaysPredict", function(pathways, measurements, class, ...)
    standardGeneric("precisionPathwaysPredict"))

#' @rdname precisionPathways
#' @export
setMethod("precisionPathwaysPredict", c("PrecisionPathways", "MultiAssayExperimentOrList"), 
          function(pathways, measurements, class)
          {
            if(is.list(measurements)) # Ensure plain list has clinical data.
            {
              # One of the tables must be named "clinical".
              if (!any(names(measurements) == "clinical"))
                stop("One of the tables must be named \"clinical\".")
            }

            prepArgs <- list(measurements, outcomeColumns = class, useFeatures = pathways[["useFeatures"]])
            measurementsAndClass <- do.call(prepareData, prepArgs)
              
            .precisionPathwaysPredict(pathways, measurementsAndClass[["measurements"]], measurementsAndClass[["outcome"]])
          })

.precisionPathwaysPredict <- function(pathways, measurements, class)
{

  # Step 1: Extract all of previously fitted models and permutations.
  modelsList <- pathways[["models"]]
  assayIDs <- lapply(pathways[["models"]], function(model) model@characteristics[model@characteristics[, 1] == "Assay Name", 2])
  assaysPermutations <- pathways[["assaysPermutations"]]
  confidenceCutoff <- pathways[["parameters"]][["confidenceCutoff"]]
  minAssaySamples <- pathways[["parameters"]][["minAssaySamples"]]
  
  # Step 2: Loop over each pathway and each assay in order to determine which samples are used at that level
  # and which are passed onwards.
  precisionPathways <- lapply(as.data.frame(assaysPermutations), function(permutation)
  {
    assaysProcessed <- character()
    samplesUsed <- character()
    individualsTableAll <- S4Vectors::DataFrame()
    tierTableAll <- S4Vectors::DataFrame()
    breakEarly = FALSE
    for(assay in permutation)
    {
      # Step 2a: Identify all samples which are consistently predicted.
      modelIndex <- match(assay, assayIDs)
      allPredictions <- predictions(modelsList[[modelIndex]])
      allSampleIDs <- sampleNames(modelsList[[modelIndex]])
      predictionsSamplesCounts <- table(allPredictions[, "sample"], allPredictions[, "class"])
      confidences <- 2 * abs(predictionsSamplesCounts[, 1] / rowSums(predictionsSamplesCounts) - 0.5)
      sampleIDsUse <- names(confidences)[confidences > confidenceCutoff]
      sampleIDsUse <- setdiff(sampleIDsUse, samplesUsed)
                
      # Check if too few samples left for next round. Include them in this round, if so.
      remainingIDs <- setdiff(allSampleIDs, c(samplesUsed, sampleIDsUse))
      if(length(remainingIDs) < minAssaySamples)
      {
        sampleIDsUse <- c(sampleIDsUse, remainingIDs)
        breakEarly = TRUE
      }
                
      predictionsSamplesCounts <- predictionsSamplesCounts[sampleIDsUse, ]
                
      # Step 2b: Individuals predictions and sample-wise accuracy, tier-wise error.
      maxVotes <- apply(predictionsSamplesCounts, 1, function(sample) which.max(sample))
      predictedClasses <- factor(colnames(predictionsSamplesCounts)[maxVotes],
                                 levels = colnames(predictionsSamplesCounts))    
      individualsTable <- S4Vectors::DataFrame(Tier = assay,
                                               `Sample ID` = sampleIDsUse,
                                               `Predicted` = predictedClasses,
                                               `Accuracy` = performance(modelsList[[modelIndex]])[["Sample Accuracy"]][sampleIDsUse],
                                                check.names = FALSE)
      knownClasses <- actualOutcome(modelsList[[modelIndex]])[match(sampleIDsUse, allSampleIDs)]
      balancedAccuracy <- calcExternalPerformance(knownClasses, predictedClasses)
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
    list(pathway = pathwayString,
         individuals = individualsTableAll, tiers = tierTableAll)
  })
  names(precisionPathways) <- sapply(precisionPathways, "[[", "pathway")
  result <- list(models = modelsList, assaysPermutations = assaysPermutations,
                 parameters = list(confidenceCutoff = confidenceCutoff, minAssaySamples = minAssaySamples),
                 pathways = precisionPathways)
  class(result) <- "PrecisionPathways"
            
  result
}

# Calculate accuracy and costs of each pathway.

#' Various Functions for Evaluating Precision Pathways
#' 
#' These functions tabulate or plot various aspects of precision pathways, such as accuracies and costs.
#' 
#' @param precisionPathways A pathway of class \code{PrecisionPathways}.
#' @param costs A named vector of assays with the cost of each one.
#' @rdname precisionPathwaysEvaluations
#' @export
calcCostsAndPerformance <- function(precisionPathways, costs = NULL)
{
  if(is.null(costs))
    stop("'costs' of each assay must be specified.")
  pathwayIDs <- names(precisionPathways[["pathways"]])
  accuraciesCosts <- do.call(rbind, lapply(precisionPathways[["pathways"]], function(pathway)
  {
    predictions <- pathway[["individuals"]][, "Predicted"]
    knownClasses <- actualOutcome(precisionPathways$models[[1]])
    allNames <- sampleNames(precisionPathways$models[[1]])
    knownClasses <- knownClasses[match(pathway[["individuals"]][, "Sample ID"], allNames)]
    balancedAccuracy <- calcExternalPerformance(knownClasses, predictions)
    costTotal <- sum(costs[na.omit(match(pathway[["individuals"]][, "Tier"], names(costs)))])
    if(is.na(costTotal)) # Only clinical used.
      costTotal <- 0
    
    data.frame(accuracy = round(balancedAccuracy, 2), cost = costTotal)
  }))

  precisionPathways$performance <- accuraciesCosts
  precisionPathways
}

# Print a summary table, including accuracy and costs.

#' @param object A set of pathways of class \code{PrecisionPathways}.
#' @param weights A numeric vector of length two specifying how to weight the predictive accuracy
#' and the cost during ranking. Must sum to 1.
#' @param ... Not used but just following the S3 requirement of the generic template.
#' @rdname precisionPathwaysEvaluations
#' @export
summary.PrecisionPathways <- function(object, weights = c(accuracy = 0.5, cost = 0.5), ...)
{
  if(!"performance" %in% names(object))
    stop("Nothing to summarise. Please run function 'calcCostsAndPerformance' first.")
  hasCost <- "cost" %in% names(object[["performance"]])   
  summaryTable <- data.frame(Pathway = rownames(object[["performance"]]),
                             `Balanced Accuracy` = object[["performance"]][, "accuracy"],
                              check.names = FALSE)
  if(hasCost) summaryTable[, "Total Cost"] <- object[["performance"]][, "cost"]
  
  if(hasCost)
    rankingScores <- list(rank(object[["performance"]][, "accuracy"]), rank(-object[["performance"]][, "cost"]))
  else
    rankingScores <- list(rank(object[["performance"]][, "accuracy"]))  
  finalScores <- mapply(function(scores, weight) scores * weight, rankingScores, as.list(weights))
  if(!is(finalScores, "tabular")) finalScores <- matrix(finalScores, nrow = 1)
  finalScores <- rowSums(finalScores)
  summaryTable <- cbind(summaryTable, Score = finalScores)
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
  ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)))    
  if(is.null(pathwayColours)) pathwayColours <- scales::hue_pal()(length(precisionPathways[["pathways"]]))
  performance <- precisionPathways[["performance"]]
  performance <- cbind(Sequence = rownames(performance), performance)
  ggplot2::ggplot(performance, aes(x = accuracy, y = cost, colour = Sequence, size = 4)) + ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = pathwayColours) + ggplot2::labs(x = "Balanced Accuracy", y = "Total Cost") + ggplot2::guides(size = "none")
}

#' @export
#' @rdname precisionPathwaysEvaluations
flowchart <- function (precisionPathways, ...) {
   UseMethod("flowchart", precisionPathways)
 }

#' @param precisionPathways A pathway of class \code{PrecisionPathways}.
#' @param pathway A character vector of length 1 specifying which pathway to plot, e.g. "clinical-mRNA".
#' @param nodeColours A named vector of colours with names being \code{"assay"}, \code{"class1"},\code{"class2"}.
#' a default colour scheme will automatically be chosen.
#' @rdname precisionPathwaysEvaluations
#' @export
flowchart.PrecisionPathways <- function(precisionPathways, pathway, nodeColours = c(assay = "#86C57C", class1 = "#ACCEE0", class2 = "#F47F72"), ...)
{
  if(!requireNamespace("data.tree", quietly = TRUE))
    stop("The package 'data.tree' could not be found. Please install it.")

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
          
  data.tree::SetGraphStyle(pathwayTree, rankdir = "LR")
  data.tree::SetEdgeStyle(pathwayTree, fontname = 'helvetica', label = .getEdgeLabel)
  data.tree::SetNodeStyle(pathwayTree, style = "filled", shape = .getNodeShape, fontcolor = "black", fillcolor = .getFillColour, fontname = 'helvetica')
  plot(pathwayTree)
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
strataPlot.PrecisionPathways <- function(precisionPathways, pathway, classColours = c(class1 = "#4DAF4A", class2 = "#984EA3"), ...)
{
  if(missing(pathway)) stop("'pathway' is not specified. Please specify one of the ones shown by printing the trained object.")    
  pathwayUse <- precisionPathways[["pathways"]][[pathway]]
  assayIDs <- pathwayUse[["tiers"]][, 1]
  possibleClasses <- levels(precisionPathways$models[[1]]@actualOutcome)
  samplesTiers <- pathwayUse$individuals
  samplesTiers$trueClass <- actualOutcome(precisionPathways$models[[1]])[match(samplesTiers[, "Sample ID"], sampleNames(precisionPathways$models[[1]]))]
  samplesTiers <- dplyr::arrange(as.data.frame(samplesTiers), Tier, trueClass, Accuracy)
  samplesTiers$ID = 1:nrow(samplesTiers)
  samplesTiers$colour = ifelse(samplesTiers$trueClass == levels(samplesTiers[, "Predicted"])[1], classColours["class1"], classColours["class2"])
  samplesTiers$Tier <- droplevels(samplesTiers$Tier)
  strataPlot <- ggplot2::ggplot(mapping = ggplot2::aes(x = ID, y = Tier), data = samplesTiers) +
                ggplot2::geom_tile(ggplot2::aes(fill = trueClass)) +
    ggplot2::scale_fill_manual(values = unname(classColours))  +
    ggplot2::labs(title = paste("Pathway:", pathway), fill = "True Class", x = "", y = "") +
    ggplot2::guides(fill = ggplot2::guide_legend(title.position = "top")) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(ggplot2::aes(fill = Accuracy)) +
    ggplot2::scale_fill_gradient(low = "#377EB8", high = "#E41A1C") +
    ggplot2::labs(fill = "Accuracy") +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top")) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          aspect.ratio = 1/4, 
          plot.title = ggplot2::element_text(face = "bold", size = 20),
          legend.title = ggplot2::element_text(face = "bold", size = 12),
          legend.text = ggplot2::element_text(size = 10),
          legend.position = "bottom",
          axis.text = ggplot2::element_text(size = 15)) +
    ggplot2::annotate("tile",
               x = samplesTiers$ID,
               y = length(levels(samplesTiers[, "Tier"])) + 0.8,
               height = 0.6,
               fill = samplesTiers$colour)  +
    ggplot2::coord_cartesian(expand = FALSE) 
    
  strataPlot
}