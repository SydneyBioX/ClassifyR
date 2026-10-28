#' Add Performance Calculations to a ClassifyResult Object or Calculate for a
#' Pair of Factor Vectors
#' 
#' If \code{calcExternalPerformance} is used, such as when having a vector of
#' known classes and a vector of predicted classes determined outside of the
#' ClassifyR package, a single metric value is calculated. If
#' \code{calcCVperformance} is used, annotates the results of calling
#' \code{\link{crossValidate}}, \code{\link{runTests}} or \code{\link{runTest}} with one of the user-specified performance measures.
#' 
#' All metrics except Matthews Correlation Coefficient are suitable for
#' evaluating classification scenarios with more than two classes and are
#' reimplementations of those available from Intel DAAL.
#' 
#' \code{\link{crossValidate}}, \code{\link{runTests}} or \code{\link{runTest}} was run in resampling mode, one performance
#' measure is produced for every resampling. Otherwise, if the leave-k-out mode was used,
#' then the predictions are concatenated, and one performance measure is
#' calculated for all classifications.
#' 
#' \code{"Balanced Error"} calculates the balanced error rate and is better
#' suited to class-imbalanced data sets than the ordinary error rate specified
#' by \code{"Error"}. \code{"Sample Error"} calculates the error rate of each
#' sample individually. This may help to identify which samples are
#' contributing the most to the overall error rate and check them for
#' confounding factors. Precision, recall and F1 score have micro and macro
#' summary versions. The macro versions are preferable because the metric will
#' not have a good score if there is substantial class imbalance and the
#' classifier predicts all samples as belonging to the majority class.
#' 
#' @aliases calcPerformance calcExternalPerformance calcCVperformance
#' calcExternalPerformance,factor,factor-method calcExternalPerformance,Surv,numeric-method
#' calcCVperformance,ClassifyResult-method
#' @param result An object of class \code{\link{ClassifyResult}}.
#' @param resultsList A list of modelling results. Each element must be of type \code{\link{ClassifyResult}}.
#' @param performanceTypes Default: \code{"auto"} A character vector. If \code{"auto"}, Balanced Accuracy will be used
#' for a classification task and C-index for a time-to-event task. If using \code{easyHard}, the default is
#' \code{"Sample Accuracy"} for a classification task and \code{"Sample C-index"} for a time-to-event task.
#' Must be one of the following options:
#' \itemize{
#' \item{\code{"Error"}: Ordinary error rate.}
#' \item{\code{"Accuracy"}: Ordinary accuracy.}
#' \item{\code{"Balanced Error"}: Balanced error rate.}
#' \item{\code{"Balanced Accuracy"}: Balanced accuracy.}
#' \item{\code{"Sample Error"}: Error rate for each sample in the data set.}
#' \item{\code{"Sample Accuracy"}: Accuracy for each sample in the data set.}
#' \item{\code{"Micro Precision"}: Sum of the number of correct predictions in each 
#'         class, divided by the sum of number of samples predicted as
#'         belonging to each class.}
#' \item{\code{"Micro Recall"}: Sum of the number of correct predictions in
#'         each class, divided by the sum of number of samples in each class.}
#' \item{\code{"Micro F1"}: F1 score obtained by calculating the
#' harmonic mean of micro precision and micro recall.}
#' \item{\code{"Macro Precision"}: Sum of the ratios of the number of correct predictions in each
#' class to the number of samples predicted to be in each class, divided by the number of classes.}
#' \item{\code{"Macro Recall"}: Sum of the ratios of the number of correct predictions
#' in each class to the number of samples in each class, divided by the number of classes.}
#' \item{\code{"Macro F1"}: F1 score obtained by calculating the harmonic mean of macro precision
#' and macro recall.}
#' \item{\code{"Matthews Correlation Coefficient"}: Matthews Correlation Coefficient (MCC). A score
#' between -1 and 1 indicating how concordant the predicted classes are to the actual classes. Only defined if
#' there are two classes.}
#' \item{\code{"AUC"}: Area Under the Curve. An area ranging from 0 to 1, under the ROC.}
#' \item{\code{"C-index"}: For survival data, the concordance index, for models which produce risk scores. Ranges from 0 to 1.}
#' \item{\code{"Sample C-index"}: Per-individual C-index.}
#' }
#' 
#' @param actualOutcome A factor vector or survival information specifying each sample's known outcome.
#' @param predictedOutcome A factor vector or survival information of the same length as \code{actualOutcome} specifying each sample's predicted outcome.
#' @param grouping Default: \code{"permutation"}. If the cross-validation was k-fold, then this determines whether the metric will be calculated for samples
#' grouped by permutation or by fold, if the value is \code{"fold"}. For small sample sizes, \code{"permutation"} would suit. But, for large sample sizes,
#' \code{"fold"} would be preferable, as class membership probabilities or risk scores are not directly comparable between folds. This setting makes
#' no difference to error or accuracy metrics, apart from their variability.
#' 
#' @return If \code{calcCVperformance} was run, an updated
#' \code{\linkS4class{ClassifyResult}} object, with new metric values in the
#' \code{performance} slot. If \code{calcExternalPerformance} was run, the
#' performance metric value itself.
#' 
#' @author Dario Strbenac
#' @examples
#' 
#'   predictTable <- DataFrame(sample = paste("A", 1:10, sep = ''),
#'                             class = factor(sample(LETTERS[1:2], 50, replace = TRUE)))
#'   actual <- factor(sample(LETTERS[1:2], 10, replace = TRUE))                             
#'   result <- ClassifyResult(DataFrame(characteristic = "Data Set", value = "Example"),
#'                            paste("A", 1:10, sep = ''), paste("Gene", 1:50), list(paste("Gene", 1:50), paste("Gene", 1:50)), list(paste("Gene", 1:5), paste("Gene", 1:10)),
#'                            list(function(oracle){}), NULL, predictTable, actual)
#'   result <- calcCVperformance(result) 
#'   performance(result)
#' 
#' @include classes.R
#' @rdname calcPerformance
#' @usage NULL
#' @export
setGeneric("calcExternalPerformance", function(actualOutcome, predictedOutcome, ...)
standardGeneric("calcExternalPerformance"))

#' @rdname calcPerformance
#' @exportMethod calcExternalPerformance
setMethod("calcExternalPerformance", c("factor", "factor"),
          function(actualOutcome, predictedOutcome, # Both are classes.
                   performanceTypes = "auto")
{
  if(length(performanceTypes) == 1 && performanceTypes == "auto") performanceTypes <- "Balanced Accuracy"
              
  if(length(levels(actualOutcome)) > 2 && "Matthews Correlation Coefficient" %in% performanceTypes)
    stop("Error: Matthews Correlation Coefficient specified but data set has more than 2 classes.")
  unknownClasses <- setdiff(as.character(predictedOutcome), c(levels(actualOutcome), NA))
  if(length(unknownClasses) > 0)
    stop("Predicted classes ", paste(unknownClasses, collapse = ", "), " are not classes of 'actualOutcome'.")
  predictedOutcome <- factor(as.character(predictedOutcome), levels = levels(actualOutcome)) # Match classes by name.
  
  sapply(performanceTypes, function(performanceType)
    .calcPerformance(list(actualOutcome), list(predictedOutcome), performanceType = performanceType)[["values"]]
  )
})

#' @rdname calcPerformance
#' @exportMethod calcExternalPerformance
setMethod("calcExternalPerformance", c("Surv", "numeric"),
          function(actualOutcome, predictedOutcome, performanceTypes = "auto")
          {
            if(length(performanceTypes) == 1 && performanceTypes == "auto") performanceTypes <- "C-index"
            
            sapply(performanceTypes, function(performanceType)
              .calcPerformance(actualOutcome, predictedOutcome, performanceType = performanceType)[["values"]]
            )
          })

#' @rdname calcPerformance
#' @exportMethod calcExternalPerformance
setMethod("calcExternalPerformance", c("factor", "tabular"), # table has class probabilities per sample.
          function(actualOutcome, predictedOutcome, performanceTypes = "auto")
          {
            if(length(performanceTypes) == 1 && performanceTypes == "auto") performanceTypes <- "AUC"
            
            sapply(performanceTypes, function(performanceType)
              .calcPerformance(actualOutcome, predictedOutcome, performanceType = performanceType)[["values"]]
            )
          })

calcCVperformance <- function(results, ...)
{
    lapply(results, calcCVperformance)
}

#' @rdname calcPerformance
#' @usage NULL
#' @export
setGeneric("calcCVperformance", function(result, ...)
    standardGeneric("calcCVperformance"))

#' @rdname calcPerformance
#' @exportMethod calcCVperformance
setMethod("calcCVperformance", "ClassifyResult",
          function(result, performanceTypes = "auto", grouping = c("permutation", "fold"))
{
  actualOutcome <- actualOutcome(result) # Extract the known outcome of each sample.
  actualOutcomeOrdered <- actualOutcome[match(result@predictions[, "sample"], sampleNames(result))]
  grouping <- match.arg(grouping)
  
  if(length(performanceTypes) == 1 && performanceTypes == "auto")
  {
      if(is.factor(actualOutcome))
      {
        performanceTypes <- "Balanced Accuracy"  
        if(all(levels(actualOutcome) %in% colnames(predictions(result))) && "fold" %in% colnames(result@predictions) && max(result@predictions[, "fold"]) < length(result@originalNames)) # Class scores present.
          performanceTypes <- c(performanceTypes, "AUC")
    } else performanceTypes <- "C-index"
  }
  
  # Allow calculation of multiple metrics at once.
  for(performanceType in performanceTypes)
  {
      ### Group by lowest level of test set.
      if(!performanceType %in% c("Sample Error", "Sample Accuracy"))
      {
        if("fold" %in% colnames(result@predictions)) # k-Fold or repeated k-Fold cross-validation.
        {
          if(grouping == "permutation")
            groupID <- result@predictions[, "permutation"]
          else
            groupID <- paste(result@predictions[, "permutation"], result@predictions[, "fold"], sep = ':')
        } else if("permutation" %in% colnames(result@predictions)) # Monte Carlo cross-validation.
        {
          groupID <- result@predictions[, "permutation"]      
        } 
        else { # Leave-k-out or independent train and test set, such as created by runTest function.
          groupID <- rep(1, nrow(result@predictions))
        }
      }
      
      ### Performance for survival data
      if(performanceType %in% c("C-index", "Sample C-index")) {
        samples <- factor(result@predictions[, "sample"], levels = sampleNames(result))
        performance <- .calcPerformance(actualOutcome = actualOutcomeOrdered,
                                        predictedOutcome = result@predictions[, "risk"], 
                                        samples = samples,
                                        performanceType = performanceType, 
                                        grouping = groupID)
        if(grepl(':', names(performance[["values"]])[1])) # Then average for each permutation.
          performance[["values"]] <- .averagePermutations(performance[["values"]])
        result@performance[[performance[["name"]]]] <- performance[["values"]]
      }
      
      if(performanceType == "AUC") {
        performance <- .calcPerformance(actualOutcomeOrdered,
                                        result@predictions[, levels(actualOutcome)],
                                        performanceType = performanceType, grouping = groupID)
        if(grepl(':', names(performance[["values"]])[1])) # Then average for each permutation.
          performance[["values"]] <- .averagePermutations(performance[["values"]])
        result@performance[[performance[["name"]]]] <- performance[["values"]]
      }
      
      ### Performance for data with classes
      if(performanceType %in% c("Balanced Accuracy", "Balanced Error", "Error", "Accuracy",
                                 "Micro Precision", "Micro Recall", "Micro F1", "Macro Precision",
                                 "Macro Recall", "Macro F1", "Matthews Correlation Coefficient",
                                 "Sample Error", "Sample Accuracy"))
      {
        if(length(levels(actualOutcome)) > 2 && performanceType == "Matthews Correlation Coefficient")
          stop("Error: Matthews Correlation Coefficient specified but data set has more than 2 classes.")
        
        classLevels <- levels(actualOutcome)
        samples <- factor(result@predictions[, "sample"], levels = sampleNames(result))
        predictedOutcome <- factor(result@predictions[, "class"], levels = classLevels)
        actualOutcome <- factor(actualOutcomeOrdered, levels = classLevels, ordered = TRUE)
        performance <- .calcPerformance(actualOutcome, predictedOutcome, samples, performanceType, groupID)
        result@performance[[performance[["name"]]]] <- performance[["values"]]
      }
  }
  result
})

# Averages values named "permutation:fold" within each permutation. Permutations are in numerical order.
.averagePermutations <- function(values)
{
  permutationIDs <- sapply(strsplit(names(values), ':'), '[', 1)
  uniqueIDs <- unique(permutationIDs)
  uniqueIDs <- uniqueIDs[order(suppressWarnings(as.numeric(uniqueIDs)), uniqueIDs)]
  sapply(split(values, factor(permutationIDs, levels = uniqueIDs)), mean)
}

#' @importFrom survival concordance
.calcPerformance <- function(actualOutcome, predictedOutcome, samples = NA, performanceType, grouping = NULL)
{
  # Make splitting by group safe.
  if(is(predictedOutcome, "DataFrame")) predictedOutcome <- as.data.frame(predictedOutcome, optional = TRUE)    
  if(performanceType %in% c("Sample Error", "Sample Accuracy"))
  {
    sampleMetricValues <- sapply(levels(samples), function(sampleID)
    {
      consider <- which(samples == sampleID)
      if(performanceType == "Sample Error")
        sum(predictedOutcome[consider] != as.character(actualOutcome[consider]))
      else
        sum(predictedOutcome[consider] == as.character(actualOutcome[consider]))
    })
    performanceValues <- as.numeric(sampleMetricValues / table(samples))
    names(performanceValues) <- levels(samples)
    return(list(name = performanceType, values = performanceValues))
  }

  if(!is.null(grouping))
  {
    actualOutcome <- split(actualOutcome, grouping)
    predictedOutcome <- split(predictedOutcome, grouping)
    if(length(samples) > 1 || !is.na(samples))
    {
      allSamples <- levels(samples)
      samples <- split(samples, grouping)
    }
  }
    
  if(performanceType == "Sample C-index")
  {
    # For each sample, count the concordant and discordant comparable pairs it belongs to, within each group.
    pairCounts <- do.call(rbind, mapply(function(iterationSurv, iterationPredictions, iterationSamples)
    {
      iterationSurv <- as.matrix(iterationSurv)
      times <- iterationSurv[, "time"]
      # Element [i, j] is TRUE if sample i has a shorter time than sample j and sample i is not censored.
      earlierEvent <- outer(times, times, '<') & iterationSurv[, "status"] == 1
      higherRisk <- outer(iterationPredictions, iterationPredictions, '>')
      lowerRisk <- outer(iterationPredictions, iterationPredictions, '<')
      # Concordant: the sample with the earlier event has the higher risk. Discordant: it has the lower risk.
      concordant <- earlierEvent & higherRisk
      discordant <- earlierEvent & lowerRisk
      data.frame(sample = iterationSamples,
                 concordant = rowSums(concordant) + colSums(concordant),
                 discordant = rowSums(discordant) + colSums(discordant))
    }, actualOutcome, predictedOutcome, samples, SIMPLIFY = FALSE))

    concordants <- rowsum(pairCounts[, "concordant"], pairCounts[, "sample"])[, 1]
    discordants <- rowsum(pairCounts[, "discordant"], pairCounts[, "sample"])[, 1]
    Cindex <- round(concordants / (concordants + discordants), 2)[allSamples]
    names(Cindex) <- allSamples
    Cindex[is.nan(Cindex)] <- NA # A censored individual with a time shorter than all events has no comparable pairs.
    return(list(name = performanceType, values = Cindex))
  }
  
  if(!is(actualOutcome, "list")) actualOutcome <- list(actualOutcome)
  if(!is(predictedOutcome, "list")) predictedOutcome <- list(predictedOutcome)
  

  if(performanceType %in% c("Accuracy", "Error")) {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      # Columns are predicted classes, rows are actual classes.
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      totalPredictions <- sum(confusionMatrix)
      correctPredictions <- sum(diag(confusionMatrix))
      diag(confusionMatrix) <- 0
      wrongPredictions <- sum(confusionMatrix)
      if(performanceType == "Accuracy")
        correctPredictions / totalPredictions
      else # It is "error".
        wrongPredictions / totalPredictions
    }, actualOutcome, predictedOutcome, SIMPLIFY = FALSE))
  } else if(performanceType %in% c("Balanced Accuracy", "Balanced Error")) {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      # Columns are predicted classes, rows are actual classes.
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      classSizes <- rowSums(confusionMatrix)
      classErrors <- classSizes - diag(confusionMatrix)
      if(performanceType == "Balanced Accuracy")
        mean(diag(confusionMatrix) / classSizes, na.rm = TRUE)
      else
        mean(classErrors / classSizes, na.rm = TRUE)
    }, actualOutcome, predictedOutcome, SIMPLIFY = FALSE))
  } else if(performanceType == "AUC") {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      classesTable <- do.call(rbind, lapply(levels(iterationClasses), function(class)
      {
        rates <- .ROCrates(iterationPredictions[, class], iterationClasses == class)
        data.frame(FPR = c(0, rates[["FPR"]]), TPR = c(0, rates[["TPR"]]), class = class)
      }))
      classesAUC <- .calcArea(classesTable, levels(actualOutcome[[1]]))
      mean(classesAUC[!duplicated(classesAUC[, c("class", "AUC")]), "AUC"]) # Average AUC in iteration.
    }, actualOutcome, predictedOutcome, SIMPLIFY = FALSE))
  } else if(performanceType %in% c("C-index")) {
    performanceValues <- unlist(mapply(function(x, y){
      y <- -y
      survival::concordance(x ~ y)$concordance
    }, actualOutcome, predictedOutcome, SIMPLIFY = FALSE))

    } else { # Metrics for which true positives, true negatives, false positives, false negatives must be calculated.
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      truePositives <- diag(confusionMatrix)
      falsePositives <- colSums(confusionMatrix) - truePositives
      falseNegatives <- rowSums(confusionMatrix) - truePositives
      trueNegatives <- sum(truePositives) - truePositives
      
      if(performanceType == "Average Accuracy")
      {
        return(sum((truePositives + trueNegatives) / (truePositives + trueNegatives + falsePositives + falseNegatives)) / nrow(confusionMatrix))
      }
      if(performanceType %in% c("Micro Precision", "Micro F1"))
      {
        microP <- sum(truePositives) / sum(truePositives + falsePositives)
        if(performanceType == "Micro Precision") return(microP)
      }
      if(performanceType %in% c("Micro Recall", "Micro F1"))
      {
        microR <- sum(truePositives) / sum(truePositives + falseNegatives)
        if(performanceType == "Micro Recall") return(microR)
      }
      if(performanceType == "Micro F1")
      {
        return(2 * microP * microR / (microP + microR))
      }
      if(performanceType %in% c("Macro Precision", "Macro F1"))
      {
        macroP <- sum(truePositives / (truePositives + falsePositives)) / nrow(confusionMatrix)
        if(performanceType == "Macro Precision") return(macroP)
      }
      if(performanceType %in% c("Macro Recall", "Macro F1"))
      {
        macroR <- sum(truePositives / (truePositives + falseNegatives)) / nrow(confusionMatrix)
        if(performanceType == "Macro Recall") return(macroR)
      }
      if(performanceType == "Macro F1")
      {
        return(2 * macroP * macroR / (macroP + macroR))
      }
      if(performanceType == "Matthews Correlation Coefficient")
      {
        return(unname((truePositives[2] * trueNegatives[2] - falsePositives[2] * falseNegatives[2]) / sqrt((truePositives[2] + falsePositives[2]) * (truePositives[2] + falseNegatives[2]) * (trueNegatives[2] + falsePositives[2]) * (trueNegatives[2] + falseNegatives[2]))))
      }
      
    }, actualOutcome, predictedOutcome, SIMPLIFY = FALSE))
  }

  list(name = performanceType, values = performanceValues)
}

#' @rdname calcPerformance
#' @param aggregate Default: \code{"median"}. Can also be \code{"mean"}. If there are multiple values, such as for repeated
#' cross-validation, then they are summarised to a single number using either mean or median.
#' @export performanceTable
#' @importFrom tidyr pivot_wider
performanceTable <- function(resultsList, performanceTypes = "auto", aggregate = c("median", "mean"))
{
  aggregate <- match.arg(aggregate)
  actualOutcome <- actualOutcome(resultsList[[1]]) # Establish outcome type.
  if(length(performanceTypes) == 1 && performanceTypes == "auto")
  {
      if(is.factor(actualOutcome)) performanceTypes <- c("AUC", "Balanced Accuracy") else performanceTypes <- "C-index"
  }

  names(performanceTypes) <- performanceTypes
  do.call(rbind, lapply(resultsList, function(result)
  {
    performances <- lapply(performanceTypes, function(performanceType)
    {
      if(!performanceType %in% names(performance(result)))
        performance <- performance(calcCVperformance(result, performanceType))[[performanceType]]
      else
        performance <- performance(result)[[performanceType]]
      if(aggregate == "median") median(performance) else mean(performance)
    })
    DataFrame(tidyr::pivot_wider(as.data.frame(result@characteristics), names_from = characteristic, values_from = value), performances, check.names = FALSE)  
  }))
}

#' @rdname calcPerformance
#' @usage NULL
#' @export
setGeneric("easyHard", function(measurements, result, ...)
    standardGeneric("easyHard"))

#' @rdname calcPerformance
#' @exportMethod easyHard
#' @importFrom broom tidy
#' @param measurements For \code{easyHard} only. Either a \code{\link{DataFrame}}, \code{\link{data.frame}}, \code{\link{matrix}}, \code{\link{MultiAssayExperiment}} 
#' or a list of the basic tabular objects containing the data.
#' @param assay For \code{easyHard} only. The assay to use to look for associations to the per-sample metric.
#' @param useFeatures For \code{easyHard} only. Default: \code{NULL} (i.e. use all provided features). A vector of features to consider of the assay specified.
#' This allows for the avoidance of variables such spike-in RNAs, sample IDs, sample acquisition dates, etc. which are not relevant for outcome prediction.
#' @param fitMode For \code{easyHard} only. Default:\code{"single"}. Either \code{"single"} or \code{"full"}. If \code{"single"},
#' an ordinary GLM model is fitted for each covariate separately. If \code{"full"}, elastic net is used to automatically tune the non-zero model coefficients.
#' @param performanceType For \code{easyHard} only. One of the valid values shown for \code{performanceType} parameter of \code{\link{calcCVperformance}}.
#' @return For \code{easyHard}, a \code{\link{DataFrame}} of logistic regression model summary.

setMethod("easyHard", "MultiAssayExperimentOrList",
          function(measurements, result, assay = "clinical", useFeatures = NULL, performanceType = "auto",
                   fitMode = c("single", "full"))
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
                            
  if(!assay %in% c(names(measurements), if(is(measurements, "MultiAssayExperiment")) "clinical"))
    stop("'assay' is not one of the names of 'measurements'.")
  fitMode  <- match.arg(fitMode)              
              
  if(is(measurements, "MultiAssayExperiment"))
  {
    if(assay == "clinical")
      assay <- colData(measurements)
    else assay <- t(as.matrix(MultiAssayExperiment::assays(measurements)[[assay]])) # Ensure that features are in columns.
  } else {assay <- measurements[[assay]]}
  if(is.matrix(assay)) assay <- S4Vectors::DataFrame(assay, check.names = FALSE)
  if(!is.null(useFeatures)) assay <- assay[, useFeatures, drop = FALSE]
  if(performanceType == "auto")
  {
      if("risk" %in% colnames(predictions(result)))
      {
          performanceType <- "Sample C-index"
      } else {performanceType <-"Sample Accuracy"}
  }
  if(!performanceType %in% names(performance(result)))
  {
    warning(paste(performanceType, "not found in result. Calculating it now."))
    result <- calcCVperformance(result, performanceType)
  }
  samplePerformance <- performance(result)[[performanceType]]
  samplePerformance <- samplePerformance[!is.na(samplePerformance)]
  if(!all(names(samplePerformance) %in% rownames(assay)))
    stop("Some samples of 'result' are not row names of the assay.")
  assay <- assay[names(samplePerformance), , drop = FALSE] # Match samples by name.
  assayOHE <- MatrixModels::model.Matrix(~ 0 + ., data = assay)
  
  if(fitMode == "single")
  {
    as(do.call(rbind, lapply(colnames(assay), function(featureID)
    {
      covariate <- assay[, featureID]
      fitted <- glm(samplePerformance ~ covariate, family = binomial, weights = rep(100, length(samplePerformance)))
      summaryDF <- broom::tidy(fitted)
      if(is.factor(covariate))
      {
        summaryDF[2:nrow(summaryDF), "term"] <- paste(featureID, levels(covariate)[2:length(levels(covariate))], sep = ": ")
        
      } else {summaryDF[, "term"] <- featureID}
      summaryDF[2:nrow(summaryDF), ]
    })), "DataFrame")
  } else { # Penalised regression.
    samplePerformanceM <- matrix(c(1 - samplePerformance, samplePerformance), ncol = 2)
    fitted <- glmnet::glmnet(assayOHE, samplePerformanceM, family = "binomial")
    lambdaConsider <- colSums(as.matrix(fitted[["beta"]])) != 0
    bestLambda <- fitted[["lambda"]][lambdaConsider][which.min(sapply(fitted[["lambda"]][lambdaConsider], function(lambda) # Largest Lambda with minimum balanced error rate.
    {
        predictions <- predict(fitted, assayOHE, s = lambda, type = "response")
        sum(abs(samplePerformanceM[, 2] - predictions))
    }))[1]]
    useVariables <- abs(fitted[["beta"]][, fitted[["lambda"]] == bestLambda]) > 0.00001
    useVariables <- colnames(assay)[unique(assayOHE@assign[useVariables])]
    dataForModel <- data.frame(assay, performance = samplePerformanceM[, 2])
    fitted <- glm(performance ~ . + 0, data = dataForModel, family = binomial(), weights = rep(100, nrow(dataForModel)))
    broom::tidy(fitted)
  }
})
# Calculates the true positive rate and false positive rate at each distinct score, from the highest
# score to the lowest. A sample is predicted positive if its score is at least the threshold.
.ROCrates <- function(scores, isPositive)
{
  totalPositives <- sum(isPositive)
  totalNegatives <- sum(!isPositive)
  scoresOrder <- order(scores, decreasing = TRUE, na.last = NA)
  scoresSorted <- scores[scoresOrder]
  truePositives <- cumsum(isPositive[scoresOrder])
  falsePositives <- cumsum(!isPositive[scoresOrder])
  # Keep the last sample of each group of tied scores, so that all tied samples are counted.
  thresholdEnds <- c(scoresSorted[-1] != scoresSorted[-length(scoresSorted)], TRUE)
  truePositives <- truePositives[thresholdEnds]
  falsePositives <- falsePositives[thresholdEnds]
  if(anyNA(scores)) # Comparisons to a missing score are unknown.
    truePositives <- falsePositives <- rep(NA_integer_, length(truePositives))
  list(FPR = falsePositives / totalNegatives, TPR = truePositives / totalPositives)
}
