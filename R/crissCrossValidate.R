#' A function to perform pairwise cross validation
#'
#' This function has been designed to perform cross-validation and model prediction on datasets in a pairwise manner.
#'
#' @param measurements A \code{list} of either \code{\link{DataFrame}}, \code{\link{data.frame}} or \code{\link{matrix}} class measurements.
#' @param outcomes A \code{list} of vectors that respectively correspond to outcomes of the samples in \code{measurements} list. /
#' Factors should be coded such that the control class is the first level.
#' @param nFeatures The number of features to be used for modelling.
#' @param selectionMethod Default: \code{"auto"}. A character keyword of the feature algorithm to be used. If \code{"auto"}, t-test (two categories) /
#' F-test (three or more categories) ranking and top \code{nFeatures} optimisation is done. Otherwise, the ranking method is per-feature Cox proportional
#' hazards p-value.
#' @param selectionOptimisation A character of "Resubstitution", "Nested CV" or "none" specifying the approach used to optimise nFeatures.
#' @param trainType Default: \code{"modelTrain"}. A keyword specifying whether a fully trained model is used to make predictions on the test
#' set or if only the feature identifiers are chosen using the training data set and a number of training-predictions are made by cross-validation
#' in the test set.
#' @param classifier Default: \code{"auto"}. A character keyword of the modelling algorithm to be used. If \code{"auto"}, then a random forest is used
#' for a classification task or Cox proportional hazards model for a survival task.
#' @param nFolds A numeric specifying the number of folds to use for cross-validation.
#' @param nRepeats A numeric specifying the number of repeats or permutations to use for cross-validation.
#' @param nCores A numeric specifying the number of cores used if the user wants to use parallelisation. 
#' @param performanceType Default: \code{"auto"}. If \code{"auto"}, then balanced accuracy for classification or C-index for survival. Otherwise, any one of the
#' options described in \code{\link{calcPerformance}} may otherwise be specified.
#' @param doRandomFeatures Default: \code{FALSE}. Whether to perform random feature selection to establish a baseline performance. Either \code{FALSE} or \code{TRUE}
#' are permitted values.
#' @param runTOP Default: \code{FALSE}. If \code{TRUE}, perform the Transferable Omics Prediction (TOP) procedure in a leave-one-dataset-out manner. 
#' @param verbose Default: 0. A number between 0 and 3 for the amount of progress messages to give.  A higher number will produce more messages.
#' 
#' @return A list with elements \code{"real"} for the matrix of pairwise performance metrics using real
#' feature selection, \code{"random"} if \code{doRandomFeatures} is \code{TRUE} for metrics of random selection, 
#' \code{"top"} if \code{runTOP} is \code{TRUE}, and \code{"params"} for a list of parameters used.
#'
#' @author Harry Robertson
#' @export
crissCrossValidate <- function(measurements, outcomes,
                               nFeatures = 20,
                               selectionMethod = "auto",
                               selectionOptimisation = "Resubstitution",
                               trainType = c("modelTrain", "modelTest"),
                               performanceType = "auto",
                               doRandomFeatures = FALSE,
                               runTOP = FALSE,
                               classifier = "auto",
                               nFolds = 5,
                               nRepeats = 20,
                               nCores = 1,
                               verbose = 0)
{
    trainType <- match.arg(trainType)
    
    if(!is.list(measurements)) stop("'measurements' is not of type list but is of type", class(measurements))
    if(is.null(names(measurements))) stop("Each element of 'measurements' must be named by the name of the data set.")
    if(!is.list(outcomes)) stop("'outcomes' is not of type list but is of type", class(outcomes))
    
    isCategorical <- is.character(outcomes[[1]]) && 
        (length(outcomes[[1]]) == 1 || length(outcomes[[1]]) == nrow(measurements[[1]])) ||
        is.factor(outcomes[[1]])
    
    # If user didn't specify performanceType, choose "Balanced Accuracy" for classification else "C-index" for survival
    if(performanceType == "auto") {
        if(isCategorical) performanceType <- "Balanced Accuracy" else performanceType <- "C-index"
    }
    # If user left selectionMethod as "auto", pick t-test for categorical or CoxPH for survival
    if(length(selectionMethod) == 1 && selectionMethod == "auto") {
        if(isCategorical) selectionMethod <- "t-test" else selectionMethod <- "CoxPH"
    }
    # If user left classifier as "auto", pick randomForest for categorical or CoxPH for survival
    if(length(classifier) == 1 && classifier == "auto") {
        if(isCategorical) classifier <- "randomForest" else classifier <- "CoxPH"
    }
    
    # Keep a copy of the original measurements/outcomes if runTOP is requested
    if (runTOP) {
        top_measurements <- measurements
        top_outcomes <- outcomes
    }
    
    dataCleaned <- mapply(function(measurementsOne, outcomesOne) {
        prepareData(measurementsOne, outcomesOne)
    }, measurements, outcomes, SIMPLIFY = FALSE)
    
    measurements <- lapply(dataCleaned, "[[", 1)
    outcomes <- lapply(dataCleaned, "[[", 2)
    
    if (trainType == "modelTrain") {
        if (verbose > 0) {
            message("Using built training models on all test data sets.")
        }
        # Build a model for each dataset
        trainedModels <- mapply(function(measurementsOne, outcomesOne) {
            train(measurementsOne, outcomesOne,
                  nFeatures            = nFeatures,
                  selectionMethod      = selectionMethod,
                  selectionOptimisation = selectionOptimisation,
                  classifier           = classifier,
                  multiViewMethod      = "none",
                  verbose              = verbose)
        }, measurements, outcomes, SIMPLIFY = FALSE)
        
        # Predict on each dataset
        performanceAllPairs <- lapply(trainedModels, function(trainedModel) {
            mapply(function(testData, testOutcomes) {
                predictions <- predict(trainedModel, testData, verbose = verbose)
                
                if (performanceType == "AUC") {
                    # Must have columns named after each factor level for multi-class AUC
                    if (is.factor(testOutcomes)) {
                        neededLevels <- levels(testOutcomes)
                        missingCols  <- setdiff(neededLevels, colnames(predictions))
                        if (length(missingCols) > 0) {
                            stop("Cannot compute AUC because the predicted probabilities\n",
                                 "do not have columns for these outcome level(s): ",
                                 paste(missingCols, collapse = ", "))
                        }
                        # Reorder columns to match factor level order
                        predictions <- predictions[, neededLevels, drop = FALSE]
                    }
                } else {
                    keepCols <- na.omit(match(c("class", "risk"), colnames(predictions)))
                    if (length(keepCols) > 0) {
                        predictions <- predictions[, keepCols]
                    }
                }
                
                # Evaluate performance
                calcExternalPerformance(testOutcomes, predictions, performanceType)
            }, measurements, outcomes)
        })
        
        realPerformance <- matrix(
            unlist(performanceAllPairs),
            ncol     = length(measurements),
            byrow    = TRUE,
            dimnames = list(paste("Select and Train", names(measurements)),
                            paste("Predict",         names(measurements)))
        )
        realPerformance <- round(realPerformance, 2)
        
    } else {
        # trainType == "modelTest"
        if (verbose > 0) {
            message("Using features chosen in training to do cross-validation in the test data sets.")
        }
        # crossValidate on each dataset to *select features* or do nested CV
        trainedModels <- mapply(function(measurementsOne, outcomesOne) {
            crossValidate(measurementsOne, outcomesOne,
                          nFeatures             = nFeatures,
                          selectionMethod       = selectionMethod,
                          selectionOptimisation = selectionOptimisation,
                          classifier            = classifier,
                          multiViewMethod       = "none",
                          nFolds                = nFolds,
                          nCores                = nCores,
                          nRepeats              = nRepeats,
                          verbose               = verbose)
        }, measurements, outcomes, SIMPLIFY = FALSE)
        
        # Build cross-validation parameters
        crossValParams <- generateCrossValParams(nRepeats, nFolds, nCores, selectionOptimisation)
        
        # Evaluate each "trainedModel" on all datasets
        performanceAllPairs <- lapply(trainedModels, function(trainedModel) {
            mapply(function(measurementsOne, outcomesOne) {
                classifierParams <- .classifierKeywordToParams(classifier, NULL)
                modellingParams <- ModellingParams(
                    selectParams = SelectParams("previousSelection",
                                                intermediate    = ".iteration",
                                                classifyResult  = trainedModel),
                    trainParams   = classifierParams$trainParams,
                    predictParams = classifierParams$predictParams
                )
                
                result <- runTests(measurementsOne, outcomesOne, crossValParams, modellingParams)
                avgPerf <- mean(performance(calcCVperformance(result, performanceType))[[performanceType]])
                avgPerf
            }, measurements, outcomes, SIMPLIFY = FALSE)
        })
        
        realPerformance <- matrix(
            unlist(performanceAllPairs),
            ncol     = length(measurements),
            byrow    = TRUE,
            dimnames = list(paste("Select", names(measurements)),
                            paste("Cross-validate", names(measurements)))
        )
        realPerformance <- round(realPerformance, 2)
    }
    
    result <- list(real = realPerformance) # Store the results
    
    if (doRandomFeatures) {
        message("Starting random feature selection procedure.")
        # For each dataset, pick random features (nFeatures)
        randomFeatures <- lapply(measurements, function(dataset) {
            sample(colnames(dataset), nFeatures)
        })
        
        performanceAllPairs <- lapply(randomFeatures, function(randomFeaturesSet) {
            mapply(function(testData, testOutcomes) {
                resultRand <- crossValidate(testData[, randomFeaturesSet, drop = FALSE],
                                            testOutcomes,
                                            nFeatures       = nFeatures,
                                            selectionMethod = "none",
                                            classifier      = classifier,
                                            multiViewMethod = "none",
                                            nFolds          = nFolds,
                                            nCores          = nCores,
                                            nRepeats        = nRepeats)
                
                mean(performance(calcCVperformance(resultRand, performanceType))[[performanceType]])
            }, measurements, outcomes)
        })
        
        randomPerformance <- matrix(
            unlist(performanceAllPairs),
            ncol     = length(measurements),
            byrow    = TRUE,
            dimnames = list(paste("Random Select", names(measurements)),
                            paste("Cross-validate", names(measurements)))
        )
        randomPerformance <- round(randomPerformance, 2)
        
        result$random <- randomPerformance
    }
    
    # TOP procedure
    if (runTOP) {
        if (verbose > 0) message("Starting TOP procedure.")
        if (!requireNamespace("TOP", quietly = TRUE)) {
            stop("The TOP package is required for runTOP.")
        }
        
        topPerformance <- numeric(length(top_measurements))
        names(topPerformance) <- names(top_measurements)
        
        allLev <- levels(factor(unlist(top_outcomes)))
        if (length(allLev) < 2) {
            stop("runTOP requires at least 2 factor levels in outcomes.")
        }
        
        for (i in seq_along(top_measurements)) {
            # Train on all but i, test on i
            topfeatures <- TOP::filterFeatures(top_measurements[-i], top_outcomes[-i],
                                               nFeatures = nFeatures,
                                               contrast = paste(allLev[2], "-", allLev[1])) 
            trainData <- lapply(top_measurements[-i], function(x) x[, topfeatures, drop = FALSE])
            trainOuts <- top_outcomes[-i]
            testData  <- top_measurements[[i]][, topfeatures, drop = FALSE]
            testOuts  <- top_outcomes[[i]]
            
            model <- TOP::TOP_model(trainData, trainOuts)
            predictedProb <- TOP::predict_TOP(model$models, testData) 
            
            if (performanceType == "AUC") {
                probDF <- data.frame(
                    Negative = 1 - predictedProb,
                    Positive = predictedProb
                )
                colnames(probDF) <- allLev[1:2] # Ensure order is correct. 
                topPerformance[i] <- calcExternalPerformance(testOuts, probDF, "AUC")
                
            } else {
                predictedClasses <- ifelse(predictedProb > 0.5, allLev[2], allLev[1])
                predictedClasses <- factor(predictedClasses, levels = allLev)
                topPerformance[i] <- calcExternalPerformance(testOuts, predictedClasses, performanceType)
            }
        }
        
        topMatrix <- matrix(topPerformance, nrow = 1,
                            dimnames = list("TOP", names(top_measurements)))
        result$top <- round(topMatrix, 2)
    }
    
    # Store the parameters
    result$params <- list(nFeatures             = nFeatures,
                          selectionMethod       = selectionMethod,
                          selectionOptimisation = selectionOptimisation,
                          classifier            = classifier,
                          nFolds                = nFolds,
                          nRepeats              = nRepeats,
                          nCores                = nCores,
                          trainType             = trainType,
                          performanceType       = performanceType,
                          doRandomFeatures      = doRandomFeatures,
                          runTOP                = runTOP)
    
    result
}

#' A function to plot the output of the crissCrossValidate function.
#'
#' This function generates a heatmap of the cross-validation results from
#' \code{\link{crissCrossValidate}}. By default, it hides the "resubstitution" diagonal
#' (where the training == test set) unless \code{showResubMetric = TRUE}.
#'
#' @param crissCrossResult The output of the \code{\link{crissCrossValidate}} function.
#' @param includeValues Logical. If \code{TRUE}, numeric values are printed on each tile.
#' @param showResubMetric Logical. If \code{FALSE}, the diagonal (resubstitution) cells
#'        are set to \code{NA} and appear grayed-out or blank. Defaults to \code{FALSE}.
#'
#' @import ggplot2
#' @import reshape2
#' @import ggpubr
#' @export
crissCrossPlot <- function(crissCrossResult,
                           includeValues    = FALSE,
                           showResubMetric  = FALSE)
{
    # We'll attach for convenience
    attach(crissCrossResult)
    on.exit(detach(crissCrossResult), add = TRUE)
    
    scalebar_title <- params$performanceType
    
    # Helper function: turn matrix into heatmap
    plotMatrix <- function(mat, xlab_text, ylab_text) {
        # Convert to matrix if needed
        mat <- as.matrix(mat)
        
        if (!showResubMetric) {
            diag(mat) <- NA
        }
        
        melted_df <- reshape2::melt(mat, na.rm = FALSE, value.name = "value")
        
        gg <- ggplot(melted_df, aes(x = Var1, y = Var2, fill = value)) +
            geom_tile(color = "white") +
            scale_fill_gradient2(
                high     = "#e25563ff",
                mid      = "white",
                low      = "#094bacff",
                midpoint = 0.5,
                limit    = c(0,1),
                space    = "Lab",
                name     = scalebar_title,
                na.value = "grey70"
            ) +
            theme_bw() +
            xlab(xlab_text) +
            ylab(ylab_text) +
            theme(
                axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1),
                axis.text.y = element_text(vjust = 1, size = 8, hjust = 1)
            ) +
            coord_fixed()
        
        if (includeValues) {
            gg <- gg + geom_text(aes(label = round(value, 3)), color = "black", size = 3, na.rm = TRUE)
        }
        gg
    }
    
    if (params$trainType == "modelTrain") {
        
        mainMatrix <- real
        if ("top" %in% names(crissCrossResult)) {  # If runTOP is TRUE and we have result$top, append it as an extra row
            mainMatrix <- rbind(mainMatrix, crissCrossResult$top)
        }
        
        heatmapObj <- plotMatrix(mainMatrix, "Training Dataset", "Testing Dataset")
        
    } else if (params$trainType == "modelTest") {
        # Real is "Features Extracted" vs. "Cross-validate"
        mainMatrix <- real
        heatmapObj1 <- plotMatrix(mainMatrix, "Features Extracted", "Dataset Tested")
        
        if (params$doRandomFeatures == TRUE && exists("random")) {
            heatmapObj2 <- plotMatrix(random, "Random Features", "Dataset Tested")
            heatmapObj <- ggarrange(
                heatmapObj1,
                heatmapObj2,
                labels        = c("A - Feature Selection", "B - Random Features"),
                ncol          = 2,
                common.legend = TRUE,
                legend        = "right"
            )
        } else {
            heatmapObj <- heatmapObj1
        }
    }
    
    print(heatmapObj)
}