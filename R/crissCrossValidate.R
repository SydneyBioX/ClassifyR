#' A function to perform pairwise cross-validation with optional logistic regression
#'
#' This function performs cross-validation and model prediction on datasets in a pairwise manner.
#' Additionally, if `runTOP` is `TRUE`, it performs logistic regression using a leave-one-dataset-out approach.
#'
#' @param measurements A \code{list} of either \code{\link{DataFrame}}, \code{\link{data.frame}}, or \code{\link{matrix}} class measurements.
#' @param outcomes A \code{list} of vectors that respectively correspond to outcomes of the samples in the \code{measurements} list.
#' @param nFeatures The number of features to be used for modelling.
#' @param selectionMethod Default: \code{"auto"}. A character keyword of the feature algorithm to be used. If \code{"auto"}, t-test (two categories) /
#' F-test (three or more categories) ranking and top \code{nFeatures} optimisation is done. Otherwise, the ranking method is per-feature Cox proportional
#' hazards p-value.
#' @param selectionOptimisation A character of "Resubstitution", "Nested CV", or "none" specifying the approach used to optimise \code{nFeatures}.
#' @param trainType Default: \code{"modelTrain"}. A keyword specifying whether a fully trained model is used to make predictions on the test
#' set or if only the feature identifiers are chosen using the training dataset and a number of training-predictions are made by cross-validation
#' in the test set.
#' @param classifier Default: \code{"auto"}. A character keyword of the modelling algorithm to be used. If \code{"auto"}, then a random forest is used
#' for a classification task or Cox proportional hazards model for a survival task.
#' @param nFolds A numeric specifying the number of folds to use for cross-validation.
#' @param nRepeats A numeric specifying the number of repeats or permutations to use for cross-validation.
#' @param nCores A numeric specifying the number of cores used if the user wants to use parallelisation.
#' @param performanceType Default: \code{"auto"}. If \code{"auto"}, then balanced accuracy for classification or C-index for survival. Otherwise, any one of the
#' options described in \code{\link{calcPerformance}} may otherwise be specified.
#' @param doRandomFeatures Default: \code{FALSE}. Whether to perform random feature selection to establish a baseline performance.
#' @param runTOP Default: \code{FALSE}. Whether to perform logistic regression using a leave-one-dataset-out approach.
#' @param verbose Default: 0. A number between 0 and 3 for the amount of progress messages to give. A higher number will produce more messages as
#' more lower-level functions print messages.
#' @return A list with elements \code{"real"} for the matrix of pairwise performance metrics using real
#' feature selection, \code{"random"} if \code{doRandomFeatures} is \code{TRUE} for metrics of random selection,
#' \code{"logisticRegression"} if \code{runTOP} is \code{TRUE} for logistic regression results, and
#' \code{"params"} for a list of parameters used during the execution of this function.
#' @author Harry Robertson
#'
#' @export

crissCrossValidate <- function(measurements, outcomes,
                               nFeatures = 20, selectionMethod = "auto",
                               selectionOptimisation = "Resubstitution",
                               trainType = c("modelTrain", "modelTest"),
                               performanceType = "auto",
                               doRandomFeatures = FALSE,
                               runTOP = FALSE,
                               classifier = "auto",
                               nFolds = 5, nRepeats = 20, nCores = 1, verbose = 0)
{
    trainType <- match.arg(trainType)
    if (!is.list(measurements)) stop("'measurements' must be a list but is of type ", class(measurements))
    if (is.null(names(measurements))) stop("Each element of 'measurements' must be named by the name of the dataset.")
    if (!is.list(outcomes)) stop("'outcomes' must be a list but is of type ", class(outcomes))
    
    # Determine if the outcome is categorical based on the first dataset
    isCategorical <- is.character(outcomes[[1]]) && (length(outcomes[[1]]) == 1 || length(outcomes[[1]]) == nrow(measurements[[1]])) || is.factor(outcomes[[1]])
    
    # Set default performanceType and selectionMethod based on the outcome type
    if (performanceType == "auto")
        if (isCategorical) performanceType <- "Balanced Accuracy" else performanceType <- "C-index"
    if (length(selectionMethod) == 1 && selectionMethod == "auto")
        if (isCategorical) selectionMethod <- "t-test" else selectionMethod <- "CoxPH"
    if (length(classifier) == 1 && classifier == "auto")
        if (isCategorical) classifier <- "randomForest" else classifier <- "CoxPH"
    
    # Clean the data for each dataset
    dataCleaned <- mapply(function(measurementsOne, outcomesOne)
    {
        prepareData(measurementsOne, outcomesOne)
    }, measurements, outcomes, SIMPLIFY = FALSE)
    measurements <- lapply(dataCleaned, "[[", 1)
    outcomes <- lapply(dataCleaned, "[[", 2)
    
    # Initialize the result list
    result <- list()
    
    # If trainType is "modelTrain", train models on each dataset and test on all datasets
    if (trainType == "modelTrain")
    {
        if (verbose > 0) message("Using built training models on all test datasets.")
        # Train a model on each dataset
        trainedModels <- mapply(function(measurementsOne, outcomesOne)
        {
            train(measurementsOne, outcomesOne,
                  nFeatures = nFeatures,
                  selectionMethod = selectionMethod, selectionOptimisation = selectionOptimisation,
                  classifier = classifier, multiViewMethod = "none", verbose = verbose)
        }, measurements, outcomes, SIMPLIFY = FALSE)
        
        # Evaluate each model on all datasets
        performanceAllPairs <- lapply(trainedModels, function(trainedModel)
        {
            mapply(function(testData, testOutcomes)
            {
                predictions <- predict(trainedModel, testData, verbose = verbose)
                if (is(predictions, "tabular")) predictions <- predictions[, na.omit(match(c("class", "risk"), colnames(predictions)))]
                calcExternalPerformance(testOutcomes, predictions, performanceType)
            }, measurements, outcomes)
        })
        
        # Compile the performance metrics into a matrix
        realPerformance <- matrix(unlist(performanceAllPairs), ncol = length(measurements), byrow = TRUE,
                                  dimnames = list(paste("Select and Train", names(measurements)), paste("Predict", names(measurements))))
        realPerformance <- round(realPerformance, 2)
        result$real <- realPerformance
    } else {
        # When trainType is "modelTest", select features using each training dataset and test them on all datasets using cross-validation
        if (verbose > 0) message("Using features chosen in training to do cross-validation in the test datasets.")
        # Perform cross-validation on each training dataset to select features
        trainedModels <- mapply(function(measurementsOne, outcomesOne)
        {
            crossValidate(measurementsOne, outcomesOne,
                          nFeatures = nFeatures,
                          selectionMethod = selectionMethod,
                          selectionOptimisation = selectionOptimisation,
                          classifier = classifier,
                          multiViewMethod = "none",
                          nFolds = nFolds,
                          nCores = nCores,
                          nRepeats = nRepeats, verbose = verbose)
        }, measurements, outcomes, SIMPLIFY = FALSE)
        
        # Generate cross-validation parameters
        crossValParams <- generateCrossValParams(nRepeats, nFolds, nCores, selectionOptimisation)
        
        # Evaluate performance on test datasets using the selected features
        performanceAllPairs <- lapply(trainedModels, function(trainedModel)
        {
            mapply(function(measurementsOne, outcomesOne)
            {
                classifierParams <- .classifierKeywordToParams(classifier)
                modellingParams <- ModellingParams(selectParams = SelectParams("previousSelection", intermediate = ".iteration", classifyResult = trainedModel),
                                                   trainParams = classifierParams$trainParams,
                                                   predictParams = classifierParams$predictParams)
                
                result <- runTests(measurementsOne, outcomesOne, crossValParams, modellingParams)
                mean(performance(calcCVperformance(result, performanceType))[[performanceType]])
            }, measurements, outcomes, SIMPLIFY = FALSE)
        })
        
        # Compile the performance metrics into a matrix
        realPerformance <- matrix(unlist(performanceAllPairs), ncol = length(measurements), byrow = TRUE,
                                  dimnames = list(paste("Select", names(measurements)), paste("Cross-validate", names(measurements))))
        realPerformance <- round(realPerformance, 2)
        result$real <- realPerformance
    }
    
    # If random features are requested, perform cross-validation using randomly selected features
    if (doRandomFeatures == TRUE){
        if (verbose > 0) message("Starting random feature selection procedure.")
        # Randomly select nFeatures from each dataset
        randomFeatures <- lapply(measurements, function(dataset) sample(colnames(dataset), nFeatures))
        performanceAllPairs <- lapply(randomFeatures, function(randomFeaturesSet)
        {
            mapply(function(testData, testOutcomes)
            {
                result <- crossValidate(testData[, randomFeaturesSet], testOutcomes,
                                        nFeatures = nFeatures,
                                        selectionMethod = "none",
                                        classifier = classifier,
                                        multiViewMethod = "none",
                                        nFolds = nFolds,
                                        nCores = nCores,
                                        nRepeats = nRepeats)
                mean(performance(calcCVperformance(result, performanceType))[[performanceType]])
            }, measurements, outcomes)
        })
        
        # Compile the performance metrics into a matrix
        randomPerformance <- matrix(unlist(performanceAllPairs), ncol = length(measurements), byrow = TRUE,
                                    dimnames = list(paste("Random Select", names(measurements)), paste("Cross-validate", names(measurements))))
        randomPerformance <- round(randomPerformance, 2)
        
        result$random <- randomPerformance
    }
    
    # If runTOP is TRUE, perform logistic regression using a leave-one-dataset-out approach
    if (runTOP == TRUE){
        if (verbose > 0) message("Starting logistic regression using leave-one-dataset-out approach.")
        logisticPerformance <- vector("numeric", length(measurements))
        names(logisticPerformance) <- names(measurements)
        for (i in seq_along(measurements)){
            # Combine all datasets except the i-th one for training
            trainMeasurements <- do.call(rbind, measurements[-i])
            trainOutcomes <- unlist(outcomes[-i])
            # Left-out test data
            testMeasurements <- measurements[[i]]
            testOutcomes <- outcomes[[i]]
            # Train logistic regression model
            logisticModel <- glm(trainOutcomes ~ ., data = as.data.frame(trainMeasurements), family = binomial())
            # Predict probabilities on the test data
            predictedProbabilities <- predict(logisticModel, newdata = as.data.frame(testMeasurements), type = "response")
            # Convert probabilities to class labels based on a threshold of 0.5
            predictedClasses <- ifelse(predictedProbabilities > 0.5, levels(factor(trainOutcomes))[2], levels(factor(trainOutcomes))[1])
            # Calculate performance
            logisticPerformance[i] <- calcExternalPerformance(testOutcomes, predictedClasses, performanceType)
        }
        # Convert logisticPerformance to a matrix to match the format of realPerformance
        logisticPerformanceMatrix <- matrix(logisticPerformance, nrow = 1,
                                            dimnames = list("Logistic Regression", names(measurements)))
        result$logisticRegression <- logisticPerformanceMatrix
    }
    
    # Store the parameters used during execution for reference
    result$params <- list(nFeatures = nFeatures, selectionMethod = selectionMethod,
                          selectionOptimisation = selectionOptimisation,
                          classifier = classifier, nFolds = nFolds, nRepeats = nRepeats, nCores = nCores,
                          trainType = trainType, performanceType = performanceType,
                          doRandomFeatures = doRandomFeatures, runTOP = runTOP)
    
    result
}

#' A function to plot the output of the crissCrossValidate function
#'
#' This function generates a heatmap visualization of the results from the \code{crissCrossValidate} function.
#' If \code{runTOP} was \code{TRUE}, the logistic regression results are included in the heatmap.
#'
#' @param crissCrossResult The output of the \code{crissCrossValidate} function.
#' @param includeValues If \code{TRUE}, then the values of the matrix will be included in the plot.
#' @author Harry Robertson
#'
#' @import ggplot2
#' @import reshape2
#' @import ggpubr
#'
#' @export

crissCrossPlot <- function(crissCrossResult, includeValues = FALSE){
    
    # Extract parameters and results for convenience
    params <- crissCrossResult$params
    scalebar_title <- params$performanceType
    
    # Depending on the trainType, plot differently
    if (params$trainType == "modelTrain"){
        # When trainType is "modelTrain", plot the performance matrix of models trained on each dataset and tested on all datasets
        real <- crissCrossResult$real
        
        # If logistic regression results are present, append them to the matrix
        if (params$runTOP == TRUE){
            real <- rbind(real, crissCrossResult$logisticRegression)
        }
        
        # Melt the performance matrix for ggplot
        melted_cormat <- reshape2::melt(real, na.rm = TRUE)
        
        # Generate the heatmap
        ggheatmap <- ggplot(melted_cormat, aes(Var1, Var2, fill = value)) +
            geom_tile(color = "white") +
            scale_fill_gradient2(high = "red", mid = "white", low = "blue",
                                 midpoint = 0.5, limit = c(0,1), space = "Lab",
                                 name=as.character(scalebar_title)) +
            theme_bw() + xlab("Training Dataset") + ylab("Testing Dataset") +
            theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
            theme(axis.text.y = element_text(vjust = 1, size = 8, hjust = 1)) +
            coord_fixed()
        
        if (includeValues == TRUE) ggheatmap <- ggheatmap + geom_text(aes(label = value), color = "black", size = 3)
    }
    
    else if (params$trainType == "modelTest"){
        # When trainType is "modelTest", plot the performance matrix using features selected from each dataset and cross-validated on test datasets
        real <- crissCrossResult$real
        
        # If logistic regression results are present, append them to the matrix
        if (params$runTOP == TRUE){
            real <- rbind(real, crissCrossResult$logisticRegression)
        }
        
        # Melt the performance matrix for ggplot
        melted_cormat_1 <- reshape2::melt(real, na.rm = TRUE)
        
        # Generate the first heatmap
        ggheatmap_1 <- ggplot(melted_cormat_1, aes(Var1, Var2, fill = value)) +
            geom_tile(color = "white") +
            scale_fill_gradient2(high = "red", mid = "white", low = "blue",
                                 midpoint = 0.5, limit = c(0,1), space = "Lab",
                                 name=as.character(scalebar_title)) +
            theme_bw() + xlab("Features Extracted") + ylab("Dataset Tested") +
            theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
            theme(axis.text.y = element_text(vjust = 1, size = 8, hjust = 1)) +
            coord_fixed()
        if (includeValues == TRUE) ggheatmap_1 <- ggheatmap_1 + geom_text(aes(label = value), color = "black", size = 3)
        
        if (params$doRandomFeatures == TRUE){
            # If random features were used, create a second heatmap for comparison
            random <- crissCrossResult$random
            melted_cormat_2 <- reshape2::melt(random, na.rm = TRUE)
            ggheatmap_2 <- ggplot(melted_cormat_2, aes(Var1, Var2, fill = value)) +
                geom_tile(color = "white") +
                scale_fill_gradient2(high = "red", mid = "white", low = "blue",
                                     midpoint = 0.5, limit = c(0,1), space = "Lab",
                                     name=as.character(scalebar_title)) +
                theme_bw() + xlab("Features Extracted") + ylab("Dataset Tested") +
                theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
                theme(axis.text.y = element_text(vjust = 1, size = 8, hjust = 1)) +
                coord_fixed()
            if (includeValues == TRUE) ggheatmap_2 <- ggheatmap_2 + geom_text(aes(label = value), color = "black", size = 3)
            
            # Arrange both heatmaps side by side for comparison
            ggheatmap <- ggarrange(ggheatmap_1, ggheatmap_2, labels = c("A - Feature Selection", "B - Random Features"),
                                   ncol = 2, common.legend = TRUE, legend = "right")
        } else {
            ggheatmap <- ggheatmap_1
        }
    }
    print(ggheatmap)
}