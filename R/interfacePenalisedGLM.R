# An Interface for glmnet Package's glmnet Function. Generalised linear models with sparsity.

penalisedGLMtrainInterface <- function(measurementsTrain, classesTrain, lambda = NULL, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting elastic net regularised GLM classifier to data.")

  # One-hot encoding needed.    
  measurementsTrain <- .encodeTrain(measurementsTrain)
  fitted <- glmnet::glmnet(measurementsTrain, classesTrain, family = "multinomial", lambda = lambda,
                           weights = as.numeric(1 / (table(classesTrain)[classesTrain] / length(classesTrain))), ...)
  # Inverse class size weighting needed to give decent predictions when class imbalance.
  
  if(is.null(lambda) || length(lambda) > 1) # fitted has numerous models for a range of lambda values.
  { # Pick one lambda based on resubstitution performance. But not the one that makes all variables excluded from model.
    lambdaConsider <- fitted[["lambda"]][colSums(as.matrix(fitted[["beta"]][[1]])) != 0]
    # Predictions for all lambda values at once. A column for each lambda.
    lambdasPredictions <- as.matrix(predict(fitted, measurementsTrain, s = lambdaConsider, type = "class"))
    balancedErrors <- apply(lambdasPredictions, 2, function(lambdaPredictions)
    {
      classPredictions <- factor(as.character(lambdaPredictions), levels = fitted[["classnames"]])
      calcExternalPerformance(classesTrain, classPredictions, "Balanced Error")
    })
    bestLambda <- lambdaConsider[which.min(balancedErrors)[1]] # Largest Lambda with minimum balanced error rate.
    attr(fitted, "tune") <- list(lambda = bestLambda)
  } else { # The user specified exactly one lambda value. Record it.
    attr(fitted, "tune") <- list(lambda = lambda)
  }
  
  attr(fitted, "featureNames") <- colnames(measurementsTrain)
  attr(fitted, "featureGroups") <- attr(measurementsTrain, "assign")
  attr(fitted, "encoding") <- attr(measurementsTrain, "encoding")
  
  fitted
}
attr(penalisedGLMtrainInterface, "name") <- "penalisedGLMtrainInterface"

# model is of class multnet
penalisedGLMpredictInterface <- function(model, measurementsTest, lambda, ..., returnType = c("both", "class", "score"), verbose = 3)
{# ... just consumes emitted tuning variables from .doTrain which are unused.
  returnType <- match.arg(returnType)

  # One-hot encoding needed.
  # Ensure that testing data has same columns names in same order as training data.
  measurementsTest <- .encodeTest(measurementsTest, model)
  
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained elastic net regularised GLM classifier.")

  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.

  classPredictions <- factor(as.character(predict(model, measurementsTest, s = lambda, type = "class")), levels = model[["classnames"]])
  classScores <- predict(model, measurementsTest, s = lambda, type = "response")[, , 1]
  
  if(is.matrix(classScores))
    classScores <- classScores[, model[["classnames"]]]
  else # Leave-one-out cross-validation likely used and glmnet doesn't have consistent return types.
    classScores <- t(classScores[model[["classnames"]]])
  
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}

################################################################################
#
# One-hot encoding of categorical features for glmnet and xgboost, which need a
# numeric matrix. The encoding of the training data is stored with the model so
# that test data are encoded into the same columns, in the same order, with the
# same factor levels, whichever levels are present in the test samples.
#
################################################################################

.encodeTrain <- function(measurementsTrain)
{
  measurementsTrain <- as(measurementsTrain, "data.frame")
  isCategorical <- sapply(measurementsTrain, function(featureValues) is.factor(featureValues) || is.character(featureValues))
  featuresLevels <- lapply(measurementsTrain[isCategorical], function(featureValues) levels(factor(featureValues)))
  trainMatrix <- model.matrix(~ 0 + ., data = measurementsTrain, xlev = featuresLevels)
  attr(trainMatrix, "encoding") <- list(features = colnames(measurementsTrain), levels = featuresLevels,
                                        columns = colnames(trainMatrix))
  trainMatrix
}

# model has an "encoding" attribute made by .encodeTrain.
.encodeTest <- function(measurementsTest, model)
{
  encoding <- attr(model, "encoding")
  # The features in the training order, so that each factor is encoded with the same contrasts.
  measurementsTest <- as(measurementsTest, "data.frame")[, encoding[["features"]], drop = FALSE]
  # Keep samples with missing values, so that each prediction stays with its sample.
  testFrame <- model.frame(~ 0 + ., data = measurementsTest, xlev = encoding[["levels"]], na.action = na.pass)
  testMatrix <- model.matrix(attr(testFrame, "terms"), testFrame)
  missingColumns <- setdiff(encoding[["columns"]], colnames(testMatrix))
  if(length(missingColumns) > 0)
    testMatrix <- cbind(testMatrix, matrix(0, nrow(testMatrix), length(missingColumns), dimnames = list(NULL, missingColumns)))
  testMatrix[, encoding[["columns"]], drop = FALSE]
}

################################################################################
#
# Get selected features (i.e. non-zero model coefficients)
#
# Note: Need to convert back to actual features when factors were expanded into
# multiple columns using indicator variables.
#
################################################################################

penalisedFeatures <- function(model)
                      {
                        # Floating point numbers test for equality.
                        whichCoefficientColumn <- which(abs(model[["lambda"]] - attr(model, "tune")[["lambda"]]) < 0.00001)[1]
                        if(is.list(model[["beta"]])) # Categorical data
                        {
                          coefficientsUsed <- sapply(model[["beta"]], function(classCoefficients) classCoefficients[, whichCoefficientColumn])
                          featureScores <- rowSums(abs(coefficientsUsed))
                        } else { # survival data
                            featureScores <- abs(model[["beta"]][, whichCoefficientColumn])
                        }
                        featureGroups <- attr(model, "featureGroups")[match(names(featureScores), attr(model, "featureNames"))]
                        groupScores <- unname(by(featureScores, featureGroups, max))
                        rankedFeaturesIndices <- order(groupScores, decreasing = TRUE)
                        selectedFeaturesIndices <- which(groupScores != 0)
                        list(rankedFeaturesIndices, selectedFeaturesIndices)
                      }