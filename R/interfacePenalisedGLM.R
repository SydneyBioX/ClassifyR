# An Interface for glmnet Package's glmnet Function. Generalised linear models with sparsity.

penalisedGLMtrainInterface <- function(measurementsTrain, classesTrain, lambda = NULL, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting elastic net regularised GLM classifier to data.")

  # One-hot encoding needed.    
  measurementsTrain <- MatrixModels::model.Matrix(~ 0 + ., data = measurementsTrain)
  fitted <- glmnet::glmnet(measurementsTrain, classesTrain, family = "multinomial", weights = as.numeric(1 / (table(classesTrain)[classesTrain] / length(classesTrain))), ...)
  # Inverse class size weighting needed to give decent predictions when class imbalance.
  
  if(is.null(lambda)) # fitted has numerous models for automatically chosen lambda values.
  { # Pick one lambda based on resubstitution performance. But not the one that makes all variables excluded from model.
    lambdaConsider <- colSums(as.matrix(fitted[["beta"]][[1]])) != 0
    bestLambda <- fitted[["lambda"]][lambdaConsider][which.min(sapply(fitted[["lambda"]][lambdaConsider], function(lambda) # Largest Lambda with minimum balanced error rate.
    {
      classPredictions <- factor(as.character(predict(fitted, measurementsTrain, s = lambda, type = "class")), levels = fitted[["classnames"]])
      calcExternalPerformance(classesTrain, classPredictions, "Balanced Error")
    }))[1]]
    attr(fitted, "tune") <- list(lambda = bestLambda)
  }
  
  attr(fitted, "featureNames") <- colnames(measurementsTrain)
  attr(fitted, "featureGroups") <- measurementsTrain@assign
  
  fitted
}
attr(penalisedGLMtrainInterface, "name") <- "penalisedGLMtrainInterface"

# model is of class multnet
penalisedGLMpredictInterface <- function(model, measurementsTest, lambda, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  # ... just consumes emitted tuning variables from .doTrain which are unused.
  returnType <- match.arg(returnType)
  # One-hot encoding needed.
  measurementsTest <- MatrixModels::model.Matrix(~ 0 + ., data = measurementsTest)
  
  # Ensure that testing data has same columns names in same order as training data.
  # Remove those annoying backquotes which glmnet adds if variables have spaces in names.
  measurementsTest <- measurementsTest[, gsub('`', '', rownames(model[["beta"]][[1]]))]
  
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained elastic net regularised GLM classifier.")

  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.

  
  measurementsTest <- measurementsTest[, rownames(model[["beta"]][[1]])]

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
                            featureScores <- model[["beta"]][, whichCoefficientColumn]
                        }
                        featureGroups <- attr(model, "featureGroups")[match(names(featureScores), attr(model, "featureNames"))]
                        groupScores <- unname(by(featureScores, featureGroups, max))
                        rankedFeaturesIndices <- order(groupScores, decreasing = TRUE)
                        selectedFeaturesIndices <- which(groupScores != 0)
                        list(rankedFeaturesIndices, selectedFeaturesIndices)
                      }