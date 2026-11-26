# An Interface for mboost Package's glmboost Function with the CoxPH Family. Componentwise boosted Cox models.

boostedCoxTrainInterface <- function(measurementsTrain, survivalTrain, mstop = 100, nu = 0.1, ..., verbose = 3)
{
  if(!requireNamespace("mboost", quietly = TRUE))
    stop("The package 'mboost' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting boosted Cox model to data.")
  measurementsMatrix <- .encodeTrain(measurementsTrain) # One-hot encoding needed.
  # A Cox model has no intercept, which glmboost warns about when the features are centred.
  fitted <- withCallingHandlers(
    mboost::glmboost(measurementsMatrix, survivalTrain, family = mboost::CoxPH(),
                     control = mboost::boost_control(mstop = mstop, nu = nu), center = TRUE, ...),
    warning = function(warning) if(grepl("does not contain intercept", conditionMessage(warning))) invokeRestart("muffleWarning"))
  attr(fitted, "featureNames") <- colnames(measurementsMatrix)
  attr(fitted, "featureGroups") <- attr(measurementsMatrix, "assign")
  attr(fitted, "encoding") <- attr(measurementsMatrix, "encoding")
  fitted
}
attr(boostedCoxTrainInterface, "name") <- "boostedCoxTrainInterface"

# model is of class glmboost. The risk is the linear predictor: higher is riskier.
boostedCoxPredictInterface <- function(model, measurementsTest, ..., verbose = 3)
{
  if(!requireNamespace("mboost", quietly = TRUE))
    stop("The package 'mboost' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting risks using boosted Cox model.")
  testMatrix <- .encodeTest(measurementsTest, model) # Same one-hot encoding as the training data.
  setNames(as.numeric(predict(model, newdata = testMatrix, type = "link")), rownames(measurementsTest))
}

# Features ranked by the size of their coefficients, and those the boosting chose.
boostedCoxFeatures <- function(model)
{
  coefficients <- coef(model) # Only the chosen columns.
  columnScores <- setNames(numeric(length(attr(model, "featureNames"))), attr(model, "featureNames"))
  columnScores[names(coefficients)] <- abs(coefficients)
  .encodedFeatures(columnScores, model)
}

# The features (columns of the data before encoding) ranked by the largest score of their encoded columns, and those
# with a score that is not zero. model has the attributes "featureGroups" and "encoding" made from .encodeTrain.
.encodedFeatures <- function(columnScores, model)
{
  featureGroups <- attr(model, "featureGroups")[match(names(columnScores), attr(model, "featureNames"))]
  nFeatures <- length(attr(model, "encoding")[["features"]])
  featureScores <- tapply(abs(columnScores), factor(featureGroups, levels = seq_len(nFeatures)), max)
  featureScores[is.na(featureScores)] <- 0
  list(order(featureScores, decreasing = TRUE), unname(which(featureScores != 0)))
}
