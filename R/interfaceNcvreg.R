# An Interface for ncvreg Package's cv.ncvreg and cv.ncvsurv Functions. Logistic or Cox regression with a
# non-convex penalty (MCP or SCAD) or the lasso.

ncvregTrainInterface <- function(measurementsTrain, outcomeTrain, penalty = c("MCP", "SCAD", "lasso"), nfolds = 5, ...,
                                 verbose = 3)
{
  if(!requireNamespace("ncvreg", quietly = TRUE))
    stop("The package 'ncvreg' could not be found. Please install it.")
  penalty <- match.arg(penalty)
  if(verbose == 3)
    message(Sys.time(), ": Fitting ", penalty, " penalised regression to data.")
  measurementsMatrix <- .encodeTrain(measurementsTrain) # One-hot encoding needed.
  # lambda is chosen by cross-validation of the training samples (lambda.min), with ncvreg's default measure: the
  # deviance.
  if(is(outcomeTrain, "Surv"))
  {
    fitted <- ncvreg::cv.ncvsurv(measurementsMatrix, outcomeTrain, penalty = penalty, nfolds = nfolds, ...)
  } else {
    if(nlevels(outcomeTrain) != 2)
      stop("ncvreg classifies two classes only, but there are ", nlevels(outcomeTrain), ".")
    fitted <- ncvreg::cv.ncvreg(measurementsMatrix, as.integer(outcomeTrain) - 1L, family = "binomial", penalty = penalty,
                                nfolds = nfolds, ...)
    attr(fitted, "classes") <- levels(outcomeTrain)
  }
  attr(fitted, "featureNames") <- colnames(measurementsMatrix)
  attr(fitted, "featureGroups") <- attr(measurementsMatrix, "assign")
  attr(fitted, "encoding") <- attr(measurementsMatrix, "encoding")
  fitted
}
attr(ncvregTrainInterface, "name") <- "ncvregTrainInterface"

# model is of class cv.ncvreg or cv.ncvsurv, fitted at lambda.min.
ncvregPredictInterface <- function(model, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  if(!requireNamespace("ncvreg", quietly = TRUE))
    stop("The package 'ncvreg' could not be found. Please install it.")
  returnType <- match.arg(returnType)
  if(verbose == 3)
    message("Predicting using penalised regression.")
  testMatrix <- .encodeTest(measurementsTest, model) # Same one-hot encoding as the training data.
  if(is(model, "cv.ncvsurv")) # The linear predictor: higher is riskier.
    return(setNames(as.numeric(predict(model, testMatrix, type = "link")), rownames(measurementsTest)))

  classes <- attr(model, "classes")
  secondClassProbability <- as.numeric(predict(model, testMatrix, type = "response"))
  classScores <- matrix(c(1 - secondClassProbability, secondClassProbability), ncol = 2,
                        dimnames = list(rownames(measurementsTest), classes))
  classPredictions <- factor(classes[(secondClassProbability > 0.5) + 1], levels = classes)
  switch(returnType, class = classPredictions,
         score = classScores,
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}

# Features ranked by the size of their coefficients at lambda.min, and those with a coefficient that is not zero.
ncvregFeatures <- function(model)
{
  coefficients <- coef(model) # At lambda.min.
  coefficients <- coefficients[names(coefficients) != "(Intercept)"]
  .encodedFeatures(coefficients, model)
}
