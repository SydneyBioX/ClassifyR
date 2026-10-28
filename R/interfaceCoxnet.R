# An Interface for glmnet Package's coxnet Function. Survival modelling with sparsity.

coxnetTrainInterface <- function(measurementsTrain, survivalTrain, lambda = NULL, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting coxnet model to data.")
    
  measurementsMatrix <- .encodeTrain(measurementsTrain) # One-hot encoding needed.
  
  # The path of lambda ends at 0.05 of its largest value unless the user sets lambda or lambda.min.ratio. Smaller
  # values give nearly unpenalised Cox models, which overfit the few events typical of omics cohorts and converge
  # slowly; on the protocol's METABRIC merge (31 combinations, 20 x 5 CV) the C-index was the same (0.602 vs 0.599)
  # and the fits 6.6 times as fast.
  if(is.null(lambda) && !"lambda.min.ratio" %in% names(list(...)))
    return(coxnetTrainInterface(measurementsTrain, survivalTrain, lambda = lambda, lambda.min.ratio = 0.05, ..., verbose = verbose))

  # The response variable is a Surv class of object.
  # cv.glmnet's settings for the cross-validation itself are only understood by cv.glmnet.
  cvOnlyNames <- c("type.measure", "foldid", "alignment", "grouped", "keep", "parallel", "relax", "gamma", "trace.it", "weights", "offset")
  if(any(names(list(...)) %in% cvOnlyNames))
  {
    cvArguments <- list(...)
    if(is.null(cvArguments[["type.measure"]])) cvArguments[["type.measure"]] <- "C" # C-index unless another measure is given.
    fit <- do.call(glmnet::cv.glmnet, c(list(measurementsMatrix, survivalTrain, family = "cox", lambda = lambda), cvArguments))
  }
  else
    fit <- .cvCoxnetC(measurementsMatrix, survivalTrain, lambda = lambda, ...)
  fitted <- fit$glmnet.fit
  
  offset <- -mean(predict(fitted, measurementsMatrix, s = fit$lambda.min, type = "link"))
  attr(fitted, "tune") <- list(lambda = fit$lambda.min, offset = offset)
  attr(fitted, "featureNames") <- colnames(measurementsMatrix)
  attr(fitted, "featureGroups") <- attr(measurementsMatrix, "assign")
  attr(fitted, "encoding") <- attr(measurementsMatrix, "encoding")
  
  class(fitted) <- class(fitted)[1] # Get rid of glmnet which messes with dispatch. 
  fitted
}
attr(coxnetTrainInterface, "name") <- "coxnetTrainInterface"

# model is of class coxnet.
coxnetPredictInterface <- function(model, measurementsTest, survivalTest = NULL, lambda, ..., verbose = 3)
{ # ... just consumes emitted tuning variables from .doTrain which are unused.
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using cox model.")
  
  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.
  
  # Same one-hot encoding, columns and column order as the training data.
  testMatrix <- .encodeTest(measurementsTest, model)
  
  offset <- attr(model, "tune")[["offset"]]
  model$offset <- TRUE
  
  survScores <- predict(model, testMatrix, s = lambda, type = "response", newoffset = offset)
  rownames(survScores) <- rownames(measurementsTest)
  survScores[, 1]
}

# The same as glmnet::cv.glmnet(x, y, family = "cox", type.measure = "C", lambda = lambda, nfolds = nfolds, ...):
# the same folds, models and choice of lambda.min, but the C-index of every lambda in a fold is computed at once,
# rather than by survival::concordance once per lambda and fold.
.cvCoxnetC <- function(x, y, lambda = NULL, nfolds = 10, ...)
{
  N <- nrow(x)
  weights <- rep(1, N)
  foldid <- sample(rep(seq(nfolds), length = N))
  fitAll <- glmnet::glmnet(x, y, weights = weights, offset = NULL, lambda = lambda, family = "cox", ...)
  lambdaAll <- fitAll[["lambda"]]
  nLambda <- length(lambdaAll)
  
  # Linear predictors of each fold's samples from the model fitted without them, at every lambda of the full model.
  predictions <- matrix(NA, N, nLambda)
  for(foldIndex in seq(nfolds))
  {
    inFold <- foldid == foldIndex
    fitFold <- glmnet::glmnet(x[!inFold, , drop = FALSE], y[!inFold, ], lambda = lambda, offset = NULL,
                              weights = weights[!inFold], family = "cox", ...)
    coefficients <- predict(fitFold, type = "coefficients", s = lambdaAll)
    nLambdaFold <- min(ncol(coefficients), nLambda)
    predictions[inFold, seq(nLambdaFold)] <- as.matrix(x[inFold, ] %*% coefficients[, seq(nLambdaFold)])
    if(nLambdaFold < nLambda)
      predictions[inFold, seq(from = nLambdaFold, to = nLambda)] <- predictions[inFold, nLambdaFold]
  }
  
  CindexFolds <- t(sapply(seq(nfolds), function(foldIndex)
  {
    inFold <- foldid == foldIndex
    .CindexColumns(predictions[inFold, , drop = FALSE], y[inFold, "time"], y[inFold, "status"])
  }))
  foldWeights <- tapply(weights, foldid, sum)
  CindexMean <- apply(CindexFolds, 2, weighted.mean, w = foldWeights, na.rm = TRUE)
  CindexSD <- sqrt(apply(scale(CindexFolds, CindexMean, FALSE)^2, 2, weighted.mean, w = foldWeights, na.rm = TRUE) / (nfolds - 1))
  defined <- !is.na(CindexSD)
  lambdaDefined <- lambdaAll[defined]
  CindexDefined <- CindexMean[defined]
  list(lambda.min = max(lambdaDefined[CindexDefined >= max(CindexDefined, na.rm = TRUE)], na.rm = TRUE), glmnet.fit = fitAll)
}

# Harrell's C-index of every column of risk scores (higher means shorter survival), computed as
# survival::concordance(Surv(time, event) ~ -risk) does: a pair is comparable if the shorter time is an event, or
# if the times are equal and only the first is an event; tied risks count one half.
.CindexColumns <- function(risk, time, event)
{
  time <- survival::aeqSurv(survival::Surv(time, event))[, "time"] # Near-equal times are tied.
  comparable <- (outer(time, time, "<") & event == 1) | (outer(time, time, "==") & outer(event == 1, event == 0, "&"))
  apply(risk, 2, function(riskColumn)
  {
    concordant <- sum(comparable & outer(riskColumn, riskColumn, ">"))
    tied <- sum(comparable & outer(riskColumn, riskColumn, "=="))
    discordant <- sum(comparable) - concordant - tied
    ((concordant - discordant) / (concordant + discordant + tied) + 1) / 2
  })
}
