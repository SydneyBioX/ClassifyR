# An Interface for LiblineaR Package's LiblineaR Function. Linear models of classes fitted by LIBLINEAR: logistic
# regression with an L2 (type 0) or L1 (type 6) penalty.

LiblineaRtrainInterface <- function(measurementsTrain, classesTrain, type = 0, cost = NULL, costs = 10^(-3:2),
                                    nFoldsCost = 5, ..., verbose = 3)
{
  if(!requireNamespace("LiblineaR", quietly = TRUE))
    stop("The package 'LiblineaR' could not be found. Please install it.")
  if(!type %in% c(0, 6, 7))
    stop("'type' must be a logistic regression type of LIBLINEAR (0, 6 or 7), which give class probabilities.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting LIBLINEAR logistic regression to data.")
  measurementsMatrix <- .encodeTrain(measurementsTrain) # One-hot encoding needed.
  # LIBLINEAR's penalty depends on the scale of the features, so they are standardised by the training samples.
  centres <- colMeans(measurementsMatrix)
  scales <- apply(measurementsMatrix, 2, sd)
  scales[is.na(scales) | scales == 0] <- 1
  scaledMatrix <- scale(measurementsMatrix, centres, scales)
  # Inverse class size weighting, so that each class contributes equally to the penalised loss.
  classWeights <- c(length(classesTrain) / (nlevels(classesTrain) * table(classesTrain)))

  if(is.null(cost)) # Choose the cost by the balanced error of cross-validation of the training samples.
  {
    nFoldsCost <- min(nFoldsCost, table(classesTrain))
    if(nFoldsCost < 3)
    {
      warning("A class has fewer than three training samples, so the cost is 1.")
      cost <- 1
    } else {
      foldsIDs <- integer(length(classesTrain)) # Folds stratified by class.
      for(classIndices in split(seq_along(classesTrain), classesTrain))
        foldsIDs[classIndices] <- sample(rep_len(seq_len(nFoldsCost), length(classIndices)))
      balancedErrors <- sapply(costs, function(aCost)
      {
        predictions <- factor(rep(NA, length(classesTrain)), levels = levels(classesTrain))
        for(fold in seq_len(nFoldsCost))
        {
          inFold <- foldsIDs == fold
          fitted <- LiblineaR::LiblineaR(scaledMatrix[!inFold, , drop = FALSE], classesTrain[!inFold], type = type,
                                         cost = aCost, wi = classWeights, ...)
          predictions[inFold] <- as.character(predict(fitted, scaledMatrix[inFold, , drop = FALSE])[["predictions"]])
        }
        calcExternalPerformance(classesTrain, predictions, "Balanced Error")
      })
      cost <- costs[which.min(balancedErrors)] # The smallest cost (strongest penalty) among ties.
    }
  }
  fitted <- LiblineaR::LiblineaR(scaledMatrix, classesTrain, type = type, cost = cost, wi = classWeights, ...)
  fitted[["classes"]] <- levels(classesTrain)
  attr(fitted, "tune") <- list(cost = cost)
  attr(fitted, "scaling") <- list(centres = centres, scales = scales)
  attr(fitted, "featureNames") <- colnames(measurementsMatrix)
  attr(fitted, "featureGroups") <- attr(measurementsMatrix, "assign")
  attr(fitted, "encoding") <- attr(measurementsMatrix, "encoding")
  fitted
}
attr(LiblineaRtrainInterface, "name") <- "LiblineaRtrainInterface"

# model is of class LiblineaR.
LiblineaRpredictInterface <- function(model, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  if(!requireNamespace("LiblineaR", quietly = TRUE))
    stop("The package 'LiblineaR' could not be found. Please install it.")
  returnType <- match.arg(returnType)
  if(verbose == 3)
    message("Predicting using LIBLINEAR logistic regression.")
  testMatrix <- .encodeTest(measurementsTest, model) # Same one-hot encoding as the training data.
  scaling <- attr(model, "scaling")
  testMatrix <- scale(testMatrix, scaling[["centres"]], scaling[["scales"]])
  predicted <- predict(model, testMatrix, proba = TRUE)
  classes <- model[["classes"]]
  classScores <- predicted[["probabilities"]][, classes, drop = FALSE]
  rownames(classScores) <- rownames(measurementsTest)
  classPredictions <- factor(as.character(predicted[["predictions"]]), levels = classes)
  switch(returnType, class = classPredictions,
         score = classScores,
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}

# Features ranked by the largest size of their weights over the classes, and those with a weight that is not zero.
LiblineaRfeatures <- function(model)
{
  weights <- model[["W"]][, colnames(model[["W"]]) != "Bias", drop = FALSE]
  columnScores <- setNames(apply(abs(weights), 2, max), attr(model, "featureNames"))
  .encodedFeatures(columnScores, model)
}
