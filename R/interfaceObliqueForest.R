# An Interface for aorsf Package's orsf Function. Oblique random forests of classes or survival.

# The outcome columns added to the training data.
.forestOutcome <- c(class = "ClassifyRoutcomeClass", time = "ClassifyRoutcomeTime", event = "ClassifyRoutcomeEvent")

obliqueForestTrainInterface <- function(measurementsTrain, outcomeTrain, mTryProportion = NULL, n_tree = 500,
                                        importance = "anova", n_thread = 1, ..., verbose = 3)
{
  if(!requireNamespace("aorsf", quietly = TRUE))
    stop("The package 'aorsf' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting oblique random forest to training data.")
  # Number of features in each linear combination. By default, aorsf's own: the square root of the number of features.
  extras <- list(...)
  if(!is.null(mTryProportion)) extras[["mtry"]] <- max(1, round(mTryProportion * ncol(measurementsTrain)))
  trainData <- .asDataFrame(measurementsTrain)
  featureNames <- colnames(trainData)
  if(any(.forestOutcome %in% featureNames))
    stop("Feature names must not be ", paste(.forestOutcome, collapse = ", "), ", which hold the outcome.")
  if(is(outcomeTrain, "Surv"))
  {
    trainData[[.forestOutcome[["time"]]]] <- outcomeTrain[, "time"]
    trainData[[.forestOutcome[["event"]]]] <- outcomeTrain[, "status"]
    formula <- as.formula(paste0("survival::Surv(", .forestOutcome[["time"]], ", ", .forestOutcome[["event"]], ") ~ ."))
  } else {
    trainData[[.forestOutcome[["class"]]]] <- outcomeTrain
    formula <- as.formula(paste(.forestOutcome[["class"]], "~ ."))
  }
  # Trees are grown from seeds drawn from R's random numbers, so set.seed makes the forest reproducible.
  forest <- do.call(aorsf::orsf, c(list(data = trainData, formula = formula, n_tree = n_tree, importance = importance,
                                        n_thread = n_thread, oobag_pred_type = "none"), extras))
  attr(forest, "featureNames") <- featureNames
  forest
}
attr(obliqueForestTrainInterface, "name") <- "obliqueForestTrainInterface"

# forest is of class ObliqueForest (an R6 object).
obliqueForestPredictInterface <- function(forest, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  if(!requireNamespace("aorsf", quietly = TRUE))
    stop("The package 'aorsf' could not be found. Please install it.")
  returnType <- match.arg(returnType)
  if(verbose == 3)
    message("Predicting using oblique random forest.")
  testData <- .asDataFrame(measurementsTest)[, attr(forest, "featureNames"), drop = FALSE]

  if(forest[["tree_type"]] == "survival") # Predicted mortality: higher is riskier.
    return(setNames(as.numeric(predict(forest, new_data = testData, pred_type = "mort", n_thread = 1)), rownames(testData)))

  classScores <- predict(forest, new_data = testData, pred_type = "prob", n_thread = 1)
  if(!is.matrix(classScores)) classScores <- matrix(classScores, nrow = 1, dimnames = list(NULL, names(classScores)))
  classes <- forest[["class_levels"]]
  colnames(classScores) <- classes
  rownames(classScores) <- rownames(testData)
  classPredictions <- factor(classes[max.col(classScores, ties.method = "first")], levels = classes)
  switch(returnType, class = classPredictions,
         score = classScores,
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}

# Features ranked by the forest's importance, and those with positive importance.
obliqueForestFeatures <- function(forest)
{
  featureNames <- attr(forest, "featureNames")
  importance <- tryCatch(aorsf::orsf_vi(forest), error = function(error) NULL)
  if(is.null(importance)) # Importance was not computed (importance = "none").
    return(list(seq_along(featureNames), seq_along(featureNames)))
  ranked <- match(names(importance), featureNames)
  list(ranked, ranked[importance > 0])
}
