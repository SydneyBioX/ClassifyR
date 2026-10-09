# An Interface for e1071 Package's Support Vector Machine Classifier.
SVMtrainInterface <- function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  if(!requireNamespace("e1071", quietly = TRUE))
    stop("The package 'e1071' could not be found. Please install it.")
  
  if(verbose == 3)
    message(Sys.time(), ": Fitting SVM classifier to data.")
  allVariables <- cbind(measurementsTrain, classesTrain) 
  trained <- e1071::svm(classesTrain ~ ., data = allVariables, probability = TRUE, ...)
  
  if(ncol(measurementsTrain) == 1) # Handle inconsistency by e1071 to not always name columns.
      colnames(trained[["SV"]]) <- colnames(measurementsTrain)
  
  trained
}
attr(SVMtrainInterface, "name") <- "SVMtrainInterface"

# model is of class svm
SVMpredictInterface <- function(model, measurementsTest, returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)

  if(!requireNamespace("e1071", quietly = TRUE))
    stop("The package 'e1071' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained SVM classifier.")
  
  # The model was fitted with a formula, so prediction on a data frame encodes the features in the same
  # way as for training, matching them by name.
  measurementsTest <- .asDataFrame(measurementsTest)
  classPredictions <- predict(model, measurementsTest, probability = TRUE)
  
  # e1071 uses attributes to pass back probabilities. Make them a standalone variable.
  classScores <- attr(classPredictions, "probabilities")[, model[["levels"]], drop = FALSE]
  attr(classPredictions, "probabilities") <- NULL
  rownames(classScores) <- names(classPredictions) <- rownames(measurementsTest)
  switch(returnType, class = classPredictions, score = classScores,
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}