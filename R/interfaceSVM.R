# An Interface for e1071 Package's Support Vector Machine Classifier.
SVMtrainInterface <- function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  if(!requireNamespace("e1071", quietly = TRUE))
    stop("The package 'e1071' could not be found. Please install it.")
    
  # Classifier requires matrix input data type.
  trainingMatrix <- as.matrix(measurementsTrain)
  
  if(verbose == 3)
    message(Sys.time(), ": Fitting SVM classifier to data.")
  
  trained <- e1071::svm(trainingMatrix, classesTrain, probability = TRUE, ...)
  
  if(ncol(trainingMatrix) == 1) # Handle inconsistency by e1071 to not always name columns.
      colnames(trained[["SV"]]) <- colnames(trainingMatrix)
  
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
  
  # Prediction function depends on test data having same set of columns in same order as
  # selected features used for training.
  colnames(measurementsTest) <- make.names(colnames(measurementsTest))
  measurementsTest <- measurementsTest[, colnames(model[["SV"]])]
  classPredictions <- predict(model, measurementsTest, probability = TRUE)
  
  # e1071 uses attributes to pass back probabilities. Make them a standalone variable.
  classScores <- attr(classPredictions, "probabilities")[, model[["levels"]], drop = FALSE]
  attr(classPredictions, "probabilities") <- NULL
  rownames(classScores) <- names(classPredictions) <- rownames(measurementsTest)
  switch(returnType, class = classPredictions, score = classScores,
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}