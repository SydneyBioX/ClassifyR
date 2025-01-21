# An Interface for Ordinary Generalised Linear Models (GLMs) with Binary Outcome
GLMtrainInterface <- function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  if(verbose == 3)
    message(Sys.time(), ": Fitting GLM classifier to data.")
  fitData <- cbind(measurementsTrain, class = classesTrain)
  glm(class ~ . + 0, family = quasibinomial, data = fitData, weights = as.numeric(1 / (table(classesTrain)[classesTrain] / length(classesTrain))), ...)
}
attr(GLMtrainInterface, "name") <- "GLMtrainInterface"

# model is of class glm.
GLMpredictInterface <- function(model, measurementsTest, returnType = c("both", "class", "score"),
                                verbose = 3)
{
  returnType <- match.arg(returnType)
  
  if(verbose == 3)
    message("Predicting classes using trained GLM classifier.")

  predictions <- predict(model, measurementsTest, type = "response")
  classes <- levels(model[["model"]][["class"]])
  classPredictions <- factor(ifelse(predictions >= 0.5, classes[2], classes[1]), classes)
  classScores <- matrix(c(1 - predictions, predictions), ncol = 2)
  colnames(classScores) <- classes
  
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}
