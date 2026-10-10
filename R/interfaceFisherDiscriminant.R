# Classification Using Fisher's LDA. Unlike ordinary LDA, Fisher's version does not have assumptions
# about the normality of the features.

fisherDiscriminant <- function(measurementsTrain, classesTrain, measurementsTest,
                               returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)
  if(nlevels(classesTrain) != 2)
    stop("Fisher discriminant is for two classes but 'classesTrain' has ", nlevels(classesTrain), " levels.")
  
  if(verbose == 3)
    message(Sys.time(), ": Fitting Fisher discriminant.")              

  trainingMatrix <- as.matrix(measurementsTrain)
  testingMatrix <- as.matrix(measurementsTest[, colnames(trainingMatrix), drop = FALSE])
  oneClassTraining <- which(classesTrain == levels(classesTrain)[1])
  otherClassTraining <- which(classesTrain == levels(classesTrain)[2])
  varOneClass <- apply(trainingMatrix[oneClassTraining, , drop = FALSE], 2, var)
  varOtherClass <- apply(trainingMatrix[otherClassTraining, , drop = FALSE], 2, var)
  # Pooled variance of each feature, weighted by the class sizes.
  varAll <- ((length(oneClassTraining) - 1) * varOneClass + (length(otherClassTraining) - 1)
             * varOtherClass) / (length(oneClassTraining) + length(otherClassTraining) - 2)
  aT <- (apply(trainingMatrix[oneClassTraining, , drop = FALSE], 2, mean) - apply(trainingMatrix[otherClassTraining, , drop = FALSE], 2, mean)) / varAll
  criticalValue <- 0.5 * aT %*% as.matrix(apply(trainingMatrix[oneClassTraining, , drop = FALSE], 2, mean) +
                                          apply(trainingMatrix[otherClassTraining, , drop = FALSE], 2, mean))
  
  if(verbose == 3)
    message("Critical value calculated.")
  
  projections <- apply(testingMatrix, 1, function(testSample) aT %*% as.matrix(testSample))
  classesPredicted <- factor(ifelse(projections >= as.vector(criticalValue), levels(classesTrain)[1], levels(classesTrain)[2]),
                             levels = levels(classesTrain))
  # In reference to the second level of 'classes'. Centred on the critical value, so a positive score
  # predicts the second class and scores from models trained on different samples are comparable.
  scores <- as.vector(criticalValue) - projections
  
  switch(returnType, class = classesPredicted,
                     score = scores,
                     both = data.frame(class = classesPredicted, score = scores, check.names = FALSE))  
}
attr(fisherDiscriminant, "name") <- "fisherDiscriminant"