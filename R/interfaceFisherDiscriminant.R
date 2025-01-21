# Classification Using Fisher's LDA. Unlike ordinary LDA, Fisher's version does not have assumptions
# about the normality of the features.

fisherDiscriminant <- function(measurementsTrain, classesTrain, measurementsTest,
                               returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)
  
  if(verbose == 3)
    message(Sys.time(), ": Fitting Fisher discriminant.")              

  oneClassTraining <- which(classesTrain == levels(classesTrain)[1])
  otherClassTraining <- which(classesTrain == levels(classesTrain)[2])
  varOneClass <- apply(trainingMatrix[oneClassTraining, ], 2, var)
  varOtherClass <- apply(trainingMatrix[otherClassTraining, ], 2, var)
  varAll <- ((length(varOneClass) - 1) * varOneClass + (length(varOtherClass) - 1)
             * varOtherClass) / (length(oneClassTraining) + length(otherClassTraining) - 2)
  aT <- (apply(trainingMatrix[oneClassTraining, ], 2, mean) - apply(trainingMatrix[otherClassTraining, ], 2, mean)) / varAll
  criticalValue <- 0.5 * aT %*% as.matrix(apply(trainingMatrix[oneClassTraining, ], 2, mean) +
                                          apply(trainingMatrix[otherClassTraining, ], 2, mean))
  
  if(verbose == 3)
    message("Critical value calculated.")
  
  classesPredicted <- factor(apply(measurementsTest, 1, function(testSample)
  {
    if(aT %*% as.matrix(testSample) >= criticalValue)
      levels(classesTrain)[1]
    else
      levels(classesTrain)[2]
  }), levels = levels(classesTrain))
  scores <- apply(measurementsTest, 1, function(testSample) -1 * (aT %*% as.matrix(testSample))) # In reference to the second level of 'classes'. 
  
  switch(returnType, class = classesPredicted,
                     score = scores,
                     both = data.frame(class = classesPredicted, score = scores, check.names = FALSE))  
}
attr(fisherDiscriminant, "name") <- "fisherDiscriminant"