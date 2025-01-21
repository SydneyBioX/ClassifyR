# An Interface for PoiClaClu Package's Classify Function. Poisson LDA for counts.

classifyInterface <- function(countsTrain, classesTrain, countsTest, ...,
                              returnType = c("both", "class", "score"), verbose = 3)
{
  if(!requireNamespace("PoiClaClu", quietly = TRUE))
    stop("The package 'PoiClaClu' could not be found. Please install it.")
  returnType <- match.arg(returnType)
  
  if(verbose == 3)
    message(Sys.time(), ": Fitting Poisson LDA classifier to training data and making predictions on test data.")

  predicted <- PoiClaClu::Classify(trainingMatrix, classesTrain, testingMatrix, ...)
  classPredictions <- predicted[["ytehat"]]
  classScores <- predicted[["discriminant"]]
  colnames(classScores) <- levels(classesTrain)
  rownames(classScores) <- names(classPredictions) <- rownames(countsTest)
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}
attr(classifyInterface, "name") <- "classifyInterface"