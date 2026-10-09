# An Interface for class Package's knn Function

kNNinterface <- function(measurementsTrain, classesTrain, measurementsTest, k = 1,
                         mode = c("unweighted", "weighted"), returnType = c("both", "class", "score"), verbose = 3)
{
  # Ensure same ordering for both tables.
  measurementsTest <- measurementsTest[, colnames(measurementsTrain), drop = FALSE]
  returnType <- match.arg(returnType)
  mode <- match.arg(mode)
  
  if(!requireNamespace("BiocNeighbors", quietly = TRUE))
    stop("The package 'BiocNeighbors' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ":Fitting k Nearest Neighbours classifier to data and predicting classes.")
  
  nearestToEach <- BiocNeighbors::queryKNN(as.matrix(measurementsTrain), as.matrix(measurementsTest), k = k)
  # Matrix, rows are test samples, columns are the k nearest training samples.
  nearestClasses <- matrix(as.character(classesTrain)[nearestToEach[["index"]]], nrow = nrow(nearestToEach[["index"]]))
  if(mode == "unweighted")
  {
      classScores <- t(apply(nearestClasses, 1, function(nearestRow) table(factor(nearestRow, levels = levels(classesTrain))) / length(nearestRow)))
      classPredictions <- levels(classesTrain)[apply(classScores, 1, which.max)]
  } else { # Mode is weighted. Euclidean distance contributes to a sample's influence on final prediction.
      classScores <- t(mapply(function(classes, distances) sapply(levels(classesTrain), function(class) sum(1/distances[classes == class])), split(nearestClasses, 1:nrow(nearestClasses)), split(nearestToEach[["distance"]], 1:nrow(nearestToEach[["distance"]]))))
      if(any(is.infinite(classScores)))
          classScores[is.infinite(classScores)] <- 999999
      classScores <- t(apply(classScores, 1, function(row) row / sum(row)))
      classPredictions <- levels(classesTrain)[apply(classScores, 1, which.max)]
  }
  
  classPredictions <- factor(classPredictions, levels = levels(classesTrain))
  rownames(classScores) <- names(classPredictions) <- rownames(measurementsTest)
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}
attr(kNNinterface, "name") <- "kNNinterface"