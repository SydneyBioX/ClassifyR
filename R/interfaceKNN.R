# An Interface for class Package's knn Function

kNNinterface <- function(measurementsTrain, classesTrain, measurementsTest, ..., verbose = 3)
{
  # Ensure same ordering for both tables.
  measurementsTest <- measurementsTest[, colnames(measurementsTrain), drop = FALSE]
  
  if(!requireNamespace("class", quietly = TRUE))
    stop("The package 'class' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ":Fitting k Nearest Neighbours classifier to data and predicting classes.")
  
  setNames(class::knn(as.matrix(measurementsTrain), as.matrix(measurementsTest), classesTrain, ...), rownames(measurementsTest))
}
attr(kNNinterface, "name") <- "kNNinterface"
