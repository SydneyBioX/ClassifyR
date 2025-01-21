# An Interface for sparsediscrim Package's dlda Function. Diagonal linear discriminant analysis.
DLDAtrainInterface <- function(measurementsTrain, classesTrain, verbose = 3)
{
  #if(!requireNamespace("sparsediscrim", quietly = TRUE))
  #stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting DLDA classifier to data.")
  
  # sparsediscrim::dlda(as.matrix(measurements), classes)
  .dlda(as.matrix(measurementsTrain), classesTrain)
}
attr(DLDAtrainInterface, "name") <- "DLDAtrainInterface"

# model is of class dlda.
DLDApredictInterface <- function(model, measurementsTest, returnType = c("both", "class", "score"),
                                 verbose = 3)
{
  isNumeric <- sapply(measurementsTest, is.numeric)
  measurementsTest <- measurementsTest[, isNumeric, drop = FALSE]
  returnType <- match.arg(returnType)
  
  # sparsediscrim doesn't match feature names to those inside trained model.
  # Ensure that there is no chance of mismatched columns.
  measurementsTest <- measurementsTest[, names(model[["var_pool"]])]
  
  #if(!requireNamespace("sparsediscrim", quietly = TRUE)) # Removed from CRAN, sadly.
  #stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained DLDA classifier.")
  
  #predict(model, as.matrix(test))
  predictions <- predict(model, as.matrix(measurementsTest)) # Copy located in utilities.R.

  switch(returnType, class = predictions[["class"]], # Factor vector.
         score = predictions[["posterior"]][, model[["groups"]]], # Numeric matrix.
         both = data.frame(class = predictions[["class"]], predictions[["posterior"]], check.names = FALSE))
}