# Interface for \code{pamr.train} Function from \code{pamr} CRAN Package
NSCtrainInterface <- function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")

  # The features are expected to be in the rows and the samples are in the columns.
  trainedModel <- pamr::pamr.train(list(x = t(as.matrix(measurementsTrain)), y = classesTrain), ...)
  
  if(verbose == 3)
    message(Sys.time(), ": Nearest shrunken centroid training completed.")
  
  trainedModel  
}
attr(NSCtrainInterface, "name") <- "NSCtrainInterface"

NSCpredictInterface <- function(model, measurementsTest, classesColumnTest = NULL, ...,
                                returnType = c("both", "class", "score"), verbose = 3)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  returnType <- match.arg(returnType)
  
  if(!is.null(classesColumnTest)) # Remove the column, since pamr uses positional matching of features.
  {
    if(is.character(classesColumnTest)) classesColumnTest <- match(classesColumnTest, colnames(measurementsTest))      
    measurementsTest <- measurementsTest[, -classesColumnTest]      
  }
  
  minError <- min(model[["errors"]])
  threshold <- model[["threshold"]][max(which(model[["errors"]] == minError))]
  
  measurementsTest <- t(as.matrix(measurementsTest))   
  classPredictions <- pamr::pamr.predict(model, measurementsTest, threshold, ...)
  classScores <- pamr::pamr.predict(model, measurementsTest, threshold, type = "posterior", ...)[, levels(model[["y"]])]
  if(!is.matrix(classScores)) # Only one sample was predicted and pamr isn't consistent with return types.
    classScores <- t(classScores)
  
  if(verbose == 3)
    message("Nearest shrunken centroid predictions made.")
  
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}


################################################################################
#
# Get selected features (pamr.listgenes)
#
################################################################################

# model is of class pamrtrained
NSCfeatures <- function(model, measurementsTrain, classesTrain)
               {
                 if(!requireNamespace("pamr", quietly = TRUE))
                   stop("The package 'pamr' could not be found. Please install it.")
            
                 minError <- min(model[["errors"]])
                 threshold <- model[["threshold"]][max(which(model[["errors"]] == minError))]
                 params <- c(list(model), list(list(x = t(as.matrix(measurementsTrain)), y = measurementsTrain, geneid = 1:ncol(measurementsTrain))), threshold)
                 chosenIndices <- as.numeric(do.call(pamr::pamr.listgenes, params)[, 1])
            
                 list(NULL, chosenIndices)
               }