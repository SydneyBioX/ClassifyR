# An Interface for survival Package's coxph Function. The standard Cox proportional hazards.
coxphTrainInterface <- function(measurementsTrain, survivalTrain, ..., verbose = 3)
{
  if(!requireNamespace("survival", quietly = TRUE))
    stop("The package 'survival' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting coxph classifier to training data and making predictions on test
            data.")
  
  # coxph doesn't like DataFrame input.
  measurementsTrain <- .asDataFrame(measurementsTrain)    
  survival::coxph(survivalTrain ~ ., measurementsTrain)
}
attr(coxphTrainInterface, "name") <- "coxphTrainInterface"

# model is of class coxph.
coxphPredictInterface <- function(model, measurementsTest, ..., verbose = 3)
{
  predict(model, .asDataFrame(measurementsTest), type = "risk")
}