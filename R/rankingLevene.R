# Selection of Differential Variability with Levene Statistic
leveneRanking <- function(measurementsTrain, classesTrain, verbose = 3)
{
  if(!requireNamespace("car", quietly = TRUE))
    stop("The package 'car' could not be found. Please install it.")            
  if(verbose == 3)
    message(Sys.time(), ": Calculating Levene statistic.")

  pValues <- apply(measurementsTrain, 2, function(featureColumn)
             car::leveneTest(featureColumn, classesTrain)[["Pr(>F)"]][1])
  
  order(pValues) # From smallest to largest.
}
attr(leveneRanking, "name") <- "leveneRanking"