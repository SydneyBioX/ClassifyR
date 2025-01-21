# Ranking of Differential Variability with Bartlett Statistic
bartlettRanking <- function(measurementsTrain, classesTrain, verbose = 3)
{
  if(verbose == 3)
    message(Sys.time(), ": Ranking features based on Bartlett statistic.")
  
  pValues <- apply(measurementsTrain, 2, function(featureColumn)
    stats::bartlett.test(featureColumn, classesTrain)[["p.value"]])
  
  order(pValues) # From smallest to largest.
}
attr(bartlettRanking, "name") <- "bartlettRanking"