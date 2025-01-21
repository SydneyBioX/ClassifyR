# Ranking of Differential Distributions with Kolmogorov-Smirnov Distance
KolmogorovSmirnovRanking <- function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  if(verbose == 3)
    message(Sys.time(), ": Ranking features by Kolmogorov Smirnov distance between classes.")

  oneClass <- classesTrain == levels(classesTrain)[1]
  otherClass <- classesTrain == levels(classesTrain)[2]
  KSdistance <- apply(measurementsTrain, 2, function(featureColumn)
                      stats::ks.test(featureColumn[oneClass], featureColumn[otherClass], ...)[["statistic"]])

  order(KSdistance, decreasing = TRUE)
}
attr(KolmogorovSmirnovRanking, "name") <- "KolmogorovSmirnovRanking"
