# Classification Using k Pairs of Features With Relative Differences Between Classes
kTSPclassifier <- function(measurementsTrain, classesTrain, measurementsTest, featurePairs = NULL,
                   difference = c("unweighted", "weighted"), minDifference = 0,
                   returnType = c("both", "class", "score"), verbose = 3)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")            
  if(verbose == 3)
    message(Sys.time(), ": Determining inequalities of feature pairs.")              
  
  difference <- match.arg(difference)
  returnType <- match.arg(returnType)
  
  classesSizes <- sapply(levels(classesTrain), function(class) sum(classesTrain == class))
  largerClass <- names(classesSizes)[which.max(classesSizes)[1]]
  secondClass <- classesTrain == levels(classesTrain)[2]

  # Order pairs so that first < second is the rule for predicting the second class, based on factor levels.
  # Effectively the classifier training.
  featurePairs <- do.call(c, lapply(featurePairs, function(pair)
  {
    isSmaller <- trainingMatrix[secondClass, S4Vectors::first(pair)] < trainingMatrix[secondClass, S4Vectors::second(pair)]
    if(sum(isSmaller) < length(isSmaller) / 2)
      Pairs(S4Vectors::second(pair), S4Vectors::first(pair))
    else
      pair
  }))

  testDataFrame <- data.frame(t(testingMatrix), check.names = FALSE)
  
  if(verbose == 3)
    message("Predicting sample classes using feature pair inequalities.")
  
  predictions <- do.call(rbind, lapply(testDataFrame, function(sampleValues)
  {
    names(sampleValues) <- rownames(testDataFrame)
    measureDifferences <- sampleValues[S4Vectors::second(featurePairs)] - sampleValues[S4Vectors::first(featurePairs)]
    useFeatures <- which(abs(measureDifferences) > minDifference)
    if(length(useFeatures) == 0) # No features have a large enough distance difference.
    {                            # Simply vote for the larger class.
      if(largerClass == levels(classesTrain)[1])
      {
        class <- levels(classesTrain)[1]
        score <- -1
      } else {
        class <- levels(classesTrain)[2]
        score <- 1
      }
    } else { # One or more features are available to vote with.
      measureDifferences <- measureDifferences[useFeatures]
      if(difference == "unweighted")
      {
        # For being in second class.
        score <- sum(measureDifferences > 0)
            
        if(score > length(measureDifferences) / 2)
          class <- levels(classesTrain)[2]
        else
          class <- levels(classesTrain)[1]
            
      } else { # Each pair contributes a score for class prediction.
               # For being in second class.
        score <- sum(measureDifferences)

        # Sum of scores is tested for being positive or negative.
        class <- levels(classesTrain)[(sum(measureDifferences) > 0) + 1]
      }
    }
    data.frame(class = factor(class, levels = levels(classesTrain)), score, check.names = FALSE)
  }))
  colnames(predictions)[2] <- levels(classesTrain)[2]

  switch(returnType, class = predictions[, "class"],
         score = predictions[, 2],
         both = predictions)
}
attr(kTSPclassifier, "name") <- "kTSPclassifier"