# Classification Using A Bayes Classifier with Kernel Density Estimates
naiveBayesKernel <- function(measurementsTrain, classesTrain, measurementsTest,
                             densityFunction = density, densityParameters = list(bw = "nrd0", n = 1024, from = expression(min(featureValues)), to = expression(max(featureValues))),
                             difference = c("unweighted", "weighted"),
                             weighting = c("height difference", "crossover distance"),
                             minDifference = 0, returnType = c("both", "class", "score"), verbose = 3)
{
  trainingMatrix <- as.matrix(measurementsTrain)
  testingMatrix <- as.matrix(measurementsTest[, colnames(trainingMatrix), drop = FALSE])
  
  difference <- match.arg(difference)
  weighting <- match.arg(weighting)
  returnType <- match.arg(returnType)
  
  classesSizes <- sapply(levels(classesTrain), function(class) sum(classesTrain == class))
  largestClass <- names(classesSizes)[which.max(classesSizes)[1]]
  
  if(verbose == 3)
    message(Sys.time(), ": Fitting densities.")
  
  featuresDensities <- lapply(measurementsTrain, function(featureValues)
  {
    densityParameters <- lapply(densityParameters, function(parameter) eval(parameter))
    lapply(levels(classesTrain), function(class)
    {
      aClassMeasurements <- featureValues[classesTrain == class]  
      do.call(densityFunction, c(list(aClassMeasurements), densityParameters))
    }) # A fitted density for each class.
  })

  classesScaleFactors <- classesSizes / nrow(trainingMatrix)
  splines <- lapply(featuresDensities, function(featureDensities) 
             {
               mapply(function(featureDensity, scaleFactor)
               {
                 splinefun(featureDensity[['x']], featureDensity[['y']] * scaleFactor, "natural")
               }, featureDensities, classesScaleFactors)
             })
  
  if(verbose == 3)
    message(Sys.time(), ": Calculating vertical distances between class densities.")

  # Needed even if horizontal distance weighting is used to determine the predicted class.
  posteriorsVertical <- mapply(function(featureSplines, testSamples)
  {
    vertical <- sapply(1:length(levels(classesTrain)), function(classIndex)
    {
      featureSplines[[classIndex]](testSamples)
    })
    if(!is.matrix(vertical)) vertical <- matrix(vertical, nrow = 1)
    vertical
  }, splines, as.data.frame(testingMatrix), SIMPLIFY = FALSE)

  classesVertical <- sapply(posteriorsVertical, function(featureVertical)
  {
      apply(featureVertical, 1, function(sampleVertical) levels(classesTrain)[which.max(sampleVertical)])
  }) # Matrix, rows are test samples, columns are features.
  if(!is.matrix(classesVertical)) classesVertical <- matrix(classesVertical, nrow = 1)
    
  distancesVertical <- sapply(posteriorsVertical, function(featureVertical)
  { # Vertical distance between highest density and second-highest, at a particular value.
    apply(featureVertical, 1, function(sampleVertical)
    {
      twoHighest <- sort(sampleVertical, decreasing = TRUE)[1:2]
      Reduce('-', twoHighest)
    })
  }) # Matrix, rows are test samples, columns are features.
  if(!is.matrix(distancesVertical)) distancesVertical <- matrix(distancesVertical, nrow = 1)
  
  if(difference == "weighted" && weighting == "crossover distance")
  {
    if(verbose == 3)
      message(Sys.time(), ": Calculating horizontal distances to crossover points of class densities.")
 
    classesVerticalIndices <- matrix(match(classesVertical, levels(classesTrain)),
                                     nrow = nrow(classesVertical), ncol = ncol(classesVertical))
    distancesHorizontal <- mapply(function(featureDensities, testSamples, predictedClasses)
    {
      classesCrosses <- .densitiesCrossover(featureDensities)
      classesDistances <- sapply(classesCrosses, function(classCrosses)
      {
        sapply(testSamples, function(testSample) min(abs(testSample - classCrosses)))
      })
      classesDistances[cbind(1:nrow(classesDistances), predictedClasses)]
    }, featuresDensities, test, as.data.frame(classesVerticalIndices)) # Matrix of horizontal distances to nearest cross-over involving the predicted class.
  }

  if(verbose == 3)
  {
    switch(returnType, class = message("Determining class labels."),
                       both = message("Calculating class scores and determining class labels."),
                       score = message("Calculating class scores."))
  }
  
  allDistances <- switch(weighting, `height difference` = distancesVertical,
                                    `crossover distance` = distancesHorizontal)

  predictions <- do.call(rbind, lapply(1:nrow(allDistances), function(sampleRow)
  {
    useFeatures <- abs(allDistances[sampleRow, ]) > minDifference
    if(all(useFeatures == FALSE)) # No features have a large enough density difference.
    {                          # Simply vote for the larger class.
      classPredicted <- largestClass
      classScores <- classesSizes / length(classesTrain)
    } else { # One or more features are available to vote with.
      distancesUsed <- allDistances[sampleRow, useFeatures]
      classPredictionsUsed <- factor(classesVertical[sampleRow, useFeatures], levels(classesTrain))
      if(difference == "unweighted")
      {
        classScores <- table(classPredictionsUsed)
        classScores <- setNames(as.vector(classScores), levels(classesTrain))
      } else { # Weighted voting.
        classScores <- tapply(distancesUsed, classPredictionsUsed, sum)
        classScores[is.na(classScores)] <- 0
      }
      classScores <- classScores / sum(classScores) # Make different feature selection sizes comparable.
      classPredicted <- names(classScores)[which.max(classScores)]
    }

    data.frame(class = factor(classPredicted, levels = levels(classesTrain)), t(classScores), check.names = FALSE)
  }))
  
  classPredictions <- predictions[, "class"]
  classScores <- predictions[, 2:ncol(predictions)]
  rownames(classScores) <- names(classPredictions) <- rownames(measurementsTest)
  switch(returnType, class = classPredictions,
         score = classScores,
         both = predictions)
}
attr(naiveBayesKernel, "name") <- "naiveBayesKernel"