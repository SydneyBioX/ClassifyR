# Classification based on Differential Distribution utilising Mixtures of Normals
mixModelsTrain <- function(measurementsTrain, classesTrain, ..., verbose = 3) # Mixed data types.
{
  if(verbose == 3)
    message(Sys.time(), ":Fitting mixtures of normals for features.")
  if(!requireNamespace("Rmixmod", quietly = TRUE))
    stop("The package 'Rmixmod' could not be found. Please install it.")

  models <- lapply(levels(classesTrain), function(class)
            {
                aClassMeasurements <- measurementsTrain[classesTrain == class, , drop = FALSE]
                apply(aClassMeasurements, 2, function(featureColumn)
                {
                   mixmodParams <- list(featureColumn)
                   mixmodParams <- append(mixmodParams, list(...))
                   do.call(Rmixmod::mixmodCluster, mixmodParams)
                })
            })

  if(verbose == 3)
    message(Sys.time(), ": Done fitting normal mixtures.")

  models <- lapply(models, function(modelSet)
            {
              class(modelSet) <- "MixModelsList"
              modelSet
            })
  names(models)[1:length(levels(classesTrain))] <- paste(levels(classesTrain), "Models", sep = '')
  models[["classSizes"]] <- setNames(as.vector(table(classesTrain)), levels(classesTrain))
  models <- MixModelsListsSet(models)
  models
}

# models is of class MixModelsListsSet
mixModelsPredict <- function(models, measurementsTest, difference = c("unweighted", "weighted"),
                             weighting = c("height difference", "crossover distance"),
                             densityXvalues = 1024, minDifference = 0,
                             returnType = c("both", "class", "score"), verbose = 3)
{
  models <- models@set
  measurementsTest <- measurementsTest[, names(models[[1]]), drop = FALSE]

  difference <- match.arg(difference)
  weighting <- match.arg(weighting)
  returnType <- match.arg(returnType)
  classesNames <- names(models[["classSizes"]])
  classesSizes <- models[["classSizes"]]
  largestClass <- names(classesSizes)[which.max(classesSizes)[1]]
  models <- models[-length(models)]

  if(verbose == 3)
    message("Predicting using normal mixtures.")

  featuresDensities <- lapply(1:ncol(measurementsTest), function(featureIndex)
  {
    featureValues <- unlist(lapply(models, function(classModels) classModels[[featureIndex]]@data))
    xValues <- seq(min(featureValues), max(featureValues), length.out = densityXvalues)
    setNames(lapply(models, function(model)
    {
      yValues <- Reduce('+', lapply(1:model[[featureIndex]]@bestResult@nbCluster, function(index)
      {
        model[[featureIndex]]@bestResult@parameters@proportions[index] * dnorm(xValues, model[[featureIndex]]@bestResult@parameters@mean[index], sqrt(as.numeric(model[[featureIndex]]@bestResult@parameters@variance[[index]])))
      }))
      list(x = xValues, y = yValues)
    }), classesNames)
  })

  splines <- lapply(featuresDensities, function(featureDensities)
             {
               lapply(featureDensities, function(classDensities)
               {
                 splinefun(classDensities[['x']], classDensities[['y']], "natural")
               })
             })

  if(verbose == 3)
    message(Sys.time(), ": Calculating vertical differences between normal mixture densities.")

  # Needed even if horizontal distance weighting is used to determine the predicted class.
  posteriorsVertical <- mapply(function(featureSplines, testSamples)
  {
    vertical <- sapply(1:length(classesNames), function(classIndex)
    {
      featureSplines[[classIndex]](testSamples)
    })
    if(!is.matrix(vertical)) vertical <- matrix(vertical, nrow = 1)
      vertical
  }, splines, measurementsTest, SIMPLIFY = FALSE)

  classesVertical <- sapply(posteriorsVertical, function(featureVertical)
  {
      apply(featureVertical, 1, function(sampleVertical) classesNames[which.max(sampleVertical)])
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

  if(difference == "crossover distance")
  {
    if(verbose == 3)
      message(Sys.time(), ": Calculating horizontal distances to crossover points of class densities.")

    classesVerticalIndices <- matrix(match(classesVertical, classesNames),
                                     nrow = nrow(classesVertical), ncol = ncol(classesVertical))
    distancesHorizontal <- mapply(function(featureDensities, testSamples, predictedClasses)
    {
      classesCrosses <- .densitiesCrossover(featureDensities)
      classesDistances <- sapply(classesCrosses, function(classCrosses)
      {
        sapply(testSamples, function(testSample) min(abs(testSample - classCrosses)))
      })
      classesDistances[cbind(1:nrow(classesDistances), predictedClasses)]
    }, featuresDensities, measurementsTest, as.data.frame(classesVerticalIndices)) # Matrix of horizontal distances to nearest cross-over involving the predicted class.
  }

  if(verbose == 3)
  {
    switch(returnType, class = message("Determining class labels."),
                       both = message("Calculating class scores and determining class labels."),
                       score = message("Calculating class scores.")
    )
  }

  allDistances <- switch(weighting, `height difference` = distancesVertical,
                         `crossover distance` = distancesHorizontal)

  predictions <- do.call(rbind, lapply(1:nrow(allDistances), function(sampleRow)
  {
    useFeatures <- abs(allDistances[sampleRow, ]) > minDifference
    if(all(useFeatures == FALSE)) # No features have a large enough density difference.
    {                          # Simply vote for the larger class.
      classPredicted <- largestClass
      classScores <- classesSizes / sum(classesSizes)
    } else { # One or more features are available to vote with.
      distancesUsed <- allDistances[sampleRow, useFeatures]
      classPredictionsUsed <- factor(classesVertical[sampleRow, useFeatures], classesNames)
      if(difference == "unweighted")
      {
        classScores <- table(classPredictionsUsed)
        classScores <- setNames(as.vector(classScores), classesNames)
      } else { # Weighted voting.
        classScores <- tapply(distancesUsed, classPredictionsUsed, sum)
        classScores[is.na(classScores)] <- 0
      }
      classScores <- classScores / sum(classScores) # Make different feature selection sizes comparable.
      classPredicted <- names(classScores)[which.max(classScores)]
    }
    data.frame(class = factor(classPredicted, levels = classesNames), t(classScores), check.names = FALSE)
  }))

  classPredictions <- predictions[, "class"]
  classScores <- predictions[, colnames(predictions) %in% classesNames]
  rownames(classScores) <- names(classPredictions) <- rownames(measurementsTest)
  switch(returnType, class = classPredictions,
         score = classScores,
         both = data.frame(class = predictions[, "class"], predictions[, colnames(predictions) %in% classesNames, drop = FALSE], check.names = FALSE)
  )
}
attr(mixModelsTrain, "name") <- "mixModelsTrain"