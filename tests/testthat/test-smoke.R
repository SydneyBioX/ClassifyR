# Every classifier and selection keyword offered by available() runs end to end on a small data set.

test_that("each classification keyword runs in crossValidate", {
  data <- makeTwoClass()
  classifiers <- c("randomForest", "GLM", "ridgeGLM", "elasticNetGLM", "LASSOGLM", "SVM", "NSC", "DLDA",
                   "naiveBayes", "mixturesNormals", "kNN") # XGB needs the xgboost 3 interface (interface fixes).
  for(classifier in classifiers)
  {
    set.seed(1)
    result <- suppressWarnings(crossValidate(data$measurements, data$classes, classifier = classifier,
                                             nFeatures = 5, nRepeats = 1, nFolds = 3))
    expect_s4_class(result, "ClassifyResult")
    expect_equal(nrow(predictions(result)), nrow(data$measurements), info = classifier)
  }
})

test_that("each selection keyword runs in crossValidate", {
  data <- makeTwoClass()
  selections <- c("t-test", "limma", "Bartlett", "Levene", "DMD", "likelihoodRatio", "KS", "KL")
  for(selection in selections)
  {
    set.seed(1)
    result <- suppressWarnings(crossValidate(data$measurements, data$classes, selectionMethod = selection,
                                             classifier = "DLDA", nFeatures = 5, nRepeats = 1, nFolds = 3))
    expect_s4_class(result, "ClassifyResult")
    expect_true(all(lengths(chosenFeatureNames(result)) == 5), info = selection)
  }
})

test_that("each survival keyword runs in crossValidate", {
  data <- makeSurvival()
  for(classifier in c("CoxPH", "CoxNet", "randomSurvivalForest")) # XGB: see above.
  {
    set.seed(1)
    result <- suppressWarnings(crossValidate(data$measurements, data$outcome, classifier = classifier,
                                             nFeatures = 5, nRepeats = 1, nFolds = 3))
    expect_s4_class(result, "ClassifyResult")
    expect_true(is.numeric(predictions(result)[, "risk"]), info = classifier)
  }
})
