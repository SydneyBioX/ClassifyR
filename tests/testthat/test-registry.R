# Every registered classifier and feature selection method, fitted and predicted through crossValidate.

classifierPackages <- list(randomForest = "ranger", kNN = "class", ridgeGLM = "glmnet", elasticNetGLM = "glmnet",
                           LASSOGLM = "glmnet", SVM = "e1071", NSC = "pamr", mixturesNormals = "Rmixmod",
                           CoxNet = "glmnet", randomSurvivalForest = "randomForestSRC", XGB = "xgboost",
                           aorsf = "aorsf", glmboost = "mboost", LiblineaR = "LiblineaR", ncvreg = "ncvreg")
selectionPackages <- list(limma = "limma", edgeR = "edgeR", Levene = "car", DMD = "robustbase")

hasPackages <- function(packages) all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))

test_that("available() lists the registered keywords with their outcomes", {
  classifiers <- available()
  expect_identical(colnames(classifiers), c("classifier Keyword", "Description", "Outcomes"))
  expect_false("previousTrained" %in% classifiers[, "classifier Keyword"])
  expect_identical(classifiers[classifiers[, 1] == "CoxPH", "Outcomes"], "survival")
  selections <- available("selectionMethod")
  expect_identical(colnames(selections), c("selectionMethod Keyword", "Description", "Outcomes"))
  expect_true(all(c("none", "t-test", "CoxPH") %in% selections[, 1]))
})

test_that("every classifier keyword fits and predicts its kinds of outcome in the expected format", {
  classesData <- makeTwoClass(nSamples = 40, nFeatures = 10, shift = 2)
  survivalData <- makeSurvival(nSamples = 60, nFeatures = 10)
  for(keyword in available()[, "classifier Keyword"])
  {
    if(!hasPackages(classifierPackages[[keyword]])) next
    outcomes <- .subset2(ClassifyR:::.ClassifyRenvir[["classifiers"]], keyword)[["outcomes"]]
    for(outcomeType in outcomes)
    {
      if(outcomeType == "classes") data <- list(classesData$measurements, classesData$classes)
      else data <- list(survivalData$measurements, survivalData$outcome)
      selection <- if(outcomeType == "classes") "t-test" else "CoxPH"
      set.seed(1)
      result <- suppressWarnings(crossValidate(data[[1]], data[[2]], nFeatures = 5, selectionMethod = selection,
                                               classifier = keyword, nFolds = 2, nRepeats = 1))
      expect_s4_class(result, "ClassifyResult")
      predicted <- predictions(result)
      if(outcomeType == "classes")
        expect_true("class" %in% colnames(predicted) && is.factor(predicted[["class"]]), label = keyword)
      else
        expect_true("risk" %in% colnames(predicted) && is.numeric(predicted[["risk"]]), label = keyword)
      expect_equal(nrow(predicted), nrow(data[[1]]), label = keyword)
    }
  }
})

test_that("every feature selection keyword ranks features for its kinds of outcome", {
  classesData <- makeTwoClass(nSamples = 40, nFeatures = 10, shift = 2)
  survivalData <- makeSurvival(nSamples = 60, nFeatures = 10)
  for(keyword in setdiff(available("selectionMethod")[, 1], "none"))
  {
    if(!hasPackages(selectionPackages[[keyword]])) next
    outcomes <- .subset2(ClassifyR:::.ClassifyRenvir[["selections"]], keyword)[["outcomes"]]
    for(outcomeType in outcomes)
    {
      if(outcomeType == "classes") data <- list(abs(classesData$measurements), classesData$classes, "DLDA")
      else data <- list(survivalData$measurements, survivalData$outcome, "CoxPH")
      set.seed(1)
      result <- suppressWarnings(crossValidate(data[[1]], data[[2]], nFeatures = 5, selectionMethod = keyword,
                                               classifier = data[[3]], nFolds = 2, nRepeats = 1))
      expect_s4_class(result, "ClassifyResult")
      expect_true(all(lengths(chosenFeatureNames(result)) == 5), label = keyword)
    }
  }
})

test_that("a method that does not accept the kind of outcome stops early with a clear error", {
  data <- makeTwoClass()
  expect_error(crossValidate(data$measurements, data$classes, classifier = "CoxPH", nFolds = 2, nRepeats = 1),
               "CoxPH can't be used with classes")
  survivalData <- makeSurvival()
  expect_error(crossValidate(survivalData$measurements, survivalData$outcome, selectionMethod = "t-test",
                             classifier = "CoxPH", nFolds = 2, nRepeats = 1), "t-test can't be used with a survival outcome")
})

test_that("keyword parameter sets match the registry, including tuning presets", {
  forest <- ClassifyR:::.classifierKeywordToParams("randomForest", "auto")
  expect_identical(forest$trainParams@tuneParams, list(mTryProportion = c(0.10, 0.25, 0.33, 0.5), num.trees = c(1, 10, 100)))
  expect_null(ClassifyR:::.classifierKeywordToParams("DLDA", "auto")$trainParams@tuneParams)
  expect_identical(ClassifyR:::.classifierKeywordToParams("ridgeGLM", NULL)$trainParams@otherParams, list(alpha = 0))
  expect_null(ClassifyR:::.classifierKeywordToParams("kNN", NULL)$predictParams)
  expect_null(ClassifyR:::.classifierKeywordToParams("notAClassifier", NULL))
})
