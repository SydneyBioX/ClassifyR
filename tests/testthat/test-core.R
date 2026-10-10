# Cross-validation driver, train/predict and data preparation.

test_that("crossValidate uses the requested selection method for list input", {
  data <- makeTwoClass()
  measurementsList <- list(a = data$measurements[, 1:15], b = data$measurements[, 16:30])
  set.seed(1)
  results <- crossValidate(measurementsList, data$classes, selectionMethod = "KS", classifier = "DLDA",
                           nFeatures = 3, nRepeats = 1, nFolds = 3)
  selectionNames <- sapply(results, function(result)
                       result@characteristics[result@characteristics[, "characteristic"] == "Selection Name", "value"])
  expect_true(all(selectionNames == "Kolmogorov-Smirnov Statistic" | selectionNames == "KS"))
  expect_false(any(selectionNames %in% c("t-test", "Difference in Means")))
})

test_that("worker pools give the same results as lapply, with forked or socket workers", {
  pool <- ClassifyR:::.makeWorkerPool(nCores = 2, nTasks = 10)
  expect_s3_class(pool, "workerPool")
  expect_equal(pool[["workers"]], 2L)
  shared <- 1:3 # Data the task function uses from its enclosing environment.
  task <- function(x) x^2 + sum(shared)
  for(type in c(if(.Platform$OS.type != "windows") "FORK", "PSOCK"))
  {
    pool[["type"]] <- type
    expect_equal(ClassifyR:::.poolApply(1:5, task, pool), lapply(1:5, task), info = type)
  }
})

test_that("socket workers give the same cross-validation results as forked ones", {
  skip_on_os("windows")
  data <- makeTwoClass()
  run <- function() { set.seed(1); crossValidate(data$measurements, data$classes, classifier = "randomForest", nFeatures = 3, nRepeats = 2, nFolds = 3, nCores = 2) }
  forked <- run()
  original <- ClassifyR:::.makeWorkerPool
  local_mocked_bindings(.makeWorkerPool = function(nCores, nTasks) { pool <- original(nCores, nTasks); pool[["type"]] <- "PSOCK"; pool }, .package = "ClassifyR")
  socket <- run()
  expect_identical(predictions(socket), predictions(forked))
})

test_that("results don't depend on the number of cores", {
  skip_on_os("windows")
  data <- makeTwoClass()
  measurementsList <- list(a = data$measurements[, 1:15], b = data$measurements[, 16:30])
  set.seed(1)
  serial <- crossValidate(measurementsList, data$classes, classifier = "randomForest", nFeatures = 3, nRepeats = 2, nFolds = 3)
  set.seed(1)
  parallel <- crossValidate(measurementsList, data$classes, classifier = "randomForest", nFeatures = 3, nRepeats = 2, nFolds = 3, nCores = 2)
  expect_identical(lapply(serial, predictions), lapply(parallel, predictions))
})

test_that("nested cross-validation tunes classifier parameters", {
  data <- makeTwoClass(nFeatures = 10)
  set.seed(1)
  result <- crossValidate(data$measurements, data$classes, classifier = "SVM", selectionMethod = "none",
                          nRepeats = 1, nFolds = 3,
                          extraParams = list(tuneCross = list(tuneMode = "Nested CV", performanceType = "Balanced Accuracy"),
                                             train = list(cost = c(0.1, 10))))
  expect_s4_class(result, "ClassifyResult")
  expect_true(all(sapply(tunedParameters(result), function(tune) "cost" %in% colnames(tune[["tuneCombinations"]]))))
})

test_that("tuning keeps the other training settings", {
  data <- makeTwoClass(nFeatures = 10)
  trainParams <- TrainParams("SVM")
  trainParams@otherParams <- c(trainParams@otherParams, list(kernel = "linear"))
  trainParams@tuneParams <- list(cost = c(0.1, 10))
  modellingParams <- ModellingParams(balancing = "none", selectParams = NULL, trainParams = trainParams, predictParams = PredictParams("SVM"))
  set.seed(1)
  result <- runTest(data$measurements[1:40, ], data$classes[1:40], data$measurements[41:60, ], data$classes[41:60],
                    crossValParams = CrossValParams(tuneMode = "Resubstitution", performanceType = "Balanced Accuracy",
                                                    parallelParams = BiocParallel::SerialParam()),
                    modellingParams = modellingParams)
  expect_equal(models(result)[[1]][["kernel"]], 0) # e1071 codes the linear kernel as 0.
})

test_that("tuning ranges without a tuning mode are an error, not ignored", {
  data <- makeTwoClass()
  set.seed(1)
  expect_error(crossValidate(data$measurements, data$classes, classifier = "SVM", nRepeats = 1, nFolds = 3,
                             extraParams = list(train = list(cost = c(0.1, 10)))), "no tuning mode")
})

test_that("preset tuning ranges are used and not passed on as a setting", {
  data <- makeTwoClass(nFeatures = 10)
  set.seed(1)
  result <- crossValidate(data$measurements, data$classes, classifier = "randomForest", nRepeats = 1, nFolds = 3,
                          extraParams = list(tuneCross = list(tuneMode = "Resubstitution", performanceType = "Balanced Accuracy"),
                                             train = list(tuneParams = "auto")))
  expect_false("tuneParams" %in% result@characteristics[, "characteristic"])
  expect_true("mTryProportion" %in% colnames(tunedParameters(result)[[1]][["tuneCombinations"]]))
})

test_that("nFeatures = 1 works", {
  data <- makeTwoClass()
  set.seed(1)
  result <- crossValidate(data$measurements, data$classes, classifier = "DLDA", nFeatures = 1, nRepeats = 1, nFolds = 3)
  expect_true(all(lengths(chosenFeatureNames(result)) == 1))
})

test_that("percentage splits use the outcome given, also for survival", {
  data <- makeTwoClass()
  if(exists("classes", envir = globalenv())) skip("A global 'classes' object would hide the bug.")
  modellingParams <- ModellingParams(balancing = "none", selectParams = SelectParams("t-test", nFeatures = 5),
                                     trainParams = TrainParams("DLDA"), predictParams = PredictParams("DLDA"))
  set.seed(1)
  result <- runTests(data$measurements, data$classes,
                     CrossValParams("Permute Percentage Split", permutations = 2, percentTest = 25, parallelParams = BiocParallel::SerialParam()),
                     modellingParams, verbose = 0)
  expect_s4_class(result, "ClassifyResult")
  survival <- makeSurvival()
  splits <- ClassifyR:::samplesSplits("Permute Percentage Split", permutations = 2, percentTest = 25, outcome = survival$outcome)
  expect_equal(length(splits$train), 2)
  expect_true(all(lengths(splits$train) + lengths(splits$test) == nrow(survival$measurements)))
})

test_that("train fits one model per assay and predict matches features by name", {
  data <- makeTwoClass()
  measurementsList <- list(a = data$measurements[, 1:15], b = data$measurements[, 16:30])
  set.seed(1)
  models <- train(measurementsList, data$classes, classifier = "DLDA", nFeatures = 3)
  expect_s3_class(models, "listOfModels")
  expect_length(models, 2)
  predicted <- predict(models, measurementsList)
  expect_length(predicted, 2)
  # Columns in another order give the same predictions.
  shuffled <- lapply(measurementsList, function(measurements) measurements[, rev(colnames(measurements))])
  expect_equal(predict(models, shuffled), predicted)
})

test_that("train combines assays by merging", {
  data <- makeTwoClass()
  measurementsList <- list(a = data$measurements[, 1:15], b = data$measurements[, 16:30])
  set.seed(1)
  model <- train(measurementsList, data$classes, classifier = "DLDA", multiViewMethod = "merge", nFeatures = 3)
  expect_s3_class(model, "trainedByClassifyR")
  predicted <- predict(model, measurementsList)
  expect_equal(nrow(as.data.frame(predicted)), nrow(data$measurements))
})

test_that("runTest accepts MultiAssayExperiment training and test sets", {
  data <- makeTwoClass()
  makeMAE <- function(samples)
  {
    MultiAssayExperiment::MultiAssayExperiment(list(RNA = t(data$measurements[samples, ])),
                                               S4Vectors::DataFrame(class = data$classes[samples], row.names = rownames(data$measurements)[samples]))
  }
  modellingParams <- ModellingParams(balancing = "none", selectParams = SelectParams("t-test", nFeatures = 5),
                                     trainParams = TrainParams("DLDA"), predictParams = PredictParams("DLDA"))
  set.seed(1)
  result <- suppressWarnings(runTest(makeMAE(1:40), makeMAE(41:60), outcomeColumns = "class", modellingParams = modellingParams))
  expect_s4_class(result, "ClassifyResult")
  expect_equal(nrow(predictions(result)), 20)
})

test_that("prepareData keeps the most variable features and drops similar ones", {
  data <- makeTwoClass()
  prepared <- prepareData(data$measurements, data$classes, topNvariance = 10)
  expect_equal(ncol(prepared$measurements), 10)
  measurements <- data$measurements
  measurements[, "g2"] <- measurements[, "g1"] + rnorm(nrow(measurements), sd = 0.01)
  prepared <- prepareData(measurements, data$classes, maxSimilarity = 0.9)
  expect_true("g1" %in% colnames(prepared$measurements))
  expect_false("g2" %in% colnames(prepared$measurements))
  expect_equal(ncol(prepared$measurements), ncol(measurements) - 1)
})

test_that("a missing or empty outcome is an error before any modelling", {
  data <- makeTwoClass()
  outcomes <- list(cohortA = data$classes)
  # A list looked up by a name it doesn't have gives NULL.
  expect_error(prepareData(data$measurements, outcomes[["cohortB"]]), "'outcome' is missing or empty")
  expect_error(prepareData(data$measurements), "'outcome' is missing or empty")
  set.seed(1)
  expect_error(crossValidate(data$measurements, outcomes[["cohortB"]], nFolds = 2, nRepeats = 1), "'outcome' is missing or empty")
  expect_error(crossValidate(as.matrix(data$measurements), factor(character()), nFolds = 2, nRepeats = 1), "'outcome' is missing or empty")
})

test_that("all assays, classifiers and combinations share the same splits", {
  data <- makeTwoClass()
  measurementsList <- list(a = data$measurements[, 1:15], b = data$measurements[, 16:30])
  foldsOf <- function(result)
  {
    predictionsTable <- as.data.frame(predictions(result))
    predictionsTable <- predictionsTable[order(predictionsTable[, "permutation"], predictionsTable[, "sample"]), ]
    paste(predictionsTable[, "permutation"], predictionsTable[, "fold"])
  }
  set.seed(1)
  separate <- crossValidate(measurementsList, data$classes, classifier = c("DLDA", "randomForest"), nFeatures = 3, nRepeats = 2, nFolds = 3)
  folds <- lapply(separate, foldsOf)
  expect_true(all(sapply(folds, identical, folds[[1]])))
  set.seed(1)
  merged <- crossValidate(measurementsList, data$classes, classifier = "DLDA", multiViewMethod = "merge", nFeatures = 3, nRepeats = 2, nFolds = 3)
  expect_true(all(sapply(lapply(merged, foldsOf), identical, folds[[1]])))
})

test_that("an unknown multiViewMethod is an error naming the choices", {
  data <- makeTwoClass()
  set.seed(1)
  expect_error(crossValidate(list(a = data$measurements[, 1:15], b = data$measurements[, 16:30]), data$classes,
                             multiViewMethod = "Merge", nRepeats = 1, nFolds = 3), "must be one of")
})

test_that("train and predict match the samples of a list of tables by name", {
  data <- makeTwoClass()
  aligned <- list(a = data$measurements[, 1:15], b = data$measurements[, 16:30])
  shuffled <- aligned
  shuffled$b <- shuffled$b[rev(rownames(shuffled$b)), ]
  set.seed(1)
  fromAligned <- train(aligned, data$classes, classifier = "DLDA", multiViewMethod = "merge", nFeatures = 3)
  set.seed(1)
  fromShuffled <- train(shuffled, data$classes, classifier = "DLDA", multiViewMethod = "merge", nFeatures = 3)
  expect_equal(predict(fromShuffled, shuffled), predict(fromAligned, aligned))
})
