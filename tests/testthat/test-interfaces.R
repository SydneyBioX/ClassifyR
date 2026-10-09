# Regression tests for the classifier interfaces, the feature rankings and their helpers.

asDataFrame <- function(measurements) S4Vectors::DataFrame(measurements, check.names = FALSE)

test_that("DLDA favours the larger class when the features carry no information", {
  set.seed(1)
  classes <- factor(rep(c("A", "B"), c(90, 10)))
  train <- asDataFrame(matrix(rnorm(100 * 5), 100, 5, dimnames = list(NULL, paste0("g", 1:5))))
  test <- asDataFrame(matrix(rnorm(100 * 5), 100, 5, dimnames = list(NULL, paste0("g", 1:5))))
  model <- ClassifyR:::DLDAtrainInterface(train, classes, verbose = 0)
  predicted <- ClassifyR:::DLDApredictInterface(model, test, verbose = 0)
  expect_gt(mean(predicted[, "class"] == "A"), 0.5)
  # The predicted class is the one with the largest posterior probability.
  expect_equal(as.character(predicted[, "class"]), c("A", "B")[apply(predicted[, c("A", "B")], 1, which.max)])
})

test_that("DLDA posterior probabilities do not underflow with many features", {
  set.seed(2)
  classes <- factor(rep(c("A", "B"), each = 20))
  train <- matrix(rnorm(40 * 2000), 40, 2000, dimnames = list(NULL, paste0("g", 1:2000)))
  train[classes == "B", 1:20] <- train[classes == "B", 1:20] + 1
  model <- ClassifyR:::DLDAtrainInterface(asDataFrame(train), classes, verbose = 0)
  scores <- ClassifyR:::DLDApredictInterface(model, asDataFrame(train), returnType = "score", verbose = 0)
  expect_false(anyNA(scores))
  expect_equal(unname(rowSums(scores)), rep(1, 40))
})

test_that("DLDA ignores non-numeric features in training and orders posterior columns by class", {
  set.seed(3)
  classes <- factor(rep(c("B", "A"), each = 15), levels = c("B", "A"))
  train <- asDataFrame(matrix(rnorm(30 * 3), 30, 3, dimnames = list(NULL, paste0("g", 1:3))))
  train[["sex"]] <- factor(rep(c("F", "M"), 15))
  model <- ClassifyR:::DLDAtrainInterface(train, classes, verbose = 0)
  both <- ClassifyR:::DLDApredictInterface(model, train, verbose = 0)
  score <- ClassifyR:::DLDApredictInterface(model, train, returnType = "score", verbose = 0)
  expect_equal(colnames(both), c("class", "B", "A"))
  expect_equal(colnames(score), c("B", "A"))
  expect_equal(as.matrix(both[, -1]), score, ignore_attr = TRUE)
})

test_that("GLM is fitted with an intercept", {
  set.seed(4)
  classes <- factor(rep(c("A", "B"), each = 50))
  train <- asDataFrame(matrix(rnorm(100), ncol = 1, dimnames = list(NULL, "g1")))
  train[["g1"]] <- train[["g1"]] + 10 + 2 * (classes == "B")
  model <- ClassifyR:::GLMtrainInterface(train, classes, verbose = 0)
  expect_true("(Intercept)" %in% names(coef(model)))
  predicted <- ClassifyR:::GLMpredictInterface(model, train, returnType = "class", verbose = 0)
  expect_gt(mean(predicted == classes), 0.75)
})

test_that("CoxNet feature ranking uses the size of protective coefficients", {
  set.seed(5)
  measurements <- matrix(rnorm(150 * 10), 150, 10, dimnames = list(NULL, paste0("g", 1:10)))
  time <- rexp(150, rate = exp(-1.5 * measurements[, 1] + 0.3 * measurements[, 2]))
  outcome <- survival::Surv(time, rep(1, 150))
  model <- suppressWarnings(ClassifyR:::coxnetTrainInterface(asDataFrame(measurements), outcome, verbose = 0))
  ranked <- ClassifyR:::penalisedFeatures(model)[[1]]
  expect_equal(ranked[1], 1)
})

# Training data with a categorical feature, and a test set in which one of its levels is absent.
makeCategorical <- function(seed = 6)
{
  set.seed(seed)
  classes <- factor(rep(c("A", "B"), each = 40))
  measurements <- data.frame(`g-1` = rnorm(80) + (classes == "B"), g2 = rnorm(80),
                             grade = factor(sample(c("low", "mid", "high"), 80, replace = TRUE), levels = c("low", "mid", "high")),
                             site = sample(c("x", "y"), 80, replace = TRUE), check.names = FALSE)
  rownames(measurements) <- paste0("s", 1:80)
  testSamples <- which(measurements[["grade"]] != "mid")[1:10]
  test <- measurements[testSamples, ]
  test[["grade"]] <- droplevels(test[["grade"]])
  list(train = asDataFrame(measurements), classes = classes, test = asDataFrame(test), testSamples = testSamples)
}

test_that("penalised GLM encodes test data with the training columns", {
  data <- makeCategorical()
  model <- ClassifyR:::penalisedGLMtrainInterface(data$train, data$classes, verbose = 0)
  allScores <- ClassifyR:::penalisedGLMpredictInterface(model, data$train, returnType = "score", verbose = 0)
  testScores <- ClassifyR:::penalisedGLMpredictInterface(model, data$test, returnType = "score", verbose = 0)
  expect_equal(testScores, allScores[data$testSamples, ], ignore_attr = TRUE)
  # Columns of the test data in a different order give the same predictions.
  reordered <- ClassifyR:::penalisedGLMpredictInterface(model, data$test[, 4:1], returnType = "score", verbose = 0)
  expect_equal(reordered, testScores)
})

test_that("CoxNet encodes test data with the training columns", {
  data <- makeCategorical()
  set.seed(7)
  outcome <- survival::Surv(rexp(80, exp(as.numeric(data$classes))), rep(1, 80))
  model <- suppressWarnings(ClassifyR:::coxnetTrainInterface(data$train, outcome, verbose = 0))
  allRisks <- ClassifyR:::coxnetPredictInterface(model, data$train, verbose = 0)
  testRisks <- ClassifyR:::coxnetPredictInterface(model, data$test, verbose = 0)
  expect_equal(unname(testRisks), unname(allRisks[data$testSamples]))
})

test_that("XGB trains and predicts with xgboost 3 and encodes test data with the training columns", {
  data <- makeCategorical()
  set.seed(8)
  model <- ClassifyR:::extremeGradientBoostingTrainInterface(data$train, data$classes, verbose = 0)
  allScores <- ClassifyR:::extremeGradientBoostingPredictInterface(model, data$train, returnType = "score", verbose = 0)
  testPredictions <- ClassifyR:::extremeGradientBoostingPredictInterface(model, data$test, verbose = 0)
  expect_equal(colnames(testPredictions), c("class", "A", "B"))
  expect_equal(as.matrix(testPredictions[, c("A", "B")]), allScores[data$testSamples, ], ignore_attr = TRUE)
  expect_s3_class(testPredictions[, "class"], "factor")
  expect_gt(length(ClassifyR:::XGBfeatures(model)[[2]]), 0)

  # A user-specified number of threads replaces the default of one.
  expect_no_error(ClassifyR:::extremeGradientBoostingTrainInterface(data$train, data$classes, nthread = 2, verbose = 0))

  outcome <- survival::Surv(rexp(80, exp(as.numeric(data$classes))), rep(1, 80))
  survivalModel <- ClassifyR:::extremeGradientBoostingTrainInterface(data$train, outcome, verbose = 0)
  risks <- ClassifyR:::extremeGradientBoostingPredictInterface(survivalModel, data$test, verbose = 0)
  expect_true(is.numeric(risks))
  expect_length(risks, 10)
})

test_that("SVM predicts with categorical features and non-syntactic feature names", {
  data <- makeCategorical()
  model <- ClassifyR:::SVMtrainInterface(data$train, data$classes, verbose = 0)
  predictions <- ClassifyR:::SVMpredictInterface(model, data$test, verbose = 0)
  expect_equal(nrow(predictions), 10)
  expect_equal(colnames(predictions), c("class", "A", "B"))
})

test_that("kNN returns a factor with the training levels, also for one test sample", {
  data <- makeTwoClass(shift = 5)
  train <- asDataFrame(data$measurements[1:50, ])
  test <- asDataFrame(data$measurements[51:60, ])
  for(mode in c("unweighted", "weighted"))
  {
    predicted <- ClassifyR:::kNNinterface(train, data$classes[1:50], test, k = 3, mode = mode, verbose = 0)
    expect_s3_class(predicted[, "class"], "factor")
    expect_equal(levels(predicted[, "class"]), c("A", "B"))
    expect_equal(as.character(predicted[, "class"]), c("A", "B")[apply(predicted[, c("A", "B")], 1, which.max)])
    expect_gt(mean(predicted[, "class"] == data$classes[51:60]), 0.5)
    one <- ClassifyR:::kNNinterface(train, data$classes[1:50], test[1, , drop = FALSE], k = 3, mode = mode, verbose = 0)
    expect_equal(nrow(one), 1)
  }
})

test_that("Fisher discriminant classifies two classes and rejects more", {
  data <- makeTwoClass(shift = 3)
  train <- asDataFrame(data$measurements[1:40, 1:5])
  test <- asDataFrame(data$measurements[41:60, 1:5])
  predicted <- ClassifyR:::fisherDiscriminant(train, data$classes[1:40], test, verbose = 0)
  expect_equal(nrow(predicted), 20)
  expect_gt(mean(predicted[, "class"] == data$classes[41:60]), 0.8)
  # The pooled variance is weighted by the class sizes.
  trainMatrix <- as.matrix(train)
  isA <- data$classes[1:40] == "A"
  pooled <- (19 * apply(trainMatrix[isA, ], 2, var) + 19 * apply(trainMatrix[!isA, ], 2, var)) / 38
  direction <- (colMeans(trainMatrix[isA, ]) - colMeans(trainMatrix[!isA, ])) / pooled
  expect_equal(unname(predicted[, "score"]), unname(-1 * as.matrix(test) %*% direction)[, 1])
  threeClasses <- factor(rep(c("A", "B", "C"), length.out = 40))
  expect_error(ClassifyR:::fisherDiscriminant(train, threeClasses, test, verbose = 0), "two classes")
})

test_that("k-TSP classifier predicts with feature pairs", {
  data <- makeTwoClass(shift = 3)
  # The pairs reverse their order between the classes.
  data$measurements[data$classes == "A", 1:3] <- data$measurements[data$classes == "A", 1:3] - 3
  train <- asDataFrame(data$measurements[1:40, ])
  test <- asDataFrame(data$measurements[41:60, ])
  pairs <- S4Vectors::Pairs(c("g1", "g2", "g3"), c("g10", "g11", "g12"))
  for(difference in c("unweighted", "weighted"))
  {
    predicted <- ClassifyR:::kTSPclassifier(train, data$classes[1:40], test, featurePairs = pairs,
                                            difference = difference, verbose = 0)
    expect_equal(nrow(predicted), 20)
    expect_gt(mean(predicted[, "class"] == data$classes[41:60]), 0.8)
  }
  threeClasses <- factor(rep(c("A", "B", "C"), length.out = 40))
  expect_error(ClassifyR:::kTSPclassifier(train, threeClasses, test, featurePairs = pairs, verbose = 0), "two classes")
})

test_that("Poisson LDA classifies counts", {
  skip_if_not_installed("PoiClaClu")
  set.seed(9)
  classes <- factor(rep(c("A", "B"), each = 20))
  counts <- matrix(rpois(40 * 20, 20), 40, 20, dimnames = list(paste0("s", 1:40), paste0("g", 1:20)))
  counts[classes == "B", 1:5] <- rpois(20 * 5, 60)
  train <- asDataFrame(counts[c(1:15, 21:35), ])
  test <- asDataFrame(counts[c(16:20, 36:40), ])
  predicted <- ClassifyR:::classifyInterface(train, classes[c(1:15, 21:35)], test, verbose = 0)
  expect_equal(levels(predicted[, "class"]), c("A", "B"))
  expect_equal(as.character(predicted[, "class"]), rep(c("A", "B"), each = 5))
})

test_that("naive Bayes and mixtures of normals run with crossover distance weighting", {
  data <- makeTwoClass(shift = 3)
  train <- asDataFrame(data$measurements[1:40, 1:4])
  test <- asDataFrame(data$measurements[41:60, 1:4])
  for(difference in c("unweighted", "weighted"))
  {
    predicted <- ClassifyR:::naiveBayesKernel(train, data$classes[1:40], test, difference = difference,
                                              weighting = "crossover distance", verbose = 0)
    expect_equal(nrow(predicted), 20)
    expect_gt(mean(predicted[, "class"] == data$classes[41:60]), 0.7)
  }
  skip_if_not_installed("Rmixmod")
  set.seed(10)
  models <- ClassifyR:::mixModelsTrain(train, data$classes[1:40], nbCluster = 1, verbose = 0)
  for(difference in c("unweighted", "weighted"))
  {
    predicted <- suppressWarnings(ClassifyR:::mixModelsPredict(models, test, difference = difference,
                                                               weighting = "crossover distance", verbose = 0))
    expect_equal(nrow(predicted), 20)
    expect_gt(mean(predicted[, "class"] == data$classes[41:60]), 0.7)
  }
})

test_that("colCoxTests handles one feature, and the slow option agrees with the fast one", {
  data <- makeSurvival()
  one <- colCoxTests(data$measurements[, 1, drop = FALSE], data$outcome)
  expect_equal(dim(one), c(1, 3))
  expect_equal(rownames(one), "g1")
  fast <- colCoxTests(data$measurements[, 1:5], data$outcome, "fast")
  slow <- colCoxTests(data$measurements[, 1:5], data$outcome, "slow")
  expect_equal(rownames(slow), paste0("g", 1:5))
  expect_equal(slow, fast, tolerance = 1e-4)
  # A two-column matrix of time and event is also accepted.
  expect_equal(colCoxTests(data$measurements[, 1:5], as.matrix(data$outcome)[, 1:2], "slow"), slow)
  # CoxPH ranking with one feature.
  expect_equal(ClassifyR:::coxphRanking(asDataFrame(data$measurements[, 1, drop = FALSE]), data$outcome, verbose = 0), 1)
})

test_that("subtractFromLocation works for one numeric feature and keeps feature names", {
  train <- asDataFrame(matrix(c(1, 2, 6), ncol = 1, dimnames = list(NULL, "g-1")))
  train[["sex"]] <- factor(c("F", "M", "F"))
  test <- asDataFrame(matrix(c(0, 10), ncol = 1, dimnames = list(NULL, "g-1")))
  test[["sex"]] <- factor(c("F", "M"))
  transformed <- ClassifyR:::subtractFromLocation(train, test, verbose = 0)
  expect_equal(colnames(transformed[[1]]), "g-1")
  expect_equal(transformed[[1]][["g-1"]], c(2, 1, 3))
  expect_equal(transformed[[2]][["g-1"]], c(3, 7))
  medians <- ClassifyR:::subtractFromLocation(train, test, location = "median", absolute = FALSE, verbose = 0)
  expect_equal(medians[[2]][["g-1"]], c(-2, 8))
})

# A clinical table and an RNA table, as prevalidation expects them.
makeMultiView <- function(nClasses, clinicalSelection = "none", seed = 11)
{
  set.seed(seed)
  classes <- factor(rep(LETTERS[seq_len(nClasses)], length.out = 60))
  clinical <- S4Vectors::DataFrame(age = rnorm(60) + as.numeric(classes), bmi = rnorm(60), row.names = paste0("s", 1:60))
  rna <- asDataFrame(matrix(rnorm(60 * 20), 60, 20, dimnames = list(paste0("s", 1:60), paste0("g", 1:20))))
  rna[, 1] <- rna[, 1] + 2 * as.numeric(classes)
  measurements <- cbind(clinical, rna)
  S4Vectors::mcols(measurements) <- S4Vectors::DataFrame(assay = rep(c("clinical", "rna"), c(2, 20)), feature = colnames(measurements))
  params <- list(clinical = ClassifyR:::generateModellingParams("clinical", clinical, 2, clinicalSelection, "DLDA", extraParams = NULL),
                 rna = ClassifyR:::generateModellingParams("rna", rna, 3, "t-test", "DLDA", extraParams = NULL))
  list(measurements = measurements, classes = classes, params = params)
}

test_that("prevalidation works for two and three classes and fits only the needed models", {
  for(nClasses in 2:3)
  {
    data <- makeMultiView(nClasses)
    model <- ClassifyR:::prevalTrainInterface(data$measurements, data$classes, data$params, verbose = 0)
    expect_equal(names(model@fullModel$prevalidationModels), "rna")
    prevalidationColumns <- if(nClasses == 2) "rna" else c("rna_B", "rna_C")
    expect_equal(model@fullModel$fullFeatures, c("age", "bmi", prevalidationColumns))
    predicted <- ClassifyR:::prevalPredictInterface(model, data$measurements, verbose = 0)
    expect_equal(nrow(predicted), 60)
    expect_gt(mean(predicted[, "class"] == data$classes), 1 / nClasses)
  }
  # A clinical table with one feature.
  data <- makeMultiView(2)
  data$measurements <- data$measurements[, -2]
  model <- ClassifyR:::prevalTrainInterface(data$measurements, data$classes, data$params, verbose = 0)
  expect_equal(model@fullModel$fullFeatures, c("age", "rna"))
})

test_that("prevalidation and PCA draw their inner cross-validation seed from the random number stream", {
  data <- makeMultiView(2, clinicalSelection = "t-test")
  seeds <- integer()
  local_mocked_bindings(SerialParam = function(...)
  {
    seeds <<- c(seeds, list(...)[["RNGseed"]])
    BiocParallel::SerialParam(...)
  }, .package = "ClassifyR")
  set.seed(2)
  ClassifyR:::prevalTrainInterface(data$measurements, data$classes, data$params, verbose = 0)
  ClassifyR:::prevalTrainInterface(data$measurements, data$classes, data$params, verbose = 0)
  ClassifyR:::pcaTrainInterface(data$measurements, data$classes, data$params["clinical"], nFeatures = c(rna = 2))
  ClassifyR:::pcaTrainInterface(data$measurements, data$classes, data$params["clinical"], nFeatures = c(rna = 2))
  expect_length(seeds, 4)
  expect_length(unique(seeds), 4)
})

test_that("limma ranking passes extra arguments to lmFit", {
  data <- makeTwoClass()
  train <- asDataFrame(data$measurements)
  unweighted <- ClassifyR:::limmaRanking(train, data$classes, verbose = 0)
  weighted <- ClassifyR:::limmaRanking(train, data$classes, weights = rep(1, nrow(train)), verbose = 0)
  expect_equal(weighted, unweighted)
})

test_that("edgesToHubNetworks accepts a matrix", {
  edges <- cbind(c("MITF", "MITF", "MITF", "KRAS"), c("HINT1", "LEF1", "PSMD14", "ARAF"))
  hubs <- edgesToHubNetworks(edges, minCardinality = 3)
  expect_equal(names(hubs@sets), "MITF")
})

test_that("previousSelection warns when few previous features are in the current data", {
  data <- makeTwoClass()
  set.seed(12)
  result <- crossValidate(data$measurements, data$classes, classifier = "DLDA", nFeatures = 5,
                          nFolds = 2, nRepeats = 1, verbose = 0)
  previous <- chosenFeatureNames(result)[[1]]
  current <- data$measurements
  colnames(current)[match(previous[1:3], colnames(current))] <- paste0("new", 1:3)
  expect_warning(selected <- ClassifyR:::previousSelection(asDataFrame(current), data$classes, result, .iteration = 1, verbose = 0),
                 "40% of the previously selected features")
  expect_equal(colnames(current)[selected], previous[4:5])
  expect_no_warning(ClassifyR:::previousSelection(asDataFrame(data$measurements), data$classes, result, .iteration = 1, verbose = 0))
})
