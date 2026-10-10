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

test_that("penalised GLM fits a logistic model for two classes and a multinomial model for more", {
  data <- makeTwoClass()
  model <- ClassifyR:::penalisedGLMtrainInterface(data$measurements, data$classes, verbose = 0)
  expect_s3_class(model, "lognet")
  predicted <- ClassifyR:::penalisedGLMpredictInterface(model, data$measurements[1:5, ], verbose = 0)
  expect_identical(colnames(predicted), c("class", "A", "B"))
  expect_equal(rowSums(predicted[, c("A", "B")]), rep(1, 5), ignore_attr = TRUE)
  expect_identical(as.character(predicted[["class"]]), c("A", "B")[(predicted[["B"]] > 0.5) + 1])
  one <- ClassifyR:::penalisedGLMpredictInterface(model, data$measurements[1, , drop = FALSE], returnType = "score", verbose = 0)
  expect_equal(unname(one), unname(as.matrix(predicted[1, c("A", "B")])))
  threeClasses <- factor(rep(c("A", "B", "C"), length.out = nrow(data$measurements)))
  model <- ClassifyR:::penalisedGLMtrainInterface(data$measurements, threeClasses, verbose = 0)
  expect_s3_class(model, "multnet")
  expect_identical(colnames(ClassifyR:::penalisedGLMpredictInterface(model, data$measurements[1:5, ], returnType = "score", verbose = 0)), c("A", "B", "C"))
})

test_that("penalised GLM chooses lambda by cross-validation unless resubstitution is asked for", {
  set.seed(8)
  classes <- factor(rep(c("A", "B"), each = 40))
  measurements <- matrix(rnorm(80 * 30), 80, 30, dimnames = list(paste0("s", 1:80), paste0("g", 1:30)))
  measurements[classes == "B", 1:3] <- measurements[classes == "B", 1:3] + 1
  measurements <- asDataFrame(measurements)
  set.seed(9)
  model <- ClassifyR:::penalisedGLMtrainInterface(measurements, classes, alpha = 0.5, verbose = 0)
  set.seed(9)
  again <- ClassifyR:::penalisedGLMtrainInterface(measurements, classes, alpha = 0.5, verbose = 0)
  expect_identical(attr(model, "tune"), attr(again, "tune"))
  expect_true(attr(model, "tune")[["lambda"]] %in% model[["lambda"]])

  # Resubstitution: the smallest balanced error of the training samples, the largest lambda among ties.
  resubstitution <- ClassifyR:::penalisedGLMtrainInterface(measurements, classes, alpha = 0.5, lambdaTuning = "resubstitution", verbose = 0)
  lambdas <- resubstitution[["lambda"]][-1]
  errors <- sapply(lambdas, function(lambda)
    calcExternalPerformance(classes, factor(predict(resubstitution, as.matrix(measurements), s = lambda, type = "class"), levels = levels(classes)), "Balanced Error"))
  expect_equal(attr(resubstitution, "tune")[["lambda"]], lambdas[which.min(errors)])
})

test_that("penalised GLM considers every lambda at which any class has a coefficient", {
  set.seed(10)
  classes <- factor(rep(c("A", "B", "C"), each = 50))
  measurements <- matrix(rnorm(150 * 40), 150, 40, dimnames = list(paste0("s", 1:150), paste0("g", 1:40)))
  measurements[classes == "B", 1:5] <- measurements[classes == "B", 1:5] + 1
  measurements[classes == "C", 6:10] <- measurements[classes == "C", 6:10] + 1
  path <- glmnet::glmnet(measurements, classes, family = "multinomial", alpha = 0.5)
  firstClassEmpty <- colSums(abs(as.matrix(path[["beta"]][["A"]]))) == 0
  anyClassUsed <- Reduce(`|`, lapply(path[["beta"]], function(coefficients) colSums(abs(as.matrix(coefficients))) != 0))
  # Lambdas at which class A, the first, has no coefficients but classes B and C do.
  lambdas <- path[["lambda"]][firstClassEmpty & anyClassUsed]
  expect_gt(length(lambdas), 1)
  model <- ClassifyR:::penalisedGLMtrainInterface(asDataFrame(measurements), classes, lambda = lambdas, alpha = 0.5,
                                                  lambdaTuning = "resubstitution", verbose = 0)
  expect_true(attr(model, "tune")[["lambda"]] %in% lambdas)
})

test_that("penalised GLM falls back to resubstitution for a class with fewer than three samples", {
  set.seed(11)
  classes <- factor(c(rep("A", 20), "B", "B"))
  measurements <- asDataFrame(matrix(rnorm(22 * 5), 22, 5, dimnames = list(paste0("s", 1:22), paste0("g", 1:5))))
  expect_warning(model <- ClassifyR:::penalisedGLMtrainInterface(measurements, classes, verbose = 0), "resubstitution")
  expect_true(attr(model, "tune")[["lambda"]] %in% model[["lambda"]])
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

test_that("weighted kNN lets training samples identical to a test sample take the vote", {
  data <- makeTwoClass(shift = 5)
  train <- asDataFrame(data$measurements[1:50, ])
  # The first test sample duplicates training sample 2 (class B); the others are new.
  test <- asDataFrame(rbind(data$measurements[2, , drop = FALSE], data$measurements[51:52, ]))
  predicted <- ClassifyR:::kNNinterface(train, data$classes[1:50], test, k = 5, mode = "weighted", verbose = 0)
  expect_equal(unlist(predicted[1, c("A", "B")]), c(A = 0, B = 1))
  expect_equal(as.character(predicted[1, "class"]), "B")
  expect_false(anyNA(predicted[, c("A", "B")]))
  expect_equal(unname(rowSums(predicted[, c("A", "B")])), rep(1, 3))
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
  # The score is centred on the critical value, so its sign gives the predicted class.
  critical <- 0.5 * sum(direction * (colMeans(trainMatrix[isA, ]) + colMeans(trainMatrix[!isA, ])))
  expect_equal(unname(predicted[, "score"]), unname(critical - as.matrix(test) %*% direction)[, 1])
  expect_equal(as.character(predicted[, "class"]), ifelse(predicted[, "score"] > 0, "B", "A"))
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

test_that("naive Bayes finds crossover points of the class densities scaled by class size", {
  # Class B is three times as common. Scaled by class size, the densities of N(0, 1) and N(1, 1) cross at
  # 0.5 + log(1/3) = -0.6; unscaled, they cross at 0.5.
  set.seed(11)
  classes <- factor(rep(c("A", "B"), c(1000, 3000)))
  values <- c(rnorm(1000, 0), rnorm(3000, 1))
  train <- asDataFrame(matrix(values, ncol = 1, dimnames = list(NULL, "g1")))
  test <- asDataFrame(matrix(c(-1.5, 0.5), ncol = 1, dimnames = list(c("t1", "t2"), "g1")))
  predicted <- ClassifyR:::naiveBayesKernel(train, classes, test, weighting = "crossover distance",
                                            minDifference = 0.3, verbose = 0)
  # Both samples are far enough from the crossover of the scaled densities to vote. Measured from the
  # unscaled crossover, t2 would be too close and get the class proportions as its scores.
  expect_equal(as.character(predicted[, "class"]), c("A", "B"))
  expect_equal(unname(unlist(predicted[2, c("A", "B")])), c(0, 1))
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

test_that("previousSelection matches non-syntactic feature names and needs no intermediate setting", {
  data <- makeTwoClass(shift = 3)
  measurements <- data$measurements
  # "g-1" and "g.1" are different features with the same syntactic name; both carry the class difference.
  colnames(measurements) <- paste0("g-", seq_len(ncol(measurements)))
  colnames(measurements)[2] <- "g.1"
  crossValParams <- CrossValParams(permutations = 1, folds = 2, parallelParams = SerialParam(RNGseed = 1))
  first <- suppressWarnings(runTests(asDataFrame(measurements), data$classes, crossValParams,
             ModellingParams(selectParams = SelectParams("t-test", nFeatures = 3), balancing = "none"), verbose = 0))
  expect_setequal(chosenFeatureNames(first)[[1]], c("g-1", "g.1", "g-3"))
  second <- suppressWarnings(runTests(asDataFrame(measurements), data$classes, crossValParams,
              ModellingParams(selectParams = SelectParams("previousSelection", classifyResult = first), balancing = "none"),
              verbose = 0))
  expect_equal(lapply(chosenFeatureNames(second), sort), lapply(chosenFeatureNames(first), sort))
})

test_that("ensemble selection keeps features ranked highly by enough of the ranking functions", {
  data <- makeTwoClass(shift = 2)
  measurements <- asDataFrame(data$measurements)
  crossValParams <- CrossValParams(permutations = 1, folds = 2, parallelParams = SerialParam(RNGseed = 1))
  # Both rankings agree on the top three, which carry the class difference.
  ensemble <- SelectParams(list("t-test", "limma"), nFeatures = 3, minPresence = 2)
  result <- runTests(measurements, data$classes, crossValParams,
                     ModellingParams(selectParams = ensemble, balancing = "none"), verbose = 0)
  expect_s4_class(result, "ClassifyResult")
  for(chosen in chosenFeatureNames(result)) expect_setequal(chosen, paste0("g", 1:3))
  # With tuning, the number of top features is chosen by resubstitution.
  ensembleTuned <- SelectParams(list("t-test", "limma"), minPresence = 2, tuneParams = list(nFeatures = c(3, 10)))
  tuned <- suppressWarnings(runTests(measurements, data$classes,
             CrossValParams(permutations = 1, folds = 2, tuneMode = "Resubstitution", parallelParams = SerialParam(RNGseed = 1)),
             ModellingParams(selectParams = ensembleTuned, balancing = "none"), verbose = 0))
  expect_true(all(lengths(chosenFeatureNames(tuned)) >= 3))
  expect_equal(colnames(tunedParameters(tuned)[[1]][["tuneCombinations"]])[1], "topN")
  # Nested-CV tuning uses the inner scheme of the training set.
  nested <- suppressWarnings(runTests(measurements, data$classes,
              CrossValParams(permutations = 1, folds = 2, tuneMode = "Nested CV", innerFolds = 2, parallelParams = SerialParam(RNGseed = 1)),
              ModellingParams(selectParams = ensembleTuned, balancing = "none"), verbose = 0))
  expect_s4_class(nested, "ClassifyResult")
  expect_true(all(lengths(chosenFeatureNames(nested)) >= 3))
})

test_that("two-class rankings stop for more than two classes", {
  data <- makeTwoClass()
  train <- asDataFrame(data$measurements)
  threeClasses <- factor(rep(c("A", "B", "C"), length.out = nrow(train)))
  pairs <- S4Vectors::Pairs(c("g1", "g2"), c("g10", "g11"))
  expect_error(ClassifyR:::KolmogorovSmirnovRanking(train, threeClasses, verbose = 0), "two classes")
  expect_error(ClassifyR:::KullbackLeiblerRanking(train, threeClasses, verbose = 0), "two classes")
  expect_error(ClassifyR:::pairsDifferencesRanking(train, threeClasses, featurePairs = pairs, verbose = 0), "two classes")
  expect_length(ClassifyR:::pairsDifferencesRanking(train, data$classes, featurePairs = pairs, verbose = 0), 2)
})

test_that("Levene ranking agrees with car::leveneTest", {
  skip_if_not_installed("car")
  set.seed(13)
  classes <- factor(rep(c("A", "B", "C"), length.out = 45))
  measurements <- matrix(rnorm(45 * 30, sd = rep(c(1, 2, 3), length.out = 45)), 45, 30, dimnames = list(NULL, paste0("g", 1:30)))
  measurements[, 1:10] <- rnorm(45 * 10)
  pValues <- apply(measurements, 2, function(featureColumn) car::leveneTest(featureColumn, classes)[["Pr(>F)"]][1])
  expect_equal(ClassifyR:::leveneRanking(asDataFrame(measurements), classes, verbose = 0), order(pValues))
  expect_equal(ClassifyR:::leveneRanking(asDataFrame(measurements[, 1, drop = FALSE]), classes, verbose = 0), 1)
})

test_that("likelihood ratio ranking agrees with the sum of normal log densities", {
  data <- makeTwoClass()
  measurements <- data$measurements
  measurements[data$classes == "B", 4:6] <- 2 * measurements[data$classes == "B", 4:6]
  logLikelihood <- function(values) sum(dnorm(values, mean(values), sd(values), log = TRUE))
  statistics <- apply(measurements, 2, logLikelihood) -
                Reduce(`+`, lapply(levels(data$classes), function(class) apply(measurements[data$classes == class, ], 2, logLikelihood)))
  expect_identical(ClassifyR:::likelihoodRatioRanking(asDataFrame(measurements), data$classes, verbose = 0), order(statistics))
  expect_identical(ClassifyR:::likelihoodRatioRanking(asDataFrame(measurements[, 1, drop = FALSE]), data$classes, verbose = 0), 1L)
})

test_that("XGB fits 100 rounds by default", {
  skip_if_not_installed("xgboost")
  data <- makeTwoClass()
  model <- ClassifyR:::extremeGradientBoostingTrainInterface(asDataFrame(data$measurements), data$classes, verbose = 0)
  expect_identical(xgboost::xgb.get.num.boosted.rounds(model), 100L)
})
