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
