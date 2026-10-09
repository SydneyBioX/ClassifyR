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
