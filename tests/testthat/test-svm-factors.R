# SVM with factor features: fitted with the formula interface and predicted on a data frame.

test_that("SVM cross-validates data with factor features", {
  data <- makeTwoClass()
  measurements <- data.frame(data$measurements[, 1:3],
                             group = factor(rep(c("x", "y", "z"), length.out = nrow(data$measurements))),
                             check.names = FALSE)
  set.seed(1)
  result <- crossValidate(measurements, data$classes, classifier = "SVM", nRepeats = 2, nFolds = 3)
  expect_s4_class(result, "ClassifyResult")
  expect_equal(nrow(predictions(result)), 2 * nrow(measurements))
})

test_that("SVM models carry their matrix features only when fitted to a matrix", {
  data <- makeTwoClass()
  numericModel <- ClassifyR:::SVMtrainInterface(as.data.frame(data$measurements[, 1:3]), data$classes, verbose = 0)
  expect_equal(attr(numericModel, "features", exact = TRUE), colnames(data$measurements)[1:3])
  mixed <- data.frame(data$measurements[, 1:2], group = factor(rep(c("x", "y"), length.out = nrow(data$measurements))))
  formulaModel <- ClassifyR:::SVMtrainInterface(mixed, data$classes, verbose = 0)
  attr(formulaModel, "featuresForTrain") <- colnames(mixed) # As runTest sets on every model.
  expect_null(attr(formulaModel, "features", exact = TRUE))
  predicted <- ClassifyR:::SVMpredictInterface(formulaModel, mixed, verbose = 0)
  expect_equal(nrow(predicted), nrow(mixed))
})
