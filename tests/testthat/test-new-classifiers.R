# Oblique random forests (aorsf), boosted Cox (mboost), LIBLINEAR logistic regression and ncvreg.

asDF <- function(measurements) S4Vectors::DataFrame(measurements, check.names = FALSE)
# Harrell's C-index of risk scores, where a higher score means a higher risk.
riskConcordance <- function(survival, risk) survival::concordance(survival ~ risk, reverse = TRUE)[["concordance"]]

test_that("oblique random forests classify, predict risks, rank features and are reproducible", {
  skip_if_not_installed("aorsf")
  data <- makeTwoClass(shift = 2)
  set.seed(1)
  forest <- ClassifyR:::obliqueForestTrainInterface(asDF(data$measurements), data$classes, n_tree = 50, verbose = 0)
  set.seed(1)
  again <- ClassifyR:::obliqueForestTrainInterface(asDF(data$measurements), data$classes, n_tree = 50, verbose = 0)
  predicted <- ClassifyR:::obliqueForestPredictInterface(forest, asDF(data$measurements[1:5, ]), verbose = 0)
  expect_identical(colnames(predicted), c("class", "A", "B"))
  expect_equal(unname(rowSums(predicted[, c("A", "B")])), rep(1, 5))
  expect_identical(predicted, ClassifyR:::obliqueForestPredictInterface(again, asDF(data$measurements[1:5, ]), verbose = 0))
  features <- ClassifyR:::obliqueForestFeatures(forest)
  expect_true(all(1:3 %in% features[[1]][1:5]))

  survivalData <- makeSurvival()
  set.seed(2)
  survivalForest <- ClassifyR:::obliqueForestTrainInterface(asDF(survivalData$measurements), survivalData$outcome, n_tree = 50, verbose = 0)
  risks <- ClassifyR:::obliqueForestPredictInterface(survivalForest, asDF(survivalData$measurements), verbose = 0)
  expect_identical(names(risks), rownames(survivalData$measurements))
  expect_gt(riskConcordance(survivalData$outcome, risks), 0.6)
})

test_that("boosted Cox models predict higher risks for earlier events and report the features they chose", {
  skip_if_not_installed("mboost")
  survivalData <- makeSurvival()
  model <- ClassifyR:::boostedCoxTrainInterface(asDF(survivalData$measurements), survivalData$outcome, verbose = 0)
  risks <- ClassifyR:::boostedCoxPredictInterface(model, asDF(survivalData$measurements), verbose = 0)
  expect_gt(riskConcordance(survivalData$outcome, risks), 0.6)
  features <- ClassifyR:::boostedCoxFeatures(model)
  expect_identical(features[[1]][1], 1L) # The first feature carries the risk.
  expect_true(all(features[[2]] %in% seq_len(ncol(survivalData$measurements))))
  expect_no_warning(ClassifyR:::boostedCoxTrainInterface(asDF(survivalData$measurements), survivalData$outcome, mstop = 10, verbose = 0))
})

test_that("LIBLINEAR logistic regression chooses its cost by cross-validation and gives class probabilities", {
  skip_if_not_installed("LiblineaR")
  data <- makeTwoClass(shift = 2)
  set.seed(3)
  model <- ClassifyR:::LiblineaRtrainInterface(asDF(data$measurements), data$classes, verbose = 0)
  expect_true(attr(model, "tune")[["cost"]] %in% 10^(-3:2))
  predicted <- ClassifyR:::LiblineaRpredictInterface(model, asDF(data$measurements[1:5, ]), verbose = 0)
  expect_identical(colnames(predicted), c("class", "A", "B"))
  expect_equal(unname(rowSums(predicted[, c("A", "B")])), rep(1, 5))
  expect_identical(as.character(predicted[["class"]]), c("A", "B")[(predicted[["B"]] > 0.5) + 1])
  fixedCost <- ClassifyR:::LiblineaRtrainInterface(asDF(data$measurements), data$classes, type = 6, cost = 0.1, verbose = 0)
  expect_identical(attr(fixedCost, "tune")[["cost"]], 0.1)
  features <- ClassifyR:::LiblineaRfeatures(fixedCost)
  expect_true(all(features[[2]] %in% 1:30) && all(1:3 %in% features[[1]][1:3]))
  expect_error(ClassifyR:::LiblineaRtrainInterface(asDF(data$measurements), data$classes, type = 2, verbose = 0), "logistic")
  threeClasses <- factor(rep(c("A", "B", "C"), length.out = nrow(data$measurements)))
  model <- ClassifyR:::LiblineaRtrainInterface(asDF(data$measurements), threeClasses, cost = 1, verbose = 0)
  expect_identical(colnames(ClassifyR:::LiblineaRpredictInterface(model, asDF(data$measurements[1:2, ]), returnType = "score", verbose = 0)), c("A", "B", "C"))
})

test_that("ncvreg fits two classes and survival, and refuses more than two classes", {
  skip_if_not_installed("ncvreg")
  data <- makeTwoClass(shift = 2)
  set.seed(4)
  model <- ClassifyR:::ncvregTrainInterface(asDF(data$measurements), data$classes, verbose = 0)
  predicted <- ClassifyR:::ncvregPredictInterface(model, asDF(data$measurements[1:5, ]), verbose = 0)
  expect_identical(colnames(predicted), c("class", "A", "B"))
  expect_true(all(ClassifyR:::ncvregFeatures(model)[[2]] %in% 1:30))
  expect_true(1L %in% ClassifyR:::ncvregFeatures(model)[[2]])
  threeClasses <- factor(rep(c("A", "B", "C"), length.out = nrow(data$measurements)))
  expect_error(ClassifyR:::ncvregTrainInterface(asDF(data$measurements), threeClasses, verbose = 0), "two classes")

  survivalData <- makeSurvival()
  set.seed(5)
  survivalModel <- ClassifyR:::ncvregTrainInterface(asDF(survivalData$measurements), survivalData$outcome, verbose = 0)
  risks <- ClassifyR:::ncvregPredictInterface(survivalModel, asDF(survivalData$measurements), verbose = 0)
  expect_gt(riskConcordance(survivalData$outcome, risks), 0.6)
})
