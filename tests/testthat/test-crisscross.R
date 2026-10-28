# Cross-cohort evaluation with crissCrossValidate and its plot.

makeCohort <- function(seed, nSamples = 40, features = paste0("g", 1:60), shift = 0)
{
  set.seed(seed)
  classes <- factor(rep(c("A", "B"), length.out = nSamples))
  measurements <- matrix(rnorm(nSamples * length(features)), nSamples,
                         dimnames = list(paste0("c", seed, "s", seq_len(nSamples)), features))
  measurements[classes == "B", 1:3] <- measurements[classes == "B", 1:3] + shift
  list(measurements = measurements, classes = classes)
}

test_that("only features measured in every data set are used", {
  first <- makeCohort(1, features = paste0("g", 1:60), shift = 1.5)
  second <- makeCohort(2, features = paste0("g", c(1:50, 61:70)), shift = 1.5)
  set.seed(1)
  expect_message(result <- crissCrossValidate(list(first = first$measurements, second = second$measurements),
                                              list(first = first$classes, second = second$classes),
                                              nFeatures = 5, classifier = "DLDA", trainType = "modelTrain"),
                 "50 features present in every data set")
  expect_equal(dim(result[["real"]]), c(2, 2))
  expect_equal(result[["params"]][["diagonal"]], "resubstitution")
})

test_that("the modelTest diagonal is cross-validation within the data set", {
  # Pure noise: features chosen with the test samples would look predictive.
  first <- makeCohort(1, nSamples = 60, features = paste0("g", 1:500))
  second <- makeCohort(2, nSamples = 60, features = paste0("g", 1:500))
  set.seed(1)
  result <- crissCrossValidate(list(first = first$measurements, second = second$measurements),
                               list(first = first$classes, second = second$classes),
                               nFeatures = 20, classifier = "DLDA", trainType = "modelTest", nRepeats = 3)
  expect_equal(result[["params"]][["diagonal"]], "cross-validation")
  expect_true(all(diag(result[["real"]]) < 0.7))

  set.seed(1)
  crossValidated <- crossValidate(first$measurements, first$classes, nFeatures = 20, selectionMethod = "t-test",
                                  classifier = "DLDA", multiViewMethod = "none", nFolds = 5, nRepeats = 3)
  expect_equal(result[["real"]][1, 1],
               round(mean(performance(calcCVperformance(crossValidated, "Balanced Accuracy"))[["Balanced Accuracy"]]), 2))
})

test_that("performance type and outcome type are resolved for every data set", {
  first <- makeCohort(1, shift = 1.5)
  survival <- survival::Surv(rexp(40), rep(1, 40))
  expect_error(crissCrossValidate(list(first = first$measurements, second = first$measurements),
                                  list(first = first$classes, second = survival), classifier = "DLDA"),
               "same kind")
})

test_that("random features work with several values of nFeatures, and the plot is returned", {
  first <- makeCohort(1, shift = 1.5)
  second <- makeCohort(2, shift = 1.5)
  set.seed(1)
  result <- suppressWarnings(crissCrossValidate(list(first = first$measurements, second = second$measurements),
                               list(first = first$classes, second = second$classes),
                               nFeatures = c(5, 10), classifier = "DLDA", trainType = "modelTest",
                               nRepeats = 2, doRandomFeatures = TRUE))
  expect_equal(dim(result[["random"]]), c(2, 2))

  searchBefore <- search()
  plotted <- crissCrossPlot(result, includeValues = TRUE)
  expect_identical(search(), searchBefore)
  expect_s3_class(plotted, "ggplot")

  result[["random"]] <- NULL
  plotted <- crissCrossPlot(result)
  expect_s3_class(plotted, "ggplot")
  expect_match(plotted[["labels"]][["caption"]], "cross-validation")
})
