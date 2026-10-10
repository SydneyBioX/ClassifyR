# Parameter classes: constructors, informative errors and show methods.

test_that("one-step classifiers give an informative error for PredictParams", {
  expect_error(PredictParams("naiveBayes"), "trains and predicts in one function")
  expect_error(PredictParams("notAClassifier"), "not a classifier keyword")
  expect_output(show(PredictParams("DLDA")), "An object of class 'PredictParams'.", fixed = TRUE)
})

test_that("SelectParams accepts 'none' and plain functions", {
  expect_null(SelectParams("none"))
  expect_s4_class(ModellingParams(selectParams = SelectParams("none")), "ModellingParams")

  fromPackage <- SelectParams(ClassifyR:::differentMeansRanking)
  expect_s4_class(fromPackage, "SelectParams")
  expect_output(show(fromPackage), "Selection Name: Difference in Means.")

  varianceRanking <- function(measurementsTrain, classesTrain, verbose = 3)
    colnames(measurementsTrain)[order(apply(measurementsTrain, 2, var), decreasing = TRUE)]
  userParams <- SelectParams(varianceRanking, nFeatures = 5)
  expect_output(show(userParams), "User-specified Ranking")
  data <- makeTwoClass()
  result <- runTests(data$measurements, data$classes, CrossValParams(permutations = 1, folds = 3, parallelParams = BiocParallel::SerialParam()),
                     ModellingParams(balancing = "none", selectParams = userParams), verbose = 0)
  expect_true(all(lengths(chosenFeatureNames(result)) == 5))
})

test_that("balancing is not accepted by TrainParams", {
  expect_error(TrainParams("DLDA", balancing = "none"), "ModellingParams")
  expect_s4_class(TrainParams("DLDA"), "TrainParams")
})

test_that("ModellingParams does not rebalance classes by default", {
  expect_identical(ModellingParams()@balancing, "none")
  expect_identical(ModellingParams(balancing = "downsample")@balancing, "downsample")
})

test_that("CrossValParams runs serially by default, reproducibly after set.seed", {
  expect_s4_class(CrossValParams()@parallelParams, "SerialParam")
  data <- makeTwoClass(shift = 1)
  measurements <- DataFrame(data$measurements, check.names = FALSE)
  runOnce <- function()
  {
    set.seed(3)
    runTests(measurements, data$classes, CrossValParams(permutations = 2, folds = 3), ModellingParams(), verbose = 0)
  }
  expect_identical(predictions(runOnce()), predictions(runOnce()))
})

test_that("show methods print complete lines", {
  ensemble <- SelectParams(list("t-test", "limma"))
  expect_output(show(ensemble), "Minimum Functions Selected By: 1.", fixed = TRUE)
  noted <- SelectParams("t-test", characteristics = S4Vectors::DataFrame(characteristic = "Note", value = "x"))
  expect_output(show(noted), "Selection Name: Difference in Means.\nNote: x.", fixed = TRUE)
  expect_output(show(TrainParams("DLDA")), "Classifier Name: Diagonal LDA.")
  expect_output(show(TransformParams("diffLoc")), "Transform Name: Subtraction From Training Set Location.")
  sets <- FeatureSetCollection(setNames(lapply(1:5, function(index) paste0("g", 1:8)), paste0("set", 1:5)))
  outputFile <- tempfile()
  sink(outputFile)
  show(sets)
  sink()
  printed <- readChar(outputFile, file.size(outputFile))
  expect_true(endsWith(printed, "set5: g1, g2, g3, g4, g5, ...\n"))
})

test_that("performance metrics record whether higher is better", {
  table <- ClassifyR:::.ClassifyRenvir[["performanceInfoTable"]]
  expect_true(all(table[, "better"] %in% c("higher", "lower")))
})
