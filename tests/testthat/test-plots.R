# Plotting functions: options which used to fail or give wrong pictures, and the user's ggplot2 theme.

makePlotResult <- function(name = "x", classifier = "DLDA", selection = "t-test", nSamples = 30, nPermutations = 3,
                           classLevels = c("No", "Yes"), seed = 1, ranked = NULL)
{
  set.seed(seed)
  sampleIDs <- paste0("s", seq_len(nSamples))
  classes <- factor(rep(classLevels, length.out = nSamples), levels = classLevels)
  predictions <- do.call(rbind, lapply(seq_len(nPermutations), function(permutation)
  {
    scores <- matrix(runif(nSamples * length(classLevels)), nSamples, dimnames = list(NULL, classLevels))
    scores[cbind(seq_len(nSamples), as.integer(classes))] <- scores[cbind(seq_len(nSamples), as.integer(classes))] + 0.4
    scores <- scores / rowSums(scores)
    data.frame(sample = sampleIDs, permutation = permutation, fold = rep(1:3, each = ceiling(nSamples / 3))[seq_len(nSamples)],
               class = factor(classLevels[max.col(scores)], levels = classLevels), scores, check.names = FALSE)
  }))
  features <- paste0("g", 1:40)
  if(is.null(ranked)) ranked <- lapply(seq_len(nPermutations * 3), function(index) sample(features))
  ClassifyResult(S4Vectors::DataFrame(characteristic = c("Assay Name", "Classifier Name", "Selection Name", "Cross-validation"),
                                      value = c(name, classifier, selection, paste(nPermutations, "Permutations, 3 Folds"))),
                 sampleIDs, features, ranked, lapply(ranked, head, 5), list(function(oracle){}), NULL,
                 S4Vectors::DataFrame(predictions, check.names = FALSE), classes)
}

# Draws a grob or ggplot to a PNG file and returns the file's checksum.
drawChecksum <- function(drawCode)
{
  fileName <- tempfile(fileext = ".png")
  png(fileName, 800, 500)
  drawCode
  dev.off()
  unname(tools::md5sum(fileName))
}

test_that("plotting functions leave the user's ggplot2 theme unchanged", {
  pdf(NULL)
  on.exit(dev.off())
  ggplot2::theme_set(ggplot2::theme_grey())
  results <- list(makePlotResult("x"), makePlotResult("y", seed = 2))
  results <- lapply(results, calcCVperformance, performanceTypes = "Balanced Accuracy")
  performancePlot(results)
  ROCplot(results)
  rankingPlot(results, topRanked = 1:10, xLabelPositions = 1:10, parallelParams = BiocParallel::SerialParam())
  selectionPlot(results, parallelParams = BiocParallel::SerialParam())
  data <- makeTwoClass()
  plotFeatureClasses(data$measurements, data$classes, useFeatures = "g1")
  expect_identical(ggplot2::theme_get(), ggplot2::theme_grey())
})

test_that("samplesMetricMap matches samples by name between results", {
  first <- makePlotResult("x")
  second <- makePlotResult("y", seed = 2)
  reordered <- second
  newOrder <- rev(seq_along(sampleNames(second)))
  reordered@originalNames <- sampleNames(second)[newOrder]
  reordered@actualOutcome <- actualOutcome(second)[newOrder]
  sameOrder <- drawChecksum(suppressWarnings(samplesMetricMap(list(first, second))))
  otherOrder <- drawChecksum(suppressWarnings(samplesMetricMap(list(first, reordered))))
  expect_identical(otherOrder, sameOrder)
})

test_that("samplesMetricMap works for a single result, more than two classes, no legends and cross-validation comparison", {
  pdf(NULL)
  on.exit(dev.off())
  single <- suppressWarnings(samplesMetricMap(makePlotResult("x")))
  expect_s3_class(single, "gtable")
  threeClasses <- list(makePlotResult("x", classLevels = c("a", "b", "c")), makePlotResult("y", classLevels = c("a", "b", "c"), seed = 2))
  expect_s3_class(suppressWarnings(samplesMetricMap(threeClasses)), "gtable")
  expect_error(suppressWarnings(samplesMetricMap(threeClasses, metricColours = list(c("white", "blue"), c("white", "red")))), "3 classes")
  twoResults <- list(makePlotResult("x"), makePlotResult("y", seed = 2))
  expect_s3_class(suppressWarnings(samplesMetricMap(twoResults, showLegends = FALSE)), "gtable")
  groups <- factor(rep(c("M", "F"), 15))
  names(groups) <- sampleNames(twoResults[[1]])
  expect_s3_class(suppressWarnings(samplesMetricMap(twoResults, showLegends = FALSE, featureValues = groups, featureName = "Sex")), "gtable")
  differentCV <- twoResults
  differentCV[[2]]@characteristics[4, "value"] <- "5 Permutations, 3 Folds"
  expect_s3_class(suppressWarnings(samplesMetricMap(differentCV, comparison = "Cross-validation")), "gtable")
})

test_that("samplesMetricMap for a matrix keeps the first and last class tiles", {
  pdf(NULL)
  on.exit(dev.off())
  metrics <- matrix(runif(20), 2, 10, dimnames = list(c("A", "B"), paste0("s", 1:10)))
  classes <- factor(rep(c("No", "Yes"), 5))
  expect_no_warning(samplesMetricMap(metrics, classes))
  expect_no_error(samplesMetricMap(metrics, classes, showLegends = FALSE))
})

test_that("ROCplot draws averaged curves", {
  result <- makePlotResult("x", classLevels = c("a", "b", "c"))
  averaged <- suppressWarnings(ROCplot(result, mode = "average"))
  expect_s3_class(averaged, "ggplot")
  expect_setequal(sub(" .*", '', unique(as.character(averaged$data[, "class"]))), c("a", "b", "c"))
  expect_true(all(averaged$data[, "lower"] <= averaged$data[, "upper"]))
})

test_that("performancePlot uses the user's yLimits when rotated", {
  results <- lapply(list(makePlotResult("x"), makePlotResult("y", seed = 2)), calcCVperformance, performanceTypes = "Balanced Accuracy")
  rotated <- performancePlot(results, yLimits = c(0.3, 0.9), rotate90 = TRUE)
  expect_equal(rotated$coordinates$limits$y, c(0.3, 0.9))
  expect_false(anyNA(rotated$layers[[2]]$data[, "Assay Name"]))
})

test_that("performancePlot draws the chance level of the metric", {
  results <- lapply(list(makePlotResult("x"), makePlotResult("y", seed = 2)), calcCVperformance, performanceTypes = "Balanced Accuracy")
  chanceLine <- function(plot) unname(unlist(lapply(plot$layers, function(layer) if(is(layer$geom, "GeomHline")) layer$data$yintercept)))
  nClasses <- length(levels(actualOutcome(results[[1]])))
  expect_equal(chanceLine(performancePlot(results)), 1 / nClasses)
  expect_equal(chanceLine(suppressWarnings(performancePlot(results, metric = "Balanced Error"))), 1 - 1 / nClasses)
  expect_equal(ClassifyR:::.chanceLevel("AUC", results[[1]]), 0.5)
  expect_equal(ClassifyR:::.chanceLevel("Matthews Correlation Coefficient", results[[1]]), 0)
  expect_true(is.na(ClassifyR:::.chanceLevel("Unknown", results[[1]])))
})

test_that("ranking overlap divides by the length of the shorter list of top features", {
  expect_equal(ClassifyR:::.topOverlap(paste0("g", 1:3), paste0("g", c(1, 2, 9, 10)), 5), 2 / 3 * 100)
  expect_equal(ClassifyR:::.topOverlap(paste0("g", 1:10), paste0("g", c(1, 2, 11:18)), 5), 2 / 5 * 100)
})

test_that("rankingPlot uses its fonts, row and column characteristics, ordering and short rankings", {
  results <- list(makePlotResult("x", classifier = "DLDA", selection = "t-test"),
                  makePlotResult("x", classifier = "SVM", selection = "t-test", seed = 2),
                  makePlotResult("x", classifier = "DLDA", selection = "limma", seed = 3),
                  makePlotResult("x", classifier = "SVM", selection = "limma", seed = 4))
  serial <- BiocParallel::SerialParam()
  plot <- rankingPlot(results, topRanked = 1:10, xLabelPositions = 1:10, comparison = "Classifier Name",
                      characteristicsList = list(row = "Selection Name"), parallelParams = serial)
  expect_equal(plot$theme$plot.title$size, 24)
  expect_setequal(unique(plot$data[, "Selection Name"]), c("t-test", "limma"))
  ordered <- rankingPlot(results, topRanked = 1:10, xLabelPositions = 1:10, characteristicsList = list(lineColour = "Classifier Name"),
                         orderingList = list(`Classifier Name` = c("SVM", "DLDA")), parallelParams = serial)
  expect_identical(levels(ordered$data[, "Classifier Name"]), c("SVM", "DLDA"))

  shortRankings <- list(paste0("g", 1:3), paste0("g", 4:6), paste0("g", 7:9))
  short <- makePlotResult("x", nPermutations = 1, ranked = shortRankings)
  overlaps <- rankingPlot(short, topRanked = 1:10, xLabelPositions = 1:10, parallelParams = serial)$data
  expect_true(all(overlaps[, "overlap"] == 0))
})

test_that("selectionPlot draws variable importance", {
  results <- lapply(1:2, function(index)
  {
    result <- makePlotResult(c("x", "y")[index], seed = index)
    result@importance <- S4Vectors::DataFrame(feature = rep(paste0("g", 1:4), 3), `Change in Balanced Error` = rnorm(12), check.names = FALSE)
    result
  })
  importance <- selectionPlot(results, comparison = "importance", characteristicsList = list(x = "Assay Name"))
  expect_s3_class(importance, "ggplot")
  expect_equal(nrow(importance$data), 24)
  expect_identical(rlang::as_label(importance$mapping$y), "Change in Balanced Error")
  facetted <- selectionPlot(results, comparison = "importance", characteristicsList = list(x = "Assay Name", row = "Assay Name"))
  expect_s3_class(facetted, "gtable")
})

test_that("featureSetSummary reports the reduction of sets in a readable message", {
  sets <- FeatureSetCollection(list(A = c("Gene 1", "Gene 2"), B = c("Gene 20", "Gene 21")))
  genes <- S4Vectors::DataFrame(matrix(rnorm(30), 3, 10, dimnames = list(NULL, paste("Gene", 1:10))), check.names = FALSE)
  expect_message(featureSetSummary(genes, featureSets = sets), "reducing 2 feature sets to 1 feature sets")
})

test_that("plotFeatureClasses fails for features which are not in the data", {
  data <- makeTwoClass()
  expect_error(plotFeatureClasses(data$measurements, data$classes, useFeatures = "notAFeature"), "useFeatures")
})

test_that("selectionPlot compares feature selections between results", {
  results <- list(makePlotResult("x", classifier = "DLDA"), makePlotResult("x", classifier = "SVM", seed = 2),
                  makePlotResult("x", classifier = "kNN", seed = 3))
  serial <- BiocParallel::SerialParam()
  common <- selectionPlot(results, comparison = "Classifier Name", characteristicsList = list(x = "Classifier Name"), parallelParams = serial)
  expect_false(any(duplicated(colnames(common$data))))
  expect_equal(nrow(common$data), 3 * 2 * 9 * 9) # Each result against two others, nine selections each.
  versusDLDA <- selectionPlot(results, comparison = "Classifier Name", referenceLevel = "DLDA",
                              characteristicsList = list(x = "Classifier Name"), parallelParams = serial)
  expect_setequal(unique(versusDLDA$data[, "Classifier Name"]), c("SVM", "kNN"))
  expect_equal(nrow(versusDLDA$data), 2 * 9 * 9)
})
