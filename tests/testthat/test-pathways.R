# Precision pathways: training, prediction of new samples, costs, summaries and plots.

makePathwaysData <- function(seed, nSamples, assays = c("RNA", "protein"))
{
  set.seed(seed)
  outcome <- factor(rep(c("Good", "Poor"), each = nSamples / 2))
  clinical <- S4Vectors::DataFrame(Outcome = outcome, Age = rnorm(nSamples) + 0.8 * (outcome == "Poor"),
                                   Stage = rnorm(nSamples), row.names = paste0("S", seed, "_", seq_len(nSamples)))
  experiments <- lapply(seq_along(assays), function(index)
  {
    values <- matrix(rnorm(30 * nSamples), 30, dimnames = list(paste0(assays[index], "f", 1:30), rownames(clinical)))
    values[1:4, outcome == "Poor"] <- values[1:4, outcome == "Poor"] + 1.2
    values
  })
  names(experiments) <- assays
  MultiAssayExperiment::MultiAssayExperiment(experiments, colData = clinical)
}
useClinical <- list(clinical = c("Age", "Stage"))

test_that(".permutations keeps a matrix when one permutation remains", {
  permutations <- ClassifyR:::.permutations(c("clinical", "RNA"), fixed = data.frame(1, "clinical"))
  expect_true(is.matrix(permutations))
  expect_equal(permutations[, 1], c("clinical", "RNA"))
  permutations <- ClassifyR:::.permutations(c("clinical", "RNA", "protein"), fixed = data.frame(1, "clinical"))
  expect_equal(dim(permutations), c(3, 2))
})

test_that("the default mode works for clinical data and one assay", {
  data <- makePathwaysData(1, 40, "RNA")
  set.seed(1)
  pathways <- precisionPathwaysTrain(data, "Outcome", useFeatures = useClinical, nRepeats = 2, nFolds = 3,
                                     minAssaySamples = 2, classifier = c(clinical = "DLDA", RNA = "DLDA"),
                                     selectionMethod = c(clinical = "none", RNA = "t-test"))
  expect_s3_class(pathways, "PrecisionPathways")
  expect_true(all(startsWith(names(pathways[["pathways"]]), "clinical")))
})

test_that("fixedAssays may be NULL or several assays, and lists without useFeatures are accepted", {
  data <- makePathwaysData(1, 40)
  set.seed(1)
  noFixed <- precisionPathwaysTrain(data, "Outcome", fixedAssays = NULL, useFeatures = useClinical, nRepeats = 2,
                                    nFolds = 3, classifier = c(clinical = "DLDA", RNA = "DLDA", protein = "DLDA"))
  expect_equal(length(noFixed[["assaysPermutations"]]), 6)
  set.seed(1)
  twoFixed <- precisionPathwaysTrain(data, "Outcome", fixedAssays = c("clinical", "protein"), useFeatures = useClinical,
                                     nRepeats = 2, nFolds = 3, classifier = c(clinical = "DLDA", RNA = "DLDA", protein = "DLDA"))
  expect_equal(twoFixed[["assaysPermutations"]], list(c("clinical", "protein", "RNA")))

  clinical <- as.data.frame(MultiAssayExperiment::colData(data))
  dataList <- list(clinical = clinical[, c("Age", "Stage")], RNA = t(MultiAssayExperiment::assay(data, "RNA")))
  set.seed(1)
  expect_warning(fromList <- precisionPathwaysTrain(dataList, clinical[, "Outcome"], nRepeats = 2, nFolds = 3,
                                                    classifier = c(clinical = "DLDA", RNA = "DLDA")), "useFeatures")
  expect_s3_class(fromList, "PrecisionPathways")
  expect_equal(fromList[["useFeatures"]][["clinical"]], c("Age", "Stage"))
})

trainSmall <- function(classifier, data = makePathwaysData(1, 60))
{
  set.seed(2)
  precisionPathwaysTrain(data, "Outcome", useFeatures = useClinical, nRepeats = 3, nFolds = 3, minAssaySamples = 3,
                         mode = "combinatorial", classifier = classifier,
                         selectionMethod = c(clinical = "none", RNA = "t-test", protein = "t-test"))
}

test_that("models are matched to classifiers and pathways by assay name, not by position", {
  newData <- makePathwaysData(2, 30)
  inOrder <- trainSmall(c(clinical = "GLM", protein = "randomForest", RNA = "DLDA"))
  expect_equal(names(inOrder[["models"]]), c("clinical", "RNA", "protein"))
  reordered <- inOrder
  reordered[["parameters"]][["classifier"]] <- inOrder[["parameters"]][["classifier"]][c("RNA", "protein", "clinical")]
  reordered[["models"]] <- inOrder[["models"]][c("protein", "clinical", "RNA")]
  set.seed(5) # Random forest breaks tied votes at random.
  predictedInOrder <- precisionPathwaysPredict(inOrder, newData, "Outcome")
  set.seed(5)
  predictedReordered <- precisionPathwaysPredict(reordered, newData, "Outcome")
  expect_equal(lapply(predictedInOrder[["pathways"]], "[[", "individuals"),
               lapply(predictedReordered[["pathways"]], "[[", "individuals"))

  # Each new sample is predicted by every fold model of every repeat.
  individuals <- predictedInOrder[["pathways"]][[1]][["individuals"]]
  expect_setequal(individuals[, "Sample ID"], rownames(MultiAssayExperiment::colData(newData)))
  # Fold models keep only what prediction needs.
  expect_null(attr(models(inOrder[["models"]][["protein"]])[[1]], "forImportance"))
})

test_that("new samples can be predicted without their classes", {
  pathways <- trainSmall(c(clinical = "DLDA", RNA = "DLDA", protein = "DLDA"))
  newData <- makePathwaysData(2, 30)
  withClasses <- precisionPathwaysPredict(pathways, newData, "Outcome")
  withoutClasses <- precisionPathwaysPredict(pathways, newData)
  expect_equal(lapply(withoutClasses[["pathways"]], function(pathway) pathway[["individuals"]][, "Predicted"]),
               lapply(withClasses[["pathways"]], function(pathway) pathway[["individuals"]][, "Predicted"]))
  costed <- calcCostsAndPerformance(withoutClasses, c(clinical = 50, RNA = 700, protein = 400))
  expect_true(all(is.na(costed[["performance"]][, "accuracy"])))
  expect_true(all(costed[["performance"]][, "cost"] > 0))
})

test_that("assay names containing '-' are kept intact", {
  data <- makePathwaysData(1, 40, c("RNA-seq", "protein"))
  set.seed(1)
  pathways <- suppressWarnings(precisionPathwaysTrain(data, "Outcome", useFeatures = useClinical, nRepeats = 2, nFolds = 3,
                                     minAssaySamples = 3, classifier = setNames(rep("DLDA", 3), c("clinical", "RNA-seq", "protein")),
                                     selectionMethod = c(clinical = "none", `RNA-seq` = "t-test", protein = "t-test")))
  expect_true(all(unlist(lapply(pathways[["pathways"]], "[[", "assays")) %in% c("clinical", "RNA-seq", "protein")))
  predicted <- suppressWarnings(precisionPathwaysPredict(pathways, makePathwaysData(2, 30, c("RNA-seq", "protein")), "Outcome"))
  expect_s3_class(predicted, "PrecisionPathways")
})

test_that("the confidence cut-off is a minimum", {
  counts <- as.table(matrix(c(18, 2, 2, 18, 17, 3), ncol = 2, byrow = TRUE,
                            dimnames = list(c("a", "b", "c"), c("Good", "Poor"))))
  expect_equal(ClassifyR:::.pathwaysConfident(counts, 0.8), c("a", "b"))
})

# A pathway built by hand: two samples predicted by clinical data and one passed on to RNA.
handPathways <- function(withTested = TRUE)
{
  individuals <- S4Vectors::DataFrame(Tier = factor(c("clinical", "clinical", "RNA"), levels = c("clinical", "RNA")),
                                      `Sample ID` = c("s1", "s2", "s3"),
                                      Predicted = factor(c("A", "B", "B"), levels = c("A", "B")),
                                      Accuracy = c(1, 1, 1), check.names = FALSE)
  pathway <- list(pathway = "clinical-RNA", assays = c("clinical", "RNA"), assaysTested = c("clinical", "RNA"),
                  individuals = individuals,
                  tiers = S4Vectors::DataFrame(Tier = c("clinical", "RNA"), `Balanced Accuracy` = c(1, 1), check.names = FALSE))
  if(!withTested) pathway[c("assays", "assaysTested")] <- NULL
  structure(list(pathways = list(`clinical-RNA` = pathway),
                 testSampleInfo = S4Vectors::DataFrame(sampleID = c("s1", "s2", "s3"), class = factor(c("A", "B", "B")))),
            class = "PrecisionPathways")
}

test_that("a sample is charged for every tier it passed through", {
  costs <- c(clinical = 50, RNA = 700)
  expect_equal(calcCostsAndPerformance(handPathways(), costs)[["performance"]][, "cost"], 50 + 50 + (50 + 700))
  expect_equal(calcCostsAndPerformance(handPathways(FALSE), costs)[["performance"]][, "cost"], 850)
  expect_equal(calcCostsAndPerformance(handPathways(), c(RNA = 700))[["performance"]][, "cost"], 700)
})

test_that("summary weights are matched by name and must sum to 1", {
  pathways <- structure(list(performance = data.frame(accuracy = c(0.9, 0.6, 0.55), cost = c(9000, 100, 50),
                                                      row.names = c("expensiveGood", "cheapA", "cheapest"))),
                        class = "PrecisionPathways")
  costFirst <- summary(pathways, weights = c(cost = 0.9, accuracy = 0.1))
  accuracyFirst <- summary(pathways, weights = c(accuracy = 0.1, cost = 0.9))
  expect_equal(costFirst, accuracyFirst)
  expect_equal(costFirst[1, "Pathway"], "cheapest")
  expect_equal(summary(pathways, weights = c(accuracy = 1, cost = 0))[1, "Pathway"], "expensiveGood")
  expect_error(summary(pathways, weights = c(accuracy = 5, cost = 5)), "sum to 1")
  expect_error(summary(pathways, weights = c(accuracy = 0.5, price = 0.5)), "named")
})

test_that("plots keep the global theme, every pathway and the true classes", {
  pathways <- handPathways()
  pathways[["performance"]] <- data.frame(accuracy = 0.3, cost = 850, row.names = "clinical-RNA")
  pathways[["models"]] <- list(ClassifyResult(S4Vectors::DataFrame(characteristic = "Classifier Name", value = "DLDA"),
                                              c("s1", "s2", "s3"), "f1", list(), list(), list(), NULL,
                                              S4Vectors::DataFrame(sample = c("s1", "s2", "s3"), class = factor(c("A", "B", "B"))),
                                              factor(c("A", "B", "B"))))
  themeBefore <- ggplot2::theme_get()
  bubbles <- bubblePlot(pathways)
  expect_identical(ggplot2::theme_get(), themeBefore)
  expect_equal(nrow(ggplot2::layer_data(bubbles, 1)), 1) # Accuracy below 0.5 is still drawn.

  strata <- strataPlot(pathways, "clinical-RNA")
  expect_identical(ggplot2::theme_get(), themeBefore)
  trueClassTiles <- ggplot2::layer_data(strata, 2)
  expect_equal(nrow(trueClassTiles), 3)
  expect_true(all(trueClassTiles[, "y"] > 2)) # Drawn above the tiers, not underneath them.
})
