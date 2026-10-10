# Performance metrics: values checked against direct implementations of their definitions.

# A ClassifyResult with class scores, from repeated cross-validation.
makeScoresResult <- function(nSamples = 40, nPermutations = 3, nFolds = 4, digits = NULL, seed = 1)
{
  set.seed(seed)
  sampleIDs <- paste0("s", seq_len(nSamples))
  classes <- factor(rep(c("No", "Yes"), length.out = nSamples), levels = c("No", "Yes"))
  predictions <- do.call(rbind, lapply(seq_len(nPermutations), function(permutation)
  {
    score <- plogis(as.integer(classes) - 1.5 + rnorm(nSamples))
    if(!is.null(digits)) score <- round(score, digits)
    data.frame(sample = sampleIDs, permutation = permutation, fold = sample(rep(seq_len(nFolds), length.out = nSamples)),
               class = factor(ifelse(score > 0.5, "Yes", "No"), levels = levels(classes)), No = 1 - score, Yes = score)
  }))
  ClassifyResult(S4Vectors::DataFrame(characteristic = "Assay Name", value = "x"), sampleIDs, "g1", list("g1"), list("g1"),
                 list(function(oracle){}), NULL, S4Vectors::DataFrame(predictions), classes)
}

# A ClassifyResult with risk scores, with ties in times and risks.
makeRiskResult <- function(nSamples = 30, nPermutations = 3, seed = 2)
{
  set.seed(seed)
  sampleIDs <- paste0("s", seq_len(nSamples))
  times <- round(rexp(nSamples, 0.1)) + 1
  times[1] <- 0.5 # Censored before every event, so it has no comparable pairs.
  events <- rbinom(nSamples, 1, 0.6)
  events[1] <- 0
  outcome <- survival::Surv(times, events)
  predictions <- do.call(rbind, lapply(seq_len(nPermutations), function(permutation)
    data.frame(sample = sample(sampleIDs), permutation = permutation, fold = rep(1:3, length.out = nSamples),
               risk = round(rnorm(nSamples), 1))))
  ClassifyResult(S4Vectors::DataFrame(characteristic = "Assay Name", value = "x"), sampleIDs, "g1", list("g1"), list("g1"),
                 list(function(oracle){}), NULL, S4Vectors::DataFrame(predictions), outcome)
}

# Sample-wise C-index computed pair by pair. Pairs with tied times are not comparable; tied risks count as half concordant.
referenceSampleCindex <- function(result)
{
  predictions <- as.data.frame(predictions(result))
  outcome <- as.matrix(actualOutcome(result))
  rownames(outcome) <- sampleNames(result)
  counts <- matrix(0, length(sampleNames(result)), 2, dimnames = list(sampleNames(result), c("concordant", "comparable")))
  for(permutationPredictions in split(predictions, predictions[, "permutation"]))
  {
    for(sampleID in permutationPredictions[, "sample"])
    {
      for(otherID in setdiff(permutationPredictions[, "sample"], sampleID))
      {
        times <- outcome[c(sampleID, otherID), "time"]
        events <- outcome[c(sampleID, otherID), "status"]
        risks <- permutationPredictions[match(c(sampleID, otherID), permutationPredictions[, "sample"]), "risk"]
        earlier <- which.min(times)
        if(times[1] == times[2] || events[earlier] != 1) next
        concordant <- if(risks[1] == risks[2]) 0.5 else as.numeric(risks[earlier] > risks[-earlier])
        counts[sampleID, ] <- counts[sampleID, ] + c(concordant, 1)
      }
    }
  }
  Cindex <- counts[, "concordant"] / counts[, "comparable"]
  Cindex[is.nan(Cindex)] <- NA
  Cindex
}

# AUC of each class by the trapezoid rule over the ROC curve, averaged over classes.
referenceAUC <- function(classes, scores)
{
  mean(sapply(levels(classes), function(class)
  {
    thresholds <- sort(unique(scores[, class]), decreasing = TRUE)
    TPR <- c(0, sapply(thresholds, function(threshold) mean(scores[classes == class, class] >= threshold)))
    FPR <- c(0, sapply(thresholds, function(threshold) mean(scores[classes != class, class] >= threshold)))
    sum(diff(FPR) * (TPR[-1] + TPR[-length(TPR)]) / 2)
  }))
}

test_that("sample-wise C-index matches the pair-by-pair definition, including samples without comparable pairs", {
  result <- makeRiskResult()
  values <- performance(calcCVperformance(result, "Sample C-index"))[["Sample C-index"]]
  expect_equal(values, referenceSampleCindex(result))
  expect_true(is.na(values[["s1"]]))
  expect_identical(names(values), sampleNames(result))
})

test_that("sample-wise C-index counts tied risks as half concordant and is not rounded", {
  outcome <- survival::Surv(c(1, 2, 3), c(1, 1, 1))
  predictions <- S4Vectors::DataFrame(sample = c("a", "b", "c"), permutation = 1, fold = 1, risk = c(2, 2, 1))
  result <- ClassifyResult(S4Vectors::DataFrame(characteristic = "Assay Name", value = "x"), c("a", "b", "c"), "g1",
                           list("g1"), list("g1"), list(function(oracle){}), NULL, predictions, outcome)
  values <- performance(calcCVperformance(result, "Sample C-index"))[["Sample C-index"]]
  # Pairs: a-b tied (1/2), a-c concordant, b-c concordant.
  expect_equal(unname(values), c(0.75, 0.75, 1))
  expect_equal(survival::concordance(outcome ~ I(-c(2, 2, 1)))$concordance, 2.5 / 3)
})

test_that("sample-wise C-index is NA for a sample which was never predicted", {
  result <- makeRiskResult()
  result@predictions <- result@predictions[result@predictions[, "sample"] != "s5", ]
  values <- performance(calcCVperformance(result, "Sample C-index"))[["Sample C-index"]]
  expect_true(is.na(values[["s5"]]))
  expect_equal(length(values), length(sampleNames(result)))
})

test_that("sample-wise C-index of a realistic result takes well under a second per hundred permutations", {
  skip_on_cran()
  result <- makeRiskResult(nSamples = 165, nPermutations = 100)
  expect_lt(system.time(calcCVperformance(result, "Sample C-index"))[["elapsed"]], 5)
})

test_that("AUC matches the trapezoid rule, with and without tied scores", {
  for(digits in list(NULL, 1))
  {
    result <- makeScoresResult(digits = digits)
    predictions <- as.data.frame(predictions(result))
    classes <- actualOutcome(result)[match(predictions[, "sample"], sampleNames(result))]
    expected <- sapply(split(seq_len(nrow(predictions)), predictions[, "permutation"]), function(rows)
      referenceAUC(classes[rows], predictions[rows, c("No", "Yes")]))
    expect_equal(performance(calcCVperformance(result, "AUC"))[["AUC"]], expected)
    expect_equal(unname(calcExternalPerformance(classes[1:40], predictions[1:40, c("No", "Yes")], "AUC")), unname(expected[1]))
  }
})

test_that("AUC with three classes is the mean of one-versus-rest AUCs", {
  set.seed(4) # No class's AUC lies on a rounding boundary.
  classes <- factor(sample(c("a", "b", "c"), 30, replace = TRUE))
  scores <- matrix(runif(90), 30, dimnames = list(NULL, levels(classes)))
  scores[cbind(1:30, as.integer(classes))] <- scores[cbind(1:30, as.integer(classes))] + 0.5
  scores <- as.data.frame(scores / rowSums(scores))
  expect_equal(unname(calcExternalPerformance(classes, scores, "AUC")), referenceAUC(classes, scores))
})

test_that("AUC is NaN with a warning if a class is absent from a group", {
  result <- makeScoresResult(nSamples = 8, nPermutations = 1, nFolds = 2)
  result@predictions <- result@predictions[order(as.character(actualOutcome(result))[match(result@predictions[, "sample"], sampleNames(result))]), ]
  result@predictions[, "fold"] <- rep(1:2, each = 4) # Each fold has samples of only one class.
  expect_warning(AUC <- performance(calcCVperformance(result, "AUC", grouping = "fold"))[["AUC"]], "2 of 2 groups")
  expect_true(is.nan(AUC[["1"]]))
  byPermutation <- performance(calcCVperformance(result, "AUC"))[["AUC"]] # Each permutation has both classes.
  expect_equal(unname(byPermutation), referenceAUC(actualOutcome(result)[match(result@predictions[, "sample"], sampleNames(result))],
                                                   as.data.frame(result@predictions[, c("No", "Yes")])))
})

test_that("AUC is not rounded", {
  set.seed(7)
  classes <- factor(rep(c("No", "Yes"), 13))
  scores <- data.frame(No = runif(26))
  scores[, "Yes"] <- 1 - scores[, "No"]
  value <- unname(calcExternalPerformance(classes, scores, "AUC"))
  expect_equal(value, referenceAUC(classes, scores))
  expect_false(isTRUE(all.equal(value, round(value, 2))))
})

test_that("macro precision and F1 count a class which is never predicted as precision 0", {
  actual <- factor(c("A", "A", "B", "B", "C", "C"))
  predicted <- factor(c("A", "A", "B", "A", "B", "B"), levels = levels(actual)) # C is never predicted.
  values <- calcExternalPerformance(actual, predicted, c("Macro Precision", "Macro Recall", "Macro F1"))
  precision <- mean(c(2 / 3, 1 / 3, 0))
  recall <- mean(c(1, 1 / 2, 0))
  expect_equal(values[["Macro Precision"]], precision)
  expect_equal(values[["Macro Recall"]], recall)
  expect_equal(values[["Macro F1"]], 2 * precision * recall / (precision + recall))
})

test_that("grouping by fold averages balanced accuracy and other class metrics within each permutation", {
  result <- makeScoresResult(nPermutations = 3)
  predictions <- as.data.frame(predictions(result))
  classes <- actualOutcome(result)[match(predictions[, "sample"], sampleNames(result))]
  for(metric in c("Balanced Accuracy", "Macro F1"))
  {
    byFold <- performance(calcCVperformance(result, metric, grouping = "fold"))[[metric]]
    expect_identical(names(byFold), as.character(1:3))
    rows <- which(predictions[, "permutation"] == 2)
    foldValues <- sapply(split(rows, predictions[rows, "fold"]), function(foldRows)
      unname(calcExternalPerformance(classes[foldRows], predictions[foldRows, "class"], metric)))
    expect_equal(byFold[["2"]], mean(foldValues))
  }
  expect_length(performance(calcCVperformance(result, "Sample Accuracy", grouping = "fold"))[["Sample Accuracy"]], 40)
})

test_that("AUC of a realistic result takes well under a second per hundred permutations", {
  skip_on_cran()
  result <- makeScoresResult(nSamples = 165, nPermutations = 100, nFolds = 5)
  expect_lt(system.time(calcCVperformance(result, "AUC"))[["elapsed"]], 3)
})

test_that("calcExternalPerformance matches predicted classes to actual classes by name", {
  actual <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  predicted <- factor(c("A", "A", "B", "B"), levels = c("B", "A"))
  expect_equal(unname(calcExternalPerformance(actual, predicted, "Accuracy")), 1)
  onlyB <- factor(c("B", "B", "B", "B"))
  expect_equal(unname(calcExternalPerformance(actual, onlyB, "Accuracy")), 0.5)
  expect_error(calcExternalPerformance(actual, factor(c("A", "C", "B", "B")), "Accuracy"), "C")
})

test_that("calcExternalPerformance calculates several metrics at once", {
  actual <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  predicted <- factor(c("A", "B", "B", "B"), levels = c("A", "B"))
  values <- calcExternalPerformance(actual, predicted, c("Accuracy", "Balanced Accuracy", "Matthews Correlation Coefficient"))
  expect_equal(values[["Accuracy"]], 0.75)
  expect_equal(values[["Balanced Accuracy"]], 0.75)
  expect_equal(values[["Matthews Correlation Coefficient"]], (2 * 1 - 1 * 0) / sqrt(3 * 2 * 1 * 2))
})

test_that("grouping by fold gives a numeric vector with permutations in numerical order", {
  result <- makeScoresResult(nPermutations = 12)
  byFold <- performance(calcCVperformance(result, "AUC", grouping = "fold"))[["AUC"]]
  expect_true(is.numeric(byFold) && is.null(dim(byFold)))
  expect_identical(names(byFold), as.character(1:12))
  predictions <- as.data.frame(predictions(result))
  classes <- actualOutcome(result)[match(predictions[, "sample"], sampleNames(result))]
  rows <- which(predictions[, "permutation"] == 10)
  foldAUCs <- sapply(split(rows, predictions[rows, "fold"]), function(foldRows)
    unname(calcExternalPerformance(classes[foldRows], predictions[foldRows, c("No", "Yes")], "AUC")))
  expect_equal(byFold[["10"]], mean(foldAUCs))

  risks <- makeRiskResult(nPermutations = 11)
  CbyFold <- performance(calcCVperformance(risks, "C-index", grouping = "fold"))[["C-index"]]
  expect_true(is.numeric(CbyFold) && is.null(dim(CbyFold)))
  expect_identical(names(CbyFold), as.character(1:11))
})

test_that("easyHard matches samples by name when the assay rows are in another order", {
  skip_if_not_installed("glmnet")
  result <- makeScoresResult(nSamples = 30, nPermutations = 4)
  result@predictions <- result@predictions[result@predictions[, "sample"] != "s3", ] # Its sample accuracy is NaN.
  result <- calcCVperformance(result, "Sample Accuracy")
  set.seed(5)
  clinical <- S4Vectors::DataFrame(age = rnorm(30), row.names = sampleNames(result))
  inOrder <- suppressWarnings(easyHard(list(clinical = clinical), result, "clinical", performanceType = "Sample Accuracy"))
  reversed <- suppressWarnings(easyHard(list(clinical = clinical[30:1, , drop = FALSE]), result, "clinical", performanceType = "Sample Accuracy"))
  expect_equal(reversed, inOrder)
})
