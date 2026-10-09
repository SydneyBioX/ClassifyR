# Faster implementations give the same results as the ones they replace.

test_that("the C-index of many score columns equals survival::concordance", {
  set.seed(3)
  for(repetition in 1:50)
  {
    n <- sample(5:40, 1)
    time <- if(repetition %% 2) round(rexp(n) * 5) + 1 else rexp(n) # Tied times in half of the cases.
    event <- rbinom(n, 1, 0.6)
    risk <- matrix(round(rnorm(n * 3), sample(0:2, 1)), n) # Tied scores too.
    expected <- apply(risk, 2, function(column)
                      survival::concordance(survival::Surv(time, event) ~ I(-column))$concordance)
    expect_identical(unname(ClassifyR:::.CindexColumns(risk, time, event)), unname(expected))
  }
})

test_that("Cox elastic net tuning chooses the same lambda and model as cv.glmnet", {
  skip_if_not_installed("glmnet")
  data <- makeSurvival(nSamples = 100, nFeatures = 15)
  for(seed in 1:5)
  {
    set.seed(seed)
    expected <- glmnet::cv.glmnet(data$measurements, data$outcome, family = "cox", type.measure = "C")
    set.seed(seed)
    fast <- ClassifyR:::.cvCoxnetC(data$measurements, data$outcome)
    expect_identical(fast$lambda.min, expected$lambda.min)
    expect_identical(fast$glmnet.fit$beta, expected$glmnet.fit$beta)
  }
})

test_that("DataFrames are converted to the same data.frame as by as.data.frame", {
  data <- makeTwoClass()
  measurements <- S4Vectors::DataFrame(data$measurements, check.names = FALSE)
  measurements$group <- factor(rep(c("x", "y"), length.out = nrow(measurements)))
  colnames(measurements)[1] <- "unsafe name-1"
  expect_identical(ClassifyR:::.asDataFrame(measurements), as.data.frame(measurements))
})

test_that("numeric design matrices equal model.matrix", {
  data <- makeTwoClass()
  measurements <- as.data.frame(data$measurements)
  colnames(measurements)[1] <- "unsafe name-1"
  expect_identical(ClassifyR:::.numericDesignMatrix(measurements), model.matrix(~ 0 + ., data = measurements))
})

test_that("an SVM fitted to a matrix of numeric features predicts as one fitted with a formula", {
  data <- makeTwoClass()
  measurements <- S4Vectors::DataFrame(data$measurements, check.names = FALSE)
  set.seed(1)
  fromMatrix <- ClassifyR:::SVMtrainInterface(measurements[1:40, ], data$classes[1:40], verbose = 0)
  set.seed(1)
  fromFormula <- e1071::svm(classes ~ ., data = data.frame(data$measurements[1:40, ], classes = data$classes[1:40]), probability = TRUE)
  expected <- predict(fromFormula, as.data.frame(data$measurements[41:60, ]), probability = TRUE)
  predicted <- ClassifyR:::SVMpredictInterface(fromMatrix, measurements[41:60, ], verbose = 0)
  expect_identical(as.character(predicted[, "class"]), as.character(expected))
  expect_identical(unname(as.matrix(predicted[, levels(data$classes)])), unname(attr(expected, "probabilities")[, levels(data$classes)]))
})
