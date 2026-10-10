################################################################################
#
# The registry of classifiers and feature selection methods.
#
# Adding a classifier or a feature selection method needs its interface file (see interfaceTemplate.R), its entry
# in Collate and one entry below. Everything else is read from the registry: the keywords of crossValidate,
# TrainParams, PredictParams and SelectParams, the tables of available(), the names shown in results and plots,
# tuning presets, the outcomes a method accepts and the fallback for too few features.
#
# The registry is built when the package is loaded, so the order of files in Collate does not matter.
#
################################################################################

# Registers a classifier.
# - keyword: the classifier keyword of crossValidate, TrainParams and PredictParams.
# - train, predict: the training and prediction functions. predict is NULL if train also predicts.
# - getFeatures: the function that extracts the features a trained model used, if the classifier selects features.
# - displayName: the name shown in results and plots. Classifiers sharing a training function share it.
# - description: one sentence for available().
# - outcomes: "classes", "survival" or both.
# - tunePresets: the ranges of parameters tuned when tuneParams is "auto". NULL if the classifier has none.
# - defaults: other parameters of the training function.
# - minFeatures, fewFeaturesClassifier: the classifier used instead when an assay has fewer features than minFeatures.
# - listed: whether available() lists it.
.registerClassifier <- function(keyword, train, predict = NULL, getFeatures = NULL, displayName, description,
                                outcomes = "classes", tunePresets = NULL, defaults = list(), minFeatures = 1,
                                fewFeaturesClassifier = NULL, listed = TRUE)
{
  .ClassifyRenvir[["classifiers"]][[keyword]] <- list(train = train, predict = predict, getFeatures = getFeatures,
    displayName = displayName, description = description, outcomes = outcomes, tunePresets = tunePresets,
    defaults = defaults, minFeatures = minFeatures, fewFeaturesClassifier = fewFeaturesClassifier, listed = listed)
}

# Registers a feature selection method. rank is the ranking (or selection) function, NULL for no selection.
.registerSelection <- function(keyword, rank, displayName = NULL, description, outcomes = c("classes", "survival"),
                               listed = TRUE)
{
  .ClassifyRenvir[["selections"]][[keyword]] <- list(rank = rank, displayName = displayName, description = description,
                                                    outcomes = outcomes, listed = listed)
}

.buildRegistry <- function()
{
  .ClassifyRenvir[["classifiers"]] <- list()
  .ClassifyRenvir[["selections"]] <- list()
  forestPresets <- list(mTryProportion = c(0.10, 0.25, 0.33, 0.5), num.trees = c(1, 10, 100))

  ##### Classifiers, in the order available() lists them #####
  .registerClassifier("randomForest", randomForestTrainInterface, randomForestPredictInterface, forestFeatures,
                      "Random Forest", "Random forest.", c("classes", "survival"), forestPresets)
  .registerClassifier("DLDA", DLDAtrainInterface, DLDApredictInterface, displayName = "Diagonal LDA",
                      description = "Diagonal Linear Discriminant Analysis.")
  .registerClassifier("kNN", kNNinterface, displayName = "k Nearest Neighbours", description = "k Nearest Neighbours.",
                      tunePresets = list(k = 1:5))
  .registerClassifier("GLM", GLMtrainInterface, GLMpredictInterface, displayName = "Logistic Regression",
                      description = "Logistic regression.")
  penalisedName <- "Penalised GLMs (Ridge, Elastic net, LASSO)"
  .registerClassifier("ridgeGLM", penalisedGLMtrainInterface, penalisedGLMpredictInterface, penalisedFeatures,
                      penalisedName, "Ridge GLM multinomial regression (alpha = 0).", defaults = list(alpha = 0),
                      minFeatures = 2, fewFeaturesClassifier = "GLM")
  .registerClassifier("elasticNetGLM", penalisedGLMtrainInterface, penalisedGLMpredictInterface, penalisedFeatures,
                      penalisedName, "Elastic net GLM multinomial regression (alpha = 0.5).", defaults = list(alpha = 0.5),
                      minFeatures = 2, fewFeaturesClassifier = "GLM")
  .registerClassifier("LASSOGLM", penalisedGLMtrainInterface, penalisedGLMpredictInterface, penalisedFeatures,
                      penalisedName, "LASSO GLM multinomial regression (alpha = 1).",
                      minFeatures = 2, fewFeaturesClassifier = "GLM")
  .registerClassifier("SVM", SVMtrainInterface, SVMpredictInterface, displayName = "Support Vector Machine",
                      description = "Support Vector Machine.",
                      tunePresets = list(kernel = c("linear", "polynomial", "radial", "sigmoid"), cost = 10^(-3:3)))
  .registerClassifier("NSC", NSCtrainInterface, NSCpredictInterface, NSCfeatures, "Nearest Shrunken Centroids",
                      "Nearest Shrunken Centroids.")
  .registerClassifier("naiveBayes", naiveBayesKernel, displayName = "Naive Bayes Kernel",
                      description = "Naive Bayes kernel feature voting classifier.",
                      tunePresets = list(difference = c("unweighted", "weighted")))
  .registerClassifier("mixturesNormals", mixModelsTrain, mixModelsPredict, displayName = "Mixtures of Normals",
                      description = "Mixture of normals feature voting classifier.", defaults = list(nbCluster = 1:2))
  .registerClassifier("CoxPH", coxphTrainInterface, coxphPredictInterface, displayName = "Cox Proportional Hazards",
                      description = "Cox proportional hazards.", outcomes = "survival")
  .registerClassifier("CoxNet", coxnetTrainInterface, coxnetPredictInterface, penalisedFeatures,
                      "Penalised Cox Proportional Hazards", "Penalised Cox proportional hazards.", outcomes = "survival")
  .registerClassifier("randomSurvivalForest", rfsrcTrainInterface, rfsrcPredictInterface, rfsrcFeatures,
                      "Random Survival Forest", "Random survival forest.", outcomes = "survival",
                      tunePresets = list(mTryProportion = c(0.10, 0.25, 0.33, 0.5), ntree = c(1, 10, 100)))
  .registerClassifier("XGB", extremeGradientBoostingTrainInterface, extremeGradientBoostingPredictInterface,
                      XGBfeatures, "Extreme Gradient Boosting", "Extreme gradient booster.", c("classes", "survival"),
                      list(mTryProportion = c(0.10, 0.25, 0.33, 0.5), nrounds = c(5, 10)))
  .registerClassifier("aorsf", obliqueForestTrainInterface, obliqueForestPredictInterface, obliqueForestFeatures,
                      "Oblique Random Forest", "Oblique random forest (aorsf).", c("classes", "survival"),
                      list(mTryProportion = c(0.10, 0.25, 0.33, 0.5)))
  .registerClassifier("glmboost", boostedCoxTrainInterface, boostedCoxPredictInterface, boostedCoxFeatures,
                      "Boosted Cox Proportional Hazards",
                      "Componentwise boosted Cox proportional hazards (mboost's glmboost).", "survival",
                      list(mstop = c(50, 100, 200)))
  .registerClassifier("LiblineaR", LiblineaRtrainInterface, LiblineaRpredictInterface, LiblineaRfeatures,
                      "LIBLINEAR Logistic Regression",
                      "Logistic regression with an L2 or L1 penalty (LiblineaR), cost chosen by cross-validation.",
                      tunePresets = list(type = c(0, 6)))
  .registerClassifier("ncvreg", ncvregTrainInterface, ncvregPredictInterface, ncvregFeatures,
                      "MCP Penalised Regression",
                      "Logistic (two classes) or Cox regression with an MCP penalty (ncvreg), lambda chosen by cross-validation.",
                      c("classes", "survival"), list(penalty = c("MCP", "SCAD", "lasso")))
  # Uses the models trained in the same iteration of a previous cross-validation.
  .registerClassifier("previousTrained", previousTrained, displayName = "Previous Trained",
                      description = "The models of a previous cross-validation.", outcomes = c("classes", "survival"),
                      listed = FALSE)

  ##### Feature selection, in the order available() lists them #####
  .registerSelection("none", NULL, description = "Skip selection procedure and use all input features.")
  .registerSelection("t-test", differentMeansRanking, "Difference in Means", "T-test.", "classes")
  .registerSelection("limma", limmaRanking, "Moderated t-test", "Moderated t-test.", "classes")
  .registerSelection("edgeR", edgeRranking, "edgeR LRT", "edgeR likelihood ratio test.", "classes")
  .registerSelection("Bartlett", bartlettRanking, "Bartlett Test", "Bartlett's test for different variance.", "classes")
  .registerSelection("Levene", leveneRanking, "Levene Test", "Levene's test for different variance.", "classes")
  .registerSelection("DMD", DMDranking, "Differences of Medians and Deviations",
                     "Differences in means/medians and/or deviations.", "classes")
  .registerSelection("likelihoodRatio", likelihoodRatioRanking, "Likelihood Ratio Test (Normal)",
                     "Likelihood ratio test (normal distribution).", "classes")
  .registerSelection("KS", KolmogorovSmirnovRanking, "Kolmogorov-Smirnov Test",
                     "Kolmogorov-Smirnov test for differences in distributions.", "classes")
  .registerSelection("KL", KullbackLeiblerRanking, "Kullback-Leibler Divergence",
                     "Kullback-Leibler divergence between distributions.", "classes")
  .registerSelection("CoxPH", coxphRanking, "Cox Proportional Hazards", "Cox proportional hazards Wald test per-feature.",
                     "survival")
  .registerSelection("randomSelection", randomSelection, "Random Selection",
                     "Randomly selects a specified number of features.")
  .registerSelection("previousSelection", previousSelection, "Previous Selection",
                     "The features chosen in the same iteration of a previous cross-validation.", listed = FALSE)
  .registerSelection("selectMulti", selectMulti, description = "Union of the selections of several assays.",
                     listed = FALSE)

  .ClassifyRenvir[["classifyKeywords"]] <- .keywordsTable(.ClassifyRenvir[["classifiers"]], "classifier Keyword")
  .ClassifyRenvir[["selectKeywords"]] <- .keywordsTable(.ClassifyRenvir[["selections"]], "selectionMethod Keyword")

  # Names shown in results and plots, by the name of a function. Functions without a keyword are listed here too.
  registered <- c(lapply(.ClassifyRenvir[["classifiers"]], function(entry) c(attr(entry[["train"]], "name"), entry[["displayName"]])),
                  lapply(.ClassifyRenvir[["selections"]], function(entry) c(attr(entry[["rank"]], "name"), entry[["displayName"]])))
  registered <- do.call(rbind, Filter(function(pair) length(pair) == 2, registered))
  others <- matrix(c("subtractFromLocation", "Subtraction From Training Set Location",
                     "classifyInterface", "Poisson LDA",
                     "fisherDiscriminant", "Fisher's LDA",
                     "kTSPclassifier", "k Top-Scoring Pairs",
                     "pairsDifferencesRanking", "Pairs Differences"), ncol = 2, byrow = TRUE)
  functionsTable <- unique(rbind(registered, others))
  dimnames(functionsTable) <- list(NULL, c("character", "name"))
  .ClassifyRenvir[["functionsTable"]] <- as.data.frame(functionsTable)
}

# The table of listed keywords that available() shows.
.keywordsTable <- function(registry, keywordColumn)
{
  listed <- Filter(function(entry) entry[["listed"]], registry)
  table <- data.frame(names(listed), sapply(listed, `[[`, "description"),
                      sapply(listed, function(entry) paste(entry[["outcomes"]], collapse = ", ")), row.names = NULL)
  setNames(table, c(keywordColumn, "Description", "Outcomes"))
}

.selectionKeywordToFunction <- function(keyword)
{
  .ClassifyRenvir[["selections"]][[keyword]][["rank"]]
}

# A classifier keyword's TrainParams and PredictParams (NULL if training also predicts), or NULL for an unknown keyword.
# tuneParams is used only by classifiers with tuning presets. "auto" chooses the presets.
.classifierKeywordToParams <- function(keyword, tuneParams)
{
  entry <- .ClassifyRenvir[["classifiers"]][[keyword]]
  if(is.null(entry)) return(NULL)
  if(is.null(entry[["tunePresets"]])) tuneParams <- NULL
  else if(is.character(tuneParams) && tuneParams == "auto") tuneParams <- entry[["tunePresets"]]
  trainParams <- do.call(TrainParams, c(list(entry[["train"]], tuneParams = tuneParams, getFeatures = entry[["getFeatures"]]),
                                        entry[["defaults"]]))
  predictParams <- if(!is.null(entry[["predict"]])) PredictParams(entry[["predict"]])
  list(trainParams = trainParams, predictParams = predictParams)
}

# Stops if a classifier or feature selection keyword does not accept the kind of outcome.
.checkOutcomeType <- function(keywords, outcome, registryName)
{
  outcomeType <- if(is(outcome, "Surv")) "survival" else "classes"
  registry <- .ClassifyRenvir[[registryName]]
  unsuitable <- Filter(function(keyword) !is.null(registry[[keyword]]) && !outcomeType %in% registry[[keyword]][["outcomes"]],
                       unique(unlist(keywords)))
  if(length(unsuitable) > 0)
    stop(paste(unsuitable, collapse = ", "), " can't be used with ", if(outcomeType == "survival") "a survival outcome" else "classes",
         ". See the Outcomes column of available(", if(registryName == "selections") "\"selectionMethod\"", ").")
}

.onLoad <- function(libname, pkgname)
{
  .buildRegistry()
}
