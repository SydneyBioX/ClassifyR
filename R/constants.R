.ClassifyRenvir <- new.env(parent = emptyenv())

# Used internally during parameter selection based on best performance.
.ClassifyRenvir[["performanceInfoTable"]] <- matrix(c("Error", "lower",
                                                      "Accuracy", "higher",
                                                      "Balanced Error", "lower",
                                                      "Balanced Accuracy", "higher",
                                                      "Micro Precision", "higher",
                                                      "Micro Recall", "higher",
                                                      "Micro F1", "higher",
                                                      "Macro Precision", "higher",
                                                      "Macro Recall", "higher",
                                                      "Macro F1", "higher",
                                                      "Matthews Correlation Coefficient", "higher",
                                                      "AUC", "higher",
                                                      "C-index", "higher"),
                                                    ncol = 2, byrow = TRUE, dimnames = list(NULL, c("type", "better"))
) |> as.data.frame()

.ClassifyRenvir[["performanceTypes"]] <- .ClassifyRenvir[["performanceInfoTable"]][, "type"]

# Nice-looking names for feature selection and classification functions, to automatically use
# in a variety of performance plots.
.ClassifyRenvir[["functionsTable"]] <- matrix(
  c("subtractFromLocation", "Subtraction From Training Set Location",
    "bartlettRanking", "Bartlett Test",
    "classifyInterface", "Poisson LDA",
    "differentMeansRanking", "Difference in Means",
    "DLDAtrainInterface", "Diagonal LDA",
    "DMDranking", "Differences of Medians and Deviations",
    "edgeRranking", "edgeR LRT",
    "GLMtrainInterface", "Logistic Regression",
    "penalisedGLMtrainInterface", "Penalised GLMs (Ridge, Elastic net, LASSO)",
    "fisherDiscriminant", "Fisher's LDA",
    "kNNinterface", "k Nearest Neighbours",
    "KolmogorovSmirnovRanking", "Kolmogorov-Smirnov Test",
    "kTSPclassifier", "k Top-Scoring Pairs",
    "KullbackLeiblerRanking", "Kullback-Leibler Divergence",
    "leveneRanking", "Levene Test",
    "likelihoodRatioRanking", "Likelihood Ratio Test (Normal)",
    "limmaRanking", "Moderated t-test",
    "mixModelsTrain", "Mixtures of Normals",
    "naiveBayesKernel", "Naive Bayes Kernel",
    "NSCtrainInterface",  "Nearest Shrunken Centroids",
    "pairsDifferencesRanking", "Pairs Differences",
    "previousSelection", "Previous Selection", 
    "previousTrained", "Previous Trained",
    "randomSelection", "Random Selection",
    "randomForestTrainInterface", "Random Forest",
    "SVMtrainInterface", "Support Vector Machine",
    "coxphTrainInterface", "Cox Proportional Hazards",
    "coxphRanking", "Cox Proportional Hazards",
    "coxnetTrainInterface", "Penalised Cox Proportional Hazards",
    #"NEMOEtrainInterface", "Nutrition-Ecotype Mixture of Experts",
    "rfsrcTrainInterface", "Random Survival Forest",
    "extremeGradientBoostingTrainInterface", "Extreme Gradient Boosting"),
  ncol = 2, byrow = TRUE, dimnames = list(NULL, c("character", "name"))
) |> as.data.frame()

.ClassifyRenvir[["selectKeywords"]] <- matrix(
  c("none", "Skip selection procedure and use all input features.",
    "t-test", "T-test.",
    "limma", "Moderated t-test.",
    "edgeR", "edgeR likelihood ratio test.",
    "Bartlett", "Bartlett's test for different variance.",
    "Levene", "Levene's test for different variance.",
    "DMD", "Differences in means/medians and/or deviations.",
    "likelihoodRatio", "Likelihood ratio test (normal distribution).",
    "KS", "Kolmogorov-Smirnov test for differences in distributions.",
    "KL", "Kullback-Leibler divergence between distributions.",
    "CoxPH", "Cox proportional hazards Wald test per-feature.",
    "randomSelection", "Randomly selects a specified number of features."
    ),
  ncol = 2, byrow = TRUE, dimnames = list(NULL, c("selectionMethod Keyword", "Description"))
) |> as.data.frame()

.ClassifyRenvir[["classifyKeywords"]] <- matrix(
  c("randomForest", "Random forest.",
    "DLDA", "Diagonal Linear Discriminant Analysis.",
    "kNN", "k Nearest Neighbours.",
    "GLM", "Logistic regression.",
    "ridgeGLM", "Ridge GLM multinomial regression (alpha = 0).",
    "elasticNetGLM", "Elastic net GLM multinomial regression (alpha = 0.5).",
    "LASSOGLM", "LASSO GLM multinomial regression (alpha = 1).",  
    "SVM", "Support Vector Machine.",
    "NSC", "Nearest Shrunken Centroids.",
    "naiveBayes", "Naive Bayes kernel feature voting classifier.",
    "mixturesNormals", "Mixture of normals feature voting classifier.",
    "CoxPH", "Cox proportional hazards.",
    "CoxNet", "Penalised Cox proportional hazards.",
    "randomSurvivalForest", "Random survival forest.",
    "XGB", "Extreme gradient booster."
    ),
  ncol = 2, byrow = TRUE, dimnames = list(NULL, c("classifier Keyword", "Description"))
) |> as.data.frame()

.ClassifyRenvir[["multiViewKeywords"]] <- matrix(
  c("none", "Keep assays separate.",
    "merge", "Concatenate all selected feaures into single table before modelling.",
    "prevalidation", "Reduce each assay into a vector and concatenate to clinical data before modelling.",
    "PCA", "Reduce each assay into a lower dimensional representation and concatenate the principal components to the clinical data before modelling."
    ),
  ncol = 2, byrow = TRUE, dimnames = list(NULL, c("multiViewMethod Keyword", "Description"))
) |> as.data.frame()

.ClassifyRenvir[["prepareDataFormals"]] <- c("useFeatures", "maxMissingProp", "topNvariance")