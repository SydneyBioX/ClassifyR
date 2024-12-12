extractPrevalidation = function(assayPreval){ #}, startingCol) {
    # seq_len(length(assayPreval)) |>
    #     lapply(
    #         function(i)
    #             assayPreval[[i]] |>
    #             tibble::column_to_rownames("sample") |>
    #             dplyr::rename_all(function(col)
    #                 paste0(names(assayPreval)[i], col))
    #     ) |>
    #     # Taking all prevalidation class vectors except the last one (needed for multiple outcome)
    #     lapply(function(x)
    #         x[, c(startingCol:(ncol(x) - 1)), drop = FALSE] |> tibble::rownames_to_column("row_names")) |>
    #     purrr::reduce(merge, by = "row_names") |>
    #     janitor::clean_names() |>
    #     tibble::column_to_rownames("row_names")
    
    
    use <- which(names(assayPreval)!="clinical")
    
    assayPreval <- sapply(assayPreval[use], function(x){
        if(is.null(ncol(x)))x <- data.frame(sample = names(x), x)
        if(!"sample"%in%colnames(x))x <- data.frame(sample = rownames(x), x)
        x[order(x$sample),]}, simplify = FALSE)
    
    vec <- sapply(assayPreval, function(x){
        x <- x[,!colnames(x) %in% c("sample", "permutation", "fold", "class"), drop = FALSE]
        if(ncol(x)>1)return(as.matrix(x[,-1, drop = FALSE]))
        as.matrix(x)
    }, simplify = TRUE)
    rownames(vec) <- assayPreval[[1]]$sample
    vec
}

setClass("prevalModel", slots = "fullModel")

prevalTrainInterface <- function(measurements, outcome, params, verbose)
          {
              ###
              # Splitting measurements into a list of each of the assays
              ###
              assayTrain <- sapply(unique(S4Vectors::mcols(measurements)[["assay"]]), function(assay) measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assay], simplify = FALSE)
              
              if(!"clinical" %in% names(assayTrain)) stop("Must have an assay called \"clinical\"")
              
              tuneMode <- "none"
              performanceType <- "N/A"
              if(!is.null(params[[1]]@selectParams@tuneParams))
              {
                  tuneMode <- "Resubstitution"
                  if(is(outcomeTrain, "Surv")) performanceType <- "C-index" else performanceType <- "Balanced Accuracy"
              }
              
              crossValParams <- CrossValParams(permutations = 1, folds = 10, parallelParams = SerialParam(RNGseed = .Random.seed[1]), tuneMode = tuneMode, performanceType = performanceType)
              if(is(outcome, "Surv")) crossValParams@performanceType <- "C-index" else crossValParams@performanceType <- "Balanced Accuracy"
              ###
              # Fit a classification model for each non-clinical data set, pulling models from "params"
              ###
              usePreval <- names(assayTrain)[names(assayTrain) != "clinical"]
              assayTests <- mapply(
                  runTests,
                  measurements = assayTrain[usePreval],
                  modellingParams = params[usePreval],
                  MoreArgs = list(
                      outcome = outcome,
                      crossValParams = crossValParams,
                      verbose = 0
                  )) |> sapply(function(result) result@predictions, simplify = FALSE)
              
              ###
              # Pull-out prevalidated vectors ie. the predictions on each of the test folds.
              ###
              prevalidationTrain <- extractPrevalidation(assayTests)
              
              # Feature select on clinical data before binding
              # selectedFeaturesClinical <- runTest(assayTrain[["clinical"]],
              #                                    outcome = outcome,
              #                                    training = seq_len(nrow(assayTrain[["clinical"]])),
              #                                    testing = seq_len(nrow(assayTrain[["clinical"]])),
              #                                    modellingParams = params[["clinical"]],
              #                                    crossValParams = CVparams,
              #                                    .iteration = 1,
              #                                    verbose = 0
              # )$selected[, "feature"]
            
              #fullTrain = cbind(assayTrain[["clinical"]][,selectedFeaturesClinical], prevalidationTrain[rownames(assayTrain[["clinical"]]), , drop = FALSE])
              
              prevalidationTrain <- S4Vectors::DataFrame(prevalidationTrain, check.names = FALSE)
              S4Vectors::mcols(prevalidationTrain)$assay = "prevalidation"
              S4Vectors::mcols(prevalidationTrain)$feature = colnames(prevalidationTrain)
              
              ###
              # Bind the prevalidated data to the clinical data
              ###
              fullTrain = cbind(assayTrain[["clinical"]], prevalidationTrain[rownames(assayTrain[["clinical"]]), , drop = FALSE])

              # Pull out clinical data
              finalModParam <- params[["clinical"]]
              #finalModParam@selectParams <- NULL
              
              # Fit classification model (from clinical in params)
              runTestOutput = runTest(
                  measurementsTrain = fullTrain,
                  outcomeTrain = outcome,
                  measurementsTest = fullTrain,
                  outcomeTest = outcome,
                  modellingParams = finalModParam,
                  crossValParams = crossValParams,
                  .iteration = 1,
                  verbose = 0
                  )
              
              
              # Extract the classification model from runTest output
              fullModel = runTestOutput$models
              fullModel$fullFeatures = colnames(fullTrain)
              
              # Fit models with each datatype for use in prevalidated prediction later..
              prevalidationModels =  mapply(
                  runTest,
                  measurementsTrain = assayTrain,
                  measurementsTest = assayTrain,               
                  modellingParams = params,
                  MoreArgs = list(
                      crossValParams = crossValParams,
                      outcomeTrain = outcome,
                      outcomeTest = outcome,
                      .iteration = 1,
                      verbose = 0
                  )
              )
              
              # Add prevalidated models and classification params for each datatype to the fullModel object
              fullModel$prevalidationModels <- prevalidationModels["models",]
              names(fullModel$prevalidationModels) <- colnames(prevalidationModels)
              fullModel$modellingParams <- params
              fullModel$prevalFeatures <- prevalidationModels["selected",]
              fullModel$prevalFeaturesRanked <- runTestOutput$ranked$feature
              fullModel$prevalFeaturesSelected <- runTestOutput$selected$feature
              
              fullModel <- new("prevalModel", fullModel = fullModel)
              fullModel
}

prevalFeatures <- function(prevalModel)
                  {
                    list(prevalModel@fullModel$prevalFeaturesRanked, prevalModel@fullModel$prevalFeaturesSelected)
                  }

prevalPredictInterface <- function(fullModel, test, returnType = "both", verbose)
          {
              fullModel <- fullModel@fullModel
              assayTest <- sapply(unique(S4Vectors::mcols(test)[["assay"]]), function(assay) test[, S4Vectors::mcols(test)[["assay"]] %in% assay], simplify = FALSE)
              
              prevalidationModels <- fullModel$prevalidationModels
              modelPredictionFunctions <- fullModel$modellingParams
              
              prevalidationPredict <- sapply(names(prevalidationModels), function(x){
                  predictParams <- modelPredictionFunctions[[x]]@predictParams
                  paramList <- list(prevalidationModels[[x]], assayTest[[x]])
                  if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
                  paramList <- c(paramList, verbose = 0)
                  prediction <- do.call(predictParams@predictor, paramList)
                  prediction}, simplify = FALSE) |>
                  extractPrevalidation()
              
              prevalidationPredict <- S4Vectors::DataFrame(prevalidationPredict)
              S4Vectors::mcols(prevalidationPredict)$assay = "prevalidation"
              S4Vectors::mcols(prevalidationPredict)$feature = colnames(prevalidationPredict)
              
              fullTest = cbind(assayTest[["clinical"]], prevalidationPredict[rownames(assayTest[["clinical"]]), , drop = FALSE])
              
              
              predictParams <- modelPredictionFunctions[["clinical"]]@predictParams
              paramList <- list(fullModel,  fullTest)
              if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
              paramList <- c(paramList, verbose = 0)
              finalPredictions <- do.call(predictParams@predictor, paramList)

              finalPredictions
          }