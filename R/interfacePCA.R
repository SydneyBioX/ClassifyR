setClass("pcaModel", slots = "fullModel")

pcaTrainInterface <- function(measurements, outcome, params, nFeatures, ...)
          {
              
              ###
              # Splitting measurements into a list of each of the datasets
              ###
              assayTrain <- sapply(unique(S4Vectors::mcols(measurements)[["assay"]]), function(assay) measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assay], simplify = FALSE)
              
              if(!"clinical" %in% names(assayTrain)) stop("Must have an assay called \"clinical\".")
              
              # Create generic crossValParams just to get things working, might be used for optimising features in runTest later???
              CVparams <- CrossValParams(permutations = 1, folds = 10, parallelParams = SerialParam(RNGseed = .Random.seed[1]), tuneMode = "Resubstitution") 
              if(is(outcome, "Surv")) CVparams@performanceType <- "C-index" else CVparams@performanceType <- "Balanced Accuracy"
              ###
              # Run PCA for all assays except clinical
              ###
              usePCA<- names(assayTrain)[names(assayTrain)!="clinical"]
              assayPCA <- sapply(assayTrain[usePCA], function(assay){
                  pca <- prcomp(assay)
              }, simplify = FALSE)
              
              
              ###
              # Pull out the top n (nFeatures) PCA loadings
              ###
              
              pcaVar <- mapply(function(pca, n){
                  data.frame(pca$x[,seq_len(n)])
              }, assayPCA, nFeatures[names(assayPCA)], SIMPLIFY = FALSE)
              
              pcaVar <- do.call("cbind", pcaVar)
              
              ###
              # Bind the PCA components to the clinical data
              ###
              
              pcaVar <- S4Vectors::DataFrame(pcaVar)
              S4Vectors::mcols(pcaVar)$assay = "PCA"
              S4Vectors::mcols(pcaVar)$feature = colnames(pcaVar)
            
              fullTrain = cbind(assayTrain[["clinical"]], pcaVar)
              
              
              ###
              # Pull out the modellingParams for the clinical data.
              ###
              
              finalModParam <- params[["clinical"]]
              
              
              ###
              # Fit model with runTest
              ### 
              
              runTestOutput = runTest(
                  fullTrain,
                  outcome,
                  fullTrain,
                  outcome,
                  modellingParams = finalModParam,
                  crossValParams = CVparams,
                  .iteration = 1,
                  verbose = 0
                  )
              
              ###
              # Pull out fitted model from runTest
              ###
              
              fullModel = runTestOutput$models
              fullModel$fullFeatures = colnames(fullTrain)
              
              # Make a class
              fullModel$pcaModels <- assayPCA # Add PCA "models" to final model.
              fullModel$nFeatures <- nFeatures
              fullModel$modellingParams <- finalModParam
              fullModel$PCAfeaturesRanked <- runTestOutput$ranked$feature
              fullModel$PCAfeaturesSelected <- runTestOutput$selected$feature
              fullModel <- new("pcaModel", fullModel = fullModel)
              fullModel
}

PCAfeatures <- function(PCAmodel)
                  {
                    list(PCAmodel@fullModel$PCAfeaturesRanked, PCAmodel@fullModel$PCAfeaturesSelected)
                  }


pcaPredictInterface <- function(fullModel, test, ..., returnType = "both", verbose = 0)
          {
              # Pull out my classification model
              fullModel <- fullModel@fullModel
              
              #Split my test data into a list of the different assays
              assayTest <- sapply(unique(S4Vectors::mcols(test)[["assay"]]), function(assay) test[, S4Vectors::mcols(test)[["assay"]] %in% assay], simplify = FALSE)
              
              # Pull out my PCA models
              pcaModels <- fullModel$pcaModels
              nFeatures <- fullModel$nFeatures
              
              # Project test data into trained PCA space
              pcaVar <- mapply(function(pca, assay, n){
                 data.frame(predict(pca, assay))[,seq_len(n)]
              }, pcaModels, assayTest[names(pcaModels)], nFeatures[names(pcaModels)], SIMPLIFY = FALSE)
              
              pcaVar <- do.call(cbind, pcaVar)
              
              pcaVar <- S4Vectors::DataFrame(pcaVar)
              S4Vectors::mcols(pcaVar)$assay = "PCA"
              S4Vectors::mcols(pcaVar)$feature = colnames(pcaVar)
              
              # Merge my PCA stuff with my clinical data
              fullTest = cbind(assayTest[["clinical"]], pcaVar)

              # Pull out my predictParams for clinical
              modelPredictionFunctions <- fullModel$modellingParams
              predictParams <- modelPredictionFunctions@predictParams
              paramList <- list(fullModel,  fullTest)
              if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
              paramList <- c(paramList, verbose = 0)
              
              # Predict with my predictor
              finalPredictions <- do.call(predictParams@predictor, paramList)

              finalPredictions
          }