# Random Forest
RFparams <- function(tuneParams) {
    if(is.character(tuneParams) && tuneParams == "auto") tuneParams <- list(mTryProportion = c(0.10, 0.25, 0.33, 0.5), num.trees = c(1, 10, 100))
    trainParams <- TrainParams(randomForestTrainInterface, tuneParams = tuneParams,
                               getFeatures = forestFeatures)
    predictParams <- PredictParams(randomForestPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Random Survival Forest
RSFparams <- function(tuneParams) {
    if(is.character(tuneParams) && tuneParams == "auto") tuneParams <- list(mTryProportion = c(0.10, 0.25, 0.33, 0.5), ntree = c(1, 10, 100))
    trainParams <- TrainParams(rfsrcTrainInterface, tuneParams = tuneParams,
                               getFeatures = rfsrcFeatures)
    predictParams <- PredictParams(rfsrcPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

XGBparams <- function(tuneParams) {
    if(is.character(tuneParams) && tuneParams == "auto") tuneParams <- list(mTryProportion = c(0.10, 0.25, 0.33, 0.5), nrounds = c(5, 10))
    trainParams <- TrainParams(extremeGradientBoostingTrainInterface, tuneParams = tuneParams,
                               getFeatures = XGBfeatures)
    predictParams <- PredictParams(extremeGradientBoostingPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))  
}

# k Nearest Neighbours
kNNparams <- function() {
    trainParams <- TrainParams(kNNinterface)
    predictParams <- NULL
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Ordinary GLM
GLMparams <- function() {
    trainParams <- TrainParams(GLMtrainInterface)
    predictParams <- PredictParams(GLMpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Ridge GLM
ridgeGLMparams <- function() {
    trainParams <- TrainParams(penalisedGLMtrainInterface, alpha = 0, getFeatures = penalisedFeatures)
    predictParams <- PredictParams(penalisedGLMpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Elastic net GLM
elasticNetGLMparams <- function() {
    trainParams <- TrainParams(penalisedGLMtrainInterface, alpha = 0.5, getFeatures = penalisedFeatures)
    predictParams <- PredictParams(penalisedGLMpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# LASSO GLM
LASSOGLMparams <- function() {
    trainParams <- TrainParams(penalisedGLMtrainInterface, getFeatures = penalisedFeatures)
    predictParams <- PredictParams(penalisedGLMpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Support Vector Machine
SVMparams = function(tuneParams) {
    if(is.character(tuneParams) && tuneParams == "auto")
        tuneParams <- list(kernel = c("linear", "polynomial", "radial", "sigmoid"), cost = 10^(-3:3))
    trainParams <- TrainParams(SVMtrainInterface, tuneParams = tuneParams)
    predictParams <- PredictParams(SVMpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Nearest Shrunken Centroid
NSCparams = function() {
    trainParams <- TrainParams(NSCtrainInterface, getFeatures = NSCfeatures)
    predictParams <- PredictParams(NSCpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Diagonal Linear Discriminant Analysis
DLDAparams = function() {
    trainParams <- TrainParams(DLDAtrainInterface)
    predictParams <- PredictParams(DLDApredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# naive Bayes Kernel
naiveBayesParams <- function(tuneParams) {
    if(is.character(tuneParams) && tuneParams == "auto") tuneParams <- list(difference = c("unweighted", "weighted")) 
    trainParams <- TrainParams(naiveBayesKernel, tuneParams = tuneParams)
    predictParams <- NULL
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Mixtures of Normals
mixModelsParams <- function() {
    trainParams <- TrainParams(mixModelsTrain, nbCluster = 1:2)
    predictParams <- PredictParams(mixModelsPredict)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Cox Proportional Hazards Model for Survival
coxphParams <- function() {
    trainParams <- TrainParams(coxphTrainInterface)
    predictParams <- PredictParams(predictor = coxphPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Cox Proportional Hazards Model with Elastic Net for Survival
coxnetParams <- function() {
    trainParams <- TrainParams(coxnetTrainInterface, getFeatures = penalisedFeatures)
    predictParams <- PredictParams(coxnetPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}