################################################################################
#
# The contract of a classifier or feature selection interface. Copy this outline for a new method, then register it
# in registry.R.
#
# Training: newTrainInterface(measurementsTrain, outcomeTrain, <parameters>, ..., verbose = 3)
# - measurementsTrain is a DataFrame of the training samples, with features in columns. Categorical features are
#   factors; methods needing a numeric matrix encode them with .encodeTrain and, at prediction, .encodeTest.
# - outcomeTrain (named classesTrain by the classifiers of classes only) is a factor of classes or a Surv object, as
#   the registry entry's outcomes allow. It is passed by position.
# - <parameters> are the method's settings. Those in the entry's tunePresets are tuned when tuneParams is "auto".
# - It returns the trained model.
#
# Prediction: newPredictInterface(model, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
# - Classes: "class" gives a factor with the training levels, "score" a samples x classes numeric matrix of class
#   scores or probabilities with the class names as column names, and "both" a data frame of a column "class"
#   followed by the score columns.
# - Survival: a numeric vector of risk scores, where higher means a higher risk of the event.
# - A classifier that trains and predicts in one function takes measurementsTest as its third argument and has no
#   prediction interface (predict = NULL in its registry entry).
#
# Features used by a model (optional), for classifiers that rank or select features themselves: newFeatures(model)
# returns a list of two vectors of column indices of measurementsTrain, the features ranked from best to worst and
# the features the model uses. Other arguments, if any, are named after variables of runTest (for example,
# measurementsTrain), which supplies them. It is the entry's getFeatures.
#
# Feature ranking: newRanking(measurementsTrain, outcomeTrain, ..., verbose = 3) returns the indices of the columns
# of measurementsTrain, from the best feature to the worst.
#
# Each function carries its name as an attribute, which results and plots use to show its display name:
#   attr(newTrainInterface, "name") <- "newTrainInterface"
#
################################################################################
