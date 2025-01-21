selectMulti <- function(measurementsTrain, outcomeTrain, params, verbose = 0)
          {                 
              assaysIndices <- lapply(unique(S4Vectors::mcols(measurementsTrain)[["assay"]]), function(assay) which(S4Vectors::mcols(measurementsTrain)[["assay"]] == assay))
              assayTrain <- lapply(assaysIndices, function(assayIndices) measurementsTrain[, assayIndices, drop = FALSE])
              tuneMode <- "none"
              performanceType <- "N/A"
              if(!is.null(params[[1]]@selectParams@tuneParams))
              {
                  tuneMode <- "Resubstitution"
                  if(is(outcomeTrain, "Surv")) performanceType <- "C-index" else performanceType <- "Balanced Accuracy"
              }
              featuresIndices <- mapply(.doSelection, 
                                         measurements = assayTrain,
                                         modellingParams = params,
                                         MoreArgs = list(outcomeTrain = outcomeTrain, 
                                                         crossValParams = CrossValParams(permutations = 1, folds = 5, tuneMode = tuneMode, performanceType = performanceType), ###### Where to get this from?
                                                         verbose = 0), SIMPLIFY = FALSE
                                        )

              unlist(mapply(function(allDataIndices, withinIndices) allDataIndices[withinIndices],
                     assaysIndices, lapply(featuresIndices, "[[", 2), SIMPLIFY = FALSE))
}
attr(selectMulti, "name") <- "Union Selection"
