# Small data sets shared by the tests. Kept small so the whole suite runs in about a minute.
makeTwoClass <- function(nSamples = 60, nFeatures = 30, shift = 1.5, seed = 1)
{
  set.seed(seed)
  classes <- factor(rep(c("A", "B"), length.out = nSamples))
  measurements <- matrix(rnorm(nSamples * nFeatures), nSamples, nFeatures,
                         dimnames = list(paste0("s", seq_len(nSamples)), paste0("g", seq_len(nFeatures))))
  measurements[classes == "B", 1:3] <- measurements[classes == "B", 1:3] + shift
  list(measurements = measurements, classes = classes)
}

makeSurvival <- function(nSamples = 80, nFeatures = 20, seed = 2)
{
  set.seed(seed)
  measurements <- matrix(rnorm(nSamples * nFeatures), nSamples, nFeatures,
                         dimnames = list(paste0("s", seq_len(nSamples)), paste0("g", seq_len(nFeatures))))
  time <- rexp(nSamples, rate = exp(0.8 * measurements[, 1]))
  censor <- rexp(nSamples, rate = 0.3)
  list(measurements = measurements, outcome = survival::Surv(pmin(time, censor), as.integer(time <= censor)))
}
