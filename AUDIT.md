# ClassifyR audit (October 2026)

Audit of ClassifyR 3.17.1 (`devel` at 378c217), carried out on branch `audit-2026-10`. It has two aims:

1. find ways to make the cross-validation machinery faster;
2. list quirks: wrong results, broken code paths and surprising behaviour.

**How it was done.** The core pipeline was read line by line:
- `crossValidate`, `runTests`, `runTest`, `utilities`, `prepareData`, `simpleParams`;
- the interfaces and rankings;
- metrics and plots;
- precision pathways, crissCrossValidate and the classes.

Suspected bugs were confirmed by running them. Timings come from the protocol's own data (METABRIC views, Procedure 2 of the draft paper) and from the package's `asthma` data (190 × 2000). Runs were on verona, single-threaded BLAS, R 4.6.1. File:line references are to this commit.

Each finding is marked:

- **[confirmed]**: reproduced by running it;
- **[read]**: from reading the code only.

---

## Part 1. Speed

### 1A. The cross-validation machinery

These change *how* the work is organised, not what is computed. Every item below should give identical
results for a fixed seed.

**M1. Subsetting `DataFrame`s is half the run time for cheap models.**
- **Profile [confirmed].** Asthma, t-test + DLDA, 20 × 5 CV, 11.2 s:
  - 52% of the time is S4 `[` (`NSBS`, `normalizeSingleBracketSubscript`);
  - feature selection is 30%;
  - training plus prediction is 5%.
- **Cost of one subset.** Taking 152 rows × 2000 columns costs **27.6 ms from a `DataFrame` against 1.1 ms from a matrix**.
- **Where the subsets happen.** Each fold subsets the data at least five times:
  - train and test rows in `runTests.R:126`;
  - selected columns in `runTest.R`;
  - `featuresForTrain` in `utilities.R:377`;
  - then most interfaces convert to `data.frame` or `matrix` again.

*Fix:* keep the numeric assays as one base `matrix` (feature metadata in a side table), and subset by integer index.
Only clinical/mixed-type columns need a `data.frame`; encode them once (`model.matrix` on the full data,
which is unsupervised) or keep them as a small separate block. Expected gain: about 2× overall for cheap
models, and less for expensive ones.

**M2. Parallelism is spread over too few tasks, so it scales poorly.**
- **How it is split now.** `crossValidate` loops *sequentially* over assays × classifiers × selection methods (`crossValidate.R:181-212`), and over every assay combination for merge, prevalidation and PCA. Each `CV()` call runs its own `bpmapply` over that one run's repeats × folds.
- **Measured [confirmed]**, asthma, DLDA, 20 × 5:

  | Run | 1 core | 8 cores | Speed-up |
  |---|---|---|---|
  | one assay | 11.0 s | 4.4 s | 2.5× (4 cores: 4.9 s) |
  | three assays, each on its own | 13.4 s | 8.9 s | 1.5× |
  | merge, 7 combinations | 46.3 s | 24.8 s | 1.9× |

*Fix:* build the full task list first, then run one `bplapply` over it, with workers started once. A task is
one (assay or combination, classifier, selection method, repeat, fold). Combined with M3 and M4, this should get
close to linear scaling.

**M3. Fixed costs paid on every `CV()` call.** These are paid in the parent process, outside the parallel section. Profile [confirmed], asthma run with 8 cores:

| Cost | Time | Where |
|---|---|---|
| Building a `MulticoreParam` (`.snowCoresMax` → `showConnections`) | 0.8 s per call | `generateCrossValParams`, `crossValidate.R:483` |
| `any(is.na(measurements))` on a `DataFrame` | 0.9 s | `runTests.R:82` |
| `prepareData` repeated: in the MAE/list method, then again inside `runTests` (`runTests.R:88`) | — | the constant-column check `apply(measurements, 2, ...)` (`prepareData.R:100`) turns the whole `DataFrame` into a matrix each time |
| First-use namespace loading (genefilter) | 1.3 s once | — |

The merge pays the per-call costs once for every combination. The paper's Procedure 2 merge (5 assays) has 31 combinations, and k assays give 2^k − 1. *Fix:* build one BiocParallel param per `crossValidate` call; check for NAs with `anyNA` on the matrix; prepare the data once.

**M4. The final full-data model runs serially, and also inside nested loops.**
- `runTests.R:193` refits on all samples after the parallel section. That is serial time on every call.
- It is also paid inside every *inner* `runTests` of nested-CV tuning, together with building a whole `ClassifyResult` that is then thrown away.
- *Fix:* send the full fit to the workers as one more task; skip it (and `ClassifyResult` construction) when `runTests` is called internally for tuning.

**M5. Merge, prevalidation and PCA redo the per-assay feature selection in every combination.**
- Each combination is cross-validated on its own folds (Q-C9). With shared folds, assay A's ranking in repeat r, fold f would be the same in every combination that contains A, and could be computed once instead of 2^(k−1) times.
- The single-assay combinations are also the same runs as `multiViewMethod = "none"`.
- *Fix:* compute each (assay, repeat, fold) ranking once and cache it. Then each combination only trains (merge), or trains the meta-model (prevalidation). In the DLDA run selection was 30% of the time; with CoxPH ranking plus CoxNet it is 4%, so the gain depends on the model.

**M6. Nested-CV tuning is far heavier than it needs to be.**
- **Same scheme inside and out.** The inner loop reuses the outer `CrossValParams` (`utilities.R:222`, `:333`). With the default 20 × 5, each tuning value costs 100 inner fits per outer fold, which is 10,000 fits per tuning value overall.
- **Nested parallelism.** The inner loop also reuses the outer `MulticoreParam`, so every worker forks again.
- **Thrown-away work.** Each inner `runTests` then pays M3 and M4.
- *Fix:* run a light inner scheme (e.g. one repeat of 5 folds, `SerialParam`). Rank once per inner training set and score every `nFeatures` value on those same inner splits.
- Note that nested tuning of classifier parameters is broken anyway (Q-B1).

**M7. Results carry every fold model.**
- `runTest` returns the fitted model for every fold; they are serialised back from the workers and kept in `@models`.
- Random forest keeps two forests per fold (S2): in a toy run (80 samples × 50 features, 10 × 5) `@models` was 42 MB of a 44 MB result.
- `PrecisionPathways` objects store every fold model per assay (`precisionPathways.R:183`) and copy them again on predict (`:308`).
- *Fix:* keep fold models only when asked (`keepModels = FALSE` by default); keep the final model.

**M8. Thread and core settings.**
- **macOS ignores `nCores` [read].** `crossValidate.R:482` tests `sysname %in% c("MacOS", "Linux")`, but macOS reports `"Darwin"`. Macs therefore fall through to `BiocParallel::bpparam()`, which ignores `nCores` *and* the `RNGseed`. The paper's timings were on a Mac, so they did not use `nCores = 8`, and parallel Mac runs are not reproducible.
- **Default worker count.** `CrossValParams()` defaults to `bpparam()`; on verona that is a `MulticoreParam` with **125 workers**.
- **Uncapped model threads.** ranger (2 threads by default) and rfsrc (all cores through OpenMP) are not limited inside BiocParallel workers (`interfaceRandomForest.R:10-11`, `interfaceRandomForestSurvival.R:12`), so N workers run 2N or more threads.
- *Fix:* test `.Platform$OS.type == "unix"`; default `num.threads = 1` and `rf.cores = 1` inside workers.

**M9. Precision pathways and criss-cross validation repeat work.**
- **`precisionPathwaysPredict`.** It runs a full `runTests` over the *test* cohort (`precisionPathways.R:239`), re-splitting and stratifying it by its labels, just to apply the stored models. It should predict with the stored models directly; that would also stop requiring test labels (Q-C6).
- **`crissCrossValidate(trainType = "modelTest")`.** It runs a full `crossValidate` per dataset only to harvest the selected features, discarding the classifier fits. It then runs k² `runTests`; the diagonal ones repeat work and are leaky (Q-A6).
- **`doRandomFeatures`.** It runs k² further `crossValidate` calls.

**M10. Evaluation after cross-validation.** These are implementation changes with identical output.
- **Sample-wise C-index** (`calcPerformance.R:270-301`): an R loop over sample pairs that also recomputes `as.matrix(Surv)` for every sample. A vectorised version with `outer()` per permutation gave identical values, **22.0 s → 0.18 s** per result (165 samples, 100 × 5).
- **AUC** (`calcPerformance.R:336-357` and `.calcArea`): one `data.frame` per unique score. The rank (Mann–Whitney) formula gave the same values, **9.9 s → 0.03 s** per result. It runs by default in `calcCVperformance(..., "auto")` and `performanceTable()`.
- **ROC curve** (`ROCplot` "merge" mode, `ROCplot.R:142-165`): the same per-score loop, about 45 s per result at 100 permutations. A cumsum-based curve would fix it.

### 1B. Faster implementations inside interfaces and rankings (same results)

- **Levene ranking** (`rankingLevene.R:9-10`): calls `car::leveneTest` once per feature, 3.2 s per 2000 features. An F-test on |x − group median| (`genefilter::colFtests`) gives an **identical ranking** in 0.8 s.
- **Likelihood-ratio ranking** (`rankingLikelihoodRatio.R:9-22`): `mapply` with `dnorm` per feature. A closed-form normal log-likelihood with `colSums` should be over 50× faster (not timed).
- **Bartlett and KS rankings**: `apply` over a `DataFrame` (which converts it to a matrix), one test per feature. Vectorise them.
- **Penalised GLM** (`interfacePenalisedGLM.R:19-23`): makes one `predict()` call per lambda to choose it. A single `predict(..., s = <all lambdas>)` gives the same choice. Training took 1.18 s, of which glmnet was 0.03 s.
- **Naive Bayes kernel and mixture models**: nested `sapply`/`apply` per feature × sample, plus a per-sample `data.frame`/`rbind`. 1.1 s at p = 200.
- **`getLocationsAndScales` / `subtractFromLocation`**: replace `apply` with `colMeans`, `matrixStats` and `sweep`.

### 1C. Statistical choices that also change speed

These change the method, so they belong to the user, not to a speed-up. They are listed so that the costs are
known when choosing defaults.

- **CoxNet lambda selection** (`interfaceCoxnet.R:14`).
  - `cv.glmnet(type = "C")` calls `survival::concordance` once per lambda per inner fold. That is 61% of a single-view CoxNet run.
  - glmnet's default path runs down to `lambda.min.ratio = 1e-4` when n > p, and near-unpenalised Cox fits converge slowly. On the Procedure 2 merge (70 features), one `cv.glmnet` took:
    - **9.2 s** with the defaults;
    - **1.1 s** with `lambda.min.ratio = 0.05`;
    - **0.19 s** with deviance as well.
  - The selected lambda stays inside the shortened path. End to end (5-assay merge plus two smaller combinations, 5 × 5 CV), the run went from **249 s to 116 s** with `lambda.min.ratio = 0.05`, and the mean C-index was about the same (0.600/0.656/0.547 → 0.608/0.637/0.547).
  - Users can already set this: `extraParams = list(train = list(lambda.min.ratio = 0.05))`.
- **Random forest** grows two forests per fit (`interfaceRandomForest.R:10-11`): one to predict, and one with `importance = "impurity_corrected"` only for feature ranking. The default `mTryProportion = 0.5` is about 12× slower than ranger's √p at p = 2000, and it is not standard random forest behaviour.
- **SVM** always fits with `probability = TRUE` (`interfaceSVM.R:10`), which adds an internal 5-fold Platt fit: about 11× slower.

---

## Part 2. Quirks and bugs

### 2A. Wrong results with no error (most serious)

- **Q-A1. `crossValidate` drops `selectionMethod` for MultiAssayExperiment and list input [confirmed].**
  - **Cause.** The `MultiAssayExperimentOrList` method calls the `DataFrame` method without passing `selectionMethod` (`crossValidate.R:336-347`).
  - **Effect.** The default is always used: t-test, or CoxPH ranking for survival. In the test, `selectionMethod = "KS"` ran as t-test, and the result's characteristics *also* say t-test.
  - **Paper.** Procedure 2 passes MAE/list input. The playbook asks for `selectionMethod = "CoxPH"`, which happens to equal the default, so the published results are unaffected; any other choice would have been ignored.
- **Q-A2. DLDA's prior has the wrong sign [confirmed].**
  - **Cause.** `predict.dlda` scores classes as `sum((x - mean)^2 / var) + log(prior)` and takes the minimum (`utilities.R:656`). The discriminant should use `- 2 log(prior)`.
  - **Effect.** Smaller classes are favoured. On pure noise with 90/10 classes, DLDA predicted 97 of 100 samples as the *minority* class.
  - **Also.** The posterior multiplies per-feature densities (`utilities.R:708`), which underflows to 0/0 = NaN with many features.
- **Q-A3. `calcExternalPerformance` relabels predicted levels by position [confirmed].**
  - **Cause.** `levels(predicted) <- levels(actual)` (`calcPerformance.R:107`).
  - **Effect.** A perfect prediction whose factor levels are in a different order scores accuracy 0. This function is used inside tuning and precision pathways.
- **Q-A4. Tuning `nFeatures` by resubstitution always picks the smallest value [confirmed].**
  - **Cause.** Passing a vector of `nFeatures` silently switches on resubstitution tuning (`crossValidate.R:124-134`). With random forest every candidate scored balanced accuracy 1, so `which.max` took the first.
  - **Effect.** In the test, `nFeatures = c(5, 10, 50, 100)` chose 5 in every fold.
  - **Same flaw elsewhere.** The same resubstitution choice is used for train-parameter tuning, for penalised GLM's lambda (`interfacePenalisedGLM.R:16-24`) and for NSC's threshold. The `"auto"` random forest grid includes `num.trees = 1` (`simpleParams.R:3`).
- **Q-A5. GLM is fitted without an intercept [confirmed].**
  - **Cause.** `glm(class ~ . + 0)` (`interfaceGLM.R:7`).
  - **Effect.** One informative feature shifted by +10 gave resubstitution accuracy 0.50, against 0.88 with an intercept. GLM is the default clinical classifier in some paths and the fallback when `nFeatures = 1`.
- **Q-A6. In criss-cross validation, the `modelTest` diagonal is leaky [confirmed].**
  - **Cause.** Features come from `crossValidate` on dataset A. They are then evaluated by a new `runTests` on A with different splits, so test samples helped choose the features.
  - **Effect.** On pure noise the diagonal gave balanced accuracy 0.85, against 0.51 for honest CV. The plot labels it "resubstitution" and hides it by default, but `result$real` returns it unlabelled.
- **Q-A7. CoxNet's feature ranking ignores protective features [confirmed].**
  - **Cause.** The survival branch of `penalisedFeatures` ranks by the signed coefficient (`interfacePenalisedGLM.R:84`), with no `abs()`.
  - **Effect.** A feature with β = −0.14 ranked below six features with zero coefficients.
- **Q-A8. Train and test encode factors differently.**
  - **CoxNet** trains on `MatrixModels::model.Matrix(~ 0 + .)` and predicts on `glmnet::makeX` (`interfaceCoxnet.R:11` vs `:38`); the two encode factors differently.
  - **Penalised GLM and XGB** rebuild the test design matrix separately. A factor level missing from a test fold errors ("non-conformable"), or, if the columns come out in another order, uses the wrong features silently.
- **Q-A9. `samplesMetricMap` aligns samples by position [confirmed].** Results with the same samples in a different order are `rbind`-ed by position and labelled with the first result's sample names, so the columns are mislabelled.
- **Q-A10. Precision pathway costs.**
  - **Cost undercount [confirmed].** `calcCostsAndPerformance` (`precisionPathways.R:345`) charges each sample only for the tier that classified it. A patient sent clinical → miRNA also had the clinical test.
    - Procedure 3 reports $8,400 for clinical-miRNA; it should be $9,300.
    - Clinical-RNA reports $13,800; it should be $14,700.
  - **Weights by position [confirmed].** `summary()` uses the weights by position: `weights = c(cost = 0.9, accuracy = 0.1)` gives accuracy 0.9.
- **Q-A11. Precision pathways match classifiers to models by position [confirmed].** At predict (`precisionPathways.R:234-246`), models come out in alphabetical assay order (`table()`), while classifiers keep the user's order. Procedure 3 works only because clinical < miRNA < RNA happens to sort that way, and `table()` order depends on the locale.
- **Q-A12. A MultiAssayExperiment without `useFeatures$clinical` uses every colData column as a predictor.**
  - **Cause.** `prepareData.R:291` substitutes all colData columns and only warns.
  - **Effect.** Columns derived from the outcome leak into the model: the METABRIC clinical table holds `LR`, `DR`, `DeathBreast`, `TLR` and `TDR` next to `timeRFS`/`eventRFS`.
  - **Paper.** Procedure 1 calls `crossValidate(ghistMAE, outcome = "subtype")` without `useFeatures`, so age and race are used as features without the user saying so.

### 2B. Code paths that fail on any input

- **Q-B1. Nested-CV tuning of classifier parameters [confirmed].** It errors with "this S4 class is not subsettable": `median(performances(result)...)` at `utilities.R:338`, where `performances` is the local variable, not the accessor `performance`.
- **Q-B2. Tuning parameters can vanish without a warning.**
  - Tuning parameters are dropped when a selection step exists and `tuneMode` is `"none"`. **[confirmed]**: `extraParams$train$tuneParams = "auto"` without `tuneCross` was ignored, and the string `"auto"` was passed on as an ordinary parameter (it shows in the characteristics).
  - After tuning, `.doTrain` *replaces* `otherParams` with the chosen combination (`utilities.R:351`), discarding the user's other training settings. **[read]**
- **Q-B3. `train()` fails for multi-assay input [confirmed].** It errors when `cleanClassifier` is called without `nFeatures` (`crossValidate.R:823`). Further problems **[read]**:
  - the merge, prevalidation and PCA branches use the undefined `modellingParams`, `measurementsUse` and `crossValParams` (`crossValidate.R:938-988`);
  - with `selectionMethod = "none"`, each assay's model is trained on *all* assays (`crossValidate.R:853`);
  - `classifierParams$trainParams@otherParams[-inTune]` is a copy-paste error (`:874`, `:895`).
- **Q-B4. The `runTest` method for MultiAssayExperiment fails [read].** It uses the undefined `extrasInputs` and `prepArgs` (`runTest.R:413`, `:419`).
- **Q-B5. `samplesSplits = "Permute Percentage Split"` only works by accident [confirmed].** It uses the undefined `classes` (`utilities.R:72`). It ran only because `data(asthma)` had put a `classes` object in the global environment, so it silently uses whatever `classes` the user has.
- **Q-B6. `prepareData` options are broken.**
  - `topNvariance` fails on the typo `unqiue` (`prepareData.R:238`) **[confirmed]**.
  - `maxSimilarity` does nothing: `if(any(pValues) < ...)` is wrong, and `dropFeatures` is reset to empty before use (`prepareData.R:260`, `:269`). **[confirmed]**: `maxSimilarity = 0.5` kept 300 of 300 features.
  - The "clinical" warning tests `"clinical" %in% is.null(...)` (`:163`), so it can never fire.
- **Q-B7. `nFeatures = 1` fails [confirmed].** It errors with "incorrect number of dimensions": `.doTest` subsets without `drop = FALSE` (`utilities.R:377`). `colCoxTests` and `subtractFromLocation` have the same problem.
- **Q-B8. Classifiers that fail on any input:**
  - `kNN`: `kNNparams()` is called without its required argument;
  - `fisherDiscriminant`, `kTSPclassifier` and Poisson LDA (`classifyInterface`): undefined `trainingMatrix`;
  - naive Bayes and mixture models with `weighting = "crossover distance"`;
  - XGB with xgboost 3.x (old `xgboost(x, y)` API);
  - SVM with factor clinical features (`model.matrix` at predict).
- **Q-B9. Plotting functions that fail:**
  - `ROCplot(mode = "average")`: `split` on a `DataFrame`, and `subset(class = ...)` uses `=`;
  - `samplesMetricMap`: `comparison = "Cross-validation"`, `showLegends = FALSE`, a single result, and more than two classes;
  - `selectionPlot` "importance" mode;
  - `calcExternalPerformance` with more than one metric.
- **Q-B10. Precision pathways fail on the simplest case [confirmed].**
  - The default mode fails for clinical plus one assay: `.permutations` drops to a vector (`utilities.R:728`).
  - `fixedAssays = NULL` errors, and so does list input.
- **Q-B11. `crissCrossValidate` hard-requires `TOP` even with `runTOP = FALSE`** (`crissCrossValidate.R:48-49`), although TOP is only in Suggests. It also does not intersect feature names across datasets, so a gene missing from one dataset gives a subscript error.
- **Q-B12. Unreachable or broken code.** The ensemble-selection branch of `.doSelection` uses the undefined `featuresLists`, `measurementsSubset` and `aResult` (`utilities.R:248-292`). `previousSelection`'s overlap check can never fire. `edgesToHubNetworks` fails on matrix input (`class(x) == "matrix"`).

### 2C. Surprising behaviour

- **Q-C1. Balancing.** `ModellingParams()` defaults to `balancing = "downsample"`, but `runTests` never applies it: only a user-level `runTest` call does (`runTest.R:112`). `crossValidate` sets `"none"`. So `runTests` users think they are downsampling, and they are not.
- **Q-C2. Variable importance.** `doImportance` always uses Balanced Error: it reads `performanceType` from `selectParams@tuneParams`, where it never is (`runTest.R:262`). For survival this is the wrong metric.
- **Q-C3. Inner-loop seeds.** Prevalidation and PCA seed their inner CV with `.Random.seed[1]`, which is the RNG *kind code* (10403), the same for every fold. They also error in a fresh session.
- **Q-C4. Sample-wise metrics.**
  - The sample-wise C-index compares risk scores from different fold models (with the default `grouping = "permutation"`), drops tied pairs (unlike the C-index), and is rounded to 2 dp.
  - The number of comparable pairs per sample ranged from 1 to 164, so some cells of the heatmap rest on one pair per permutation.
  - NA arises for any censored patient censored before the earliest event, not only "the shortest censored patient".
  - AUC is rounded to 2 dp before averaging, and is the unweighted mean of one-vs-rest AUCs.
- **Q-C5. `grouping = "fold"`.** It averages C-index and AUC per permutation (returning a `by` object), but not Balanced Accuracy, which returns 500 per-fold values. Its documentation says grouping makes no difference to accuracy metrics.
- **Q-C6. Precision pathways: how samples are routed.**
  - "Confidence" is agreement between resampled models, not accuracy.
  - The cut-off is strict `>`: with 20 repeats and cut-off 0.8, 18/20 agreement fails, so in practice the rule is ≥ 19/20.
  - At prediction time, `minAssaySamples` makes one patient's routing depend on the rest of the test batch.
  - Prediction requires test labels.
  - Pathways assume two classes.
- **Q-C7. Plots.**
  - `performancePlot`, `bubblePlot` and others call `ggplot2::theme_set()`, changing the user's global theme.
  - `performancePlot` overrides user `yLimits` when `rotate90 = TRUE`, and draws a chance line at 0.5 whatever the metric or number of classes.
  - `bubblePlot` silently drops pathways below 0.5 accuracy.
- **Q-C8. Two-class-only methods.** KS, KL, pairs-differences, Fisher and kTSP silently use only the first two classes. Prevalidation fails for more than two classes.
- **Q-C9. Seed handling, and folds that differ between assays [confirmed].**
  - `crossValidate` refuses to run without `set.seed`. It derives the BiocParallel seed by reading `.Random.seed`
    directly (`crossValidate.R:471`), which assumes Mersenne-Twister.
  - Each assay, classifier and combination is cross-validated on **different folds**: the splits are drawn from the
    global random stream, which earlier cross-validations have already advanced. In the harness, two assays'
    fold assignments agreed for 19% of samples, which is chance for 5 folds.
  - So comparisons between assays or combinations (`performancePlot`, `samplesMetricMap`) are not paired by fold, and
    they include split-to-split noise. With shared folds, the comparisons would be paired, and the per-assay
    rankings in merge, prevalidation and PCA could be computed once and reused (M5).

### 2D. Packaging and documentation

- **DESCRIPTION.**
  - `Packaged: 2014-10-18` is stale.
  - `survival`, `generics`, `S4Vectors`, `MultiAssayExperiment` and `BiocParallel` are in Depends; Imports is the current Bioconductor recommendation.
  - `ggpubr`, `dcanr`, `reshape2`, `broom` and `ggupset` are each imported for one or two calls.
  - The Title still advertises "differential variability and differential distribution testing".
- **NAMESPACE.** `import(ggplot2, ggpubr, reshape2, grid)` brings in whole namespaces. `S4Vectors::DataFrame`, `combn`, `setNames` and `na.omit` are used without being imported; they only work because the packages are attached.
- **Tests.** There is no `tests/` folder. Most of the bugs above would be caught by a smoke test calling each `available()` keyword once.
- **Documentation errors.**
  - Precision and recall are swapped in `calcPerformance`'s documentation.
  - The vignette says the C-index "is simply the number of pairs".
  - `performancePlot` says the default is a violin plot; it is a box plot.
  - `METABRICclinical` says "eight features" (two of them are the outcome) and cites the wrong paper's URL.
  - `precisionPathways`, `crissCrossValidate` and `crissCrossPlot` have no examples.

---

## Relevance to the rejected protocol paper

1. **Procedure 2 merge (27 min).** About 80% of the time is the glmnet Cox solver walking the lambda path to near-zero penalty. The rest is the machinery (M2–M5), and parallel scaling was poor (about 2× on 8 cores). The machinery fixes alone should give several-fold; the lambda path is a separate, statistical choice (1C).
2. **Mac timings.** Every timing in the paper was taken on a Mac, where `nCores = 8` was ignored and all cores were used (M8).
3. **Results that may need re-checking:**
   - Procedure 1: clinical columns used without `useFeatures` (Q-A12);
   - Procedure 3: pathway costs (Q-A10); classifier order, correct only by luck (Q-A11);
   - Procedure 4: the leaky diagonal, if shown (Q-A6);
   - any DLDA result (Q-A2).
4. **Definitions the paper should state:** the sample-wise C-index and its NA rule (Q-C4), and what "confidence" means in precision pathways (Q-C6).

## Reproducing

Benchmark and confirmation scripts are in the ProjectClassifyR repository (`scripts/audit_2026-10/`).
