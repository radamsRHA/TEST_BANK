---
title: "Discriminating Models of Trait Evolution with EvoDA"
author: "Jenny Roa Lozano, Michael DeGiorgio, Raquel Assis, and Rich Adams"
format:
  html:
    toc: true
    toc-title: "🔍 EvoDA Tutorial Overview"
    toc-location: left
    code-fold: false
    theme: cosmo
    highlight-style: github
    df-print: paged
---

#### Author affiliations

Jenny Roa Lozano<sup>1,2</sup>, Michael DeGiorgio<sup>3</sup>, Raquel Assis<sup>3,4</sup>, and Rich Adams<sup>1,2\*</sup>\
<sup>1</sup> Center for Agricultural Data Analytics, University of Arkansas, Fayetteville, AR\
<sup>2</sup> Department of Entomology and Plant Pathology, University of Arkansas, Fayetteville, AR\
<sup>3</sup> Department of Electrical Engineering and Computer Science, Florida Atlantic University, Boca Raton, FL\
<sup>4</sup> Institute for Human Health and Disease Intervention, Florida Atlantic University, Boca Raton, FL

# Introduction

Phylogenetic comparative methods (PCMs) have revolutionized how we study trait evolution across species. Traditional approaches such as AIC-based model selection have proven powerful, but they may struggle in the face of model complexity, measurement error, and noisy empirical data. In contrast, supervised learning, especially discriminant analysis, offers a promising alternative for classifying traits into evolutionary models using labeled training data.

In this tutorial, we introduce **Evolutionary Discriminant Analysis (EvoDA)**, a supervised learning framework designed to predict the evolutionary processes shaping trait variation. EvoDA leverages training and testing datasets simulated under known models, and performs classification with high accuracy—even when conventional methods fall short. This guide walks you through every step using real fungal phylogeny case studies.

::: {.callout-tip collapse="true" title="💡 Tip: Tutorial Manual and where find it"}
See the file https://github.com/radamsRHA/TraitTrainR/blob/main/TraitTrainR_CommandManual.pdf for detailed instructions on functions NOTE See the file https://github.com/radamsRHA/TraitTrainR/blob/main/TraitTrainR_TutorialManual.pdf for detailed example Run
:::

::: {.callout-tip collapse="true" title="💡 Tip: Before to start README"}
This README serves as a quick start guide. Please see the above tutorial files for more indepth instructions and details of options and implementation
:::

# Tutorial Outline

Jump to the section you need:

-   [Installing Dependencies](#installing-dependencies)
-   [Quick Example using TraitTrainR](#quick-example-using-traittrainr)
-   [Training Data: Fungal Phylogeny - Case Study II](#case-study-ii-training-data-fungal-phylogeny)
-   [Testing Data: Fungal Phylogeny](#testing-data)
-   [Rejection Sampling](#rejection-sampling-coming-soon)

------------------------------------------------------------------------

# Installing Dependencies {#installing-dependencies}

Before using EvoDA, install the required R packages. This includes the EvoDA helper package `TraitTrainR`:

::: {.callout-tip collapse="true" title="💡 Tip: Note: R version"}
TraitTrainR was written in R 4.4.0 ("Puppy cup") and we recommend that version or later for installing `TraitTrainR`
:::

::: cell
``` r
install.packages("devtools")
library(devtools)
install_github("radamsRHA/TraitTrainR")

install.packages("phytools")   # version >= 2.1-1
install.packages("geiger")     # version >= 2.0.11

library(TraitTrainR)
library(phytools)
library(geiger)
```
:::

------------------------------------------------------------------------

# TraitTrainR Evolutionary Models Overview

Before exploring the example, it's helpful to understand the evolutionary models available in `TraitTrainR`. The table below summarizes the parameters, evolutionary processes, and distributions used for training and testing simulations. These distributions align with default assumptions in `fitContinuous` from the geiger package.

```{r, echo=FALSE}
library(knitr)
kable(
  data.frame(
    Model = c("BM", "OU", "EB", "Lambda", "Delta", "Kappa", "Trend"),
    Parameters = c("z0; σ²", "z0; σ²; α", "z0; σ²; a", "z0; σ²; λ", "z0; σ²; δ", "z0; σ²; κ", "z0; σ²; slope"),
    Evolutionary_Process = c(
      "Random-walk model of trait change (Felsenstein 1973)",
      "Stabilizing selection (Butler and King 2004)",
      "Adaptive radiation: rates decline over time (Harmon et al. 2010)",
      "Phylogenetic signal scaling (Pagel 1999)",
      "Time-dependent evolutionary rate (Pagel 1999)",
      "Branch length scaling in phylogeny (Pagel 1999)",
      "Linear change in trait evolution over time (Pennell et al. 2014)"
    ),
    Distributions = c(
      "σ² ~ Exp(1); z0 ~ N(0, 1)",
      "σ² ~ Exp(1); z0 ~ N(0,1); α ~ U(exp(-500), exp(1))",
      "σ² ~ Exp(1); z0 ~ N(0, 1); a ~ U(-5/depth, -1e-6)",
      "σ² ~ Exp(1); z0 ~ N(0, 1); λ ~ U(exp(-500), 1)",
      "σ² ~ Exp(1); z0 ~ N(0, 1); δ ~ U(exp(-500), 3)",
      "σ² ~ Exp(1); z0 ~ N(0, 1); κ ~ U(exp(-500), 1)",
      "σ² ~ Exp(1); z0 ~ N(0, 1); slope ~ U(-100, 100)"
    )
  ),
  caption = "**Table 1. Trait models, parameters, evolutionary processes, and sampling distributions used in TraitTrainR.**",
  align = 'llll'
)
```

------------------------------------------------------------------------

# QUICK START CODE: simulate under Brownian Motion and Ornstein–Uhlenbeck models (see next section for more detailed information): {#quick-example-using-traittrainr}

To demonstrate `TraitTrainR` functionality, here is an example using three basic models (BM, OU, EB) and a simple tree.

::: {.callout-tip collapse="true" title="💡 Tip: Why start simple?"}
This quick example serves as a sandbox to understand the structure of simulations in `TraitTrainR`. It ensures you’re comfortable before scaling up to realistic phylogenies and complex parameterizations.
:::

::: cell
``` r
library(TraitTrainR)
library(phytools)
library(geiger)

MyTree <- read.tree(text = "((A:1, B:1):1, C:2);") 
plot(MyTree)

list.SimulationModelSettings <- list()
NumberOfReps <- 5
list.Rmatrix <- lapply(1:NumberOfReps, function(i) matrix(1, 1, 1))

list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "BM", 
                                          vector.Sig2 = rexp(NumberOfReps, rate = 1), 
                                          vector.AncestralState = rep(1, NumberOfReps), 
                                          list.Rmatrix = list.Rmatrix)

list.SimulationModelSettings[[2]] <- list(string.SimulationModel = "OU", 
                            vector.Sig2 = rexp(NumberOfReps, rate = 1), 
                            vector.AncestralState = rnorm(NumberOfReps),
                            vector.Alpha = runif(NumberOfReps, min = exp(-500), max = exp(1)),
                            list.Rmatrix = list.Rmatrix)

list.SimulationModelSettings[[3]] <- list(string.SimulationModel = "EB", 
                                     vector.Sig2 = rexp(NumberOfReps, rate = 1), 
                                     vector.AncestralState = rnorm(NumberOfReps), 
                                     vector.A = runif(NumberOfReps, min = log(10^-5)/310, max = -0.000001),
                                     list.Rmatrix = list.Rmatrix)

MySimulationResults <- TraitTrain(handle.Phylogeny = MyTree,
                       list.SimulationModelSettings = list.SimulationModelSettings,
                       logical.PIC = TRUE, logical.PROJECT = TRUE)
```
:::

------------------------------------------------------------------------

Now that we've explored a basic example using a toy tree, let's take things to the next level.

We’re now ready to dive into the real application: a **fungal phylogeny**. This case study will provide the training data needed to power EvoDA’s discriminant learning approach for model classification.

# Case Study II: Training and Testing data with Fungal phylogeny {#case-study-ii-training-data-fungal-phylogeny}

In this case study, we use a realistic fungal phylogeny composed of 17 species, including major taxa like *Aspergillus*, *Neurospora*, *Saccharomyces*, and *Cryptococcus*. This deep and diverse tree reflects a broad evolutionary timespan, ideal for exploring evolutionary model discrimination. Here, we generate trait data under three models—Brownian Motion (BM), Ornstein-Uhlenbeck (OU), and Early Burst (EB)—and divide the workflow into training and testing phases to support EvoDA analysis.

------------------------------------------------------------------------

## 🧪 Training Data Simulation

```{r}
# Load dependencies
rm(list = ls())
library(TraitTrainR)
library(ape)
library(geiger)
library(phytools)
library(caret)
library(corrplot)
library(MASS)
library(mda)
```

------------------------------------------------------------------------

### Load fungal phylogeny

::: cell
```{r}

handle.TargetTree <- read.tree(text = "((((A.fumigatus:138.308934,A.nidulans:138.308934):191.332495,((M.oryzae:175.091151,(N.discreta:14.736545,(N.tetrasperma:5.443496,N.crassa:5.443496):9.293049):160.354606):33.404309,F.graminearum:208.49546):121.145969):232.724291,((((N.castellii:145.870505,(S.bayanus:60.615136,(((S.paradoxus:25.739496,S.cerevisiae:25.739496):12.037587,S.mikatae:37.777084):13.200718,S.kudriavzevii:50.977802):9.637334):85.255369):21.321276,C.glabrata:167.191781):29.950882,L.kluyveri:197.142664):111.428464,(C.albicans:135.842772,C.parapsilosis:135.842772):172.728356):253.794591):281.634281,C.neoformans:844);")
```
:::

### Define simulation parameters for training

::: cell
```{r}

numeric.MeasurementError <- 0
numeric.NumberTrainingReps <- 100   #1e5
list.Rmatrix <- lapply(1:numeric.NumberTrainingReps, function(i) matrix(1, 1, 1))

# Define model settings for training
list.SimulationModelSettings <- list()
list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "BM", vector.Sig2 = rexp(numeric.NumberTrainingReps, rate = 1), vector.AncestralState = rnorm(numeric.NumberTrainingReps), list.Rmatrix = list.Rmatrix)

list.SimulationModelSettings[[2]] <- list(string.SimulationModel = "OU", vector.Sig2 = rexp(numeric.NumberTrainingReps, rate = 1), vector.AncestralState = rnorm(numeric.NumberTrainingReps), vector.Alpha = runif(numeric.NumberTrainingReps, min = exp(-500), max = exp(1)), list.Rmatrix = list.Rmatrix)

list.SimulationModelSettings[[3]] <- list(string.SimulationModel = "EB", vector.Sig2 = rexp(numeric.NumberTrainingReps, rate = 1), vector.AncestralState = rnorm(numeric.NumberTrainingReps), vector.A = runif(numeric.NumberTrainingReps, min = log(10^-5)/373.1363, max = -0.000001), list.Rmatrix = list.Rmatrix)

```
:::

------------------------------------------------------------------------

### Simulate training traits

::: cell
```{r}
handle.RESULTS_TRAIN <- TraitTrain(handle.Phylogeny = handle.TargetTree, list.SimulationModelSettings = list.SimulationModelSettings, logical.PIC = TRUE, logical.PROJECT = TRUE, numeric.MeasurementError = numeric.MeasurementError)

```
:::

------------------------------------------------------------------------

## 🧪 Testing Data Simulation

### Define simulation parameters for testing

::: cell
```{r}
numeric.NumberTestingReps <- 50   #1e3
list.Rmatrix <- lapply(1:numeric.NumberTestingReps, function(i) matrix(1, 1, 1))

# Define model settings for testing
list.SimulationModelSettings <- list()
list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "BM", vector.Sig2 = rexp(numeric.NumberTestingReps, rate = 1), vector.AncestralState = rnorm(numeric.NumberTestingReps), list.Rmatrix = list.Rmatrix)

list.SimulationModelSettings[[2]] <- list(string.SimulationModel = "OU", vector.Sig2 = rexp(numeric.NumberTestingReps, rate = 1), vector.AncestralState = rnorm(numeric.NumberTestingReps), vector.Alpha = runif(numeric.NumberTestingReps, min = exp(-500), max = exp(1)), list.Rmatrix = list.Rmatrix)

list.SimulationModelSettings[[3]] <- list(string.SimulationModel = "EB", vector.Sig2 = rexp(numeric.NumberTestingReps, rate = 1), vector.AncestralState = rnorm(numeric.NumberTestingReps), vector.A = runif(numeric.NumberTestingReps, min = log(10^-5)/310, max = -0.000001), list.Rmatrix = list.Rmatrix)

```
:::

## 🧬 Simulate test traits

```{r}
handle.RESULTS_TEST <- TraitTrain(
  handle.Phylogeny = handle.TargetTree,
  list.SimulationModelSettings = list.SimulationModelSettings,
  logical.PIC = TRUE,
  logical.PROJECT = TRUE,
  numeric.MeasurementError = numeric.MeasurementError)
```


# Prepare training and testing datasets for EvoDA

::: cell
```{r}

DATA_TRAIN <- handle.RESULTS_TRAIN$RESULTS_PIC
DATA_TRAIN$SimulationModelNumber <- factor(DATA_TRAIN$SimulationModelNumber)

DATA_TEST <- handle.RESULTS_TEST$RESULTS_PIC
DATA_TEST$SimulationModelNumber <- factor(DATA_TEST$SimulationModelNumber)

```
:::

------------------------------------------------------------------------

# 🎯 Rejection Sampling: {#rejection-sampling:-Filtering-Simulted-Traits-with-Empirical-Data}

Now that we’ve generated both training and testing datasets under three evolutionary models (BM, OU, and EB), it’s time to bring the **real data** into play.

To make our simulations biologically meaningful, we **filter** them using **rejection sampling** based on **Fungal gene expression data**.

This ensures that simulated traits (used to train EvoDA) exhibit **variation comparable to empirical gene expression traits**, making downstream model classification more robust and interpretable.

::: {.callout-tip title="💡 What is Rejection Sampling?"}
Rejection sampling is a filtering method where only data points that meet specific empirical criteria are retained. In our case, we compute the variance of phylogenetically independent contrasts (PICs) from real gene expression data and use that variance as a benchmark to select comparable simulations.
:::

::: {.callout-tip title="🧠 Why PICs?"}
PICs remove phylogenetic autocorrelation from trait data, allowing you to assess trait variation independently of shared ancestry. This is key for accurate comparison between empirical and simulated data.
:::

------------------------------------------------------------------------

::: cell
```{r}

# Load depends 
library(dplyr)
library(tidyr)

# Load input GEX data 
handle.GEX_Data <- as_tibble(read.table(file = '~/Downloads/EvoDA_Tutorial/EvoDA_Tutorial/gene_expression_tpm_matrix_updated_Standard.LogNorm.tsv', header = TRUE))

# Filter the GEX data based on missing values 
handle.GEX_Data <- handle.GEX_Data %>% drop_na()  # remove any genes with NA

# Check overlap of species 
vector.KeepNames <- c(handle.TargetTree$tip.label, colnames(handle.GEX_Data)[colnames(handle.GEX_Data) != "Protein"])
vector.KeepNames <- names(table(vector.KeepNames))[table(vector.KeepNames) == 2]

# Only keep species in both tree and traits 
handle.GEX_Data <- handle.GEX_Data[, vector.KeepNames]
handle.TargetTree <- drop.tip(handle.TargetTree, handle.TargetTree$tip.label[-match(vector.KeepNames, handle.TargetTree$tip.label)])

# Compute PICs 
handle.GEX_PICs <- apply(t(handle.GEX_Data), MARGIN = 2, FUN = function(x) pic(x = x, phy = handle.TargetTree))
handle.GEX_PICs <- t(handle.GEX_PICs)

# Rename columns to match DATA_TEST
colnames(handle.GEX_PICs) <- colnames(DATA_TEST)[-ncol(DATA_TEST)]


```
:::

------------------------------------------------------------------------

📌 *Next step*: Use `handle.GEX_PICs` to define empirical variance bounds, and apply rejection sampling to both training and testing simulations to ensure they fall within biologically plausible ranges.

------------------------------------------------------------------------

------------------------------------------------------------------------

# Rejection Sampling (Coming Soon) {#rejection-sampling-coming-soon}

------------------------------------------------------------------------

# Rejection Sampling (Coming Soon) {#rejection-sampling-coming-soon}

------------------------------------------------------------------------

This step ensures only biologically plausible test data are retained by comparing trait variance against empirical bounds.

------------------------------------------------------------------------

# References (Coming Soon) {#References}
