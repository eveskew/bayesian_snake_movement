# Snake Movement and Activity Analyses using Bayesian Modeling

This repository contains code for in-progress work on snake movement and activity patterns. Data were collected in the southeastern United States, and the general objective is to identify environmental factors influencing snake activity. Hopefully, this code may be of use to those trying to implement Bayesian models in R, as code for all Stan models are included. Count data representing snake observations are modeled using a parameterization of the negative binomial distribution that is described by a mean and dispersion parameter. Models are fit using the `rstan` package, while model comparison is conducted using the `loo` package. Various functions from `rethinking` are also used to summarize and visualize model results. 

### Repository Contents

- The .R script file contains my in-progress code and some commentary
- All the raw data files necessary to conduct the analysis
