# Snake Movement Analyses using Bayesian Modeling

This repository contains code for in-progress work on snake movement/activity patterns. Data were collected in the southeastern United States, and the general objective is to identify environmental factors influencing snake movement. Potentially (hopefully?) this code may be of use to those trying to implement Bayesian models in R, as code for all Stan models are included. Count data representing snake observations are modeled using a parameterization of the negative binomial distribution that is described by a mean and dispersion parameter. Models are fit using the `rstan` package, while model comparison is conducted using the `loo` package. Various functions from `rethinking` are also used to summarize and visualize model results. Currently I am not sharing the raw data, but as this project gets closer to publication that may change. 

### Repository Contents

- The .R script file contains my in-progress code and some commentary
