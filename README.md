# Snake Movement and Activity using Bayesian Modeling

This repository contains code for work that will appear in publication as:

Eskew, E.A., and B.D. Todd. *In press*. Too cold, too wet, too bright, or just right? Environmental predictors of snake movement and activity. Copeia.

Data were collected in the southeastern United States, and the general objective is to identify environmental factors influencing snake activity. Hopefully, this code may be of use to those trying to implement Bayesian models in R, as code for all Stan models are included. Count data representing snake observations are modeled using a parameterization of the negative binomial distribution that is described by a mean and dispersion parameter. Models are fit using the `rstan` package, while model comparison is conducted using the `loo` package. Various functions from `rethinking` are also used to summarize and visualize model results.
 

### Repository Contents

- `bayesian_snake_movement.R` is the primary script file containing code and commentary
- The `R` directory contains code for functions and the Stan models that are sourced in the main script
- The `data` directory contains the raw data files necessary to conduct the analysis
- The `outputs` directory contains figures in both PNG and TIFF formats
