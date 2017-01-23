

# Model m1
# Main effects: trapping effort, Julian day, Julian day squared
# Varying intercepts by: species, bay

m1.nbinom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Count[N];
real TrapEffort[N];
real JulianDay[N];
real JulianDaySquared[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
int N_Bay;
int Bay[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bJulianDay;
real bJulianDaySquared;
real<lower=0> phi;

// declare parameters related to varying effects
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
vector[N_Bay] a_Bay;
real<lower=0> sigma_Bay;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + a_Bay[Bay[i]] + bTrapEffort*TrapEffort[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i];
}

Count ~ neg_binomial_2_log(eta, phi);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);
phi ~ cauchy(0, 1);

// priors for varying effects
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
a_Bay ~ normal(0, sigma_Bay);
sigma_Bay ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + a_Bay[Bay[i]] + bTrapEffort*TrapEffort[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i];

log_lik[i] <- neg_binomial_2_log_log(Count[i], eta[i], phi);
}

}"

#==============================================================================


# Model m2
# Main effects: trapping effort, Julian day, Julian day squared, precipitation
# Varying intercepts by: species, bay

m2.nbinom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Count[N];
real TrapEffort[N];
real JulianDay[N];
real JulianDaySquared[N];
real Precip[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
int N_Bay;
int Bay[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bJulianDay;
real bJulianDaySquared;
real bPrecip;
real<lower=0> phi;

// declare parameters related to varying effects
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
vector[N_Bay] a_Bay;
real<lower=0> sigma_Bay;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + a_Bay[Bay[i]] + bTrapEffort*TrapEffort[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bPrecip*Precip[i];
}

Count ~ neg_binomial_2_log(eta, phi);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);
bPrecip ~ normal(0, 10);
phi ~ cauchy(0, 1);

// priors for varying effects
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
a_Bay ~ normal(0, sigma_Bay);
sigma_Bay ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + a_Bay[Bay[i]] + bTrapEffort*TrapEffort[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bPrecip*Precip[i];

log_lik[i] <- neg_binomial_2_log_log(Count[i], eta[i], phi);
}

}"

#==============================================================================


# Model m3
# Main effects: trapping effort, Julian day, Julian day squared, temperature
# Varying intercepts by: species, bay

m3.nbinom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Count[N];
real TrapEffort[N];
real JulianDay[N];
real JulianDaySquared[N];
real Tmin[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
int N_Bay;
int Bay[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bJulianDay;
real bJulianDaySquared;
real bTmin;
real<lower=0> phi;

// declare parameters related to varying effects
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
vector[N_Bay] a_Bay;
real<lower=0> sigma_Bay;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + a_Bay[Bay[i]] + bTrapEffort*TrapEffort[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bTmin*Tmin[i];
}

Count ~ neg_binomial_2_log(eta, phi);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);
bTmin ~ normal(0, 10);
phi ~ cauchy(0, 1);

// priors for varying effects
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
a_Bay ~ normal(0, sigma_Bay);
sigma_Bay ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + a_Bay[Bay[i]] + bTrapEffort*TrapEffort[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bTmin*Tmin[i];

log_lik[i] <- neg_binomial_2_log_log(Count[i], eta[i], phi);
}

}"

#==============================================================================


# Model m4
# Main effects: trapping effort, Julian day, Julian day squared, 
# lunar brightness
# Varying intercepts by: species, bay

m4.nbinom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Count[N];
real TrapEffort[N];
real JulianDay[N];
real JulianDaySquared[N];
real LunarBright[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
int N_Bay;
int Bay[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bJulianDay;
real bJulianDaySquared;
real bLunarBright;
real<lower=0> phi;

// declare parameters related to varying effects
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
vector[N_Bay] a_Bay;
real<lower=0> sigma_Bay;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + a_Bay[Bay[i]] + bTrapEffort*TrapEffort[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bLunarBright*LunarBright[i];
}

Count ~ neg_binomial_2_log(eta, phi);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);
bLunarBright ~ normal(0, 10);
phi ~ cauchy(0, 1);

// priors for varying effects
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
a_Bay ~ normal(0, sigma_Bay);
sigma_Bay ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + a_Bay[Bay[i]] + bTrapEffort*TrapEffort[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bLunarBright*LunarBright[i];

log_lik[i] <- neg_binomial_2_log_log(Count[i], eta[i], phi);
}

}"

#==============================================================================


# Model m5
# Main effects: trapping effort, Julian day, Julian day squared, precipitation,
# temperature, lunar brightness
# Varying intercepts by: species, bay

m5.nbinom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Count[N];
real TrapEffort[N];
real JulianDay[N];
real JulianDaySquared[N];
real Precip[N];
real Tmin[N];
real LunarBright[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
int N_Bay;
int Bay[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bJulianDay;
real bJulianDaySquared;
real bPrecip;
real bTmin;
real bLunarBright;
real<lower=0> phi;

// declare parameters related to varying effects
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
vector[N_Bay] a_Bay;
real<lower=0> sigma_Bay;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + a_Bay[Bay[i]] + bTrapEffort*TrapEffort[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bPrecip*Precip[i] + bTmin*Tmin[i] + bLunarBright*LunarBright[i];
}

Count ~ neg_binomial_2_log(eta, phi);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);
bPrecip ~ normal(0, 10);
bTmin ~ normal(0, 10);
bLunarBright ~ normal(0, 10);
phi ~ cauchy(0, 1);

// priors for varying effects
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
a_Bay ~ normal(0, sigma_Bay);
sigma_Bay ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + a_Bay[Bay[i]] + bTrapEffort*TrapEffort[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bPrecip*Precip[i] + bTmin*Tmin[i] + bLunarBright*LunarBright[i];

log_lik[i] <- neg_binomial_2_log_log(Count[i], eta[i], phi);
}

}"

#==============================================================================


# Model m6
# Main effects: trapping effort, trapped during day only?, Julian day, 
# Julian day squared
# Varying intercepts by: species

m6.nbinom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Count[N];
real TrapEffort[N];
int DayOnly[N];
real JulianDay[N];
real JulianDaySquared[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bDayOnly;
real bJulianDay;
real bJulianDaySquared;
real<lower=0> phi;

// declare parameters related to varying effects
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i];
}

Count ~ neg_binomial_2_log(eta, phi);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bDayOnly ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);
phi ~ cauchy(0, 1);

// priors for varying effects
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i];

log_lik[i] <- neg_binomial_2_log_log(Count[i], eta[i], phi);
}

}"

#==============================================================================


# Model m7
# Main effects: trapping effort, trapped during day only?, Julian day, 
# Julian day squared, precipitation
# Varying intercepts by: species

m7.nbinom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Count[N];
real TrapEffort[N];
int DayOnly[N];
real JulianDay[N];
real JulianDaySquared[N];
real Precip[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bDayOnly;
real bJulianDay;
real bJulianDaySquared;
real bPrecip;
real<lower=0> phi;

// declare parameters related to varying effects
vector[N_Species] bPrecip_Species;
real<lower=0> sigma_bPrecip_Species;
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
(bPrecip + bPrecip_Species[Species[i]])*Precip[i];
}

Count ~ neg_binomial_2_log(eta, phi);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bDayOnly ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);
bPrecip ~ normal(0, 10);
phi ~ cauchy(0, 1);

// priors for varying effects
bPrecip_Species ~ normal(0, sigma_bPrecip_Species);
sigma_bPrecip_Species ~ cauchy(0, 1);
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
(bPrecip + bPrecip_Species[Species[i]])*Precip[i];

log_lik[i] <- neg_binomial_2_log_log(Count[i], eta[i], phi);
}

}"

#==============================================================================


# Model m8
# Main effects: trapping effort, trapped during day only?, Julian day, 
# Julian day squared, temperature
# Varying intercepts by: species

m8.nbinom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Count[N];
real TrapEffort[N];
int DayOnly[N];
real JulianDay[N];
real JulianDaySquared[N];
real Tmin[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bDayOnly;
real bJulianDay;
real bJulianDaySquared;
real bTmin;
real<lower=0> phi;

// declare parameters related to varying effects
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bTmin*Tmin[i];
}

Count ~ neg_binomial_2_log(eta, phi);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bDayOnly ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);
bTmin ~ normal(0, 10);
phi ~ cauchy(0, 1);

// priors for varying effects
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bTmin*Tmin[i];

log_lik[i] <- neg_binomial_2_log_log(Count[i], eta[i], phi);
}

}"

#==============================================================================


# Model m9
# Main effects: trapping effort, trapped during day only?, Julian day, 
# Julian day squared, lunar brightness (plus diurnality)
# Varying intercepts by: species

m9.nbinom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Count[N];
real TrapEffort[N];
int DayOnly[N];
real JulianDay[N];
real JulianDaySquared[N];
real LunarBright[N];
int Diurnal[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bDayOnly;
real bJulianDay;
real bJulianDaySquared;
real bLunarBright;
real bDiurnal;
real bLunarBrightDiurnal;
real<lower=0> phi;

// declare parameters related to varying effects
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bLunarBright*LunarBright[i] + bDiurnal*Diurnal[i] +
bLunarBrightDiurnal*LunarBright[i]*Diurnal[i];
}

Count ~ neg_binomial_2_log(eta, phi);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bDayOnly ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);
bLunarBright ~ normal(0, 10);
bDiurnal ~ normal(0, 10);
bLunarBrightDiurnal ~ normal(0, 10);
phi ~ cauchy(0, 1);

// priors for varying effects
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
bLunarBright*LunarBright[i] + bDiurnal*Diurnal[i] +
bLunarBrightDiurnal*LunarBright[i]*Diurnal[i];

log_lik[i] <- neg_binomial_2_log_log(Count[i], eta[i], phi);
}

}"

#==============================================================================


# Model m10
# Main effects: trapping effort, trapped during day only?, Julian day, 
# Julian day squared, precipitation, temperature, lunar brightness (plus
# diurnality)
# Varying intercepts by: species

m10.nbinom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Count[N];
real TrapEffort[N];
int DayOnly[N];
real JulianDay[N];
real JulianDaySquared[N];
real Precip[N];
real Tmin[N];
real LunarBright[N];
int Diurnal[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bDayOnly;
real bJulianDay;
real bJulianDaySquared;
real bPrecip;
real bTmin;
real bLunarBright;
real bDiurnal;
real bLunarBrightDiurnal;
real<lower=0> phi;

// declare parameters related to varying effects
vector[N_Species] bPrecip_Species;
real<lower=0> sigma_bPrecip_Species;
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
(bPrecip + bPrecip_Species[Species[i]])*Precip[i] +
bTmin*Tmin[i] + bLunarBright*LunarBright[i] + bDiurnal*Diurnal[i] +
bLunarBrightDiurnal*LunarBright[i]*Diurnal[i];
}

Count ~ neg_binomial_2_log(eta, phi);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bDayOnly ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);
bPrecip ~ normal(0, 10);
bTmin ~ normal(0, 10);
bLunarBright ~ normal(0, 10);
bDiurnal ~ normal(0, 10);
bLunarBrightDiurnal ~ normal(0, 10);
phi ~ cauchy(0, 1);

// priors for varying effects
bPrecip_Species ~ normal(0, sigma_bPrecip_Species);
sigma_bPrecip_Species ~ cauchy(0, 1);
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] eta;

for (i in 1:N) {

eta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i] +
(bPrecip + bPrecip_Species[Species[i]])*Precip[i] +
bTmin*Tmin[i] + bLunarBright*LunarBright[i] + bDiurnal*Diurnal[i] +
bLunarBrightDiurnal*LunarBright[i]*Diurnal[i];

log_lik[i] <- neg_binomial_2_log_log(Count[i], eta[i], phi);
}

}"

#==============================================================================


# Example global binomial regression model
# Same model as m10.nbinom, expect fit with a binomial distribution

m10.binom.modelcode <- 
  
  "data {

// number of observations
int<lower=1> N; 

// declare outcome and predictor variables
int<lower=0> Success[N];
real TrapEffort[N];
int DayOnly[N];
real Precip[N];
real Tmin[N];
real LunarBright[N];
int Diurnal[N];
real JulianDay[N];
real JulianDaySquared[N];

// declare variables related to clustering units
int N_Species; 
int Species[N];
}

///////////////////////////////////////////////////////////////////////////////
parameters {

// declare parameters related to main effects
real a;
real bTrapEffort;
real bDayOnly;
real bPrecip;
real bTmin;
real bLunarBright;
real bDiurnal;
real bLunarBrightDiurnal;
real bJulianDay;
real bJulianDaySquared;

// declare parameters related to varying effects
vector[N_Species] bPrecip_Species;
real<lower=0> sigma_bPrecip_Species;
vector[N_Species] a_Species; 
real<lower=0> sigma_Species;
}

///////////////////////////////////////////////////////////////////////////////
model {

vector[N] theta;

for (i in 1:N) {

theta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
(bPrecip + bPrecip_Species[Species[i]])*Precip[i] +
bTmin*Tmin[i] + bLunarBright*LunarBright[i] + bDiurnal*Diurnal[i] +
bLunarBrightDiurnal*LunarBright[i]*Diurnal[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i];

}

Success ~ bernoulli_logit(theta);

// priors

// priors for main effects
a ~ normal(0, 10);
bTrapEffort ~ normal(0, 10);
bDayOnly ~ normal(0, 10);
bPrecip ~ normal(0, 10);
bTmin ~ normal(0, 10);
bLunarBright ~ normal(0, 10);
bDiurnal ~ normal(0, 10);
bLunarBrightDiurnal ~ normal(0, 10);
bJulianDay ~ normal(0, 10);
bJulianDaySquared ~ normal(0, 10);

// priors for varying effects
bPrecip_Species ~ normal(0, sigma_bPrecip_Species);
sigma_bPrecip_Species ~ cauchy(0, 1);
a_Species ~ normal(0, sigma_Species); 
sigma_Species ~ cauchy(0, 1);
}

///////////////////////////////////////////////////////////////////////////////
generated quantities {

vector[N] log_lik;
vector[N] theta;

for (i in 1:N) {

theta[i] <- 
a + a_Species[Species[i]] + 
bTrapEffort*TrapEffort[i] + bDayOnly*DayOnly[i] +
(bPrecip + bPrecip_Species[Species[i]])*Precip[i] +
bTmin*Tmin[i] + bLunarBright*LunarBright[i] + bDiurnal*Diurnal[i] +
bLunarBrightDiurnal*LunarBright[i]*Diurnal[i] +
bJulianDay*JulianDay[i] + bJulianDaySquared*JulianDaySquared[i];

log_lik[i] <- bernoulli_logit_log(Success[i], theta[i]);
}

}"