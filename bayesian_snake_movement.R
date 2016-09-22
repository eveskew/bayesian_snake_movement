# Snake movement analyses

# Evan Eskew
# Original: 15 June 2016
# Modified: 06 September 2016, 12 September 2016, 16 September 2016,
# 22 September 2016

# Data on snake movement (activity) are drawn from two long-term datasets:
# the Land-use Effects on Amphibian Populations (LEAP) study and monitoring
# of Ellenton Bay (EBay) conducted at the Savannah River Site.

# Code for the Bayesian models are organized at the END of this script. The
# preceding code calls fit model objects saved as .Rdata files.


library(dplyr)
library(ggplot2)
library(cowplot)
library(lubridate)
library(lme4)
library(rethinking)
library(rstan)
library(loo) 

#==============================================================================


# Define a custom precis plot function that allows you to alter x axis limits

custom_precis_plot <- 
  function (x, y, pars, cex, col.ci = "black", xlab = "Value", xlim...) {
    
    x <- x@output
    if (!missing(pars)) {
      x <- x[pars, ]
    }
    n <- nrow(x)
    mu <- x[n:1, 1]
    left <- x[[3]][n:1]
    right <- x[[4]][n:1]
    set_nice_margins()
    dotchart(mu, labels = rownames(x)[n:1], cex = cex, 
             xlab = xlab, xlim = xlim...)
    for (i in 1:length(mu)) lines(c(left[i], right[i]), c(i, i), 
                                  lwd = 2, col = col.ci)
    abline(v = 0, lty = 1, col = col.alpha("black", 0.15))
  }

#==============================================================================


# Import raw LEAP and EBay data


# Import LEAP effort and environmental data
leap_effort <- read.csv("data/LEAP effort and environment 2016-03-22.csv")

# Import LEAP capture data
leap_capture <- read.csv("data/LEAP snake fence captures 2016-03-09.csv")


# Import EBay effort and enviromental data
ebay_effort <- read.csv("data/ebay effort and environment 2016-03-22.csv")

# Import EBay capture data
ebay_capture <- read.csv("data/ebay snake fence captures 2016-03-07.csv")


# Import diurnal dataframe (binary variable indicating diurnal species)
diurnal <- read.csv("Data/diurnal.csv")

#==============================================================================


# Clean/modify LEAP data


# Add Julian day to the LEAP effort dataframe
tmp_dates <- as.POSIXlt(leap_effort$Date, format = "%m/%d/%y")
leap_effort$JulianDay <- tmp_dates$yday


# How many LEAP captures were recorded?
nrow(leap_capture)

# Create new LEAP capture dataframe summarizing species observations by counts
leap_capture2 <- summarise(group_by(leap_capture, Date, Bay, Species), n())
colnames(leap_capture2) <- c("Date", "Bay", "Species", "Count")

# Check to make sure captures are correct
sum(leap_capture2$Count)


# How many rows are in the LEAP effort dataframe?
nrow(leap_effort)

# Eliminate days that had no trapping effort
leap_effort2 <- filter(leap_effort, TrapEffort > 0)
nrow(leap_effort2)

# Eliminate days with no precipitation data
leap_effort2 <- filter(leap_effort2, !is.na(Precip))
nrow(leap_effort2)

# Replicate the effort dataframe for as many species as exist in the
# capture dataset
leap_effort2 <- merge(leap_effort2, levels(leap_capture2$Species), all = T)
colnames(leap_effort2)[9] <- "Species"
nrow(leap_effort2)


# Create a final LEAP dataframe by merging the effort and capture dataframes
d_leap <- merge(leap_effort2, leap_capture2, 
           by = c("Date", "Bay", "Species"), all = T)
d_leap$Count <- ifelse(is.na(d_leap$Count), 0, d_leap$Count)

# Does it have the correct number of total captures?
sum(d_leap$Count)


# Add a CapturesPerEffort and Success column
d_leap$CapturesPerEffort <- d_leap$Count/d_leap$TrapEffort
d_leap$Success <- ifelse(d_leap$Count == 0, 0, 1)

# Add standardized predictor columns
d_leap$TrapEffort.s <- as.numeric(scale(d_leap$TrapEffort))
d_leap$Precip.s <- as.numeric(scale(d_leap$Precip))
d_leap$Tmin.s <- as.numeric(scale(d_leap$Tmin))
d_leap$LunarBright.s <- as.numeric(scale(d_leap$LunarBright))
d_leap$JulianDay.s <- as.numeric(scale(d_leap$JulianDay))
d_leap$JulianDay.s.Squared <- d_leap$JulianDay.s*d_leap$JulianDay.s

# Add on a Year column and sort by Year, then JulianDay, then Species
tmp_dates <- as.POSIXlt(d_leap$Date, format = "%m/%d/%y")
d_leap$Year <- as.integer(format(tmp_dates, "%Y"))
d_leap <- arrange(d_leap, Year, JulianDay, Species)

# Drop unused levels
d_leap <- droplevels(d_leap)

#==============================================================================


# Clean/modify EBay data 


# Add Julian day to the EBay effort dataframe and clean up the Date column
tmp_dates <- as.POSIXlt(ebay_effort$Date, format = "%m/%d/%y")
ebay_effort$JulianDay <- tmp_dates$yday
ebay_effort$Date <- format(tmp_dates, format = "%m/%d/%y")

# Round temperatures to whole numbers
ebay_effort$Tmin <- round(ebay_effort$Tmin)
ebay_effort$Tmax <- round(ebay_effort$Tmax)


# Modify the ebay_capture Date column so that it matches with the ebay_effort
# dataframe
tmp_dates <- as.POSIXlt(ebay_capture$Date, format = "%d-%b-%y") 
ebay_capture$Date <- format(tmp_dates, format = "%m/%d/%y")


# Modify ebay_capture to rename the TrapType, Species, and Count columns
colnames(ebay_capture)[3] <- "TrapType"
colnames(ebay_capture)[4] <- "Species"
colnames(ebay_capture)[5] <- "Count"


# How many EBay captures were recorded?  
sum(ebay_capture$Count)

# Modify the ebay_capture data to combine different captures within days into
# one count
ebay_capture2 <- summarise(group_by(ebay_capture, Date, Species), sum(Count))
colnames(ebay_capture2) <- c("Date", "Species", "Count")

# Check to make sure captures are still correct
sum(ebay_capture2$Count)


# How many rows are in the EBay effort dataframe?
nrow(ebay_effort)

# Eliminate days with no precipitation data
ebay_effort2 <- filter(ebay_effort, !is.na(Precip))
nrow(ebay_effort2)

# Replicate the effort dataframe for as many species as exist in the
# capture dataset
ebay_effort2 <- merge(ebay_effort2, levels(ebay_capture$Species), all = T)
colnames(ebay_effort2)[14] <- "Species"
nrow(ebay_effort2)


# Create a final EBay dataframe by merging the effort and capture dataframes
d_ebay <- merge(ebay_effort2, ebay_capture2, 
                by = c("Date", "Species"), all = T)
d_ebay$Count <- ifelse(is.na(d_ebay$Count), 0, d_ebay$Count)

# Does it have the correct number of total captures?
sum(d_ebay$Count)


# Add a TrapEffort, CapturesPerEffort, and Success column
d_ebay$TrapEffort <- d_ebay$Buckets + d_ebay$CoffeeCan + d_ebay$SnakeTraps
d_ebay$CapturesPerEffort <- d_ebay$Count/d_ebay$TrapEffort
d_ebay$Success <- ifelse(d_ebay$Count == 0, 0, 1)

# Add a Diurnal column
d_ebay$Diurnal <- rep(NA, nrow(d_ebay))
for (i in 1:nrow(d_ebay)) {
  d_ebay$Diurnal[i] <- ifelse(filter(
    diurnal, Species == as.character(d_ebay$Species[i]))$Diurnal == 1, 1, 0)
}

# Add standardized predictor columns
d_ebay$TrapEffort.s <- as.numeric(scale(d_ebay$TrapEffort))
d_ebay$Precip.s <- as.numeric(scale(d_ebay$Precip))
d_ebay$Tmin.s <- as.numeric(scale(d_ebay$Tmin))
d_ebay$LunarBright.s <- as.numeric(scale(d_ebay$LunarBright))
d_ebay$JulianDay.s <- as.numeric(scale(d_ebay$JulianDay))
d_ebay$JulianDay.s.Squared <- d_ebay$JulianDay.s*d_ebay$JulianDay.s

# Add on a Year column and sort by Year, then JulianDay, then Species
tmp_dates <- as.POSIXlt(d_ebay$Date, format = "%m/%d/%y")
d_ebay$Year <- as.integer(format(tmp_dates, "%Y"))
d_ebay <- arrange(d_ebay, Year, JulianDay, Species)

# Drop "seminatrix" as a Species (because it has such an unusual life history
# compared to other species)
d_ebay <- filter(d_ebay, Species != "seminatrix")

# Drop unused levels
d_ebay <- droplevels(d_ebay)

#==============================================================================


# Define data for all LEAP Stan models
dat.leap = list(N = nrow(d_leap), 
                Count = d_leap$Count, 
                Success = d_leap$Success,
                TrapEffort = d_leap$TrapEffort.s,
                Precip = d_leap$Precip.s,
                Tmin = d_leap$Tmin.s,
                LunarBright = d_leap$LunarBright.s,
                JulianDay = d_leap$JulianDay.s,
                JulianDaySquared = d_leap$JulianDay.s.Squared,
                N_Species = max(as.integer(d_leap$Species)),
                Species = as.integer(d_leap$Species),
                N_Bay = max(as.integer(d_leap$Bay)),
                Bay = as.integer(d_leap$Bay))

# Define data for all EBay Stan models
dat.ebay = list(N = nrow(d_ebay), 
                Count = d_ebay$Count,
                Success = d_ebay$Success,
                TrapEffort = d_ebay$TrapEffort.s,
                DayOnly = d_ebay$DayOnly,
                Precip = d_ebay$Precip.s,
                Tmin = d_ebay$Tmin.s,
                LunarBright = d_ebay$LunarBright.s,
                JulianDay = d_ebay$JulianDay.s,
                JulianDaySquared = d_ebay$JulianDay.s.Squared,
                Diurnal = d_ebay$Diurnal,
                N_Species = max(as.integer(d_ebay$Species)), 
                Species = as.integer(d_ebay$Species))

#==============================================================================


# Rough plotting of raw count data

my_theme <- theme_minimal() +
  theme(text = element_text(size = 20),
        #plot.background = 
          #element_rect(fill = adjustcolor("floralwhite", alpha.f = 0.4)),
        panel.grid.major = element_line(size = 0.25, color = "gray77"),
        panel.grid.minor = element_line(size = 0.15, color = "gray77"),
        panel.grid.minor.x = element_line(size = 0),
        legend.title.align = 0.5)

plot1 <- filter(d_leap) %>%
ggplot(aes(x = JulianDay, y = Count)) + 
  scale_y_continuous(limits = c(0, 8), 
                     minor_breaks = seq(0, 10, 1), breaks = seq(0, 8, 2)) +
  xlab("Julian Day") + xlim(0, 365) + 
  ggtitle("LEAP Snake Captures") +
  geom_point(size = 3) +
  #geom_smooth(col = "red", size = 1, span = 1, method = "loess", se = F) +
  #facet_wrap(~ Species) +
  my_theme

plot2 <- filter(d_ebay) %>%
ggplot(aes(x = JulianDay, y = Count)) +
  scale_y_continuous(limits = c(0, 8), 
                     minor_breaks = seq(0, 10, 1), breaks = seq(0, 8, 2)) +
  xlab("Julian Day") + xlim(0, 365) +
  ggtitle("Ellenton Bay Snake Captures") +
  geom_point(size = 3) +
  #geom_smooth(col = "red", size = 1, span = 1, method = "loess", se = F) +
  #facet_wrap(~ Species) +
  my_theme

png("outputs/Fig1.png", width = 800, height = 1200)
plot_grid(plot1, plot2, nrow = 2, scale = 0.95,  
          labels = c("A", "B"), label_size = 30, vjust = 0.5)
dev.off()

#==============================================================================


# Load fit Bayesian models for LEAP data and generate summaries (tabular and
# visual) of parameter estimates

load("m1.nbinom.Rdata")
load("m2.nbinom.Rdata")
load("m3.nbinom.Rdata")
load("m4.nbinom.Rdata")
load("m5.nbinom.Rdata")


m1.nbinom.e <- extract(m1.nbinom, permuted = TRUE)
pars.trim <- names(m1.nbinom)[1:(length(names(m1.nbinom))-(dat.leap$N*2 + 1))]
m1.nbinom.df <- as.data.frame(m1.nbinom, pars = pars.trim)
precis(m1.nbinom.df, prob = 0.95)
plot(precis(m1.nbinom.df, prob = 0.95))

m2.nbinom.e <- extract(m2.nbinom, permuted = TRUE)
pars.trim <- names(m2.nbinom)[1:(length(names(m2.nbinom))-(dat.leap$N*2 + 1))]
m2.nbinom.df <- as.data.frame(m2.nbinom, pars = pars.trim)
precis(m2.nbinom.df, prob = 0.95)
plot(precis(m2.nbinom.df, prob = 0.95))

m3.nbinom.e <- extract(m3.nbinom, permuted = TRUE)
pars.trim <- names(m3.nbinom)[1:(length(names(m3.nbinom))-(dat.leap$N*2 + 1))]
m3.nbinom.df <- as.data.frame(m3.nbinom, pars = pars.trim)
precis(m3.nbinom.df, prob = 0.95)
plot(precis(m3.nbinom.df, prob = 0.95))

m4.nbinom.e <- extract(m4.nbinom, permuted = TRUE)
pars.trim <- names(m4.nbinom)[1:(length(names(m4.nbinom))-(dat.leap$N*2 + 1))]
m4.nbinom.df <- as.data.frame(m4.nbinom, pars = pars.trim)
precis(m4.nbinom.df, prob = 0.95)
plot(precis(m4.nbinom.df, prob = 0.95))

m5.nbinom.e <- extract(m5.nbinom, permuted = TRUE)
pars.trim <- names(m5.nbinom)[1:(length(names(m5.nbinom))-(dat.leap$N*2 + 1))]
m5.nbinom.df <- as.data.frame(m5.nbinom, pars = pars.trim)
precis(m5.nbinom.df, prob = 0.95)
plot(precis(m5.nbinom.df, prob = 0.95))


# Create a nice dotplot for m5 model parameter estimates

m5.precis <- precis(m5.nbinom.df, prob = 0.95)
# Remove the intercept parameter
m5.precis@output <- m5.precis@output[2:nrow(m5.precis@output), ]
rownames(m5.precis@output) <- c("Trap Effort", "Julian Day",
                               "Julian Day Squared", "Precipitation", 
                               "Temperature", "Moon Fraction", "Φ",
                               "σ (Var. Intercept by Species)", 
                               "σ (Var. Intercept by Location)")
custom_precis_plot(m5.precis, xlab = "Parameter Estimate", 
                   xlim = c(-2.5, 3.5), cex = 1.5)


# Generate WAIC estimates for all models and perform model comparison

m1.nbinom.log_lik <- extract_log_lik(m1.nbinom)
m1.nbinom.waic <- waic(m1.nbinom.log_lik)
print(m1.nbinom.waic)

m2.nbinom.log_lik <- extract_log_lik(m2.nbinom)
m2.nbinom.waic <- waic(m2.nbinom.log_lik)
print(m2.nbinom.waic)

m3.nbinom.log_lik <- extract_log_lik(m3.nbinom)
m3.nbinom.waic <- waic(m3.nbinom.log_lik)
print(m3.nbinom.waic)

m4.nbinom.log_lik <- extract_log_lik(m4.nbinom)
m4.nbinom.waic <- waic(m4.nbinom.log_lik)
print(m4.nbinom.waic)

m5.nbinom.log_lik <- extract_log_lik(m5.nbinom)
m5.nbinom.waic <- waic(m5.nbinom.log_lik)
print(m5.nbinom.waic)


m.comp <- as.data.frame(compare(m1.nbinom.waic, m2.nbinom.waic, 
                                m3.nbinom.waic, m4.nbinom.waic, 
                                m5.nbinom.waic))
m.comp$dWAIC <- round(m.comp$waic - min(m.comp$waic), digits = 5)
m.comp$wWAIC <- round(ICweights(m.comp$waic), digits = 5)
m.comp


# Load fit Bayesian models for EBay data and generate summaries (tabular and
# visual) of parameter estimates

load("m6.nbinom.Rdata")
load("m7.nbinom.Rdata")
load("m8.nbinom.Rdata")
load("m9.nbinom.Rdata")
load("m10.nbinom.Rdata")


m6.nbinom.e <- extract(m6.nbinom, permuted = TRUE)
pars.trim <- 
  names(m6.nbinom)[1:(length(names(m6.nbinom))-(dat.ebay$N*2 + 1))]
m6.nbinom.df <- as.data.frame(m6.nbinom, pars = pars.trim)
precis(m6.nbinom.df, prob = 0.95)
plot(precis(m6.nbinom.df, prob = 0.95))

m7.nbinom.e <- extract(m7.nbinom, permuted = TRUE)
pars.trim <- 
  names(m7.nbinom)[1:(length(names(m7.nbinom))-(dat.ebay$N*2 + 1))]
m7.nbinom.df <- as.data.frame(m7.nbinom, pars = pars.trim)
precis(m7.nbinom.df, prob = 0.95)
plot(precis(m7.nbinom.df, prob = 0.95))

m8.nbinom.e <- extract(m8.nbinom, permuted = TRUE)
pars.trim <- 
  names(m8.nbinom)[1:(length(names(m8.nbinom))-(dat.ebay$N*2 + 1))]
m8.nbinom.df <- as.data.frame(m8.nbinom, pars = pars.trim)
precis(m8.nbinom.df, prob = 0.95)
plot(precis(m8.nbinom.df, prob = 0.95))

m9.nbinom.e <- extract(m9.nbinom, permuted = TRUE)
pars.trim <- 
  names(m9.nbinom)[1:(length(names(m9.nbinom))-(dat.ebay$N*2 + 1))]
m9.nbinom.df <- as.data.frame(m9.nbinom, pars = pars.trim)
precis(m9.nbinom.df, prob = 0.95)
plot(precis(m9.nbinom.df, prob = 0.95))

m10.nbinom.e <- extract(m10.nbinom, permuted = TRUE)
pars.trim <- 
  names(m10.nbinom)[1:(length(names(m10.nbinom))-(dat.ebay$N*2 + 1))]
m10.nbinom.df <- as.data.frame(m10.nbinom, pars = pars.trim)
precis(m10.nbinom.df, prob = 0.95)
plot(precis(m10.nbinom.df, prob = 0.95))


# Create a nice dotplot for m10 model parameter estimates

m10.precis <- precis(m10.nbinom.df, prob = 0.95)
# Remove the intercept parameter
m10.precis@output <- m10.precis@output[2:nrow(m10.precis@output), ]
rownames(m10.precis@output) <- c("Trap Effort", "Traps Open Day Only?",
                                "Julian Day", "Julian Day Squared", 
                                "Precipitation", "Temperature", 
                                "Moon Fraction", "Diurnal Species?", 
                                "Moon-Diurnal Interaction", "Φ",
                                "σ (Var. Precip. Slope by Species)", 
                                "σ (Var. Intercept by Species)")
custom_precis_plot(m10.precis, xlab = "Parameter Estimate", 
                   xlim = c(-2.5, 3.5), cex = 1.5)


# Create a dotplot for m10 model species-specific varying effects
# of precipitation

m10.precis.precip <- precis(m10.nbinom.df[12:31], depth = 2, prob = 0.95)
m10.precis.precip@output$names <- 
  c("Farancia abacura", "Opheodrys aestivus", 
    "Pantherophis alleghaniensis", "Cemophora coccinea", 
    "Coluber constrictor", "Agkistrodon contortrix", 
    "Storeria dekayi", "Diadophis punctatus", 
    "Nerodia erythrogaster", "Farancia erytrogramma", 
    "Nerodia fasciata", "Masticophis flagellum", 
    "Nerodia floridana", "Pantherophis guttatus", 
    "Crotalus horridus", "Storeria occipitomaculata",
    "Agkistrodon piscivorus", "Heterodon platirhinos", 
    "Thamnophis sauritus", "Thamnophis sirtalis")
m10.precis.precip@output <- arrange(m10.precis.precip@output, Mean)
rownames(m10.precis.precip@output) <- m10.precis.precip@output$names
custom_precis_plot(m10.precis.precip, xlab = "Parameter Estimate", 
                   cex = 1.5, xlim = c(-1.5, 1.5))


# Create a dotplot for m10 model species-specific preciptation intercepts
# (i.e., add together the overall precipitation estimate and the 
# species-specific varying effects to get a realized species-specific effect)

for (i in 1:20) {
  m10.nbinom.df[paste0("new_bPrecip_", i)] <-
    m10.nbinom.df$bPrecip + m10.nbinom.df[11 + i]
}

m10.precis.precip <- precis(m10.nbinom.df[54:73], depth = 2, prob = 0.95)
m10.precis.precip@output$names <- 
  c("Farancia abacura", "Opheodrys aestivus", 
    "Pantherophis alleghaniensis", "Cemophora coccinea", 
    "Coluber constrictor", "Agkistrodon contortrix", 
    "Storeria dekayi", "Diadophis punctatus", 
    "Nerodia erythrogaster", "Farancia erytrogramma", 
    "Nerodia fasciata", "Masticophis flagellum", 
    "Nerodia floridana", "Pantherophis guttatus", 
    "Crotalus horridus", "Storeria occipitomaculata",
    "Agkistrodon piscivorus", "Heterodon platirhinos", 
    "Thamnophis sauritus", "Thamnophis sirtalis")
m10.precis.precip@output <- arrange(m10.precis.precip@output, Mean)
rownames(m10.precis.precip@output) <- m10.precis.precip@output$names
custom_precis_plot(m10.precis.precip, xlab = "Parameter Estimate", 
                   cex = 1.5, xlim = c(-2, 1))
abline(v = mean(m10.nbinom.df$bPrecip), lty = 3)


# Generate WAIC estimates for all models and perform model comparison

m6.nbinom.log_lik <- extract_log_lik(m6.nbinom)
m6.nbinom.waic <- waic(m6.nbinom.log_lik)
print(m6.nbinom.waic)

m7.nbinom.log_lik <- extract_log_lik(m7.nbinom)
m7.nbinom.waic <- waic(m7.nbinom.log_lik)
print(m7.nbinom.waic)

m8.nbinom.log_lik <- extract_log_lik(m8.nbinom)
m8.nbinom.waic <- waic(m8.nbinom.log_lik)
print(m8.nbinom.waic)

m9.nbinom.log_lik <- extract_log_lik(m9.nbinom)
m9.nbinom.waic <- waic(m9.nbinom.log_lik)
print(m9.nbinom.waic)

m10.nbinom.log_lik <- extract_log_lik(m10.nbinom)
m10.nbinom.waic <- waic(m10.nbinom.log_lik)
print(m10.nbinom.waic)


m.comp <- as.data.frame(compare(m6.nbinom.waic, m7.nbinom.waic, 
                                m8.nbinom.waic, m9.nbinom.waic, 
                                m10.nbinom.waic))
m.comp$dWAIC <- round(m.comp$waic - min(m.comp$waic), digits = 5)
m.comp$wWAIC <- round(ICweights(m.comp$waic), digits = 5)
m.comp

#==============================================================================


# Predictive checks


# m5 results with varying temperature
set.seed(4)
n.samples <- 15000

# Predictions for Tmin = 0
arbitrary.temp <- (0 - mean(d_leap$Tmin))/(sd(d_leap$Tmin))

preds.Tmin0 <- 
  rnbinom(n.samples, 
          mu = exp(m5.nbinom.e$a + 
                     m5.nbinom.e$bTmin*arbitrary.temp), 
          size = m5.nbinom.e$phi)

# Predictions for Tmin = 12
arbitrary.temp <- (12 - mean(d_leap$Tmin))/(sd(d_leap$Tmin))

preds.Tmin12 <- 
  rnbinom(n.samples, 
          mu = exp(m5.nbinom.e$a +
                     m5.nbinom.e$bTmin*arbitrary.temp), 
          size = m5.nbinom.e$phi)

# Predictions for Tmin = 24
arbitrary.temp <- (24 - mean(d_leap$Tmin))/(sd(d_leap$Tmin))

preds.Tmin24 <-
  rnbinom(n.samples, 
          mu = exp(m5.nbinom.e$a +
                     m5.nbinom.e$bTmin*arbitrary.temp), 
          size = m5.nbinom.e$phi)

# Package predictions into a dataframe
preds.Temp <- as.data.frame(c(preds.Tmin0, preds.Tmin12, preds.Tmin24))
colnames(preds.Temp) <- "Count"
preds.Temp$Treatment <- rep(c("0", "12", "24"), each = n.samples)

# Plot
ggplot(preds.Temp, aes(x = Count, color = Treatment)) +
  geom_density() + coord_cartesian(xlim = c(0, 3))

ggplot(preds.Temp, aes(x = Treatment, y = Count)) +
  geom_jitter() +
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", 
               size= 0.5, color = "red", geom = "crossbar")


# m5 results with varying Julian day
set.seed(4)
n.samples <- 15000

# Predictions for Julian Day 50
arbitrary.day <- (50 - mean(d_leap$JulianDay))/(sd(d_leap$JulianDay))

preds.JulianDay50 <-
    rnbinom(n.samples, 
            mu = exp(m5.nbinom.e$a +
                       m5.nbinom.e$bJulianDay*arbitrary.day +
                       m5.nbinom.e$bJulianDaySquared*(arbitrary.day^2)), 
            size = m5.nbinom.e$phi)

# Predictions for Julian Day 150
arbitrary.day <- (150 - mean(d_leap$JulianDay))/(sd(d_leap$JulianDay))

preds.JulianDay150 <-
    rnbinom(n.samples, 
            mu = exp(m5.nbinom.e$a +
                       m5.nbinom.e$bJulianDay*arbitrary.day +
                       m5.nbinom.e$bJulianDaySquared*(arbitrary.day^2)), 
            size = m5.nbinom.e$phi)

# Predictions for Julian Day 250
arbitrary.day <- (250 - mean(d_leap$JulianDay))/(sd(d_leap$JulianDay))

preds.JulianDay250 <- 
    rnbinom(n.samples, 
            mu = exp(m5.nbinom.e$a +
                       m5.nbinom.e$bJulianDay*arbitrary.day +
                       m5.nbinom.e$bJulianDaySquared*(arbitrary.day^2)), 
            size = m5.nbinom.e$phi)

# Package predictions into a dataframe
preds.Julian <- 
  as.data.frame(c(preds.JulianDay50, preds.JulianDay150, preds.JulianDay250))
colnames(preds.Julian) <- "Count"
preds.Julian$Treatment <- rep(c("50", "150", "250"), each = n.samples)

# Plot
ggplot(preds.Julian, aes(x = Count, color = Treatment)) +
  geom_density() + coord_cartesian(xlim = c(0, 3))

ggplot(preds.Julian, aes(x = Treatment, y = Count)) +
  geom_jitter() +
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", 
               size= 0.5, color = "red", geom = "crossbar")


# m5 results with varying temperature and varying precipitation
set.seed(4)
n.samples <- 15000

# Predictions for Tmin = -3, precip = 24
arbitrary.temp <- (-3 - mean(d_leap$Tmin))/(sd(d_leap$Tmin))
arbitrary.precip <- (24 - mean(d_leap$Precip))/(sd(d_leap$Precip))

preds.low <-
    rnbinom(n.samples, 
            mu = exp(m5.nbinom.e$a +
                       m5.nbinom.e$bTmin*arbitrary.temp +
                       m5.nbinom.e$bPrecip*arbitrary.precip), 
            size = m5.nbinom.e$phi)

# Predictions for Tmin = 12, precip = 4
arbitrary.temp <- (12 - mean(d_leap$Tmin))/(sd(d_leap$Tmin))
arbitrary.precip <- (4 - mean(d_leap$Precip))/(sd(d_leap$Precip))

preds.med <- 
    rnbinom(n.samples, 
            mu = exp(m5.nbinom.e$a +
                       m5.nbinom.e$bTmin*arbitrary.temp +
                       m5.nbinom.e$bPrecip*arbitrary.precip), 
            size = m5.nbinom.e$phi)

# Predictions for Tmin = 24, precip = 0
arbitrary.temp <- (24 - mean(d_leap$Tmin))/(sd(d_leap$Tmin))
arbitrary.precip <- (0 - mean(d_leap$Precip))/(sd(d_leap$Precip))

preds.high <- 
    rnbinom(n.samples, 
            mu = exp(m5.nbinom.e$a +
                       m5.nbinom.e$bTmin*arbitrary.temp +
                       m5.nbinom.e$bPrecip*arbitrary.precip), 
            size = m5.nbinom.e$phi)

# Package predictions into a dataframe
preds.Temp.Precip <- as.data.frame(c(preds.low, preds.med, preds.high))
colnames(preds.Temp.Precip) <- "Count"
preds.Temp.Precip$Treatment <- rep(c("low", "med", "high"), each = n.samples)
preds.Temp.Precip$Treatment <-
  factor(preds.Temp.Precip$Treatment, 
         levels = c("low", "med", "high"))

# Plot
ggplot(preds.Temp.Precip, aes(x = Count, color = Treatment)) +
  geom_density() + coord_cartesian(xlim = c(0, 3))

ggplot(preds.Temp.Precip, aes(x = Treatment, y = Count)) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 24, 4)) +
  xlab("Predictive Scenario") +
  geom_jitter(width = 0.7, height = 0.2) +
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", 
               size= 0.6, color = "grey", geom = "crossbar") +
  scale_x_discrete(labels = c("low" = "-3°C, 24 mm",
                              "med" = "12°C, 4 mm",
                              "high" = "24°C, 0 mm")) +
  my_theme


# m10 results with varying lunar brightness and diurnality
set.seed(4)
n.samples <- 15000

# Predictions for Lunar = 0, Diurnal = 0
arbitrary.lunar <- (0 - mean(d_ebay$LunarBright))/(sd(d_leap$LunarBright))
arbitrary.diurnal <- 0

preds.1 <-
  rnbinom(n.samples, 
          mu = exp(m10.nbinom.e$a +
                     m10.nbinom.e$bLunarBright*arbitrary.lunar +
                     m10.nbinom.e$bDiurnal*arbitrary.diurnal +
                     m10.nbinom.e$bLunarBrightDiurnal*
                     arbitrary.lunar*arbitrary.diurnal), 
          size = m10.nbinom.e$phi)

# Predictions for Lunar = 0.5, Diurnal = 0
arbitrary.lunar <- (0.5 - mean(d_ebay$LunarBright))/(sd(d_leap$LunarBright))
arbitrary.diurnal <- 0

preds.2 <-
  rnbinom(n.samples, 
          mu = exp(m10.nbinom.e$a +
                     m10.nbinom.e$bLunarBright*arbitrary.lunar +
                     m10.nbinom.e$bDiurnal*arbitrary.diurnal +
                     m10.nbinom.e$bLunarBrightDiurnal*
                     arbitrary.lunar*arbitrary.diurnal), 
          size = m10.nbinom.e$phi)

# Predictions for Lunar = 1, Diurnal = 0
arbitrary.lunar <- (1 - mean(d_ebay$LunarBright))/(sd(d_leap$LunarBright))
arbitrary.diurnal <- 0

preds.3 <-
  rnbinom(n.samples, 
          mu = exp(m10.nbinom.e$a +
                     m10.nbinom.e$bLunarBright*arbitrary.lunar +
                     m10.nbinom.e$bDiurnal*arbitrary.diurnal +
                     m10.nbinom.e$bLunarBrightDiurnal*
                     arbitrary.lunar*arbitrary.diurnal), 
          size = m10.nbinom.e$phi)

# Predictions for Lunar = 0, Diurnal = 1
arbitrary.lunar <- (0 - mean(d_ebay$LunarBright))/(sd(d_leap$LunarBright))
arbitrary.diurnal <- 1

preds.4 <-
  rnbinom(n.samples, 
          mu = exp(m10.nbinom.e$a +
                     m10.nbinom.e$bLunarBright*arbitrary.lunar +
                     m10.nbinom.e$bDiurnal*arbitrary.diurnal +
                     m10.nbinom.e$bLunarBrightDiurnal*
                     arbitrary.lunar*arbitrary.diurnal), 
          size = m10.nbinom.e$phi)

# Predictions for Lunar = 0.5, Diurnal = 1
arbitrary.lunar <- (0.5 - mean(d_ebay$LunarBright))/(sd(d_leap$LunarBright))
arbitrary.diurnal <- 1

preds.5 <-
  rnbinom(n.samples, 
          mu = exp(m10.nbinom.e$a +
                     m10.nbinom.e$bLunarBright*arbitrary.lunar +
                     m10.nbinom.e$bDiurnal*arbitrary.diurnal +
                     m10.nbinom.e$bLunarBrightDiurnal*
                     arbitrary.lunar*arbitrary.diurnal), 
          size = m10.nbinom.e$phi)

# Predictions for Lunar = 1, Diurnal = 1
arbitrary.lunar <- (1 - mean(d_ebay$LunarBright))/(sd(d_leap$LunarBright))
arbitrary.diurnal <- 1

preds.6 <-
  rnbinom(n.samples, 
          mu = exp(m10.nbinom.e$a +
                     m10.nbinom.e$bLunarBright*arbitrary.lunar +
                     m10.nbinom.e$bDiurnal*arbitrary.diurnal +
                     m10.nbinom.e$bLunarBrightDiurnal*
                     arbitrary.lunar*arbitrary.diurnal), 
          size = m10.nbinom.e$phi)

# Package predictions into a dataframe
preds.Lunar <- 
  as.data.frame(c(preds.1, preds.2, preds.3, preds.4, preds.5, preds.6))
colnames(preds.Lunar) <- "Count"
preds.Lunar$Treatment <- 
  factor(rep(c("New Moon, Nocturnal", "Halfmoon, Nocturnal",
               "Full moon, Nocturnal", "New moon, Diurnal",
               "Halfmoon, Diurnal", "Full moon, Diurnal"),
             each = n.samples),
         levels = c("New Moon, Nocturnal", "Halfmoon, Nocturnal",
                    "Full moon, Nocturnal", "New moon, Diurnal",
                    "Halfmoon, Diurnal", "Full moon, Diurnal"))

# Plot
ggplot(preds.Lunar, aes(x = Count, color = Treatment)) +
  geom_density() + coord_cartesian(xlim = c(0, 3))

ggplot(preds.Lunar, aes(x = Treatment, y = Count)) +
  geom_jitter() +
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", 
               size= 0.5, color = "red", geom = "crossbar")

#==============================================================================


# Linear modeling with lme4


# Here I'm fitting m1.nbinom and m5.nbinom using "glmer()" from lme4 with
# a Poisson family outcome. Even though these actually have a different 
# outcome distribution and different method of fitting, in general 
# parameter estimates from this method and Stan are quite similar, suggesting 
# our inferences are robust.

m1.pois <- glmer(Count ~ TrapEffort.s + JulianDay.s + JulianDay.s.Squared + 
              (1|Species) + (1|Bay), 
            data = d_leap, family = poisson)

m5.pois <- glmer(Count ~ TrapEffort.s + JulianDay.s + JulianDay.s.Squared +
              Precip.s + Tmin.s + LunarBright.s +
              (1|Species) + (1|Bay), 
            data = d_leap, family = poisson)


precis(m1.pois, prob = 0.95)
ranef(m1.pois)
precis(m1.nbinom.df, prob = 0.95, depth = 2)

precis(m5.pois, prob = 0.95)
ranef(m5.pois)
precis(m5.nbinom.df, prob = 0.95, depth = 2)

#==============================================================================


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

m1.nbinom <- stan(model_code = m1.nbinom.modelcode, 
                  model_name = "m1.nbinom", data = dat.leap, 
                  iter = 2000, warmup = 1000, chains = 3, 
                  cores = 3, sample_file = "m1.nbinom.csv", verbose = TRUE)
save("m1.nbinom", file = "m1.nbinom.Rdata")

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

m2.nbinom <- stan(model_code = m2.nbinom.modelcode, 
              model_name = "m2.nbinom", data = dat.leap, 
              iter = 2000, warmup = 1000, chains = 3, 
              cores = 3, sample_file = "m2.nbinom.csv", verbose = TRUE)
save("m2.nbinom", file = "m2.nbinom.Rdata")

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

m3.nbinom <- stan(model_code = m3.nbinom.modelcode, 
              model_name = "m3.nbinom", data = dat.leap, 
              iter = 2000, warmup = 1000, chains = 3, 
              cores = 3, sample_file = "m3.nbinom.csv", verbose = TRUE)
save("m3.nbinom", file = "m3.nbinom.Rdata")

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

m4.nbinom <- stan(model_code = m4.nbinom.modelcode, 
              model_name = "m4.nbinom", data = dat.leap, 
              iter = 2000, warmup = 1000, chains = 3, 
              cores = 3, sample_file = "m4.nbinom.csv", verbose = TRUE)
save("m4.nbinom", file = "m4.nbinom.Rdata")

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

m5.nbinom <- stan(model_code = m5.nbinom.modelcode, 
              model_name = "m5.nbinom", data = dat.leap, 
              iter = 2000, warmup = 1000, chains = 3, 
              cores = 3, sample_file = "m5.nbinom.csv", verbose = TRUE)
save("m5.nbinom", file = "m5.nbinom.Rdata")

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

m6.nbinom <- stan(model_code = m6.nbinom.modelcode, 
                   model_name = "m6.nbinom", data = dat.ebay, 
                   iter = 2000, warmup = 1000, chains = 3, 
                   cores = 3, sample_file = "m6.nbinom.csv", verbose = TRUE)
save("m6.nbinom", file = "m6.nbinom.Rdata")

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

m7.nbinom <- stan(model_code = m7.nbinom.modelcode, 
              model_name = "m7.nbinom", data = dat.ebay, 
              iter = 2000, warmup = 1000, chains = 3, 
              cores = 3, sample_file = "m7.nbinom.csv", verbose = TRUE)
save("m7.nbinom", file = "m7.nbinom.Rdata")

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

m8.nbinom <- stan(model_code = m8.nbinom.modelcode, 
              model_name = "m8.nbinom", data = dat.ebay, 
              iter = 2000, warmup = 1000, chains = 3, 
              cores = 3, sample_file = "m8.nbinom.csv", verbose = TRUE)
save("m8.nbinom", file = "m8.nbinom.Rdata")

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

m9.nbinom <- stan(model_code = m9.nbinom.modelcode, 
              model_name = "m9.nbinom", data = dat.ebay, 
              iter = 2000, warmup = 1000, chains = 3, 
              cores = 3, sample_file = "m9.nbinom.csv", verbose = TRUE)
save("m9.nbinom", file = "m9.nbinom.Rdata")

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

m10.nbinom <- stan(model_code = m10.nbinom.modelcode, 
              model_name = "m10.nbinom", data = dat.ebay, 
              iter = 2000, warmup = 1000, chains = 3, 
              cores = 3, sample_file = "m10.nbinom.csv", verbose = TRUE)
save("m10.nbinom", file = "m10.nbinom.Rdata")

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

m10.binom <- stan(model_code = m10.binom.modelcode, 
              model_name = "m10.binom", data = dat.ebay, 
              iter = 2000, warmup = 1000, chains = 3, 
              cores = 3, sample_file = "m10.binom.csv", verbose = TRUE)
save("m10.binom", file = "m10.binom.Rdata")

load("m10.binom.Rdata")

m10.binom.e <- extract(m10.binom, permuted = TRUE)
pars.trim <- names(m10.binom)[1:(length(names(m10.binom))-(dat.ebay$N*2 + 1))]
m10.binom.df <- as.data.frame(m10.binom, pars = pars.trim)
precis(m10.binom.df, prob = 0.95)
plot(precis(m10.binom.df, prob = 0.95))