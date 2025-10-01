# Exploring WNV data 
# Seroprevalence study on humans in Algeria, 2017 - 2018

library(dplyr)
library(tidyverse)
library(ggplot2)
library(Rsero)
library(patchwork)

wnv <- read.csv("Data/WNV seroprevalence study_Aissam_IP DZ samp dates.csv")

table(wnv$sex)
summary(wnv$Age..y.)
table(wnv$Region, exclude = F)
table(wnv$District)
table(wnv$Residency.type)
table(wnv$YF.vaccine)

rural <- filter(wnv, Residency.type == "Rural")
table(rural$Region)

# Outcomes
table(wnv$IgG.anti.WNV..EIA., exclude = F)
table(wnv$PRNT90, exclude = F)
table(wnv$PRNT90.final, exclude = F)

table(wnv$PRNT90.final, exclude = F)[2]/table(wnv$PRNT90.final, exclude = F)[1]

wnv$res <- ifelse(wnv$PRNT90.final == "Positif", TRUE, FALSE)

# Remove obs with missing age
wnv <- wnv |>
  filter(!is.na(Age..y.)) |>
  rename(age = Age..y.)

# Sample times
wnv$sample_year <- as.numeric(gsub("\\D", "", wnv$Collection.date))

table(wnv$sample_year)

xtabs(~District+sample_year, data=wnv)


# Separate by region
south <- filter(wnv, Region == "Sud")
north <- filter(wnv, Region == "Nord")
oran <- filter(wnv, District == "Oran")
jijel <- filter(wnv, District == "Jijel")
tizi <- filter(wnv, District == "Alger")
biskra <- filter(wnv, District == "Biskra")
tim <- filter(wnv, District == "Adrar")

table(north$District)
table(south$District)

#Use Rsero R package to analyze data (change "biskra" to whichever province or region of interest)
wnv_survey  = SeroData(age_at_sampling = biskra$age,
                       Y               = biskra$res, # seropositive or seronegative
                       sampling_year   = biskra$sample_year)

seroprevalence(wnv_survey)
seroprevalence.plot(wnv_survey)

# Estimate the FOI using different models
model = list()
model[[1]] = FOImodel(type='constant')
model[[2]] = FOImodel(type='constant', seroreversion = 1)
model[[3]] = FOImodel(type='constantoutbreak')
model[[4]] = FOImodel(type='outbreak', K=1)
model[[5]] = FOImodel(type='emergence')
#model[[6]] = FOImodel(type='outbreak', prior_distribution_T = "exponential", priorT1 = 0.0001, priorT2 = 10)
model[[6]] = FOImodel(type='outbreak', prior_distribution_T = "exponential", priorT1 = 0.1)
model[[7]] = FOImodel(type='outbreak', K=1, priorT1 = 35, priorT2 = 10)

#model[[6]] = FOImodel(type='constantoutbreak', sp=0.99)
#model[[7]] = FOImodel(type='constantoutbreak', sp=0.95)
#model[[8]] = FOImodel(type='outbreak', K=1, sp=0.99)
#model[[9]] = FOImodel(type='outbreak', K=1, sp=0.95)
#model[[6]] = FOImodel(type='constantoutbreak', K=2,sp=0.99)
#model[[7]] = FOImodel(type='outbreak', K=1, sp = 0.99,priorT1 = 20,priorT2=1)
#model[[8]] = FOImodel(type='constantoutbreak', K=1, sp = 0.99,priorT1 = 20,priorT2=1)

#model_age_dependent=list()
#model_age_dependent[[1]] = FOImodel(type='constant',age_dependent_foi = 1)
#model_age_dependent[[2]] = FOImodel(type='constant', seroreversion = 1,age_dependent_foi = 1)
#model_age_dependent[[3]] = FOImodel(type='constantoutbreak',age_dependent_foi = 1)
#model_age_dependent[[4]] = FOImodel(type='outbreak', K=1, sp = 0.99,age_dependent_foi = 1)

#model[[5]] = FOImodel(type= "piecewise", K=2, seroreversion=0) 
#model[[6]] = FOImodel(type='outbreak', K=1, background = 1, priorT1 = 20, priorT2 = 40)
#model[[7]] = FOImodel(type='outbreak', K=1, sp = 0.98)
#model[[8]] = FOImodel(type='outbreak', K=2)
#model[[9]] = FOImodel(type='outbreak', K=4)
#model[[10]] = FOImodel(type='emergence')

tictoc::tic()
Fit = fit(data = wnv_survey,
          model = model[[2]], # Choose the model from the list here or create a loop to run all of them
          chain=4,
          cores=4, 
          iter=5000, 
          warmup = 2500,  
          thin = 1)
#control = list(adapt_delta = 0.99)) # 0.95
# max_treedepth = 15))

tictoc::toc()
beepr::beep()

Chains=rstan::extract(Fit$fit)
plot(Chains$T_raw)
#saveRDS(Fit, "Results/chains new Rsero updated/oran_outbreak.RDS")
#saveRDS(Fit, "Results/chains new Rsero updated/tizi_outbreak.RDS")
#saveRDS(Fit, "Results/chains new Rsero updated/jijel_outbreak.RDS")
#saveRDS(Fit, "Results/chains new Rsero updated/biskra_constant.RDS")
#saveRDS(Fit, "Results/chains new Rsero updated/tim_constant.RDS")

#saveRDS(Fit, "Results/chains new Rsero updated/oran_outbreak_T_exp_prior.RDS")
#saveRDS(Fit, "Results/chains new Rsero updated/tizi_outbreak_T_exp_prior.RDS")


#library("shinystan") 
#launch_shinystan(Fit$fit)

C = compute_information_criteria(Fit) 
#C$DIC # save the DIC and effective number of parameters 
#C$pD
C$WAIC
C$PSIS_LOO

# If some Pareto k values are too high, use this modified function from the package
# to do moment matching to increase the accuracy of the sampling
# Note that you need to do this in the same session in which the model was run
# not on saved results
source("compute_info_criteria2_function.R")
C = compute_information_criteria2(Fit) 
C$PSIS_LOO

p = seroprevalence.fit(Fit)

# Plot the model fit
i=1
p[[1]] + 
  ylim(0,1)

Fit$model # print parameters and their priors
print(Fit$fit, digits = 4) #
print(Fit$fit, digits = 4, pars = c("annual_foi")) # Need to print more decimals
print(Fit$fit, digits = 4, pars = c("annual_foi","rho")) # Need to print more decimals
print(Fit$fit, digits = 4, pars = c("T[1]","alpha[1]","annual_foi")) # Need to print more decimals
print(Fit$fit, digits = 4, pars = c("T[1]","alpha[1]")) # Need to print more decimals
print(Fit$fit, digits = 4, pars = c("annual_foi","T")) # Need to print more decimals

plot_posterior(Fit)

rstan::traceplot(Fit$fit, pars = c("annual_foi"))
rstan::traceplot(Fit$fit, pars = c("annual_foi","rho"))
rstan::traceplot(Fit$fit, pars = c("T[1]","alpha[1]","annual_foi"))
rstan::traceplot(Fit$fit, pars = c("T[1]","alpha[1]"))
rstan::traceplot(Fit$fit, pars = c("T","annual_foi"))

parameters_credible_intervals(Fit)# Other diagnostics
stan_ess(Fit$fit) # show the ratio of effective sample size over total sample size
stan_diag(Fit$fit, info = 'sample') # shows diagnostic (loglikelihood, acceptance rate)

pairs(Fit$fit, pars = c("T[1]", "alpha[1]"))




######## NEW SIMPLE MODEL ##############
library(rstan)
# 
# rstan_options(auto_write = TRUE)
# options(mc.cores = max(1L, parallel::detectCores() - 1L))
# 
# # ---- Load data ----
# # Expect a CSV file "data.csv" with a column named Y of 0/1 values.
# # If you already have a vector/data.frame in R, just set `Y <- as.integer(your_vector)` and skip the read.csv.
# df <- read.csv("data.csv", stringsAsFactors = FALSE)
# if (!"Y" %in% names(df)) stop("The CSV must contain a column named 'Y' with 0/1 values.")
# 
# Y <- as.integer(df$Y)
# if (any(!Y %in% c(0L, 1L))) stop("Column 'Y' must contain only 0/1 values.")


N = wnv_survey$N
Y = wnv_survey$Y
stan_data <- list(N = length(Y), Y = Y, alpha=1, beta=1)



# ---- Compile & Sample ----
stan_file <- "simple_model.stan"
sm <- stan_model(file = stan_file)

fit <- sampling(
  sm,
  data = stan_data,
  seed = 1234,
  chains = 4,
  iter = 5000,
  warmup = 2500,
  refresh = 200
)

print(fit, pars = c("annual_foi_raw"))
rstan::traceplot(fit, pars = c("annual_foi_raw"))

#saveRDS(fit, "Results/chains new Rsero updated/tizi_simple.RDS")
#saveRDS(fit, "Results/chains new Rsero updated/oran_simple.RDS")
#saveRDS(fit, "Results/chains new Rsero updated/jijel_simple.RDS")
#saveRDS(fit, "Results/chains new Rsero updated/biskra_simple.RDS")
#saveRDS(fit, "Results/chains new Rsero updated/tim_simple.RDS")


# ---- Compute DIC ----
# DIC = 2 * E[D(theta) | y] - D(E[theta | y]),
# where D(theta) = -2 * log p(y | theta) (constants in the likelihood cancel).
post <- rstan::extract(fit)
theta_draws <- as.numeric(post$annual_foi_raw)

# Per-draw log-likelihood:
# log p(y | theta) = sum_j [ y_j * log(theta) + (1 - y_j) * log(1 - theta) ]
loglik_draw <- vapply(theta_draws, function(th) {
  if (th <= 0 || th >= 1) return(-Inf) # safety, should be rare due to bounds
  sum(Y * log(th) + (1 - Y) * log1p(-th))
}, numeric(1))

D_draw <- -2 * loglik_draw
Dbar <- mean(D_draw, na.rm = TRUE)

# Posterior mean of theta
theta_hat <- mean(theta_draws, na.rm = TRUE)

# Deviance at posterior mean
loglik_hat <- sum(Y * log(theta_hat) + (1 - Y) * log1p(-theta_hat))
D_hat <- -2 * loglik_hat

DIC <- 2 * Dbar - D_hat

cat("\n----- DIC Summary -----\n")
cat(sprintf("E[Deviance | y] (Dbar): %.3f\n", Dbar))
cat(sprintf("Deviance at E[theta|y] (Dhat): %.3f\n", D_hat))
cat(sprintf("DIC: %.3f\n", DIC))
cat("-----------------------\n")

# Optional: effective number of parameters (pD)
pD <- Dbar - D_hat
cat(sprintf("Effective number of parameters (pD): %.3f\n", pD))

log_lik1 <- extract_log_lik(fit)
waic1 <- loo::waic(log_lik1)
waic1

loo1 <- loo::loo(log_lik1)
print(loo1)
plot(loo1)

p = seroprevalence.fit(fit)

# Plot the model fit
i=1
p[[1]] + 
  ylim(0,1)



