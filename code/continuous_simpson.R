# set working directory to file location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(brms)       # Bayesian regression models
library(HDInterval) # for calculating highest density intervals
library(tidybayes)  # for tidy output of Bayesian models
# for data wrangling 
library(tibble)     # dataframes
library(tidyr)      # for tidying dataframes
library(dplyr)      # intuitive data manipulation
# for plotting
library(ggplot2)    # for plots

source("utils.R")

# seeding for reproducibility
set.seed(123)

### global variables
sample_size = 700 # sample size of simulated continuous data
n_iter = 5000     # number of iterations 


##### DATA SIMULATION

### Experiment 1: Gender as confound

testosterone_sample <- function(n, pi = 0.5) {
  # Generate mixing component: 1 = N1, 0 = N2
  component <- rbinom(n, 1, pi)
  
  # Generate values based on component
  x <- numeric(n)
  x[component == 1] <- rnorm(sum(component == 1), mean = 1.2, sd = 0.4)
  x[component == 0] <- rnorm(sum(component == 0), mean = 19.85, sd = 5.525)
  
  return(x)
}

G <- testosterone_sample(sample_size)

drug_sample <- function(n, G) {
  # Set parameters (determined by trying to replicate the probabilities in the
  # original Simpson's paradox data)
  alpha <- 1.30101
  beta_G <- -0.123612
  
  p <- 1 / (1 + exp(-(alpha + beta_G * G)))
  return(rbinom(n, 1, p))
}

D_g <- drug_sample(sample_size, G)

# Simulate recovery
recovery_sample <- function(n, D, G) {
  
  # Set parameters (determined by trying to replicate the probabilities in the
  # original Simpson's paradox data)
  alpha <- 0.679958
  beta_D <- 0.5
  beta_G <- 0.0615114
  
  p <- 1 / (1 + exp(-(alpha + beta_D * D + beta_G * G)))
  return(rbinom(n, 1, p))
}

R_g <- recovery_sample(sample_size, D_g, G)

data_gender <- tibble(
  gender = G,
  drug = ifelse(D_g == 1, "Take", "Refuse"),
  recover = R_g == 1
)

# summarize data to compare fractions to original Simpson's paradox data
gender_summary <- data_gender %>%
  mutate(
    gender_group = ifelse(gender < 3.0, "female", "male")
  ) %>%
  group_by(drug, gender_group) %>%
  summarise(
    recovery_rate = mean(recover),
    .groups = "drop"
  )

# standardize data to make it easier to model
data_gender <- data_gender %>%
  mutate(gender_std = (gender - mean(gender)) / sd(gender))

### Experiment 2: Blood pressure as mediator

# Simulate drug assignment
drug_sample_b <- function(n, pi=0.5) {
  return(rbinom(n, 1, 0.5))
}

D_b <- drug_sample_b(sample_size)

# Simulate blood pressure
bloodpressure_sample <- function(n, D) {
  #return(rnorm(n, mean = 71 + 7 * D, sd = 15))
  x <- numeric(n)
  x[D == 1] <- rnorm(sum(D == 1), mean = 78, sd = 15)
  x[D == 0] <- rnorm(sum(D == 0), mean = 71, sd = 15)
  
  return(x)
}

B <- bloodpressure_sample(sample_size, D_b)

# Simulate recovery
recovery_sample_b <- function(n, D, B) {
  
  # Set parameters (determined by trying to replicate the probabilities in the
  # original Simpson's paradox data)
  alpha <- 4.2
  beta_D <- 0.34
  beta_B <- -0.038
  
  p <- 1 / (1 + exp(-(alpha + beta_D * D + beta_B * B)))
  return(rbinom(n, 1, p))
}

R_b <- recovery_sample_b(sample_size, D_b, B)

# Return tibble with formatted columns
data_bloodpressure <- tibble(
  bloodP = B,
  drug = ifelse(D_b == 1, "Take", "Refuse"),
  recover = R_b == 1
)

bloodpressure_summary <- data_bloodpressure %>%
  mutate(
    bp_group = ifelse(bloodP < 75.0, "low", "high")
  ) %>%
  group_by(drug, bp_group) %>%
  summarise(
    recovery_rate = mean(recover),
    .groups = "drop"
  )


##### ESTIMATION

### Experiment 1: Gender as confound

# estimate the distribution of gender 
#   (by intercept-only regression with a mixture model)
fit_gender <- brm(
  formula = gender_std ~ 1,
  data = data_gender,
  family = mixture(gaussian, gaussian),
  iter = n_iter*2,
  prior = c(
    prior(normal(0, 1), Intercept, dpar = mu1),
    prior(normal(4, 1), Intercept, dpar = mu2))
)

# estimate the distribution of recovery rates as predicted
#   by gender and treatment
fit_recovery_gender <- brm(
  formula = recover ~ gender_std * drug,
  data = data_gender,
  family = bernoulli(link = "logit"),
  iter = n_iter*2
)

# draw from estimated testosterone distribution
posterior_gender <- tidybayes::predicted_draws(
  object = fit_gender,
  newdata = tibble(Intercept = 1),
  value = "testosterone",
  ndraws = 100
) |> ungroup() |> 
  select(.draw, testosterone)

# assign each of the new draws the Take and the Refuse condition
newdata <- tibble(gender_std = rep(posterior_gender$testosterone, 2),
                  drug = append(rep("Take", length(posterior_gender$testosterone)),
                                rep("Refuse", length(posterior_gender$testosterone))))


# generate posterior recovery chances for the created sample
posterior_g <- tidybayes::epred_draws(
  object  = fit_recovery_gender,
  newdata = newdata,
  value   = "recovery",
  ndraws  = 100
) |> ungroup() |> 
  select(.draw, gender_std, drug, recovery)|>
  pivot_wider(names_from = drug, values_from = recovery)|>
  mutate(causal_effect = Take - Refuse)

# samples of posterior predictive samples when the drug is 
#   refused and taken
posterior_DrugRefused_g <- posterior_g$Refuse
posterior_DrugTaken_g <- posterior_g$Take

# summarizing posterior estimates
post_sum_g <- summarize_posterior(posterior_DrugRefused_g, 
                                  posterior_DrugTaken_g)




#plot
# plot posteriors
ggplot() + 
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             color = "darkgrey") +
  # samples of CE for gender
  geom_density(data=posterior_g, 
               mapping = aes(x=causal_effect,
                             color = "darkorange"), 
               kernel="gaussian",
               linewidth = 1,
               key_glyph = "path") +
  # 95% CIs 
  geom_segment(aes(x = post_sum_g$CI_lower[3],
                   xend = post_sum_g$CI_upper[3],
                   y = 1.6, yend = 1.6),
               linewidth = 1) +
  geom_point(aes(x = post_sum_g$mean[3], y=1.6)) +
  theme_classic() +
  theme(legend.position.inside = c(0.8, 0.8)) +
  labs(x = "causal effect",
       y = "density") +
  guides(color=guide_legend(title=NULL,
                            position = "inside")) +
  scale_color_identity(guide = "legend",
                       labels = c("Gender as \n confound")) 

