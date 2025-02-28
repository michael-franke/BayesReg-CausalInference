#### packages, options, seed ----

# for Bayesian regression modelling
library(brms)       # Bayesian regression models
library(HDInterval) # for calculating highest density intervals
library(tidybayes)  # for tidy output of Bayesian models
# for data wrangling 
library(tibble)     # dataframes
library(tidyr)      # for tidying dataframes
library(dplyr)      # intuitive data manipulation
# for plotting
library(ggplot2)    # for plots

# option for Bayesian regression models: 
# use all available cores for parallel computing
options(mc.cores = parallel::detectCores())

# seeding for reproducibility
set.seed(123)


#### global parameters ----

# number of iterations / samples
n_iter = 5000


#### helper functions ----

#' function to compare posteriors of two conditions 
#' returns a dataframe with the means and CIs of both conditions, as well as
#'   their difference
summarize_posterior <- function(post_do0, post_do1, 
                                label_do0 = "refuse drug",
                                label_do1 = "take drug") {
  
  # Ensure the inputs are tibbles and the columns are named "value"  
  # (this is ugly but necessary to ensure that the 
  #   function can handle tibbles and normal lists as input)
  post_do0 <- tibble(value = pull(as_tibble(post_do0)))
  post_do1 <- tibble(value = pull(as_tibble(post_do1)))
  
  #calculate difference between conditions
  diff <- post_do1 - post_do0
  
  # combine data into single tibble
  combined_data <- bind_rows(
    post_do1 |> mutate(condition = label_do1),
    post_do0 |> mutate(condition = label_do0),
    diff |> mutate(condition = "causal effect")
  )
  
  # calculate means and confidence intervals
  post_sum <- combined_data |>
    group_by(condition) |>
    summarize(
      CI_lower = HDInterval::hdi(value)[1],
      mean = mean(value),
      CI_upper = HDInterval::hdi(value)[2]
    )
  
  # reorder columns
  order <- c(label_do1, label_do0, "causal effect")
  post_sum <- post_sum |>
    mutate(condition = factor(condition, levels = order)) |>
    arrange(condition)
  
  return(post_sum)
  
}

#### preparing the data set ----

# generate dataset with Simpsons paradox
data_simpsons_paradox <- tibble(
  gender = c("Male", "Male", "Female", "Female"),
  bloodP = c("Low", "Low", "High", "High"),
  drug = c("Take", "Refuse", "Take", "Refuse"),
  k = c(81, 234, 192, 55),
  N = c(87, 270, 263, 80),
  proportion = k/N
)

# cast into long format (one observation per row)
data_SP_long <- rbind(
  data_simpsons_paradox |> 
    uncount(k) |>
    mutate(recover = TRUE) |>
    select(-N, -proportion), 
  data_simpsons_paradox |>
    uncount(N-k) |>
    mutate(recover=FALSE) |>
    select(-k, -N, -proportion)
  ) |> # this part is new
  mutate(
    gender = factor(gender), 
    drug = factor(drug)
  )


#### scenario 1: causal inference with gender as confound ----

# estimate the distribution of gender 
#   (by intercept-only regression)
fit_gender <- brm(
  formula = gender ~ 1,
  data = data_SP_long,
  family = bernoulli(link = "logit"),
  iter = n_iter
)

# estimate the distribution of recovery rates as predicted
#   by gender and treatment
fit_recovery_gender <- brm(
  formula = recover ~ gender * drug,
  data = data_SP_long,
  family = bernoulli(link = "logit"),
  iter = n_iter
)

#  estimate gender distribution by sampling from the expected value
#    of the posterior predictive distribution
#  (estimate the fraction of men in the sample)
posterior_gender_proportion <- tidybayes::epred_draws(
  object = fit_gender,
  newdata = tibble(Intercept = 1),
  value = "maleProp", #proportion of men in the sample
  ndraws = n_iter * 2
) |> ungroup() |> 
  select(.draw, maleProp)

# make a data frame with all compbinations of values for variables
#   of drug and gender
newdata <- tibble(
  gender = factor(c("Male", "Male", "Female", "Female"), levels = levels(data_SP_long$gender)),
  drug = factor(c("Take", "Refuse", "Take", "Refuse"), levels = levels(data_SP_long$drug))
)

# get posterior predictive samples for each combination of values
#  for drug and genders
posterior_g <- tidybayes::epred_draws(
  object  = fit_recovery_gender,
  newdata = newdata,
  value   = "recovery",
  ndraws  = n_iter * 2
) |> ungroup() |> 
  select(.draw, gender, drug, recovery) |>
  full_join(posterior_gender_proportion)

# obtain estimates of causal effect by
#   calculating the weighted average of recovery rates for each 
#   value of drug, weighted by the est. proportions of males
posterior_g <- posterior_g |> 
  mutate(weights = ifelse(gender == "Male", 
                          maleProp, 1-maleProp)) |> 
  group_by(`.draw`, drug) |> 
  summarize(predRecover = sum(recovery * weights)) |> 
  pivot_wider(names_from = drug, values_from = predRecover) |> 
  # causal effect estimates are obtained from by subtracting 
  #   the predicted recovery rates when the drug is refused from
  #   the predicted recovery rates when the drug is taken
  mutate(causal_effect = Take - Refuse) 

# samples of posterior predictive samples when the drug is 
#   refused and taken
posterior_DrugRefused_g <- posterior_g$Refuse
posterior_DrugTaken_g <- posterior_g$Take

# summarizing posterior estimates
post_sum_g <- summarize_posterior(posterior_DrugRefused_g, 
                                  posterior_DrugTaken_g)
                                  

#### scenario 2: causal inference with blood pressure as mediator ----

# fitting a regression model to predict recovery based on drug use alone
fit_recovery_bp <- brm(
  formula = recover ~ drug,
  data = data_SP_long,
  family = bernoulli(link = "logit"),
  iter = n_iter
)

# posterior predictive samples for the recovery rate
#   when the drug is refused and taken
posterior_bp <- 
  tidybayes::epred_draws(
    object = fit_recovery_bp,
    newdata = tibble(drug = c("Refuse", "Take")),
    value = "predRecover", # prob. of recovery
    ndraws = n_iter * 2
  ) |> ungroup() |> 
  select(.draw, drug, predRecover) |> 
  pivot_wider(names_from = drug, values_from = predRecover) |> 
  # causal effect estimates are obtained from by subtracting 
  #   the predicted recovery rates when the drug is refused from
  #   the predicted recovery rates when the drug is taken
  mutate(causal_effect = Take - Refuse) 

# samples of posterior predictive samples when the drug is 
#   refused and taken
posterior_DrugRefused_bp <- posterior_bp$Refuse
posterior_DrugTaken_bp   <- posterior_bp$Take

# summarizing posterior estimates
post_sum_bp <- summarize_posterior(posterior_DrugRefused_bp, 
                                   posterior_DrugTaken_bp)

#### plotting ----

# plot posteriors
ggplot() + 
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             color = "darkgrey") +
  # samples of CE for gender
  geom_density(data=posterior_g, 
               mapping = aes(x=causal_effect,
                             color = "darkred"), 
               kernel="gaussian",
               linewidth = 1,
               key_glyph = "path") +
  # samples of CE for gender for blood pressure
  geom_density(data=posterior_bp, 
               mapping = aes(x=causal_effect,
                             color ="darkorange"), 
               kernel="gaussian",
               linewidth = 1,
               key_glyph = "path") +
  # 95% CIs 
  geom_segment(aes(x = post_sum_g$CI_lower[3],
                   xend = post_sum_g$CI_upper[3],
                   y = 1.6, yend = 1.6),
               linewidth = 1) +
  geom_segment(aes(x = post_sum_bp$CI_lower[3],
                   xend = post_sum_bp$CI_upper[3],
                   y = 2.9, yend = 2.9),
               linewidth = 1) +
  # mean points
  geom_point(aes(x = post_sum_g$mean[3], y=1.6)) +
  geom_point(aes(x = post_sum_bp$mean[3], y=2.9)) +
  # design
  theme_classic() +
  theme(legend.position.inside = c(0.8, 0.8)) +
  labs(x = "causal effect",
       y = "density") +
  guides(color=guide_legend(title=NULL,
                            position = "inside")) +
  scale_color_identity(guide = "legend",
                       labels = c("Blood pressure \n as mediator",
                                  "Gender as \n confound")) 

# save figure
ggsave(plot = last_plot(), filename = "posterior_causal_effect_smooth.png",
       width = 6, height = 4)



