# set working directory to file location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

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
source("utils.R")

# option for Bayesian regression models: 
# use all available cores for parallel computing
options(mc.cores = parallel::detectCores())

# seeding for reproducibility
set.seed(123)

# color definition
mygreen <- rgb(34,178,34,maxColorValue = 255)

#### global parameters ----

# number of iterations / samples
n_iter = 5000


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

# make a data frame with all combinations of values for variables
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
  select(.draw, gender, drug, recovery) 

# obtain estimates of causal effect by
#   calculating the weighted average of recovery rates for each 
#   value of drug, weighted by the est. proportions of males
posterior_g <- posterior_g |>
  # add the previous samples of P(G=1)
  #  and flip, where necessary, to P(G = 0)
  full_join(posterior_gender_proportion) |> 
  mutate(weights = ifelse(gender == "Male", 
                          maleProp, 1-maleProp)) |> 
  # calculate P(G=g) * P(R|G=g, D=d) for each combination of G and D
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

#Combine data of both experiments in one data format
posterior_combined <- bind_rows(
  tibble(causal_effect = posterior_g$causal_effect) %>% 
    mutate(mediator_label = "Gender as \nconfound"),
  tibble(causal_effect = posterior_bp$causal_effect) %>% 
    mutate(mediator_label = "Blood pressure \nas mediator")
)

# Combine summary data for means and CIs
summary_combined <- bind_rows(
  post_sum_g[3, ] %>% mutate(mediator_label = "Gender as \nconfound"),
  post_sum_bp[3, ] %>% mutate(mediator_label = "Blood pressure \nas mediator")
) %>%
  select(mediator_label, mean, CI_lower, CI_upper)

# Plot
ggplot(posterior_combined, aes(x = causal_effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  
  # Density plot with fill and color
  geom_density(aes(fill = mediator_label, color = mediator_label),
               alpha = 0.4, linewidth = 1, kernel = "gaussian",
               key_glyph = "path") +
  
  # Confidence intervals
  geom_segment(data = summary_combined,
               aes(x = CI_lower, xend = CI_upper, y = 0.1, yend = 0.1),
               inherit.aes = FALSE,
               linewidth = 1) +
  
  # Mean points
  geom_point(data = summary_combined,
             aes(x = mean, y = 0.1),
             inherit.aes = FALSE,
             size = 3) +
  
  # Facet by mediator
  facet_wrap(~ mediator_label, ncol = 1, scales = "free_y") +
  
  # Clean design
  theme_classic() +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position.inside = c(0.8, 0.8),
    legend.spacing.y = unit(5, "lines"),
    legend.margin = margin(t = 5, b = 5)
  ) +
  
  labs(x = "Causal effect") +
  scale_color_manual(values = c("Gender as \nconfound" = "darkorange",
                                "Blood pressure \nas mediator" = mygreen),
                     labels = c( #this is necessary to force spacing
                       "\nGender as \nconfound",
                       " \nBlood pressure \nas mediator"  
                     )) +
  scale_fill_manual(values = c("Gender as \nconfound" = "darkorange",
                               "Blood pressure \nas mediator" = mygreen)) +
  guides(color = guide_legend(title = NULL,
                              position = "inside",
                              override.aes = list(alpha = 1)),
         fill = "none")






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
ggsave(plot = last_plot(), filename = "posterior_causal_effect.png",
       width = 6, height = 4)



