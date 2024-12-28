# packages for Bayesian regression modelling
library(brms)      
# install faintr from Github (easy evaluation of brm models)
#devtools::install_github("michael-franke/faintr")

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

#### helper functions

#' function to compare posteriors of two conditions 
#' returns a dataframe with the means and CIs of both conditions, as well as
#' their difference
summarize_posterior <- function(post_do0, post_do1, 
                                label_do0 = "refuse drug",
                                label_do1 = "take drug") {
  
  # Ensure the inputs are tibbles and the columns are named "value"  
  # (this is ugly but necessary to ensure the 
  # function can handle tibbles and normal lists as input)
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
      CI_lower = quantile(value, 0.025),
      mean = mean(value),
      CI_upper = quantile(value, 0.975)
    )
  
  # reorder columns
  order <- c(label_do1, label_do0, "causal effect")
  post_sum <- post_sum |>
    mutate(condition = factor(condition, levels = order)) |>
    arrange(condition)
  
  return(post_sum)
  
}



# generate dataset with Simpsons paradox
data_simpsons_paradox <- tibble(
  gender = c("Male", "Male", "Female", "Female"),
  bloodP = c("Low", "Low", "High", "High"),
  drug = c("Take", "Refuse", "Take", "Refuse"),
  k = c(81, 234, 192, 55),
  N = c(87, 270, 263, 80),
  proportion = k/N
)

# cast into long format
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


#### Causal inference with gender as confound
n_iter = 2000

# fist, estimate the distribution of gender (by intercept-only regression)
# we need this later when we marginalize over gender
fit_gender <- brm(
  formula = gender ~ 1,
  data = data_SP_long,
  family = bernoulli(link = "logit"),
  iter = n_iter
)

# then estimate the distribution of recovery rates as predicted
# by gender and treatment (this is before the do-calculus step)
fit_recovery_gender <- brm(
  formula = recover ~ gender * drug,
  data = data_SP_long,
  family = bernoulli(link = "logit"),
  iter = n_iter
)



# sample from estimated gender dist
# why? because for the do(treatment) operation, we make gender 
# independent from treatment (remove causal connection)
# in the formula: we no longer compute P(G=g|D=d), but only P(G=g)

# Option A): estimate gender distribution by sampling data points
#            (build a virtual participant pool)

posterior_gender_sample <- tidybayes::predicted_draws(
  # predicted_draws draws from posterior predictive dist
  # P(y_new | x_new, y_obs), because we want data points
  object = fit_gender,
  newdata = tibble(Intercept = 1),
  value = "gender",
  ndraws = n_iter*2
) |> 
  ungroup() |>
  mutate(gender = ifelse(gender, "Male", "Female")) |>
  select(gender)


# posterior predictive samples for D=1 (do(D=1))
posterior_DrugTaken_g <- tidybayes::epred_draws(
  # epred_draws draws from the expectation of the post pred dist 
  # E(y_new | x_new, y_obs), because we want the mean
  object = fit_recovery_gender,
  newdata = posterior_gender_sample |> mutate(drug="Take"),
  value = "taken", #can we just leave this argument out and call them all "value"?
  ndraws = n_iter * 2
) |> ungroup() |>
  select(taken)

# posterior predictive samples for D=0
posterior_DrugRefused_g <- tidybayes::epred_draws(
  object = fit_recovery_gender,
  newdata = posterior_gender_sample |> mutate(drug="Refuse"),
  value = "refused",
  ndraws = n_iter * 2
) |> ungroup() |>
  select(refused)

# summarize results 
post_sum_g <- summarize_posterior(posterior_DrugRefused_g, posterior_DrugTaken_g)

# Option B): estimate gender distribution by sampling from the expected value
#            of the posterior predictive distribution
#            (estimate the fraction of men in the sample)
posterior_gender_proportion <- tidybayes::epred_draws(
  object = fit_gender,
  newdata = tibble(Intercept = 1),
  value = "maleProp", #proportion of men in the sample
  ndraws = n_iter * 2
) |> ungroup() |> 
  select(.draw, maleProp)

# posterior predictive samples for all combinations of drug and gender
newdata <- tibble(
  gender = factor(c("Male", "Male", "Female", "Female"), levels = levels(data_SP_long$gender)),
  drug = factor(c("Take", "Refuse", "Take", "Refuse"), levels = levels(data_SP_long$drug))
)

posterior_drug_g_smooth <- tidybayes::epred_draws(
  object  = fit_recovery_gender,
  newdata = newdata,
  value   = "recovery",
  ndraws  = n_iter * 2
) |> ungroup() |> 
  select(.draw, gender, drug, recovery) |>
  full_join(posterior_gender_proportion)  |> 
  mutate(weights = ifelse(gender == "Male", maleProp, 1-maleProp)) |> 
  group_by(`.draw`, drug) |> 
  summarize(predRecover = sum(recovery * weights)) |> 
  pivot_wider(names_from = drug, values_from = predRecover) |> 
  mutate(causal_effect = Take - Refuse) 

post_sum_g_smooth <- summarize_posterior(posterior_drug_g_smooth$Refuse, 
                                         posterior_drug_g_smooth$Take)





#### Causal inference with blood pressure as mediator

#' in this formula, we ignore bp alltogether, because we only care about the 
#' total causal effect that the drug has on recovery, independent of whether 
#' this effect ocurrs directly or via the change in bp
fit_recovery_bp <- brm(
  formula = recover ~ drug,
  data = data_SP_long,
  family = bernoulli(link = "logit"),
  iter = n_iter
)

posterior_DrugTaken_bp <-
  faintr::filter_cell_draws(fit_recovery_bp, drug == "Take") |>
  pull(draws) |>
  psych::logistic() |>
  # transform to tibble
  (\(.) tibble(taken=.))()


posterior_DrugRefused_bp <-
  faintr::filter_cell_draws(fit_recovery_bp, drug == "Refuse") |>
  pull(draws) |>
  psych::logistic() |>
  # transform to tibble
  (\(.) tibble(refused=.))()


# summarize results 
post_sum_bp <- summarize_posterior(posterior_DrugRefused_bp, posterior_DrugTaken_bp)


# compute posteriors of causal effect by subtracting DrugRefused from DrugTaken
posterior_effect_g <- posterior_DrugTaken_g - posterior_DrugRefused_g
posterior_effect_bp <- posterior_DrugTaken_bp - posterior_DrugRefused_bp

# plot posteriors
p <- ggplot() + 
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             color = "darkgrey") +
  # mean lines
  #geom_vline(aes(xintercept = post_sum_g$mean[3]),
  #           color = "darkgrey") +
  #geom_vline(aes(xintercept = post_sum_bp$mean[3]),
  #           color = "darkgrey") +
  # gender
  #geom_density(data=posterior_effect_g, 
  #             mapping = aes(x=taken,
  #                           color = "darkred"), 
  #             kernel="gaussian",
  #             linewidth = 1,
  #             key_glyph = "path") +
  geom_density(data=posterior_drug_g_smooth, 
               mapping = aes(x=causal_effect,
                             color = "darkred"), 
               kernel="gaussian",
               linewidth = 1,
               key_glyph = "path") +
  # blood pressure
  geom_density(data=posterior_effect_bp, 
               mapping = aes(x=taken,
                             color ="darkorange"), 
               kernel="gaussian",
               linewidth = 1,
               key_glyph = "path") +
  # CIs
  # ragged gender
  #geom_segment(aes(x = post_sum_g$CI_lower[3],
  #                 xend = post_sum_g$CI_upper[3],
  #                 y = 1.5, yend = 1.5),
  #             linewidth = 1) +
  # smooth gender
  geom_segment(aes(x = post_sum_g_smooth$CI_lower[3],
                   xend = post_sum_g_smooth$CI_upper[3],
                   y = 1.5, yend = 1.5),
               linewidth = 1) +
  geom_segment(aes(x = post_sum_bp$CI_lower[3],
                   xend = post_sum_bp$CI_upper[3],
                   y = 3, yend = 3),
               linewidth = 1) +
  # mean points
  #geom_point(aes(x = post_sum_g$mean[3], y=1.5)) +
  geom_point(aes(x = post_sum_g_smooth$mean[3], y=1.5)) +
  geom_point(aes(x = post_sum_bp$mean[3], y=3)) +
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



