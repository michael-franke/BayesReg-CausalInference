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