# Linear models - List 2
library(ggplot2)
options(warn = -1)
set.seed(100)

# Task 1

copiers_df <- read.table('CH01PR20', col.names = c("service_time", "machines_count"))
ggplot(copiers_df) +
  geom_point(aes(x = machines_count, y = service_time), shape = 19, size = 2.2, color="darkblue") +
  xlab('Number of machines serviced') +
  ylab('Service time') +
  ggtitle("The impact of the number of machines serviced on the total service time") +
  theme(plot.title = element_text(hjust = 0.5))


# Task 2

get_linear_model <- function (x_values, y_values) {
  linear_model <- data.frame(matrix(nrow = 1, ncol = 7))
  colnames(linear_model) <- c("beta0", "beta1", "var_of_errors", "t_statistic", "f_statistic", "p_value", "rse_squared")
  x_mean <- mean(x_values)
  y_mean <- mean(y_values)
  beta1 <- sum((x_values - x_mean) * (y_values - y_mean)) / sum((x_values - x_mean) ^ 2)
  beta0 <- y_mean - beta1 * x_mean
  predicted_values <- beta0 + beta1 * x_values
  errors <-  y_values - (beta0  + beta1 * x_values)
  var_of_errors <- var(errors) * (length(x_values) - 1) / (length(x_values) - 2)
  beta1_se <- sqrt(var_of_errors / sum((x_values - mean(x_values)) ^ 2))
  t_statistic <-  beta1 / beta1_se
  ssm <- sum((predicted_values - y_mean) ^ 2)
  sse <- sum((predicted_values - y_values) ^ 2)
  mse <- sse / (length(y_values) - 2)
  f_statistic <- ssm / mse
  p_value <- pf(f_statistic, df1 = 1, df2 = length(y_values) - 2, lower.tail = FALSE)
  rse_squared <- ssm / (ssm + sse)
  linear_model["beta0"] <- beta0
  linear_model["beta1"] <- beta1
  linear_model["var_of_errors"] <- var_of_errors
  linear_model["t_statistic"] <- t_statistic
  linear_model["f_statistic"] <- f_statistic
  linear_model["p_value"] <- p_value
  linear_model["rse_squared"] <- rse_squared
  return(linear_model)
}

round(get_linear_model(copiers_df$machines_count, copiers_df$service_time), digits = 4)

# Task 3

fit_and_estimate_value <- function (x_values, y_values, estimation_x, significance_level = 0.05) {
  linear_model <- get_linear_model(x_values, y_values)
  beta0 <- linear_model$beta0
  beta1 <- linear_model$beta1
  estimation_y <- beta0 + beta1 * estimation_x
  values_count <- length(x_values)
  x_mean <- mean(x_values)
  errors <- y_values - (beta0 + beta1 * x_values)
  ssx <- sum((x_values - x_mean) ^ 2)
  standard_error <- sqrt(var(errors) * (1 / values_count + ((estimation_x - x_mean) ^ 2) / ssx))
  confidence_interval_length <- qt(1 - significance_level / 2, values_count - 2) * standard_error
  min_interval <- estimation_y - confidence_interval_length
  max_interval <- estimation_y + confidence_interval_length
  return(data.frame(
    estimation_y=estimation_y, standard_error=standard_error, min_interval=min_interval, max_interval=max_interval
  ))
}

round(fit_and_estimate_value(copiers_df$machines_count, copiers_df$service_time, estimation_x = 11.0), digits = 4)


# Task 4

fit_and_predict_value <- function (x_values, y_values, estimation_x, significance_level = 0.05) {
  linear_model <- get_linear_model(x_values, y_values)
  beta0 <- linear_model$beta0
  beta1 <- linear_model$beta1
  estimation_y <- beta0 + beta1 * estimation_x
  values_count <- length(x_values)
  x_mean <- mean(x_values)
  errors <- y_values - (beta0 + beta1 * x_values)
  ssx <- sum((x_values - x_mean) ^ 2)
  standard_error <- sqrt(var(errors) * (1 + 1 / values_count + ((estimation_x - x_mean) ^ 2) / ssx))
  confidence_interval_length <- qt(1 - significance_level / 2, values_count - 2) * standard_error
  min_interval <- estimation_y - confidence_interval_length
  max_interval <- estimation_y + confidence_interval_length
  return(data.frame(
    estimation_y=estimation_y, standard_error=standard_error, min_interval=min_interval, max_interval=max_interval
  ))
}

round(fit_and_predict_value(copiers_df$machines_count, copiers_df$service_time, estimation_x = 11.0), digits = 4)


# Task 5

plot_prediction_interval <- function (x_values, y_values, significance_level = 0.05) {
  linear_model <- get_linear_model(x_values, y_values)
  beta0 <- linear_model$beta0
  beta1 <- linear_model$beta1
  predicted_values <- beta0 + beta1 * x_values
  values_count <- length(x_values)
  x_mean <- mean(x_values)
  errors <- y_values - (beta0 + beta1 * x_values)
  ssx <- sum((x_values - x_mean) ^ 2)
  standard_errors <- sqrt(var(errors) * (1 + 1 / values_count + ((x_values - x_mean) ^ 2) / ssx))
  fisher_constant <- sqrt(2 * qf(1 - significance_level, 2, values_count - 2))
  confidence_interval_lengths <- fisher_constant * standard_errors
  min_values <- predicted_values - confidence_interval_lengths
  max_values <- predicted_values + confidence_interval_lengths
  data_df <- data.frame(x_value = x_values, min_value = min_values, max_value = max_values)
  ggplot(data_df) +
    geom_point(aes(x = x_value, y = min_value), shape = 19, size = 2.2, color = "red") +
    geom_point(aes(x = x_value, y = max_value), shape = 19, size = 2.2, color = "red") +
    xlab('Number of machines serviced') +
    ylab('Service time') +
    ggtitle("Confidence intervals for individual observations") +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_prediction_interval(copiers_df$machines_count, copiers_df$service_time)


# Task 6

calculate_test_power <- function(elements_count, error_var, ssx, beta_value, significance_level = 0.05) {
  beta1_sd <- error_var / ssx
  t_mean <- beta_value / sqrt(beta1_sd)
  quantile <- qt(1 - significance_level / 2, df = elements_count - 2)
  lower_pbb <- pt(-quantile, df = elements_count - 2, ncp = t_mean)
  upper_pbb <- pt(quantile, df = elements_count - 2, ncp = t_mean)
  return (1 - upper_pbb + lower_pbb)
}

calculate_test_power(elements_count = 40, error_var = 120, ssx = 1000, beta_value = 1.0)
beta_values <- seq(from = -2, to = 2, by = 0.05)
test_powers <- c()
for (index in seq_along(beta_values)) {
  power <- calculate_test_power(elements_count = 40, error_var = 120, ssx = 1000, beta_value = beta_values[index])
  test_powers[[index]] <- power
}
powers_df <- data.frame(beta_value = beta_values, test_power = test_powers)
ggplot(powers_df) +
  geom_line(aes(x = beta_value, y = test_power), color="darkblue") +
  xlab('True slope value') +
  ylab('Test power') +
  ggtitle("The power for rejecting the null hypothesis that the regression slope is 0") +
  theme(plot.title = element_text(hjust = 0.5))

# Task 7

test_sample_rejection_frequency <- function (
  x_generator, error_generator, beta, experiments_count = 1000, samples_per_experiment = 200, significance_level = 0.05
) {
  x_samples <- x_generator(samples_count = samples_per_experiment)
  x_mean <- mean(x_samples)
  rejections <- 0
  for (i in 1:experiments_count) {
    error_samples <- error_generator(samples_count = samples_per_experiment)
    y_samples <- 5 + beta * x_samples + error_samples
    y_mean <- mean(y_samples)
    beta1_estimator <- sum((x_samples - x_mean) * (y_samples - y_mean)) / sum((x_samples - x_mean) ^ 2)
    beta0_estimator <- y_mean - beta1_estimator * x_mean
    errors <- y_samples - (beta0_estimator + beta1_estimator * x_samples)
    standard_error <- sqrt(var(errors) / sum((x_samples - x_mean) ^ 2))
    test_statistic <- beta1_estimator / standard_error
    is_rejected <- abs(test_statistic) >= qt(1 - significance_level / 2, samples_per_experiment - 2)
    if (is_rejected) {
      rejections <- rejections + 1
    }
  }
  return(rejections / experiments_count)
}

x_generator <- function (samples_count) rnorm(samples_count, mean = 0, sd = 5e-3)
error_generator1 <- function (samples_count) rnorm(samples_count, mean = 0, sd = 1.0)
error_generator2 <- function (samples_count) rexp(samples_count, rate = 1.0)

test_sample_rejection_frequency(x_generator, error_generator1, beta = 0)
test_sample_rejection_frequency(x_generator, error_generator2, beta = 0)
test_sample_rejection_frequency(x_generator, error_generator1, beta = 1.5)
test_sample_rejection_frequency(x_generator, error_generator2, beta = 1.5)
