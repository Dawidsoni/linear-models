# Linear models - List 2
library(ggplot2)
options(warn = -1)
set.seed(1)

# Task 1

copiers_df <- read.table('CH01PR20', col.names = c("service_time", "machines_count"))
ggplot(copiers_df) +
  geom_point(aes(x = machines_count, y = service_time), shape = 19, size = 2.2, color="darkblue") +
  xlab('Number of machines serviced') +
  ylab('Service time') +
  ggtitle("The impact of the number of machines serviced on the total service time") +
  theme(plot.title = element_text(hjust = 0.5))


# Task 2

linear_model <- lm(service_time ~ machines_count, copiers_df)
summary(linear_model)
confint(linear_model)

# Task 3

# TODO


# Task 4

# TODO


# Task 5

# TODO


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


# Task 8

# TODO
