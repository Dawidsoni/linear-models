# Linear models - List 3
library(ggplot2)
library(MASS)
set.seed(100)

# Task 1

qt(1 - 0.05 / 2, df = 10) ^ 2
qf(1 - 0.05, df1 = 1, df2 = 10)


# Task 2

perform_f_test <- function (ssm, sse, dfE, significance_level = 0.05) {
  test_df <- data.frame(matrix(nrow = 1, ncol = 7))
  colnames(test_df) <- c("samples_count", "error_std", "f_df2", "f_value", "f_qunatile", "explained_var", "xy_corr")
  test_df$samples_count <- dfE + 2
  test_df$error_std <- sqrt(sse / dfE)
  test_df$f_df2 <- dfE
  test_df$f_value <- ssm / (sse / dfE)
  test_df$f_qunatile <- qf(1 - significance_level, df1 = 1, df2 = dfE)
  test_df$explained_var <- ssm / (ssm + sse)
  test_df$xy_corr <- sqrt(test_df$explained_var)
  return(test_df)
}

round(perform_f_test(ssm = 100, sse = 400, dfE = 20), digits = 3)


# Task 3

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

create_prediction_intervals <- function (x_values, y_values, linear_model, x_inputs, significance_level = 0.05) {
  predictions_df <- data.frame(matrix(nrow = length(x_inputs), ncol = 3))
  colnames(predictions_df) <- c("x_input", "y_min", "y_max")
  predictions_df$x_input <- x_inputs
  y_outputs <- linear_model$beta0 + linear_model$beta1 * x_inputs
  var_of_errors <- linear_model$var_of_errors
  ssx <- sum((x_values - mean(x_values)) ^ 2)
  standard_errors <- sqrt(var_of_errors * (1 + 1 / length(x_values) + ((x_inputs - mean(x_values)) ^ 2) / ssx))
  confidence_interval_length <- qt(1 - significance_level / 2, length(x_values) - 2) * standard_errors
  predictions_df$y_min = y_outputs - confidence_interval_length
  predictions_df$y_max = y_outputs + confidence_interval_length
  return(predictions_df)
}

plot_data_with_confidence_intervals <- function (
  data_df, linear_model, title, x_title, y_title, x_transform = identity, y_transform = identity,
  curve_x_transform = identity, curve_margin = 1.0
) {
  transformed_data_df <- data_df
  transformed_data_df$x_input <- x_transform(transformed_data_df$x_input)
  transformed_data_df$y_min <- y_transform(transformed_data_df$y_min)
  transformed_data_df$y_max <- y_transform(transformed_data_df$y_max)
  transformed_data_df$y_value <- y_transform(transformed_data_df$y_value)
  model_curve_func <- function (x) y_transform(linear_model$beta0 + linear_model$beta1 * curve_x_transform(x))
  ggplot(transformed_data_df) +
    geom_point(aes(x = x_input, y = y_value), shape = 19, size = 2.2, color = "blue") +
    geom_point(aes(x = x_input, y = y_min), shape = 19, size = 2.2, color = "darkorange") +
    geom_point(aes(x = x_input, y = y_max), shape = 19, size = 2.2, color = "darkorange") +
    geom_segment(aes(x = x_input, y = y_min, xend = x_input, yend = y_max), color = "black") +
    stat_function(fun = model_curve_func, color = "springgreen4", size = 1.2, n = 1000) +
    xlim(min(transformed_data_df$x_input) - curve_margin, max(transformed_data_df$x_input) + curve_margin) +
    xlab(x_title) +
    ylab(y_title) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
}


students_df <- read.table('table1_6.csv', col.names = c("id", "gpa", "iq", "gender", "piers_test"))
iq_linear_model <- get_linear_model(students_df$iq, students_df$gpa)
round(iq_linear_model, digits = 4)
create_prediction_intervals(
  students_df$iq, students_df$gpa, iq_linear_model, x_inputs = c(100), significance_level = 0.1
)
iq_joint_intervals <- create_prediction_intervals(
  students_df$iq, students_df$gpa, iq_linear_model, students_df$iq, significance_level = 0.05
)
iq_joint_intervals$y_value <- students_df$gpa
plot_data_with_confidence_intervals(
  iq_joint_intervals, iq_linear_model, title = "Prediction intervals of GPA ~ IQ linear model", x_title = "IQ",
  y_title = "GPA"
)


# Task 4

piers_test_linear_model <- get_linear_model(students_df$piers_test, students_df$gpa)
round(piers_test_linear_model, digits = 4)
create_prediction_intervals(
  students_df$piers_test, students_df$gpa, piers_test_linear_model, x_inputs = c(60), significance_level = 0.1
)
piers_test_joint_intervals <- create_prediction_intervals(
  students_df$piers_test, students_df$gpa, piers_test_linear_model, students_df$piers_test, significance_level = 0.05
)
piers_test_joint_intervals$y_value <- students_df$gpa
plot_data_with_confidence_intervals(
  piers_test_joint_intervals, piers_test_linear_model,
  title = "Prediction intervals of GPA ~ Piers-Harris score linear model", x_title = "Piers-Harris score",
  y_title = "GPA"
)


# Task 5

analyse_residuals <- function (x_values, y_values, linear_model, binwidth, min_x_lim, max_x_lim) {
  residuals <- y_values - (linear_model$beta0 + linear_model$beta1 * x_values)
  print(paste0("Sum of residuals: ", round(sum(residuals), digits = 5)))
  residual_df <- data.frame(matrix(nrow = length(residuals), ncol = 3))
  colnames(residual_df) <- c("x_value", "residual", "position")
  residual_df$x_value <- x_values
  residual_df$residual <- residuals
  residual_df$position <- seq_along(residuals)
  print(ggplot(residual_df) +
    geom_point(aes(x = x_value, y = residual), shape = 19, size = 2.2, color = "darkorange") +
    xlab("Machines count") +
    ylab("Residual") +
    ggtitle("Residuals depeneding on machines counts") +
    theme(plot.title = element_text(hjust = 0.5)))
  print(ggplot(residual_df) +
    geom_point(aes(x = position, y = residual), shape = 19, size = 2.2, color = "darkorange") +
    xlab("Measurement order") +
    ylab("Residual") +
    ggtitle("Residuals depending on measurement orders") +
    theme(plot.title = element_text(hjust = 0.5)))
  print(ggplot(residual_df, aes(x = residual)) +
    geom_histogram(aes(y = ..density..), fill = "blue", color = "black", binwidth = binwidth) +
    stat_function(fun = dnorm, args = list(mean = mean(residuals), sd = sqrt(var(residuals))), size = 2.0, n = 1000) +
    xlim(min_x_lim, max_x_lim) +
    ggtitle("Histogram of residuals") +
    xlab("Residuals") +
    theme(plot.title = element_text(hjust = 0.5)))
}

copiers_df <- read.table('ch01pr20.csv', col.names = c("service_time", "machines_count"))
copiers_linear_model <- get_linear_model(copiers_df$machines_count, copiers_df$service_time)
round(copiers_linear_model, digits = 4)
analyse_residuals(
  copiers_df$machines_count, copiers_df$service_time, copiers_linear_model, binwidth = 4.0, min_x_lim = -30.0,
  max_x_lim = 30.0
)


# Task 6

noisy_copiers_df <- copiers_df
noisy_copiers_df$service_time[1] <- 2000
noisy_copiers_linear_model <- get_linear_model(noisy_copiers_df$machines_count, noisy_copiers_df$service_time)
round(noisy_copiers_linear_model, digits = 4)
analyse_residuals(
  noisy_copiers_df$machines_count, noisy_copiers_df$service_time, noisy_copiers_linear_model, binwidth = 50.0,
  min_x_lim = -2000.0, max_x_lim = 2000.0
)


# Task 7

solutions_df <- read.table('ch03pr15.csv', col.names = c("concentration", "time"))
solutions_linear_model <- get_linear_model(solutions_df$time, solutions_df$concentration)
round(solutions_linear_model, digits = 4)

# Task 8

solutions_joint_intervals <- create_prediction_intervals(
  solutions_df$time, solutions_df$concentration, solutions_linear_model, solutions_df$time, significance_level = 0.05
)
solutions_joint_intervals$y_value <- solutions_df$concentration
plot_data_with_confidence_intervals(
  solutions_joint_intervals, solutions_linear_model,
  title = "Prediction intervals of solution time ~ solution concetration linear model", x_title = "Solution time",
  y_title = "Solution concentration"
)


# Task 9

boxcox(solutions_df$concentration ~ solutions_df$time)


# Task 10

log_solutions_df <- solutions_df
log_solutions_df$concentration <- log(log_solutions_df$concentration)
log_solutions_linear_model <- get_linear_model(log_solutions_df$time, log_solutions_df$concentration)
round(log_solutions_linear_model, digits = 4)
log_solutions_joint_intervals <- create_prediction_intervals(
  log_solutions_df$time, log_solutions_df$concentration, log_solutions_linear_model, log_solutions_df$time,
  significance_level = 0.05
)
log_solutions_joint_intervals$y_value <- log_solutions_df$concentration
plot_data_with_confidence_intervals(
  log_solutions_joint_intervals, log_solutions_linear_model,
  title = "Prediction intervals of solution time ~ log(solution concetration) linear model", x_title = "Solution time",
  y_title = "Log(solution concentration)"
)

# Task 11

plot_data_with_confidence_intervals(
  log_solutions_joint_intervals, log_solutions_linear_model, y_transform = exp,
  title = "Prediction intervals of solution time ~ log(solution concetration) linear model", x_title = "Solution time",
  y_title = "Solution concentration", curve_margin = 0.3
)


# Task 12

sqrt_solutions_df <- solutions_df
sqrt_solutions_df$time <- sqrt_solutions_df$time ^ (-0.5)
sqrt_solutions_linear_model <- get_linear_model(sqrt_solutions_df$time, sqrt_solutions_df$concentration)
round(sqrt_solutions_linear_model, digits = 4)
sqrt_solutions_joint_intervals <- create_prediction_intervals(
  sqrt_solutions_df$time, sqrt_solutions_df$concentration, sqrt_solutions_linear_model, sqrt_solutions_df$time,
  significance_level = 0.05
)
sqrt_solutions_joint_intervals$y_value <- sqrt_solutions_df$concentration
plot_data_with_confidence_intervals(
  sqrt_solutions_joint_intervals, sqrt_solutions_linear_model,
  title = "Prediction intervals of (solution time)^(-1/2) ~ solution concetration linear model",
  x_title = "(solution time)^(-1/2)", y_title = "Solution concentration", curve_margin = 0.05
)
sqrt_beta0 <- sqrt_solutions_linear_model$beta0
sqrt_beta1 <- sqrt_solutions_linear_model$beta1
plot_data_with_confidence_intervals(
  sqrt_solutions_joint_intervals, sqrt_solutions_linear_model, x_transform = function (x) x^(-2),
  curve_x_transform = function (x) x^(-1/2), curve_margin = 0.1,
  title = "Prediction intervals of (solution time)^(-1/2) ~ solution concetration linear model",
  x_title = "Solution time", y_title = "Solution concentration"
)
