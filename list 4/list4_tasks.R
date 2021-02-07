# Linear models - List 4
library(ggplot2)
set.seed(100)

# Task 1

round(qt(p = 0.975, 17), digits = 4)


# Task 2

round(pf(1.5, df1 = 1, df2 = 20, lower.tail = FALSE), digits = 4)
round(pf(1.5, df1 = 2, df2 = 20, lower.tail = FALSE), digits = 4)
round(pf(6, df1 = 3, df2 = 20, lower.tail = FALSE), digits = 4)
round(pf(14.3478, df1 = 1, df2 = 22, lower.tail = FALSE), digits = 4)
round(sqrt(300 / 760), digits = 4)


# Task 3

transform_vectors_from_standard_norm <- function(vectors, mean_matrix, covariance_matrix) {
  a_matrix <- t(chol(covariance_matrix))
  b_matrix <- matrix(rep(mean_matrix, dim(vectors)[1]), nrow=dim(vectors)[2])
  transformed_vectors <- t(a_matrix %*% t(vectors) + b_matrix)
  return(transformed_vectors)
}


fit_and_explore_model_manually <- function (features_matrix, real_outputs, samples_count, df) {
  design_matrix <- cbind(rep(1, samples_count), features_matrix)
  estimated_betas <- solve(t(design_matrix) %*% design_matrix) %*% t(design_matrix) %*% real_outputs
  errors <- real_outputs - design_matrix %*% estimated_betas
  mse <- sum(errors ^ 2) / df
  beta1_var <- mse * solve(t(design_matrix) %*% design_matrix)[2, 2]
  t_statistic <- estimated_betas[2] / sqrt(beta1_var)
  estimated_p_value <- 2 * (1 - pt(abs(t_statistic), df))
  print(paste0("Estimated standard deviation: ", round(sqrt(beta1_var) ,digits = 4)))
  print(paste0("Estimated t-statistic: ", round(t_statistic, digits = 4)))
  print(paste0("Estimated p-value: ", round(estimated_p_value, digits = 4)))
}


explore_estimator_variance <- function (
  data_matrix, features_matrix, beta1, samples_count, df, simulations_count = 1000, significance_level = 0.05
) {
  design_matrix <- cbind(rep(1, samples_count), features_matrix)
  beta1_values <- rep(NA, simulations_count)
  null_hypothesis_rejection_count <- 0
  for (i in 1:simulations_count) {
    real_outputs <- beta1 * data_matrix[, 1]  + rnorm(samples_count)
    estimated_betas <- solve(t(design_matrix) %*% design_matrix) %*% t(design_matrix) %*% real_outputs
    errors <- real_outputs - design_matrix %*% estimated_betas
    mse <- sum(errors ^ 2) / df
    beta1_var <- mse * solve(t(design_matrix) %*% design_matrix)[2, 2]
    t_statistic <- estimated_betas[2] / sqrt(beta1_var)
    estimated_p_value <- 2 * (1 - pt(abs(t_statistic), df))
    if (estimated_p_value <= significance_level) {
      null_hypothesis_rejection_count <- null_hypothesis_rejection_count + 1
    }
    beta1_values[[i]] <- estimated_betas[2]
  }
  print(paste0("Standard deviation of beta1 estimators: ", round(sd(beta1_values), digits = 4)))
  print(paste0("Power of a test: ", round(null_hypothesis_rejection_count / simulations_count, digits = 4)))
}


explore_correlation_influence <- function (beta1) {
  features_matrix <- transform_vectors_from_standard_norm(
    vectors = matrix(rnorm(200), ncol = 2),
    mean_matrix = matrix(c(0.0, 0.0), nrow=2),
    covariance_matrix = matrix(c(1.0, 0.9, 0.9, 1.0) / 100.0, nrow=2)
  )
  feature1_vector <- features_matrix[, 1]
  feature2_vector <- features_matrix[, 2]
  real_outputs <- beta1 * feature1_vector  + rnorm(100)
  real_model <- lm(real_outputs ~ feature1_vector)
  extended_model <- lm(real_outputs ~ feature1_vector + feature2_vector)
  ci_of_real_model <- round(confint(real_model, 'feature1_vector', level = 0.95), digits = 4)
  ci_of_extended_model <- round(confint(extended_model, 'feature1_vector', level = 0.95), digits = 4)
  print(paste0("Real model confidence interval: [", ci_of_real_model[1], ", ", ci_of_real_model[2], "]"))
  print(paste0("Extended model confidence interval: [", ci_of_extended_model[1], ", ", ci_of_extended_model[2], "]"))
  pvalue_of_real_model <- round(summary(real_model)[["coefficients"]][2, 4], digits = 4)
  pvalue_of_extended_model <- round(summary(extended_model)[["coefficients"]][2, 4], digits = 4)
  print(paste0("P-value for real model: ", pvalue_of_real_model))
  print(paste0("P-value for extended model: ", pvalue_of_extended_model))
  print("Exploring properties of real model")
  fit_and_explore_model_manually(feature1_vector, real_outputs, samples_count = 100, df = 98)
  print("Exploring properties of extended model")
  fit_and_explore_model_manually(features_matrix, real_outputs, samples_count = 100, df = 98)
  print("Exploring variance of real model")
  explore_estimator_variance(features_matrix, feature1_vector, beta1, samples_count = 100, df = 98)
  print("Exploring variance of extended model")
  explore_estimator_variance(features_matrix, features_matrix, beta1, samples_count = 100, df = 98)
}


explore_correlation_influence(beta1 = 3)


# Task 4

calculate_aic_statistic <- function (model) {
  samples_count <- length(model$residuals)
  variables_count <- length(model$coefficients)
  sse <- sum(model$residuals ^ 2)
  return(samples_count * log(sse / samples_count) + 2 * variables_count)
}


construct_and_explore_model <- function (design_matrix, real_coefficients, real_outputs) {
  output_df <- data.frame(matrix(nrow = 1, ncol = 7))
  colnames(output_df) <- c("cols_count", "sse", "mse", "aic", "p_value1", "p_value2", "false_discoveries_count")
  samples_df <- data.frame(x = design_matrix, y = real_outputs)
  model <- lm(y ~ . -1, data = samples_df)
  output_df$cols_count <- dim(design_matrix)[2]
  output_df$sse <- sum(model$residuals ^ 2)
  output_df$mse <- sum((design_matrix %*% (model$coefficients - real_coefficients)) ^ 2)
  output_df$aic <- calculate_aic_statistic(model)
  p_values <- as.matrix(summary(model)$coefficients[, 4])
  output_df$p_value1 <- p_values[1]
  output_df$p_value2 <- p_values[2]
  is_null_hypothesis_true <- real_coefficients != 0
  output_df$false_discoveries_count <- sum((p_values < 0.05) != is_null_hypothesis_true)
  return(output_df)
}


explore_impact_of_multiple_dimensions <- function (
  samples_count, explanatory_variables_count, dimensions, reverse_columns = FALSE
) {
  output_df <- data.frame()
  design_matrix <- matrix(rnorm(samples_count * explanatory_variables_count, mean = 0.0, sd = 0.1), nrow = 1000)
  beta_variables <- c(rep(3.0, 5), rep(0.0, explanatory_variables_count - 5))
  random_errors <- rnorm(samples_count)
  real_outputs <- design_matrix %*% beta_variables + random_errors
  max_cols_count <- dim(design_matrix)[2]
  for (cols_count in dimensions) {
    columns_subset <- if (reverse_columns) max_cols_count:(max_cols_count - cols_count + 1) else 1:cols_count
    output_df <- rbind(output_df, construct_and_explore_model(
      as.matrix(design_matrix[, columns_subset]), beta_variables[columns_subset], real_outputs
    ))
  }
  return(output_df)
}


explored_dimensions <- c(1, 2, 5, 10, 50, 100, 500, 950)
explored_dimensions_df <- explore_impact_of_multiple_dimensions(
  samples_count = 1000, explanatory_variables_count = 950, dimensions = explored_dimensions
)
round(explored_dimensions_df, digits = 4)
explored_dimensions_df <- explore_impact_of_multiple_dimensions(
  samples_count = 1000, explanatory_variables_count = 950, dimensions = explored_dimensions, reverse_columns = TRUE
)
round(explored_dimensions_df, digits = 4)


# Task 5

patients_df <- read.csv("CH06PR15.txt", header = FALSE, sep = "")
colnames(patients_df) <- c("age", "severity_level", "anxiety_level", "satisfaction_level")
reduced_patients_model <- lm(satisfaction_level ~ 1, data = patients_df)
full_patients_model <- lm(satisfaction_level ~ age + severity_level + anxiety_level, data = patients_df)
summary(full_patients_model)
round(anova(reduced_patients_model, full_patients_model), digits = 4)


# Task 6

print(round(confint(full_patients_model, 'age', level = 0.95), digits = 3))
print(round(confint(full_patients_model, 'severity_level', level = 0.95), digits = 3))
print(round(confint(full_patients_model, 'anxiety_level', level = 0.95), digits = 3))


# Task 7

satisfactions_df <- data.frame(
  residual = full_patients_model$residuals,
  prediction = full_patients_model$fitted.values,
  age = patients_df$age,
  severity_level = patients_df$severity_level,
  anxiety_level = patients_df$anxiety_level
)
print(ggplot(satisfactions_df) +
    geom_point(aes(x = prediction, y = residual), shape = 19, size = 2.2, color = "dark orange") +
    xlab("Predicted value") +
    ylab("Residual") +
    ggtitle("Plot of predicted values versus residuals for the satisfaction model") +
    theme(plot.title = element_text(hjust = 0.5)))
print(ggplot(satisfactions_df) +
    geom_point(aes(x = age, y = residual), shape = 19, size = 2.2, color = "dark orange") +
    xlab("Age") +
    ylab("Residual") +
    ggtitle("Plot of peoples' ages versus residuals for the satisfaction model") +
    theme(plot.title = element_text(hjust = 0.5)))
print(ggplot(satisfactions_df) +
    geom_point(aes(x = severity_level, y = residual), shape = 19, size = 2.2, color = "dark orange") +
    xlab("Severity level") +
    ylab("Residual") +
    ggtitle("Plot of severity levels versus residuals for the satisfaction model") +
    theme(plot.title = element_text(hjust = 0.5)))
print(ggplot(satisfactions_df) +
    geom_point(aes(x = anxiety_level, y = residual), shape = 19, size = 2.2, color = "dark orange") +
    xlab("Anxiety level") +
    ylab("Residual") +
    ggtitle("Plot of anxiety levels versus residuals for the satisfaction model") +
    theme(plot.title = element_text(hjust = 0.5)))


# Task 8

print(ggplot(satisfactions_df, aes(sample = residual)) +
  ggtitle("QQ-plot of residuals versus standard normal distribution") +
  xlab("Theoretical value") +
  ylab("Residual") +
  geom_qq(colour="dark orange") +
  geom_qq_line(colour="blue") +
  theme(plot.title = element_text(hjust = 0.5)))
shapiro.test(satisfactions_df$residual)


# Task 9

cs_df <- read.csv("csdata.txt", header = FALSE, sep = "")
colnames(cs_df) <- c("id", "gpa", "hsm", "hss", "hse", "satm", "satv", "sex")
reduced_gpa_model <- lm(gpa ~ hsm + hss + hse, data = cs_df)
full_gpa_model <- lm(gpa ~ hsm + hss + hse + satm + satv, data = cs_df)
reduced_gpa_sse <- sum(reduced_gpa_model$residuals ^ 2)
full_gpa_sse <- sum(full_gpa_model$residuals ^ 2)
f_statistic <- ((reduced_gpa_sse - full_gpa_sse) / 2) / (full_gpa_sse / (dim(cs_df)[1] - 6))
round(f_statistic, digits = 4)
round(anova(reduced_gpa_model, full_gpa_model), digits = 4)


# Task 10


print_type2_sum_of_squares <- function (df, explanatory_variables, response_variable) {
  formula_of_full_model <- as.formula(
    paste(response_variable, paste(explanatory_variables, collapse=" + "), sep=" ~ ")
  )
  full_model <- lm(formula_of_full_model, data = df)
  for (variable_index in seq_along(explanatory_variables)) {
    chosen_explanatory_variables <- explanatory_variables[-variable_index]
    formula_of_reduced_model <- as.formula(
      paste(response_variable, paste(chosen_explanatory_variables, collapse=" + "), sep=" ~ ")
    )
    reduced_model <- lm(formula_of_reduced_model, data = df)
    print(round(anova(reduced_model, full_model), digits = 4))
  }
}


verify_hsm_type1_sum_of_squares <- function (df) {
  basic_model <- lm(gpa ~ satm + satv, data = df)
  extended_model <- lm(gpa ~ satm + satv + hsm, data = df)
  sse_of_basic_model <- sum(basic_model$residuals ^ 2)
  sse_of_extended_model <- sum(extended_model$residuals ^ 2)
  print(sse_of_basic_model - sse_of_extended_model)
}


sums_gpa_model <- lm(gpa ~ satm + satv + hsm + hse + hss, data = cs_df)
round(anova(sums_gpa_model), digits = 4)
print_type2_sum_of_squares(
  cs_df, explanatory_variables = c("satm", "satv", "hsm", "hse", "hss"), response_variable = "gpa"
)
verify_hsm_type1_sum_of_squares(cs_df)


# Task 11

cs_df$sat <- cs_df$satm + cs_df$satv
sat_gpa_model <- lm(gpa ~ satm + satv + sat, data = cs_df)
summary(sat_gpa_model)


# Task 12

plot_partial_regression <- function (df, explanatory_variables, examined_variable, response_variable) {
  response_variable_formula <- as.formula(
    paste(response_variable, paste(explanatory_variables, collapse=" + "), sep=" ~ ")
  )
  response_variable_model <- lm(response_variable_formula, data = df)
  examined_variable_formula <- as.formula(
    paste(examined_variable, paste(explanatory_variables, collapse=" + "), sep=" ~ ")
  )
  examined_variable_model <- lm(examined_variable_formula, data = df)
  residuals_df <- data.frame(
    response_variable = response_variable_model$residuals, examined_variable = examined_variable_model$residuals
  )
  print(ggplot(residuals_df, aes(x = response_variable, y = examined_variable)) +
    geom_point(shape = 19, size = 2.2, color = "dark orange") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    xlab("Residuals of a model predicting 'GPA'") +
    ylab(paste0("Residuals of a model predicting '", examined_variable, "'")) +
    ggtitle(paste0("Partial regression plot of the '", examined_variable, "' variable")) +
    theme(plot.title = element_text(hjust = 0.5)))
}


create_partial_regressions_plots <- function (df, explanatory_variables, response_variable) {
  for (variable_index in seq_along(explanatory_variables)) {
    chosen_explanatory_variables <- explanatory_variables[-variable_index]
    examined_variable <- explanatory_variables[variable_index]
    plot_partial_regression(df, chosen_explanatory_variables, examined_variable, response_variable)
  }
}


create_partial_regressions_plots(
  cs_df,
  explanatory_variables = c("hsm", "hss", "hse", "satm", "satv", "sex"),
  response_variable = "gpa"
)


# Task 13

plot_deleted_residuals <- function (df, explanatory_variables, response_variable) {
  df$intercept <- 1.0
  model_formula <- as.formula(
    paste(response_variable, paste(explanatory_variables, collapse=" + "), sep=" ~ ")
  )
  data_matrix <- data.matrix(df[c("intercept", c(explanatory_variables))])
  h_matrix <- data_matrix %*% solve(t(data_matrix) %*% data_matrix) %*% t(data_matrix)
  studentized_residuals <- rep(NA, dim(df)[1])
  for (index in seq(1, dim(df)[1])) {
    sample_discarded_df <- df[-index, ]
    sample_discarded_model <- lm(model_formula, data = sample_discarded_df)
    enumerator_value <- df[index, "gpa"] - predict(sample_discarded_model, df[index, ])
    denominator_value <- sqrt(sum((sample_discarded_model$residuals) ^ 2) * (1 - h_matrix[index, index]))
    studentized_residuals[index] <- enumerator_value / denominator_value
  }
  studentized_residuals_df <- data.frame(index = seq(1, length(studentized_residuals)), value = studentized_residuals)
  print(ggplot(studentized_residuals_df, aes(x = index, y = value)) +
    geom_point(shape = 19, size = 2.2, color = "dark orange") +
    xlab("Sample index") +
    ylab("Studentized deleted residual") +
    ggtitle("Plot of studentized deleted residuals of the GPA model") +
    theme(plot.title = element_text(hjust = 0.5)))
}

plot_deleted_residuals(
  cs_df,
  explanatory_variables = c("hsm", "hss", "hse", "satm", "satv", "sex"),
  response_variable = "gpa"
)


# Task 14

examine_dffits <- function (df, explanatory_variables, response_variable) {
  df$intercept <- 1.0
  model_formula <- as.formula(
    paste(response_variable, paste(explanatory_variables, collapse=" + "), sep=" ~ ")
  )
  full_model <- lm(model_formula, data = df)
  data_matrix <- data.matrix(df[c("intercept", c(explanatory_variables))])
  h_matrix <- data_matrix %*% solve(t(data_matrix) %*% data_matrix) %*% t(data_matrix)
  dffits_values <- rep(NA, dim(data_matrix)[1])
  for (index in seq(1, dim(data_matrix)[1])) {
    sample_discarded_df <- df[-index, ]
    sample_discarded_model <- lm(model_formula, data = sample_discarded_df)
    enumerator_value <- predict(full_model, df[index, ]) - predict(sample_discarded_model, df[index, ])
    sample_discarded_sse <- sum((sample_discarded_model$residuals) ^ 2)
    sample_discarded_mse <- sample_discarded_sse / (dim(data_matrix)[1] - dim(data_matrix)[2])
    denominator_value <- sqrt(sample_discarded_mse * h_matrix[index, index])
    dffits_values[index] <- enumerator_value / denominator_value
  }
  dffits_df <- data.frame(index = seq(1, length(dffits_values)), value = dffits_values)
  print(ggplot(dffits_df, aes(x = index, y = value)) +
    geom_point(shape = 19, size = 2.2, color = "dark orange") +
    xlab("Sample index") +
    ylab("DFFITS value") +
    ggtitle("Plot of DFFITS values of the GPA model") +
    theme(plot.title = element_text(hjust = 0.5)))
}


examine_dffits(
  cs_df,
  explanatory_variables = c("hsm", "hss", "hse", "satm", "satv", "sex"),
  response_variable = "gpa"
)


# Task 15

plot_tolerance_statistic <- function (df, explanatory_variables) {
  tolerances <- rep(NA, length(explanatory_variables))
  for (variable_index in seq_along(explanatory_variables)) {
    chosen_explanatory_variables <- explanatory_variables[-variable_index]
    response_variable <- explanatory_variables[variable_index]
    variable_formula <- as.formula(
      paste(response_variable, paste(chosen_explanatory_variables, collapse=" + "), sep=" ~ ")
    )
    variable_model <- lm(variable_formula, data = df)
    model_sse <- sum((variable_model$residuals) ^ 2)
    model_sst <- sum((df[[response_variable]] - mean(df[[response_variable]])) ^ 2)
    model_r_squared <- 1 - model_sse / model_sst
    tolerances[[variable_index]] <- 1 - model_r_squared
  }
  tolerances_df <- data.frame(variable = explanatory_variables, tolerance = tolerances)
  ggplot(tolerances_df, aes(x = variable, y = tolerance)) +
    geom_bar(stat = "identity", fill = "blue") +
    xlab("Variable name") +
    ylab("Tolerance") +
    ggtitle("Plot of the tolerance statistc for explanatory variables of the GPA model") +
    theme(plot.title = element_text(hjust = 0.5))
}



plot_tolerance_statistic(
  cs_df, explanatory_variables = c("hsm", "hss", "hse", "satm", "satv", "sex")
)


# Task 16


calculate_bic_statistic <- function (model) {
  samples_count <- length(model$residuals)
  variables_count <- length(model$coefficients)
  sse <- sum(model$residuals ^ 2)
  return(samples_count * log(sse / samples_count) + log(samples_count) * variables_count)
}


create_and_rate_submodels <- function (df, explanatory_variables, response_variable) {
  ratings_df <- data.frame(matrix(nrow = 2 ^ length(explanatory_variables) - 1, ncol = 3))
  colnames(ratings_df) <- c("variables", "aic", "bic")
  index_of_model <- 1
  for (variables_count in seq_along(explanatory_variables)) {
    variables_subsets <- combn(explanatory_variables, variables_count, simplify = FALSE)
    for (chosen_explanatory_variables in variables_subsets) {
      formula_of_model <- as.formula(
        paste(response_variable, paste(chosen_explanatory_variables, collapse=" + "), sep=" ~ ")
      )
      model <- lm(formula_of_model, data = df)
      ratings_df$variables[index_of_model] <- paste(chosen_explanatory_variables, collapse=", ")
      ratings_df$aic[index_of_model] <- round(calculate_aic_statistic(model), digits = 2)
      ratings_df$bic[index_of_model] <- round(calculate_bic_statistic(model), digits = 2)
      index_of_model <- index_of_model + 1
    }
  }
  return(ratings_df)
}

cs_ratings_df <- create_and_rate_submodels(
  cs_df, explanatory_variables = c("hsm", "hss", "hse", "satm", "satv", "sex"), response_variable = "gpa"
)
cs_ratings_df
