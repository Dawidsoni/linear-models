# Linear models - List 1

# Task 1

plot_vectors <- function(vectors, title) {
  plot(vectors, xlab='X', ylab='Y', pch=16, col="blue", main=title)
}

vectors <- matrix(rnorm(n=200), nrow=100)
plot_vectors(vectors, title="Samples from N(0, I)")


# Task 2

transform_vectors_from_standard_norm <- function(vectors, mean_matrix, covariance_matrix) {
  a_matrix <- t(chol(covariance_matrix))
  b_matrix <- matrix(rep(mean_matrix, dim(vectors)[1]), nrow=dim(vectors)[2])
  transformed_vectors <- t(a_matrix %*% t(vectors) + b_matrix)
  return(transformed_vectors)
}

mean_matrix <- matrix(c(4.0, 2.0), nrow=2)
covariance_matrix1 <- matrix(c(1.0, 0.9, 0.9, 1.0), nrow=2)
transformed_vectors1 <- transform_vectors_from_standard_norm(vectors, mean_matrix, covariance_matrix1)
plot_vectors(transformed_vectors1, title="Samples from N((4, 2), (1, 0.9, 0.9, 1))")
covariance_matrix2 <- matrix(c(1.0, -0.9, -0.9, 1.0), nrow=2)
transformed_vectors2 <- transform_vectors_from_standard_norm(vectors, mean_matrix, covariance_matrix2)
plot_vectors(transformed_vectors2, title="Samples from N((4, 2), (1, -0.9, -0.9, 1))")
covariance_matrix3 <- matrix(c(9.0, 0.0, 0.0, 1.0), nrow=2)
transformed_vectors3 <- transform_vectors_from_standard_norm(vectors, mean_matrix, covariance_matrix3)
plot_vectors(transformed_vectors3, title="Samples from N((4, 2), (9, 0, 0, 1))")


# Task 3

plot_samples_variances_covariances <- function(vectors, samples_count) {
  sample_covariance_matrix <- cov(vectors)
  variances <- diag(sample_covariance_matrix)
  hist(
    variances, main = paste("Plot of variances of", samples_count, "samples"),
    xlab = "Variances", breaks = 7, col = "darkorange"
  )
  print(paste0("Mean of variances: ", round(mean(variances), digits=2)))
  covariances <- sample_covariance_matrix[lower.tri(sample_covariance_matrix, diag=FALSE)]
  hist(
    covariances, main = paste("Plot of covariances of", samples_count, "samples"),
    xlab = "Covariances", breaks = 7, col = "darkorange"
  )
  print(paste0("Mean of covariances: ", round(mean(covariances), digits=2)))
}

analyse_vectors_transform <- function(dimension, samples_count) {
  vectors <- matrix(rnorm(n=dimension * samples_count), nrow=samples_count)
  mean_vector <- matrix(rep(0, dimension), nrow=dimension)
  covariance_matrix <- matrix(rep(0.9, dimension * dimension), nrow=dimension)
  diag(covariance_matrix) <- 1
  transformed_vectors <- transform_vectors_from_standard_norm(
    vectors, mean_vector, covariance_matrix
  )
  plot_samples_variances_covariances(transformed_vectors, samples_count)
}

analyse_vectors_transform(dimension = 100, samples_count = 200)
analyse_vectors_transform(dimension = 100, samples_count = 20000)
