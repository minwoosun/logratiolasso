test_that("gaussian two-stage runs", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  lambda_1 <- exp(0:-10)

  k_max <- 5
  model <- two_stage(x, y, k_max = k_max, lambda_1 = lambda_1, second.stage = "y")

  expect_equal(model$lambda, lambda_1)
  expect_equal(dim(model$selected_vars), c(p, length(lambda_1)))
  expect_equal(nrow(model$coef[[length(lambda_1)]]$theta_ind), k_max)

  expect_equal(length(model$betas), length(lambda_1))
  expect_equal(dim(model$betas[[length(lambda_1)]]), c(p, k_max))
})

test_that("two-stage recovers noiseless case", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] - x[, 2]
  y <- y - mean(y)

  lambda_1 <- exp(0:-10)
  k_max <- 1
  model <- two_stage(x, y, k_max = k_max, lambda_1 = lambda_1, second.stage = "y")

  expect_equal(model$coef[[length(lambda_1)]]$theta_vals[[1]], c("X" = 1))
})

test_that("gaussian two-stage conservative runs", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  lambda_1 <- exp(0:-10)

  k_max <- 5
  model <- two_stage(x, y, k_max = k_max, lambda_1 = lambda_1, second.stage = "yhat")

  expect_equal(model$lambda, lambda_1)
  expect_equal(dim(model$selected_vars), c(p, length(lambda_1)))
  expect_equal(nrow(model$coef[[length(lambda_1)]]$theta_ind), k_max)
})

test_that("out_to_beta is dimensionally correct", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  lambda_1 <- exp(0:-10)

  k_max <- 5
  model <- two_stage(x, y, k_max = k_max, lambda_1 = lambda_1, second.stage = "y")
  output <- out_to_beta(model$coef[[length(lambda_1)]], k_max, p)

  expect_equal(dim(output), c(p, k_max))
  expect_lte(sum(output[, 1] != 0), 2)
  expect_lte(sum(output[, 2] != 0), 4)
  expect_lte(sum(output[, 3] != 0), 6)
  expect_lte(sum(output[, 4] != 0), 8)
  expect_lte(sum(output[, 5] != 0), 10)
})

test_that("predict_two_stage is dimensionally correct", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  lambda_1 <- exp(0:-10)

  k_max <- 5
  model <- predict_two_stage(x, y, x,
                             k_max = k_max, lambda_1 = lambda_1, second.stage = "y")

  expect_equal(dim(model$y_pred[[length(lambda_1)]]), c(n, k_max))
})

test_that("cv_two_stage is dimensionally correct", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  lambda_1 <- exp(0:-10)

  k_max <- 5
  model <- cv_two_stage(x, y, k_max = k_max, lambda_1 = lambda_1, second.stage = "y")

  expect_equal(dim(model$mse), c(k_max, length(lambda_1)))
  expect_equal(length(model$beta_min), p)
})

test_that("cv_two_stage conservative is dimensionally correct", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  lambda_1 <- exp(0:-10)

  k_max <- 5
  model <- cv_two_stage(x, y, k_max = k_max, lambda_1 = lambda_1, second.stage = "yhat")

  expect_equal(dim(model$mse), c(k_max, length(lambda_1)))
  expect_equal(length(model$beta_min), p)
})


