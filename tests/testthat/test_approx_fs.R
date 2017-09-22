test_that("approximate_fs runs", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  k_max <- 5
  model <- approximate_fs(x, y, k_max)

  expect_equal(dim(model$beta), c(p, k_max))
  expect_equal(length(model$intercept), k_max)
})

test_that("predict.approximate_fs works in simple case", {
  n <- 5
  p <- 3
  k <- 2
  x <- matrix(2, n, p)

  fs_model <- list(beta = matrix(1, p, k), intercept = c(10, 12))
  predictions <- predict_approximate_fs(fs_model, x)

  expect_equal(dim(predictions), c(n, k))
  expect_equal(predictions[1, ], c(16, 18))
  expect_equal(predictions[2, ], c(16, 18))
})

test_that("cv.approximate_fs runs", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  k_max <- 5
  model <- cv_approximate_fs(x, y, k_max)

  expect_equal(length(model$cvm), k_max)
  expect_equal(length(model$beta), p)
})

test_that("cv.approximate_fs recovers noiseless case", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] - x[, 2]
  y <- y - mean(y)

  k_max <- 1
  model <- cv_approximate_fs(x, y, k_max)

  expect_equal(model$beta, c(1,-1, rep(0,p-2)))
})
