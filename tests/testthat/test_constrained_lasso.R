test_that("gaussian constrained lasso runs", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  model <- glmnet.constr(x, y, family = "gaussian")
})

test_that("gaussian cv constrained lasso runs", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  model <- glmnet.constr(x, y, family = "gaussian")
  cv_model <- cv.glmnet.constr(model, x, y)
})

test_that("binomial constrained lasso runs", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- rbinom(n, 1, .5)

  model2 <- glmnet.constr(x, y, family = "binomial")
})

test_that("binomial cv constrained lasso runs", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- rbinom(n, 1, .5)

  model2 <- glmnet.constr(x, y, family = "binomial")
  cv_model2 <- cv.glmnet.constr(model2, x, y)
})

test_that("predict.glmnet.constr gives expected output", {
  fit <- list(a0 = 11, family = "gaussian",
              b = matrix( c(1, 2, 3,
                            1, 2, 3,
                            1, 2, 3),
                         nrow = 3, byrow= TRUE))
  x <- matrix(1, 3, 3)

  out <- predict.glmnet.constr(fit, x)

  expect_equal(out[1, ], c(14, 17, 20))
  expect_equal(out[1, ], out[2, ])
  expect_equal(out[1, ], out[3, ])
})


