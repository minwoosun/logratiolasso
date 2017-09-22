test_that("myfs runs and produces dimensionally correct output", {
  set.seed(1)
  n <- 100
  p <- 30

  x <- abs(matrix(rnorm(n*p),nrow = n))
  y <- x[, 1] + x[, 2] + .1 * rnorm(n)
  y <- y - mean(y)

  model <- myfs(x,y, nsteps = 5)

  expect_equal(nrow(model$ind), 5)
  expect_equal(length(model$beta[[5]]), 5)
  expect_equal(length(model$rss), 5)
})
