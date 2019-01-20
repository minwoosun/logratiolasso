## ------------------------------------------------------------------------
library(logratiolasso)
devtools::load_all()

set.seed(10)
n <- 100 #number of observations
p <- 20 #number of features

x <- abs(matrix(rnorm(n*p), nrow = n)) #positive raw features
w <- log(x) #logarithmically transformed features
y <- w[, 1] - w[, 2] + rnorm(n) #response

## ------------------------------------------------------------------------
centered_w <- scale(w, center = TRUE, scale = FALSE)
centered_y <- y - mean(y)

model_fit <- glmnet.constr(centered_w, y, family = "gaussian")

## ------------------------------------------------------------------------
dim(model_fit$beta)

## ------------------------------------------------------------------------
cv_model_fit <- cv.glmnet.constr(model_fit, centered_w, y)

## ------------------------------------------------------------------------
cv_model_fit$cvm #CV estimate of error

## ------------------------------------------------------------------------
cv_model_fit$beta #best beta value

## ---- warning = FALSE----------------------------------------------------
ts_model <- two_stage(centered_w, centered_y, k_max = 5)

## ------------------------------------------------------------------------
ts_model$betas[[10]]

## ------------------------------------------------------------------------
cv_ts_model <- cv_two_stage(centered_w, centered_y, k_max = 5)

## ------------------------------------------------------------------------
cv_ts_model$lambda_min #index of best lambda
cv_ts_model$k_min #number of ratios

## ------------------------------------------------------------------------
cv_ts_model$beta_min

## ---- warning = FALSE----------------------------------------------------
cv_ts_model2 <- cv_two_stage(centered_w, centered_y, k_max = 5, second.stage = "yhat")
cv_ts_model2$beta_min

## ------------------------------------------------------------------------
afs_model <- approximate_fs(w, y, k_max = 5)
afs_model$beta

## ------------------------------------------------------------------------
afs_cv <- cv_approximate_fs(x, y, k_max = 5, n_folds = 10)
afs_cv$cvm

