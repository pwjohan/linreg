my_lm <- function(formula, data) {
  cl <- match.call()
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, mf)
  xtx_inv <- solve(t(X) %*% X)
  b_hat <- drop(xtx_inv %*% t(X) %*% y)
  fitted_values <- drop(X %*% b_hat)
  resid <- y - fitted_values
  df <- nrow(X) - ncol(X)
  values <- list(
    coefficients = b_hat,
    fitted.values = fitted_values,
    residuals = resid,
    df.residual = df,
    call = cl,
    xtx_inv = xtx_inv
  )
  class(values) <- "my_lm"
  return(values)
}

fit
