#' Fit a linear regression model
#'
#' Fits a linear model by ordinary least squares, computing
#' \eqn{\hat\beta = (X^TX)^{-1}X^Ty} directly instead of calling `lm()`.
#'
#' @param formula A model formula of the form `y ~ x1 + x2`.
#' @param data A data frame containing the variables in the formula.
#'
#' @returns An object of class `"my_lm"`: a list containing
#' \describe{
#'   \item{y}{The response vector.}
#'   \item{coefficients}{Named vector of estimated coefficients.}
#'   \item{fitted.values}{Fitted values \eqn{X\hat\beta}.}
#'   \item{residuals}{Residuals \eqn{y - X\hat\beta}.}
#'   \item{df.residual}{Residual degrees of freedom, \eqn{n - p}.}
#'   \item{call}{The matched call.}
#'   \item{xtx_inv}{The matrix \eqn{(X^TX)^{-1}}.}
#'   \item{terms}{The model terms, used by `predict()`.}
#' }
#' @export
#'
#' @examples
#' fit <- my_lm(mpg ~ wt + hp, data = mtcars)
#' fit
#' summary(fit)
#' predict(fit, newdata = data.frame(wt = 3, hp = 150))
my_lm <- function(formula, data) {
  cl <- match.call()
  if (is.data.frame(data) == FALSE) {
    stop("data must be a data frame")
  }
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, mf)
  if (qr(X)$rank < ncol(X)) {
    stop("Design matrix is rank deficient")
  }
  xtx_inv <- solve(t(X) %*% X)
  b_hat <- drop(xtx_inv %*% t(X) %*% y)
  fitted_values <- drop(X %*% b_hat)
  resid <- y - fitted_values
  df <- nrow(X) - ncol(X)
  values <- list(
    y = y,
    coefficients = b_hat,
    fitted.values = fitted_values,
    residuals = resid,
    df.residual = df,
    call = cl,
    xtx_inv = xtx_inv,
    terms = attr(mf, "terms")
  )
  class(values) <- "my_lm"
  return(values)
}

#' Summarize a my_lm fit
#'
#' Computes standard errors, t-tests, \eqn{R^2}, adjusted \eqn{R^2}, and the
#' overall F-test for a fitted `my_lm` model.
#'
#' @param object A fitted model of class `"my_lm"`.
#' @param ... Further arguments (currently ignored).
#'
#' @returns An object of class `"summary.my_lm"`: a list containing
#' \describe{
#'   \item{call}{The call used to fit the model.}
#'   \item{coefficients}{Matrix with columns for the estimate, standard
#'     error, t value, and two-sided p-value of each coefficient.}
#'   \item{residuals}{The model residuals.}
#'   \item{sigma}{Residual standard error, \eqn{\hat\sigma}.}
#'   \item{df}{Residual degrees of freedom, \eqn{n - p}.}
#'   \item{r.squared}{\eqn{R^2}.}
#'   \item{adj.r.squared}{Adjusted \eqn{R^2}.}
#'   \item{fstatistic}{Named vector: the F value and its numerator and
#'     denominator degrees of freedom.}
#' }
#' @export
#'
#' @examples
#' fit <- my_lm(mpg ~ wt + hp, data = mtcars)
#' summary(fit)
#' summary(fit)$r.squared
summary.my_lm <- function(object, ...) {
  n <- length(object$y)
  p <- length(object$coefficients)
  num_df <- p - 1
  den_df <- n - p
  sse <- sum(object$residuals^2)
  sst <- sum((object$y - mean(object$y))^2)
  sig2hat <- sse / den_df
  se <- sqrt(diag(sig2hat * object$xtx_inv))
  t_vals <- object$coefficients / se
  p_vals <- 2 * pt(abs(t_vals), den_df, lower.tail = FALSE)
  coef_table <- cbind(object$coefficients, se, t_vals, p_vals)
  colnames(coef_table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  r_squared <- 1 - (sse/sst)
  adj_r_squared <- 1 - ((1 - r_squared) * (n - 1)) / den_df
  f <- ((sst - sse) / (p - 1)) / (sse / den_df)
  f_stat <- c(f, num_df, den_df)
  names(f_stat) <- c("value", "numdf", "dendf")
  values <- list(
    call = object$call,
    coefficients = coef_table,
    residuals = object$residuals,
    sigma = sqrt(sig2hat),
    df = den_df,
    r.squared = r_squared,
    adj.r.squared = adj_r_squared,
    fstatistic = f_stat
  )
  class(values) <- "summary.my_lm"
  return(values)
}

#' Print a my_lm fit
#'
#' Displays the call and estimated coefficients of a fitted `my_lm` model.
#'
#' @param x An object of class `"my_lm"`.
#' @param ... Further arguments (currently ignored).
#'
#' @returns `x` invisibly
#' @export
print.my_lm <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Coefficients:\n")
  print(x$coefficients)
  invisible(x)
}

#' Print a summary of a my_lm fit
#'
#' Displays the call, residual quantiles, coefficient table, residual standard error, \eqn{R^2}, and F-statistic.
#'
#' @param x An object of class `"summary.my_lm"`
#' @param ... Further arguments (currently ignored)
#'
#' @returns `x` invisibly
#' @export
print.summary.my_lm <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Residuals:\n")
  res <- quantile(x$residuals)
  names(res) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(signif(res, digits = 4))
  cat("\n")
  cat("Coefficients:\n")
  printCoefmat(x$coefficients)
  cat("\n")
  cat("Residual standard error:", signif(x$sigma, 4), "on", x$df, "degrees of freedom\n")
  cat("Multiple R-squared:", signif(x$r.squared, 4), "Adjusted R-squared:", signif(x$adj.r.squared, 4), "\n")
  f_pval <- pf(x$fstatistic["value"], x$fstatistic["numdf"], x$fstatistic["dendf"], lower.tail = FALSE)
  cat("F-statistic:", signif(x$fstatistic["value"], 4), "on", x$fstatistic["numdf"], "and", x$fstatistic["dendf"], "DF, p-value:", signif(f_pval, 4), "\n")
  invisible(x)
}

#' Predict from a my_lm fit
#'
#' Computes predicted values from a fitted `my_lm` model, for the original
#' data or for new observations.
#'
#' @param object A fitted model of class `"my_lm"`.
#' @param newdata Optional data frame containing values of the predictor
#'   variables. If omitted, the fitted values are returned.
#' @param ... Further arguments (currently ignored).
#'
#' @returns A named numeric vector of predicted values.
#' @export
#'
#' @examples
#' fit <- my_lm(mpg ~ wt + hp, data = mtcars)
#' predict(fit)
#' predict(fit, newdata = data.frame(wt = c(2.5, 3.5), hp = c(100, 200)))
predict.my_lm <- function(object, newdata, ...) {
  if (missing(newdata)) {
    return(object$fitted.values)
  }
  else {
    tt <- delete.response(object$terms)
    mf <- model.frame(tt, newdata)
    X <- model.matrix(tt, mf)
    preds <- drop(X %*% object$coefficients)
    return(preds)
  }
}

?my_lm
