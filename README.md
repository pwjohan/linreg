# linreg

linreg fits linear regression models by ordinary least squares, computing
the estimates directly from (XᵀX)⁻¹Xᵀy instead of calling `lm()`. It was
written for BIOSTAT 625 at the University of Michigan.

## Installation

```r
# install.packages("remotes")
remotes::install_github("pwjohan/linreg", build_vignettes = TRUE)
```

## Usage

```r
library(linreg)

fit <- my_lm(mpg ~ wt + hp, data = mtcars)
fit
#> Call:
#> my_lm(formula = mpg ~ wt + hp, data = mtcars)
#>
#> Coefficients:
#> (Intercept)          wt          hp
#> 37.22727012 -3.87783074 -0.03177295

summary(fit)
predict(fit, newdata = data.frame(wt = 3, hp = 150))
```

## Functions

- `my_lm()`: fits a linear model from a formula and data frame
- `summary()`: standard errors, t-tests, R², adjusted R², and the F-test
- `print()`: displays a fit or its summary in the same layout as `lm()`
- `predict()`: fitted values, or predictions for new data

## Comparison with lm()

The vignette shows that `my_lm()` matches `lm()` exactly on real data,
simulated data, and a model with a categorical predictor, and benchmarks
the two with `bench::mark()`. After installing, open it with:

```r
vignette("my_lm-intro", package = "linreg")
```

## License

MIT
