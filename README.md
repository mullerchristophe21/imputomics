

``` r
# LOAD LIBRARY
devtools::install_github("mullerchristophe21/imputomics")
library(imputomics)

# GENERATE DATA
n <- 100
sd <- 1

data_random <- data.frame(x = rnorm(n=n,mean=0,sd=sd), y=rnorm(n=n,mean=0,sd=sd), z = rnorm(n=n,mean=0,sd=sd))
data_random$y <- data_random$x + data_random$y
data_random$z <- data_random$x + data_random$y + data_random$z

data_wiht_na <- data_random
data_wiht_na[1:10, 1] <- NA

# IMPUTATION
impute_mice_rf(missdf = data_wiht_na)[1:10,]
impute_mice_drf(missdf = data_wiht_na)[1:10,]
impute_mice_cart(missdf = data_wiht_na)[1:10,]
```
