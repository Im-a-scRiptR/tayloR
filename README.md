
# tayloR

The goal of tayloR is to perform the Taylor Series expanion for a single varibaled
function.

## Installation

You can install the development version of tayloR like so:

``` r
# To install

devtools::install_github("Im-a-scRiptR/tayloR")

```

## Example

This is a basic demo usage of the make_taylor_series function:

``` r
library(tayloR)

f <- \(x) {
  exp(x^2)
}

taylor <-
  make_taylor_series(
    f,
    num_derivs   = 3,
    center       = 0,
    intr_start   = 0,
    intr_end     = 0.1,
    var_name     = "x",
    open_closed  = "closed"
  )

taylor

```

