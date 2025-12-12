# Software installation

Please use the instructions below to install and check your
installation. If you run into issues, you can post a question with
information about what you tried and what didn’t work, on the [course
github discussion
page](https://github.com/inlabru-org/inlabruCourseMay2025/discussions)
and Finn will reply when he’s able. Since 29 April 2025, the latest INLA
package is built for R 4.5, so if you’re able to upgrade your R
installation, please do so to avoid unnecessary issues. The package will
in many cases also work with older R versions, but compatibility is
sometimes difficult.

## Installing `INLA` and `inlabru`

Due to the work involved in building the binaries for the INLA package C
software for different architectures, the INLA package is not on CRAN,
but it can be installed from its own distribution repository. The
easiest approach is to add the repository to your local R options,
e.g. in `~/.Rprofile`, so that `install.packages("INLA")` can work as
usual. Alternatively, the same repository information can be supplied
via `install.packages("INLA", repos = ...)`.

In `~/.Rprofile`, add the following lines, or run the code in your R
session:

``` r
local({
  r <- c(INLA = "https://inla.r-inla-download.org/R/testing",
         CRAN = "https://cloud.r-project.org/",
         inlabru_universe = "https://inlabru-org.r-universe.dev")
  options(repos = r)
})
```

Note: if you have installed `INLA` from the `"stable"` repository
instead of `"testing"`, you may need to upgrade to the `"testing"`
version, which typically contains bug fixes and improvements, and is the
version `inlabru` is tested against.

The third repository is optional, and gives you easy access to the
development versions of `inlabru` (2.13.0.9024 or later) and `fmesher`
(0.6.0.9000 or later, which is a package dependency of `INLA` and
`inlabru`). Most of the tutorials in this course will work with the CRAN
versions of `inlabru` and `fmesher`, `2.12.0`, but installing the
development version from `r-universe.dev` is recommended, and required
for a few new features and contains bug fixes, see [the inlabru
changelog](https://inlabru-org.github.io/inlabru/news/).

To install the packages, including optional dependencies many of which
are required for data analysis, run the following code in your R
session:

``` r
install.packages(c("INLA", "inlabru"), dependencies = TRUE)
```

The warning about the unavailable `HKprocess` package can safely be
ignored. `INLA` also suggests two Bioconductor packages, `graph` and
`Rgraphviz`, which are not essential for this course, but can be
installed with

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("graph", "Rgraphviz"), dep = TRUE)
```

## Installation check

Please check your installation using the basic model runs below. If you
run into issues, you can post a question with information about what you
tried and what didn’t work, on the [course github discussion
page](https://github.com/inlabru-org/inlabruCourseMay2025/discussions)
and Finn will reply when he’s able. Since 29 April 2025, the latest INLA
package is built for R 4.5, so if you’re able to upgrade your R
installation, please do so to avoid unnecessary issues. The package will
in many cases also work with older R versions, but compatibility is
sometimes difficult.

You can check that `INLA` is correctly installed by running

``` r
df <- data.frame(y = rnorm(100) + 10)
fit <- INLA::inla(
  y ~ 1,
  data = df
)
summary(fit)
#> Time used:
#>     Pre = 0.382, Running = 0.121, Post = 0.0141, Total = 0.517 
#> Fixed effects:
#>               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
#> (Intercept) 10.071 0.105      9.865   10.071     10.277 10.071   0
#> 
#> Model hyperparameters:
#>                                          mean   sd 0.025quant 0.5quant
#> Precision for the Gaussian observations 0.924 0.13      0.687    0.918
#>                                         0.975quant  mode
#> Precision for the Gaussian observations       1.20 0.905
#> 
#> Marginal log-Likelihood:  -158.73 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
```

If the simple [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html) call
fails with a crash, you may need to install different `inla` binaries
for your hardware/software combination, with
[`INLA::inla.binary.install()`](https://rdrr.io/pkg/INLA/man/binary.install.html).

When [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html) works, you can
check that `inlabru` is installed correctly by running the same model in
`inlabru`:

``` r
fit <- inlabru::bru(
  y ~ Intercept(1, prec.linear = exp(-7)),
  data = df
)
summary(fit)
#> inlabru version: 2.13.0.9024 
#> INLA version: 25.12.12 
#> Latent components:
#> Intercept: main = linear(1)
#> Observation models:
#>   Model tag: <No tag>
#>     Family: 'gaussian'
#>     Data class: 'data.frame'
#>     Response class: 'numeric'
#>     Predictor: y ~ Intercept
#>     Additive/Linear/Rowwise: TRUE/TRUE/TRUE
#>     Used components: effect[Intercept], latent[] 
#> Time used:
#>     Pre = 0.226, Running = 0.12, Post = 0.033, Total = 0.379 
#> Fixed effects:
#>             mean    sd 0.025quant 0.5quant 0.975quant   mode kld
#> Intercept 10.071 0.105      9.865   10.071     10.277 10.071   0
#> 
#> Model hyperparameters:
#>                                          mean   sd 0.025quant 0.5quant
#> Precision for the Gaussian observations 0.924 0.13      0.687    0.918
#>                                         0.975quant  mode
#> Precision for the Gaussian observations       1.20 0.905
#> 
#> Marginal log-Likelihood:  -163.19 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
```
