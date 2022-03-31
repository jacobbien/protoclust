
<!-- README.md is generated from README.Rmd. Please edit that file -->

# protoclust

This R package implements hierarchical clustering with minimax linkage
(a.k.a., prototype clustering). See [Bien, J. and Tibshirani, R. (2011)
“Hierarchical Clustering with Prototypes via Minimax
Linkage”](http://faculty.marshall.usc.edu/Jacob-Bien/papers/jasa2011minimax.pdf)
for more details.

Since `protoclust` is on CRAN, it can be easily installed using

``` r
install.packages("protoclust")
```

in R. Occasionally, this github version of the package will be more
up-to-date. The easiest way to install this version is by using the
[devtools](https://CRAN.R-project.org/package=devtools) R package (if
not already installed, open R and type `install.packages("devtools")`).
To install `protoclust`, type

``` r
devtools::install_github("jacobbien/protoclust")
```

in R.

## Simple example of prototype clustering

Suppose we have some data. ![](README_figs/data-1.png)<!-- -->

After computing an n-by-n matrix of dissimilarities `d`, we can apply
prototype clustering.

``` r
# perform minimax linkage clustering:
library(protoclust)
hc <- protoclust(d)

# cut the tree at height 1:
h <- 1
cut <- protocut(hc, h = h)
# plot dendrogram (and show cut):
plot(hc, imerge = cut$imerge, col=2)
abline(h = h, lty = 2)
```

![](README_figs/dendrogram-1.png)<!-- -->

This gives us 9 prototypes with the guarantee that each of the n points
is no more than 1 unit from one of the prototypes.

![](README_figs/dataclusters-1.png)<!-- -->
