# Find the K in K-means by Parametric Bootstrap

Nina Zumel, [Win-Vector LLC](http://www.win-vector.com)

Code to explore the use of parametric bootstrap simulations to determine the appropriate *k* for K-means clustering a data set. The approach was inspired by the `boot.comp` function in [`mixtools`](https://cran.r-project.org/web/packages/mixtools/index.html), an R package for fitting and analyzing finite mixture models (see [this article](http://exploringdatablog.blogspot.com/2011/08/fitting-mixture-distributions-with-r.html) by Ron Pearson for an example of using `mixtools` to fit mixtures of gaussians).

The primary code is in `kcomp_functions.R`. The R markdown file `kcomp.Rmd` shows an example of running this code. The R markdown file `stepthrough.Rmd` steps through one iteration of the bootstrap simulation, and produces the graphs used in our blog post.

## Links and References

* Our blog post on this approach is [here](http://www.win-vector.com/blog/2016/02/finding-the-k-in-k-means-by-parametric-bootstrap/).

* A Shiny app that demonstrates this code interactively is [here](https://win-vector.shinyapps.io/kcompshiny/).

* We cover other approaches to estimating the number of clusters in Chapter 8 of our book [*Practical Data Science in R*](https://www.manning.com/books/practical-data-science-with-r) (Manning, 2014). This chapter is available as a free sample chapter, [here](https://manning-content.s3.amazonaws.com/download/e/dc31390-3cb7-49dd-ab02-937c1af1c2e1/PDSwR_CH08.pdf).



