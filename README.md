weierstrass
===========

An R implementation of Weierstrass rejection sampler for combining posterior samples from multiple subsets.

#Description
The Weierstrass sampler is a "divide-conquer-combine" type parallel MCMC sampler. The implementation in this package is the rejection sampler which directly works on the subset posterior samples as a post-processing sampler. The method make use of rejection sampling or importance sampling to combine the subset posterior samples for approximating the posterior obtained on full data set.

To combine subset posterior samples, the algorithm adopts the 'pairwise-combining' strategy, i.e, it first combines the subset pairwisely to obtain half numbers of new subsets and then repeats the procedure until obtaining the final one. More detailed information is provided in the R-help documentation accompanied with functions.

Besides the main function, two testing models are also included in the pacakge. Uers can build logistic model and binomial model with user-specified features (such as number of predictors, predictor correlations, etc) for testing the performance of the weierstrass rejection sampler. The testing functions integrate full functionalities of data generating, inference and combining. The logistic model requires the package `BayesLogit` (available on R-CRAN) for posterior sampling. More detailed information can be found in the R-help documentation.

#Installation
There are several simple ways to install the package for R. Building the package requires `devtools`. Just type the following lines in your R console (Suppose your current working directory is where you want to keep the package),
```
library(devtools)
install_github('weierstrass','wwrechard')
```
Now the package is fully functional.

Alternatively, one can download the whole repository manually (in one folder and assume the folder is named "weierstrass") into the target directory, or obtain via `git`,
```
cd your_target_directory
git clone https://github.com/wwrechard/weierstrass.git
```
and then in the R console (Suppose your current working directory is where you want to keep the package),
```
library(devtools)
install('weierstrass')
```
However, this way might require a restart of R to correctly view the R-help documentation for the functions.

Finally, if you are not interested at all in the R-help documentation, you can simply download the "weierstrass.r" under "./R" on this repository. This is the main function for weierstrass rejection sampling. "logitTest.r" and "BinTest.r" are the two testing functions comming along with this packages. All of the three can be directly `source` and used in R.

#Usage
To call up the package, type in the R console,
```
setwd(directory_that_weierstrass_is_installed)
library(weierstrass)
```
The detail of functions are documented in the R-help docs for different functions. Type
```
?weierstrass
```
will call up the help document for the function `weierstrass`. One simple linear regression example is provided with this function, to view that,
```
example('weierstrass')
```
`logitTest` and `BinTest` are two testing functions for Weierstrass rejection sampling. The details can be found with their help documents as well.

#References
Xiangyu Wang, David B. Dunson. Parallelizing MCMC via Weierstrass Sampler. [http://arxiv.org/abs/1312.4605]
