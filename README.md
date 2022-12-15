# *herbivar*
<img src="man/figures/logo.png" align="right" height="200"/>

### Quantitative Tools For Variance Explicit Ecology and Herbivory



## Overview 

This package makes available several methods used for analyzing variability in 
    ecological data. Currently, it supports fitting a neutral herbivory model to data using 
    Maximum Likelihood Estimation (MLE). It also contains several helper functions for investigating 
    distributions and wrangling the HerbVar dataset. In addition, using the 'imager' infrastructure, 
    the package also supports automatic processing of leaf images, mostly in the context of leaf 
    scans, that can retrieve rich within-leaf herbivory data. The package is integrated with 
    the 'spatstat' package that allows for spatial analysis of these data. 
    
Currently work in progress, mostly for internal use only. Proceed with caution :)


## Installation

To install the package, you'll need to first use the `install_github()` function from the *devtools* package to install the latest version of *herbivar* from github. 


```{r}
  #install.packages("devtools")
  devtools::install_github("vsbpan/herbivar")
```
## Vignette

You can access the vignette via this command:

```{r}
browseVignettes("herbivar")
```


