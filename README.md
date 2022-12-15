# *herbivar*
<img src="man/figures/logo.png" align="right" height="200"/>

### Quantitative Tools For Variance Explicit Ecology and Herbivory



## Overview 

This package makes several methods used for analyzing variability inecological data available. Currently, it supports fitting a neutral herbivory model to data using Maximum Likelihood Estimation (MLE). It also contains several helper functions for investigating distributions and wrangling the HerbVar dataset. In addition, using the *imager* infrastructure, the package also supports automatic processing of leaf images, mostly in the context of leaf scans, that can retrieve rich within-leaf herbivory data. The package is integrated with the *spatstat* package that allows for spatial analysis of these data. 
    
Currently work in progress, mostly for internal use only. Proceed with caution :)

Report any bugs to vsbpan@gmail.com Thanks!

## Installation

To install the package, you'll need to first use the `install_github()` function from the *devtools* package to install the latest version of *herbivar* from github. The package depends on *EBImage* and some of the packages associated with *spatstat*, so you'll need to install those too. Installing them can be kind of tricky, so I have provided the code below. If you not using Windows, you might have a problem with the file `libX11.6.dylib`, which the dependent package *imager* uses. You'll need to install XQuartz (https://www.xquartz.org/) to fix the error. 


```{r}
  #install.packages("devtools")
  #install.packages("BiocManager") Needed for BiocManager::install
  #BiocManager::install("EBImage")
  #install.packages("spatstat")
  devtools::install_github("vsbpan/herbivar", build_opts = c("--no-resave-data", "--no-manual"))
```
## Vignette

You can access the vignette via this command:

```{r}
browseVignettes("herbivar")
```


