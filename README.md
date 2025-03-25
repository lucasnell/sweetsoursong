
# A Song of Sweet and Sour (Nectar)


[![DOI](https://zenodo.org/badge/768350141.svg)](https://zenodo.org/badge/latestdoi/768350141)


This repo contains the R package `sweetsoursong` that simulates metacommunities
of nectar microbes with effects on pollinator visitations to plants.
It also contains scripts (many of which use the R package) that are used to
create figures in the manuscript *Dispersal–community feedback can promote 
regional coexistence despite priority effects*.


# Organization

This repo contains the following top-level files/folders:

```
├── _data
├── _figures
├── _mathematica
├── _scripts
├── DESCRIPTION
├── LICENSE.md
├── man
├── NAMESPACE
├── R
├── README.md
├── src
└── sweetsoursong.Rproj
```

Here are descriptions of each:


* `_data`: Folder containing `*.rds` and `*.mx`files that store output 
    from R and Mathematica simulations that take a while to run.
    None are present in this repo, but they can be created using the
    provided scripts.
* `_figures`: Folder used to store figures created using R and Mathematica
    scripts. None are present in this repo, but they can be created using
    the provided scripts.
* `_mathematica`: Folder containing Mathematica scripts (and rendered PDFs)
   that create some of the plots used in the manuscript.
   Their file name indicates what it's used for, as does the description at
   the top of each document.
* `_scripts`: Folder containing R scripts that create most of the plots used
    in the manuscript. They are ordered by use in the paper, and
    each file contains a description of what it contains.
* `DESCRIPTION`: Description file for the `sweetsoursong` package.
* `LICENSE.md`: License file for the `sweetsoursong` package.
* `man`: Folder containing the documentation files for the `sweetsoursong` package.
* `NAMESPACE`: File defining imports and export for the `sweetsoursong` package.
* `R`: Folder containing the R files for the `sweetsoursong` package.
* `README.md`: This file, the top-level README.
* `src`: Folder containing the C++ files for the `sweetsoursong` package.
* `sweetsoursong.Rproj`: RStudio project file for this project.






# Replicating R environment

I used R version 4.4.3 (platform: aarch64-apple-darwin20) for all my scripts.

This project uses the `renv` package, so if you want to use this, you must
first install it:

```r
install.packages("renv")
```

Then to install all the packages I used for these analyses, you can simply run
the following while having this project's main directory as your working
directory:

```r
renv::restore()
```


If you'd rather avoid `renv`, then you can install all the packages 
(in the versions I used) this way:

```r
pkgs <- c("BH@1.87.0.1", "dplyr@1.1.4", "dqrng@0.4.1", "ggplot2@3.5.1", 
          "ggtext@0.1.2", "patchwork@1.3.0", "Rcpp@1.0.14", 
          "RcppArmadillo@14.4.0.1", "RcppParallel@5.1.10", 
          "tibble@3.2.1", "tidyverse@2.0.0")
install.packages(pkgs)
install.packages("remotes")
remotes::install_github("lucasnell/sweetsoursong@v1.0.0")
```

Note that these only install the proper versions of the packages I manually 
installed (or were dependencies of `sweetsoursong`),
so dependencies might vary from what I used.


Similarly, but without version numbers at all:

```r
pkgs <- c("BH", "dplyr", "dqrng", "ggplot2", 
          "ggtext", "patchwork", "Rcpp", 
          "RcppArmadillo", "RcppParallel", 
          "tibble", "tidyverse")
install.packages(pkgs)
install.packages("remotes")
remotes::install_github("lucasnell/sweetsoursong")
```
