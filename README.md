
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


