[![DOI](https://zenodo.org/badge/465693977.svg)](https://zenodo.org/badge/latestdoi/465693977)

# SIFT-seq data for patient MR4050

R data package related to 

>Chandran et al. *Immunogenicity and therapeutic targeting of a public neoantigen derived from mutated PIK3CA* Nat Medicine (2022).

This package contains access to the processed TCR-seq and scRNA-seq data of T cells obtained from patient MR4050.

## How to use it

```
## install
devtools::install_github("abcwcm/KlebanoffMR4050")
```

Upon installation, the processed data, e.g. in the form of `SingleCellExperiment` objects, can be loaded thus:

```
## load SingleCellExperiment objects
sce.mr4050 <- KlebanoffMR4050::load_MR4050shared()

## load results of differential gene expression comparisons
## see 
KlebanoffMR4050::load_DE_results() # loads an object named `delist.both`
de.mr4050 <- delist.both; rm(delist.both)
```

For more details, see the [code repository](https://github.com/abcwcm/Chandran2021) and the `Rmd` files detailing how the data underlying the figures were obtained.

