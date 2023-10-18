# GGoutlieR

`GGoutlieR` (**G**eo-**G**enetic **outlieR**) is an R package for identifying and visualizing unusual geo-genetic patterns for landscape genomic studies.
It provides a summary table with heuristic p-values allowing users to easily identify outliers from thousands of biological samples.
Moreover, its visualization tool enables users to present unusual geo-genetic association patterns on a geographical map.

#### Details about GGoutlieR framework

To find the detailed algorithm of `GGoutlieR`, please check the supplementary material of our manuscript [HERE](https://github.com/kjschmidlab/GGoutlieR/blob/master/paper/suppinfo.pdf) (`./paper/suppinfo.pdf`).

## Installation

You can install `GGoutlieR` either from CRAN or from this Github repository with the R commands below.

```
##-----------------------
## install from CRAN
##-----------------------
install.packages("GGoutlieR")

##-----------------------
## install from GitHub
##-----------------------

## install dependencies
install.packages(c("stats4","FastKNN","foreach","doParallel","parallel","scales", "RColorBrewer","ggforce", "rlang", "tidyr", "utils", "rnaturalearth", "rnaturalearthdata", "sf", "ggplot2","cowplot"))

## install GGoutlieR from Github
library(devtools)
install_git("https://github.com/kjschmidlab/GGoutlieR.git")
```

## Get started

You can find a more detailed tutorial with an example of `GGoutlieR` using a global barley landrace collection [HERE](https://github.com/kjschmidlab/GGoutlieR/blob/master/vignettes/outlier_detection.pdf) (`./vignettes/outlier_detection.pdf`). 
The VCF file of the barley landrace collection is available [HERE](https://github.com/kjschmidlab/GGoutlieR/blob/master/supp_data/IPKbarley_GBS1661landraces_LDpruned.vcf.gz) (`./supp_data/IPKbarley_GBS1661landraces_LDpruned.vcf.gz`).

### Simple example

```
# example of barley landrace data 
library(GGoutlieR)
data("ipk_anc_coef") # get ancestry coefficients
data("ipk_geo_coord") # get geographical coordinates

# run ggoutlier (takes about 40 sec)
ggoutlier_example <- 
                ggoutlier(geo_coord = ipk_geo_coord,
                          gen_coord = ipk_anc_coef,
                          make_fig = FALSE,
                          p_thres = 0.005, 
                          cpu = 4, 
                          klim = c(3,15), 
                          method = "composite",
                          verbose = FALSE,
                          multi_stages = FALSE)) # switch off multi-stage test to reduce computational time
                          
# check the first few rows of the summary table
head(summary_ggoutlier(ggoutlier_example))     

# visualize GGoutlieR results (set boundaries using `plot_xlim` and `plot_ylim` to focus on Eurasia and North Africa)

plot_ggoutlier(ggoutlier_res = ggoutlier_example,
               gen_coord = ipk_anc_coef,
               geo_coord = ipk_geo_coord,
               p_thres = 0.005,
               map_type = "both",
               plot_xlim = c(-20,140),
               plot_ylim = c(10,62),
               pie_r_scale = 1.8,
               map_resolution = "medium")
```

## Reference manual

The reference manual is at  
[https://cran.r-project.org/web/packages/GGoutlieR/GGoutlieR.pdf](https://cran.r-project.org/web/packages/GGoutlieR/GGoutlieR.pdf)

## Contributing

We appreciate your interest in using `GGoutlieR` in your study. 
You can contribute to the improvement of the package via following ways.

### Questions and bug reports

Please open issues if you notice an issue with the `GGoutlieR` repository. 
You can also contact [cheweichang92@gmail.com](cheweichang92@gmail.com) or [karl.schmid@uni-hohenheim.de](karl.schmid@uni-hohenheim.de) for your questions.

### Pull request

If you would like to contribute to the code:

* Fork the `GGoutlieR` repository
* Contribute to your forked repository.
* Create a pull request.

Your changes or additions will be merged in the master branch of `GGoutlieR` repository if they pass required checks.

## Preprint
Che-Wei Chang and Karl Schmid. 2023. GGoutlieR: an R package to identify and visualize unusual geo-genetic patterns of biological samples. bioRxiv. DOI: [https://doi.org/10.1101/2023.04.06.535838](https://doi.org/10.1101/2023.04.06.535838)



