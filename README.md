# GGoutlieR

`GGoutlieR` is an R package for identifying and visualizing unusual geo-genetic patterns for landscape genomic studies.
It provides a summary table with heuristic p-values allowing users to easily identify outliers from thousands of biological samples.
Moreover, its visualization tool enables users to present unusual geo-genetic association patterns on a geographical map.

## Installation

To install `GGoutlieR` v1.0.0, please use the R commands below:

```
## install dependencies
install.packages(c("geosphere","stats4","FastKNN","foreach","doSNOW","iterators","parallel","rworldmap","scales","plotrix","mapplots", "RColorBrewer","rworldxtra","dichromat","sp"))

## install GGoutlieR
library(devtools)
install_git("https://github.com/kjschmidlab/GGoutlieR.git")
```

## Tutorial

You can find an example of `GGoutlieR` analysis using a global barley landrace collection [HERE](https://github.com/kjschmidlab/GGoutlieR/blob/master/vignettes/outlier_detection.pdf) (`./vignettes/outlier_detection.pdf`). 

## Details about GGoutlieR framework

To find the detailed algorithm of `GGoutlieR`, please check the supplementary material of our manuscript [HERE](https://github.com/kjschmidlab/GGoutlieR/blob/master/paper/suppinfo.pdf) (`./paper/suppinfo.pdf`).

## Preprint
Che-Wei Chang and Karl Schmid. 2023. GGoutlieR: an R package to identify and visualize unusual geo-genetic patterns of biological samples. bioRxiv. DOI: [https://doi.org/10.1101/2023.04.06.535838](https://doi.org/10.1101/2023.04.06.535838)
