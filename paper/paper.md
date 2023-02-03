---
title: 'GGoutlieR: an R package to identify and visualize unusual geo-genetic patterns of biological samples'
tags:
  - R
  - landscape genomics
  - K nearest neighbors
  - geo-genetic patterns

authors:
  - name: Che-Wei Chang
    orcid: 0000-0001-6704-3732
    affiliation: 1
  - name: Karl Schmid
    orcid: 0000-0001-5129-895X
    affiliation: 1

affiliations:
 - name: University of Hohenheim, Stuttgart, Germany
   index: 1
 

date: 09 Jan 2023
bibliography: paper.bib
---

# Summary

Landscape genomics is a rising research field integrating genomic and environmental information to explore driving forces of evolution.
Reliable geographical origin data of biological samples are a prerequisite for landscape genomics studies.
Conventionally, researchers discover potentially questionable samples with visualization-based tools.
However, such approaches are infeasible to handle large sample sizes due to overlapping data points on a graph and may encumber reproducible research.
To address this shortage, we developed **G**eo-**G**enetic **outlier** (`GGoutlieR`), an R package of a heuristic framework to reveal and visualize samples with unusual geo-genetic patterns.
`GGoutlieR` can calculate empirical p values for every sample, allowing researchers to easily spot outliers from thousands of samples.
Furthermore, `GGoutlieR` provides a plotting function to display the geo-genetic patterns of outliers on a geographical map.
`GGoutlieR` could greatly reduce the researcher's effort for data cleaning before conducting further landscape genomics analyses.

# Statement of need

Landscape genomics is a thriving field in ecological conservation and evolutionary genetics [@aguirre2021evolutionary; @lasky2023genotype], which provides insights into associations between genetic variation and environmental factors.
This methodology requires reliable geographical and genomic information of biological samples.
To recognize whether data are reliable, researchers may scrutinize associations between genetic similarities and geographical origins of biological samples before carrying out further studies.
The pairwise genetic similarities of samples are expected to decline as geographical distances between origin habitats increase, so-called isolation-by-distance assumption.
This assumption could be violated due to long-distance migration or artificial factors, such as human activities or mistakes in data management.

Visualization-based tools, such as `SPA` [@yang2012model], `SpaceMix` [@bradburd2016spatial], `unPC` [@house2018evaluating], help researchers to unveil samples with geo-genetic patterns opposing the isolation-by-distance assumption, but those tools do not provide statistics to simply pinpoint outliers.
This shortage could be detrimental to the reproducibility of research.
Moreover, with the advances in genome sequencing technologies, researchers nowadays work on increasing sample sizes, for example, genebank collection studies of rice [@gutaker2020genomic; @wang2018genomic], barley [@milner2019genebank], wheat [@schulthess2022genomics], soybean [@liu2020pan] and maize [@li2019identifying].
Visualization-based approaches may have difficulty in presenting unusual geo-genetic patterns in a large data set because numerous data points on a graph tend to overlap.
Therefore, a new approach is needed to facilitate the detection of unusual geo-genetic associations in biological samples.
We developed a heuristic statistic framework to detect **G**eo-**G**enetic **outlier**s, named `GGoutlieR`.
Our package `GGoutlieR` computes empirical p-values of violating the isolation-by-distance assumption for individual samples according to geographical origins and genotypic data.
This feature enables users to easily select outliers from thousands of samples for further investigation.
Furthermore, `GGoutlieR` visualizes the geo-genetic patterns of outliers in network fashion on a geographical map, providing insights into the relationships of geography and genetic clusters.

# Concept of `GGoutlieR`

Under the isolation-by-distance assumption, the geographical origins are predictable from genetic variations [@battey2020predicting;@guillot2016accurate], and vice versa.
With this respect, prediction models should result in large prediction errors for samples that oppose the isolation-by-distance assumption.
We developed the `GGoutlieR` framework with this idea to model anomalous geo-genetic patterns.

In brief, `GGoutlierR` uses KNN regression to predict genetic components with the K nearest geographical neighbors and also does prediction contrariwise.
Next, prediction errors are transformed to distance-based statistics following a Gamma distribution with unknown parameters.
An empirical Gamma distribution is obtained as the null distribution by searching optimal parameters with maximum likelihood estimation.
With the null Gamma distribution, `GGoutlieR` tests the null hypothesis that the geo-genetic pattern of a given sample agrees with the isolation-by-distance assumption.
Finally, p values for every sample are computed with the empirical null distribution and statistics computed from prediction errors.
The details of the `GGoutlieR` framework are described step-by-step in the supplementary material (GITHUB_LINK).

# Example

### Outlier identification

For demonstration, we used the genotypic data and passport data of the global barley landrace collection with 1,661 accessions from the IPK genebank [@milner2019genebank; @konig2020bridge].
The full analysis of the barley data set with `GGoutlieR` is available in the vignette (GITLAB_LINK).
The outlier analysis was done with the function `ggoutlier`.
The function `summary_ggoutlier` was then used to obtain a summary table of outliers by taking the output of `ggoutlier`.

```R
library(GGoutlieR)
data("ipk_anc_coef") # get ancestry coefficients
data("ipk_geo_coord") # get geographical coordinates

pthres = 0.025 # set a p-value threshold

## run GGoutlieR
ggoutlier_result <- ggoutlier(geo_coord = ipk_geo_coord,
                              gen_coord = ipk_anc_coef,
                              plot_dir = "./fig", 
                              p_thres = pthres, 
                              cpu = 4, 
                              klim = c(3,50), 
                              method = "composite",
                              verbose = F,
                              min_nn_dist = 1000)

## print out outliers
head(summary_ggoutlier(ggoutlier_result))

#>                 ID     method      p.value
#> 1  BRIDGE_HOR_2827     geoKNN 0.0002533251
#> 2 BRIDGE_HOR_12795     geoKNN 0.0002871882
#> 3    BRIDGE_BCC_37     geoKNN 0.0003011807
#> 4 BRIDGE_HOR_10557     geoKNN 0.0003500990
#> 5 BRIDGE_HOR_10555     geoKNN 0.0003697789
#> 6        BTR_FT519 geneticKNN 0.0003816026
```

### Visualization of unusual geo-genetic patterns

The unusual geo-genetic patterns detected by `GGoutlieR` can be presented on a geographical map with the function `plot_ggoutlier`.

```R
plot_ggoutlier(ggoutlier_res = ggoutlier_result,
               gen_coord = ipk_anc_coef,
               geo_coord = ipk_geo_coord,
               p_thres = pthres,
               map_type = "both",
               select_xlim = c(-20,140), 
               select_ylim = c(10,62),
               plot_xlim = c(-20,140),
               plot_ylim = c(10,62),
               pie_r_scale = 1.2,
               map_resolution = "course",
               adjust_p_value_projection = F)
```

![Visualization example of GGoutlieR with IPK barley landrace data. The red lines show the individual pairs with unusual genetic similarities across long geographical distances. The blue lines indicate the unusual genetic differences between geographical neighbors. Pie charts present the ancestry coefficients of outliers identified by GGoutlieR.](../fig/IPK_ggoutlier_for_paper.jpg)

# References
