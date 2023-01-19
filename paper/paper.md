---
title: 'GGoutlieR: an R package to identify and visualize abnormal geo-genetic patterns of biological samples'
tags:
  - R
  - landscape genomics
  - K nearest neighbors
  - geo-genetic patterns

authors:
  - name: Che-Wei Chang
    orcid: 0000-0001-6704-3732
    equal-contrib: true
    affiliation: 1
  - name: Karl Schmid
    orcid: 0000-0001-5129-895X
    equal-contrib: true 
    affiliation: 1

affiliations:
 - name: University of Hohenheim, Stuttgart, Germany
   index: 1
 - name: Institution Name, Country

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

Landscape genomics is a thriving field in ecological conservation (CITATION) and evolutionary genetics (CITATION), which provides insights into associations between genetic variation and environmental factors.
This methodology requires reliable geographical and genomic information of biological samples.
To recognize whether data are reliable, researchers may scrutinize associations between genetic similarities and geographical origins of biological samples before carrying out further studies.
The pairwise genetic similarities of samples are expected to decline as geographical distances between origin habitats increase, so-called isolation-by-distance assumption.
This assumption could be violated due to long-distance migration or artificial factors, such as human activities or mistakes in data management.

Visualization-based tools, such as `SPA` [@yang2012model], `SpaceMix` [@bradburd2016spatial], `unPC` [@house2018evaluating], help researchers to unveil samples with geo-genetic patterns opposing the isolation-by-distance assumption, but those tools do not provide statistics to simply pinpoint outliers.
This shortage could be detrimental to the reproducibility of research.
Moreover, with the advances in genome sequencing technologies, researchers nowadays work on much larger sample sizes than in the past, for example, the studies of genebank collections (CITATION).
Visualization-based approaches may have difficulty in presenting unusual geo-genetic patterns because thousands of samples projected on a graph tend to overlap.
Therefore, a new approach is needed to facilitate the detection of unusual geo-genetic associations in biological samples.
We developed a heuristic statistic framework to detect **G**eo-**G**enetic **outlier**s, named `GGoutlieR`.
Our package `GGoutlieR` computes empirical p-values of violating the isolation-by-distance assumption for individual samples according to geographical origins and genotypic data.
This feature enables users to easily select outliers from thousands of samples for further investigation.
Furthermore, `GGoutlieR` visualizes the geo-genetic patterns of outliers in network fashion on a geographical map, providing insights into the relationships of geography and genetic clusters.

# Identification of outliers with abnormal geo-genetic patterns

Under the isolation-by-distance assumption, the geographical origins are predictable from genetic variations (CITATION), and vice versa.
With this respect, prediction models should result in large prediction errors for samples that oppose the isolation-by-distance assumption.
We developed the `GGoutlieR` framework with this idea to model abnormal geo-genetic patterns.

In brief, `GGoutlierR` uses KNN regression to predict genetic components with the K nearest geographical neighbors and also does prediction contrariwise.
Next, prediction errors are transformed to distance-based statistics following a Gamma distribution with unknown parameters.
An empirical Gamma distribution is obtained as the null distribution by searching optimal parameters with maximum likelihood estimation.
With the null Gamma distribution, `GGoutlieR` tests the null hypothesis that the geo-genetic pattern of a given sample agrees with the isolation-by-distance assumption.
Finally, p values for every sample are computed with the empirical null distribution and statistics computed from prediction errors.
The details of the `GGoutlieR` framework are described in the supplementary material (GITHUB_LINK).

# Visualization on a geographical map

# References
