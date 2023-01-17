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

# Statement of need

Landscape genomics is a thriving field in ecological conservation (CITATION) and evolutionary genetics (CITATION).
It provides insights into associations between genetic variation and environmental factors.
This methodology requires reliable geographical and genomic information for biological samples.
The pairwise genetic similarities of biological samples are expected to decline as geographical distances between origin habitats increase, so-called isolation-by-distance.
This assumption could be violated due to long-distance migration or artificial factors, such as human transportation or mistakes in data management.
Thus, researchers may wish to scrutinize samples with the isolation-by-distance assumptions before carrying out further studies.

Visualization-based tools help researchers to study unusual associations between genetic similarities and geographical origins, such as `SPA` [@yang2012model], `SpaceMix` [@bradburd2016spatial], `unPC` [@house2018evaluating], but they do not provide simple statistics to point out outliers. 
Moreover, with the advances in genome sequencing technologies, researchers nowadays work on much larger sample sizes than in the past, for example, the studies of genebank collections (CITATION).
Visualization-based approaches may have difficulty in presenting unusual geo-genetic patterns because thousands of samples projected on a figure tend to overlap.
Therefore, a new approach is needed to assist researchers to detect unusual geo-genetic associations in biological samples.
We developed a heuristic statistic framework to detect **G**eo-**G**enetic **outlier**s, named `GGoutlieR`.
Our package `GGoutlieR` computes empirical p-values of violating the isolation-by-distance assumption for individual samples according to geographical origins and genotypic data.
`GGoutlieR` can also visualize the geo-genetic patterns of outliers on a geographical map.

# Identification of outliers with abnormal geo-genetic patterns

Under the isolation-by-distance assumption, the geographical origins are predictable from genetic variations (CITATION), and vice versa.
With this respect, prediction models should result in large prediction errors for samples that oppose the isolation-by-distance assumption.
In brief, `GGoutlierR` uses KNN regression to predict genetic components with geographical coordinates and also do prediction contrariwise.
Next, two empirical null Gamma distributions are obtained by searching optimal parameters with maximum likelihood estimation.
Finally, p values for every sample are computed with the empirical null distributions and statistics computed from prediction errors.
The details of `GGoutlieR` framework are described in the supplementary material (GITHUB_LINK).




# Visualization on a geographical map

# References
