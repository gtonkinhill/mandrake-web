---
title: "Algorithm"
date: 2018-01-27T15:42:17+01:00
weight: 20
---

<br></br>
<br></br>
<br></br>

### Mandrake

Mandrake combines the fast pairwise genetic distance methods used in **[pairsnp](https://github.com/gtonkinhill/pairsnp-cpp)** and **[poppunk](https://github.com/johnlees/pp-sketchlib)** with the Stochastic Cluster Embedding algorithm to provide fast and informative dimension reduction plots of genetic population structure.

#### Stochastic Cluster Embedding

Neighbor Embedding (NE) that aims to preserve pairwise similarities between data items has been shown to yield an effective principle for data visualization. However, even the currently best NE methods such as Stochastic Neighbor Embedding (SNE) may leave large-scale patterns such as clusters hidden despite strong signals being present in the data. Stochastic Cluster Embedding (SCE) generalizes SNE by using non-normalized Kullback-Leibler divergence with a scale parameter. In this family, much better cluster visualizations often appear with a scale parameter value different from the one corresponding to SNE. Experimental results have demonstrated that this method consistently improves the visualization of data clusters compared with the state-of-the-art NE approaches. The original code for the SCE algorithm can be found at https://github.com/rozyangno/sce.

