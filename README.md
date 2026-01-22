# Neighborjoining_phylogenetic_tree

```r
# install.packages(c("tidyverse","adegenet","ape","poppr","ggtree","ggplot2"))
library(tidyverse)
library(ape)
library(phangorn)
library(poppr)
library(ggtree)
library(proxy)
library(parallel)
library(viridisLite)
library(parallelDist)
```
```r
# upload the markers_data.csv file
dos <- read_csv("markers_data.csv", col_types = cols(.default = "c"))
```

```r
# Compute NA‐aware Euclidean distance & build NJ tree
dm_euc <- as.dist(proxy::dist(
  doseMat,
  method   = "Euclidean",
  pairwise = TRUE    # ignores any NAs
))

nj_euc <- NJ(dm_euc)
nj_mid <- midpoint(nj_euc)  # midpoint‐root for symmetry

# write tree with support labels for MEGA/R
write.tree(nj_mid, file = "roses_NJ_midpoint_bootstrap.nwk")
```
```r
# Follow this for making colorful phylogenetic tree
https://github.com/RezwanCAAS/Circular_Colorful_Phylogenetic_Tree
```
