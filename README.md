### Neighbor-joining phylogenetic tree

```r
# 1) install.packages(c("tidyverse","ape","phangorn","proxy","parallelDist"))
library(tidyverse)
library(ape)
library(phangorn)
library(proxy)
library(parallelDist)
```
```r
# 2) upload the markers_data.csv file
dos <- read_csv("markers_data.csv", col_types = cols(.default = "c")) # makers_data.csv is attached in the repository
```

```r
# 3) Compute NA‐aware Euclidean distance & build NJ tree
dm_euc <- as.dist(proxy::dist(
  doseMat,
  method   = "Euclidean",
  pairwise = TRUE    # ignores any NAs
))
# for neighbor-joining tree
nj_euc <- NJ(dm_euc)
nj_mid <- midpoint(nj_euc)  # midpoint‐root for symmetry

# write tree with support labels for MEGA/R
write.tree(nj_mid, file = "roses_NJ_midpoint_bootstrap.nwk")
```
```r
# 4) Follow this for making colorful phylogenetic tree
https://github.com/RezwanCAAS/Circular_Colorful_Phylogenetic_Tree
```
