### Colorful Phylogenetic Tree
This R based code used the ".nwk" tree as input to make the large dataset phylogenetic tree

# Packages
# install.packages(c("tidyverse","ggtree","ape"))  # run once if needed
library(tidyverse)
library(ggtree)
library(ape)
library(phangorn) 


# 1) Inputs file
tree_file   <- "phylogenetic_tree.nwk"   # your Newick tree
ploidy_file <- "genotype_ploidy.csv"  # columns: ploidy, genotypes


## example of genotype_ploidy.csv
# ploidy,genotypes
# 2,genotype_1
# 2,genotype_2
# 2,genotype_3
# 2,genotype_4



# 2) Read tree and midpoint-root (if not already)
tr <- read.tree(tree_file)


# 3) If the tree isn't rooted or you want consistent display
tr <- midpoint(tr)


# 4) Read ploidy map
pld <- readr::read_csv(ploidy_file, show_col_types = FALSE)


# 5) Normalize columns
stopifnot(all(c("ploidy","genotypes") %in% names(pld)))
pld <- pld %>%
  mutate(
    ploidy = str_trim(as.character(ploidy)),  # remove spaces
    ploidy = as.integer(ploidy)  # convert to integers
  ) %>%
  mutate(
    PloidyClass = case_when(
      ploidy == 2 ~ "Diploid (2x)",
      ploidy == 3 ~ "Triploid (3x)",
      ploidy == 4 ~ "Tetraploid (4x)",
      ploidy == 5 ~ "Pentaploid (5x)",
      TRUE        ~ "Unknown"
    )
  )


# 6) Build a data frame of tree tips
tips_df <- tibble(label = tr$tip.label)


# 7) Join ploidy onto tips by label <-> genotypes
annot <- tips_df %>%
  left_join(pld, by = c("label" = "genotypes")) %>%
  mutate(
    PloidyClass = replace_na(PloidyClass, "Unknown"),
    ploidy      = replace_na(ploidy, NA_integer_)
  )


# 8) Quick match report in console
cat("Total tips in tree: ", nrow(tips_df), "\n")
cat("Matched ploidy rows:", sum(!is.na(annot$ploidy)), "\n")
if (any(is.na(annot$ploidy))) {
  cat("Unmatched tips (no ploidy):\n")
  print(annot$label[is.na(annot$ploidy)])
}


# 9) Plotting colors (set however you like)
ploidy_cols <- c(
  "Diploid (2x)"    = "#1f77b4",
  "Triploid (3x)"   = "#ff7f0e",
  "Tetraploid (4x)" = "#d62728",
  "Pentaploid (5x)" = "purple",
  "Unknown"         = "grey70"
)


# 10a) option1: set the dimension of tree for clean circular tree
tree <- ggtree(tr, layout = "circular", aes(color = PloidyClass), size = 0.10) %<+% annot +
  xlim_tree(1.4) +  # increase distance to make room
  geom_tippoint(aes(fill = PloidyClass),
                shape = 21, size = 1.2, stroke = 0.15, color = "black", alpha = 0.8) +
  geom_tiplab2(aes(angle = angle), align = TRUE, offset = 0.08, size = 0.8) +  # smaller & farther
  scale_fill_manual(values = ploidy_cols, name = "Ploidy") +
  scale_color_manual(values = ploidy_cols, name = "Ploidy", na.value = "grey80") +
  guides(fill = guide_legend(override.aes = list(size = 2)),
         color = guide_legend(override.aes = list(size = 1.5))) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 6),
    plot.margin  = margin(10, 16, 10, 10)
  )


# 10b) option2: cicular tree for tight branch to genotype readable tree
p <- ggtree(
  tr,
  layout = "circular",
  branch.length = "branch.length",   # ensure real branch lengths
  aes(color = PloidyClass),
  size = 0.10
) %<+% annot +
  
  # give more radial space
  xlim_tree(1.8) +
  
  # branch tip points
  geom_tippoint(
    aes(fill = PloidyClass),
    shape = 21,
    size = 1.0,
    stroke = 0.15,
    color = "black",
    alpha = 0.9
  ) +
  
  # genotypes labeling to branch tip
  geom_tiplab(
    size = 0.55,
    offset = 0.02,          # label just outside the tip
    linesize = 0.15,        # thin leader line
    align = FALSE
  ) +
  
  scale_fill_manual(values = ploidy_cols, name = "Ploidy") +
  scale_color_manual(values = ploidy_cols, name = "Ploidy", na.value = "grey80") +
  
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    plot.margin  = margin(10, 30, 10, 10)
  )


# 11) save the figures in pdf and jpg format
ggsave("phylogeny_circular_dense.png", tree, width = 24, height = 24, dpi = 600)
ggsave("phylogeny_circular_dense.pdf", tree, width = 24, height = 24)

