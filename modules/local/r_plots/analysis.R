#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(ggnewscale)
})

# arguments: tree, metadata CSV, output prefix
args <- commandArgs(trailingOnly = TRUE)
tree_file     <- args[1]
tree_metadata <- args[2]
out_prefix    <- args[3]
rgi_files  <- args[-(1:3)]

#
# PLOTTING PHYLOGENETIC TREE
#

# column in the metadata that matches the tree tip names
# (change to "sample" if your genome files are named by sample, e.g. cp1032.fa)
key_col <- "accession"

# read tree and strip file extensions from tip labels
tree <- read.tree(tree_file)
tree$tip.label <- sub("\\.(fa|fna|fasta|fas|fsa)$", "", tree$tip.label)

# read metadata and match it to the tips
meta <- read.csv(tree_metadata, stringsAsFactors = FALSE)
idx  <- match(tree$tip.label, meta[[key_col]])

anno <- data.frame(
  Country = meta$country[idx],
  Year    = as.numeric(meta$date_of_collection[idx]),
  row.names = tree$tip.label
)

# base tree
p <- ggtree(tree) +
  geom_tiplab(size = 2, align = TRUE, linetype = "dotted", linesize = 0.2) +
  geom_treescale() +
  theme_tree2()

# spacing is in tree-length units, so scale it to the tree
depth  <- max(node.depth.edgelength(tree))
off1   <- depth * 0.35   # gap between tree tips and first strip (room for labels)
colw   <- 0.05           # strip width as a fraction of tree width
off2   <- off1 + depth * colw + depth * 0.01

# strip 1: country (categorical)
p <- gheatmap(p, anno[, "Country", drop = FALSE],
              offset = off1, width = colw,
              colnames = FALSE, color = NA) +
  scale_fill_discrete(name = "Country", na.value = "grey85")

# strip 2: year (continuous), with its own scale
p <- p + new_scale_fill()
p <- gheatmap(p, anno[, "Year", drop = FALSE],
              offset = off2, width = colw,
              colnames = FALSE, color = NA) +
  scale_fill_viridis_c(name = "Year", na.value = "grey85")

# save
h <- max(4, length(tree$tip.label) * 0.15)
ggsave(paste0(out_prefix, "_tree", ".pdf"), p, width = 10, height = h, limitsize = FALSE)
ggsave(paste0(out_prefix, "_tree", ".png"), p, width = 10, height = h, dpi = 300, limitsize = FALSE)



#
# PLOTTING AMR PRESENCE/ABSENCE HEATMAP
#

sample_of <- function(f) sub("\\.txt$", "", sub("\\.rgi\\.txt$", "", basename(f)))

# read all RGI outputs (tab-delimited)
hits <- do.call(rbind, lapply(rgi_files, function(f) {
  d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(NULL)
  data.frame(sample = sample_of(f),
             gene   = d$Best_Hit_ARO,
             cutoff = d$Cut_Off,
             stringsAsFactors = FALSE)
}))

# if a sample has several hits for the same gene, keep the best one
rank <- c(Perfect = 3, Strict = 2, Loose = 1)
hits <- hits[order(-rank[hits$cutoff]), ]
hits <- hits[!duplicated(hits[, c("sample", "gene")]), ]

samples <- sapply(rgi_files, sample_of)
genes   <- names(sort(table(hits$gene), decreasing = TRUE))

# matrix of calls: Absent / Loose / Strict / Perfect
mat <- matrix("Absent", nrow = length(samples), ncol = length(genes),
              dimnames = list(samples, genes))
mat[cbind(match(hits$sample, samples), match(hits$gene, genes))] <- hits$cutoff

fill_scale <- scale_fill_manual(
  values = c(Absent = "grey92", Loose = "#fdd49e",
             Strict = "#fc8d59", Perfect = "#b30000"),
  breaks = c("Perfect", "Strict", "Loose", "Absent"),
  name   = "RGI hit")

#plotting heatmap
long <- expand.grid(sample = rownames(mat), gene = colnames(mat),
                    stringsAsFactors = FALSE)
long$call <- factor(mat[cbind(long$sample, long$gene)],
                    levels = c("Perfect", "Strict", "Loose", "Absent"))
long$gene   <- factor(long$gene, levels = genes)
long$sample <- factor(long$sample, levels = rev(samples))

p_amr <- ggplot(long, aes(gene, sample, fill = call)) +
  geom_tile(color = "white") + fill_scale +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(), panel.grid = element_blank())

h <- max(4, nrow(mat) * 0.15) + 2
w <- 3 + ncol(mat) * 0.25

ggsave(paste0(out_prefix, "_AMR.pdf"), p_amr, width = w, height = h, limitsize = FALSE)
ggsave(paste0(out_prefix, "_AMR.png"), p_amr, width = w, height = h, dpi = 300, limitsize = FALSE)