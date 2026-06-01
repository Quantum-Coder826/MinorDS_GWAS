library(tidyverse)
library(circlize)
library(ComplexHeatmap)

qtl.table <- read_csv("outputs/qtl_table.csv")
chrom_lenghts <- read_csv("data/potato_chrom_lengths_DMv6.csv")

###################################
### generate genome coordinates ###
###################################
qtl.table$coord <- paste0(qtl.table$Chrom, ":", qtl.table$Position, "..", qtl.table$Position)


circos.genomicInitialize(chrom_lenghts)

circos.trackPoints(
  factors = qtl.table$Chrom,
  x = qtl.table$Position,
  y = rep(0, nrow(qtl.table)),
  pch = 3, cex = 1.2
)

circos.trackText(
  factors = qtl.table$Chrom,
  x = qtl.table$Position,
  y = rep(-1, nrow(qtl.table)),
  labels = qtl.table$Marker,
  facing = "clockwise",
  cex = 0.7
)

circos.clear()
