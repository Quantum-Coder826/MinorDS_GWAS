library(tidyverse)
library(circlize)
library(ComplexHeatmap)

qtl.table <- read_rds("./outputs/qtl_table.rds")
qtl.ans <- read_rds("./outputs/fit_ans.rds")
chrom_lenghts <- read_csv("data/potato_chrom_lengths_DMv6.csv")

###################################
### generate genome coordinates ###
###################################
qtl.table$coord <- paste0(qtl.table$Chrom, ":", qtl.table$Position, "..", qtl.table$Position)

### maak per trait een tabel bepaal of ze signifikant zijn
trait.no = 1
fit.ans[[trait.no]] %>%
  filter(Trait == traits[trait.no], R2 >= 0.05) %>%
  arrange(desc(R2)) %>%
  knitr::kable(.,digits = 3)

qtl.table <- qtl.table %>%
  filter(Marker %in% c("SNP23778", "SNP34378", "SNP28271", "SNP33303", "SNP24920", "SNP22740"))


circos.genomicInitialize(chrom_lenghts,
                         sector.width = 2,
                         labels.cex = 1*par("cex"),
                         axis.labels.cex = 0.6*par("cex"))

model.shape <- c(
  "additive" = 15,
  "1-dom-alt" = 16,
  "1-dom-ref" = 17 
)

trait.colors <- c(
  "KOOKSCORE" = "red",
  "NABAKKEN" = "orange",
  "RIJPTIJD" = "green",
  "VERKLEURING.KOKEN" = "blue",
  "VLEESKLEUR.NA.KOKEN" = "purple"
)

circos.trackPoints(
  factors = qtl.table$Chrom,
  x = qtl.table$Position,
  y = rep(0, nrow(qtl.table)),
  cex = 1.2,
  pch = model.shape[qtl.table$Model],
  col = trait.colors[qtl.table$Trait]
)

circos.trackText(
  factors = qtl.table$Chrom,
  x = qtl.table$Position,
  y = rep(-1.3, nrow(qtl.table)),
  labels = qtl.table$Marker,
  facing = "reverse.clockwise",
  cex = 0.7
)

legend.traits <- Legend(at = c("NABAKKEN","RIJPTIJD","VERKLEURING.KOKEN"),
       title = "Trait", legend_gp = gpar(fill=c("orange", "green", "blue")))

legend.models <- Legend(at = c("additive", "1-dom-alt", "1-dom-ref"),
                        type = "points", pch = c(15,16,17), title = "Model")

text(0, 0, "Signifikante SNPs", cex = 1.5) # De "title"
text(0, -0.06, "p >= 0.05", cex = 0.8)

legends <- packLegend(legend.models, legend.traits)
draw(legends, x = unit(5, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))

circos.clear()
