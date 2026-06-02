library(tidyverse)
library(circlize)
library(ComplexHeatmap)

qtl.table <- read_csv("outputs/qtl_table.csv")
qtl.ans <- read_rds("./outputs/fit_ans.rds")
chrom_lenghts <- read_csv("data/potato_chrom_lengths_DMv6.csv")

###################################
### generate genome coordinates ###
###################################
qtl.table$coord <- paste0(qtl.table$Chrom, ":", qtl.table$Position, "..", qtl.table$Position)
View(qtl.table)

qtl.table <- qtl.table %>%
  arrange(desc(Effect)) %>%
  top_n(n = 10, wt = Effect)

### Genereer alle trait signifikante ###
trait.no = 6
fit.ans[[trait.no]] %>%
  filter(Trait == traits[trait.no]) %>%
  arrange(desc(R2)) %>%
  knitr::kable(.,digits = 3)


circos.genomicInitialize(chrom_lenghts)

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
  y = rep(-1.5, nrow(qtl.table)),
  labels = qtl.table$Marker,
  facing = "clockwise",
  cex = 0.7
)

legend.traits <- Legend(at = c("KOOKSCORE","NABAKKEN","RIJPTIJD","VERKLEURING.KOKEN","VLEESKLEUR.NA.KOKEN"),
       title = "Trait", legend_gp = gpar(fill=c("red", "orange", "green", "blue", "purple")))

legend.models <- Legend(at = c("additive", "1-dom-alt", "1-dom-ref"),
                        type = "points", pch = 15:17, title = "Model")

#text(0, 0, "Trait kookscore", cex = 1.5) # De "title"


legends <- packLegend(legend.models, legend.traits)
draw(legends, x = unit(5, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))

circos.clear()

