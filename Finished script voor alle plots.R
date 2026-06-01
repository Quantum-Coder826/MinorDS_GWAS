library(GWASpoly)
library(tidyverse)

Saveplot <- function(x, y) {
  name <- paste0("Plots_Brent/Manhattan/Single_Chrom/ManhattanPlot_", x, "_", y, ".png")
  ggsave(
    filename = name,
    plot = manhattan.plot(data2, traits=x, chrom=y) + 
      theme(axis.text.x = element_text(angle=270,vjust=0.5)),
    scale = 1,
    limitsize = TRUE,
    create.dir = TRUE,
  )
}

redo = TRUE
if (redo == FALSE) {
  data.loco.scan <- readRDS("./data_Brent/data_loco_scan.rds")
  data.original.scan <- readRDS("./data_Brent/data_original_scan.rds")
} 
if(redo == TRUE){
Pheno <- "data_Brent/potato_phenos.csv"
Geno <- "data_Brent/potato_genos.csv"


data <- read.GWASpoly(
  ploidy     = 4,
  pheno.file = Pheno,
  geno.file  = Geno,
  format     = "numeric",
  n.traits   = 7, 
  delim      = ";"
)
######################################################
data.loco <- set.K(data,LOCO=TRUE,n.core=2)
data.original <- set.K(data,LOCO=FALSE,n.core=2)
######################################################
N <- 260
traits <- c("RIJPTIJD","KOOKSCORE","VLEESKLEUR.NA.KOKEN","VERKLEURING.KOKEN","NABAKKEN")
params <- set.params(geno.freq = 1 - 5/N, fixed="Negatief", fixed.type = "numeric")
######################################################
data.loco.scan <- GWASpoly(data=data.loco,
                           models=c("additive","1-dom"),
                           traits=traits,
                           params=params,
                           n.core=2)


data.original.scan <- GWASpoly(data.original,
                               models=c("additive","1-dom"),
                               traits=traits,
                               params=params,
                               n.core=2)

saveRDS(data.loco.scan, "./data_Brent/data_loco_scan.rds")
saveRDS(data.original.scan, "./data_Brent/data_original_scan.rds")
#Save hoeft niet elke run opnieuw
}
######################################################
for(tr in traits) {
  name <- paste0("Plots_Brent/QQ//QQplot_", tr, ".png")
  
  ggsave(
    filename = name,
    plot = qq.plot(data.original.scan, trait=tr) + ggtitle(label=tr),
    scale = 1,
    limitsize = TRUE,
    create.dir = TRUE,
  )
}

######################################################
#data2 <- set.threshold(data.loco.scan,method="Bonferroni",level=0.05) #niet nodig
data2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05)
######################################################
ggsave(
  filename = paste0("Plots_Brent/Manhattan/ManhattanPlot_All.png"),
  plot = manhattan.plot(data2,traits=traits) + 
    theme(axis.text.x = element_text(angle=270,vjust=0.5)),
  scale = 1,
  limitsize = TRUE,
  create.dir = TRUE,
)

for(tr in traits) {
  name <- paste0("Plots_Brent/Manhattan/ManhattanPlot_", tr, ".png")
  ggsave(
    filename = name,
    plot = manhattan.plot(data2, traits=tr) + 
      theme(axis.text.x = element_text(angle=270,vjust=0.5)),
    scale = 1,
    limitsize = TRUE,
    create.dir = TRUE,
  )
    if (tr == "RIJPTIJD"){
      Saveplot(tr, "chr03")
      Saveplot(tr, "chr04")
      Saveplot(tr, "chr05")
    }
    if (tr == "KOOKSCORE"){
      Saveplot(tr, "chr04")
    }
    if (tr == "VLEESKLEUR.NA.KOKEN"){
      Saveplot(tr, "chr02")
      Saveplot(tr, "chr03")
    }
    if (tr == "VERKLEURING.KOKEN"){
      Saveplot(tr, "chr01")
      Saveplot(tr, "chr03")
      Saveplot(tr, "chr06")
    }
    if (tr == "NABAKKEN"){
      Saveplot(tr, "chr05")
      Saveplot(tr, "chr08")
    }
}
######################################################
ggsave(
  filename = paste0("Plots_Brent/LD/LD_Plot_All.png"),
  plot = LD.plot(data2, max.loci=1000) + 
    xlim(0,30),
  scale = 1,
  limitsize = TRUE,
  create.dir = TRUE,
)
####################################################
qtl <- get.QTL(data=data2,
               traits=traits,
               models="additive",
               bp.window=5e6)

knitr::kable(qtl) 
#####################################################
for(tr in traits) {
  x <- fit.QTL(data=data2,
               trait=tr,
               qtl=qtl[,c("Marker","Model")])
  knitr::kable(x, digits=3)
}