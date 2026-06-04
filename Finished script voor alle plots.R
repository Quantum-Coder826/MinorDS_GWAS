library(GWASpoly)
library(tidyverse)
library(gridExtra)
All_man_tr_chr <<- NULL
All_man_tr_chr <- character()

Saveplot <- function(x, y) {
  name <- paste0("Plots_Brent/Manhattan/Single_Chrom/ManhattanPlot_", x, "_", y, ".png")
  name2 <- paste0("ManhattanPlot_", x, "_", y)
  name3 <- paste0(x, "_", y)
  
  Man_tr_chr <- manhattan.plot(data2, traits=x, chrom=y) + 
    theme(axis.text.x = element_text(angle=270,vjust=0.5)) +
    ggtitle(label = name2)
  
  
  All_man_tr_chr[[name3]] <<- Man_tr_chr
  
  
  ggsave(
    filename = name,
    plot = Man_tr_chr,
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
      #Saveplot(tr, "chr02")
      Saveplot(tr, "chr03")
    }
    if (tr == "VERKLEURING.KOKEN"){
      Saveplot(tr, "chr01")
      Saveplot(tr, "chr03")
      #Saveplot(tr, "chr06")
    }
    if (tr == "NABAKKEN"){
     # Saveplot(tr, "chr05")
      #Saveplot(tr, "chr08")
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

#####################################################
AllelDosage <- function(marker, traitname, data, chr){
  # marker = bv. "SNP24920"
  #trait = "trait"
GenoMarker <- data.frame(
  ID = rownames(data@geno),
  dosage = data@geno[, marker]
)

 MarkerData<- merge(GenoMarker,
            data@pheno[, c("ID", traitname)],
            by = "ID")

 ggsave(
   filename = paste0("Plots_Brent/Allele_Dosage/AD_", marker,"_", traitname,"_", chr, ".png"),
   plot = ggplot(MarkerData, aes(x = factor(dosage), y = .data[[traitname]])) +
     geom_boxplot() +
     geom_jitter(width = 0.15, alpha = 0.5) +
     labs(
       x = "Allele dosage",
       y = traitname,
       title = marker
     ) +
     theme_bw(),
   scale = 1,
   limitsize = TRUE,
   create.dir = TRUE,
 )
}
#################################################
AllelDosage("SNP20248", "RIJPTIJD", data, "chr04")
AllelDosage("SNP17347", "RIJPTIJD", data, "chr03")
AllelDosage("SNP23778", "RIJPTIJD", data, "chr05")
AllelDosage("SNP24824", "RIJPTIJD", data, "chr05")

AllelDosage("SNP20811", "KOOKSCORE", data, "chr04")
AllelDosage("SNP21172", "KOOKSCORE", data, "chr04")
AllelDosage("SNP21737", "KOOKSCORE", data, "chr04")

AllelDosage("SNP15706", "VLEESKLEUR.NA.KOKEN", data, "chr03")
AllelDosage("SNP16489", "VLEESKLEUR.NA.KOKEN", data, "chr03")
AllelDosage("SNP18001", "VLEESKLEUR.NA.KOKEN", data, "chr03")

AllelDosage("SNP1655", "VERKLEURING.KOKEN", data, "chr01")
AllelDosage("SNP16716", "VERKLEURING.KOKEN", data, "chr03")
AllelDosage("SNP18001", "VERKLEURING.KOKEN", data, "chr03")
AllelDosage("SNP19842", "VERKLEURING.KOKEN", data, "chr03")

###################################################


