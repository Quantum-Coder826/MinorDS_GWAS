library(GWASpoly)
library(tidyverse)
source("./qq.R")

# Ik zet mijn wd hadmatig
#setwd("C:/Users/Hp/OneDrive - NHL Stenden/minor data science/gwasie")
NCORE = 16 #Hoeveel cores ik wil gebuiken, zet ik graag in een var om het snel te veranderen.

# In bestandspaden betekend `./` deze folder als in de wd van R.
geno <- "./data/potato_genos.csv"
pheno <- "./data/potato_phenos.csv"

data <- read.GWASpoly(ploidy=4, pheno.file=pheno, geno.file=geno,
                      format="numeric", n.traits=7, delim=";")

############################
### Analyseer de samples ###
############################
data.loco <- set.K(data,LOCO=TRUE,n.core=NCORE)

N <- 260 #Population size
params <- set.params(geno.freq = 1 - 5/N, fixed="Negatief", fixed.type = "numeric")

traits <- c("RIJPTIJD","KOOKSCORE","VLEESKLEUR.NA.KOKEN","VERKLEURING.KOKEN","NABAKKEN")
data.loco.scan <- GWASpoly(data=data.loco,models=c("additive","1-dom"), 
                           traits=traits,params=params,n.core=NCORE)

saveRDS(data.loco.scan, "./outputs/data_loco_scan.rds") # Save hoeft niet elke run opnieuw

##########################
#### Visualiseren data ###
##########################
# Start hier met runnen waneer alle analyse onnodig is.
if (!exists("data.loco.scan")) {
  data.loco.scan <- readRDS("./outputs/data_loco_scan.rds")
}

# qqplots
for (trait in traits) {
  qq.plot(data.loco.scan, trait = trait)
  p + ggtitle(trait)
  ggsave(paste0("./Plots/QQ/", trait, ".png"))
  print(paste("Saved:", trait))
}

data.m.eff <- set.threshold(data.loco.scan,method="M.eff",level=0.05, n.core = NCORE)

# manhattan plots
for (trait in traits) {
  p <- manhattan.plot(data.m.eff,traits=trait)
  p + theme(axis.text.x = element_text(angle=90,vjust=0.5))
  p + ggtitle(trait)
  ggsave(paste0("./Plots/Manhattan/", trait, ".png"))
  print(paste("Saved:", trait))
}

p <- LD.plot(data.m.eff, max.loci=1000)
p + xlim(0,40) 
ggsave("./Plots/LDplot.png")

qtl <- get.QTL(data=data.m.eff,traits=traits,bp.window=15e6) #Uit DL kwam ~15mBp
write.csv(qtl, "./outputs/qtl_table.csv", sep = ",")
View(qtl)