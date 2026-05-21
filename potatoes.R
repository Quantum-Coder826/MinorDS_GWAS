library(GWASpoly)
library(ggplot2)

# Ik zet mijn wd hadmatig
#setwd("C:/Users/Hp/OneDrive - NHL Stenden/minor data science/gwasie")
NCORE = 16 #Hoeveel cores ik wil gebuiken, zet ik graag in een var om het snel te veranderen.

# In bestandspaden betekend `./` deze folder als in de wd van R.
geno <- "./data/potato_genos.csv"
pheno <- "./data/potato_phenos.csv"

data <- read.GWASpoly(ploidy=4, pheno.file=pheno, geno.file=geno,
                      format="numeric", n.traits=7, delim=";")

data.loco <- set.K(data,LOCO=TRUE,n.core=NCORE)

N <- 260 #Population size
params <- set.params(geno.freq = 1 - 5/N, fixed= "Negatief", fixed.type = "numeric")


traits <- c("RIJPTIJD","KIEMRUST","KOOKSCORE","VLEESKLEUR.NA.KOKEN","VERKLEURING.KOKEN","NABAKKEN")
data.loco.scan <- GWASpoly(data=data.loco,models=c("additive","1-dom"),
                           traits=traits,params=params,n.core=NCORE)

saveRDS(data.loco.scan, "./outputs/data_loco_scan.rds")