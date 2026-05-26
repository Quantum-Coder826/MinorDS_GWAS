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
params <- set.params(geno.freq = 0.90, fixed="Negatief", fixed.type = "numeric")

traits <- c("RIJPTIJD","KIEMRUST","KOOKSCORE","VLEESKLEUR.NA.KOKEN","VERKLEURING.KOKEN","NABAKKEN")
data.loco.scan <- GWASpoly(data=data.loco,models=c("additive","1-dom"), 
                           traits=traits,params=params,n.core=NCORE)

saveRDS(data.loco.scan, "./outputs/data_loco_scan.rds") # Save hoeft niet elke run opnieuw

##########################
#### Visualiseren data ###
##########################
# Remember me?
if (!exists("data.loco.scan")) {
  data.loco.scan <- readRDS("./outputs/data_loco_scan.rds")
}

qq.plot(data.loco.scan, trait = "KOOKSCORE")

qq(data.loco.scan, trait="RIJPTIJD") %>%
  filter(Chrom == c("chr03", "chr04", "chr05")) %>%
    ggplot(aes(x, y, colour = model)) + 
    facet_wrap(~Chrom) + 
    geom_point() + 
    theme_bw() + 
    xlab(expression(paste("Expected -log"[10], "(p)", sep = ""))) + 
    ylab(expression(paste("Observed -log"[10], "(p)", sep = ""))) + 
    scale_colour_brewer(palette = "Set1") + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) + 
    theme(text = element_text(size = 15))

#m.eff gaf niks, te streng?
data2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05) # errors???

p <- manhattan.plot(data2,traits="RIJPTIJD")
p + theme(axis.text.x = element_text(angle=90,vjust=0.5))
manhattan.plot(data2,traits="RIJPTIJD",chrom="chr03") #NOTE: chromesomen zijn chrXX geformateerd

#wont work
p <- LD.plot(data2, max.loci=1000)
p + xlim(0,30) 

qtl <- get.QTL(data=data2,traits="KOOKTYPE",models="additive",bp.window=5e6)
knitr::kable(qtl)

fit.ans <- fit.QTL(data=data2,trait="KOOKTYPE_num",
                   qtl=qtl[,c("Marker","Model")])
knitr::kable(fit.ans,digits=3)
