library(GWASpoly)
library(ggplot2)

# Ik zet mijn wd hadmatig
#setwd("C:/Users/Hp/OneDrive - NHL Stenden/minor data science/gwasie")
NCORE = 16 #Hoeveel cores ik wil gebuiken, zet ik graag in een var om het snel te veranderen.

# In bestandspaden betekend `./` deze folder als in de wd van R.
geno <- "./data/potato_genos.csv"
pheno <- "./data/potato_phenos.csv"

ld[is.na(ld)] <- 0

Heatmap(ld, name = "r2")

geno <- "./data/potato_genos.csv"
pheno <- "./data/potato_phenos.csv"

data <- read.GWASpoly(ploidy=4, pheno.file=pheno, geno.file=geno,
                      format="numeric", n.traits=8, delim=";")

data.loco <- set.K(data,LOCO=TRUE,n.core=NCORE)

N <- 260 #Population size
params <- set.params(geno.freq = 0.90, fixed="Negatief", fixed.type = "numeric") # haal fixed weg

traits <- c("RIJPTIJD","KIEMRUST","KOOKSCORE","VLEESKLEUR.NA.KOKEN","VERKLEURING.KOKEN","NABAKKEN")
data.loco.scan <- GWASpoly(data=data.loco,models=c("additive","1-dom"), 
                           traits=traits,params=params,n.core=3)

data.original.scan <- GWASpoly(data.original,models=c("additive","1-dom"),
                               traits=traits,params=params,n.core=3)

#kooktype ========= werkt voor nu=======
data.loco@pheno$KOOKTYPE_num <- as.numeric(as.factor(trimws(data.loco@pheno$KOOKTYPE)))

data.original@pheno$KOOKTYPE <- trimws(data.original@pheno$KOOKTYPE)
data.original@pheno$KOOKTYPE_num <- as.numeric(as.factor(data.original@pheno$KOOKTYPE))

data.loco.scan2 <- GWASpoly(data=data.loco,models=c("additive","1-dom"), 
                           traits="KOOKTYPE_num",params=params,n.core=3)

data.original.scan2 <- GWASpoly(data.original,models=c("additive","1-dom"),
                               traits="KOOKTYPE_num",params=params,n.core=3)


#alles behalve kooktype
qq.plot(data.original.scan,trait="RIJPTIJD") + ggtitle(label="Original")

#kooktype
qq.plot(data.original.scan2,trait="KOOKTYPE_num") + ggtitle(label="kooktype")

#m.eff gaf niks, te streng?
data2 <- set.threshold(data.loco.scan2,method="FDR",level=0.05)

p <- manhattan.plot(data2,traits="RIJPTIJD")
p + theme(axis.text.x = element_text(angle=90,vjust=0.5))
manhattan.plot(data2,traits="RIJPTIJD",chrom="3")

#wont work
p <- LD.plot(data2, max.loci=1000)
p + xlim(0,30) 

qtl <- get.QTL(data=data2,traits="KOOKTYPE_num",models="additive",bp.window=5e6)
knitr::kable(qtl)

fit.ans <- fit.QTL(data=data2,trait="KOOKTYPE_num",
                   qtl=qtl[,c("Marker","Model")])
knitr::kable(fit.ans,digits=3)

saveRDS(data.loco.scan, "./outputs/data_loco_scan.rds")