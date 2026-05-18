#test run
genofile <- system.file("extdata", "new_potato_geno.csv", package = "GWASpoly")
phenofile <- system.file("extdata", "new_potato_pheno.csv", package = "GWASpoly")

library(GWASpoly)
data <- read.GWASpoly(ploidy=4, pheno.file=phenofile, geno.file=genofile,
                      format="numeric", n.traits=1, delim=",")

data.loco <- set.K(data,LOCO=TRUE,n.core=3)
data.original <- set.K(data,LOCO=FALSE,n.core=3)

N <- 957 #Population size
params <- set.params(geno.freq = 1 - 5/N, fixed = "env", fixed.type = "factor")

data.loco.scan <- GWASpoly(data=data.loco,models=c("additive","1-dom"),
                           traits=c("vine.maturity"),params=params,n.core=3)

data.original.scan <- GWASpoly(data.original,models=c("additive","1-dom"),
                               traits=c("vine.maturity"),params=params,n.core=3)

library(ggplot2)
qq.plot(data.original.scan,trait="vine.maturity") + ggtitle(label="Original")

data2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05)

p <- manhattan.plot(data2,traits="vine.maturity")
p + theme(axis.text.x = element_text(angle=90,vjust=0.5))
manhattan.plot(data2,traits="vine.maturity",chrom="chr05")


p <- LD.plot(data2, max.loci=1000)
p + xlim(0,30) 

qtl <- get.QTL(data=data2,traits="vine.maturity",models="additive",bp.window=5e6)
knitr::kable(qtl)

fit.ans <- fit.QTL(data=data2,trait="vine.maturity",
                   qtl=qtl[,c("Marker","Model")],
                   fixed=data.frame(Effect="env",Type="factor"))
knitr::kable(fit.ans,digits=3)

install.packages("ComplexHeatmap")
library(ComplexHeatmap)

geno <- data@geno[c(1:10), c(1:10)]

#je wil R2. cor geeft alleen richting en sterkte, r2 geeft koppeling tussen 2 snp
ld <- cor(geno, use = "pairwise.complete.obs")^2

ld[is.na(ld)] <- 0

Heatmap(ld, name = "r2")

#for real
setwd("C:/Users/Hp/OneDrive - NHL Stenden/minor data science/gwasie")
library(readxl)

geno <- "C:/Users/Hp/OneDrive - NHL Stenden/minor data science/gwasie/potato_genos.csv"
pheno <- "C:/Users/Hp/OneDrive - NHL Stenden/minor data science/gwasie/potato_phenos.csv"

library(GWASpoly)
data <- read.GWASpoly(ploidy=4, pheno.file=pheno, geno.file=geno,
                      format="numeric", n.traits=9, delim=";")

data.loco <- set.K(data,LOCO=TRUE,n.core=3)
data.original <- set.K(data,LOCO=FALSE,n.core=3)

N <- 260 #Population size
params <- set.params(geno.freq = 1 - 5/N, fixed= "Negatief", fixed.type = "numeric")










