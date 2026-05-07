#test run
genofile <- system.file("extdata", "new_potato_geno.csv", package = "GWASpoly")
phenofile <- system.file("extdata", "new_potato_pheno.csv", package = "GWASpoly")

library(GWASpoly)
data <- read.GWASpoly(ploidy=4, pheno.file=phenofile, geno.file=genofile,
                      format="numeric", n.traits=1, delim=",")

data.loco <- set.K(data,LOCO=TRUE,n.core=2)
data.original <- set.K(data,LOCO=FALSE,n.core=2)

N <- 957 #Population size
params <- set.params(geno.freq = 1 - 5/N, fixed = "env", fixed.type = "factor")

data.loco.scan <- GWASpoly(data=data.loco,models=c("additive","1-dom"),
                           traits=c("vine.maturity"),params=params,n.core=2)

data.original.scan <- GWASpoly(data.original,models=c("additive","1-dom"),
                               traits=c("vine.maturity"),params=params,n.core=2)

library(ggplot2)
qq.plot(data.original.scan,trait="vine.maturity") + ggtitle(label="Original")
