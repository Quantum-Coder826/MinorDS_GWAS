library(GWASpoly)
data <- read.GWASpoly(ploidy=4, pheno.file=phenofile, geno.file=genofile,
                      format="numeric", n.traits=1, delim=",")

data.loco <- set.K(data,LOCO=TRUE,n.core=3)
data.original <- set.K(data,LOCO=FALSE,n.core=3)

N <- 957 #Population size
params <- set.params(geno.freq = 1 - 5/N, fixed = "env", fixed.type = "factor")

data.loco.scan <- GWASpoly(data=data.loco,models=c("additive","1-dom"),
                           traits=c("Negatief"),params=params,n.core=3)

data.original.scan <- GWASpoly(data.original,models=c("additive","1-dom"),
                               traits=c("vine.maturity"),params=params,n.core=3)

library(ggplot2)

# Ik zet mijn wd hadmatig
#setwd("C:/Users/Hp/OneDrive - NHL Stenden/minor data science/gwasie")
NCORE = 16 #Hoeveel cores ik wil gebuiken, zet ik graag in een var om het snel te veranderen.

# In bestandspaden betekend `./` deze folder als in de wd van R.
geno <- "./data/potato_genos.csv"
pheno <- "./data/potato_phenos.csv"

ld[is.na(ld)] <- 0

Heatmap(ld, name = "r2")

#for real#########################################################################
setwd("C:/Users/Hp/OneDrive - NHL Stenden/minor data science/gwasie")
library(readxl)

geno <- "C:/Users/Hp/OneDrive - NHL Stenden/minor data science/gwasie/potato_genos.csv"
pheno <- "C:/Users/Hp/OneDrive - NHL Stenden/minor data science/gwasie/potato_phenos.csv"

library(GWASpoly)
data <- read.GWASpoly(ploidy=4, pheno.file=pheno, geno.file=geno,
                      format="numeric", n.traits=8, delim=";")

data.loco <- set.K(data,LOCO=TRUE,n.core=3)
data.original <- set.K(data,LOCO=FALSE,n.core=3)

N <- 260 #Population size
params <- set.params(geno.freq = 0.87, fixed= "Negatief", fixed.type = "numeric") # haal fixed weg

traits <- c("RIJPTIJD","KIEMRUST","KOOKSCORE","VLEESKLEUR.NA.KOKEN","VERKLEURING.KOKEN","NABAKKEN")
data.loco.scan <- GWASpoly(data=data.loco,models=c("additive","1-dom"), 
                           traits= traits,params=params,n.core=3)

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


library(ggplot2)
#alles behalve kooktype
qq.plot(data.original.scan,trait="VLEESKLEUR.NA.KOKEN") + ggtitle(label="VLEESKLEUR.NA.KOKEN")

#kooktype
install.packages("tidyverse")
library(tidyverse)

qq <- function(data, trait, models = NULL){
  stopifnot(inherits(data, "GWASpoly.fitted"))
  stopifnot(is.element(trait, names(data@scores)))
  
  all.models <- colnames(data@scores[[trait]])
  
  if (is.null(models)) {
    models <- all.models
  } else {
    dom.models <- models[grep("dom", models, fixed = TRUE)]
    models <- setdiff(models, dom.models)
    
    if (length(dom.models) > 0) {
      dom.models <- unlist(lapply(as.list(dom.models), 
                                  function(x) {
                                    all.models[grep(x, all.models, fixed = TRUE)]
                                  }))
    }
    
    models <- union(models, dom.models)
    stopifnot(all(is.element(models, all.models)))
  }
  
  scores <- as.data.frame(data@scores[[trait]][, models])
  colnames(scores) <- models
  scores$Chrom <- data@map$Chrom
  
  tmp <- pivot_longer(
    data = scores,
    cols = 1:length(models),
    names_to = "model",
    values_to = "y",
    values_drop_na = TRUE
  )
  
  tmp <- as.data.frame(tmp)
  tmp$model <- factor(tmp$model, levels = models, ordered = TRUE)
  
  tmp <- tmp[order(tmp$model, tmp$Chrom, tmp$y,
                   decreasing = c(FALSE, FALSE, TRUE)), ]
  
  tmp2 <- tapply(tmp$y, list(tmp$Chrom, tmp$model), function(x) {
    n <- length(x)
    -log10(ppoints(n))
  })
  
  tmp$x <- unlist(tmp2)
  
  return(tmp)
}

rijptijd <- qq(data.loco.scan, trait="RIJPTIJD") %>%
  filter(Chrom %in% c("chr03", "chr04", "chr05", "chr06")) %>%
  ggplot(aes(x, y, colour = model)) + 
  facet_wrap(~Chrom ) + 
  geom_point() + 
  theme_bw() + 
  xlab(expression(paste("Expected -log"[10], "(p)", sep = ""))) + 
  ylab(expression(paste("Observed -log"[10], "(p)", sep = ""))) + 
  scale_colour_brewer(palette = "PRGn") + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  theme(text = element_text(size = 15))  + ggtitle(label="RIJPTIJD")

rijptijd

#combine
library(patchwork)
library(gridGraphics)
install.packages("gridExtra")
library(gridExtra)

#maak de combo plot
plot1 <- (rijptijd / vleeskleurna / nabak) + 
  plot_annotation(title = "alles") &
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 6.3))

#zie plot
plot1
rijptijd
grid.arrange(rijptijd, vleeskleurna, nabak, ncol = 3)
grid.arrange(kiemrust, verkleurkoken, kookscore, ncol=3)

#m.eff gaf niks, te streng?
data2 <- set.threshold(data.loco.scan,method="Bonferroni",level=0.05)

p <- manhattan.plot(data2,traits="RIJPTIJD")
p + theme(axis.text.x = element_text(angle=90,vjust=0.5))
manhattan.plot(data2,traits="KIEMRUST",chrom="chr03")

#wont work
p <- LD.plot(data2, max.loci=1000)
p + xlim(0,60) 

#haal de diploweg
qtl <- get.QTL(data=data2,traits="RIJPTIJD",models=c("additive","1-dom","general"),bp.window=5e6)
knitr::kable(qtl)

fit.ans <- fit.QTL(data=data2,trait="RIJPTIJD",
                   qtl=qtl[,c("Marker","Model")])
knitr::kable(fit.ans,digits=3)


saveRDS(data.loco.scan, "./outputs/data_loco_scan.rds")

str(data)
#circel lmao
install.packages("circlize")


library(circlize)

# --------------------------
# 1. Definieer chromosomen
# --------------------------
# Stel: 3 chromosomen (bijv. aardappel)
chrom_lengths <- data.frame(
  chr = c("chr3", "chr4", "chr5", "chr6"),
  start = c(0, 0, 0, 0),
  end = c(60e6, 69e6, 55e6, 59e6)  # lengtes in bp (denkbeeldig)
)

# Initialiseer cirkel
circos.initialize(
  factors = chrom_lengths$chr,
  xlim = chrom_lengths[, c("start", "end")]
)

# --------------------------
# 2. Plot chromosoom-banden
# --------------------------

circos.trackPlotRegion(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    
    # label boven de chromosoomring
    circos.text(
      x = CELL_META$xcenter,
      y = 1.3,
      labels = chr,
      cex = 0.8
    )
  },
  bg.col = c("lightgrey"),   # 🔥 dit maakt de chromosomen zichtbaar
  bg.border = "black"
)
``


# --------------------------
# 3. Definieer SNPs
# --------------------------
# 3 SNPs met 3 eigenschappen (kleuren)
snps <- data.frame(
  chr = c("chr3", "chr4", "chr5", "chr6", "chr4"),
  pos = c(27758989, 1669280, 4809697, 52223818, 6422158),
  trait = c("RIJPTIJD", "RIJPTIJD", "RIJPTIJD", "RIJPTIJD", "KOOKSCORE")
)

# Kleur mapping
trait_colors <- c(
  "RIJPTIJD" = "red",
  "KOOKSCORE" = "blue",
  "RIJPTIJD" = "green"
)

# --------------------------
# 4. Voeg SNP punten toe
# --------------------------
circos.trackPoints(
  factors = snps$chr,
  x = snps$pos,
  y = rep(0.5, nrow(snps)),
  col = trait_colors[snps$trait],
  pch = 16,
  cex = 1.2
)

# --------------------------
# 5. Voeg labels toe
# --------------------------
circos.trackText(
  factors = snps$chr,
  x = snps$pos,
  y = rep(0.8, nrow(snps)),
  labels = snps$trait,
  col = trait_colors[snps$trait],
  cex = 0.7
)

# --------------------------
# 6. Reset plot
# --------------------------
circos.clear()
``

#dosage=============================================================

qtl <- get.QTL(data2, trait = "NABAKKEN")
qtl

library(ggplot2)

marker <- "SNP24920"

df <- data.frame(
  ID = rownames(data@geno),
  dosage = data@geno[, marker]
)

# Fenotype toevoegen
df <- merge(df,
            data@pheno[, c("ID", "NABAKKEN")],
            by = "ID")

ggplot(df, aes(x = factor(dosage), y = NABAKKEN)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +
  labs(
    x = "Allele dosage",
    y = "NABAKKEN",
    title = marker
  ) +
  theme_bw()
