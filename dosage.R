library(tidyverse)
library(GWASpoly)
library(gridExtra)

if (!exists("data.m.eff")) {data.m.eff <- readRDS("./outputs/data_m_eff.rds")}
traits <- c("RIJPTIJD","KOOKSCORE","VLEESKLEUR.NA.KOKEN","VERKLEURING.KOKEN","NABAKKEN")

dosage.plot <- function(snp, trait) {
  qtl <- get.QTL(data.m.eff, trait = trait)
  
  df <- data.frame(
    ID = rownames(data@geno),
    dosage = data@geno[, snp]
  )
  
  # Fenotype toevoegen
  df <- merge(df,
              data@pheno[, c("ID", trait)],
              by = "ID")
  
  plot <- ggplot(df, aes(x = factor(dosage), y = .data[[trait]])) +
    geom_boxplot() +
    geom_jitter(width = 0.15, alpha = 0.5) +
    labs(
      x = "Allele dosage",
      y = "Score",
      title = paste0(snp,"\n", trait)
    ) +
    theme_bw()
  
  return(plot)
}

plots_gorb <- list()
#plots_gorb <- append(plots_gorb, dosage.plot("", ""))
plots_gorb <- append(plots_gorb, dosage.plot("SNP23778", "RIJPTIJD")) # dosage 2 & 3
#plots_gorb <- append(plots_gorb, dosage.plot("SNP34378", "RIJPTIJD"))

plots_gorb <- append(plots_gorb, dosage.plot("SNP23047", "VLEESKLEUR.NA.KOKEN"))# dosage 3
plots_gorb <- append(plots_gorb, dosage.plot("SNP16489", "VLEESKLEUR.NA.KOKEN"))# dosage 1 & 2

#lots_gorb <- append(plots_gorb, dosage.plot("SNP28271", "VERKLEURING.KOKEN"))
plots_gorb <- append(plots_gorb, dosage.plot("SNP16782", "VERKLEURING.KOKEN")) # dosage 2 & 3
#plots_gorb <- append(plots_gorb, dosage.plot("SNP1655", "VERKLEURING.KOKEN")) # All dosage Honarable mention

plots_gorb <- append(plots_gorb, dosage.plot("SNP33303", "NABAKKEN")) # dosage 2
plots_gorb <- append(plots_gorb, dosage.plot("SNP24920", "NABAKKEN")) # all exept dosage 4
#plots_gorb <- append(plots_gorb, dosage.plot("SNP22740", "NABAKKEN"))

grid.arrange(grobs = plots_gorb, ncol = 3)
