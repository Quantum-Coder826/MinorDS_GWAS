library(tidyverse)
library(GWASpoly)

if (!exists("data.m.eff")) {data.m.eff <- readRDS("./outputs/data_m_eff.rds")}
traits <- c("RIJPTIJD","KOOKSCORE","VLEESKLEUR.NA.KOKEN","VERKLEURING.KOKEN","NABAKKEN")

snp <- "SNP24920"
trait <- "NABAKKEN"

qtl <- get.QTL(data.m.eff, trait = trait)

df <- data.frame(
  ID = rownames(data@geno),
  dosage = data@geno[, snp]
)

# Fenotype toevoegen
df <- merge(df,
            data@pheno[, c("ID", trait)],
            by = "ID")

ggplot(df, aes(x = factor(dosage), y = .data[[trait]])) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +
  labs(
    x = "Allele dosage",
    y = trait,
    title = marker
  ) +
  theme_bw()
