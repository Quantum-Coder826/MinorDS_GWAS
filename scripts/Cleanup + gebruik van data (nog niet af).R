library(readxl)
library(tidyverse)
library(GWASpoly)
library(ggplot2)

# ================================================
# STEP 1: Find the row where geno starts
# ================================================
raw <- read_excel("data/DatasetGWAS.xlsx", sheet = 2, col_names = FALSE)

# Find which row contains "Chromosome" - this is the geno header row
geno_header_row <- which(raw[[1]] == "Chromosome")
cat("Geno starts at row:", geno_header_row, "\n")

# How many pheno rows are there (excluding the header row)
pheno_data_rows <- geno_header_row - 2  # -1 for header, -1 for blank/gap row

# ================================================
# STEP 2: Read pheno block
# ================================================
pheno_raw <- read_excel("data/DatasetGWAS.xlsx", 
                        sheet = 2,
                        n_max = pheno_data_rows,
                        col_names = TRUE)

# Drop Allele.A, Allele.B
pheno_raw <- pheno_raw[, !names(pheno_raw) %in% c("Allele.A", "Allele.B", "...2", "...3")]
names(pheno_raw)[1] <- "Trait"

# Remove uninformative/Useless traits
pheno_raw <- pheno_raw[!pheno_raw$Trait %in% c("KIEMRUST"), ]
is_letter_row <- apply(pheno_raw[, -1], 1, 
                       function(x) { any( grepl("[A-Za-z]", na.omit(x)) ) }
                       )
pheno_raw <- pheno_raw[!is_letter_row, ]

# Convert all trait values to numeric (Just in case)
pheno_raw[, -1] <- lapply(pheno_raw[, -1], as.numeric)

# Transpose: individuals become rows, traits become columns
pheno_t           <- as.data.frame(t(pheno_raw[, -1]))
colnames(pheno_t) <- pheno_raw$Trait
pheno_t$Name      <- rownames(pheno_t)

# Name column first, then traits
pheno_final <- pheno_t[, c("Name", pheno_raw$Trait)]

# ================================================
# STEP 3: Read geno block
# ================================================
geno_raw <- read_excel("data/DatasetGWAS.xlsx", 
                       sheet = 2,
                       skip = geno_header_row - 1,
                       col_names = TRUE)

# Grab individual names from the pheno header row in the raw sheet
ind_names <- as.character(raw[1, -(1:3)])  # skip first 3 cols (Phenotype, Allele.A, Allele.B)
ind_names <- ind_names[!is.na(ind_names) & ind_names != "NA"]

cat("Individual names found:", length(ind_names), "\n")
cat("First few:", head(ind_names), "\n")

# Remove annotation columns
geno <- geno_raw[, !names(geno_raw) %in% c("Ref", "Alt", "Allele.A", "Allele.B")]

# Rename marker metadata columns
geno <- geno %>% rename(
  Marker     = `SNP number`,
  Chromosome = Chromosome,
  Position   = Position
)
# Assign the individual names to the empty columns (columns 4 onward)
colnames(geno)[4:ncol(geno)] <- ind_names

# Format chromosome as chr01, chr02 etc. & remove Unknown
geno <- geno %>% 
  filter(Chromosome %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
print(geno)
names(geno)[duplicated(names(geno))]
geno$Chromosome <- paste0("chr", sprintf("%02d", as.integer(geno$Chromosome)))

# Enforce column order: Marker, Chromosome, Position, then individuals
geno <- geno[, c("Marker","Chromosome", "Position",
                 setdiff(names(geno), c("Marker", "Chromosome", "Position")))]

# ================================================
# STEP 4: Verify sample names match
# ================================================
pheno_names <- pheno_final$Name
geno_names  <- colnames(geno)[4:ncol(geno)]

missing_from_geno  <- setdiff(pheno_names, geno_names)
missing_from_pheno <- setdiff(geno_names, pheno_names)

cat("In pheno but not geno:", length(missing_from_geno),  "\n")
cat("In geno but not pheno:", length(missing_from_pheno), "\n")

if(length(missing_from_geno) > 0)  print(missing_from_geno)
if(length(missing_from_pheno) > 0) print(missing_from_pheno)

# ================================================
# STEP 5: Clean geno markers
# ================================================
# Remove markers with >10% missing data
missing_rate <- rowMeans(is.na(geno[, 4:ncol(geno)]))
geno <- geno[missing_rate < 0.10, ]

# Remove monomorphic markers (no variation = no information)
is_mono <- apply(geno[, 4:ncol(geno)], 1,
                 function(x) length(unique(na.omit(x))) == 1)
geno <- geno[!is_mono, ]

cat("Markers remaining after cleaning:", nrow(geno), "\n")
cat("Individuals in pheno:",            nrow(pheno_final), "\n")

# ================================================
# STEP 6: Save CSVs
# ================================================
write.csv(pheno_final, "data/pheno_clean.csv", row.names = FALSE)
write.csv(geno,        "data/geno_clean.csv",  row.names = FALSE)

# ================================================
# STEP 7: Load into GWASpoly
# ================================================
data <- read.GWASpoly(
  ploidy     = 4,
  pheno.file = "data/pheno_clean.csv",
  geno.file  = "data/geno_clean.csv",
  format     = "numeric",
  n.traits   = ncol(pheno_final) - 1,  # all columns except Name
  delim      = ","
)
?read.GWASpoly
######################################################
data.loco <- set.K(data,LOCO=TRUE,n.core=2)
data.original <- set.K(data,LOCO=FALSE,n.core=2)
######################################################
N <- nrow(data@pheno)
Traits <- colnames(data@pheno)[-1]
params <- set.params(geno.freq = 1 - 5/N)
######################################################
data.loco.scan <- GWASpoly(data=data.loco,
                           models=c("additive","1-dom"),
                           traits=Traits,
                           params=params,
                           n.core=2)

data.original.scan <- GWASpoly(data.original,
                               models=c("additive","1-dom"),
                               traits=Traits,
                               params=params,
                               n.core=2)
######################################################
print(data@map$Chrom)
for(tr in Traits) {
  name <- paste0("Plots/QQplot_", tr, ".png")

  ggsave(
    filename = name,
    plot = qq.plot(data.original.scan, trait=tr) + ggtitle(label=tr),
    scale = 1,
    limitsize = TRUE,
    create.dir = TRUE,
  )
}
#print(qq.plot(data.original.scan, trait=tr) + ggtitle(label=tr))
######################################################
data2 <- set.threshold(data.loco.scan,method="Bonferroni",level=0.05)
#data2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05)

#M.eff is beter maar werkt niet?
######################################################
p <- manhattan.plot(data2,traits=Traits)
p + theme(axis.text.x = element_text(angle=270,vjust=0.5))

#te ingezoomd?

#manhattan.plot(data2, traits=Traits, chrom="chr04")

######################################################
p <- LD.plot(data2, max.loci=1000)
p + xlim(0,30) 
####################################################
qtl <- get.QTL(data=data2,
               traits=Traits,
               models="additive",
               bp.window=5e6)
knitr::kable(qtl) 

#nog niet gefixt

#####################################################
fit.ans <- fit.QTL(data=data2,
                   traits=Traits,
                   qtl=qtl[,c("Marker","Model")])
knitr::kable(fit.ans,digits=3)

#nog niet gefixt


#maak van alle plots files, het laad veeeeel te sloom
