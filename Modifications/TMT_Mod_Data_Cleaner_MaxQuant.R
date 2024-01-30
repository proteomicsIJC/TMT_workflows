### TMT Standard Workflow

### Libraries and WD----
library(dplyr)
library(tidyr)
library(limma)
library(sva)
library(janitor)
library(seqinr)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(reshape2)
library(mice)
library(gplots)
library(pheatmap)
library(openxlsx)
library(stringr)
library(rstudioapi)
library(data.table)
setwd(dirname(getActiveDocumentContext()$path))
#----

### Function deffinition----
source("./functions/general")
source("./functions/TMT_MaxQuant")
#----

### Get the data----
maxquant <- read.table("./phospho_data/Phospho (STY)Sites.txt", header = T, check.names = F, sep = "\t", quote = "", dec = ".")
maxquant2 <- maxquant
#----

maxquant_initalizer <- function(tmt,n_plex,maxquant_data){
  old_names <- colnames(maxquant_data)[grep(x = colnames(maxquant_data), pattern = "^Reporter intensity corrected")]
  reactives6 <- c("126","127n","128c","129n","130c","131")
  reactives10 <- c("126","127n","127c","128n","128c","129n","129c","130n","130c","131n")
  reactives11 <- c("126","127n","127c","128n","128c","129n","129c","130n","130c","131n","131c")
  reactives16 <- c("126","127n","127c","128n","128c","129n","129c","130n","130c","131n","131c","132n","132c","133n","133c","134n")
  if (tmt == "tmt6"){
    # tmt6
    reactives <- reactives6
  }
  
  if (tmt == "tmt10"){
    # tmt10
    reactives <- reactives10
  }
  
  if (tmt == "tmt11"){
    # tmt11
    reactives <- reactives11
  }
  
  if (tmt == "tmt16"){
    # tmt16
    reactives <- reactives16
  }
  
  
  number_of_plexes <- c(1:n_plex)
  plex <- c()
  for (i in 1:length(number_of_plexes)){
    plex[i] <- paste("plex",number_of_plexes[i],sep = "")
  }
  
  plex <- rep(plex, each = length(reactives))
  reactives <- rep(reactives, times = length(unique(plex)))
  
  new_names <- c()
  for (i in 1:length(plex)){
    new_names[i] <- paste("Intensity ",plex[i]," : ",toupper(tmt),"-",toupper(reactives[i]),sep = "")
  }
  
  name_changer <- data.frame(old = old_names,
                             new = new_names)
  
  setnames(maxquant_data, name_changer$old, name_changer$new)
  
  return(maxquant_data)
}

#maxquant <- maxquant_initalizer(maxquant_data = maxquant, tmt = "tmt16", n_plex = 3)

### Get the correct columns----
# Retrive only intensities, and protein annotation columns
intensities <- grep(pattern = "^Intensity plex[0-9] \\: TMT", colnames(maxquant))
old_names <- colnames(maxquant)[33:128]
reverse <- grep(pattern = "^Reverse", colnames(maxquant))
potential_cont <- grep(pattern = "^Potential contaminant", colnames(maxquant))
only_site <- grep(pattern = "^Only identified by site", colnames(maxquant))
maxquant <- maxquant[,c(c(1,2,6,7),
                        c(8,24,30,10,11,25,28),
                  c(33:128),
                  c(reverse,potential_cont,only_site))]
colnames(maxquant)
# Clear dataset colnames
maxquant <- clean_names(maxquant)
colnames(maxquant)

# Rename majority to protein group
colnames(maxquant)[1] <- "protein_group"
# Transform all zero intensity values to NA
maxquant <- zero_to_NA(patterns = "plex", dataset = maxquant)
#----

### Remove contaminants----
maxquant_clean <- proteinGroupsCleanner(ds = maxquant)
maxquant_clean <- subset(maxquant_clean, select = -c(reverse, only_identified_by_site))
#----

#### Deconvoluting multiplicity from columns to rows...
# Keep intensity and meta-data apart
intensityColumns <- grep("reporter_intensity_corrected", colnames(maxquant_clean))
intensityColumns <- colnames(maxquant_clean)[intensityColumns]
nonIntensityColumns <- grep("reporter_intensity_corrected", colnames(maxquant_clean), invert = T)
nonIntensityColumns <- colnames(maxquant_clean)[nonIntensityColumns]

# Repair the colnames a bit
for (i in 1:length(colnames(maxquant_clean))){
  if (colnames(maxquant_clean)[i] %in% intensityColumns){
    colnames(maxquant_clean)[i] <- sub(pattern = "_([^_]*)$", replacement = ".\\1",x = colnames(maxquant_clean[i]))
  }
}

maxquant_clean <- maxquant_clean %>%
  tidyr::pivot_longer(
    cols = -all_of(nonIntensityColumns),
    names_to = c("sample","multiplicity"),
    names_pattern = "(reporter_intensity_corrected_[0-9]{1,2}_plex[1-9])\\.([1-3])",
    )

maxquant_clean <- maxquant_clean %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = value
  )

# In multiplicity put the "___" (there's 3)
for (i in 1:length(maxquant_clean$multiplicity)){
  maxquant_clean$multiplicity[i] <- paste("___",maxquant_clean$multiplicity[i],sep = "")
}


### NAs removal
row_has_na <- apply(maxquant_clean[,13:44], 1, function(x){all(is.na(x))})
maxquant_clean <- maxquant_clean[!row_has_na,]

### Loc probability
maxquant_clean <- maxquant_clean %>% 
  filter(localization_prob >= 0.5)
  
### Get the old names again
for (i in 1:length(old_names)){
  old_names[i] <- gsub(x = old_names[i],pattern = "___[1-3]",replacement = "")
}

old_names <- unique(old_names)
new_names <- colnames(maxquant_clean)[startsWith(x = colnames(maxquant_clean), prefix = "reporter")]

####### CREUAT A POSTA !!!######################################################
name_changer <- data.frame(old = new_names,
                           new = old_names)

setnames(maxquant_clean, name_changer$old, name_changer$new)

### Create a protein group column
colnames(maxquant_clean)[1] <- "protein"

maxquant_clean$protein_group <- "Not_computed"
for (i in 1:length(maxquant_clean$protein)){
  maxquant_clean$protein_group[i] <- paste(maxquant_clean$protein[i],"_",maxquant_clean$positions_within_proteins[i],maxquant_clean$multiplicity[i],sep = "")
}

maxquant_clean <- maxquant_clean %>% 
  relocate(protein_group, .before = protein)





