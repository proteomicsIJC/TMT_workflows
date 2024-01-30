#################################
### TMT Mod Worflow with PEAKS ##
#################################

#### This script is based on 029

### libraries and WD
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
setwd(dirname(getActiveDocumentContext()$path))

### functions
# general
source("./functions/general/")

# TMT peaks
source("./functions/TMT_PEAKS/")

### Data Importation
# Raw peaks data
peaks <- read.csv2("./raw_data/peptide_ADPR.csv", header = T, check.names = F, sep = ",", dec = ".")

# Meta data
to_get_vect <- read.csv2("./raw_data/to_get_vector.csv", sep = ",", header = T,row.names = 1, check.names = F)
to_get_groups <- read.csv2("./raw_data/to_get_name.csv", sep = ",", header = T, check.names = T)

# Contaminants
cont <- readLines("./raw_data/contaminants.fasta")

# Check if folders exist, if NO create them
wd <- getwd()
dir.create(file.path(wd,"./raw_data"))
dir.create(file.path(wd,"./results"))
dir.create(file.path(wd,"./plots"))
file.remove(file.path(wd,"./results/used_parameters.txt"))
file.create(file.path(wd, "./results/used_parameters.txt"))


### Retrieve only the columns for the analysis and clean the data----
# Retrive only intensities, and protein annotation columns
intensities <- grep(pattern = "^Intensity", colnames(peaks))
peaks <- peaks[,c(c(1,12),
                  c(intensities))]

# Remove intensity columns consisting on the combination of samples
peaks <- subset(peaks, select = -c(grep("\\;",colnames(peaks)),
                                   grep("\\(",colnames(peaks))))

# Clear dataset colnames
peaks <- clean_names(peaks)
colnames(peaks) <-gsub("intensity_","",colnames(peaks))
for (i in 1:length(colnames(peaks))){
  colnames(peaks)[i] <- gsub(pattern = "tm_tset", replacement = "plex", colnames(peaks)[i])
}

# Transform all zero intensity values to NA
peaks <- zero_to_NA(patterns = "^plex", dataset = peaks)

### Restore annotation missing values
##### !!!!!!!!!! HEYY CUIDAO AQUI !!!!!!!!!!!!!!!!
##peaks <- peaks[peaks$accession != "",]

# Get db without /ns 
# Specify the path to your FASTA file
fasta_file <- "./raw_data/SP_Mouse_20210426.fasta"

# Read the file using readLines
fasta_lines <- readLines(fasta_file)

# Initialize empty vectors to store headers and sequences
headers <- character(0)
sequences <- character(0)

# Iterate through the lines and separate headers and sequences
current_sequence <- NULL
for (line in fasta_lines) {
  if (substr(line, 1, 1) == ">") {
    # This line is a header
    headers <- c(headers, line)
    if (!is.null(current_sequence)) {
      sequences <- c(sequences, paste(current_sequence, collapse = ""))
    }
    current_sequence <- character(0)
  } else {
    # This line is part of a sequence
    current_sequence <- c(current_sequence, line)
  }
}

# Don't forget to add the last sequence
if (!is.null(current_sequence)) {
  sequences <- c(sequences, paste(current_sequence, collapse = ""))
}

# Create a data frame with headers and sequences
fasta_data <- data.frame(Header = headers, Sequence = sequences)

# Print the first few records to verify the results
head(fasta_data)

# Gene names and description
fasta_data$Description <- "Not done"
fasta_data$Gene_name <- "Not done"
fasta_data$Accession <- "Not done"
for (i in 1:length(rownames(fasta_data))){
  # Gene name
  fasta_data$Gene_name[i] <- sub(".*GN=([^ ]*).*","\\1", fasta_data$Header[i])
  # Description
  fasta_data$Description[i] <- sub(".*MOUSE(.*?)OS.*","\\1", fasta_data$Header[i])
  fasta_data$Description[i] <- trimws(fasta_data$Description[i])
  
  # Accession
  fasta_data$Accession[i] <- sub("^[^|]*\\|([^|]*).*", "\\1", fasta_data$Header[i])
  
  cat("\r", paste(i, (length(fasta_data$Header)), sep = "/"))
}

# Clean peptide names to do the search
for (i in 1:length(peaks$peptide)){
  peaks$peptide_seq[i] <- gsub("\\s*\\([^\\)]+\\)","",peaks$peptide[i])
  cat("\r", paste(i, (length(peaks$peptide)), sep = "/"))
}

peaks <- relocate(.data = peaks, peptide_seq,.after = peptide)

# Do the search and complete the databaset
## peptide accessions
for (i in 1:length(peaks$peptide_seq)){
  if (peaks$accession[i] == ""){
    n <- grep(peaks$peptide_seq[i], x = fasta_data$Sequence)
    peaks$accession[i] <- fasta_data$Accession[n]}
  cat("\r", paste(i, (length(peaks$peptide_seq)), sep = "/"))
}

## Now time for the accession extraction
peaks$good_accession <- "Not done"
for (i in 1:length(peaks$accession)){
  peaks$good_accession[i] <- sub("^([^|]*)\\|.*", "\\1", peaks$accession[i])
}
peaks <- peaks %>%
  relocate(good_accession, .before = accession)

## Peptide Gene_names and Descriptions
annotation_later <- peaks %>% 
  subset(select = c(good_accession)) %>% 
  distinct()
annotation_later$Gene_name <- "Not done"
annotation_later$Description <- "Not done"

for (i in 1:length(annotation_later$good_accession)){
  n <- grep(pattern = annotation_later$good_accession[i], x = fasta_data$Accession)
  annotation_later$Gene_name[i] <- fasta_data$Gene_name[n]
  annotation_later$Description[i] <- fasta_data$Description[n]
  cat("\r", paste(i, (length(annotation_later$good_accession)), sep = "/"))
}

# merge annotation_later with peaks
peaks <- merge(peaks, annotation_later,by = "good_accession")
peaks <- peaks %>% 
  relocate(all_of(colnames(annotation_later)), .after = peptide_seq)
colnames(peaks)[3] <- "Accession"

# create the Uniprot_ID column
colnames(peaks)[6] <- "Uniprot_ID"
for (i in 1:length(peaks$Uniprot_ID)){
  peaks$Uniprot_ID[i] <- sub("^[^|]*\\|(.*)", "\\1", peaks$Uniprot_ID[i])
  cat("\r", paste(i, (length(peaks$Uniprot_ID)), sep = "/"))
}

# remove annotation_later
remove(annotation_later)

##peaks$accession <- gsub(">sp\\|","",peaks$accession)
##peaks$accession <- gsub("OS\\=.*","",peaks$accession)
##peaks$accession <- gsub(" .*","",peaks$accession)

### write.csv2(x = peaks, "./raw_data/curated_peaks.csv", row.names = F)
### peaks <- read.csv2("./raw_data/curated_peaks.csv")

### Create a protein group column 
# This will be our grouping column, if it is required it will be the fusion of two columns, 
# if not it will be only a copy of the peptide column 
length(unique(peaks$peptide)) == length(peaks$peptide)

prot_group <- peaks$peptide
peaks$protein_group <- prot_group

peaks <- relocate(.data = peaks, 
                  protein_group, .after = Accession)

### Remove contaminants
peaks_clean <- REMOVALcontaminantsPEAKS(dataset = peaks, contaminants = "./raw_data/contaminants.fasta", accession_name = "Accession")

### peaks_clean <- peaks
### Peaks reshaping
# Collapse redundant rows
intensities <- grep(pattern = "^plex", colnames(peaks_clean))
redundant_peaks <- peaks_clean[,c(4,c(intensities))]
redundant_peaks <- redundant_peaks %>% 
  group_by(protein_group) %>%
  distinct(redundant_peaks[,1:ncol(redundant_peaks)])

# Concatenate rows with information to keep
concatenate_peaks <- peaks_clean[,c("peptide","peptide_seq","Accession","protein_group","Gene_name","Description","Uniprot_ID")]

concatenate_peaks <- concatenate_peaks %>%
  group_by(protein_group) %>%
  mutate(accession = paste(Accession, collapse = ";"),
         peptide = paste(peptide, collapse = ";"))

concatenate_peaks <- concatenate_peaks %>%
  distinct(concatenate_peaks[1:ncol(concatenate_peaks)])

### Change Peaks to ProteinGroup
peaks_PG_clean <- merge(x = concatenate_peaks, y = redundant_peaks,
                        by = "protein_group")

### Add meta-data
# Clean the group assigner matrix
to_get_groups <- to_get_groups %>%
  dplyr::rename(sample_number = 1, group_number = 2, exp_group = 3)
to_get_groups$sample_number <- as.character(to_get_groups$sample_number) 

for (k in 1:length(to_get_groups$sample_number)){
  if (!startsWith(x = to_get_groups$sample_number[k], prefix = "POOL")){
    to_get_groups$sample_number[k] <- paste("sample",to_get_groups$sample_number[k],sep = "_")  
  }
}

# Get sample_numbers and names correspondence
name_number <- create_meta_data(group_matrix_dataset = to_get_vect, what_is_your_tmt = "tmt16")

# Crete meta-data object
meta_data <- merge(name_number, to_get_groups, by="sample_number")
meta_data$sample_name <- tolower(meta_data$sample_name)

# Save meta-data
write.table(meta_data, "./results/meta_data.tsv", row.names = F, sep = "\t", dec = ".")

### log2 transformation
# log2 transformation 
peaks_PG_clean <- log2_to_pattern(patterns = "^plex", dataset = peaks_PG_clean)

# Save median intensity value
samples <- grep(x = colnames(peaks_PG_clean), pattern = "plex")
x_num <- as.numeric(unlist(peaks_PG_clean[,samples]))
median_all <- median(x_num, na.rm = T)

### Quality graphs
# Change data to long format
long_format <- peaks_PG_clean %>%
  pivot_longer(cols = starts_with("plex"),
               names_to = "sample_name",
               values_to = "intens")

# Add median calculation (by sample)
long_format_median <- long_format %>%
  group_by(sample_name) %>%
  summarise(MED = median(intens, na.rm=T))

# Assign if median values are high
long_format_median <- long_format_median %>%
  mutate(high = ifelse(abs(MED-mean(MED)) >= 1.5, "YES","NO"))

long_format <- merge(x = long_format, y = long_format_median,
                     by = "sample_name")

long_format <- long_format  %>% 
  group_by(sample_name)%>%
  mutate(is_max = 
           ifelse(intens == max(intens, na.rm = T), "YES","NO")) %>%
  ungroup()

# Assign experimental group (name and number)
long_format <- merge(long_format, meta_data, by = "sample_name")

# Assign TMT Plex
long_format <- long_format %>%
  mutate(plex = 
           substr(sample_name, 1,5))

# Assign TMT reactive
long_format <- long_format %>%
  mutate(tmt = 
           substr(sample_name, 13,16))

# Assign TMT sample or pool
long_format <- long_format %>%
  mutate(sample_or_pool = 
           ifelse(startsWith(group_number, "POOL"),"POOL","sample"))

### General Phospho state----
# Prepare the data
#peptides <- peaks_PG_clean$peptide
#total_peptides <- length(peaks_clean$peptide)

# Quesitos
# POSAR DINS DEL DATASET !!!!!1

#phosho_counter = 0
#for (i in 1:length(peptides)){
#  if (grepl("+79.97", peptides[i], fixed = TRUE)){
#    phosho_counter = phosho_counter + 1
#  }}
#phosho_counter
#no_phospho <- total_peptides - phosho_counter
## https://r-coder.com/pie-chart-r/?utm_content=cmp-true

# Bars
# FER COM AQUI !!!
## https://stackoverflow.com/questions/12427385/how-to-calculate-the-number-of-occurrence-of-a-given-character-in-each-row-of-a
#-----

### Graphs
# Set a colour palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
          "#00B159", "#FCD612", "#FF2F03", "#03D3FF",
          "#2FF923", "#FF03B4", "#f08a0c", "#0d407f",
          "#037A68", "#510840", "#F70D1A", "#e0218a")

# Boxplot
intensity_boxplots <- ggplot(long_format, mapping = aes(x = sample_name, y = intens))+
  geom_boxplot(data = long_format, mapping = aes(x = sample_name, y = intens, fill = plex))+
  #  geom_text_repel(data = long_format[which((long_format$high == "YES") &
  #                                          (long_format$is_max == "YES")),], aes(label = sample_name, y = max(intens+2), x = sample_name))+
  theme_bw()+
  scale_fill_manual(values = cbp1)+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

intensity_boxplots  

# Completenes
completeness <- long_format %>%
  group_by(sample_name) %>%
  summarize(count_na = 100 - (((sum(is.na(intens))/nrow(peaks_PG_clean))))*100)
completeness <- merge(completeness, meta_data, by = "sample_name")
completeness <- completeness %>%
  mutate(plex = 
           substr(sample_name, 1,5))

completeness_barplot <- ggplot(completeness, mapping = aes(y = count_na, x = sample_name))+
  geom_bar(stat = "identity", aes(fill = plex))+
  geom_hline(data = completeness ,mapping = aes(yintercept = mean(count_na), col = "red"))+
  theme_bw()+
  scale_fill_manual(values = cbp1)+
  xlab("Sample") +
  ylab("% of completeness")+
  ggtitle("Completeness of the data per sample")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
completeness_barplot

# Number of proteins per sample
nprot <- long_format %>%
  group_by(sample_name) %>%
  summarize(count_prot = sum(!is.na(intens)))
nprot <- merge(nprot, meta_data, by = "sample_name")
nprot <- nprot %>%
  mutate(plex = 
           substr(sample_name, 1,5))

number_of_prot <- ggplot(nprot, mapping = aes(y = count_prot, x = sample_name))+
  geom_bar(stat = "identity", aes(fill = plex))+
  geom_hline(data = nprot ,mapping = aes(yintercept = mean(count_prot), col = "red"))+
  theme_bw()+
  scale_fill_manual(values = cbp1)+
  xlab("Sample") +
  ylab("# of proteins")+
  ggtitle("Proteins per sample")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
number_of_prot

# NA per protein
naprot <- long_format %>%
  group_by(protein_group) %>%
  summarize(count_prot = sum(is.na(intens)))

# NAs density per protein
na_density <- ggplot(data = naprot, mapping = aes(x = count_prot))+
  geom_histogram(binwidth = 1, bins = 1)+
  theme_bw()+
  ggtitle("NAs density per protein groups")+
  xlab("# NAs")+
  ylab("# proteins")+
  theme(legend.position= "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
na_density

### Remove samples
long_format <- remove_samp(dataset = long_format)

### Normalize the data
# Use long format data to work and remove two non-used columns
peaks_PG_clean_median <- long_format
peaks_PG_clean_median <- subset(peaks_PG_clean_median, select = -c(high,is_max))

# Do the calculation
peaks_PG_clean_median$normalized_intensity <- (peaks_PG_clean_median$intens - peaks_PG_clean_median$MED) + median_all 

# Boxplot graph for normalized values
intensity_boxplots_norm <- ggplot(peaks_PG_clean_median, mapping = aes(x = sample_name, y = normalized_intensity))+
  geom_boxplot(data = peaks_PG_clean_median, mapping = aes(x = sample_name, y = normalized_intensity, fill = plex))+
  theme_bw()+
  scale_fill_manual(values = cbp1)+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection (Normalized intensities)")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

intensity_boxplots_norm

### Imputation of missing values
peaks_PG_clean_median_imp <- tim(impute = "no", dataset = peaks_PG_clean_median, NAs_prop = 0.7)

### PCA 1 DO THE PCA THE NEW WAYY !!!
# Change long to wide format
peaks_PG_clean_median_imp_topca <- reshape2::dcast(peaks_PG_clean_median_imp, 
                                                   protein_group
                                                   ~ sample_name,value.var="normalized_intensity",
                                                   fun.aggregate = median)

peak_mat <- t(peaks_PG_clean_median_imp_topca[,c(2:ncol(peaks_PG_clean_median_imp_topca))])


# PCA
pca1 <- prcomp(peak_mat, scale. = TRUE, center = TRUE)

# PCA colouring
to_colour <- as.data.frame(peak_mat)
to_colour <- as.data.frame(peak_mat[,c(1), drop = F])
colnames(to_colour)[1] <- "sample_name"
to_colour$sample_name <- rownames(to_colour)

# Plex sets and tmt sets
plexes <- substr(rownames(to_colour), 1,6)
tmts <- substr(rownames(to_colour), 13,16)

# Construct the dataset
to_colour$plex <- plexes
to_colour$tmt <- tmts
to_colour <- merge(to_colour, meta_data, by = "sample_name")
rownames(to_colour) <- to_colour$sample_name
to_colour <- to_colour[rownames(peak_mat),]

# PCA data construction
## positional data
pca1x <- as.data.frame(pca1$x)

## eigenvalues
eigenvalues <- pca1$sdev^2
pca1x$eigenvalues <- eigenvalues

## construct the data for colouring
pca1x <- pca1x[to_colour$sample_name,]
pca1xx <- merge(pca1x, to_colour,by = "row.names")
rownames(pca1xx) <- pca1xx$Row.names
pca1xx <- pca1xx[,-1]

# PCA graph
pca_plot <- ggplot(pca1xx, 
                   aes(x = PC1, y = PC2, label = exp_group)) +
  
  geom_point(size = 4, aes(color = plex)) +
  scale_color_manual(values = c("plex1_" = "blue", "plex2_" = "red"))+
  guides(color=guide_legend(title = "Tmt plex"))+
  geom_text_repel(nudge_y = 0.01, 
                  size = 6)+
  xlab(paste("PC1(", round(100*(pca1$sdev[1]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca1$sdev[2]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA)"))+
  geom_polygon(aes(group = plex, fill = plex), alpha = 0.2, show.legend = FALSE)+
  scale_fill_manual(values = c("plex1_" = "blue", "plex2_" = "red"))
pca_plot

### Remove batch effect
# Check the order of your columns, POOL samples are not plotted
columns_checker(peaks_PG_clean_median_imp)

# Remove batch effect
peaks_PG_clean_median_imp_unbatch <- remove_batch(dataset = peaks_PG_clean_median_imp,remove = "yes", use_combat = T,
                                                  where_is_the_batch1 = c(rep("batch1",12),
                                                                          rep("batch2",12)))

write.table(file = "./results/processed_PEAKS_output.tsv", x = peaks_PG_clean_median_imp_unbatch, sep = "\t", dec = ".")

### PCA 2 DO THE PCA THE NEW WAYY !!!
# Cahnge long to wide format
peaks_PG_clean_median_imp_topca2 <- reshape2::dcast(peaks_PG_clean_median_imp_unbatch, 
                                                    protein_group
                                                    ~ sample_name,value.var="unbatched_intensity",
                                                    fun.aggregate = median)

peak_mat2 <- t(peaks_PG_clean_median_imp_topca2[,c(2:ncol(peaks_PG_clean_median_imp_topca2))])

# PCA
pca2 <- prcomp(peak_mat2, scale. = TRUE, center = TRUE)

# to_colour
samples_afer_batch <- intersect(rownames(peak_mat2), to_colour$sample_name)
meta_data_tracker <- subset(to_colour, sample_name %in% c(samples_afer_batch))

# PCA data construction
## positional data
pca2x <- as.data.frame(pca2$x)

## eigenvalues
eigenvalues <- pca2$sdev^2
pca2x$eigenvalues <- eigenvalues

## data construction
pca2xx <- merge(pca2x, to_colour, by = "row.names")
pca2xx <- pca2xx[,-1]

## PCA representation
pca_plot2 <- ggplot(pca2xx, 
                    aes(x = PC1, y = PC2, label = exp_group)) +
  
  geom_point(size = 4, aes(color = exp_group)) +
  ###scale_color_manual(values = c("plex1_" = "blue", "plex2_" = "red"))+
  guides(color=guide_legend(title = "Tmt plex"))+
  geom_text_repel(nudge_y = 0.01, 
                  size = 6)+
  xlab(paste("PC1(", round(100*(pca1$sdev[1]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca1$sdev[2]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA)"))+
  geom_polygon(aes(group = exp_group, fill = exp_group), alpha = 0.2, show.legend = FALSE)
###scale_fill_manual(values = c("plex1_" = "blue", "plex2_" = "red"))
pca_plot2


### limma analysis
# Retrieve expression matrix
expression_matrix <- as.data.frame((reshape2::dcast(peaks_PG_clean_median_imp_unbatch, 
                                                    protein_group ~ sample_name,value.var="unbatched_intensity", fun.aggregate = median)))

rownames(expression_matrix) <- expression_matrix[,1]
expression_matrix <- expression_matrix[,-1]
expression_matrix <- expression_matrix %>% 
  mutate_if(is.character, as.numeric)
expression_matrix <- expression_matrix[,meta_data_tracker$sample_name]

###### TRANSFORM SAMPLE NAME TO A LEGIBLE NAME
## legible name
all(colnames(expression_matrix) == meta_data_tracker$sample_name)
meta_data_tracker$legible_name <- meta_data_tracker$exp_group

rownames(meta_data_tracker) <- NULL
for (i in 1:length(meta_data_tracker$legible_name)){
  ## CONTROL
  meta_data_tracker$legible_name[i] <- gsub(pattern = "Control", replacement = "Ct", x = meta_data_tracker$legible_name[i])
  meta_data_tracker$exp_group[i] <- gsub(pattern = "Control", replacement = "Ct", x = meta_data_tracker$exp_group[i])
  ## SIR7
  meta_data_tracker$legible_name[i] <- gsub(pattern = "Sir7KO", replacement = "S7KO", x = meta_data_tracker$legible_name[i])
  meta_data_tracker$exp_group[i] <- gsub(pattern = "Sir7KO", replacement = "S7KO", x = meta_data_tracker$exp_group[i])
  ## SIRT7
  meta_data_tracker$legible_name[i] <- gsub(pattern = "SirT7KO", replacement = "S7KO", x = meta_data_tracker$legible_name[i])
  meta_data_tracker$exp_group[i] <- gsub(pattern = "SirT7KO", replacement = "S7KO", x = meta_data_tracker$exp_group[i])
}

meta_data_tracker$legible_name <- with(meta_data_tracker, paste(legible_name, ave(rep(1, length(exp_group)), exp_group, FUN = seq_along), sep = "_"))
rownames(meta_data_tracker) <- meta_data_tracker$sample_name
meta_data_tracker <- meta_data_tracker[colnames(expression_matrix),]
all(rownames(meta_data_tracker) == colnames(expression_matrix))
colnames(expression_matrix) <- meta_data_tracker$legible_name

### Annotation extraction
annotation <- peaks_PG_clean_median_imp_unbatch %>% 
  subset(select = c(protein_group, peptide, peptide_seq, Accession, Gene_name, Description, Uniprot_ID)) %>% 
  distinct()
colnames(annotation)[4] <- "ID"
colnames(annotation)[1] <- "Accession"
rownames(annotation) <- annotation$Accession

# Create design and contrast matrix
groups <- meta_data_tracker$exp_group
design <- model.matrix(~0 + groups)
colnames(design) <- gsub("^groups", "", colnames(design))
colnames(design) <- gsub(" ","_", colnames(design))
design

# Fit the model
fit <- lmFit(expression_matrix, design = design)

# Contrast matrix
cnt <- makeContrasts("S7KO_CR_old vs. S7KO_Ct_old" = S7KO_CR_old - S7KO_Ct_old,
                     "WT_Ct_old vs. WT_Ct_young" = WT_Ct_old - WT_Ct_young, 
                     levels = design)

# Fit contrast to the model
fit1 <- contrasts.fit(fit = fit, contrasts = cnt)
fit1 <- eBayes(fit = fit1)

writting <- xlsx_tt(fit__1 = fit1, meta_data = meta_data_tracker, 
                    meta_sample_column = "legible_name", meta_data_column = "exp_group", 
                    annotation = annotation, expression_matrix = expression_matrix, filename = "./results/ADPR_conntrasts_results_ep2.xlsx",
                    color_samples = c("S7KO_CR_old" = "red",
                                      "S7KO_Ct_old" = "blue",
                                      "WT_Ct_old" = "yellow",
                                      "WT_Ct_young" = "green"), differentation_element = " vs. ")

tt <- topTable(fit = fit1, number = Inf)

### limma contrasts
### HERE PUT THE CODE TO PRODUCE A CLEAN VERSION OF THE RESULTS !!!


ggsave("./plots/intensity_detection_raw.tiff", intensity_boxplots)

ggsave("./plots/completeness.tiff", plot = completeness_barplot)

ggsave("./plots/number_of_proteins.tiff", plot = number_of_prot)

ggsave("./plots/na_density.tiff", plot = na_density)

ggsave("./plots/intensity_detection_normalized.tiff", plot = intensity_boxplots)

ggsave("./plots/pca1.tiff", plot = pca_plot)

ggsave("./plots/pca2.tiff", plot = pca_plot2)

