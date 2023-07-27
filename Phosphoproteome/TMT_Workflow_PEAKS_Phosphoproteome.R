## TMT Workflow

### libraries and WD----
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
#----

### Data Importation----
# Raw peaks data
peaks <- read.csv2("./raw_data/peptide.csv", header = T, check.names = F, sep = ",", dec = ".")

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
#----

### Function deffinition----
source("./functions/general")
source("./functions/TMT_PEAKS")
#----

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

# Transform all zero intensity values to NA
peaks <- zero_to_NA(patterns = "^plex", dataset = peaks)

##### !!!!!!!!!! HEYY CUIDAO AQUI !!!!!!!!!!!!!!!!
#peaks <- peaks[peaks$accession != "",]

# Restore annotation missing values
# Get db without /n 
db <- readLines("./raw_data/SP_Homo sapiens_Canonical_20230315.fasta")

#db.headers <- grep(">", db)
#db.f <- c()
#for (i in 1:(length(db.headers) - 1)) {
  #first <- db.headers[i] + 1
  #last <- db.headers[i + 1] - 1
  #seqs <- paste(db[first:last], collapse = "")
  #db.f <- c(db.f, db[db.headers[i]], seqs)
  #cat("\r", paste(i, (length(db.headers) - 1), sep = "/"))
  #}

#db.f <- c(db.f, db[db.headers[length(db.headers)]], paste(db[(db.headers[length(db.headers)] + 1):length(db)], collapse = ""))

# Get headers and seq apart
#db_heads <- db.f[grepl(">", db.f)]
#db_seqs <- db.f[!grepl(">", db.f)]

# Clean peptide names to do the search
#for (i in 1:length(peaks$peptide)){
  #  peaks$peptide_seq[i] <- gsub("\\s*\\([^\\)]+\\)","",peaks$peptide[i])
  #}
#peaks <- relocate(.data = peaks, peptide_seq,.after = peptide)

# Do the search and complete the databaset
#for (i in 1:length(peaks$peptide_seq)){
  #if (peaks$accession[i] == ""){
    #n <- grep(peaks$peptide_seq[i], x = db_seqs)
    #peaks$accession[i] <- db_heads[n]
    #}
  #cat("\r", paste(i, (length(db.headers) - 1), sep = "/"))
  #}

#peaks$accession <- gsub(">sp\\|","",peaks$accession)
#peaks$accession <- gsub("OS\\=.*","",peaks$accession)
#peaks$accession <- gsub(" .*","",peaks$accession)

write.csv2(x = peaks, "./raw_data/curated_peaks.csv", row.names = F)
peaks <- read.csv2("./raw_data/curated_peaks.csv")


# Create a protein group column, this will be our grouping column, if it is required it will be the fusion of two columns, if not it will be only a 
length(unique(peaks$peptide)) == length(peaks$peptide)

prot_group <- peaks$peptide
peaks$protein_group <- prot_group

peaks <- relocate(.data = peaks, 
                        protein_group, .after = accession)
#----

### Remove contaminants----
peaks_clean <- REMOVALcontaminantsPEAKS(dataset = peaks, contaminants = "./raw_data/contaminants.fasta", accession_name = "accession")
#----

### Peaks reshaping----
# Collapse redundant rows
intensities <- grep(pattern = "^plex", colnames(peaks_clean))
redundant_peaks <- peaks_clean[,c(4,c(intensities))]
redundant_peaks <- redundant_peaks %>% 
  group_by(protein_group) %>%
  distinct(redundant_peaks[,1:ncol(redundant_peaks)])

# Concatenate rows with information to keep
concatenate_peaks <- peaks_clean[,c("peptide","peptide_seq","accession","protein_group")]

concatenate_peaks <- concatenate_peaks %>%
  group_by(protein_group) %>%
  mutate(accession = paste(accession, collapse = ";"),
         peptide = paste(peptide, collapse = ";"))

concatenate_peaks <- concatenate_peaks %>%
  distinct(concatenate_peaks[1:ncol(concatenate_peaks)])

### Change Peaks to ProteinGroup
peaks_PG_clean <- merge(x = concatenate_peaks, y = redundant_peaks,
                        by = "protein_group")
#----

### Add meta-data----
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
name_number <- get_good_df(group_matrix_dataset = to_get_vect, what_is_your_tmt = "tmt16")

# Crete meta-data object
meta_data <- merge(name_number, to_get_groups, by="sample_number")

# Save meta-data
write.table(meta_data, "./results/meta_data.tsv", row.names = F, sep = "\t", dec = ".")
#----

### log2 transformation-----
# log2 transformation 
peaks_PG_clean <- log2_to_pattern(patterns = "^plex", dataset = peaks_PG_clean)

# Save median intensity value
samples <- grep(x = colnames(peaks_PG_clean), pattern = "plex")
x_num <- as.numeric(unlist(peaks_PG_clean[,samples]))
median_all <- median(x_num, na.rm = T)

# Check for pool ratios
peaks_PG_clean$pool_ratios <- peaks_PG_clean$plex1_tmt16_133n - peaks_PG_clean$plex2_tmt16_133n

#peaks_PG_clean$pool_ratios[is.na(peaks_PG_clean$pool_ratios)] <- 0

mean <- mean(peaks_PG_clean$pool_ratios, na.rm = T)
sd <- sd(peaks_PG_clean$pool_ratios, na.rm = T)

hist(peaks_PG_clean$pool_ratios, breaks = 50)
abline(v = mean + (3*sd))
abline(v = mean - (3*sd))
text(paste0("Ratio mean:",mean), y = 1500, x = 15)
#----

### Quality graph----
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
           ifelse(startsWith(group_number, "POOL"),"POOL",
                  ifelse(startsWith(group_number, "Group"), "sample",NA)))

#----

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

### Graphing----
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
#----

### Remove samples----
#----

### Normalize the data----
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
#----

### Imputation of missing values----
peaks_PG_clean_median_imp <- tim(impute = "no", dataset = peaks_PG_clean_median, NAs_prop = 0.7)
#----

### PCA 1----
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

# Plex sets and tmt sets
plexes <- substr(rownames(to_colour), 1,6)
tmts <- substr(rownames(to_colour), 13,16)

# Construct the dataset
to_colour$plex <- plexes
to_colour$tmt <- tmts
to_colour$sample_name <- rownames(to_colour)
to_colour <- merge(to_colour, meta_data, by = "sample_name")

pca1_graph <- autoplot(pca1, data = to_colour, colour = "plex",
                       frame = T)+
  scale_fill_manual(values = cbp1) +
  scale_color_manual(values = rep("black",9))
pca1_graph
#----

### Remove batch effect----
# Check the order of your columns, POOL samples are not plotted
columns_checker(peaks_PG_clean_median_imp)

# Remove batch effect
peaks_PG_clean_median_imp_unbatch <- remove_batch(dataset = peaks_PG_clean_median_imp,remove = "no", use_combat = T,
                                                  where_is_the_batch1 = c(rep("batch1",12),
                                                                          rep("batch2",12)))

write.table(file = "./results/processed_PEAKS_output.tsv", x = peaks_PG_clean_median_imp_unbatch, sep = "\t", dec = ".")
#----

### PCA 2----
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

pca2_graph <- autoplot(pca2, data = meta_data_tracker, colour = "exp_group",
                       frame = T)+
  scale_fill_manual(values = cbp1) +
  scale_color_manual(values = rep("black",9))
pca2_graph
#----

### Pair-wise----
# Check the name of your groups
unique(peaks_PG_clean_median_imp_unbatch$exp_group)

# Plot PCAs
pair_wise_pca(dataset = peaks_PG_clean_median_imp_unbatch, group1 = "Combo_R",
              group2 = "ON123300_R", column = "exp_group", to_colour = to_colour)
#----

### limma----
# Retrieve expression matrix
expression_matrix <- as.data.frame((reshape2::dcast(peaks_PG_clean_median_imp_unbatch, 
                                          protein_group ~ sample_name,value.var="unbatched_intensity", fun.aggregate = median)))

rownames(expression_matrix) <- expression_matrix[,1]
expression_matrix <- expression_matrix[,-1]
expression_matrix <- expression_matrix %>% 
  mutate_if(is.character, as.numeric)

# Create design and contrast matrix
meta_data_tracker <- subset(to_colour, sample_name %in% colnames(expression_matrix))
meta_data_tracker <- meta_data_tracker[match(colnames(expression_matrix), meta_data_tracker$sample_name)]
meta_data_tracker <- subset(meta_data_tracker, select = sample_name)
meta_data_tracker <- merge(meta_data_tracker, meta_data)

groups <- meta_data_tracker$exp_group
design <- model.matrix(~0 + groups)
colnames(design) <- gsub("^groups", "", colnames(design))
colnames(design) <- gsub(" ","_", colnames(design))
design

# Annotation
annotation <- subset(peaks_PG_clean, select = c(peptide,peptide_seq,accession,protein_group))
annotation <- distinct(annotation)
rownames(annotation) <- annotation$protein_group

# Fit the model
fit <- lmFit(expression_matrix, design = design)
#----

### limma contrasts----
# Prepare all possile contrasts
contrasts_all <- make_all_contrasts(design)

# Check if contrasts are well made
colnames(contrasts_all)

# Reverse the desired contrasts
contrasts_all <- upsidedown(contrasts_all,
                            comparisons_to_change = c("Combo_R_vs_Ibrutinib_R",
                                                      "Combo_R_vs_ON123300_R",
                                                      "Control_R_vs_Ibrutinib_R",
                                                      "Control_R_vs_ON123300_R"))

# This should show the desired contrasts
contrasts_all
#----

### Fit the model----
fit1 <- contrasts.fit(fit = fit, contrasts = contrasts_all)
fit1 <- eBayes(fit = fit1)
#----

### Extract TopTables----
# For all contrasts
tt <- topTable(fit1 ,number = Inf)

tt <- merge(tt, annotation, by = "row.names")
colnames(tt)[1] <- "protein_group"
tt <- tt[, !duplicated(colnames(tt))]
tt <- relocate(.data = tt, c(peptide,accession), .after = protein_group)
write.table(tt, "./results/TopTable_all.tsv", row.names = F, sep = "\t", dec = ".")

# Extract all Top Tables
tt_all <- tt_extractor(fit1 = fit1, annotation = annotation)

# Clean TT to share
tt_all <- tt_list_cleaner(list = tt_all, meta_data = meta_data)

write.xlsx(tt_all, file = "./results/final_toptables.xlsx")
#----

### Heatmap----
# Filter the data
ttplot <- tt
ttplot <- ttplot %>% 
  filter(P.Value < 0.05)

# Create rownnames matrix for the plot
heat_matrix <- as.matrix(ttplot[,grep("vs",colnames(ttplot))])
rownames(heat_matrix) <- ttplot$accession

# Plot the heatmap
all_heatmap <- heatmap.2(x = heat_matrix,
                         trace = "none", density.info = "none",
                         main = "Differential protein expression", scale = "row", )

#----

ggsave("./plots/intensity_detection_raw.tiff", intensity_boxplots)

ggsave("./plots/completeness.tiff", plot = completeness_barplot)

ggsave("./plots/number_of_proteins.tiff", plot = number_of_prot)

ggsave("./plots/na_density.tiff", plot = na_density)

ggsave("./plots/intensity_detection_normalized.tiff", plot = intensity_boxplots)

ggsave("./plots/pca1.tiff", plot = pca1_graph)

ggsave("./plots/pca2.tiff", plot = pca2_graph)

ggsave("./plots/pheatmap_all.tiff", plot = all_heatmap)
