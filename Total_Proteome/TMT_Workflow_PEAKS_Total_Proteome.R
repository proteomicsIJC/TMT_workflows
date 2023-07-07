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
source("./functions/TMT_PEAKS")
#----

### Get the data----
peaks <- read.csv2("./raw_data/proteins.csv", header = T, check.names = F, sep = ",", dec = ".")

# Meta data
to_get_vect <- read.csv2("./raw_data/to_get_vector.csv", sep = ",", header = T,row.names = 1, check.names = F)
to_get_groups <- read.csv2("./raw_data/to_get_name.csv", sep = ",", header = T)

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


### Get the correct columns----
# Retrive only intensities, and protein annotation columns
intensities <- grep(pattern = "^Intensity", colnames(peaks))
peaks <- peaks[,c(c(1,3),
                  c(intensities),
                  c(ncol(peaks)))]
colnames(peaks)
# Remove intensity columns consisting on the combination of samples
peaks <- subset(peaks, select = -c(grep("\\;",colnames(peaks)),
                                   grep("\\(",colnames(peaks))))
colnames(peaks)
# Clear dataset colnames
peaks <- clean_names(peaks)
colnames(peaks) <-gsub("intensity_","",colnames(peaks))
colnames(peaks)

# Make protein_group as character
peaks$protein_group <- as.character(peaks$protein_group)

# Transform all zero intensity values to NA
peaks <- zero_to_NA(patterns = "^plex", dataset = peaks)

# Filter base intensity
## peaks <- base_to_NA(patterns = "^plex", dataset = peaks)
#----

### Remove contaminants----
peaks_clean <- REMOVALcontaminantsPEAKS(dataset = peaks, contaminants = "./raw_data/contaminants.fasta", accession_name = "accession")
#----

### PEAKS reshape----
# Obtain a list consisting of all the accession names in the dataset
g_acs <- strsplit(peaks_clean$accession, "\\|")

# Extract correct accession names
good_names <- beforelast(g_acs)

# Introduce the data to the dataset
peaks_clean$good_accession <- good_names
peaks_clean <- peaks_clean %>%
  relocate(good_accession, .before = accession)

# Visualize the generation of the new protein groups
head(peaks_clean[,c("accession","good_accession")])

# Drop accession column and rename good_accession as accession
peaks_clean <- (subset(peaks_clean, select = -c(accession)))
colnames(peaks_clean)[colnames(peaks_clean) == "good_accession"] <- "accession"
colnames(peaks_clean)
head(peaks_clean$accession)

# Collapse redundant rows
intensities <- grep(pattern = "^plex", colnames(peaks_clean))
redundant_peaks <- peaks_clean[,c(1,c(intensities))]
redundant_peaks <- redundant_peaks %>% 
  group_by(protein_group) %>%
  distinct(redundant_peaks[,1:ncol(redundant_peaks)])


# Concatenate rows with information to keep
concatenate_peaks <- peaks_clean[,c("protein_group","accession","description")]

concatenate_peaks <- concatenate_peaks %>%
  group_by(protein_group) %>%
  mutate(accession = paste(accession, collapse = ";"),
         description = paste(description, collapse = ";"))

concatenate_peaks <- concatenate_peaks %>%
  distinct(concatenate_peaks[1:ncol(concatenate_peaks)])

### Change Peaks to ProteinGroup
peaks_PG_clean <- merge(x = concatenate_peaks, y = redundant_peaks,
                        by = "protein_group")
#----

### Meta-data----
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

### log2 transformation----
# log2 transformation
peaks_PG_clean <- log2_to_pattern(patterns = "^plex", dataset = peaks_PG_clean)

# save median value
samples <- grep(x = colnames(peaks_PG_clean), pattern = "plex")
samples <- intersect(colnames(peaks_PG_clean[samples]), meta_data$sample_name)
x_num <- as.numeric(unlist(peaks_PG_clean[,samples]))
median_all <- median(x_num, na.rm = T)
#----

### Data quality----
# Change data to long format
long_format <- peaks_PG_clean %>%
  pivot_longer(cols = starts_with("plex"),
               names_to = "sample_name",
               values_to = "intens")

# Remove colnames description as won't be used for the graphs' sake
long_format <- subset(long_format, select=-c(description))

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

### Remove dataset----
long_format <- remove_samp(dataset = long_format)
#----


### Normalization by median----
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

### Missing values iputation----
peaks_PG_clean_median_imp <- tim(impute = "no", dataset = peaks_PG_clean_median, NAs_prop = 0.5)
#----

### PCA 1----
# Change long to wide format
peaks_PG_clean_median_imp_topca <- data.table::dcast(peaks_PG_clean_median_imp, 
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

pca1_graph <- autoplot(pca1, data = to_colour, colour = "exp_group",
                       frame = T)+
  scale_fill_manual(values = cbp1) +
  scale_color_manual(values = rep("black",9))
pca1_graph
#----

### Batch effect removal----
# Check the order of your columns, POOL samples are not plotted
columns_checker(peaks_PG_clean_median_imp)

# Remove batch effect
peaks_PG_clean_median_imp_unbatch <- remove_batch(dataset = peaks_PG_clean_median_imp,remove = "no")

write.table(file = "./results/processed_PEAKS_output.tsv", x = peaks_PG_clean_median_imp_unbatch, sep = "\t", dec = ".")
#----


### PCA 2----
# Cahnge long to wide format
peaks_PG_clean_median_imp_topca2 <- data.table::dcast(peaks_PG_clean_median_imp_unbatch, 
                                          protein_group
                                          ~ sample_name,value.var="unbatched_intensity",
                                          fun.aggregate = median)

peak_mat2 <- t(peaks_PG_clean_median_imp_topca2[,c(2:ncol(peaks_PG_clean_median_imp_topca2))])

# PCA
pca2 <- prcomp(peak_mat2, scale. = TRUE, center = TRUE)

# to_colour
samples_afer_batch <- intersect(rownames(peak_mat2), to_colour$sample_name)
to_colour2 <- subset(to_colour, sample_name %in% c(samples_afer_batch))

pca2_graph <- autoplot(pca2, data = to_colour2, colour = "exp_group",
                       frame = T)+
  scale_fill_manual(values = cbp1) +
  scale_color_manual(values = rep("black",9))
pca2_graph
#----


##########
# limma #
#########

### Enter data to limma---- 
# Retrive expression matrix
expression_matrix <- as.data.frame((data.table::dcast(peaks_PG_clean_median_imp_unbatch, 
                                          protein_group ~ sample_name,value.var="unbatched_intensity", fun.aggregate = median)))

rownames(expression_matrix) <- expression_matrix[,1]
expression_matrix <- expression_matrix[,-1]
expression_matrix <- expression_matrix %>% 
  mutate_if(is.character, as.numeric)

# Create design and contrast matrix
to_colour2 <- subset(to_colour, sample_name %in% colnames(expression_matrix))
to_colour2 <- to_colour2[match(colnames(expression_matrix), to_colour2$sample_name)]
to_colour2 <- subset(to_colour2, select = sample_name)
to_colour2 <- merge(to_colour2, meta_data)

groups <- to_colour2$exp_group
design <- model.matrix(~0 + groups)
colnames(design) <- gsub("^groups", "", colnames(design))
colnames(design) <- gsub(" ","_", colnames(design))
design

# Annotation
annotation <- subset(peaks_PG_clean, select = c(protein_group,accession,description))
annotation <- distinct(annotation)
rownames(annotation) <- annotation$protein_group
annotation <- annotation[,-1]

# Fit the model
fit <- lmFit(expression_matrix, design = design)
#----

### Perform all possible samples----
# Prepare all possile contrasts
contrasts_all <- make_all_contrasts(design)

# Reverse the desired contrasts
contrasts_all <- upsidedown(contrasts_all, comparisons_to_change = c("Combo_vs_Ibrutinib",
                                                                     "Combo_vs_ON123300",
                                                                     "Control_vs_Ibrutinib",
                                                                     "Control_vs_ON123300"))
contrasts_all
#----


### Fit contrasts to the model----
fit1 <- contrasts.fit(fit = fit, contrasts = contrasts_all)
fit1 <- eBayes(fit = fit1)
#----


### Extract limma data----
tt <- topTable(fit1, number = Inf)
tt <- merge(tt, annotation, by = "row.names")
colnames(tt)[1] <- "protein_group"
tt <- relocate(.data = tt, c(accession,description), .after = protein_group)
write.table(tt, "./results/TopTable_all.tsv", row.names = F, sep = "\t", dec = ".")

# Extract all Top Tables
tt_all <- tt_extractor(fit1 = fit1, annotation = annotation)

# Clean TT to share
tt_all <- tt_list_cleaner(list = tt_all, meta_data = meta_data)

# Extract xlsx finalTopTables
xlsx_tt(tt_cleaned_list = tt_all)
#----

#### Heatmap----
cat("Data has been filtered folowing a P.Value threshold of 0.01")

# Filter the data
ttplot <- tt
ttplot <- ttplot %>% 
  filter(P.Value < 0.01)

ttplot <- merge(ttplot, batching, by.x = "protein_group", by.y = "row.names")

heat_matrix <- as.matrix(ttplot[,grep("plex",colnames(ttplot))])
heat_matrix <- as.data.frame(heat_matrix)

# set colnames
old_names <- colnames(heat_matrix)[colnames(heat_matrix) %in% to_colour2$sample_name]
new_names <- c()

for ( k in 1:length(old_names)){
  grup <- meta_data$exp_group[meta_data$sample_name == old_names[k]]
  sample_number <- meta_data$sample_number[meta_data$sample_name == old_names[k]]
  sample_number <- gsub(x = sample_number, "sample_","")
  new_names[k] <- paste(grup,sample_number,sep = ".")
}

name_changer <- data.frame(old = old_names,
                           new = new_names)

setnames(heat_matrix, old =  name_changer$old, name_changer$new)

heat_matrix <- as.matrix(heat_matrix)

# Plot the heatmap
all_heatmap <- heatmap.2(x = heat_matrix,
                         trace = "none", density.info = "none",
                         main = "Differential Protein Group expression", scale = "row",cexCol = 0.75, cexRow = 0.75,
                         margins = c(7,2), col = colorRampPalette(c("red", "white", "blue"))(n = 256))
#-----



ggsave("./plots/intensity_detection_raw.tiff", intensity_boxplots, width = 10, height = 10)

ggsave("./plots/completeness.tiff", plot = completeness_barplot, width = 10, height = 10)

ggsave("./plots/number_of_proteins.tiff", plot = number_of_prot, width = 10, height = 10)

ggsave("./plots/na_density.tiff", plot = na_density, width = 10, height = 10)

ggsave("./plots/intensity_detection_normalized.tiff", plot = intensity_boxplots, width = 10, height = 10)

ggsave("./plots/pca1.tiff", plot = pca1_graph, width = 10, height = 10)

ggsave("./plots/pca2.tiff", plot = pca2_graph, width = 10, height = 10)

png("./plots/pheatmap_all.png",all_heatmap)


