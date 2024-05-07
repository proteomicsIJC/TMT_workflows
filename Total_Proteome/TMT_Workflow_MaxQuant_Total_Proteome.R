#############################
### TMT Standard Workflow ###
#############################

### Libraries, WD and functions
library(dplyr)
library(tidyr)
library(data.table)
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
library(reshape2)
library(readr)
setwd(dirname(getActiveDocumentContext()$path))

# Functions deffinition
## General proteomics functions
source("./functions/general/columns_checker.R")
source("./functions/general/create_meta_data.R")
source("./functions/general/log2_to_pattern.R")
source("./functions/general/make_all_contrasts.R")
source("./functions/general/presence_vs_no_presence.R")
source("./functions/general/remove_batch.R")
source("./functions/general/remove_samp.R")
source("./functions/general/tim.R")
source("./functions/general/upsidedown.R")
source("./functions/general/zero_to_na.R")
source("./functions/general/xlsx_tt.R")

## TMT_Maxquant_functions
source("./functions/TMT_MaxQuant/maxquantinitializer.R")
source("./functions/TMT_MaxQuant/proteinGroupsCleaner.R")

## Folder system 
wd <- getwd()
dir.create(file.path(wd,"./raw_data"))
dir.create(file.path(wd,"./results"))
dir.create(file.path(wd,"./plots"))
file.remove(file.path(wd,"./results/used_parameters.txt"))
file.create(file.path(wd, "./results/used_parameters.txt"))

### Data importation
# Raw data
maxquant <- read.table("./raw_data/proteinGroups.txt", header = T, check.names = F, sep = "\t", dec = ".")

# Meta data
to_get_vect <- openxlsx::read.xlsx("./raw_data/to_get_vector.xlsx", rowNames = T)
to_get_groups <- openxlsx::read.xlsx("./raw_data/to_get_name_common_cells.xlsx")
##to_get_groups <- openxlsx::read.xlsx("./raw_data/to_get_name.xlsx")

# Contaminants
cont <- readLines("./raw_data/contaminants.fasta")

### Work the data tables
# Clean the maxquant data
maxquant <- maxquant_initalizer(maxquant_data = maxquant, tmt = "tmt10", n_plex = 3)

# Retrive only intensities, and protein annotation columns
intensities <- grep(pattern = "^Intensity plex[0-9] \\: TMT", colnames(maxquant))
reverse <- grep(pattern = "^Reverse", colnames(maxquant))
potential_cont <- grep(pattern = "^Potential contaminant", colnames(maxquant))
only_site <- grep(pattern = "^Only identified by site", colnames(maxquant))
maxquant <- maxquant[,c(c(2,6,7),
                        c(intensities),
                        c(reverse,potential_cont,only_site))]
colnames(maxquant)

# Clear dataset colnames
maxquant <- clean_names(maxquant)
colnames(maxquant) <-gsub("intensity_","",colnames(maxquant))
colnames(maxquant)

# Rename majority to protein group
colnames(maxquant)[1] <- "protein_group"
# Transform all zero intensity values to NA
maxquant <- zero_to_NA(patterns = "^plex", dataset = maxquant)

### Remove contaminants
maxquant_clean <- proteinGroupsCleanner(ds = maxquant)

### Meta-data work
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
name_number <- create_meta_data(group_matrix_dataset = to_get_vect, what_is_your_tmt = "tmt10")

# Crete meta-data object
meta_data <- merge(name_number, to_get_groups, by="sample_number")
meta_data <- meta_data %>% 
  arrange(group_number)

# Save meta-data
write.table(meta_data, "./results/meta_data.tsv", row.names = F, sep = "\t", dec = ".")

### log2 transformation
# log2 transformation
maxquant_clean <- log2_to_pattern(patterns = "^plex", dataset = maxquant_clean)

# save median value
samples <- grep(x = colnames(maxquant_clean), pattern = "plex")
samples <- intersect(colnames(maxquant_clean[samples]), meta_data$sample_name)
x_num <- as.numeric(unlist(maxquant_clean[,samples]))
median_all <- median(x_num, na.rm = T)

### Data quality graphs
# Change data to long format
long_format <- maxquant_clean %>%
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
  summarize(count_na = 100 - (((sum(is.na(intens))/nrow(maxquant_clean))))*100)
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

### Normalization by median
# Use long format data to work and remove two non-used columns
maxquant_clean_median <- long_format
maxquant_clean_median <- subset(maxquant_clean_median, select = -c(high,is_max))

# Do the calculation
maxquant_clean_median$normalized_intensity <- (maxquant_clean_median$intens - maxquant_clean_median$MED) + median_all 

# Boxplot graph for normalized values
intensity_boxplots_norm <- ggplot(maxquant_clean_median, mapping = aes(x = sample_name, y = normalized_intensity))+
  geom_boxplot(data = maxquant_clean_median, mapping = aes(x = sample_name, y = normalized_intensity, fill = plex))+
  theme_bw()+
  scale_fill_manual(values = cbp1)+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection (Normalized intensities)")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

intensity_boxplots_norm

### Missing values iputation
maxquant_clean_median_imp <- tim(impute = "no", dataset = maxquant_clean_median, NAs_prop = 0, intensity_to_impute = "normalized_intensity")

### PCA 1
# Change long to wide format
maxquant_clean_median_imp_to_pca <- reshape2::dcast(maxquant_clean_median_imp, 
                                                    protein_group
                                                    ~ sample_name,value.var="normalized_intensity",
                                                    fun.aggregate = median)

peak_mat <- t(maxquant_clean_median_imp_to_pca[,c(2:ncol(maxquant_clean_median_imp_to_pca))])

# PCA calculation
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
                   aes(x = PC1, y = PC2, 
                       ####label = plex 
                   )) +
  
  geom_point(size = 4, aes(color = plex)) +
  scale_color_manual(values = c("plex1_" = "#999999",
                                "plex2_" = "#E69F00",
                                "plex3_" = "#56B4E9"))+
  
  
  guides(color=guide_legend(title = "Experimental group"))+
  ##geom_text_repel(nudge_y = 0.01, 
  ##                size = 6)+
  xlab(paste("PC1(", round(100*(pca1$sdev[1]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca1$sdev[2]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA)"))+
  geom_polygon(aes(group = plex, fill = plex), alpha = 0.2, show.legend = FALSE)+
  scale_fill_manual(values = c("plex1_" = "#999999",
                               "plex2_" = "#E69F00",
                               "plex3_" = "#56B4E9"))
pca_plot

### Batch effect removal
# Check the order of your columns, POOL samples are not plotted
columns_checker(maxquant_clean_median_imp)

# Remove batch effect
## POOL CANNOT BE USED !!!
maxquant_clean_median_imp_unbatch <- remove_batch(dataset = maxquant_clean_median_imp, 
                                                  remove = "yes", use_pool = T)

##maxquant_clean_median_imp_unbatch <- remove_batch(dataset = maxquant_clean_median_imp, 
##                                                  remove = "yes",use_combat =  T, 
##                                                  where_is_the_batch1 = c(rep("batch1", times = 9),
##                                                                          rep("batch2", times = 9),
##                                                                          rep("batch3", times = 9)))

write.table(file = "./results/processed_maxquant_output.tsv", x = maxquant_clean_median_imp_unbatch, sep = "\t", dec = ".", row.names = F)

# Boxplot
intensity_boxplots_unbatch <- ggplot(maxquant_clean_median_imp_unbatch, mapping = aes(x = sample_name, y = unbatched_intensity))+
  geom_boxplot(data = maxquant_clean_median_imp_unbatch, mapping = aes(x = sample_name, y = unbatched_intensity, fill = plex))+
  #  geom_text_repel(data = long_format[which((long_format$high == "YES") &
  #                                          (long_format$is_max == "YES")),], aes(label = sample_name, y = max(intens+2), x = sample_name))+
  theme_bw()+
  scale_fill_manual(values = cbp1)+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection after batch correction")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

intensity_boxplots_unbatch  

# Median normalization
## median_all calculation
median_all2 <- median(maxquant_clean_median_imp_unbatch$unbatched_intensity)

## MEDIAN Calculation
MED_data <- maxquant_clean_median_imp_unbatch %>%
  subset(select = c(sample_name, unbatched_intensity)) %>% 
  group_by(sample_name) %>%
  summarise(MEDIAN = median(unbatched_intensity))

## merge the data
maxquant_clean_median_imp_unbatch <- merge(maxquant_clean_median_imp_unbatch, MED_data, 
                                           by = "sample_name")

## do the calculation 
maxquant_clean_median_imp_unbatch$unbatched_intensity <- (maxquant_clean_median_imp_unbatch$unbatched_intensity - 
                                                            maxquant_clean_median_imp_unbatch$MEDIAN) + median_all2 

# Boxplot
intensity_boxplots_unbatch_norm <- ggplot(maxquant_clean_median_imp_unbatch, mapping = aes(x = sample_name, y = unbatched_intensity))+
  geom_boxplot(data = maxquant_clean_median_imp_unbatch, mapping = aes(x = sample_name, y = unbatched_intensity, fill = plex))+
  #  geom_text_repel(data = long_format[which((long_format$high == "YES") &
  #                                          (long_format$is_max == "YES")),], aes(label = sample_name, y = max(intens+2), x = sample_name))+
  theme_bw()+
  scale_fill_manual(values = cbp1)+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("log2 Detection Intensity after batch correction normalization")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

intensity_boxplots_unbatch_norm  


### PCA 2
# Cahnge long to wide format
maxquant_clean_median_imp_to_pca2 <- reshape2::dcast(maxquant_clean_median_imp_unbatch, 
                                                     protein_group
                                                     ~ sample_name,value.var="unbatched_intensity",
                                                     fun.aggregate = median)

peak_mat2 <- t(maxquant_clean_median_imp_to_pca2[,c(2:ncol(maxquant_clean_median_imp_to_pca2))])

# PCA
pca2 <- prcomp(peak_mat2, scale. = TRUE, center = TRUE)

# PCA data construction
## positional data
pca2x <- as.data.frame(pca2$x)

## to_colour
samples_afer_batch <- intersect(rownames(peak_mat2), to_colour$sample_name)
meta_data_tracker <- subset(to_colour, sample_name %in% c(samples_afer_batch))

## construct the data for colouring
pca2x <- pca2x[meta_data_tracker$sample_name,]
pca2xx <- merge(pca2x, to_colour,by = "row.names")
rownames(pca2xx) <- pca2xx$Row.names
pca2xx <- pca2xx[,-1]

# PCA graph
pca_plot2 <- ggplot(pca2xx, 
                    aes(x = PC1, y = PC2, label = exp_group)) +
  
  geom_point(size = 4, aes(color = exp_group))+
  ###scale_color_manual(values = c("NCI_7_ref" = "red","pool4" = "green",
  ###                              "A549" = "yellow","CCRF_CEM" = "orange",
  ###                              "COLO_205" = "blue","H226" = "gray",
  ###                              "H23" = "black","KARPAS" = "pink",
  ###                              "RPM1_8226" = "brown","T47D" = "cyan"))+
  
  scale_color_manual(values = c("NCI_7_ref" = "red","pool4" = "green",
                                "A549" = "yellow","CCRF_CEM" = "orange",
                                "COLO_205" = "blue","H226" = "gray",
                                "H23" = "pink","RPM1_8226" = "brown","T47D" = "cyan"))+
  
  
  guides(color=guide_legend(title = "Experimental group"))+
  xlab(paste("PC1(", round(100*(pca2$sdev[1]^2/sum(pca2$sdev^2)), 2), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca2$sdev[2]^2/sum(pca2$sdev^2)), 2), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA) Un-batched data"))+
  geom_polygon(aes(group = exp_group, fill = exp_group), alpha = 0.2, show.legend = FALSE)+
  ###scale_fill_manual(values = c("NCI_7_ref" = "red","pool4" = "green",
  ###                                "A549" = "yellow","CCRF_CEM" = "orange",
  ###                            "COLO_205" = "blue","H226" = "gray",
  ###                            "H23" = "black","KARPAS" = "pink",
  ###                            "RPM1_8226" = "brown","T47D" = "cyan"))
  scale_fill_manual(values = c("NCI_7_ref" = "red","pool4" = "green",
                               "A549" = "yellow","CCRF_CEM" = "orange",
                               "COLO_205" = "blue","H226" = "gray",
                               "H23" = "pink","RPM1_8226" = "brown","T47D" = "cyan"))

pca_plot2


### UMAP
umap_res <- umap::umap(peak_mat2, n_neighbors = 3, n_components = 2, metric = "euclidean")
umap_df <- umap_res$layout
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df <- merge(umap_df, to_colour, by = "row.names")
rownames(umap_df) <- umap_df$Row.names
umap_df <- umap_df[,-1]

umap_plot <- ggplot(umap_df, 
                    aes(x = UMAP1, y = UMAP2, label = exp_group)) +
  geom_point(size = 4, aes(color = exp_group))+
  ###scale_color_manual(values = c("NCI_7_ref" = "red","pool4" = "green",
  ###                              "A549" = "yellow","CCRF_CEM" = "orange",
  ###                              "COLO_205" = "blue","H226" = "gray",
  ###                              "H23" = "black","KARPAS" = "pink",
  ###                              "RPM1_8226" = "brown","T47D" = "cyan"))+
  scale_color_manual(values = c("NCI_7_ref" = "red","pool4" = "green",
                                "A549" = "yellow","CCRF_CEM" = "orange",
                                "COLO_205" = "blue","H226" = "gray",
                                "H23" = "pink","RPM1_8226" = "brown","T47D" = "cyan"))+
  
  guides(color=guide_legend(title = "Experimental group"))+
  ggtitle(("UMAP"))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))+
  geom_polygon(aes(group = exp_group, fill = exp_group), alpha = 0.2, show.legend = FALSE)+
  ###scale_fill_manual(values = c("NCI_7_ref" = "red","pool4" = "green",
  ###                             "A549" = "yellow","CCRF_CEM" = "orange",
  ###                             "COLO_205" = "blue","H226" = "gray",
  ###                             "H23" = "black","KARPAS" = "pink",
  ###                             "RPM1_8226" = "brown","T47D" = "cyan"))
  scale_fill_manual(values = c("NCI_7_ref" = "red","pool4" = "green",
                               "A549" = "yellow","CCRF_CEM" = "orange",
                               "COLO_205" = "blue","H226" = "gray",
                               "H23" = "pink","RPM1_8226" = "brown","T47D" = "cyan"))

umap_plot


### Save plots 
ggsave("./plots/1_intensity_detection_raw.tiff", intensity_boxplots, width = 10, height = 10)
ggsave("./plots/2_completeness.tiff", plot = completeness_barplot, width = 10, height = 10)
ggsave("./plots/3_number_of_proteins.tiff", plot = number_of_prot, width = 10, height = 10)
ggsave("./plots/4_na_density.tiff", plot = na_density, width = 10, height = 10)
ggsave("./plots/5_intensity_detection_normalized.tiff", plot = intensity_boxplots_norm, width = 10, height = 10)
ggsave("./plots/6_pca1.tiff", plot = pca_plot, width = 10, height = 10)
ggsave("./plots/7_intensity_detection_unbatched.tiff", plot = intensity_boxplots_unbatch, width = 10, height = 10)
ggsave("./plots/8_intensity_detection_unbatched_normalized.tiff", plot = intensity_boxplots_unbatch_norm, width = 10, height = 10)
ggsave("./plots/9_pca2.tiff", plot = pca_plot2, width = 10, height = 10)
ggsave("./plots/10_umap.tiff", plot = umap_plot, width = 10, height = 10)

#### limma
# Enter data to limma
## Retrive expression matrix
expression_matrix <- as.data.frame((reshape2::dcast(maxquant_clean_median_imp_unbatch, 
                                                    protein_group ~ sample_name,value.var="unbatched_intensity", fun.aggregate = median)))

rownames(expression_matrix) <- expression_matrix[,1]
expression_matrix <- expression_matrix[,-1]
expression_matrix <- expression_matrix %>% 
  mutate_if(is.character, as.numeric)

## Create design and contrast matrix
meta_data_tracker <- subset(to_colour, sample_name %in% colnames(expression_matrix))
meta_data_tracker <- meta_data_tracker[match(colnames(expression_matrix), meta_data_tracker$sample_name),]
meta_data_tracker <- subset(meta_data_tracker, select = sample_name)
meta_data_tracker <- merge(meta_data_tracker, meta_data)
rownames(meta_data_tracker) <- meta_data_tracker$sample_name

## Expression matrix same order as meta_data_tracker
expression_matrix <- expression_matrix[,rownames(meta_data_tracker)]

## Check a correct order 
all(rownames(meta_data_tracker) == colnames(expression_matrix))

# Sample name correspondance
## meta_data
for (i in 1:length(rownames(meta_data_tracker))){
  meta_data_tracker$sample_name[i] <- paste0(meta_data_tracker$exp_group[i], "_",
                                             sub(".+_(.+)","\\1", meta_data_tracker$sample_number[i]))
}
meta_data_tracker$plex_name <- rownames(meta_data_tracker)
rownames(meta_data_tracker) <- meta_data_tracker$sample_name
all(meta_data_tracker$plex_name == colnames(expression_matrix)) 

## expression_matrix
colnames(expression_matrix) <- rownames(meta_data_tracker)

## Check a correct order 
all(rownames(meta_data_tracker) == colnames(expression_matrix))

# Export expression matrix
expression_matrix_to_export <- expression_matrix
expression_matrix_to_export <- tibble::rownames_to_column(expression_matrix_to_export, var = "Protein_group")
write_tsv(file = "./results/expression_matrix.tsv", x = expression_matrix_to_export)

## Create desing matrix
groups <- meta_data_tracker$exp_group
design <- model.matrix(~0 + groups)
colnames(design) <- gsub("^groups", "", colnames(design))
colnames(design) <- gsub(" ","_", colnames(design))
design

## Annotation
annotation <- subset(maxquant_clean, select = c(protein_group,protein_names,gene_names))
annotation <- distinct(annotation)
rownames(annotation) <- annotation$protein_group
annotation <- annotation[,-1]
annotation$Accession <- rownames(annotation)

## Fit the model
fit <- lmFit(expression_matrix, design = design)

# Contrast matrix
## Prepare all possile contrasts IF ALL CONTRASTS ARE THE ONES TO PERFORM
###contrasts_all <- make_all_contrasts(design)

# Reverse the desired contrasts
##contrasts_all <- upsidedown(contrasts_all, comparisons_to_change = c("Combo_R_vs_Ibrutinib_R",
##       "Combo_R_vs_ON123300_R",
##       "Control_R_vs_Ibrutinib_R",
##       "Control_R_vs_ON123300_R"))
contrasts_all <- makeContrasts("NCI_7_ref vs. rest" = NCI_7_ref  - (A549 + CCRF_CEM + COLO_205 + H226 + H23 + pool4 + RPM1_8226 + T47D)/8,  ## NCI_7_ref vs.rest
                               "NCI_7_ref vs. A549"  = NCI_7_ref - A549,
                               "NCI_7_ref vs. CCRF_CEM"  = NCI_7_ref - CCRF_CEM,
                               "NCI_7_ref vs. COLO_205"  = NCI_7_ref - COLO_205,
                               "NCI_7_ref vs. H226"  = NCI_7_ref - H226,
                               "NCI_7_ref vs. H23"  = NCI_7_ref - H23,
                               "NCI_7_ref vs. pool4"  = NCI_7_ref - pool4,
                               "NCI_7_ref vs. RPM1_8226"  = NCI_7_ref - RPM1_8226,
                               "NCI_7_ref vs. T47D"  = NCI_7_ref - T47D,
                               "A549 vs. COLO_205" = A549 - COLO_205,
                               levels = design)


# Fit contrasts to the model
fit1 <- contrasts.fit(fit = fit, contrasts = contrasts_all)
fit1 <- eBayes(fit = fit1)

### Extract limma data
tt <- topTable(fit1, number = Inf)
tt <- merge(tt, annotation, by = "row.names")
colnames(tt)[1] <- "protein_group"
tt <- relocate(.data = tt, c("protein_names","gene_names"), .after = "protein_group")
write.table(tt, "./results/TopTable_all.tsv", row.names = F, sep = "\t", dec = ".")

## write xlsx results
results_xlsx <- xlsx_tt(fit__1 = fit1, meta_data = meta_data_tracker, meta_sample_column = "sample_name", meta_data_column = "exp_group",
                        annotation = annotation, expression_matrix = expression_matrix, 
                        ###color_samples = c(
                        ###    "NCI_7_ref" = "red","pool4" = "green",
                        ###    "A549" = "yellow","CCRF_CEM" = "orange",
                        ###    "COLO_205" = "blue","H226" = "gray",
                        ###    "H23" = "black","KARPAS" = "pink",
                        ###    "RPM1_8226" = "brown","T47D" = "cyan"),
                        color_samples = c(
                          "NCI_7_ref" = "red","pool4" = "green",
                          "A549" = "yellow","CCRF_CEM" = "orange",
                          "COLO_205" = "blue","H226" = "gray",
                          "H23" = "pink","RPM1_8226" = "brown","T47D" = "cyan"),
                        differentation_element = " vs. ",
                        filename = "./results/TopTable_results_USA_ep.xlsx")


