library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

# Read the new_serf_df.csv data
serf_data <- read.csv("new_serf_df.csv", header = TRUE)

# Create a wide format matrix for diversity calculations
serf_wide <- serf_data %>%
  select(Sample_base, ASV, Abundance) %>%
  pivot_wider(names_from = ASV, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("Sample_base")

# Calculate Shannon and Simpson diversity indices
shannon_diversity <- diversity(serf_wide, index = "shannon")
simpson_diversity <- diversity(serf_wide, index = "simpson")
invsimpson_diversity <- diversity(serf_wide, index = "invsimpson")

# Create a dataframe with diversity indices and sample metadata
diversity_df <- data.frame(
  Sample_base = names(shannon_diversity),
  Shannon = shannon_diversity,
  Simpson = simpson_diversity,
  InvSimpson = invsimpson_diversity
) %>%
  # Extract treatment information from sample names
  mutate(
    # Split sample base to extract components
    OilType = case_when(
      grepl("^RAW_", Sample_base) ~ "Unweathered Diesel",
      grepl("^WTD1_", Sample_base) ~ "Weathered Diesel",
      grepl("^CRD_", Sample_base) ~ "Crude Oil",
      TRUE ~ "Unknown"
    ),
    # Extract nutrient/dispersant treatment
    NutDisp = case_when(
      grepl("_ND$", Sample_base) ~ "ND",
      grepl("_D$", Sample_base) ~ "D",
      grepl("_N$", Sample_base) ~ "N", 
      TRUE ~ "Unknown"
    ),
    # Extract concentration
    Conc_pct = case_when(
      grepl("_0.1_", Sample_base) ~ 0.1,
      grepl("_1_", Sample_base) ~ 1.0,
      grepl("_0_", Sample_base) ~ 0.0,
      TRUE ~ NA_real_
    ),
    # Create treatment group labels
    Treatment_Group = case_when(
      NutDisp == "ND" ~ "Nutri+Disp",
      NutDisp == "D" ~ "Dispersants", 
      NutDisp == "N" ~ "Nutrients",
      TRUE ~ as.character(NutDisp)
    ),
    # Create concentration labels
    Conc_Label = case_when(
      is.na(Conc_pct) ~ "Unknown",
      TRUE ~ paste0(Conc_pct, "%")
    )
  )

diversity_df_clean <- diversity_df %>%
  filter(OilType != "Unknown" & Treatment_Group != "Unknown")


# Print summary of diversity data

cat("Shannon diversity range:", round(min(diversity_df$Shannon), 3), "to", round(max(diversity_df$Shannon), 3), "\n")
cat("Simpson diversity range:", round(min(diversity_df$Simpson), 3), "to", round(max(diversity_df$Simpson), 3), "\n")

print(head(diversity_df))

# Create violin plots for Shannon diversity
p_shannon <- ggplot(diversity_df_clean, aes(x = Treatment_Group, y = Shannon, fill = Treatment_Group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  facet_wrap(~ OilType, scales = "free_x") +
  scale_fill_brewer(palette="Purples")+
  labs(
    x = "Treatment", 
    y = "Shannon Diversity Index"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(family= "Arial", size = 15, face= "bold"),
    axis.text.x = element_text(family= "Arial", hjust = 0.5, size = 15),
    strip.text = element_text(size = 15, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(family= "Arial", size= 14),
    legend.title = element_text(family= "Arial", size=14)
    #plot.title = element_text(size = 14, face = "bold")
    #plot.subtitle = element_text(size = 11, style = "italic")
  )

print(p_shannon)

ggsave("p_shannon.png", plot = p_shannon, 
       width = 16, height = 8, dpi = 300)


# Create violin plots for Simpson diversity
p_simpson <- ggplot(diversity_df_clean, aes(x = Treatment_Group, y = Simpson, fill = Treatment_Group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  facet_wrap(~ OilType, scales = "free_x") +
  scale_fill_brewer(palette="Purples")+
  labs(
    x = "Treatment", 
    y = "Simpson Diversity Index"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(family= "Arial", size = 15, face= "bold"),
    axis.text.x = element_text(family= "Arial", hjust = 0.5, size = 15),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(family= "Arial", size= 14),
    legend.title = element_text(family= "Arial", size=14)
    #plot.title = element_text(size = 14, face = "bold"),
    #plot.subtitle = element_text(size = 11, style = "italic")
  )

print(p_simpson)

ggsave("p_simpson.png", plot = p_simpson, 
       width = 16, height = 8, dpi = 300)


# Combined plot showing both diversity indices
diversity_long <- diversity_df_clean %>%
  select(Sample_base, OilType, Treatment_Group, Conc_Label, Shannon, Simpson) %>%
  pivot_longer(cols = c(Shannon, Simpson), names_to = "Diversity_Index", values_to = "Value")

p_combined <- ggplot(diversity_long, aes(x = Treatment_Group, y = Value, fill = Treatment_Group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  facet_grid(Diversity_Index ~ OilType, scales = "free") +
  scale_fill_brewer(palette="Purples")+
  labs(
    x = "Treatment", 
    y = "Diversity Index Value"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(family= "Arial", size= 15, face= "bold"),
    axis.text.x = element_text(family= "Arial", hjust = 0.5, size = 15),
    strip.text = element_text(size = 15, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(family= "Arial", size= 14),
    legend.title = element_text(family= "Arial", size=14)
    #plot.title = element_text(size = 14, face = "bold")
  )

print(p_combined)

ggsave("p_combined.png", plot = p_combined, 
       width = 16, height = 8, dpi = 300)
...............................................................................................................

serf_wide_nmds <- serf_wide[rowSums(serf_wide) > 0, ]

cat("NMDS input data dimensions:", nrow(serf_wide_nmds), "samples x", ncol(serf_wide_nmds), "ASVs\n")

# Perform NMDS ordination using Bray-Curtis dissimilarity
set.seed(123)  # For reproducible results
nmds_result <- metaMDS(serf_wide_nmds, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)

# Extract NMDS scores
nmds_scores <- as.data.frame(scores(nmds_result, display="sites"))
nmds_scores$Sample_base <- rownames(nmds_scores)

nrow(scores(nmds_result))   


# Add sample metadata to NMDS scores
nmds_data <- nmds_scores %>%
  mutate(
    # Extract treatment information from sample names
    OilType = case_when(
      grepl("^RAW_", Sample_base) ~ "RAW",
      grepl("^WTD1_", Sample_base) ~ "WTD1",
      grepl("^CRD_", Sample_base) ~ "CRD",
      grepl("^PRE_", Sample_base) ~ "PRE",
      grepl("^CTL_", Sample_base) ~ "CTL",
      TRUE ~ "Unknown"
    ),
    NutDisp = case_when(
      grepl("_ND$", Sample_base) ~ "ND",
      grepl("_D$", Sample_base) ~ "D",
      grepl("_N$", Sample_base) ~ "N",
      grepl("^PRE_", Sample_base) ~ "PRE",
      grepl("^CTL_.*_$", Sample_base) ~ "CTL", 
      TRUE ~ "Unknown"
    ),
    Conc_pct = case_when(
      grepl("_0.1_", Sample_base) ~ 0.1,
      grepl("_1_", Sample_base) ~ 1.0,
      grepl("^PRE_0_", Sample_base) ~ 0.0,
      grepl("^CTL_0.1_", Sample_base) ~ 0.1,
      TRUE ~ NA_real_
    ),
    Treatment_Group = case_when(
      NutDisp == "ND" ~ "Nutrients+Dispersants",
      NutDisp == "D" ~ "Dispersant Only", 
      NutDisp == "N" ~ "Nutrients Only",
      NutDisp == "PRE" ~ "Pre-treatment",
      NutDisp == "CTL" ~ "Control",
      TRUE ~ as.character(NutDisp)
    ),
    Conc_Label = case_when(
      is.na(Conc_pct) ~ "Unknown",
      TRUE ~ paste0(Conc_pct, "%")
    )
  ) %>%
  filter(OilType != "Unknown" & Treatment_Group != "Unknown") %>%
  mutate(
    OilType = factor(OilType, levels = c("RAW", "WTD1", "CRD", "PRE", "CTL")),
    Treatment_Group = factor(Treatment_Group)
  )

unique(nmds_data)

nmds_scores$Sample_base[!(nmds_scores$Sample_base %in% serf_data$Sample_base)]
head(nmds_scores)

# Stress interpretation
#stress_interpretation <- case_when(
  #nmds_result$stress < 0.05 ~ "Excellent representation",
  #nmds_result$stress < 0.1 ~ "Good representation", 
  #nmds_result$stress < 0.2 ~ "Fair representation",
  #nmds_result$stress < 0.3 ~ "Poor representation",
  #TRUE ~ "Very poor representation"
#)
#cat("Stress interpretation:", stress_interpretation, "\n")


# Create NMDS plot colored by Treatment Group
p_nmds_treatment <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, color = Treatment_Group, shape = OilType)) +
  stat_ellipse(aes(color = Treatment_Group, group = Treatment_Group), geom = "polygon", alpha = 0, level= 0.95,type = "norm", linetype = 1) +
  
  geom_point(size = 4, alpha = 1) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8), 
                     labels = c("Raw Diesel", "Weathered Diesel", "Crude Oil", "Pre-treatment", "Control"),
                     name = "Oil Type") +
  scale_color_brewer(palette = "Set2", name = "Treatment Group") +
  scale_fill_brewer(palette = "Set2", name = "Treatment Group") +
  labs(
    subtitle = paste0("Stress = ", round(nmds_result$stress, 3)),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size=15),
    axis.text = element_text(size=15),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size=10),
    legend.title = element_text(size=15)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    shape = guide_legend(override.aes = list(size = 3))
  )
 print(p_nmds_treatment)
 
 ggsave("p_nmds_treatment.png", plot = p_nmds_treatment, 
        width = 10, height = 8, dpi = 300)

.............................................................................................................................
 #PERMANOVA test between oil types
 
 str(seqtab.serf.new)
 seqtab.serf.numeric <- seqtab.serf.new %>%
   select(where(is.numeric))
 #seqtab.numeric <- seqtab.numeric[match(meta.filtered$SampleID, rownames(seqtab.numeric)), ]
 #all(rownames(seqtab.numeric) == meta.filtered$SampleID)
 permanova_serf<- adonis2(seqtab.serf.numeric ~ OilType, data = meta.serf.new, method = "bray", permutations = 999)
 print(permanova_serf)
 
 
 #PERMANOVA test between treatmnet types 
 
 permanova_serf <- adonis2(seqtab.serf.numeric ~ NutDisp,
                      data = meta.serf.new,
                      method = "bray",
                      permutations = 999)
 
 print(permanova_serf)
 
 # ANOSIM test for oil types 
 ano_s <- anosim(seqtab.serf.numeric, grouping = meta.serf.new$OilType, distance = "bray")
 print(ano_s)
 
 
 #counting taxonmoic levels in the dataset
 
 length(unique(serf_data$ASV))
 length(unique(serf_data$Genus))
 length(unique(serf_data$Phylum))
 length(unique(serf_data$Family))
 length(unique(serf_data$Class))
 length(unique(serf_data$Order))
 
 #See if the CMO and SERF data are different or identical
 table(serf_data$Genus == working_data_avg$Genus)
 identical(serf_data, working_data_avg)
 
 #Counting  the taxonomic levels - option 2
 count_unique <- function(df) {
   sapply(df[, c("ASV","Genus","Phylum","Family","Class","Order")], function(x) length(unique(x)))
 }
 
 count_unique(serf_data)
 count_unique(working_data_avg)
 
 
 
 meta.serf.new$Treatment <- with(meta.serf.new, paste(Conc_pct, NutDisp, sep = "_"))
 adonis2(seqtab.serf.numeric ~ Treatment,
         data = meta.serf.new,
         method = "bray")
 
 table(meta.serf.new$Conc_pct, meta.serf.new$NutDisp)
 
 ado <- adonis2(seqtab.serf.numeric ~ Conc_pct + NutDisp,
                data = meta.serf.new,
                method = "bray")
ado 

