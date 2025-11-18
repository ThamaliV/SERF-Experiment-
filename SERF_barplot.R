library(gplots)
library(dendextend)
library(Rtsne)
library(vegan)
library(fastcluster)
library(ggplot2)
library(tidyverse)
library(scales)
library(viridis)
library(tibble)
library(tidyr)
library(forcats)
library(dplyr)



seqtab <- readRDS("E:/Thami-Uni/Oil Microcosm Data/seqtab_final.rds")

rownames(seqtab) = gsub(pattern = "716S",replacement = "7-16S",x = rownames(seqtab))

  
  
  
  
  sepnames.16S = separate(data=data.frame(x=rownames(seqtab)),col=x, into=c("project","primer","plate","well"),sep="[|]", fill="right")
  sepnames.16S$plate.geno = toupper(paste(sepnames.16S$platenum$p2,sepnames.16S$well,sep="_"))
  sepnames.16S$platenum = separate(data=sepnames.16S,col=plate, into=c("p1","p2","p3","p4","p5","p6"), sep="-", fill="right")
  sepnames.16S$plate.geno = toupper(paste(sepnames.16S$platenum$p2,sepnames.16S$well,sep="_"))
  seqtab.new = seqtab[grep(x=rownames(seqtab),"20250811"),]
  
  head(rownames(seqtab.new))
  seqtab.new = seqtab.new[,colSums(seqtab.new)>0]
  rownames(seqtab.new) = sepnames.16S$platenum$well[grep("20250811",rownames(seqtab))]
  
  meta.in = read.csv("E:/Thami-Uni/SERF_MOS analysis/Genohub IDs Project 2599682, 3694692.csv", 
                     header = TRUE)
  rownames(meta.in) = meta.in$uniq.id
  meta.in = meta.in[rownames(seqtab.new),]
  
  taxout.16S = as.data.frame(readRDS("E:/Thami-Uni/Oil Microcosm Data/taxout_edit.rds")$`ref_dada2_silva_v138-2_train_set_uniq.fasta.gz`)
  taxout.16S$ASV = paste0("bASV",1:nrow(taxout.16S))
  taxout.16S$short = paste(taxout.16S$ASV,taxout.16S$Domain, taxout.16S$Order, taxout.16S$Genus)
  
  taxout.all = taxout.16S
  seqtab.new.small.serf = seqtab.new[,colMeans(seqtab.new)>10]
  
  dim(seqtab.new.small.serf)
  head(seqtab.new.small.serf)
  taxnames = taxout.all[colnames(seqtab.new.small.serf),"short"]
  
  
  rn_serf <- rownames(seqtab.new.small.serf)
  
  # get the index positions of start and end wells
  start <- which(rn_serf == "plate2-C10")
  end   <- which(rn_serf == "plate2-E05")
  
  
  
  # subset rows between them
  
  seqtab.serf <- seqtab.new.small.serf[start:end, ]   # subset rows first
  
  extra <- seqtab.new.small.serf[which(rn_serf == "plate2-F02"), , drop = FALSE]
  
  seqtab.serf.new <- rbind(seqtab.serf, extra)
  
  
  seqtab.serf.new <- as.data.frame(seqtab.serf.new)
  
  seqtab.serf.new <- cbind("uniq.id" = rownames(seqtab.serf.new), seqtab.serf.new)
  rownames(seqtab.serf.new) <- NULL 
  
 
  meta.serf <- meta.in[start:end, ]
  
  # Add the metadata row for plate2-F02
  extra.meta <- meta.in[rownames(meta.in) == "plate2-F02", , drop = FALSE]
  
  meta.serf.new <- rbind(meta.serf, extra.meta)
  
  
  
  seqtab_serf_long <- seqtab.serf.new %>%
    pivot_longer(
      cols = -uniq.id,             # all columns except SampleID
      names_to = "Sequence",             # new column for ASV IDs
      values_to = "Abundance" )  # new column for counts  
  
  
  merged_serf_seqtab <- seqtab_serf_long %>%
    left_join(meta.in, by = "uniq.id")
  
  taxout.all <- rownames_to_column(taxout.all, var = "Sequence")
  
  working_serf_data <- merged_serf_seqtab %>%
    left_join(taxout.all, by = "Sequence") 
  
  working_serf_data <- working_serf_data %>%
    mutate(Sample_base = paste(OilType,Conc_pct,NutDisp, sep = "_"))
  
  new_serf_df <- working_serf_data[, c("Sample_base", "Site", "Time", "NutDisp", "OilType", "Conc_pct","ASV", "Domain","Phylum","Class","Order","Family","Genus", "Abundance")]
  
  
  #write.csv(new_serf_df, "new_serf_df.csv", row.names = FALSE)
  
  
  genus_serf <- new_serf_df %>%
    group_by(Sample_base) %>%
    mutate(RelAbundance = Abundance / sum(Abundance),
        
   Genus = ifelse(RelAbundance < 0.01, "Other", Genus))%>%
    ungroup()
  
  
  #write.csv(genus_serf, "genus_serf.csv", row.names = FALSE)
  
  my_colors <- c(
    
    "#003f5c",
    "#2f4b7c",
    "#665191",
    "#a05195",
    "#d45087",
    "#7b3c66",
    "#b94c5b",
    "#da743b",
    "#d0ae14",
    "#cc5c8f",
    "#565983",
    "#7a5d93",
    "#a35d97",
    "#ed607c",
    "#ff6e60",
    "#ff873e",
    "#ffa600",
    "#718380",
    "#669081",
    "#005b76",
    "#00788a",
    "#009694",
    "#22b397",
    "#68cf92",
    "#a7e98c",
    "#eaff89",
    "#393540"
    
    
  )
  set.seed(123)  
  genus_colors <- sample(my_colors)
  
  genera <- unique(genus_serf$Genus)
  print(genera)
  
  
  diesel_df <- genus_serf %>%
    filter(OilType %in% c("PRE","CTL","RAW","WTD1")) %>%
    group_by(Sample_base, Site, Time, NutDisp, OilType, Conc_pct, Genus) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Sample_base) %>%
    mutate(RelAbundance = Abundance / sum(Abundance)) %>%
    ungroup()
  
  
  diesel <- unique(diesel_df$OilType)
  print(diesel)
  
  gen <- unique(genus_serf$OilType)
  print(gen)
  
  crude_df <- genus_serf %>%
    filter(OilType %in% c("PRE","CRD")) %>%
    group_by(Sample_base, Site, Time, NutDisp, OilType, Conc_pct, Genus) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Sample_base) %>%
    mutate(RelAbundance = Abundance / sum(Abundance)) %>%
    ungroup()
  
  
  crude <- unique(crude_df$OilType)
  crude
  

  
  # Create better labels for diesel oil types and treatment groups
  diesel_df <- diesel_df %>%
    mutate(
      OilType_Label = case_when(
        OilType == "PRE" ~ "Pre",
        OilType == "CTL" ~ "Control", 
        OilType == "RAW" ~ "Raw Diesel",
        OilType == "WTD1" ~ "Weathered Diesel",
        TRUE ~ as.character(OilType)
      ),
      # Create clean treatment labels for X-axis
      Treatment_Group = case_when(
        NutDisp == "ND" ~ "Nutri+Disp",
        NutDisp == "D" ~ "Dispersants", 
        NutDisp == "N" ~ "Nutrients",
        is.na(NutDisp) | NutDisp == "" ~ "Control",
        TRUE ~ as.character(NutDisp)
      ),
      # Create concentration labels
      Conc_Label = case_when(
        is.na(Conc_pct) ~ "",
        TRUE ~ paste0(Conc_pct, "%")
      ),
      # Create clean X-axis labels combining treatment and concentration
      Clean_X_Label = case_when(
        Conc_Label == "" ~ Treatment_Group,
        TRUE ~ paste(Treatment_Group, Conc_Label, sep = "\n")
      )
    )
  
  diesel_df$OilType <- factor(
    diesel_df$OilType,
    levels = c("PRE", "CTL", "RAW", "WTD1")
  )
  
  # Factor the labels in the same order
  diesel_df$OilType_Label <- factor(
    diesel_df$OilType_Label,
    levels = c("Pre", "Control", "Raw Diesel", "Weathered Diesel")
  )
  
  # Reorder Genus factor to put "Other" at the bottom (last level)
  diesel_df$Genus <- factor(diesel_df$Genus, 
                            levels = c(setdiff(unique(diesel_df$Genus), "Other"), "Other"))
  
  
  # Create a combined grouping variable
  #genus_serf$Group <- paste(genus_serf$OilType, genus_serf$NutDisp, genus_serf$Conc_pct, sep = "_")
  
  p_diesel <- ggplot(diesel_df, aes(x = Clean_X_Label, y = RelAbundance, fill = Genus)) +
    #geom_bar(stat = "identity", position = "stack", width = 0.8) +
    geom_bar(stat = "identity", linewidth = 0.1, colour= "black", linewidth = 0.2, position = "stack" ) +
    facet_grid(~OilType_Label, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = my_colors) +
    labs(x = "Treatment", y = "Relative Abundance", fill = "Genus") +
    theme_bw() +
    theme(
      axis.title = element_text(family= "Arial", size = 25, face= "bold"),
      axis.text.x = element_text(family = "Arial",hjust = 0.5, size = 13),
      strip.text = element_text(family = "Arial", size = 22, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size= 15),
      legend.title = element_text(size=15),
      panel.spacing = unit(0.3, "lines")
    )
print(p_diesel)  
  

taxa_colors <- c(
    "Colwellia"= "#003f5c",
    "Erythrobacter"= "#2f4b7c",
    "Hyphomonas"=   "#665191",
    "Incertae Sedis"=  "#a05195",
    "Legionella"= "#d45087",
    "Lentimonas"= "#7b3c66",
    "Marinobacter"=  "#ff6e60",
    "Marinomonas"=  "#ffa600",
    "Marivita"=   "#b94c5b",
    "Methylophaga"=  "#da743b",
    "Parasedimentitalea"=  "#d0ae14",
    "Pseudoalteromonas"=  "#718380",
    "Pseudohongiella"= "#cc5c8f",
    "Salinibacterium"=  "#565983",
    "Seohaeicola"= "#7a5d93",
    "Tropicimonas-Psychromarinibacter" =  "#a35d97",
    "Yoonia"=  "#ed607c",
    "Other"=   "#393540")
    

crude_df <- crude_df %>%
  mutate(
    Treatment_Group = case_when(
      NutDisp == "ND" ~ "Nutri+Disp",
      NutDisp == "D" ~ "Dispersant", 
      NutDisp == "N" ~ "Nutrients",
      NutDisp == "PRE" ~ "Pre-treatment",
      is.na(NutDisp) | NutDisp == "" ~ "Control",
      TRUE ~ as.character(NutDisp)
    ),
    # Rename OilType values for better panel labels
    OilType_Label = case_when(
      OilType == "PRE" ~ "Pre", 
      OilType == "CRD" ~ "Weathered Crude Oil  ",
      is.na(OilType) | OilType == "" ~ "Control",
      TRUE ~ as.character(OilType)
    ),
    # Create concentration labels
    Conc_Label = case_when(
      is.na(Conc_pct) ~ "",
      TRUE ~ paste0(Conc_pct, "%")
    ),
    # Combine treatment and concentration for detailed labeling
    Detailed_Treatment = case_when(
      Conc_Label == "" ~ Treatment_Group,
      TRUE ~ paste(Treatment_Group, Conc_Label, sep = "\n")
    )
  )


# Reorder Genus factor to put "Other" at the bottom (last level)
crude_df$Genus <- factor(crude_df$Genus, 
                          levels = c(setdiff(unique(crude_df$Genus), "Other"), "Other"))


p_crude <- ggplot(crude_df, aes(x = Detailed_Treatment, y = RelAbundance, fill = Genus)) +
  
  geom_bar(stat = "identity", linewidth = 0.1, colour= "black", linewidth = 0.2, position = "stack" ) +
  facet_grid(~OilType_Label, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = taxa_colors) +
  labs(x = "Treatment", y = "Relative Abundance", fill = "Genus") +
  theme_bw() +
  theme(
    axis.title = element_text(family= "Arial", size = 25, face= "bold"),
    axis.text.x = element_text(family= "Arial", hjust = 0.5, size = 13),
    strip.text = element_text(family= "Arial", size = 22, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size= 15),
    legend.title = element_text(size=15),
    panel.spacing = unit(0.3, "lines")
  )
    
print(p_crude)

ggsave("p_crude.png", plot = p_crude, 
       width = 16, height = 8, dpi = 300)

ggsave("p_diesel.png", plot = p_diesel,
       width =16, height =8, dpi=300)




........................................................................................................


# creating the bubble plot

genus_serf <- new_serf_df %>%
  group_by(Sample_base) %>%
  mutate(RelAbundance = Abundance / sum(Abundance))%>%
           ungroup()
# Sum abundance across all samples for each genus
top_genera <- genus_serf %>%
  group_by(Genus) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  #slice_head(n =30) %>%
  pull(Genus)  # Get the genus names

# Filter the main data
bubble_df <- genus_serf %>%
  filter(Genus %in% top_genera, OilType %in% c("CRD","WTD1"))



# Determine which genera are abundant in each OilType
genus_presence <- bubble_df %>%
  group_by(Genus, OilType) %>%
  summarise(max_rel = max(RelAbundance), .groups = "drop") %>%
  pivot_wider(names_from = OilType, values_from = max_rel, values_fill = 0) %>%
  mutate(
    in_diesel = if_else(WTD1 > 0.05, TRUE, FALSE),
    in_crude  = if_else(CRD > 0.05, TRUE, FALSE),
    Group = case_when(
      in_diesel & in_crude ~ "Both",
      in_diesel & !in_crude ~ "Diesel only",
      !in_diesel & in_crude ~ "Crude only",
      TRUE ~ "None"
    )
  )

# Check results
table(genus_presence$Group)
head(genus_presence)

#Join group info and filter
bubble_df <- bubble_df %>%
  left_join(genus_presence %>% select(Genus, Group), by = "Genus") %>%
  filter(Group != "None")

# Define factor levels for group order
bubble_df <- bubble_df %>%
  mutate(Group = factor(Group, levels = c("Diesel only", "Crude only", "Both")))

# Create ordered genus list by group
ordered_genera <- bubble_df %>%
  distinct(Genus, Group) %>%
  arrange(Group) %>%
  pull(Genus)

# Apply that ordering to the Genus factor
bubble_df <- bubble_df %>%
  mutate(Genus = factor(Genus, levels = ordered_genera))


bubble_df$Treatment_Conc <- paste(bubble_df$NutDisp, bubble_df$Conc_pct, sep = "_")
bubble_df$Treatment_Conc <- factor(bubble_df$Treatment_Conc)


facet_labels <- c(
  "CRD"  = "Weathered Crude Oil",
  "WTD1" = "Weathered Diesel"
)


group_positions <- bubble_df %>%
  distinct(Genus, Group) %>%
  group_by(Group) %>%
  summarise(ymax = max(as.numeric(Genus)), .groups = "drop") %>%
  arrange(ymax)



# Create the bubble plot without facets, single Y-axis
bubble_plot3 <- ggplot(bubble_df, aes(x = Sample_base, y = Genus)) +
 
  geom_point(aes(size = RelAbundance, color = Treatment_Conc), 
             alpha = 0.6) +
  facet_wrap(~ OilType, scales = "free_x", labeller = labeller(OilType = facet_labels)) +
  scale_size_continuous(range = c(1, 20), 
                        name = "Relative Abundance") +
  scale_color_brewer(palette = "Set2", name = "Treatment Group",
                     labels = c(
      
                       "N_0.1"   = "Nutrients 0.1%",
                       "N_1"   = "Nutrients 1%",
                       "ND_0.1"  = "Nutri + Disp 0.1%",
                       "ND_1" = "Nutri+Disp 1%",
                       "D_0.1"= "Dispersant 0.1%",
                       "D_1" = "Dispersant 1%"
                     )) +
  labs(
    
    x = "Time",
    y = "Genus") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 50,size=20, face= "bold"),
    strip.text = element_text(size = 20, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size= 12),
    legend.title = element_text(size= 15, face="bold"),
        legend.spacing.x = unit(1, "cm"),
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(3, "lines"),
    legend.title.align = 0.5, # center-align title above items
      
    
  )+
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5,   # center the title
      override.aes = list(size = 7)
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,   # center the title
      override.aes = list(shape = 21, fill = "white", color = "black", stroke = 1)
    )
  )


  for (pos in group_positions$ymax[-nrow(group_positions)]) {  # skip last group
    bubble_plot3 <- bubble_plot3 +
      geom_hline(yintercept = pos + 0.50, linetype = "dashed", color = "grey50")
  }



print(bubble_plot3)

#ggsave("crd v dsl.png", plot = bubble_plot, 
 #      width = 15, height = 7, dpi = 300)

ggsave("crdvdsl.png", plot = bubble_plot3, 
       width = 15, height = 7, dpi = 450)

dim(genus_presence)
head(genus_presence)




#bubbble plot of WTD vs Raw 

bubble_df2 <- genus_serf %>%
  filter(Genus %in% top_genera, OilType %in% c("RAW","WTD1"))



# Sum abundance across all samples for each genus
top_genera2 <- genus_serf %>%
  group_by(Genus) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  #slice_head(n =30) %>%
  pull(Genus)  # Get the genus names

# Determine which genera are abundant (>0.15) in each OilType
genus_presence2 <- bubble_df2 %>%
  group_by(Genus, OilType) %>%
  summarise(max_rel = max(RelAbundance), .groups = "drop") %>%
  pivot_wider(names_from = OilType, values_from = max_rel, values_fill = 0) %>%
  mutate(
    in_wtd = if_else(WTD1 > 0.05, TRUE, FALSE),
    in_raw  = if_else(RAW > 0.05, TRUE, FALSE),
    Group = case_when(
      in_wtd & in_raw ~ "Both",
      in_wtd & !in_raw ~ "Weathered only",
      !in_wtd & in_raw ~ "Raw only",
      TRUE ~ "None"
    )
  )

# Check results
table(genus_presence2$Group)
head(genus_presence2)

# Join group info and filter
bubble_df2 <- bubble_df2%>%
  left_join(genus_presence2 %>% select(Genus, Group), by = "Genus") %>%
  filter(Group != "None")

# Define factor levels for group order
bubble_df2 <- bubble_df2 %>%
  mutate(Group = factor(Group, levels = c("Weathered only", "Raw only", "Both")))

# Create ordered genus list by group
ordered_genera2 <- bubble_df2 %>%
  distinct(Genus, Group) %>%
  arrange(Group) %>%
  pull(Genus)

# Apply that ordering to the Genus factor
bubble_df2 <- bubble_df2 %>%
  mutate(Genus = factor(Genus, levels = ordered_genera2))



facet_labels <- c(
  "RAW"  = "Raw Diesel",
  "WTD1" = "Weathered Diesel"
)


group_positions2 <- bubble_df2 %>%
  distinct(Genus, Group) %>%
  group_by(Group) %>%
  summarise(ymax = max(as.numeric(Genus)), .groups = "drop") %>%
  arrange(ymax)




bubble_df2$Treatment_Conc2 <- paste(bubble_df2$NutDisp, bubble_df2$Conc_pct, sep = "_")
bubble_df2$Treatment_Conc2 <- factor(bubble_df2$Treatment_Conc2)






# Create the bubble plot without facets, single Y-axis
bubble_plot2 <- ggplot(bubble_df2, aes(x = Sample_base, y = Genus)) +
  
  geom_point(aes(size = RelAbundance, color = Treatment_Conc2), 
             alpha = 0.6) +
  facet_wrap(~ OilType, scales = "free_x", labeller = labeller(OilType = facet_labels)) +
  scale_size_continuous(range = c(1, 20), 
                        name = "Relative Abundance") +
  scale_color_brewer(palette = "Set2", name = "Treatment Group",
                     labels = c(
                       
                       "N_0.1"   = "Nutrients 0.1%",
                       "N_1"   = "Nutrients 1%",
                       "ND_0.1"  = "Nutri + Disp 0.1%",
                       "ND_1" = "Nutri+Disp 1%",
                       "D_0.1"= "Dispersant 0.1%",
                       "D_1" = "Dispersant 1%"
                     )) +
  labs(
    
    y = "Genus") +
   
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 50, size=20, face= "bold"),
    strip.text = element_text(size = 20, face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size= 12),
    legend.title = element_text(size= 15, face="bold"),
    legend.spacing.x = unit(1, "cm"),
    legend.key.width = unit(2, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.title.align = 0.5, # center-align title above items
    
    
  )+
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5,  # center the title
      order= 2,
      override.aes = list(size = 7)
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5, # center the title
      order =1, 
      override.aes = list(shape = 21, fill = "white", color = "black", stroke = 1)
    )
  )

for (pos in group_positions2$ymax[-nrow(group_positions2)]) {  # skip last group
  bubble_plot2 <- bubble_plot2 +
    geom_hline(yintercept = pos +0.5, linetype = "dashed", color = "grey50")
}

  


print(bubble_plot2)

ggsave("Raw v wtd.png", plot = bubble_plot2, 
       width = 16, height = 7, dpi = 300)


table(bubble_df$Sample_base, bubble_df$Conc_pct)





