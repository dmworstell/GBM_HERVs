library(tidyverse)
library(ggpubr)
library(dplyr)


calculate_p_value <- function(df1_PDUI, df2_PDUI) {
  result <- tryCatch(
    wilcox.test(df1_PDUI, df2_PDUI)$p.value,
    error = function(e) {
      print("Error in wilcox.test:")
      print(e)
      print("Data passed to wilcox.test:")
      print(df1_PDUI)
      print(df2_PDUI)
      return(NA)
    }
  )
  return(result)
}

# Read in the data
#HERVs_Dapars_OL <- read_tsv("/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/alternative_polyA/DaPars2_results/output/simplified_sorted_ALLprvs.tsv")
HERVs_Dapars_OL = read_tsv("C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/alternative_polyA/DaPars2_results/output/simplified_sorted_ALLprvs.tsv")

colnames(HERVs_Dapars_OL)[colnames(HERVs_Dapars_OL) == 'Provirus_chr'] <- 'Chromosome'
HERVs_Dapars_OL$Overlaps_HERV = TRUE

#HERVs_Dapars_NO_OL <- read_tsv("/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/alternative_polyA/DaPars2_results/output/simplified_sorted_NOprvs.tsv")
HERVs_Dapars_NO_OL <- read_tsv("C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/alternative_polyA/DaPars2_results/output/simplified_sorted_NOprvs.tsv")

colnames(HERVs_Dapars_NO_OL)[colnames(HERVs_Dapars_NO_OL) == 'DaPars_chr'] <- 'Chromosome'
HERVs_Dapars_NO_OL$Overlaps_HERV = FALSE

data <- bind_rows(HERVs_Dapars_OL, HERVs_Dapars_NO_OL)

# Add a new column to indicate if a gene has provirus overlaps in multiple features
data <- data %>%
  group_by(Gene_name) %>%
  mutate(Multiple_Features = ( n_distinct(Provirus_ID) > 1  | n_distinct(Overlapping_exon_number) > 1 )) %>%
  ungroup()



# Define the filter conditions
conditions <- list(
  "Overlap with 3' terminal exon vs no HERV overlap" = "(Overlaps_HERV == FALSE) | Terminal_exon_overlap == 'last_exon'",
#  "Overlap with 5' terminal exon vs no HERV overlap" = "(Overlaps_HERV == FALSE) | Terminal_exon_overlap == 'first_exon'",
#  "Overlap with single exon vs no HERV overlap" = "(Overlaps_HERV == FALSE) | Terminal_exon_overlap == 'single_exon'",
  "Overlap with other exon vs no HERV overlap" =  "(Overlaps_HERV == FALSE) | Terminal_exon_overlap == 'middle_exon'",
  "Overlap with intron vs no HERV overlap" =  "(Overlaps_HERV == FALSE) | Terminal_exon_overlap == 'intron'"
  )


# Define gene types
gene_types <- c("ProteinCoding", "lncRNA")


# Loop through the gene types and generate a plot for each
#for (gene_type in gene_types) {
filtered_data_list <- list()

for (name in names(conditions)) {
  filtered_data <- data %>%
    dplyr::filter(eval(parse(text = conditions[[name]]))) %>%
    dplyr::group_by(Gene_name, Overlaps_HERV) %>%
    #     filter(PC_or_lncRNA == gene_type) %>%
    #      dplyr::group_by(Gene_name, Overlaps_HERV, PC_or_lncRNA) %>%
    dplyr::summarise(Control_median_PDUI = median(Control_median_PDUI,na.rm=TRUE),
              Disease_median_PDUI = median(Disease_median_PDUI,na.rm=TRUE),
              PDUI_difference = median(Control_median_PDUI,na.rm=TRUE) - median(Disease_median_PDUI,na.rm=TRUE),
              .groups = "drop") %>% dplyr::ungroup()
  
  filtered_data$Condition <- name
  filtered_data_list[[name]] <- filtered_data
}

# Combine all of the filtered data into one data frame
all_filtered_data <- bind_rows(filtered_data_list)

# Remove rows with NA or Inf values in PDUI_difference
#all_filtered_data <- all_filtered_data %>%
#  filter(!is.na(PDUI_difference) & PDUI_difference != Inf & PDUI_difference != -Inf)

# Remove conditions with fewer than 2 genes
all_filtered_data <- all_filtered_data %>% 
  dplyr::group_by(Condition) %>% 
  #dplyr::group_by(Condition, PC_or_lncRNA) %>% 
  dplyr::filter(n() > 1)

# For each group (defined by Condition and PC_or_lncRNA), calculate the count and the position for the label
count_data <- all_filtered_data %>%
  dplyr::group_by(Condition, Overlaps_HERV) %>%
  dplyr::summarise(n = n(),
            .groups = "drop") %>%
  mutate(label_y = case_when(
    Overlaps_HERV == "TRUE" ~ 1.05,
    Overlaps_HERV == "FALSE" ~ 1.0
  ))

# Convert Condition into a factor within each Overlaps_HERV
all_filtered_data$Condition <- factor(all_filtered_data$Condition)

# Drop unused levels in Condition
all_filtered_data$Condition <- droplevels(all_filtered_data$Condition)

# Convert Condition into a factor within each Overlaps_HERV
all_filtered_data$Condition <- factor(all_filtered_data$Condition)

# Define the order of the levels
ordered_conditions <- c(
  "Overlap with 3' terminal exon vs no HERV overlap",
#  "Overlap with 5' terminal exon vs no HERV overlap",
#  "Overlap with single exon vs no HERV overlap",
  "Overlap with other exon vs no HERV overlap",
  "Overlap with intron vs no HERV overlap"
)

#detach(package:plyr)
library(dplyr)

# Set the order of the levels
all_filtered_data$Condition <- factor(all_filtered_data$Condition, levels = ordered_conditions)

# Drop unused factor levels
all_filtered_data$Condition <- factor(all_filtered_data$Condition)

#Drop empty tibbles from the list
filtered_data_list <- Filter(nrow, filtered_data_list)

p_values_GBM <- lapply(names(filtered_data_list), function(name) {
  df <- filtered_data_list[[name]]
  HERV_OL_df <- df %>% filter(Overlaps_HERV == "TRUE")
  NO_HERV_OL_df <- df %>% filter(Overlaps_HERV == "FALSE")
  
  HERV_OL_p_value <- calculate_p_value(HERV_OL_df$Disease_median_PDUI, NO_HERV_OL_df$Disease_median_PDUI)
  
  data.frame(
    Condition = name,
    Overlaps_HERV = c("TRUE", "FALSE"),
    p_value = HERV_OL_p_value
  )
})

p_values_Control <- lapply(names(filtered_data_list), function(name) {
  df <- filtered_data_list[[name]]
  HERV_OL_df <- df %>% filter(Overlaps_HERV == "TRUE")
  NO_HERV_OL_df <- df %>% filter(Overlaps_HERV == "FALSE")
  
  HERV_OL_p_value <- calculate_p_value(HERV_OL_df$Control_median_PDUI, NO_HERV_OL_df$Control_median_PDUI)
  
  data.frame(
    Condition = name,
    Overlaps_HERV = c("TRUE", "FALSE"),
    p_value = HERV_OL_p_value
  )
})

p_values_Diff <- lapply(names(filtered_data_list), function(name) {
  df <- filtered_data_list[[name]]
  HERV_OL_df <- df %>% filter(Overlaps_HERV == "TRUE")
  NO_HERV_OL_df <- df %>% filter(Overlaps_HERV == "FALSE")
  
  HERV_OL_p_value <- calculate_p_value(HERV_OL_df$PDUI_difference, NO_HERV_OL_df$PDUI_difference)
  
  data.frame(
    Condition = name,
    Overlaps_HERV = c("TRUE", "FALSE"),
    p_value = HERV_OL_p_value
  )
})

# Combine the list of dataframes into one dataframe
p_values_GBM <- do.call(rbind, p_values_GBM)
p_values_Control <- do.call(rbind, p_values_Control)
p_values_Diff <- do.call(rbind, p_values_Diff)

# Order the Overlaps_HERV in p_values
p_values_GBM <- p_values_GBM %>%
  mutate(Overlaps_HERV = factor(Overlaps_HERV, levels = c("TRUE", "FALSE")))

p_values_Control <- p_values_Control %>%
  mutate(Overlaps_HERV = factor(Overlaps_HERV, levels = c("TRUE", "FALSE")))

p_values_Diff <- p_values_Diff %>%
  mutate(Overlaps_HERV = factor(Overlaps_HERV, levels = c("TRUE", "FALSE")))


# Filter the data to remove rows with non-finite or missing values
#all_filtered_data <- all_filtered_data %>%
#  filter(!is.na(PDUI_difference) & PDUI_difference != Inf & PDUI_difference != -Inf)

combined_data <- all_filtered_data %>% 
  #semi_join(count_data, by = c("Condition", "PC_or_lncRNA", "Overlaps_HERV"))
  semi_join(count_data, by = c("Condition"))

# Get the maximum PDUI values for each Condition
max_PDUI <- combined_data %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarize(max_PDUI = max(Disease_median_PDUI, na.rm = TRUE))

# Assign label_y positions for each Condition  combination
p_value_labels_GBM <- p_values_GBM %>%
  left_join(max_PDUI, by = "Condition") %>%
  dplyr::arrange(Condition, Overlaps_HERV) %>%
  dplyr::group_by(Condition) %>%
  dplyr::mutate(label_y = 1 + 0.05 * (n() - row_number() + 1),  # Add 0.05 for each subsequent label within the same condition
         p_value_label = case_when(
           p_value < 0.001 ~ "***",
           p_value < 0.01 ~ "**",
           p_value < 0.05 ~ "*",
           TRUE ~ "ns"
         ))

p_value_labels_Control <- p_values_Control %>%
  left_join(max_PDUI, by = "Condition") %>%
  arrange(Condition, Overlaps_HERV) %>%
  group_by(Condition) %>%
  mutate(label_y = 1 + 0.05 * (n() - row_number() + 1),  # Add 0.05 for each subsequent label within the same condition
         p_value_label = case_when(
           p_value < 0.001 ~ "***",
           p_value < 0.01 ~ "**",
           p_value < 0.05 ~ "*",
           TRUE ~ "ns"
         ))

p_value_labels_Diff <- p_values_Diff %>%
  left_join(max_PDUI, by = "Condition") %>%
  arrange(Condition, Overlaps_HERV) %>%
  group_by(Condition) %>%
  mutate(label_y = 1 + 0.05 * (n() - row_number() + 1),  # Add 0.05 for each subsequent label within the same condition
         p_value_label = case_when(
           p_value < 0.001 ~ "***",
           p_value < 0.01 ~ "**",
           p_value < 0.05 ~ "*",
           TRUE ~ "ns"
         ))


### GBM PDUI ###
# Filter the data
filtered_data <- combined_data %>%
  group_by(Condition) %>%
  filter(n_distinct(Overlaps_HERV) > 1) %>%
  ungroup()

# Create a summary data frame
summary_data <- filtered_data %>%
  group_by(Condition, Overlaps_HERV) %>%
  dplyr::summarise(median_PDUI = median(Disease_median_PDUI, na.rm=TRUE), 
            max_PDUI = max(Disease_median_PDUI, na.rm=TRUE), 
            n = n())

# Create a new variable for the x-coordinate
summary_data$x_position <- ifelse(summary_data$Overlaps_HERV == "TRUE", as.numeric(summary_data$Condition) + 0.2, as.numeric(summary_data$Condition) - 0.2)

library(ggplot2)
library(grid)

# Create the plot
combined_plot <- ggplot(filtered_data, aes(x = Condition, y = Disease_median_PDUI, fill = Overlaps_HERV)) +
  geom_violin() +
  stat_summary(fun = median, geom = "crossbar", width = 0.2, aes(color = Overlaps_HERV), show.legend = FALSE) +
  stat_compare_means(aes(group = Overlaps_HERV), label="p.signif", label.x.npc = "left", method = "wilcox.test", size = 12, vjust = -0.5) +
  labs(
  title = "GBM Percent Distal Utilization",  # Add the title here
  x = NULL,
  y = NULL
  ) +
  theme_minimal() +
  scale_x_discrete(drop = TRUE, labels = c("Overlap with 3' terminal exon vs no HERV overlap" = "3’ Exon",
                                           "Overlap with other exon vs no HERV overlap" = "Middle Exon",
                                           "Overlap with single exon vs no HERV overlap" = "Single Exon",
                                           "Overlap with intron vs no HERV overlap" = "Intron")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)), breaks = seq(0, 1, by = 0.25), labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), guide = FALSE) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "Overlaps HERV") +
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5, size=28, face = "bold"),
        axis.text.x = element_text(size=26, hjust = 0.5, face = "bold", color = "black"),
        axis.text.y = element_text(size=26, face = "bold", color = "black"),  # Increase y-axis text size
        axis.title.y = element_blank(),  # Remove y-axis title
        legend.title = element_text(size=26, face = "bold"),
        legend.text = element_text(size=26, face = "bold"),
        plot.margin = margin(t = 100, r = 10, b = 100, l = 150))  # Increase left margin to create space for labels

# Print the plot to a grid object
grid.newpage()
plot_grob <- ggplotGrob(combined_plot)
grid.draw(plot_grob)

# Get the width of the plot panel
panel_index <- which(plot_grob$layout$name == "panel")
panel <- plot_grob$grobs[[panel_index]]

panel_width <- plot_grob$widths[plot_grob$layout[panel_index, "l"]]

# Calculate the center position of the plot panel
center_position <- unit(0.5, "npc") - panel_width / 2

# Calculate the width of the "HERV-Overlapping Feature" text
text_label <- "HERV-Overlapping Feature"
text_grob <- textGrob(text_label, gp = gpar(fontsize = 26, fontface = "bold"))
text_width <- convertWidth(grobWidth(text_grob), "npc", valueOnly = TRUE)

# Adjust the center position to account for the text width
adjusted_center_position <- center_position + (panel_width - unit(text_width, "npc")) / 2

# Add the "HERV-Overlapping Feature" label within the margin below the plot
grid.text(text_label, x = unit(0.5, "npc"), y = unit(0.08, "npc"), just = "center", gp = gpar(fontsize = 24, fontface="bold"))

# Draw the solid black line above the "HERV-Overlapping Feature" label
grid.lines(x = unit.c(adjusted_center_position, adjusted_center_position + unit(text_width, "npc")), y = unit(0.12, "npc"), gp = gpar(col = "black", lwd = 4))

## CTRL PDUI ###
# Create a summary data frame
summary_data <- filtered_data %>%
  group_by(Condition, Overlaps_HERV) %>%
  dplyr::summarise(median_PDUI = median(Control_median_PDUI, na.rm=TRUE), 
            max_PDUI = max(Control_median_PDUI, na.rm=TRUE), 
            n = n())

# Create a new variable for the x-coordinate
summary_data$x_position <- ifelse(summary_data$Overlaps_HERV == "TRUE", as.numeric(summary_data$Condition) + 0.2, as.numeric(summary_data$Condition) - 0.2)

combined_plot <- ggplot(filtered_data, aes(x = Condition, y = Control_median_PDUI, fill = Overlaps_HERV)) +
  geom_violin() +
  stat_summary(fun = median, geom = "crossbar", width = 0.2, aes(color = Overlaps_HERV), show.legend = FALSE) +
  stat_compare_means(aes(group = Overlaps_HERV), label="p.signif", label.x.npc = "left", method = "wilcox.test", size = 12, vjust = -0.5) +
  labs(
    title = "Control Percent Distal Utilization",  # Add the title here
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  scale_x_discrete(drop = TRUE, labels = c("Overlap with 3' terminal exon vs no HERV overlap" = "3’ Exon",
                                           "Overlap with other exon vs no HERV overlap" = "Middle Exon",
                                           "Overlap with single exon vs no HERV overlap" = "Single Exon",
                                           "Overlap with intron vs no HERV overlap" = "Intron")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)), breaks = seq(0, 1, by = 0.25), labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), guide = FALSE) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "Overlaps HERV") +
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5, size=28, face = "bold"),
        axis.text.x = element_text(size=26, hjust = 0.5, face = "bold", color = "black"),
        axis.text.y = element_text(size=26, face = "bold", color = "black"),  # Increase y-axis text size
        axis.title.y = element_blank(),  # Remove y-axis title
        legend.title = element_text(size=26, face = "bold"),
        legend.text = element_text(size=26, face = "bold"),
        plot.margin = margin(t = 100, r = 10, b = 100, l = 150))  # Increase left margin to create space for labels

# Print the plot to a grid object
grid.newpage()
plot_grob <- ggplotGrob(combined_plot)
grid.draw(plot_grob)

# Get the width of the plot panel
panel_index <- which(plot_grob$layout$name == "panel")
panel <- plot_grob$grobs[[panel_index]]

panel_width <- plot_grob$widths[plot_grob$layout[panel_index, "l"]]

# Calculate the center position of the plot panel
center_position <- unit(0.5, "npc") - panel_width / 2

# Calculate the width of the "HERV-Overlapping Feature" text
text_label <- "HERV-Overlapping Feature"
text_grob <- textGrob(text_label, gp = gpar(fontsize = 26, fontface = "bold"))
text_width <- convertWidth(grobWidth(text_grob), "npc", valueOnly = TRUE)

# Adjust the center position to account for the text width
adjusted_center_position <- center_position + (panel_width - unit(text_width, "npc")) / 2

# Add the "HERV-Overlapping Feature" label within the margin below the plot
grid.text(text_label, x = unit(0.5, "npc"), y = unit(0.08, "npc"), just = "center", gp = gpar(fontsize = 24, fontface = "bold"))

# Draw the solid black line above the "HERV-Overlapping Feature" label
grid.lines(x = unit.c(adjusted_center_position, adjusted_center_position + unit(text_width, "npc")), y = unit(0.12, "npc"), gp = gpar(col = "black", lwd = 4))

###

### PDUI DIFF ###
# Create a summary data frame
summary_data <- filtered_data %>%
  group_by(Condition, Overlaps_HERV) %>%
  dplyr::summarise(median_PDUI = median(PDUI_difference, na.rm=TRUE), 
            max_PDUI = max(PDUI_difference, na.rm=TRUE), 
            n = n())

# Create a new variable for the x-coordinate
summary_data$x_position <- ifelse(summary_data$Overlaps_HERV == "TRUE", as.numeric(summary_data$Condition) + 0.2, as.numeric(summary_data$Condition) - 0.2)

combined_plot <- ggplot(filtered_data, aes(x = Condition, y = PDUI_difference, fill = Overlaps_HERV)) +
  geom_violin() +
  stat_summary(fun = median, geom = "crossbar", width = 0.2, aes(color = Overlaps_HERV), show.legend = FALSE) +
  stat_compare_means(aes(group = Overlaps_HERV), label="p.signif", label.x.npc = "left", method = "wilcox.test", size = 12, vjust = -0.5) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  scale_x_discrete(drop = TRUE, labels = c("Overlap with 3' terminal exon vs no HERV overlap" = "3’ Exon",
                                           "Overlap with other exon vs no HERV overlap" = "Middle Exon",
                                           "Overlap with single exon vs no HERV overlap" = "Single Exon",
                                           "Overlap with intron vs no HERV overlap" = "Intron")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)), breaks = seq(-1.00, 1, by = 0.25), labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), guide = FALSE) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "Overlaps HERV") +
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5, size=28, face = "bold"),
        axis.text.x = element_text(size=26, hjust = 0.5, face = "bold", color = "black"),
        axis.text.y = element_text(size=26, face = "bold", color = "black"),  # Increase y-axis text size
        axis.title.y = element_blank(),  # Remove y-axis title
        legend.title = element_text(size=26, face = "bold"),
        legend.text = element_text(size=26, face = "bold"),
        plot.margin = margin(t = 100, r = 10, b = 100, l = 150))  # Increase left margin to create space for labels

# Print the plot to a grid object
grid.newpage()
plot_grob <- ggplotGrob(combined_plot)
grid.draw(plot_grob)

# Get the width of the plot panel
panel_index <- which(plot_grob$layout$name == "panel")
panel <- plot_grob$grobs[[panel_index]]

panel_width <- plot_grob$widths[plot_grob$layout[panel_index, "l"]]

# Calculate the center position of the plot panel
center_position <- unit(0.5, "npc") - panel_width / 2

# Calculate the width of the "HERV-Overlapping Feature" text
text_label <- "HERV-Overlapping Feature"
text_grob <- textGrob(text_label, gp = gpar(fontsize = 26, fontface = "bold"))
text_width <- convertWidth(grobWidth(text_grob), "npc", valueOnly = TRUE)

# Adjust the center position to account for the text width
adjusted_center_position <- center_position + (panel_width - unit(text_width, "npc")) / 2

# Add the "HERV-Overlapping Feature" label within the margin below the plot
grid.text(text_label, x = unit(0.5, "npc"), y = unit(0.08, "npc"), just = "center", gp = gpar(fontsize = 24, fontface = "bold"))

# Draw the solid black line above the "HERV-Overlapping Feature" label
grid.lines(x = unit.c(adjusted_center_position, adjusted_center_position + unit(text_width, "npc")), y = unit(0.12, "npc"), gp = gpar(col = "black", lwd = 4))
