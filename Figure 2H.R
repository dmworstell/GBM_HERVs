suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

#base_dir <- "/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/results/bootstrapping_results"
base_dir = 'C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/results/bootstrapping_results'

lncRNA_DE_prvs <- read.csv(file.path(base_dir, "lncRNA_proviruses_1000_bootstraps_7146_MT_80percentExon_DEHervs.csv"), row.names = 1, check.names=FALSE)
lncRNA_All_Proviruses <- read.csv(file.path(base_dir, "lncRNA_proviruses_1000_bootstraps_7146_MT_80percentExon_AllHervs.csv"), row.names = 1, check.names=FALSE)
mRNA_All_prvs <- read.csv(file.path(base_dir, "mRNA_proviruses_1000_bootstraps_7146_MT_80percentExon_AllHervs.csv"), row.names = 1, check.names=FALSE)
mRNA_DE_Proviruses <- read.csv(file.path(base_dir, "mRNA_proviruses_1000_bootstraps_7146_MT_80percentExon_DEHervs.csv"), row.names = 1, check.names=FALSE)

Raw_Bootstrap_Data <- bind_rows(
  lncRNA_DE_prvs %>% mutate(`DE or All HERVs` = "DE HERVs", `Transcript Type` = "lncRNA"),
  lncRNA_All_Proviruses %>% mutate(`DE or All HERVs` = "All HERVs", `Transcript Type` = "lncRNA"),
  mRNA_All_prvs %>% mutate(`DE or All HERVs` = "All HERVs", `Transcript Type` = "mRNA"),
  mRNA_DE_Proviruses %>% mutate(`DE or All HERVs` = "DE HERVs", `Transcript Type` = "mRNA")
)

# Function to process the data
process_data <- function(data) {
  summary_data <- data.frame(
    Category = colnames(data),
    Median = apply(data, 2, median),
    Q1 = apply(data, 2, quantile, probs = 0.25),
    Q3 = apply(data, 2, quantile, probs = 0.75),
    Lower_Whisker = apply(data, 2, function(x) max(min(x), quantile(x, 0.25) - 1.5 * IQR(x))),
    Upper_Whisker = apply(data, 2, function(x) min(max(x), quantile(x, 0.75) + 1.5 * IQR(x)))
  )
  
  return(summary_data)
}

# Process lncRNA DE Proviruses data
lncRNA_DE_summary <- process_data(lncRNA_DE_prvs) %>%
  mutate(Type = "DE Proviruses", RNA = "lncRNA")

# Process lncRNA All Proviruses data
lncRNA_All_summary <- process_data(lncRNA_All_Proviruses) %>%
  mutate(Type = "All Proviruses", RNA = "lncRNA")

# Process mRNA DE Proviruses data
mRNA_DE_summary <- process_data(mRNA_DE_Proviruses) %>%
  mutate(Type = "DE Proviruses", RNA = "mRNA")

# Process mRNA All Proviruses data
mRNA_All_summary <- process_data(mRNA_All_prvs) %>%
  mutate(Type = "All Proviruses", RNA = "mRNA")

# Combine all summaries into a single data frame
combined_data <- bind_rows(lncRNA_DE_summary, lncRNA_All_summary, mRNA_DE_summary, mRNA_All_summary)

# Define the new order for the categories
new_order <- c("+ | First", "- | First", "+ | Middle", "- | Middle", "+ | Last", "- | Last", "+ | Intron", "- | Intron")

# Filter out the "Single" categories
combined_data <- combined_data %>%
  filter(!Category %in% c("+ | Single", "- | Single"))

# Update the Category factor levels
combined_data$Category <- factor(combined_data$Category, levels = new_order)

# Create the interaction term and set its levels to the desired order
combined_data$Interaction <- interaction(combined_data$RNA, combined_data$Type)
combined_data$Interaction <- factor(combined_data$Interaction, levels = c("lncRNA.DE Proviruses", "lncRNA.All Proviruses", "mRNA.DE Proviruses", "mRNA.All Proviruses"))

# Plot the combined data
ggplot(combined_data, aes(x = Category, y = Median, fill = Interaction)) +
  geom_boxplot(aes(
    lower = Q1, upper = Q3, middle = Median, 
    ymin = Lower_Whisker, ymax = Upper_Whisker
  ), stat = "identity", position = position_dodge2(preserve = "single",width=1),width=0.7) +
  #scale_x_discrete(expand = c(0, 0)) + 
  scale_y_log10() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 32),
        axis.title.y = element_text(size = 26, color = "black"),
        axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24, color = "black"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  labs(title = "HERV / RNA Overlaps", y = "Median Counts (1000 bootstraps)", x = NULL, fill = "RNA Type") +
  scale_fill_manual(values = c("lncRNA.DE Proviruses" = "blue", "lncRNA.All Proviruses" = "lightblue",
                               "mRNA.DE Proviruses" = "red", "mRNA.All Proviruses" = "pink"),
                    labels = c("lncRNA DE Proviruses", "lncRNA All Proviruses", "mRNA DE Proviruses", "mRNA All Proviruses"))

library(openxlsx)

# Define the file path
#file_path <- '/Users/Daniel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Manuscripts/2024/Supplemental/Supporting_Data.xlsx'
file_path <- 'C:/Users/worst/iCloudDrive/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Manuscripts/2024/Supplemental/Supporting_Data.xlsx'

# Define the sheet name
sheet_name <- "Figure 2C"

# Check if the file exists
if (file.exists(file_path)) {
  # Load the existing workbook
  wb <- loadWorkbook(file_path)
} else {
  # Create a new workbook
  wb <- createWorkbook()
}

# Add a new sheet with the specified name
addWorksheet(wb, sheet_name)

# Write the data to the new sheet
writeData(wb, sheet_name, Raw_Bootstrap_Data)

# Save the workbook
saveWorkbook(wb, file_path, overwrite = TRUE)



## FOR INSET

# Filter for intron-overlapping categories only
introns_data <- combined_data %>%
  filter(Category %in% c("+ | Intron", "- | Intron"))

# Create a plot on a log scale, restricted to y values from 500 to 5000
ggplot(introns_data, aes(x = Category, y = Median, fill = Interaction)) +
  geom_boxplot(aes(
    lower = Q1, upper = Q3, middle = Median, 
    ymin = Lower_Whisker, ymax = Upper_Whisker
  ), stat = "identity", position = position_dodge2(preserve = "single",width=1),width=0.7) +
  #scale_x_discrete(expand = c(0, 0)) +
  # Use a log10 scale with limits = c(500, 5000)
  scale_y_log10(limits = c(700, 3000)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 32),
        axis.title.y = element_text(size = 26, color = "black"),
        axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  labs(title = NULL, 
       y = NULL, 
       x = NULL, 
       fill = "RNA Type") +
  scale_fill_manual(values = c("lncRNA.DE Proviruses" = "blue", 
                               "lncRNA.All Proviruses" = "lightblue",
                               "mRNA.DE Proviruses" = "red", 
                               "mRNA.All Proviruses" = "pink"),
                    labels = c("lncRNA DE Proviruses", 
                               "lncRNA All Proviruses", 
                               "mRNA DE Proviruses", 
                               "mRNA All Proviruses"))