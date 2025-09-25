#====Methylation Discordant Site Analysis================================
#This analysis looks at the number of concordant versus discordant sites
#in all 66 pairwise comparisons of UT9728 ORA assemblies.
#It is looking at the difference in fraction of reads methylated at
#every matching site in the pairwise assemblies.
#========================================================================


#load libraries
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(tidyr)


#all pairwise comparisons
# -------- 1. Input files --------
files <- c(
  "path/to/ORA/assembly/UT9728_HWa_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_HWb_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_IPa_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_IPb_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_IPc_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_MLa_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_MLb_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_MLc_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_OHa_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_OHb_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_TSa_meth_sites_Ont.csv",
  "path/to/ORA/assembly/UT9728_TSb_meth_sites_Ont.csv")

# -------- 2. Labels --------
labels <- c("UT9728_HWa", "UT9728_HWb", "UT9728_IPa", "UT9728_IPb", "UT9728_IPc","UT9728_MLa", 
            "UT9728_MLb", "UT9728_MLc", "UT9728_OHa", "UT9728_OHb", "UT9728_TSa", "UT9728_TSb")

# -------- 3. Load all assemblies into a named list --------
assemblies <- map2(files, labels, ~ {
  read.csv(.x) %>%
    filter(motif == "GATC") %>%
    mutate(Assembly = .y)
}) %>% set_names(labels)

# -------- 4. Comparison function --------
#create function to join tables based on defined site key and perform discordant site analysis
compare_assemblies <- function(df1, df2, label1, label2) {
  joined <- inner_join(
    df1, df2,
    by = "Sequence",
    suffix = c(paste0("_", label1), paste0("_", label2)),
    relationship = "many-to-many"
  ) %>%
    filter(
      .data[[paste0("Strand_", label1)]] == .data[[paste0("Strand_", label2)]] |
        abs(.data[[paste0("Position_", label1)]] - .data[[paste0("Position_", label2)]]) <= 5
    ) %>%
    mutate(
      Diff_methylation = abs(.data[[paste0("Percent_modified_", label1)]] -
                               .data[[paste0("Percent_modified_", label2)]]),
      Discordant = Diff_methylation > 0.15,
      Comparison = paste(label1, "vs", label2, sep = "_")
    )
  
  # keep only relevant columns
  result <- joined %>%
    select(
      !!paste0("Contig_", label1) := paste0("Contig_", label1),
      !!paste0("Contig_", label2) := paste0("Contig_", label2),
      Strand = paste0("Strand_", label1),
      !!paste0("Position_", label1) := paste0("Position_", label1),
      !!paste0("Position_", label2) := paste0("Position_", label2),
      !!paste0("Coverage_", label1) := paste0("Total_coverage_", label1),
      !!paste0("Coverage_", label2) := paste0("Total_coverage_", label2),
      !!paste0("Percent_modified_", label1) := paste0("Percent_modified_", label1),
      !!paste0("Percent_modified_", label2) := paste0("Percent_modified_", label2),
      Sequence,
      !!paste0("motif_", label1) := paste0("motif_", label1),
      Diff_methylation,
      Discordant,
      Comparison)
  
  return(result)}

# -------- 5. Run and save all 66 pairwise comparisons --------
outdir <- "path/to/output/directory"
dir.create(outdir, showWarnings = FALSE)

pairwise_results <- list()

combinations <- combn(labels, 2, simplify = FALSE)

for (pair in combinations) {
  label1 <- pair[1]
  label2 <- pair[2]
  
  res <- compare_assemblies(assemblies[[label1]], assemblies[[label2]], label1, label2)
  
  # Save individual CSV
  outfile <- file.path(outdir, paste0("pairwise_", label1, "_vs_", label2, ".csv"))
  write.csv(res, outfile, row.names = FALSE)
  
  # Store in list
  pairwise_results[[paste(label1, label2, sep = "_vs_")]] <- res}

# -------- 6. Summary table --------
#create a summary table containing comparison, total methylated sites,
#total discordant sites, and percentage of discordant sites
summary_table <- bind_rows(pairwise_results) %>%
  group_by(Comparison) %>%
  summarize(
    Total_sites = n(),
    Discordant_sites = sum(Discordant),
    Percent_discordant = mean(Discordant) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(Percent_discordant))

# Save summary table too
write.csv(summary_table, file.path(outdir, "summary_table.csv"), row.names = FALSE)


# Output directory for plots
plot_dir <- "path/to/plot/directory"
dir.create(plot_dir, showWarnings = FALSE)

# Loop over all 66 comparisons and plot discordant sites heatmap
for (comp_name in names(pairwise_results)) {
  df <- pairwise_results[[comp_name]]
  
  # extract the two assembly labels from the Comparison column
  labels <- strsplit(comp_name, "_vs_")[[1]]
  label1 <- labels[1]
  label2 <- labels[2]
  
  # dynamically get the Percent_modified columns
  x_col <- paste0("Percent_modified_", label1)
  y_col <- paste0("Percent_modified_", label2)
  
  # dynamically build axis labels
  x_label <- paste0(label1, ": Fraction methylated")
  y_label <- paste0(label2, ": Fraction methylated")
  
  # make plot
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(color = Discordant), alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("FALSE" = "skyblue", "TRUE" = "tomato"),
                       labels = c("FALSE" = "Concordant", "TRUE" = "Discordant")) +
    labs(x = x_label,y = y_label,color = "Discordant (>0.15)") +
    theme_classic(base_size = 14)
  
  # save plot
  ggsave(filename = file.path(plot_dir, paste0("plot_", comp_name, ".jpeg")),
         plot = p, width = 6, height = 6)}


#plot heatmap of discordant sites percentage
summary_table <- read.csv("path/to/summary_table.csv")

# 1. Split into Assembly1 / Assembly2
heatmap_data <- summary_table %>%
  separate(Comparison, into = c("Assembly1", "Assembly2"), sep = "_vs_")

# 2. Duplicate rows with swapped assemblies so heatmap is not one sided
heatmap_data_full <- bind_rows( heatmap_data,heatmap_data %>% 
                                  rename(Assembly1 = Assembly2, Assembly2 = Assembly1))

# 3. Add self-comparisons (0% discordant)
assemblies <- unique(c(heatmap_data$Assembly1, heatmap_data$Assembly2))
self_comparisons <- expand.grid(Assembly1 = assemblies, Assembly2 = assemblies) %>%
  filter(Assembly1 == Assembly2) %>%
  mutate(Percent_discordant = 0)

heatmap_data_full <- bind_rows(heatmap_data_full, self_comparisons)

# 4. Plot full heatmap
p <- ggplot(heatmap_data_full, aes(x = Assembly1, y = Assembly2, fill = Percent_discordant)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "viridis", name = "% Discordant Sites") +
  labs(x = "Assembly 1", y = "Assembly 2") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 0)) +
  ggtitle("Heatmap of Discordant Site Percentage Per Pairwise Comparison")

print(p)

# 5. Save
ggsave("path/to/output/directory",
       plot = p, width = 8, height = 6, dpi = 600)

#pairwise discordance vs coverage
# -------- Output directory for plots --------
plot_dir <- file.path(outdir, "coverage_scatterplots")
dir.create(plot_dir, showWarnings = FALSE)

# -------- Loop over all pairwise results --------
#plot all sites in a scatterplot with discordant sites colored red and assembly coverages on axes
for (comp_name in names(pairwise_results)) {
  
  res <- pairwise_results[[comp_name]]
  
  # Extract assembly names from the column headers
  cov_cols <- grep("^Coverage_", colnames(res), value = TRUE)
  asm1 <- sub("Coverage_", "", cov_cols[1])
  asm2 <- sub("Coverage_", "", cov_cols[2])
  
  # Make scatter plot with discordant plotted last
  p <- ggplot(res, aes_string(x = cov_cols[1], y = cov_cols[2])) +
    geom_point(data = subset(res, Discordant == FALSE), aes(color = "FALSE"), alpha = 0.6) +
    geom_point(data = subset(res, Discordant == TRUE), aes(color = "TRUE"), alpha = 0.8) +
    scale_color_manual(values = c("FALSE" = "skyblue", "TRUE" = "tomato"),
                       labels = c("FALSE" = "Concordant", "TRUE" = "Discordant")) +
    labs(title = paste0("Coverage correlation: ", asm1, " vs ", asm2),
         x = paste0(asm1, " Coverage"), y = paste0(asm2, " Coverage"), color = "Discordant") +
    theme_bw(base_size = 14)
  
  # Save plot
  ggsave(filename = file.path(plot_dir, paste0("scatter_", asm1, "_vs_", asm2, ".jpeg")),
         plot = p, width = 6, height = 5, dpi = 600)}
