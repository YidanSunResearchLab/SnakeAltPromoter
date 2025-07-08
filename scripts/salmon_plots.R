#!/usr/bin/env Rscript

# ---------------------
# Salmon plots script
# Usage:
# Rscript salmon_plots.R <output_dir> <promoter_rds> <activity> <gene_expr> <cell_lines> <condition>
# ---------------------


# --------------------
# Load Libraries
# --------------------
library(proActiv)
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
library(Rtsne)
library(future.apply)
library(ggpubr)

# --------------------
# Parse Inputs
# --------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript salmon_plots.R <output_dir> <promoter_rds> <activity> <gene_expr> <cell_lines> <condition>")
}

output_dir   <- args[1]
promoter_rds_path <- args[2]
activity_path     <- args[3]
gene_expr_path    <- args[4]
cell_lines   <- strsplit(args[5], " ")[[1]]
condition_all <- strsplit(args[6], ",")[[1]]

# --------------------
# Validate Inputs

stopifnot(file.exists(promoter_rds_path),
          file.exists(activity_path),
          file.exists(gene_expr_path))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --------------------
# Load Data and filter outliers
# --------------------
activity      <- readRDS(activity_path)
samples_all <- colnames(activity)
condition <- condition_all
activity <- readRDS(activity_path)
gene_expr <- readRDS(gene_expr_path)

print(cell_lines)
print(condition)
print(head(activity))
print(head(gene_expr))
stopifnot(length(condition) == ncol(activity))

# If gene_expr uses geneClean it aligns with promoterId
if (!identical(rownames(gene_expr), rownames(activity))) {
  if (all(rownames(activity) %in% rownames(gene_expr))) {
    gene_expr <- gene_expr[rownames(activity), , drop = FALSE]
  } else {
    stop("gene_expr and activity have mismatched promoterId sets")
  }
}

if (!identical(rownames(gene_expr), rownames(activity))) {
  if (nrow(gene_expr) == nrow(activity) && all(rownames(activity) %in% rownames(gene_expr))) {
    gene_expr <- gene_expr[rownames(activity), , drop = FALSE] 
  } else {
    stop("gene_expr and activity have different row sets; please supply gene_expr with promoterId rows or equal length to activity")
  }
}

# --------------------
# Load promoter annotations
# --------------------
promoter_anno <- readRDS(promoter_rds_path)
#Extract the core mapping between promoterId and geneId
#Sort promoters by position: on + strand ascending start while on - strand descending start
rowData <- promoterIdMapping(promoter_anno)
# Add geneClean for aggregation
rowData$geneClean <- sub("\\..*$", "", rowData$geneId)

coords  <- as.data.frame(promoter_anno@promoterCoordinates)[ ,
            c("promoterId","seqnames","start","strand")] %>%
  left_join(rowData, by = "promoterId") %>%
  group_by(geneClean) %>%
  arrange(strand, if_else(strand == "-", desc(start), start)) %>%
  mutate(promoterPosition = row_number()) %>%
  ungroup() %>%
  select(promoterId, promoterPosition)

# Ensure consistent datatype before joining
rowData$promoterId <- as.character(rowData$promoterId)   # int → char
coords$promoterId  <- as.character(coords$promoterId)

# Merge promoter position back to rowData table
# Warning and deduplicate if multiple annotations found for the same promoter
dup_promoters <- coords$promoterId[duplicated(coords$promoterId)]
if (length(dup_promoters) > 0) {
  coords <- coords[!duplicated(coords$promoterId), ]  
  warning("Deduplicated. Multiple entries found for some promoterIds in coords:\n",
          paste(head(dup_promoters), collapse = ", "))
}  
rowData <- left_join(rowData, coords, by = "promoterId")
# Keep same order as activity
rowData <- rowData[match(rownames(activity), rowData$promoterId), ]
# Ensure promoter IDs match between activity and rowData
# Make sure every promoter has assigned class
stopifnot(identical(rownames(activity), rowData$promoterId),
          !any(is.na(rowData$class)))

# --------------------
# Final checks
# --------------------

cat("Final dimensions:\n")
cat("rowData: ", nrow(rowData), "\n")
cat("activity: ", nrow(activity), "\n")
cat("gene_expr: ", nrow(gene_expr), "\n")

stopifnot(nrow(rowData) == nrow(activity))
stopifnot(nrow(rowData) == nrow(gene_expr))
stopifnot(!any(is.na(rowData$promoterId)))

# ----- Theme and Color for all plots -------------------------------------------
plot_theme <- theme_classic() + theme(
  plot.title = element_text(hjust = .5, family = "Helvetica"),
  axis.text  = element_text(colour = "black", family = "Helvetica"),
  axis.title = element_text(colour = "black", family = "Helvetica"))
cols_type <- c("Major" = "#ff1e56", "Minor" = "#ffac41", "Inactive" = "#323232")

make_plots <- function(cl) {
  message("Processing cell line: ", cl)

  cl_dir <- file.path(output_dir, cl)
  dir.create(cl_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Generate column names for the cell line
  col_class <- paste0(cl, ".class")
  col_gene  <- paste0(cl, ".gene.mean")
  col_mean  <- paste0(cl, ".mean")
  # Find the column corresponding to the cell line in activity matrix
  #idx <- which(startsWith(colnames(activity), cl))
  # Use grep because sample name does not start with cell line
  idx <- grepl(paste0("_", cl, "_"), colnames(activity))
  message("[", cl, "] length(idx) = ", length(idx))
  message("First few matching columns: ", paste(head(colnames(activity)[idx]), collapse = ", "))

  # Calculate the average promoter activity for the cell line
  avg_cl <- rowMeans(activity[ , idx, drop = FALSE], na.rm = TRUE)
  # Ensure IDs are are contained in both activity and rowData
  rowData$promoterId <- as.character(rowData$promoterId)
  rownames(activity) <- as.character(rownames(activity))
  if (!all(as.character(rowData$promoterId) %in% names(avg_cl))) {
  stop("Some promoter IDs in rowData are missing from activity matrix.")
}
# Store promoter-level average activity and gene level average activity in rowData
# Uses the order of rowData$promoterId to extract names from avg_cl, thus able to match even if rowData and avg_cl are in different orders
# Returns na if Id not found in avg_cl
# Safety check to ensure all Ids found and are matched in the correct order
  rowData[[col_mean]] <- avg_cl[as.character(rowData$promoterId)]
  if (any(is.na(avg_cl[as.character(rowData$promoterId)]))) {
  stop("Some promoter IDs from rowData not found in activity rownames.")
}
if (!identical(as.character(rowData$promoterId), names(avg_cl))) {
  stop("Mismatch in promoterId order between rowData and activity")
}
  rowData[[col_gene]] <- rowMeans(gene_expr[ , idx, drop = FALSE], na.rm = TRUE)
  # Cell-line specific classification for an activity that is considered as active
  threshold_minor <- 0.25
  
  classify_one_gene <- function(indices) {
  # Retrieve the average promoter activity values for all promoters of this gene
  v <- rowData[[col_mean]][indices]
  # Initialize all promoters as "Inactive" by default
  labs <- rep("Inactive", length(v))
  # If no activity is above the threshold, maintain default as "Inactive" and exit function
  if (all(is.na(v)) || max(v, na.rm = TRUE) < threshold_minor) {
    return(labs)
  }
  # Only one active promoter per gene with the highest activity is classify as "Major"
  major_idx <- which.max(v)
  labs[major_idx] <- "Major"
  # The rest of the promoters with activity above the threshold are classified as "Minor"
  minor_idx <- which(v >= threshold_minor & seq_along(v) != major_idx)
  labs[minor_idx] <- "Minor"
  
  return(labs)
}
  # Split all promoter rows into groups by geneId. Each group contains indices of promoters for that gene.
  gene_groups <- split(seq_len(nrow(rowData)), rowData$geneId)
  # Apply classification to each group of promoters sharing the same gene identifier 
  # Classify all promoters from activity of its corresponding gene
  class_labels <- unlist(lapply(gene_groups, classify_one_gene))
  # Store the classification in a new column, maintaining the order of classification.
  rowData[[col_class]] <- factor(class_labels, levels = c("Major", "Minor", "Inactive"))
  # Extract a list of genes with at least one active promoter for downstream alternative promoter analysis.
  minor_genes <- rowData$geneId[rowData[[col_class]] == "Minor"]

# Plot 1: promoter_activity_category_percentage --------------------

# A bar plot showing the percentage of promoters in each category(level): Major, Minor, Inactive, for a specific cell line.
# y axis is the number of promoters
# Bar is divided to segments of level (Major, Minor, Inactive), each with promoter counts labelled.

  plot_cat_percentage <- function() {
    df_cat <- rowData %>%
      # Filter out NA values to avoid counting them.
      filter(!is.na(.data[[col_class]])) %>%
      # Select column specific to one cell line, count the number of promoters in each category, and set fixed order for factor levels.
      # Plot will be stacked bar with Major, Minor, Inactive from top to bottom.  
      count(Type = factor(.data[[col_class]], levels = c("Major","Minor","Inactive")))
    # Create a single stacked bar plot with the promoter counts per category
    ggplot(df_cat, aes(x = 1, y = n, fill = Type)) +
      # Stacked bar with black border
      geom_col(width = .6, colour = "black") +
      # Add count labels on each stack segment; text adjusted to be inside each segment
      geom_text(aes(label = comma(n)), position = position_stack(.8),
                colour = "white", size = 1.8, fontface = "bold") +
      # Apply defined colors and theme
      scale_fill_manual(values = cols_type) +
      plot_theme +
      # Remove x-axis related because only one bar: one graph for each cell line
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), legend.position = "right") +
      labs(y = "Number of promoters", x = NULL) 
  }

    # Plot 2: promoter_activity_category_comparison------------------

    # A boxplot comparing promoter activity across different promoter categories (Major / Minor / Inactive) for the given cell line.
    # y axis is the average promoter activity; x axis is Major, Minor, Inactive.
    # 3 bars with different colors for each promoter category.
    plot_activity_comparison <- function(){ 
        # Filter so promoter activity >0
        df2 <- rowData[rowData$geneId %in% minor_genes, ]
        df2[[paste0(cl, "_mean")]] <- rowData[[paste0(cl, "_mean")]][match(df2$promoterId, rowData$promoterId)]
        df2[[col_class]] <- rowData[[col_class]][match(df2$promoterId, rowData$promoterId)]
        df2 <- df2[!is.na(df2[[col_class]]), ]
        # Plot filtered data with x as promoter category, y as average promoter activity, and fill accordingly by promoter category.

        # Plot object assigned to variable
        ggplot(df2, aes(x = .data[[col_class]], y = .data[[col_mean]], fill = .data[[col_class]])) +
            geom_boxplot(outlier.shape = NA, notch = TRUE, width = .5, colour = "black") +
            scale_fill_manual(values = cols_type, name = "Type") +
            stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Major",
                            na.rm = TRUE, label.y.npc = "top", tip.length = 0.01) +
            plot_theme +
            theme(
                legend.position = "right",
                axis.text.x = element_text(size = 6)
            ) +
            labs(y = "Average promoter activity", x = NULL)
    }

# Plot 3: Promoter position category-----------------------------

# A continuous bar plot showing the proportion of Major and Minor promoters at different promoter positions.
# y axis is the proportion of Major and Minor promoters at each position; x axis is the promoter position (1 to 5, with >=5 grouped).
# Number of promoters at each position is labelled on the top of the bar.
# The percentage of total promoter positions larger than 1 out of all promoter positions is labelled at the top of the bar.

  plot_position_category <- function() {
    # Cap promoter positions greater than 5 as 5. Further promoter counts are usually not large so could be grouped.
    rowData$promoterPosition <- ifelse(rowData$promoterPosition > 5, 5, rowData$promoterPosition)
    # Filter to include only active promoters (Major and Minor)
    filtered <- rowData[rowData[[col_class]] %in% c("Major", "Minor"), ]
      if (nrow(filtered) == 0) {
    message(sprintf("No Major or Minor promoters found for cell line %s. Skipping position plot.", cl))
    return(NULL)
  }
    # Count promoters per position and class
    # position_counts has three columns: Position (1 to 5), Type (Major or Minor), and Freq (number of promoters).
    position_counts <- as.data.frame(table(
      Position = filtered$promoterPosition,
      Type = filtered[[col_class]]))
    # Set factor levels for plotting order
    # Bar will be stacked with 1 to 5 from left to right; descending from major to minor.
    position_counts$Position <- factor(position_counts$Position, levels = 1:5)
    position_counts$Type <- factor(position_counts$Type, levels = c("Major", "Minor"))
    # Convert counts to wide format (row = Type, columns = Position) to calculate total counts per position and bar widths.
    # plot_data has two columns: Type (Major or Minor); Position (1 to 5) as column names; counts as values.
    # Fill missing positions with 0 counts to prevent NA.
    plot_data <- tidyr::pivot_wider(position_counts, names_from = Position, values_from = Freq, values_fill = 0)
    # may return a tibble with character columns, and ggplot2 might fail to match character values to factor levels in scale_fill_manual().
    # This can lead to incorrect coloring (e.g., all gray bars) and a legend labeled "NA".
    # Converting to data.frame ensures compatibility with factor-based fill mappings.
    plot_data <- as.data.frame(plot_data)
    #####rownames(plot_data) <- plot_data$Type
    # Remove the Type column, only keeping values.
    #####plot_data <- plot_data[, -1]
    # Add Type as a normal column for ggplot, because it needs an explicit column.
    #plot_data$Type <- rownames(plot_data)
    #plot_data$Type <- as.character(plot_data$Type) 
    # Ensure order of Type and Position for plotting
    plot_data$Type <- factor(plot_data$Type, levels = c("Major", "Minor"))
    # Reshape back to long format for ggplot.
    # Type as grouping variable, Position as variable, Count as value.
    plot_data_melt <- reshape2::melt(plot_data, id.vars = "Type", variable.name = "Position", value.name = "Count")
    plot_data_melt$Position <- factor(plot_data_melt$Position, levels = as.character(1:5))
    # Compute bar widths and center positions based on total counts.
    # This step is necessary to ensure bar labels match with bar positions, 
    # given the bar widths are proportional to the number of promoters, which differs for each position.
    plot_data_melt$Count <- as.numeric(plot_data_melt$Count)
    barwidth <- plot_data_melt %>% group_by(Position) %>% summarise(total = sum(Count)) %>%
    # Calculate the right side of bar and minus half to get the center position.
      mutate(width = log10(total / sum(total) * 100), center = cumsum(width) - width / 2)
    # Merge back to main data for plotting
    # plot_data_melt has position, type, and count; barwidth has total, width, center for each position
    # Every merge stack will have count, width, and center position.
    plot_data_final <- merge(plot_data_melt, barwidth, by = "Position") %>%
      # Sum at each position; stack from the descending sequence of major at the bottom.
      group_by(Position) %>% arrange(desc(Type)) %>%
      # prop is the proportion of count out of all, which corresponds to the height of the bar segment since max is 1.
      # ypos centerizes text label, which is the cumulative proportion minus half of the proportion of the segment.
      mutate(prop = Count / sum(Count), ypos = cumsum(prop) - prop / 2)
    # Set levels for the factor Type to ensure Major and minor are colored/
    plot_data_final$Type <- factor(plot_data_final$Type, levels = c("Major", "Minor"))
    # Compute percentage out of all promoters for promoter counts after position 1
    pct <- round(sum(barwidth$total[2:5]) / sum(barwidth$total) * 100)

    # Generate stacked bar plot
    ggplot(plot_data_final,
           aes(x = center, y = prop, fill = Type, width = width)) +
      geom_col(position = "stack", colour = "black") +
      # Add white count labels on top of each stacked bar
      # Adjusted y position to be slightly below the bar top
      geom_text(data = barwidth, inherit.aes = FALSE,
                aes(x = center, y = 0.97, label = comma(total)),
                colour = "white", size = 1.8, fontface = "bold", vjust = 0) +
      # Draw a horizontal segment over positions 2-5, adjust y position to be slightly above the bar top.
      annotate("segment", x = barwidth$center[2] - barwidth$width[2]/2,
               xend = barwidth$center[5] + barwidth$width[5]/2,
               y = 1.07, yend = 1.07, size = 0.3) +
      # Add percentage label above the segment, adjust y position to be slightly above the segment.
      annotate("text", x = mean(barwidth$center[2:5]), y = 1.11,
               label = paste0(pct, "%"), fontface = "bold") +
      # Format y axis with % labels, expand to let the graph not touch y axis.
      scale_y_continuous(breaks = seq(0, 1, 0.2),
                         labels = paste0(seq(0, 100, 20), "%"),
                         expand = c(0, 0.18)) +
      # Format x axis with custom widths and center labels, mult to let graph not touch x axis.
      scale_x_continuous(breaks = barwidth$center,
                         labels = c("1", "2", "3", "4", ">=5"),
                         expand = expansion(mult = c(0.05, 0.05))) +
      # Assign fill colors for promoter types
      scale_fill_manual(values = c("Major" = "#ff1e56", "Minor" = "#ffac41"), na.translate = FALSE) +
      # Set coordinate limits and axis clipping, so that the segment and percentage labels above the bars are not clipped.
      coord_cartesian(ylim = c(0, 1.12), clip = "off") +
      labs(x = expression("Promoter position (5'" %->% "3')"),
           y = "Proportion of promoter types") +
      plot_theme +
      theme(legend.title = element_blank(), plot.margin = margin(6, 6, 6, 6))
  }

  # Plot 4: promoter_activity_geneexpression_correlation---------------------------------------

  # A scatter plot showing the correlation between average gene expression and average promoter activity for Major promoters.
  # The closer from the slope = 1 line, the more contribution the major promoter has.

  plot_gene_corr <- function() {
    # Only keep Major promoters for correlation analysis.
    #### Filter 0 value
    majors <- rowData %>% filter(.data[[col_class]] == "Major") %>% filter(.data[[col_mean]] > 0, .data[[paste0(cl, ".gene.mean")]] > 0)
    # x axis is the average gene expression; y axis is the average promoter activity.
    ggplot(majors, aes(x = .data[[paste0(cl, ".gene.mean")]], y = .data[[col_mean]])) +
      # Points are colored black, with a coordinate limit of 0 to 15 for both axes.
      geom_point(size = .8, alpha = .8, colour = "black") +
      coord_cartesian(xlim = c(0, 15), ylim = c(0, 15)) +
      plot_theme +
      labs(x = "Average gene expression", y = "Average promoter activity")
  }

  # Plot 5: promoter_activity_category_percentage_genewisegenewise-----------------------------------------

  # A bar plot showing the number of genes with Major, Minor, and Inactive promoters.
  # y axis is the number of genes; bar stacks counts from each promoter category (Major, Minor, Inactive).
  # Differs from Plot 1 in that inactive counts in genes that include active promoters are not counted, so inactive counts are lower.

  plot_genewise_percentage <- function() {
    # Group all promoters by geneId. Each gene is now a list of promoter classes its promoters belong to.
    gene_split <- split(rowData[[col_class]], rowData$geneId)
    # Count the number of genes with at least one Major or Minor promoter.
    # Classify each gene based on the promoters it has
    classify_gene <- function(x) {
      if ("Major" %in% x && !"Minor" %in% x) {
        return("Major")
      } else if ("Minor" %in% x) {
        return("Minor")
      } else {
        return("Inactive")
      }
    }

    # Apply to each gene group
    gene_category <- sapply(gene_split, classify_gene)

    # Count each category
    counts <- table(factor(gene_category, levels = c("Major", "Minor", "Inactive")))
    # Create a data frame for plotting and set levels so order in the plot is Major, Minor, Inactive.
    #df5 <- data.frame(Type = factor(names(counts), levels = c("Major","Minor","Inactive")), Freq = counts)
    df5 <- data.frame(Type = factor(names(counts), levels = c("Major","Minor","Inactive")),
                  Freq = as.numeric(counts))
    # Create a bar plot with one bar stacked of gene counts for each promoter type.
    ggplot(df5, aes(x = 1, y = Freq, fill = Type)) +
      # Bar plot with black border.
      geom_col(width = .6, colour = "black") +
      # Add count labels on each stack segment; text position adjusted to be inside the corresponding segment
      geom_text(aes(label = comma(Freq)), position = position_stack(.8),
                colour = "white", size = 1.8, fontface = "bold") +
      scale_fill_manual(values = cols_type) +
      plot_theme +
      # Disable x axis related because only one bar: one graph for each cell line.
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), legend.position = "right") +
      labs(y = "Number of genes", x = NULL)
  }
  
  # Plot 6: Single/multiple promoter activity category -----------------------------------

  # A horizontally stacked bar plot showing the percentage of genes for each category.
  # x-axis is the percentage of active promoters.
  
    active_promoter_category <- function(x) {
      # Filter out inactive promoters.
      if (any(x %in% "Major")) {
        if (any(x %in% "Minor")) {
          # Multiple active promoters if has major and at least one minor.
          return("Multipromoter.Multiactive")
        } else if (any(x %in% "Inactive")) {
          # Single active promoter if has major but no minor. Multipromoter because it has at least one inactive.
          return("Multipromoter.Singleactive")
        } else {
          # Single active promoter if has major but no minor and no inactive.
          return("Singlepromoter.Singleactive")
        }
      }
      return("Inactive")
    }
    # Group all promoters by geneId. Each gene is now a list of promoter classes its promoters belong to.
    gene_split  <- split(rowData[[col_class]], rowData$geneId)
    # Apply the custom classification function to each gene’s promoter class vector, returning one of several activity category labels per gene.
    category_v  <- sapply(gene_split, active_promoter_category)
    # Convert the result into a frequency table of number of genes per category.
    df6 <- as.data.frame(table(Category = category_v), stringsAsFactors = FALSE)
    # Exclude Inactive category from the plot, as only active promoters are of interest.
    df6 <- df6[df6$Category != "Inactive", ]
    # Categorize the categories into factors with plot labels clearer for interpretation.
    df6$Category <- factor(df6$Category,
                          levels = c("Singlepromoter.Singleactive", "Multipromoter.Singleactive", "Multipromoter.Multiactive"),
                          labels  = c("Single active promoters\n(Single promoter genes)", "Single active promoters\n(Multi promoter genes)", "Multiple active promoters\n(Multi promoter genes)"))
    # The bar is 100% wide, partitioned into colored segments corresponding to different categories.
    p6 <- ggplot(df6, aes(x = 1, y = Freq / sum(Freq) * 100, fill = Category)) +
            # Stacked bar with black border; bar is horizontal.
            geom_col(colour = "black") + coord_flip() +
            # Apply theme and colors.
            scale_fill_manual(values = c("#851d41", "#db3056", "#ff6464")) +
            plot_theme +
            # Move legend to top and remove y axis related since there is only one bar and coordinates are flipped.
            theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), legend.position = "top",
              legend.title = element_blank(), legend.text = element_text(size = 7),
              legend.box = "horizontal"
            ) +
            # Legend is arranged in 3 rows for better visibility.
            guides(fill = guide_legend(nrow = 3)) +
            labs(y = "% of active promoters", x = NULL)
    
  # Graph from function of each graph and save---------------------------------------------------
  # Name file according to its cell line; save file to the directory named by its cell line.
  ggsave(plot_cat_percentage(), file = file.path(cl_dir, paste0(cl, "_promoter_activity_category_percentage.pdf")), width = 2.7, height = 2.9, units = "in", dpi = 300)
  ggsave(plot_activity_comparison(), file = file.path(cl_dir, paste0(cl, "_promoter_activity_category_comparison.pdf")), width = 2.5, height = 2.5, units = "in", dpi = 300)
  ggsave(plot_position_category(), file = file.path(cl_dir, paste0(cl, "_promoter_activity_position_category.pdf")), width = 3, height = 3, units = "in", dpi = 300)
  ggsave(plot_gene_corr(), file = file.path(cl_dir, paste0(cl, "_promoter_activity_geneexpression_correlation.pdf")), width = 3, height = 3, units = "in", dpi = 300)
  ggsave(plot_genewise_percentage(), file = file.path(cl_dir, paste0(cl, "_promoter_activity_category_percentage_genewise.pdf")), width = 2.7, height = 3, units = "in", dpi = 300)
  ggsave(p6, file = file.path(cl_dir, paste0(cl, "_promoter_activity_single_multiple_category.pdf")), width = 2.5, height = 2.5, units = "in", dpi = 300)

  message("[", cl, "] finished")
}

# Parallel execution ---------------------------------------
# Use multiple cores for different cell lines to speed up the plotting process.
# Set up parallel processing with one less than total cores to avoid overloading.
plan(multisession, workers = min(length(cell_lines), availableCores() - 1))
# apply function to each cell line
future_lapply(cell_lines, make_plots)
#lapply(cell_lines, make_plots)
#future_lapply(cell_lines, plot_activity_comparison)


  message("Plots 1-6 completed for all cell lines.")

# =======================
# Generate Plots that do not need to iterate through each cell line
# =======================

# Plot 7: promoter_activity_number_hist_all

# A bar plot showing the distribution of the number of promoters per gene across all genes.
# x axis is the number of promoters; y axis is the number of genes.
# Plot 7 includes all genes, while Plot 8 excludes genes with only one promoter.

promoterAnnotationData <- promoter_anno
# Get promoter-to-gene mapping table: each row links a promoterId to its geneId.
promoter_id_mapping <- promoterIdMapping(promoterAnnotationData)
# Remove duplicate promoter-gene pairs, keeping only unique pairs.
# Table count frequency each geneId appears, representing the number of promoters associated with each geneId.
# Create frequency table with each row representing geneId and its frequency (promoter counts).
gene_prom_table <- as.data.frame(table(unique(promoter_id_mapping[, c("promoterId", "geneId")])$geneId))
# Cap the frequency of promoters per gene at 11.
gene_prom_table$Freq[gene_prom_table$Freq > 11] <- 11
# Count how many genes fall into each promoter number category
plot_hist <- as.data.frame(table(gene_prom_table$Freq))

# Plot 7: promoter_activity_number_hist_all -----------------------------------

p7 <- ggplot(plot_hist, aes(x = Var1, y = Freq)) +
        geom_col(fill = "white", colour = "black") +
        # Add count labels on top of each bar
        geom_text(aes(label = Freq), vjust = -0.25, size = 2.5) +
        # Replace capped 11 as ">=11" in x axis labels
        scale_x_discrete(labels = c(1:10, ">=11")) +
        plot_theme +
        labs(x = "Promoter number", y = "Number of genes",
             title = "Promoter number distribution (all genes)") +
        # Adjust title size and centering
        theme(
        plot.title = element_text(size = 8, hjust = 0.5)
        )

# Plot 8: promoter_activity_number_hist_without1 -----------------------------------

# Remove first row (corresponds to promoter count == 1)
plot_hist2 <- plot_hist[-1, ]
p8 <- ggplot(plot_hist2, aes(x = Var1, y = Freq)) +
        geom_col(fill = "white", colour = "black") +
        # Add count labels on top of each bar
        geom_text(aes(label = Freq), vjust = -0.25, size = 2.5) +
        # Replace capped 11 as ">=11" in x axis labels
        scale_x_discrete(labels = c(2:10, ">=11")) +
        plot_theme +
        labs(x = "Promoter number", y = "Number of genes",
             title = "Promoter number distribution (>=2 promoters)") +
        # Adjust title size and centering
        theme(
        plot.title = element_text(size = 8, hjust = 0.5)
        )

# Run the plots and save them to the output directory
ggsave(p7, file = file.path(output_dir, "promoter_activity_number_hist_all.pdf"),      width = 3.5, height = 3, units = "in", dpi = 300)
ggsave(p8, file = file.path(output_dir, "promoter_activity_number_hist_without1.pdf"), width = 3.5, height = 3, units = "in", dpi = 300)

message("Plots 6-8 completed. All figures generated successfully.")


# Plot 9: t-SNE -----------------------------------

# A t-SNE plot visualizing promoter activity across all samples (every point is a sample).
# X axis is tSNE1; y axis is tSNE2.
# Ideally, samples from the same cell line cluster together, indicating similar promoter activity patterns.

message("Generating global t-SNE plot …")

# Extract promoter activity matrix from the assay
mat <- activity
# Filter out promoters with no activity across all samples.
#mat <- mat[rowSums(mat) > 0, ]
mat <- mat[rowSums(mat, na.rm = TRUE) > 0, ]
mat <- mat[complete.cases(mat), ]
# Transpose the matrix: now rows = samples (observation), columns = promoters (feature).
# tsne expects one observation per row.
tsne_df <- data.frame(t(mat))
# Use condition to match samples and convert to factor with specified levels for later coloring.
tsne_df$Sample <- factor(condition, levels = cell_lines)

# Set random seed for reproducibility
set.seed(40)
# Run t-SNE on expression matrix, excluding the Sample column.
# Y is a 2D matrix with t-SNE coordinates for each sample
#Y <- Rtsne(as.matrix(tsne_df[ , !names(tsne_df) %in% "Sample"]), perplexity = 1)$Y
#pca_res <- prcomp(tsne_df[ , !names(tsne_df) %in% "Sample"], scale.=TRUE)
#Y <- Rtsne(pca_res$x[, 1:2], perplexity=1)$Y
Y <- Rtsne(as.matrix(tsne_df[ , !names(tsne_df) %in% "Sample"]), perplexity = 1)$Y

# Define colors for each cell line and match by name.
cell_cols <- setNames(c("#ff1e56", "#ffac41", "#323232")[seq_along(cell_lines)],
                      cell_lines)

pdf(file.path(output_dir, "promoter_activity_tsne_plot.pdf"), width = 12/2.54, height = 10/2.54)
# Allow drawing outside the plot region for legend.
par(xpd = NA)
# Set shape, filled color, border color, size, and aspect ratio for the t-SNE plot.
# ,1 is tSNE1, x axis, the first dimension after dimentionally reduced to 2D.
plot(Y[,1], Y[,2], pch = 24, bg = cell_cols[as.character(tsne_df$Sample)], col = "black", cex = 1.4, asp = 1,
     xlab = "tSNE1", ylab = "tSNE2", main = "t-SNE plot (promoters active >=1 sample)")
# legend matches point shape, filled colors, and border color of cell lines.
legend("topright", title = "Cell Lines", legend = levels(tsne_df$Sample), pch = 24,
       pt.bg = cell_cols[levels(tsne_df$Sample)], col = "black", bty = "n", cex = 0.9)
dev.off()
message("t-SNE plot saved as promoter_activity_tsne_plot.pdf")

message("All plots (1-9) generated successfully.")

# Print warnings -----------------------------------
if (length(warnings()) > 0) {
  cat("\n=== WARNINGS ===\n")
  print(warnings())
}


# -- End of script ------------------------------------------------