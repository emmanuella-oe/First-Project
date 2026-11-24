
# ============================================================
# BINF*6210 - Assignment 2
# Trichinella spiralis vs. Trichinella britovi:
# Transposable Element Sequence Divergence
#
# Primary author: Emmanuella
# Secondary contributor: Anusha
#
# Original analysis:
# TE filtering, IQR outlier assessment, Wilcoxon tests,
# and comparison of divergence across shared TE classes.
#
# proposed improvements:
# 1. Correct filtering and add a reusable data import function.
# 2. Implement class-specific outlier detection.
# 3. Improve statistical summaries and visualizations.
# ============================================================


# 01 Load libraries and datasets -----------------------------

library(tidyverse)

# Run the script from the First-Project RStudio project root.
data_dir <- "2026_11_02_Assignment_1/data"

# Preserve original non-TE exclusions.
ls_non_te <- c(
  "Simple_repeat",
  "Low_complexity",
  "Satellite",
  "snRNA",
  "rRNA",
  "ARTEFACT",
  "Unknown",
  "simple_repeat"
)


# EDIT 1: Reusable data preparation function -----------------

read_te_data <- function(filename, species_name) {
  
  df <- read_csv(
    filename,
    show_col_types = FALSE
  )
  
  required_columns <- c("matching_class", "perc_div")
  
  if (!all(required_columns %in% names(df))) {
    stop("Missing required columns in: ", filename)
  }
  
  df |>
    mutate(
      annotation_id = row_number(),
      species = species_name
    ) |>
    filter(
      !is.na(matching_class),
      !matching_class %in% ls_non_te,
      is.finite(perc_div)
    )
}

# Import both species with the same function.
t_spiralis <- read_te_data(
  file.path(data_dir, "GCA_001447595.1.csv"),
  "T. spiralis"
)

t_britovi <- read_te_data(
  file.path(data_dir, "GCA_001447585.1.csv"),
  "T. britovi"
)

# Combine data for comparative analyses.
te_combined <- bind_rows(
  t_spiralis,
  t_britovi
)

# Check annotation counts after filtering.
te_combined |>
  count(species, name = "n_annotations") |>
  print()


# 02 Class-specific outlier analysis -------------------------

# EDIT 2: Reusable class-specific IQR function.
# The IQR is calculated separately for each species and
# TE matching class.
#
# Classes with fewer than 5 observations are not flagged.
# Outliers are retained rather than automatically removed.

flag_class_outliers <- function(df) {
  
  df |>
    group_by(species, matching_class) |>
    mutate(
      class_n = n(),
      q1 = quantile(perc_div, 0.25),
      q3 = quantile(perc_div, 0.75),
      class_iqr = q3 - q1,
      lower_bound = q1 - 1.5 * class_iqr,
      upper_bound = q3 + 1.5 * class_iqr,
      iqr_outlier = class_n >= 5 &
        (perc_div < lower_bound |
           perc_div > upper_bound)
    ) |>
    ungroup()
}

# Apply to combined TE data.
te_combined <- flag_class_outliers(te_combined)

# Summarize outlier flags by species and matching class.
outlier_summary <- te_combined |>
  group_by(species, matching_class) |>
  summarise(
    n_annotations = n(),
    n_outliers = sum(iqr_outlier),
    percent_outliers = 100 * mean(iqr_outlier),
    .groups = "drop"
  ) |>
  arrange(desc(n_outliers))

print(outlier_summary, n = Inf)

# Overall summary of outlier flags.
te_combined |>
  group_by(species) |>
  summarise(
    total_annotations = n(),
    flagged_outliers = sum(iqr_outlier),
    .groups = "drop"
  ) |>
  print()

# Visualize distributions with class-specific IQR boxplots.
# For readability, include classes with at least
# five annotations in both species.

shared_class_counts <- te_combined |>
  count(matching_class, species) |>
  pivot_wider(
    names_from = species,
    values_from = n,
    values_fill = 0
  )

eligible_classes <- shared_class_counts |>
  filter(
    `T. spiralis` >= 5,
    `T. britovi` >= 5
  ) |>
  pull(matching_class)

plot_data <- te_combined |>
  filter(matching_class %in% eligible_classes)

outlier_plot <- ggplot(
  plot_data,
  aes(
    x = matching_class,
    y = perc_div,
    fill = species
  )
) +
  geom_boxplot(
    position = position_dodge2(width = 0.8),
    outlier.alpha = 0.5
  ) +
  coord_flip() +
  labs(
    title = "TE Divergence Across Shared Matching Classes",
    subtitle = "Class-specific distributions in two Trichinella species",
    x = "TE Matching Class",
    y = "Percent Sequence Divergence",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

print(outlier_plot)


# 03 Wilcoxon tests by matching class ------------------------

# EDIT 3: Reusable statistical summary function.
#
# The primary analysis retains all valid annotations.
# Tests are only conducted for classes with at least five
# annotations per species and variable divergence values.

compare_te_classes <- function(df) {
  
  df |>
    group_by(matching_class) |>
    group_modify(~ {
      
      x <- .x$perc_div[
        .x$species == "T. spiralis"
      ]
      
      y <- .x$perc_div[
        .x$species == "T. britovi"
      ]
      
      can_test <- (
        length(x) >= 5 &&
          length(y) >= 5 &&
          n_distinct(c(x, y)) > 1
      )
      
      tibble(
        n_spiralis = length(x),
        n_britovi = length(y),
        
        median_spiralis = if (length(x) > 0) {
          median(x)
        } else {
          NA_real_
        },
        
        median_britovi = if (length(y) > 0) {
          median(y)
        } else {
          NA_real_
        },
        
        median_difference = if (
          length(x) > 0 && length(y) > 0
        ) {
          median(x) - median(y)
        } else {
          NA_real_
        },
        
        p_value = if (can_test) {
          wilcox.test(
            x, y,
            alternative = "two.sided",
            paired = FALSE,
            exact = FALSE
          )$p.value
        } else {
          NA_real_
        }
      )
    }) |>
    ungroup() |>
    mutate(
      p_adjusted = p.adjust(
        p_value,
        method = "BH"
      ),
      significant = !is.na(p_adjusted) &
        p_adjusted < 0.05
    ) |>
    arrange(p_adjusted)
}

# Primary analysis: include all valid TE annotations.
wilcox_results <- compare_te_classes(
  te_combined
)

print(wilcox_results, n = Inf)

# Show only classes eligible for statistical testing.
tested_results <- wilcox_results |>
  filter(!is.na(p_value))

print(tested_results, n = Inf)


# 04 Visualize statistical comparisons ----------------------

# Original-style comparison, now limited to tested classes.
comparison_plot <- ggplot(
  te_combined |>
    filter(matching_class %in%
             tested_results$matching_class),
  aes(
    x = matching_class,
    y = perc_div,
    fill = species
  )
) +
  geom_boxplot(
    outlier.alpha = 0.3
  ) +
  coord_flip() +
  labs(
    title = "Percent Divergence in Shared TE Classes",
    subtitle = "T. spiralis and T. britovi",
    x = "TE Matching Class",
    y = "Percent Sequence Divergence",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

print(comparison_plot)


# New visualization: direction and magnitude of
# the difference between species' median divergence.
#
# Positive values: higher median in T. spiralis.
# Negative values: higher median in T. britovi.
#
# Significance is based on BH-adjusted Wilcoxon p-values.
# The median difference is descriptive, not the
# Wilcoxon test statistic.

median_difference_plot <- ggplot(
  tested_results,
  aes(
    x = median_difference,
    y = reorder(
      matching_class,
      median_difference
    ),
    color = significant
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(size = 3) +
  labs(
    title = "Differences in Median TE Divergence",
    subtitle = paste(
      "T. spiralis minus T. britovi;",
      "BH-adjusted Wilcoxon results"
    ),
    x = "Difference in Median Percent Divergence",
    y = "TE Matching Class",
    color = "Adjusted p < 0.05"
  ) +
  theme_minimal()

print(median_difference_plot)


# 05 Sensitivity analysis -----------------------------------

# Compare the main analysis with one that excludes
# IQR-flagged observations.
#
# This is a sensitivity check, not automatic evidence
# that the flagged TE annotations are erroneous.

wilcox_without_outliers <- te_combined |>
  filter(!iqr_outlier) |>
  compare_te_classes()

sensitivity_summary <- wilcox_results |>
  select(
    matching_class,
    p_adjusted,
    significant
  ) |>
  rename(
    adjusted_p_all = p_adjusted,
    significant_all = significant
  ) |>
  left_join(
    wilcox_without_outliers |>
      select(
        matching_class,
        p_adjusted,
        significant
      ) |>
      rename(
        adjusted_p_filtered = p_adjusted,
        significant_filtered = significant
      ),
    by = "matching_class"
  )

print(sensitivity_summary, n = Inf)


# 06 Final summary ------------------------------------------

cat("\n--- Analysis summary ---\n")

cat(
  "Total TE annotations:",
  nrow(te_combined), "\n"
)

cat(
  "Classes tested:",
  nrow(tested_results), "\n"
)

cat(
  "Significant classes after BH correction:",
  sum(tested_results$significant), "\n"
)

cat(
  "Class-specific IQR flags:",
  sum(te_combined$iqr_outlier), "\n"
)

# End of revised analysis.

