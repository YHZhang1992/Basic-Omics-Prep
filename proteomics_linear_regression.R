# Install needed packages if not already installed
install.packages(c("tidyverse", "broom", "pheatmap", "ggplot2"))
BiocManager::install("vsn")

# Load libraries
library(tidyverse)
library(broom)
library(vsn)
library(pheatmap)
library(ggplot2)

# Load proteomics data (rows = samples, columns = proteins)
proteins_raw <- read.csv("somascan_proteins.csv", row.names = 1)

# Load sample metadata
metadata <- read.csv("metadata.csv", row.names = 1)

# Ensure matching sample order
proteins_raw <- proteins_raw[rownames(metadata), ]
stopifnot(identical(rownames(proteins_raw), rownames(metadata)))

# Convert to matrix and apply VSN
vsn_fit <- vsn2(as.matrix(proteins_raw))
proteins_norm <- predict(vsn_fit, newdata = as.matrix(proteins_raw))

# Convert normalized data to data frame
proteins <- as.data.frame(proteins_norm)

# Example phenotype: binary or continuous
metadata$disease_status <- as.factor(metadata$disease_status)  # categorical
# metadata$FEV1 <- as.numeric(metadata$FEV1)                   # continuous

# Function to run regression for each protein
run_lm <- function(protein, outcome, covariates = NULL) {
  df <- cbind(metadata, protein_expr = proteins[[protein]])
  formula_str <- paste("protein_expr ~", paste(c(outcome, covariates), collapse = " + "))
  fit <- lm(as.formula(formula_str), data = df)
  tidy(fit) %>% filter(term == outcome) %>% mutate(protein = protein)
}

# Define outcome and covariates
outcome_var <- "disease_status"
covariates <- c("age", "sex", "smoking_status")  # Optional

# Apply to all proteins
results_list <- lapply(colnames(proteins), run_lm, outcome = outcome_var, covariates = covariates)
results_df <- bind_rows(results_list)

results_df <- results_df %>%
  arrange(p.value) %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

sig_results <- filter(results_df, p_adj < 0.05)

annotations <- read.csv("somascan_annotations.csv")  # includes SOMAmer ID and gene symbol
sig_results_annot <- left_join(sig_results, annotations, by = c("protein" = "SOMAmer_ID"))

ggplot(results_df, aes(x = estimate, y = -log10(p_adj))) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = filter(results_df, p_adj < 0.01 & abs(estimate) > 1),
            aes(label = protein), vjust = 1.5, size = 3) +
  theme_minimal() +
  labs(x = "Effect Size (Beta)", y = "-log10 Adjusted p-value")

top_proteins <- sig_results %>% arrange(p_adj) %>% head(30) %>% pull(protein)
pheatmap(scale(proteins[, top_proteins]),
         annotation_col = metadata[, c(outcome_var, covariates)],
         main = "Top 30 Differential Proteins")

write.csv(results_df, "all_linear_regression_results.csv", row.names = FALSE)
write.csv(sig_results_annot, "significant_proteins_annotated.csv", row.names = FALSE)
