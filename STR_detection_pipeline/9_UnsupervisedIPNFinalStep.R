# Required libraries
library(tidyverse)
library(cluster)
library(factoextra)
library(mclust)
library(NbClust)
library(dendextend)
library(readxl)

set.seed(1234)

filepath <- /file/path
sheet_name <- "sheet_name"

data <- read_excel(filepath, sheet = sheet_name)

# Exclude controls with dx
data <- data[is.na(data$exclude_reason) | data$exclude_reason == "", ]

# Select and scale features
features <- c("percent_reads", "max_str_length", "mean_str_length")
features <- c("percent_reads", "mean_str_length")
data_subset <- data[features]
scaled_data <- scale(data_subset)

data$casecontrol <- ifelse(grepl("PNRR", data$sample_name), "Case", "Control")
print(table(data$casecontrol))

# ------------------------------
# 1. K-means Clustering
# ------------------------------

# Determine optimal k
# Elbow method
wss <- sapply(1:10, function(k) {
  kmeans(scaled_data, centers = k, nstart = 25)$tot.withinss
})
plot(1:10, wss, type = "b",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares",
     main = "K-means Elbow Method")

# Silhouette method
fviz_nbclust(scaled_data, kmeans, method = "silhouette")

# Perform k-means with chosen k
k <- 2
kmeans_result <- kmeans(scaled_data, centers = k, nstart = 25)

# Show k-means grouping results with case/control status
kmeans_results_df <- data.frame(
  Sample = data$sample_name,
  Status = data$casecontrol,
  Cluster = kmeans_result$cluster
) %>% arrange(Cluster)
print("\nK-means Clustering Results:")
print(kmeans_results_df)

# Print cluster composition
print("\nCluster composition by case/control status:")
print(table(Cluster = kmeans_result$cluster, Status = data$casecontrol))

# Visualize k-means clusters
fviz_cluster(kmeans_result, data = scaled_data,
             main = "K-means Cluster Plot",
             geom = "point",
             ggtheme = theme_minimal()) +
  geom_text(aes(label = paste0(data$casecontrol), 
                color = data$casecontrol),
            hjust = -0.1, size = 3) +
  scale_color_manual(values = c("Case" = "#EE7733", "Control" = "#0077BB")) +
  guides(color = guide_legend(title = "Status"))

# ------------------------------
# 2. Hierarchical Clustering
# ------------------------------

# Calculate distance matrix
dist_matrix <- dist(scaled_data)
hc <- hclust(dist_matrix, method = "ward.D2")

# Plot full dendrogram with colored labels for case/control
labels_colors <- ifelse(data$casecontrol == "Case", "#EE7733", "#0077BB")
plot(hc, labels = paste0(data$casecontrol),
     main = "Hierarchical Clustering Dendrogram",
     cex = 0.6, hang = -1)

# Cut tree and show groups
k_hc <- k  # Use same k as k-means
hc_clusters <- cutree(hc, k = k_hc)
hc_results_df <- data.frame(
  Sample = data$sample_name,
  Status = data$casecontrol,
  Cluster = hc_clusters
) %>% arrange(Cluster)
print("\nHierarchical Clustering Results:")
print(hc_results_df)

print("\nHierarchical cluster composition by case/control status:")
print(table(Cluster = hc_clusters, Status = data$casecontrol))

# Visualize hierarchical clusters in 2D
fviz_cluster(list(data = scaled_data, cluster = hc_clusters),
             main = "Hierarchical Cluster Plot",
             geom = "point",
             ggtheme = theme_minimal()) +
  geom_text(aes(label = paste0(data$casecontrol),
                color = data$casecontrol),
            hjust = -0.1, size = 3) +
  scale_color_manual(values = c("Case" = "#EE7733", "Control" = "#0077BB")) +
  guides(color = guide_legend(title = "Status"))

# ------------------------------
# 3. GMM Clustering
# ------------------------------

# Perform GMM
gmm_result <- Mclust(scaled_data, G=1:5)

# Show GMM results with case/control status
gmm_results_df <- data.frame(
  Sample = data$sample_name,
  Status = data$casecontrol,
  Cluster = gmm_result$classification,
  Uncertainty = round(1 - apply(gmm_result$z, 1, max), 3)
) %>% arrange(Cluster)
print("\nGMM Clustering Results (with uncertainty):")
print(gmm_results_df)

print("\nGMM cluster composition by case/control status:")
print(table(Cluster = gmm_result$classification, Status = data$casecontrol))

# Visualize GMM clusters
fviz_cluster(list(data = scaled_data, cluster = gmm_result$classification),
             main = "GMM Cluster Plot",
             geom = "point",
             ggtheme = theme_minimal()) +
  geom_text(aes(label = paste0(data$casecontrol),
                color = data$casecontrol),
            hjust = -0.1, size = 3) +
  scale_color_manual(values = c("Case" = "#EE7733", "Control" = "#0077BB")) +
  guides(color = guide_legend(title = "Status"))

# ------------------------------
# 4. PCA Visualization
# ------------------------------

# Perform PCA
pca_result <- prcomp(scaled_data)

# Create PCA plot data frame with case/control status
pca_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Status = data$casecontrol  # Add the case/control status
)

# Create the plot with colored points and labels
ggplot(pca_data, aes(x = PC1, y = PC2, color = Status)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Case" = "#EE7733", "Control" = "#0077BB")) +  # Using the same colors as your clustering plots
  theme_minimal() +
  labs(title = "PCA Plot of Samples",
       x = paste0("First Principal Component (", var_explained[1], "%)"),
       y = paste0("Second Principal Component (", var_explained[2], "%)"))

# Show variance explained by each PC
var_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 2)
print("\nVariance explained by each PC:")
print(paste("PC1:", var_explained[1], "%"))
print(paste("PC2:", var_explained[2], "%"))

kmeans_all <- data.frame(
  Sample = data$sample_name,
  Status = data$casecontrol,
  KMeans = kmeans_result$cluster
) %>%
  left_join(data[c("sample_name", features)], by = c("Sample" = "sample_name"))

hc_all <- data.frame(
  Sample = data$sample_name,
  Status = data$casecontrol,
  Hierarchical = hc_clusters
) %>%
  left_join(data[c("sample_name", features)], by = c("Sample" = "sample_name"))

gmm_all <- data.frame(
  Sample = data$sample_name,
  Status = data$casecontrol,
  GMM = gmm_result$classification
) %>%
  left_join(data[c("sample_name", features)], by = c("Sample" = "sample_name"))

kmeans_metrics <- kmeans_all %>%
  group_by(KMeans) %>%
  summarize(
    avg_percent_reads = mean(percent_reads),
    avg_mean_str_length = mean(mean_str_length)
  )

kmeans_biallelic_cluster <- kmeans_metrics %>%
  mutate(total_score = avg_percent_reads + avg_mean_str_length) %>%
  slice_max(total_score) %>%
  pull(KMeans)

hc_metrics <- hc_all %>%
  group_by(Hierarchical) %>%
  summarize(
    avg_percent_reads = mean(percent_reads),
    avg_mean_str_length = mean(mean_str_length)
  )

hc_biallelic_cluster <- hc_metrics %>%
  mutate(total_score = avg_percent_reads + avg_mean_str_length) %>%
  slice_max(total_score) %>%
  pull(Hierarchical)

gmm_metrics <- gmm_all %>%
  group_by(GMM) %>%
  summarize(
    avg_percent_reads = mean(percent_reads),
    avg_mean_str_length = mean(mean_str_length)
  )

gmm_biallelic_cluster <- gmm_metrics %>%
  mutate(total_score = avg_percent_reads + avg_mean_str_length) %>%
  slice_max(total_score) %>%
  pull(GMM)

final_results <- data.frame(
  sample_name = data$sample_name,
  kmeans = ifelse(kmeans_result$cluster == kmeans_biallelic_cluster, "Biallelic", "Monoallelic"),
  hierarchical = ifelse(hc_clusters == hc_biallelic_cluster, "Biallelic", "Monoallelic"),
  GMM = ifelse(gmm_result$classification == gmm_biallelic_cluster, "Biallelic", "Monoallelic"),
  other_detected_patterns = data$other_detected_patterns
)

final_results_2 <- final_results
final_results_2$WGS_new <- ifelse(final_results$hierarchical == "Biallelic" | 
                                    final_results$GMM == "Biallelic", 
                                  "Biallelic", "Monoallelic")

write.csv(clustering_comparison, "clustering_comparison.csv", row.names = FALSE)
