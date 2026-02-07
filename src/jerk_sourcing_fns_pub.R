######################################################################
######################################################################
# Jerk tag - sourcing functions 
# by Jacey Van Wert
# Date: Dec 12 2025
######################################################################
######################################################################



##############################################################
# Load libraries
##############################################################

if (!require("pacman")) install.packages("pacman")

# Load packages and install if not present
pacman::p_load(
  here,
  tidyverse,
  janitor,
  bestNormalize,
  dtw,
  ape,
  brms,
  bayesplot,
  tidybayes,
  ggplot2,
  ggdist,
  ggridges,
  xts,
  circlize,
  dendextend,
  doParallel,
  viridis,
  corrplot,
  foreach,
  zoo,
  ggh4x,
  patchwork,
  emmeans,
  vegan
)

##############################################################
# Set path
##############################################################

# set baseline path
path <- "~/Documents/Github_repositories/fightr/"
setwd(path)


# where things will be saved 
savepath <- "~/Documents/Github_repositories/fightr/"



##############################################################
# Load functions
##############################################################
source(paste0(path, "functions/clusGap.ZS.R")) #for finding k - optimal number of clusters
source(paste0(path, "functions/jerk_fns.R")) #for suite of functions
source(paste0(path, "functions/compute_dtw_fn.R")) #for dtw function


##############################################################
# Load data
##############################################################

# Load jerk data
jerk <- read.csv(paste0(path,"data/all_jerk_2024_pub.csv")) 
jerk$id <- str_sub(jerk$fish_id, -3, -1)


# Load metadata 
meta <- read.csv(paste0(path, "data/salmon_2024_meta_pub.csv")) %>% 
  mutate(weight_kg_estim =  as.numeric(fl_cm) * as.numeric(girth_cm) * as.numeric(girth_cm) / 27120)
meta$id <- as.character(meta$id)


##############################################################
# Transform variables of interest
# transform() 
##############################################################

# select names of variables to normalize
totrans <- meta %>% 
  select(weight_kg_estim, blood_ph_post_angling, lactate) %>% 
  names()

# transform 
transformations <- transform(
  data = meta,
  vars_to_transform = totrans,
  save = TRUE
)


#add normalized variables back to meta
transformed_cols <- transformations$data %>%
  select(fish_id, ends_with(".trans"))
meta <- meta %>%
  left_join(transformed_cols, by = "fish_id")
meta$fish_id <- as.factor(meta$fish_id)

##############################################################
# Run DTW on raw, differential, high pass, multivariate models (12 methods)
# compute_dtw()
##############################################################

dtw_results <- compute_dtw(
  data = jerk,
  fish_id_col = "fish_id",
  x_col = "act_x",
  y_col = "act_y",
  z_col = "act_z",
  duration_col = "duration",
  event_id_col = "event_id",
  methods = c("raw", "differential", "multivariate", "highpass"),
  fig_output_dir = paste0(path,"figs/traces")
)

View(dtw_results)
saveRDS(dtw_results, paste0(path, "dtw_results/dtw_results.rds"))


##############################################################
# Run PCoA on the 12 DTW methods
# run_pcoa()
##############################################################

pcoa_results <- lapply(names(dtw_results$dtw_matrices), function(method_name) {
  message("Running PCoA for: ", method_name)
  run_pcoa(
    dtw_matrix = dtw_results$dtw_matrices[[method_name]],
    method_name = method_name,
    meta_data = meta,
    n_dims = 3
  )
})
names(pcoa_results) <- names(dtw_results$dtw_matrices)





##############################################################
# Clustering for 12 DTW methods
# hclust_analysis()
##############################################################

linkage_methods <- c("ward.D2")
hclust_results <- list()

for(linkage in linkage_methods) {
  for(method_name in names(dtw_results$dtw_matrices)) {
    
    result_name <- paste0(method_name, "_", linkage)
    message("Running hierarchical clustering: ", result_name)
    
    hclust_results[[result_name]] <- hclust_analysis(
      matrix = dtw_results$dtw_matrices[[method_name]],
      method_name = method_name,
      linkage_method = linkage
    )
    
    # show results for each  
    cat(sprintf("  -> k=%d, cophenetic=%.3f\n", 
                hclust_results[[result_name]]$k,
                hclust_results[[result_name]]$cophenetic_correlation))
  }
}

# clustering summary
# higher cophenetic correlation = better 
cluster_summary <- data.frame(
  method = sapply(hclust_results, `[[`, "method_name"),
  cophenetic_correlation = sapply(hclust_results, `[[`, "cophenetic_correlation"),
  clusters = sapply(hclust_results, `[[`, "k")
) %>%
  arrange(desc(cophenetic_correlation))
cat("\n=== Final Summary ===\n")
print(cluster_summary)



##############################################################
# Plot PCoA for 12 DTW methods
# plot_pcoa()
##############################################################

for(method_name in names(pcoa_results)) {
  message("Plotting PCoA for: ", method_name)
  
  plot_pcoa(
    pcoa_obj = pcoa_results[[method_name]],
    shape_by = "sp",
    cluster_obj = hclust_results[[paste0(method_name, "_ward.D2")]], #including clusters
    show_ellipses = TRUE,
    label_points = FALSE,
    save_plot = TRUE,
    pc_a = 1,
    pc_b = 2
  )
}

##############################################################
# Jerk summary metrics- calculate
# summarize_jerk()
##############################################################

jerk_results <- summarize_jerk(jerk, 
                               burst_threshold = 1.25, save = TRUE)


View(jerk_results)


##############################################################
# Jerk summary metrics- PCA
# run_jerk_pca()
##############################################################

jerk_matrix <- jerk_results$z_matrix #get matrix for PCA
jerk_pca_results <- run_jerk_pca(jerk_matrix, 
                                 n_dims = 3, save = TRUE)




##############################################################
# Jerk summary metrics- Clustering
# hclust_analysis()
# plot_dendro_with_covs()
##############################################################

# pc matrix
View(jerk_pca_results$pc_matrix)

# dist matrix from pc scores
pc_dist_matrix_jerk <- dist(jerk_pca_results$pc_matrix, method = "euclidean")
pc_dist_matrix_jerk <- as.matrix(pc_dist_matrix_jerk)

# run clustering fxn
hc_jerk_results <- hclust_analysis(
  matrix = pc_dist_matrix_jerk,
  method_name = "PC_scores",
  linkage_method = "ward.D2",
  return_gap_stats = TRUE)


# plot dendrogram with trait heatmap
jerk_summary <- jerk_results$summary

plot_dendro_with_covs(
  hc_obj = hc_jerk_results,
  meta_data = meta, 
  additional_data = jerk_summary,
  numb_labels =  c("fight_duration"),
  trait_vars = c("mean_actvsum", "cv_jerk", "early_intensity", "late_intensity", "fatigue", "total_effort", 
                 "fight_duration"),
  show_clusters = FALSE,
  color_by = NULL,
  save_plot = TRUE,
  circle_size = 0.14,
  width=6
)


 
# summary for each cluster
jerk_summary$cluster <- hc_jerk_results$clusters[jerk_summary$fish_id]
cluster_summaries <- jerk_summary %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_fight_duration = mean(fight_duration),
    sd_fight_duration = sd(fight_duration),
    mean_total_effort = mean(total_effort),
    sd_total_effort = sd(total_effort)
  )




##############################################################
# Visualize pearson correlations 
# correlo()
##############################################################

correlo(phys_data = meta,
        jerk_data = jerk_summary,
        vars = c("fl_cm","n_bursts", "fight_duration", "mean_actvsum", "total_effort",
                 "water_temp_c",
                 "lactate", "blood_ph_post_angling"))





##############################################################
# BRMS
# analyze_fish_fight()
##############################################################

# FIRST: check sample size to param ratio
n_fish <- 14
n_responses <- 2
n_predictors <- 3 + 3 # 3 PCos + 3 covs (4 predictors for jerk_summary)
n_params_per_response <- n_predictors + 1  # +1 for intercept
total_params <- n_responses * n_params_per_response
ratio <- n_fish / total_params
cat("Total parameters:", total_params, "\n")
cat("Ratio:", round(ratio, 2), "\n")
cat("Recommendation:", ifelse(ratio < 1, "Underpowered", "Adequate"), "\n")


# Run brms 
brms_results <- analyze_fish_fight(
  jerk_data = jerk_results$z_scores,
  dtw_matrices = dtw_results$dtw_matrices,
  meta_data = meta,
  response_vars = c("blood_ph_post_angling.trans", "lactate.trans"),
  covariates =  c("water_temp_c", "sp", "weight_kg_estim"),
  n_dims = 3
)

# All results
comparison_table <- brms_results$comparison
all_models <- brms_results$models
best_model <- brms_results$best_model
summary(brms_results$best_model$model) #beta and 95% CI 

saveRDS(brms_results, paste0(path, "/brms_outputs/brms_results.rds"))



##############################################################
# Forest plot
##############################################################

library(bayesplot)

# quick view regression coefficients for best model
mcmc_plot(results_temp$best_model$model, type = "intervals", 
          pars = "^b_",
          prob = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")

# models with sig.
mcmc_plot(results_temp$models$x_highpass$model, type = "intervals", 
          pars = "^b_",
          prob = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")
summary(results_temp$models$x_highpass$model) #beta and 95% CI 

# models with sig.
mcmc_plot(results_temp$models$y_diff$model, type = "intervals", 
          pars = "^b_",
          prob = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")
summary(results_temp$models$y_diff$model) #beta and 95% CI 


##############################################################
# Best output - final plot 
##############################################################


# plot best output:
plot_dendro_with_covs(
  hclust_results$y_highpass_ward.D2 , 
  meta_data = meta,
  additional_data = jerk_summary,
  trait_vars = c(
    "lactate","blood_ph_post_angling",
    "fight_duration"),
  color_by = NULL,
  numb_labels = c("fight_duration"),
  save_plot = TRUE,
  height = 7, 
  width = 4
)





