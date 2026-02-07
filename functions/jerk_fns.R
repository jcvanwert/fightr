################################################################################
# DTW FUNCTIONS
# Author: Jacey Van Wert 
# Date: Dec 12 2025
# This file has a variety of functions useful for the jerk/accel data analysis
# These functions will be integrated into a new R package (In progress)
################################################################################


# Load packages and install if not present
if (!require("pacman")) install.packages("pacman")

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

# Function for transformation
transform <- function(data,
                      vars_to_transform,
                      save = TRUE,
                      loo=TRUE,
                      suffix = ".trans"){
  
  require(bestNormalize)
  
  transformations <- list()
  
  # create save directory if it doesn't exist
  if(save) {
    transform_dir <- paste0(savepath, "transformations")
    if(!dir.exists(transform_dir)) {
      dir.create(transform_dir, recursive = TRUE)
    }
  }
  
  # loop through variables
  for(var in vars_to_transform) {
    cat("Transforming:", var, "\n")
    bn <- bestNormalize(data[[var]], loo = loo, na.rm = TRUE) # apply bestNormalize
    transformations[[var]] <- bn # store 
    data[[paste0(var, suffix)]] <- bn$x.t  # add transformed variable to dataframe

    
    # save transformation object
    if(save) {
      save_path <- file.path(transform_dir, paste0(var, "_transform.rds"))
      saveRDS(bn, save_path)
      cat("  Transformation:", class(bn$chosen_transform)[1], "\n")
      message("  Saved: ", save_path)
    } else {
      cat("  Transformation:", class(bn$chosen_transform)[1], "\n")
    }
    cat("\n")
  }
  
  # return list with updated data and transformations
  return(list(
    data = data,
    transformations = transformations
  ))
  
}


# Function for correlogram
correlo <- function(phys_data = NULL, jerk_data = NULL, 
                    vars = NULL, 
                    colors = c("#baae00", "#ded471", "#FFFFFF", "#dab8e2", "#a962b6"),
                    save_plot = TRUE,
                    width = 8, height = 8,
                    output_dir = paste0(savepath,"figs/correlogram")){
  
  require(corrplot)
  
  # handle data input
  if(!is.null(phys_data) && !is.null(jerk_data)) {
    final_data <- left_join(phys_data, jerk_data, by = "fish_id")
  } else if(!is.null(phys_data)) {
    final_data <- phys_data
  } else if(!is.null(jerk_data)) {
    final_data <- jerk_data
  } else {
    stop("Must provide 'phys_data', 'jerk_data', or both")
  }
  
  if(!is.null(vars)) {
    final_data <- final_data[, vars, drop = FALSE]  
  } else {
    # use all numeric columns if vars not specified
    final_data <- final_data %>% select(where(is.numeric))
  }
  

  
  # calculate correlation
  corr <- cor(final_data, use = "complete.obs")
  
  # rename labels
  new_names <- colnames(corr)
  new_names <- tools::toTitleCase(new_names)  # capitalize all
  new_names[new_names == "Blood_ph_post_angling"] <- "Blood pH"
  new_names[new_names == "Blood_ph"] <- "Blood pH"
  new_names[new_names == "Weight_kg_estim"] <- "Weight"
  new_names[new_names == "Fl_cm"] <- "Fork length"
  new_names[new_names == "AUC"] <- "Recovery\ncost"
  new_names[new_names == "Mean_actvsum"] <- "Mean ActVSum"
  new_names[new_names == "Water_temp_c"] <- "Water temperature"
  new_names <- gsub("_", " ", new_names)  # Replace underscores with spaces
  
  colnames(corr) <- new_names
  rownames(corr) <- new_names
  
  # create color palette
  col <- colorRampPalette(colors)
  
  if(save_plot) {
    if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- "correlogram.png"
    filepath <- file.path(output_dir, filename)
    png(filepath, width = width, height = height, units = "in", res = 300)
    
    corrplot(corr, 
             method="color", 
             # order="hclust", #automatic ordering by similarities
             col=col(200),  
             type="upper",  
             addCoef.col = "black", 
             tl.col="black", 
             tl.srt=45, 
             diag=FALSE)
  
    
    dev.off()
    message("Saved: ", filepath)
    
  } else {
    if(dev.cur() > 1) dev.off()
    dev.new()
    
    corrplot(corr, 
             method = "color",
             col = col(200),
             type = "upper",
             addCoef.col = "black",
             tl.col = "black", 
             tl.srt = 45,
             diag = FALSE)
  }
  
  # return correlation matrix invisibly
  invisible(corr)
}





# Functions for hierarchical analysis and plotting
hclust_analysis <- function(matrix, method_name, linkage_method = NULL, k =NULL,
                            return_gap_stats = FALSE #if want to calculate k, otherwise will provide
) {
  
  UseCores <- detectCores() - 4
  options('mc.cores' = UseCores)
  registerDoParallel(UseCores)
  
  # convert to distance object
  dist <- as.dist(matrix)
  
  # if no k provided, find optimal k using gap stat
  if(is.null(k)){
    cat("Finding optimal k with gap statistic...\n")
    
    hcluster<- function(x, k_val){
      cutree(hclust(as.dist(x), method = linkage_method), k = k_val)
    }
    
    k_max <- 15
    n_samples <- attr(dist, "Size")
    K.max <- min(n_samples - 1, k_max)
    
    if(K.max < 2){
      stop("K of 2 or less")
    }
    
    #calc gap statistic (need to source clusGap.ZS fxn)
    cG <- clusGap.ZS(as.matrix(dist),
                     FUNcluster = hcluster,
                     K.max = K.max,
                     B = 100,
                     cores = UseCores)
    
    # optimal k
    k <- .maxSE(f = cG$Tab[, "gap"],
                SE.f = cG$Tab[, "SE.sim"],
                method = "Tibs2001SEmax",
                SE.factor = 1)
  }
  
  
  # hierarchical clustering
  hc <- hclust(dist, method = linkage_method)
  
  clusters <- cutree(hc, k = k)
  
  # calculate cophenetic correlation
  cophenetic_cor <- cor(dist, cophenetic(hc))
  
  # results
  results <- list(
    hclust_obj = hc,
    distance_matrix = as.matrix(dist),
    clusters = clusters,
    k = k,
    method_name = method_name,
    linkage = linkage_method,
    cophenetic_correlation = cophenetic_cor
  )
  
  if(return_gap_stats && exists("cG")) {
    results$gap_stats <- cG
  }
  
  return(results)
}


plot_dendro_with_covs <- function(hc_obj, 
                                  meta_data, 
                                  additional_data = NULL,
                                  trait_vars = NULL,
                                  show_clusters = TRUE,
                                  color_by = "cluster",
                                  cluster_colors = c("#62A8BF", "#8A6BB0", "#6B8E9F", "#67A38A", 
                                                              "#B7D3F2", "#AFAFDC", "#DCAB6B", "#F4B9B2"), 
                                  numb_labels = TRUE, #TRUE, FALSE, or vars c()
                                  circle_size = 0.12,
                                  filename = NULL,
                                  shorten_labels = TRUE,  #if want to shorten fish id labels
                                  label_start = 6,    #if do want to shorten, then by where
                                  output_dir = paste0(savepath,"figs/clusters"),
                                  save_plot = TRUE, 
                                  height = 7, 
                                  width = 4) {
  
  require(dendextend)
  require(circlize)
  require(viridis)
  
  # reset graphics 
  if(!save_plot) {
    if(dev.cur() > 1) dev.off()
    dev.new()
  }
  
  # create dendrogram 
  dend <- as.dendrogram(hc_obj$hclust_obj)
  if(shorten_labels) {
    labels(dend) <- substring(labels(dend), label_start)
  }
  
  # match metadata to dendrogram order
  fish_order <- hc_obj$hclust_obj$labels[hc_obj$hclust_obj$order]
  meta_subset <- meta_data[match(fish_order, meta_data$fish_id), ]
  
  # additional data if provided
  if(!is.null(additional_data)) {
    meta_subset <- left_join(meta_subset, additional_data, by = "fish_id")
  }
  
  #add cluster info to meta
  if(!is.null(hc_obj$clusters)) {
    meta_subset$cluster <- hc_obj$clusters[match(meta_subset$fish_id, names(hc_obj$clusters))]
  }
  
  if(!is.null(color_by)) {
    if(color_by == "cluster" && !is.null(hc_obj$clusters)) {
      color_vals <- as.factor(meta_subset$cluster)
      n_groups <- nlevels(color_vals)
      label_colors <- cluster_colors[as.numeric(color_vals)]
    } else {
      color_vals <- as.factor(meta_subset[[color_by]])
      n_groups <- nlevels(color_vals)
      label_colors <- colorRampPalette(c("black", "black"))(n_groups)[as.numeric(color_vals)]
    }
    labels_colors(dend) <- label_colors
  }
  
  
  # open device
  if(save_plot) {
    if(is.null(filename)) {
      filename <- paste0(hc_obj$method_name, "_", hc_obj$linkage, "_dendrogram_traits.png")
    } else if(!grepl("\\.(png|pdf|jpg|jpeg)$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".png")  # add extension if missing
    }
    
    if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    filepath <- file.path(output_dir, filename)
    png(filepath, width = width, height = height, units = "in", res = 300)
  }
  
  # layout
  graphics::layout(matrix(c(1, 2,    # row 1: dendrogram, traits
                  3, 4),   # row 2: blank, legend
                nrow = 2, byrow = TRUE), 
         widths = c(0.35, 0.65), 
         heights = c(0.9, 0.08))
  
  
  ##### panel 1: dendrogram #####
  #############################  
  par(mar = c(0.1, 4, 6, 2))
  plot(dend, 
       horiz = TRUE,
       main = NULL,
       xlab = "Distance",
       ylab = "",
       yaxt = "n")
  
  ##### panel 2: heatmaps  #####
  #############################  
  par(mar = c(0.1, 5, 6, 2))
  
  n_fish <- nrow(meta_subset)
  n_traits <- length(trait_vars)
  
  # standardize trait vals
  trait_matrix <- sapply(trait_vars, function(var) {
    vals <- meta_subset[[var]]
    (vals - min(vals, na.rm = TRUE)) / (max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE))
  })
  
  # color palette
  mako_colors <- rev(mako(100))[1:70] 
  col_fun <- colorRamp2(seq(0, 1, length.out = 70), mako_colors)
  
  
  # positioning of circles
  spacing <- 0.3  #smaller = closer together
  x_positions <- (1:n_traits) * spacing
  
  # empty plot
  plot(1, type = "n", 
       xlim = c(0, max(x_positions) + spacing),
       ylim = c(0.5, n_fish + 0.5),
       xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n",
       main = NULL,
       bty = "n")
  
  
  
  trait_labels <- c(
    blood_ph_post_angling = "Blood pH",
    fat1        = "Fat",
    weight_kg_estim       = "Weight",
    glucose = "Glucose",
    lactate        = "Lactate",
    potassium = "Potassium",
    sodium = "Sodium",
    testosterone = "Testosterone",
    cortisol = "Cortisol",
    chloride = "Chloride",
    osmolality = "Osmolality",
    estradiol = "Estradiol",
    ht_avg = "Hematocrit",
    n_bursts = "N bursts",
    water_temp_c = "Water temp (°C)",
    fight_duration = "Duration (s)",
    total_effort = "Total effort",
    mean_actvsum = "Mean ActVSum",
    sd_jerk = "SD ActVSum",
    cv_jerk = "CV ActVSum",
    mean_sv = "Mean SV",
    sd_sv = "SD SV",
    cv_sv = "CV SV",
    early_intensity = "Early intensity",
    late_intensity = "Late intensity",
    fatigue = "Fatigue"
  )
  
  
  # add trait column labels
  text(x = x_positions, y = n_fish + 0.8, 
       labels = trait_labels[trait_vars], 
       srt = 45, adj = 0, xpd = TRUE, cex = 1.0)
  
  # draw circles
  for(i in 1:n_fish) {
    for(j in 1:n_traits) {
      std_val <- trait_matrix[i, j]
      actual_val <- meta_subset[[trait_vars[j]]][i]
      
      if(!is.na(std_val)) {
        symbols(x = x_positions[j], y = i, 
                circles = circle_size,
                inches = FALSE, 
                bg = col_fun(std_val),
                fg = "black",
                add = TRUE)
        
        show_number <- if(is.logical(numb_labels)) {
          numb_labels  # if TRUE or FALSE for all
        } else if(is.character(numb_labels)) {
          trait_vars[j] %in% numb_labels  #only if variable name is in the vector
        } else {
          FALSE
        }
        
        if(show_number) {
          text_color <- ifelse(std_val > 0.6, "white", "black")
          
          label_text <- ifelse(actual_val > 99,
                               round(actual_val, 0),
                               round(actual_val, 1))
          
          text(x = x_positions[j], y = i,
               labels = label_text,
               cex = 0.8,
               col = text_color)
        }
        
        
      }
    }
  }
  
  
  
  ##### panel 3: blank    #####
  #############################  
  par(mar = c(0,0,0,0))
  plot.new()
  
  
  ##### panel 4: legend  #####
  #############################  
  
  par(mar = c(0.5, 0.5, 0.01, 2))
  
  legend_vals <- seq(0, 1, length.out = 70)
  plot(1, type = "n", 
       xlim = c(0, 70), ylim = c(0, 2),
       xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n",
       bty = "n")
  
  bar_start <- 10  
  bar_end <- 60    
  
  for(i in bar_start:(bar_end-1)) {
    rect(i-1, 0.5, i+1, 1.5,  
         col = col_fun(legend_vals[round((i - bar_start) / (bar_end - bar_start) * 70)]), 
         border = NA)
  }
  
  # black border
  rect(bar_start, 0.5, bar_end, 1.5, border = "black", lwd = 1)
  
  # labels 
  text(x = c(bar_start, (bar_start + bar_end)/2, bar_end), y = 0.2, 
       labels = c("Low", "Med", "High"),
       adj = 0.5, xpd = TRUE, cex = 0.9)
  
  
  
  if(save_plot) {
    dev.off()
    message("Saved: ", filepath)
  }
  
  # reset layout after plotting
  graphics::layout(1)
}


# Functions for PCoA 
run_pcoa <- function(dtw_matrix, method_name, meta_data, n_dims=3,
                     output_dir = paste0(savepath,"figs/raw")) {
  
  # convert to distance object
  dtw_dist <- as.dist(dtw_matrix)
  
  # run PCoA
  pcoa_result <- pcoa(dtw_dist)
  
  # extract coordinates
  pco_scores <- as.data.frame(pcoa_result$vectors[, 1:n_dims])
  colnames(pco_scores) <- paste0("PCo", 1:n_dims)
  pco_scores$fish_id <- rownames(pco_scores)
  
  # merge with metadata
  model_data <- pco_scores %>%
    left_join(meta_data, by = "fish_id")
  
  # variance explained
  var_exp <- pcoa_result$values$Relative_eig[1:n_dims]
  
  # return results
  list(
    pcoa_result = pcoa_result,
    scores = pco_scores,
    model_data = model_data,
    variance_explained = var_exp,
    method_name = method_name
  )
}



plot_pcoa <- function(pcoa_obj, shape_by = NULL, label_points = FALSE,
                      meta_data = NULL, cluster_obj = NULL, 
                      show_ellipses = FALSE, ellipse_type = "norm", 
                      show_shape_legend = TRUE, show_cluster_legend = TRUE, 
                      fit_vectors = FALSE, vector_vars = NULL,
                      vector_label = NULL, show_vector_pvals = FALSE, 
                      vector_colors = NULL,  
                      vector_scale = 3, vector_label_size = 4, 
                      pc_a = 1, pc_b = 2, title = NULL, 
                      output_dir = paste0(savepath,"figs/pcoa"), 
                      save_plot = FALSE) {
  
  require(ggrepel)
  require(vegan)
  
  # use model_data from pcoa_obj
  plot_data <- pcoa_obj$model_data
  
  # add cluster info if provided
  if(!is.null(cluster_obj)){
    cluster_df <- data.frame(
      fish_id = names(cluster_obj$clusters),
      cluster = as.factor(cluster_obj$clusters)
    )
    plot_data <- plot_data %>% left_join(cluster_df, by = "fish_id")
  }
  
  # variance explained
  var_exp <- round(pcoa_obj$variance_explained * 100, 1)
  var_a <- var_exp[pc_a]
  var_b <- var_exp[pc_b]
  
  # get PC column names
  pc_a_name <- paste0("PCo", pc_a)
  pc_b_name <- paste0("PCo", pc_b)
  
  vector_data <- NULL
  if(fit_vectors && !is.null(vector_vars)) {
    # ordination scores for the specified axes
    ord_scores <- plot_data %>% 
      select(all_of(c(pc_a_name, pc_b_name))) %>% 
      as.data.frame()
    
    # environmental variables
    env_data <- plot_data %>% 
      select(all_of(vector_vars)) %>% 
      as.data.frame()
    
    # remove rows with NAs
    complete_cases <- complete.cases(ord_scores, env_data)
    ord_scores <- ord_scores[complete_cases, ]
    env_data <- env_data[complete_cases, , drop = FALSE]
    
    # check if we have data
    has_data <- !is.null(env_data) && !is.null(nrow(env_data)) && 
      !is.null(ncol(env_data)) && nrow(env_data) > 0 && ncol(env_data) > 0
    
    # fit vectors using vegan envfit
    if(has_data) {
      vector_fit <- envfit(ord_scores, env_data, permutations = 999, na.rm = TRUE)
      
      # extract vector scores
      vector_scores <- as.data.frame(scores(vector_fit, display = "vectors"))
      colnames(vector_scores) <- c(pc_a_name, pc_b_name)
      vector_scores$variable <- rownames(vector_scores)
      vector_scores$pval <- vector_fit$vectors$pvals
      vector_scores$r2 <- vector_fit$vectors$r
      
      # add colors to vector_data
      if(!is.null(vector_colors)) {
        vector_scores$color <- vector_colors[vector_scores$variable]
        # use black as default for any missing colors
        vector_scores$color[is.na(vector_scores$color)] <- "black"
      } else {
        vector_scores$color <- "black"
      }
      
      # create label_text - use custom labels if provided
      if(!is.null(vector_label) && length(vector_label) == nrow(vector_scores)) {
        vector_scores$display_label <- vector_label
      } else if(!is.null(vector_label) && is.list(vector_label)) {
        # if vector_label is a named list/vector, map it
        vector_scores$display_label <- sapply(vector_scores$variable, function(v) {
          ifelse(v %in% names(vector_label), vector_label[[v]], v)
        })
      } else {
        # use variable names as default
        vector_scores$display_label <- vector_scores$variable
      }
      
      vector_data <- vector_scores # filter(pval < .05) # filter by p-value
      
      # print summary
      if(nrow(vector_data) > 0) {
        message("\nSignificant vectors (p < ", .05, "):")
        print(vector_data %>% select(variable, r2, pval) %>% arrange(pval))
      } else {
        message("\nNo significant vectors found at p < ", .05)
        vector_data <- NULL
      }
    } else {
      message("\nInsufficient data for vector fitting after removing NAs")
    }
  }
  
  # base plot - add color here if cluster_obj exists
  if(!is.null(cluster_obj)){
    p <- ggplot(plot_data, aes(x = .data[[pc_a_name]], y = .data[[pc_b_name]], 
                               color = cluster, fill = cluster))
  } else {
    p <- ggplot(plot_data, aes(x = .data[[pc_a_name]], y = .data[[pc_b_name]]))
  }
  
  # add shape
  if (!is.null(shape_by) && shape_by %in% colnames(plot_data)) {
    p <- p + geom_point(aes(shape = .data[[shape_by]]), size = 3)+
      scale_shape_manual(values = c(21, 22), name = "Species",
                         guide = if(show_shape_legend) "legend" else "none")
  } else {
    p <- p + geom_point(size = 3, alpha = 0.7)
  }
  
  # add cluster ellipses
  if(!is.null(cluster_obj) && show_ellipses){
    p <- p + stat_ellipse(aes(color = cluster, fill = cluster), 
                          type = ellipse_type, level = 0.95, geom = "polygon", 
                          alpha = 0.2, show.legend = FALSE)
  }
  
  # color for cluster ellipses
  if (!is.null(cluster_obj)) {
    n_clusters <- length(unique(plot_data$cluster))
    cluster_colors <- c("#62A8BF", "#8A6BB0", "#6B8E9F", "#67A38A", 
                                 "#B7D3F2", "#AFAFDC", "#DCAB6B", "#F4B9B2")
                                 p <- p + 
                                   scale_color_manual(values = cluster_colors, name = "Cluster",
                                                      guide = if(show_cluster_legend) "legend" else "none") +
                                   scale_fill_manual(values = cluster_colors, name = "Cluster",
                                                     guide = if(show_cluster_legend) "legend" else "none")
  }
  
  if(!is.null(vector_data) && nrow(vector_data) > 0) {
    # add vector arrows with custom colors
    p <- p + 
      geom_segment(data = vector_data,
                   aes(x = 0, y = 0,
                       xend = .data[[pc_a_name]] * vector_scale,
                       yend = .data[[pc_b_name]] * vector_scale),
                   arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                   color = vector_data$color,
                   alpha = 1, linewidth = 1, inherit.aes = FALSE)
    
    # create label text with or without p-values (for arrows - can be disabled)
    if(show_vector_pvals) {
      vector_data$label_text <- paste0(vector_data$display_label, 
                                       "\n(p=", round(vector_data$pval, 3), ")")
    } else {
      vector_data$label_text <- vector_data$display_label
    }
    
    # add vector labels on arrows (optional - can comment out if only want corner legend)
    # p <- p + 
    #   geom_text_repel(data = vector_data,
    #                  aes(x = .data[[pc_a_name]] * vector_scale * 1.1,
    #                      y = .data[[pc_b_name]] * vector_scale * 1.1,
    #                      label = label_text),
    #                  color = "black", size = vector_label_size,
    #                  fontface = "bold", inherit.aes = FALSE,
    #                  box.padding = 1, segment.color = NA)
    
    # vector legend in top corner
    # plot limits for positioning
    x_range <- range(plot_data[[pc_a_name]], na.rm = TRUE)
    y_range <- range(plot_data[[pc_b_name]], na.rm = TRUE)
    
    # create legend data frame with significance flag
    n_vectors <- nrow(vector_data)
    
    # calculate consistent spacing
    y_increment <- diff(y_range) * 0.07  # 7% of y-axis range per item - consistent spacing
    y_start <- y_range[1] + diff(y_range) * 0.02  # start 2% from bottom
    
    legend_data <- data.frame(
      x = rep(x_range[1] + diff(x_range) * -0.1, n_vectors),  # left side (.1% from left)
      y = y_start + (0:(n_vectors-1)) * y_increment,  
      label = paste0(vector_data$display_label, ": p = ", 
                     format(round(vector_data$pval, 3), nsmall = 2)),
      color = vector_data$color,
      is_significant = vector_data$pval < 0.05
    )
    
    # split into significant and non-significant
    legend_sig <- legend_data[legend_data$is_significant, ]
    legend_nonsig <- legend_data[!legend_data$is_significant, ]
    
    # add significant vectors (bold)
    if(nrow(legend_sig) > 0) {
      p <- p + 
        geom_text(data = legend_sig,
                  aes(x = x, y = y, label = label),
                  color = legend_sig$color,
                  hjust = 0,  # left-align text
                  vjust = 0,  # bottom-align text
                  size = 3.5,
                  fontface = "bold",
                  inherit.aes = FALSE)
    }
    
    # add non-significant vectors 
    if(nrow(legend_nonsig) > 0) {
      p <- p + 
        geom_text(data = legend_nonsig,
                  aes(x = x, y = y, label = label),
                  color = legend_nonsig$color,
                  hjust = 0,  # left-align text
                  vjust = 0,  # bottom-align text
                  size = 3.5,
                  fontface = "plain",
                  inherit.aes = FALSE)
    }
  }
  
  # add fish labels (only if label_points = TRUE)
  if (label_points) {
    p <- p + geom_text_repel(aes(label = substring(fish_id, 6)),
                             min.segment.length = Inf, size = 3)
  }
  
  p <- p + 
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(title = title %||% paste("PCoA -", pcoa_obj$method_name),
         x = paste0("PCo", pc_a, " (", var_a, "%)"),
         y = paste0("PCo", pc_b, " (", var_b, "%)")) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = c(.3, .99),
      legend.justification = c(1, 1),
      legend.box.just = "right",
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
      legend.key.size = unit(0.4, "cm"),  
      legend.text = element_text(size = 8), 
      legend.title = element_text(size = 9), 
      legend.spacing.y = unit(0.1, "cm"), 
      legend.box.spacing = unit(0.2, "cm") 
    )
  
  # save plot
  if(save_plot) {
    filename <- paste0(pcoa_obj$method_name, "_pcoa_", pc_a, "v", pc_b, ".png")
    filepath <- file.path(output_dir, filename)
    ggsave(filepath, plot = p, width = 8, height = 7)
    message("Saved: ", filepath)
  }
  
  return(p)
}




# Functions for PCA
run_jerk_pca <- function(jerk_matrix,
                         n_dims = 3, save = TRUE) {
  
  # check number of dims are available
  max_dims <- min(ncol(jerk_matrix), nrow(jerk_matrix))
  if(n_dims > max_dims) {
    warning(sprintf("Requested %d dimensions but only %d available. Using %d.", 
                    n_dims, max_dims, max_dims))
    n_dims <- max_dims
  }
  
  
  # PCA
  cat("Running PCA...\n")
  pca_result <- prcomp(jerk_matrix, center = FALSE, scale. = FALSE)
  
  # PC scores
  pc_scores <- as.data.frame(pca_result$x[, 1:n_dims, drop = FALSE])
  colnames(pc_scores) <- paste0("PC", 1:n_dims)
  pc_scores$fish_id <- rownames(pc_scores)
  
  # PC matrix
  pc_matrix <- pc_scores %>% 
    select(starts_with("PC")) %>% 
    as.matrix()
  rownames(pc_matrix) <- pc_scores$fish_id
  
  # variance explained
  var_explained <- summary(pca_result)$importance[2, 1:n_dims]
  cat("\nVariance explained by PC axes:\n")
  for(i in 1:n_dims) {
    cat(sprintf("PC%d: %.2f%%\n", i, var_explained[i] * 100))
  }
  
  # loadings 
  loadings <- as.data.frame(pca_result$rotation[, 1:n_dims, drop = FALSE])
  loadings$variable <- rownames(loadings)
  
  
  output <- list(
    pca_result = pca_result,
    pc_scores = pc_scores,
    pc_matrix = pc_matrix,
    variance_explained = var_explained,
    loadings = loadings
  )
  
  # save 
  if(save) {
    if(!dir.exists(paste0(savepath,"jerk_results/"))) dir.create(paste0(savepath,"jerk_results/"), recursive = TRUE)
    saveRDS(output, file.path(paste0(savepath,"jerk_results/jerk_pca_results.rds")))
    message("Saved: ", file.path(paste0(savepath,"jerk_results/jerk_pca_results.rds")))
    
  }  
  return(output)
  
}



# Function for jerk summary
summarize_jerk <- function(jerk_data, burst_threshold = 1.25,  #assumes act_v_sum already present
                           return_table = TRUE, return_z_scores = TRUE,
                           save = TRUE) {
  
  # main summary
  jerk_summary <- jerk_data %>%
    group_by(fish_id) %>%
    summarize(
      
      # timing
      peak_timing_proportion = which.max(act_v_sum) / n(),  # 0-1 scale; #n being total numb obs for that fish
      
      # fatigue 
      fatigue = {
        peak_idx <- which.max(act_v_sum)
        if(peak_idx < n() - 2) { # need at least 3 points after peak
          # use slope of regression line from peak to end
          post_peak_time <- seq(peak_idx, n())
          post_peak_vals <- act_v_sum[peak_idx:n()]
          coef(lm(post_peak_vals ~ post_peak_time))[2] * -1 # returns slope
        } else NA
      },
      
      # burst (using provided threshold)
      n_bursts = {
        rle_obj <- rle(act_v_sum > burst_threshold) #run length encoding (groups burst or no burst in order)
        sum(rle_obj$values) #count how many bursts there were 
      },
      
      fight_duration = max(duration),
      total_effort = sum(act_v_sum, na.rm = TRUE),
      mean_actvsum = mean(act_v_sum, na.rm = TRUE),
      sd_jerk = sd(act_v_sum, na.rm = TRUE),
      time_at_max = sum(act_v_sum == max(act_v_sum), na.rm = TRUE),
      early_intensity = mean(act_v_sum[1:floor(n()/3)], na.rm = TRUE),
      late_intensity = mean(act_v_sum[ceiling(2*n()/3):n()], na.rm = TRUE),
      cv_jerk = sd_jerk / abs(mean_actvsum)
    )
  
  # output list
  output <- list(summary = jerk_summary)
  
  # table
  if(return_table) {
    output$table <- jerk_summary %>%
      select(fish_id, n_bursts, mean_actvsum, cv_jerk, 
             early_intensity, late_intensity, fatigue, total_effort) %>%
      mutate(across(where(is.numeric), ~round(., 2)))
  }
  
  # z-standardize 
  if(return_z_scores) {
    output$z_scores <- jerk_summary %>%
      mutate(across(c(n_bursts, mean_actvsum, cv_jerk, 
                      early_intensity, late_intensity, fatigue, 
                      total_effort),
                    ~as.numeric(scale(.)), #converts to numeric
                    .names = "{.col}_z")) %>%
      select(fish_id, ends_with("_z"))
  

    # make matrix
    jerk_matrix <- output$z_scores %>% 
      select(ends_with("_z")) %>%
      as.matrix()
    rownames(jerk_matrix) <- output$z_scores$fish_id
    
    output$z_matrix <- jerk_matrix
  }
  
  # save 
  if(save) {
    if(!dir.exists(paste0(savepath,"jerk_results/"))) dir.create(paste0(savepath,"jerk_results/"), recursive = TRUE)
    saveRDS(output, file.path(paste0(savepath,"jerk_results/jerk_results.rds")))
    message("Saved: ", file.path(paste0(savepath,"jerk_results/jerk_results.rds")))
    
  }  
  return(output)
}

# Function for accel summary
summarize_accel <- function(accel_data, 
                            x_col = "X",
                            y_col = "Y",
                            z_col = "Z",
                            burst_threshold = 1.25, 
                            return_table = TRUE, return_z_scores = TRUE,
                            save = TRUE) {
  
  
  # main summary
  accel_summary <- accel_data %>%
    group_by(fish_id) %>%
    summarize(
      
      # timing
      peak_timing_proportion = which.max(sv) / n(),  # 0-1 scale; #n being total numb obs for that fish
      
      # fatigue 
      fatigue = {
        peak_idx <- which.max(sv)
        if(peak_idx < n() - 2) { # need at least 3 points after peak
          # use slope of regression line from peak to end
          post_peak_time <- seq(peak_idx, n())
          post_peak_vals <- sv[peak_idx:n()]
          coef(lm(post_peak_vals ~ post_peak_time))[2] * -1 # returns slope
        } else NA
      },
      
      # burst (using provided threshold)
      n_bursts = {
        rle_obj <- rle(sv > burst_threshold) #run length encoding (groups burst or no burst in order)
        sum(rle_obj$values) #count how many bursts there were 
      },
      
      fight_duration = max(fight_duration_s),
      total_effort = sum(sv, na.rm = TRUE),
      mean_sv = mean(sv, na.rm = TRUE),
      sd_sv = sd(sv, na.rm = TRUE),
      time_at_max = sum(sv == max(sv), na.rm = TRUE),
      early_intensity = mean(sv[1:floor(n()/3)], na.rm = TRUE),
      late_intensity = mean(sv[ceiling(2*n()/3):n()], na.rm = TRUE),
      cv_sv = sd_sv / abs(mean_sv)
    )
  
  # output list
  output <- list(summary = accel_summary)
  
  # table
  if(return_table) {
    output$table <- accel_summary %>%
      select(fish_id, n_bursts, mean_sv, cv_sv, 
             early_intensity, late_intensity, fatigue, total_effort) %>%
      mutate(across(where(is.numeric), ~round(., 2)))
  }
  
  # z-standardize 
  if(return_z_scores) {
    output$z_scores <- accel_summary %>%
      mutate(across(c(n_bursts, mean_sv, cv_sv, 
                      early_intensity, late_intensity, fatigue, 
                      total_effort),
                    ~as.numeric(scale(.)), #converts to numeric
                    .names = "{.col}_z")) %>%
      select(fish_id, ends_with("_z"))
    
    
    # make matrix
    accel_matrix <- output$z_scores %>% 
      select(ends_with("_z")) %>%
      as.matrix()
    rownames(accel_matrix) <- output$z_scores$fish_id
    
    output$z_matrix <- accel_matrix
  }
  
  # save 
  if(save) {
    if(!dir.exists(paste0(savepath,"accel_results/"))) dir.create(paste0(savepath,"accel_results/"), recursive = TRUE)
    saveRDS(output, file.path(paste0(savepath,"accel_results/accel_results.rds")))
    message("Saved: ", file.path(paste0(savepath,"accel_results/accel_results.rds")))
    
  }  
  return(output)
}

# Function for complete brms analysis
analyze_fish_fight <- function(jerk_data,
                               dtw_matrices = NULL, 
                               meta_data,
                               response_vars = NULL,
                               covariates = NULL,
                               n_dims = 3,
                               iter = 4000,
                               warmup = 3000,
                               cores = 4,
                               save_model = TRUE,
                               output_dir = paste0(savepath,"models/brms_outputs")) {
  
  require(ape)
  require(brms)
  require(bayesplot)
  require(tidyverse)
  
  cat("\n========================================\n")
  cat("║   FISH FIGHT brms ANALYSIS SUITE   ║\n")
  cat("\n========================================\n")
  
  # input validation
  cat("Validating inputs...\n")
  if(!("fish_id" %in% colnames(jerk_data))) stop("jerk_data must have 'fish_id' column")
  if(!("fish_id" %in% colnames(meta_data))) stop("meta_data must have 'fish_id' column")
  
  # DTW validation only if provided
  if(!is.null(dtw_matrices)) {
    if(!is.list(dtw_matrices)) stop("dtw_matrices must be a named list")
  }
  
  if(any(!response_vars %in% colnames(meta_data))) {
    stop("Response variables not found in meta_data: ",
         paste(response_vars[!response_vars %in% colnames(meta_data)], collapse=", "))
  }
  
  cat("Inputs validated\n")
  cat("  - Jerk data:", nrow(jerk_data), "fish\n")
  cat("  - Meta data:", nrow(meta_data), "fish\n")
  
  if(!is.null(dtw_matrices)) {
    cat("  - DTW methods:", length(dtw_matrices), "\n")
  } else {
    cat("  - DTW methods: NONE (summary metrics only)\n")
  }
  
  cat("  - Response vars:", paste(response_vars, collapse=", "), "\n\n")
  
  # initialize storage
  all_results <- list()
  
  # ============================================================
  # NESTED FUNCTION: Fit single model
  # ============================================================
  fit_single_model <- function(input_data, data_type, method_name) {
    
    cat("\n ANALYZING:", method_name, "(", toupper(data_type), ")\n")
    
    # --- PCoA ORDINATION (for DTW) ---
    if(data_type == "dtw") {
      dtw_dist <- as.dist(input_data)
      ord_result <- pcoa(dtw_dist)
      dim_scores <- as.data.frame(ord_result$vectors[, 1:n_dims])
      colnames(dim_scores) <- paste0("Dim", 1:n_dims)
      dim_scores$fish_id <- rownames(dim_scores)
      var_explained <- ord_result$values$Relative_eig[1:n_dims]
      analysis_type <- "PCoA"
      
      # --- PCA ORDINATION (for jerk/accel summary metrics) ---
    } else {  # jerk
      predictor_vars <- colnames(input_data)[grepl("_z$", colnames(input_data))]
      jerk_matrix <- input_data %>% 
        select(all_of(predictor_vars)) %>% 
        as.matrix()
      rownames(jerk_matrix) <- input_data$fish_id
      
      ord_result <- prcomp(jerk_matrix, center=FALSE, scale.=FALSE)
      dim_scores <- as.data.frame(ord_result$x[, 1:n_dims])
      colnames(dim_scores) <- paste0("Dim", 1:n_dims)
      dim_scores$fish_id <- rownames(dim_scores)
      var_explained <- summary(ord_result)$importance[2, 1:n_dims]
      analysis_type <- "PCA"
    }
    
    cat("  Variance explained:", paste0(round(var_explained*100, 1), "%", collapse=", "), "\n")
    
    # --- MERGE WITH METADATA ---
    model_data <- dim_scores %>% left_join(meta_data, by = "fish_id")
    
    # prep factors
    if("dna_sex" %in% covariates) {
      model_data$dna_sex <- factor(model_data$dna_sex, levels = c("male", "presumed_female"))
    }
    if("sp" %in% covariates) {
      model_data$sp <- factor(model_data$sp, levels = c("CHINOOK", "COHO"))
    }
    
    # standardize numeric covariates
    cov_final <- covariates
    for(cov in covariates) {
      if(is.numeric(model_data[[cov]]) && !grepl("_z$", cov)) {
        z_name <- paste0(cov, "_z")
        model_data[[z_name]] <- as.vector(scale(model_data[[cov]]))
        cov_final[cov_final == cov] <- z_name
      }
    }
    
    # --- BUILD FORMULA ---
    dim_terms <- paste0("Dim", 1:n_dims, collapse = " + ")
    cov_terms <- paste(cov_final, collapse = " + ")
    
    if(data_type == "jerk") {
      # fight_time_s is in model_data from the meta_data join
      if("fight_time_s" %in% colnames(model_data)) {
        model_data$fight_time_s_z <- as.vector(scale(model_data$fight_time_s))
        cov_terms <- paste(cov_terms, "fight_time_s_z", sep = " + ")
      }
    }
    
    response_terms <- paste0("mvbind(", paste(response_vars, collapse = ", "), ")")
    formula_str <- paste0(response_terms, " ~ ", dim_terms, " + ", cov_terms)
    
    cat("  Formula:", formula_str, "\n")
    
    bf_formula <- bf(as.formula(formula_str)) + set_rescor(TRUE)
    
    # --- FIT MODEL ---
    cat("  Fitting Bayesian model...\n")
    
    priors <- c(prior(lkj(2), class = "rescor"))
    
    fit <- brm(bf_formula,
               data = model_data,
               prior = priors,
               control = list(adapt_delta = 0.95, max_treedepth = 15),
               iter = iter,
               warmup = warmup,
               cores = cores,
               silent = 2,
               refresh = 0,
               backend = "cmdstanr")
    
    cat("Model fitted.\n")
    
    # --- EXTRACT RESULTS ---
    all_params <- posterior_summary(fit)
    
    # diagnostics
    neff <- neff_ratio(fit)
    max_rhat <- max(summary(fit)$fixed[, "Rhat"], na.rm = TRUE)
    
    # save 
    if(save_model) {
      if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      saveRDS(fit, file.path(output_dir, paste0(method_name, "_model.rds")))
    }
    
    # return results
    return(list(
      method_name = method_name,
      model_data = model_data,
      data_type = data_type,
      analysis_type = analysis_type,
      model = fit,
      variance_explained = var_explained,
      all_params = all_params,
      max_rhat = max_rhat,
      min_neff = min(neff)
    ))
  }
  
  # ============================================================
  # MAIN ANALYSIS LOOP
  # ============================================================
  
  # ALWAYS run jerk/accel summary analysis
  cat("\n>>> Running SUMMARY METRICS analysis <<<\n")
  all_results[["summary_metrics"]] <- fit_single_model(
    input_data = jerk_data,
    data_type = "jerk",
    method_name = "summary_metrics"
  )
  
  # ONLY run DTW analyses if matrices provided
  if(!is.null(dtw_matrices)) {
    cat("\n>>> Running DTW analyses <<<\n")
    for(method_name in names(dtw_matrices)) {
      all_results[[method_name]] <- fit_single_model(
        input_data = dtw_matrices[[method_name]],
        data_type = "dtw",
        method_name = method_name
      )
    }
  }
  
  
  # ============================================================
  # CREATE COMPARISON TABLE
  # ============================================================
  cat("\n Creating comparison table...\n")
  
  comparison <- data.frame()
  
  for(method_name in names(all_results)) {
    result <- all_results[[method_name]]
    
    # basics
    data_type <- result$data_type
    analysis_type <- result$analysis_type
    var_explained <- sum(result$variance_explained) * 100
    
    # model fit
    ic_score <- tryCatch({
      loo_result <- loo(result$model)
      loo_result$estimates["looic", "Estimate"]
    }, error = function(e) {
      cat("    Warning: LOO failed:", conditionMessage(e), "\n")
      NA
    })
    
    bayes_r2 <- tryCatch({
      median(bayes_R2(result$model))
    }, error = function(e) NA)
    
    # EFFECTS
    all_params <- result$all_params
    
    # effects: dimensions
    dim_params <- all_params[grepl("Dim[1-3]", rownames(all_params)), ]
    n_sig_dim <- sum(dim_params[, "Q2.5"] * dim_params[, "Q97.5"] > 0)
    sig_dims <- rownames(dim_params)[dim_params[, "Q2.5"] * dim_params[, "Q97.5"] > 0]
    sig_dims <- unique(gsub("^b_[^_]+_", "", sig_dims))
    sig_dims_str <- if(length(sig_dims) > 0) paste(sig_dims, collapse = ", ") else "none"
    
    # effects: covariates
    covariate_pat <- paste0(covariates, collapse = "|")
    covariate_params <- all_params[grepl(covariate_pat, rownames(all_params)), ]
    n_sig_covs <- sum(covariate_params[, "Q2.5"] * covariate_params[, "Q97.5"] > 0)
    sig_covs <- rownames(covariate_params)[covariate_params[, "Q2.5"] * covariate_params[, "Q97.5"] > 0]
    sig_covs <- unique(gsub("^b_[^_]+_", "", sig_covs))
    sig_covs_str <- if(length(sig_covs) > 0) paste(sig_covs, collapse = ", ") else "none"
    
    # add to table
    comparison <- rbind(comparison, data.frame(
      method = method_name,
      data_type = data_type,
      analysis_type = analysis_type,
      var_explained_pct = round(var_explained, 1),
      ic_score = round(ic_score, 1),
      bayes_r2 = round(bayes_r2, 3),
      n_sig_dimensions = n_sig_dim,
      sig_dimensions = sig_dims_str,
      n_sig_covariates = n_sig_covs,
      sig_covariates = sig_covs_str,
      max_rhat = round(result$max_rhat, 3),
      min_neff = round(result$min_neff, 3),
      stringsAsFactors = FALSE
    ))
  }
  
  # sort by performance
  comparison <- comparison %>% arrange(ic_score)
  
  # ============================================================
  # PRINT SUMMARY
  # ============================================================
  cat("\n========================================\n")
  cat("║          ANALYSIS COMPLETE!            ║\n")
  cat("\n========================================\n")
  
  cat("BEST MODEL:", comparison$method[1], "\n")
  cat("  - Data type:", comparison$data_type[1], "\n")
  cat("  -IC:", comparison$ic_score[1], "\n")
  cat("  - Bayes R²:", comparison$bayes_r2[1], "\n")
  cat("  - Significant dims:", comparison$n_sig_dimensions[1], "\n\n")
  
  cat("TOP 3 MODELS:\n")
  print(comparison[1:min(3, nrow(comparison)), c("method", "ic_score", "bayes_r2")], 
        row.names = FALSE)
  
  cat("\n Full comparison table available in: results$comparison\n")
  cat("All models available in: results$models\n")
  if(save_model) {
    cat("Models saved to:", output_dir, "\n")
  }
  
  # ============================================================
  # RETURN RESULTS
  # ============================================================
  return(list(
    models = all_results,
    comparison = comparison,
    best_model = all_results[[comparison$method[1]]],
    summary = list(
      n_methods = length(all_results),
      best_method = comparison$method[1],
      best_ic_score = comparison$ic_score[1]
    )
  ))
  
}

  

