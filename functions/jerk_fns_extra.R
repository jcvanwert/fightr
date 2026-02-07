################################################################################
# DTW FUNCTIONS - EXTRA
# Author: Jacey Van Wert 
# Date: Jan 02 2026
# This file has a variety of functions useful for the accel data analysis but NOT necessary
################################################################################



# Set path accordingly - this is the overview folder; everything will save and create folders within
path <- "~/UFL Dropbox/Jacey Van Wert/JVW_Collab/Salmon_accel/jerk_2024"




# Function creating a PCA on the physiological data
run_physio_pca <- function(meta_data,
                           fish_id_col = "fish_id",
                           physio_vars = c("lactate", "glucose"),
                           n_dims = 3) {
  
  # extract and prepare data
  physio_matrix <- meta_data %>%
    select(all_of(physio_vars)) %>%
    as.matrix()
  
  rownames(physio_matrix) <- meta_data[[fish_id_col]]
  
  # run PCA (center=TRUE, scale=TRUE standardizes the variables)
  cat("Running PCA...\n")
  pca_result <- prcomp(physio_matrix, center = TRUE, scale. = TRUE)
  
  # PC scores
  pc_scores <- as.data.frame(pca_result$x[, 1:n_dims])
  colnames(pc_scores) <- paste0("PC", 1:n_dims)
  pc_scores$fish_id <- rownames(pc_scores)
  
  # variance explained
  var_explained <- summary(pca_result)$importance[2, 1:n_dims]
  cat("\nVariance explained by PC axes:\n")
  for(i in 1:n_dims) {
    cat(paste0("PC", i, ": "), round(var_explained[i] * 100, 2), "%\n")
  }
  
  # loadings
  loadings <- as.data.frame(pca_result$rotation[, 1:n_dims])
  loadings$variable <- rownames(loadings)
  
  cat("\nLoadings (variable contributions):\n")
  print(round(loadings[, 1:n_dims], 3))
  
  return(list(
    pca_result = pca_result,
    pc_scores = pc_scores,
    variance_explained = var_explained,
    loadings = loadings
  ))
}



# Function for plotting PCA 
plot_pca <- function(pca_res, meta_data = NULL, shape_by = NULL, 
                     label_points = TRUE, show_vectors=TRUE,
                     cluster_obj = NULL,
                     show_ellipses = FALSE,
                     ellipse_type = "norm",
                     show_shape_legend = TRUE,
                     show_cluster_legend=TRUE,
                     pc_a = 1, pc_b = 2,
                     title = NULL) {
  
  require(ggrepel)
  
  # merge with metadata
  if (!is.null(meta_data)) {
    plot_data <- pca_res$pc_scores %>%
      left_join(meta_data, by = "fish_id")
  } else {
    plot_data <- pca_res$pc_scores
  }
  
  # add cluster info if provided
  if(!is.null(cluster_obj)){
    cluster_df <- data.frame(
      fish_id = names(cluster_obj$clusters),
      cluster = as.factor(cluster_obj$clusters)
    )
    plot_data <- plot_data %>% 
      left_join(cluster_df, by = "fish_id")
    
  }
  
  # var 
  var_exp <- round(pca_res$variance_explained * 100, 1)
  var_a <- var_exp[pc_a]
  var_b <- var_exp[pc_b]
  
  
  # get PC column names
  pc_a_name <- paste0("PC", pc_a)
  pc_b_name <- paste0("PC", pc_b)
  
  
  # base plot- add color here if cluster_obj exists
  
  if(!is.null(cluster_obj)){
    p <- ggplot(plot_data, aes(x = .data[[pc_a_name]], y = .data[[pc_b_name]],
                               color = cluster, fill = cluster))
  } else {
    p <- ggplot(plot_data, aes(x = .data[[pc_a_name]], y = .data[[pc_b_name]]))
    
  }
  
  
  
  # add shape 
  if (!is.null(shape_by) && shape_by %in% colnames(plot_data)) {
    p <- p + geom_point(aes(shape = .data[[shape_by]]), size = 4)+
      scale_shape_manual(values=c(16,17), name = "Species",
                         guide = if(show_shape_legend) "legend" else "none")
    # scale_fill_manual(values = c("#3670A0", "#4FA4E5"), name = "Species")
  } else {
    p <- p + geom_point(size = 4, alpha = 0.7)
  }
  
  
  # add cluster ellipses 
  if(!is.null(cluster_obj) && show_ellipses){
    p <- p + stat_ellipse(aes(color = cluster, fill = cluster),
                          type=ellipse_type,
                          level=0.95,
                          geom= "polygon",
                          alpha = 0.2,
                          show.legend = FALSE)
    
  }
  
  # color for cluster ellipses
  if (!is.null(cluster_obj)) {
    n_clusters <- length(unique(plot_data$cluster))
    cluster_colors <- c("#62A8BF", "#8A6BB0", "#6B8E9F", "#67A38A", "#B7D3F2", "#AFAFDC", "#DCAB6B", "#F4B9B2")
    
    p <- p + 
      scale_color_manual(values = cluster_colors, name = "Cluster",
                         guide = if(show_cluster_legend) "legend" else "none") +
      scale_fill_manual(values = cluster_colors, name = "Cluster",
                        guide = if(show_cluster_legend) "legend" else "none")
  }
  
  
  # add labels
  if (label_points) {
    p <- p + geom_text_repel(aes(label = substring(fish_id,6)),
                             min.segment.length = Inf, size = 3)
  }
  
  
  # add biplot vectors
  if(show_vectors){
    scale_factor <- 2.24
    
    arrow_data <- pca_res$loadings %>% 
      mutate(
        PC_a_arrow = .data[[pc_a_name]] * scale_factor,
        PC_b_arrow = .data[[pc_b_name]] * scale_factor,
        variable_label = gsub("_z$", "", variable) %>% 
          gsub("_", " ", .),
        variable_num = row_number(), 
        arrow_length = sqrt(PC_a_arrow^2 + PC_b_arrow^2),
        shorten_by = 0.15,
        PC_a_arrow_end = PC_a_arrow * (1 - shorten_by/arrow_length),
        PC_b_arrow_end = PC_b_arrow * (1 - shorten_by/arrow_length)
      )
    
    
    # arrows as diamond
    p <- p + 
      geom_segment(data = arrow_data,
                   aes(x = 0, y = 0, xend = PC_a_arrow, yend = PC_b_arrow),
                   color = "black", linewidth = 0.6, inherit.aes = FALSE)+
      geom_point(data = arrow_data,
                 aes(x = PC_a_arrow, y = PC_b_arrow),
                 shape = 23,  # diamond shape
                 size = 6,
                 fill = "black",
                 color = "black",
                 inherit.aes = FALSE) +
      geom_text(data = arrow_data,
                aes(x = PC_a_arrow, y = PC_b_arrow, label = variable_num),
                color = "white",
                size = 3,
                fontface = "bold",
                inherit.aes = FALSE)
    
  }
  
  
  p <- p +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(title = title,
         x = paste0("PC", pc_a, " (", var_a, "%)"),
         y = paste0("PC", pc_b, " (", var_b, "%)"))+
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = c(0.98, 0.98),   
      legend.justification = c(1, 1),
      legend.box.just = "right",
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
    )
  
  return(p)
}

# for pca plot
create_variable_key <- function(pca_res, pc_a = 1, pc_b = 2, col_spacing = 0.15,
                                custom_labels = NULL) {
  
  # get PC column names
  pc_a_name <- paste0("PC", pc_a)
  pc_b_name <- paste0("PC", pc_b)
  
  
  arrow_data <- pca_res$loadings %>% 
    mutate(
      PC_a_arrow = .data[[pc_a_name]],
      PC_b_arrow = .data[[pc_b_name]],
      variable_label = gsub("_z$", "", variable) %>% 
        gsub("_", " ", .),
      variable_num = row_number()
    )
  
  
  #if custom labels are provided
  if(!is.null(custom_labels)) {
    cat("Variables in PCA:\n")
    print(arrow_data$variable)
    cat("\nCustom label keys:\n")
    print(names(custom_labels))
    
    for(i in 1:nrow(arrow_data)) {
      orig_var <- gsub("_z$", "", arrow_data$variable[i])
      cat("\nTrying to match:", orig_var)
      if(orig_var %in% names(custom_labels)) {
        cat(" -> FOUND:", custom_labels[[orig_var]])
        arrow_data$variable_label[i] <- custom_labels[[orig_var]]
      } else {
        cat(" -> NOT FOUND")
      }
    }
  }
  
  # create variable key
  variable_key <- data.frame(
    num = arrow_data$variable_num,
    label = arrow_data$variable_label
  )
  
  # numb rows for 2 col
  n_vars <- nrow(variable_key)
  n_rows = ceiling(n_vars/2)
  
  # split into two columns
  variable_key$col <- rep(1:2, length.out = n_vars)
  variable_key$row <- rep(1:n_rows, each = 2)[1:n_vars]
  
  # format labels
  variable_key$full_label <- paste0(variable_key$num, ": ", variable_key$label)
  
  # calculate positions 
  variable_key$x <- ifelse(variable_key$col == 1, 0, col_spacing)
  variable_key$y <- 1 - (variable_key$row) / (n_rows + 1)
  
  # box dimensions
  box_width <- col_spacing * 1.5 + 0.05
  
  # text grob for the key
  key_plot <- ggplot(variable_key, aes(x = x, y = y, label = full_label)) + 
    annotate("rect", 
             xmin = -0.02, xmax = box_width, 
             ymin = -0.05, ymax = 1.05,
             fill = "white", color = "black", linewidth = 0.4) +
    geom_text(hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(-0.1, box_width + 0.1), 
                    ylim = c(-0.1, 1.1), 
                    clip = "off") +
    theme_void() +
    theme(plot.margin = margin(1, 2, 2, 2))
  
  return(key_plot)
}

