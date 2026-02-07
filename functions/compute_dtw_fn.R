################################################################################
# COMPREHENSIVE DTW ANALYSIS FUNCTION
# Author: Jacey Van Wert 
# Date: Dec 12 2025
# This function computes multiple DTW distance matrices and creates visualizations
################################################################################

#' Compute Dynamic Time Warping Distance Matrices
#'
#' @param data Dataframe containing accelerometer/jerk data
#' @param fish_id_col Name of column containing individual fish/animal IDs
#' @param x_col Name of column containing X-axis acceleration
#' @param y_col Name of column containing Y-axis acceleration  
#' @param z_col Name of column containing Z-axis acceleration
#' @param duration_col Name of column containing time/duration values
#' @param event_id_col Name of column containing event IDs (for high-pass filter)
#' @param timestamp_col Name of column containing timestamps (for high-pass filter)
#' @param methods Vector of methods to compute. Options: "raw", "differential", "multivariate", "highpass"
#' @param normalize Logical. Use normalized DTW distance? (accounts for different sequence lengths)
#' @param highpass_smooth Integer. Smoothing window for high-pass filter (number of observations)
#' @param n_cores Integer. Number of cores for parallel processing (high-pass filter only)
#' @param create_plots Logical. Create visualization plots?
#' @param fig_output_dir Directory to save plots (if create_plots = TRUE)
#' @param plot_width Width of plots in inches
#' @param plot_height Height of plots in inches
#'
#' @return List containing:
#'   - dtw_matrices: List of distance matrices for each method
#'   - plots: List of ggplot objects (if create_plots = TRUE) in fig_output_dir
#'   - data_processed: Processed data including high-pass filtered values (if applicable)
#'
#' @export


################################################################################
# COMPREHENSIVE DTW FUNCTION 
################################################################################

compute_dtw <- function(
    data,
    fish_id_col = "fish_id",
    x_col = "act_x",
    y_col = "act_y",
    z_col = "act_z",
    duration_col = "duration",
    event_id_col = "event_id",
    timestamp_col = "date_time_utc",
    methods = c("raw", "differential", "multivariate", "highpass"),
    normalize = TRUE,
    highpass_smooth = 4,
    highpass_parallel = FALSE,  
    n_cores = NULL,
    create_plots = TRUE,
    fig_output_dir = NULL,
    plot_width = 8,
    plot_height = 7,
    axis_colors = c(x = "#21DEE8", y = "#109FBC", z = "#13505B")
) {
  
  # ============================================================================
  # SETUP & VALIDATION
  # ============================================================================
  
  cat("\n")
  cat("================================================================================\n")
  cat("                                  DTW ANALYSIS                                  \n")
  cat("================================================================================\n\n")
  
  # check required packages
  required_packages <- c("tidyverse", "dtw", "xts", "zoo")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if(length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
  }
  
  # load libraries
  suppressPackageStartupMessages({
    library(tidyverse)
    library(dtw)
    library(xts)
    library(zoo)
  })
  
  # validate inputs
  required_cols <- c(fish_id_col, x_col, y_col, z_col, duration_col)
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # create output directory
  if(create_plots && !is.null(fig_output_dir)) {
    if(!dir.exists(fig_output_dir)) {
      dir.create(fig_output_dir, recursive = TRUE)
      cat("✓ Created output directory:", fig_output_dir, "\n")
    }
  }
  
  cat("Configuration:\n")
  cat("  Individuals:", length(unique(data[[fish_id_col]])), "\n")
  cat("  Observations:", nrow(data), "\n")
  cat("  Methods:", paste(methods, collapse = ", "), "\n")
  cat("  Normalize DTW:", normalize, "\n")
  if("highpass" %in% methods) {
    cat("  High-pass smooth:", highpass_smooth, "\n")
    cat("  Parallel processing:", highpass_parallel, "\n")
  }
  cat("\n")
  
  # ============================================================================
  # PREPARE DATA STRUCTURES
  # ============================================================================
  
  cat("=== PREPARING DATA ===\n")
  
  # extract fish IDs
  fish_ids <- unique(data[[fish_id_col]])
  n_fish <- length(fish_ids)
  
  # create data lists for each axis
  jerk_list_x <- data %>%
    split(.[[fish_id_col]]) %>%
    lapply(function(x) x[[x_col]])
  
  jerk_list_y <- data %>%
    split(.[[fish_id_col]]) %>%
    lapply(function(x) x[[y_col]])
  
  jerk_list_z <- data %>%
    split(.[[fish_id_col]]) %>%
    lapply(function(x) x[[z_col]])
  
  cat("✓ Created data lists for", n_fish, "individuals\n")
  
  # create template distance matrix
  dtw_matrix_template <- matrix(0, nrow = n_fish, ncol = n_fish)
  rownames(dtw_matrix_template) <- colnames(dtw_matrix_template) <- fish_ids
  
  # initialize storage
  dtw_matrices <- list()
  plots <- list()
  data_processed <- list()
  
  # store raw data
  data_processed$raw <- data
  
  # ============================================================================
  # HELPER FUNCTION: Compute DTW for a single axis
  # ============================================================================
  
  compute_dtw_axis <- function(jerk_list_axis, axis_name, method_suffix = "") {
    dtw_matrix <- dtw_matrix_template
    
    pb <- txtProgressBar(min = 0, max = n_fish * (n_fish - 1) / 2, style = 3)
    counter <- 0
    
    for(i in 1:n_fish) {
      for(j in i:n_fish) {
        if(i != j) {
          alignment <- dtw(jerk_list_axis[[i]], jerk_list_axis[[j]])
          
          dist_val <- if(normalize) {
            alignment$normalizedDistance
          } else {
            alignment$distance
          }
          
          dtw_matrix[i, j] <- dist_val
          dtw_matrix[j, i] <- dist_val
          
          counter <- counter + 1
          setTxtProgressBar(pb, counter)
        }
      }
    }
    close(pb)
    
    return(dtw_matrix)
  }
  
  # ============================================================================
  # METHOD 1: RAW AXES DTW
  # ============================================================================
  
  if("raw" %in% methods) {
    cat("\n=== METHOD 1: RAW AXES DTW ===\n")
    
    # X axis
    cat("Computing X-axis DTW...\n")
    dtw_matrices$x_axis <- compute_dtw_axis(jerk_list_x, "X")
    
    # Y axis
    cat("Computing Y-axis DTW...\n")
    dtw_matrices$y_axis <- compute_dtw_axis(jerk_list_y, "Y")
    
    # Z axis
    cat("Computing Z-axis DTW...\n")
    dtw_matrices$z_axis <- compute_dtw_axis(jerk_list_z, "Z")
    
    cat("Raw axes DTW complete!\n")
  }
  
  # ============================================================================
  # METHOD 2: MULTIVARIATE DTW (RAW)
  # ============================================================================
  
  if("multivariate" %in% methods) {
    cat("\n=== METHOD 2: MULTIVARIATE DTW ===\n")
    
    dtw_matrix_mv <- dtw_matrix_template
    
    pb <- txtProgressBar(min = 0, max = n_fish * (n_fish - 1) / 2, style = 3)
    counter <- 0
    
    for(i in 1:n_fish) {
      for(j in i:n_fish) {
        if(i != j) {
          # combine axes into matrix (time x dimensions)
          s1 <- cbind(jerk_list_x[[i]], jerk_list_y[[i]], jerk_list_z[[i]])
          s2 <- cbind(jerk_list_x[[j]], jerk_list_y[[j]], jerk_list_z[[j]])
          
          # compute pairwise distances
          cost_matrix <- proxy::dist(s1, s2, method = "Euclidean")
          
          # DTW on multivariate series
          alignment <- dtw(cost_matrix)
          
          dist_val <- if(normalize) {
            alignment$normalizedDistance
          } else {
            alignment$distance
          }
          
          dtw_matrix_mv[i, j] <- dist_val
          dtw_matrix_mv[j, i] <- dist_val
          
          counter <- counter + 1
          setTxtProgressBar(pb, counter)
        }
      }
    }
    close(pb)
    
    dtw_matrices$multivariate <- dtw_matrix_mv
    cat("Multivariate DTW complete!\n")
  }
  
  # ============================================================================
  # METHOD 3: DIFFERENTIAL DTW
  # ============================================================================
  
  if("differential" %in% methods) {
    cat("\n=== METHOD 3: DIFFERENTIAL DTW ===\n")
    
    # compute differences
    cat("Computing first differences...\n")
    jerk_list_x_diff <- lapply(jerk_list_x, diff)
    jerk_list_y_diff <- lapply(jerk_list_y, diff)
    jerk_list_z_diff <- lapply(jerk_list_z, diff)
    
    # X axis differential
    cat("Computing X-axis differential DTW...\n")
    dtw_matrices$x_diff <- compute_dtw_axis(jerk_list_x_diff, "X_diff")
    
    # Y axis differential
    cat("Computing Y-axis differential DTW...\n")
    dtw_matrices$y_diff <- compute_dtw_axis(jerk_list_y_diff, "Y_diff")
    
    # Z axis differential
    cat("Computing Z-axis differential DTW...\n")
    dtw_matrices$z_diff <- compute_dtw_axis(jerk_list_z_diff, "Z_diff")
    
    # multivariate differential
    if("multivariate" %in% methods) {
      cat("Computing multivariate differential DTW...\n")
      
      dtw_matrix_mvd <- dtw_matrix_template
      
      pb <- txtProgressBar(min = 0, max = n_fish * (n_fish - 1) / 2, style = 3)
      counter <- 0
      
      for(i in 1:n_fish) {
        for(j in i:n_fish) {
          if(i != j) {
            s1 <- cbind(jerk_list_x_diff[[i]], jerk_list_y_diff[[i]], jerk_list_z_diff[[i]])
            s2 <- cbind(jerk_list_x_diff[[j]], jerk_list_y_diff[[j]], jerk_list_z_diff[[j]])
            
            cost_matrix <- proxy::dist(s1, s2, method = "Euclidean")
            alignment <- dtw(cost_matrix)
            
            dist_val <- if(normalize) {
              alignment$normalizedDistance
            } else {
              alignment$distance
            }
            
            dtw_matrix_mvd[i, j] <- dist_val
            dtw_matrix_mvd[j, i] <- dist_val
            
            counter <- counter + 1
            setTxtProgressBar(pb, counter)
          }
        }
      }
      close(pb)
      
      dtw_matrices$multivariate_diff <- dtw_matrix_mvd
    }
    
    cat("Differential DTW complete!\n")
  }
  
  # ============================================================================
  # METHOD 4: HIGH-PASS FILTERED DTW 
  # ============================================================================
  
  if("highpass" %in% methods) {
    cat("\n=== METHOD 4: HIGH-PASS FILTERED DTW ===\n")
    
    # check for required columns
    if(!event_id_col %in% colnames(data)) {
      warning("event_id_col '", event_id_col, "' not found. Skipping high-pass filter.")
    } else {
      
      # calculate sampling rate
      sampling_rate <- mean(diff(data[[duration_col]]), na.rm = TRUE)
      cat("  Sampling rate:", round(sampling_rate, 3), "seconds\n")
      cat("  Smoothing window:", highpass_smooth * sampling_rate, "seconds\n")
      
      # split data by event
      fishjerk_list <- split(data, data[[event_id_col]])
      cat("  Processing", length(fishjerk_list), "events\n")
      
      # processing method
      if(highpass_parallel) {
        cat("  Using parallel processing...\n")
        library(foreach)
        library(doParallel)
        
        if(is.null(n_cores)) {
          n_cores <- max(1, parallel::detectCores() - 2)
        }
        
        core.cl <- makeCluster(n_cores)
        registerDoParallel(core.cl)
        
        fishjerk_hp <- foreach(
          i = 1:length(fishjerk_list),
          .packages = c('xts', 'zoo'),
          .errorhandling = 'pass'
        ) %dopar% {

          current_data <- fishjerk_list[[i]]
          
          # create timestamps
          if(timestamp_col %in% colnames(current_data)) {
            base_time <- as.POSIXct(current_data[[timestamp_col]][1], format = "%Y-%m-%d %H:%M")
            current_data$date_time <- base_time + current_data[[duration_col]]
          } else {
            current_data$date_time <- as.POSIXct("2020-01-01") + current_data[[duration_col]]
          }
          
          # skip if too few observations
          if(nrow(current_data) < highpass_smooth) {
            return(NULL)
          }
          
          # create xts object
          fishjerk_xts <- xts(current_data[, c(x_col, y_col, z_col)], 
                              order.by = current_data$date_time)
          
          # low-pass filter (smoothed version)
          fishjerk_lp <- rollmean(fishjerk_xts, k = highpass_smooth, 
                                  align = 'center', fill = NA)
          
          # interpolate edges
          fishjerk_lp <- na.approx(fishjerk_lp, rule = 2)
          
          # high-pass = raw - low-pass
          fishjerk_hp_vals <- fishjerk_xts - fishjerk_lp
          
          # add to df
          current_data$act_x_lowpass <- as.numeric(fishjerk_lp[, x_col])
          current_data$act_y_lowpass <- as.numeric(fishjerk_lp[, y_col])
          current_data$act_z_lowpass <- as.numeric(fishjerk_lp[, z_col])
          
          current_data$act_x_highpass <- as.numeric(fishjerk_hp_vals[, x_col])
          current_data$act_y_highpass <- as.numeric(fishjerk_hp_vals[, y_col])
          current_data$act_z_highpass <- as.numeric(fishjerk_hp_vals[, z_col])
          
          # return processed data
          current_data
      }
        
        stopCluster(core.cl)
        
      } else {
        # SEQUENTIAL PROCESSING 
        cat("  Using sequential processing...\n")
        fishjerk_hp <- list()
        
        for(i in 1:length(fishjerk_list)) {
          cat("  Processing event", i, "of", length(fishjerk_list), "\r")
          
          current_data <- fishjerk_list[[i]]
          
          # create timestamps
          if(timestamp_col %in% colnames(current_data)) {
            base_time <- as.POSIXct(current_data[[timestamp_col]][1],
                                    format = "%Y-%m-%d %H:%M")
            current_data$date_time <- base_time + current_data[[duration_col]]
          } else {
            current_data$date_time <- as.POSIXct("2020-01-01") + current_data[[duration_col]]
          }
          
          # skip if too few observations
          if(nrow(current_data) < highpass_smooth) {
            warning("Event ", i, " has too few observations (", nrow(current_data), 
                    "). Skipping.")
            next
          }
          
          # create xts object
          fishjerk_xts <- xts(current_data[, c(x_col, y_col, z_col)],
                              order.by = current_data$date_time)
          
          # low-pass filter (smoothed version)
          fishjerk_lp <- rollmean(fishjerk_xts, k = highpass_smooth, 
                                  align = 'center', fill = NA)
          
          # interpolate edges
          fishjerk_lp <- na.approx(fishjerk_lp, rule = 2)
          
          # high-pass = raw - low-pass
          fishjerk_hp_vals <- fishjerk_xts - fishjerk_lp
          
          # add to df
          current_data$act_x_lowpass <- as.numeric(fishjerk_lp[, x_col])
          current_data$act_y_lowpass <- as.numeric(fishjerk_lp[, y_col])
          current_data$act_z_lowpass <- as.numeric(fishjerk_lp[, z_col])
          current_data$act_x_highpass <- as.numeric(fishjerk_hp_vals[, x_col])
          current_data$act_y_highpass <- as.numeric(fishjerk_hp_vals[, y_col])
          current_data$act_z_highpass <- as.numeric(fishjerk_hp_vals[, z_col])
          
          fishjerk_hp[[i]] <- current_data
        }
        
        cat("\n")
      }
      
      # remove NULL entries
      fishjerk_hp <- fishjerk_hp[!sapply(fishjerk_hp, is.null)]
      
      # combine into df
      if(length(fishjerk_hp) > 0) {
        data_hp <- do.call(rbind, fishjerk_hp)
        rownames(data_hp) <- NULL
        
        # store processed data
        data_processed$highpass <- data_hp
        
        cat("High-pass filter applied to", nrow(data_hp), "observations\n")
        
        # create high-pass data lists
        jerk_list_x_hp <- data_hp %>%
          split(.[[fish_id_col]]) %>%
          lapply(function(x) x$act_x_highpass)
        
        jerk_list_y_hp <- data_hp %>%
          split(.[[fish_id_col]]) %>%
          lapply(function(x) x$act_y_highpass)
        
        jerk_list_z_hp <- data_hp %>%
          split(.[[fish_id_col]]) %>%
          lapply(function(x) x$act_z_highpass)
        
        # X axis high-pass
        cat("Computing X-axis high-pass DTW...\n")
        dtw_matrices$x_highpass <- compute_dtw_axis(jerk_list_x_hp, "X_HP")
        
        # Y axis high-pass
        cat("Computing Y-axis high-pass DTW...\n")
        dtw_matrices$y_highpass <- compute_dtw_axis(jerk_list_y_hp, "Y_HP")
        
        # Z axis high-pass
        cat("Computing Z-axis high-pass DTW...\n")
        dtw_matrices$z_highpass <- compute_dtw_axis(jerk_list_z_hp, "Z_HP")
        
        # multivariate high-pass
        if("multivariate" %in% methods) {
          cat("Computing multivariate high-pass DTW...\n")
          
          dtw_matrix_mv_hp <- dtw_matrix_template
          
          pb <- txtProgressBar(min = 0, max = n_fish * (n_fish - 1) / 2, style = 3)
          counter <- 0
          
          for(i in 1:n_fish) {
            for(j in i:n_fish) {
              if(i != j) {
                s1 <- cbind(jerk_list_x_hp[[i]], jerk_list_y_hp[[i]], jerk_list_z_hp[[i]])
                s2 <- cbind(jerk_list_x_hp[[j]], jerk_list_y_hp[[j]], jerk_list_z_hp[[j]])
                
                cost_matrix <- proxy::dist(s1, s2, method = "Euclidean")
                alignment <- dtw(cost_matrix)
                
                dist_val <- if(normalize) {
                  alignment$normalizedDistance
                } else {
                  alignment$distance
                }
                
                dtw_matrix_mv_hp[i, j] <- dist_val
                dtw_matrix_mv_hp[j, i] <- dist_val
                
                counter <- counter + 1
                setTxtProgressBar(pb, counter)
              }
            }
          }
          close(pb)
          
          dtw_matrices$multivariate_highpass <- dtw_matrix_mv_hp
        }
        
        cat("High-pass DTW complete!\n")
        
      } else {
        warning("High-pass filter failed for all events.")
      }
    }
  }
  
  # ============================================================================
  # VISUALIZATIONS
  # ============================================================================
  
  if(create_plots) {
    cat("\n=== CREATING VISUALIZATIONS ===\n")
    
    # Plot 1: raw time series
    cat("Creating raw time series plot...\n")
    
    data_long <- data %>%
      select(all_of(c(fish_id_col, duration_col, x_col, y_col, z_col))) %>%
      pivot_longer(cols = all_of(c(x_col, y_col, z_col)),
                   names_to = "axis",
                   values_to = "value")
    
    p1 <- ggplot(data_long, aes(x = .data[[duration_col]], y = value, color = axis)) +
      geom_line(alpha = 0.8, linewidth = 0.5) +
      facet_wrap(as.formula(paste("~", fish_id_col)), ncol=3) +
      scale_color_manual(
        values = c(setNames(axis_colors["x"], x_col),
                   setNames(axis_colors["y"], y_col),
                   setNames(axis_colors["z"], z_col)),
        labels = c("X-axis", "Y-axis", "Z-axis"),
        name = ""
      ) +
      labs(
        x = "Time (s)",
        y = expression("Mean change in acceleration units")
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        legend.box.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank() 

      )
    
    plots$raw_timeseries <- p1
    
    if(!is.null(fig_output_dir)) {
      ggsave(file.path(fig_output_dir, "raw_timeseries_all.png"),
             plot = p1, width = plot_width, height = plot_height, dpi = 300)
    }
    
    # Plot 2: differential time series 
    if("differential" %in% methods) {
      cat("Creating differential time series plot...\n")
      
      # create differential data
      data_diff_list <- list()
      
      for(fish in fish_ids) {
        fish_data <- data %>% filter(.data[[fish_id_col]] == fish)
        
        fish_data_diff <- fish_data %>%
          mutate(
            act_x_diff = c(NA, diff(.data[[x_col]])),
            act_y_diff = c(NA, diff(.data[[y_col]])),
            act_z_diff = c(NA, diff(.data[[z_col]]))
          ) %>%
          filter(!is.na(act_x_diff))
        
        data_diff_list[[fish]] <- fish_data_diff
      }
      
      data_diff <- bind_rows(data_diff_list)
      data_processed$differential <- data_diff
      
      # create long format for plotting
      data_diff_long <- data_diff %>%
        select(all_of(c(fish_id_col, duration_col)),
               act_x_diff, act_y_diff, act_z_diff) %>%
        pivot_longer(cols = c(act_x_diff, act_y_diff, act_z_diff),
                     names_to = "axis",
                     values_to = "value")
      
      p2 <- ggplot(data_diff_long, aes(x = .data[[duration_col]], y = value, color = axis)) +
        geom_line(alpha = 0.8, linewidth = 0.5) +
        facet_wrap(as.formula(paste("~", fish_id_col)), ncol=3) +
        scale_color_manual(
          values = c(setNames(axis_colors["x"], "act_x_diff"),
                     setNames(axis_colors["y"], "act_y_diff"),
                     setNames(axis_colors["z"], "act_z_diff")),
          labels = c("X-axis", "Y-axis", "Z-axis"),
          name = ""
        ) +
        labs(
          x = "Time (s)",
          y = expression("Mean change in jerk units")
        ) +
        theme_bw() +
        theme(
          legend.position = "bottom",
          strip.text = element_text(face = "bold"),
          legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
          legend.box.background = element_rect(fill = "white", color = NA),
          panel.grid.minor = element_blank() 
          
        )
      
      plots$differential_timeseries <- p2
      
      if(!is.null(fig_output_dir)) {
        ggsave(file.path(fig_output_dir, "differential_timeseries_all.png"),
               plot = p2, width = plot_width, height = plot_height, dpi = 300)
      }
    }
    
    
    # Plot 3: high-pass comparison 
    if("highpass" %in% methods && !is.null(data_processed$highpass)) {
      
      cat("Creating high-pass filtered time series plot...\n")
      
      data_hp <- data_processed$highpass
      
      # create long format for plotting
      data_hp_long <- data_hp %>%
        select(all_of(c(fish_id_col, duration_col)),
               act_x_highpass, act_y_highpass, act_z_highpass) %>%
        pivot_longer(cols = c(act_x_highpass, act_y_highpass, act_z_highpass),
                     names_to = "axis",
                     values_to = "value")
      
      p4 <- ggplot(data_hp_long, aes(x = .data[[duration_col]], y = value, color = axis)) +
        geom_line(alpha = 0.8, linewidth = 0.5) +
        facet_wrap(as.formula(paste("~", fish_id_col)), ncol=3) +
        scale_color_manual(
          values = c(setNames(axis_colors["x"], "act_x_highpass"),
                     setNames(axis_colors["y"], "act_y_highpass"),
                     setNames(axis_colors["z"], "act_z_highpass")),
          labels = c("X-axis", "Y-axis", "Z-axis"),
          name = ""
        ) +
        labs(
          x = "Time (s)",
          y = expression("Mean change in acceleration units")
        ) +
        theme_bw() +
        theme(
          legend.position = "bottom",
          strip.text = element_text(face = "bold"),
          legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
          legend.box.background = element_rect(fill = "white", color = NA),
          panel.grid.minor = element_blank()
        )
      
      plots$highpass_timeseries <- p4
      
      if(!is.null(fig_output_dir)) {
        ggsave(file.path(fig_output_dir, "highpass_timeseries_all.png"),
               plot = p4, width = plot_width, height = plot_height, dpi = 300)
        
        
      cat("Creating high-pass comparison plot...\n")
      
      # sample one fish for comparison
      sample_fish <- fish_ids[1]
      sample_raw <- data %>% filter(.data[[fish_id_col]] == sample_fish)
      sample_hp <- data_processed$highpass %>% filter(.data[[fish_id_col]] == sample_fish)
      
      p3 <- ggplot() +
        # raw
        geom_line(data = sample_raw, aes(x = .data[[duration_col]], 
                                         y = .data[[x_col]], color = "Raw"),
                  alpha = 0.6) +
        # low-pass
        geom_line(data = sample_hp, aes(x = .data[[duration_col]], 
                                        y = act_x_lowpass, color = "Low-pass"),
                  linewidth = 1.2) +
        # high-pass
        geom_line(data = sample_hp, aes(x = .data[[duration_col]], 
                                        y = act_x_highpass, color = "High-pass")) +
        scale_color_manual(
          values = c("Raw" = "gray50", "Low-pass" = "blue", "High-pass" = "red"),
          name = ""
        ) +
        labs(
          title = paste("High-pass Filter Comparison - Example Fish:", sample_fish),
          x = "Time (s)",
          y = "Mean change in acceleration units"
        ) +
        theme_bw() +
        theme(legend.position = "bottom")
      
      plots$highpass_comparison <- p3
      
      if(!is.null(fig_output_dir)) {
        ggsave(file.path(fig_output_dir, "highpass_comparison.png"),
               plot = p3, width = 8, height = 5, dpi = 300)
      }
    }
    
    }
  }
  
  # ============================================================================
  # FINAL OUTPUT
  # ============================================================================
  
  cat("\n")
  cat("================================================================================\n")
  cat("                    DTW ANALYSIS COMPLETE!                                      \n")
  cat("================================================================================\n")
  cat("\nResults:\n")
  cat("  DTW matrices computed:", length(dtw_matrices), "\n")
  if(create_plots) {
    cat("  Plots created:", length(plots), "\n")
  }
  if(!is.null(fig_output_dir)) {
    cat("  Outputs saved to:", fig_output_dir, "\n")
  }
  cat("\n")
  
  # return results
  result <- list(
    dtw_matrices = dtw_matrices,
    data_processed = data_processed,
    plots = if(create_plots) plots else NULL,
    config = list(
      fish_id_col = fish_id_col,
      x_col = x_col,
      y_col = y_col,
      z_col = z_col,
      methods = methods,
      normalize = normalize,
      highpass_smooth = highpass_smooth,
      n_fish = n_fish
    )
  )
  
  class(result) <- c("dtw_analysis", "list")
  
  return(result)
}



################################################################################
# USAGE EXAMPLE
################################################################################

# # # Run the analysis
# dtw_results <- compute_dtw(
#   data = jerk,
#   fish_id_col = "fish_id",
#   x_col = "act_x",
#   y_col = "act_y",
#   z_col = "act_z",
#   duration_col = "duration",
#   event_id_col = "event_id",
#   methods = c("raw", "differential", "multivariate", "highpass"),
#   fig_output_dir = paste0(path,"dtw_results")
# )
# 
# View(dtw_results)
