#create a function to get the cell abundances
getCellAbundances <- function(spe_obj, CN_column, imageID = "imageID", cellType = "cellType") {
  
  cell_abundances <- colData(spe) %>% 
    as.data.frame() %>%
    group_by(!!sym(imageID), !!sym(CN_column), !!sym(cellType)) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(col_name = paste(!!sym(CN_column), !!sym(cellType), sep = "_")) %>%
    pivot_wider(names_from = col_name, values_from = count, values_fill = 0) %>%
    select(-!!sym(CN_column), -!!sym(cellType)) %>% # Ensure region and cellType columns are not included
    group_by(!!sym(imageID)) %>%
    summarise(across(everything(), sum, na.rm = TRUE)) 
  
  
  # Convert counts to proportions for each imageID
  cell_abundances <- cell_abundances %>%
    rowwise() %>%
    mutate(across(where(is.numeric), ~ . / sum(c_across(where(is.numeric)), na.rm = TRUE))) %>%
    ungroup()  %>%
    column_to_rownames(imageID)
  
  return(cell_abundances)
}



calculate_NND <- function(pp_list, take_mean = FALSE) {
  #Then we will calculate the NND. We have to use two functions. 
  #One that finds the distance and one that tells us the index of the cell that corresponds to the first NN. 
  #We do this so we can map the cell types to the first NN.
  #We can aso specify if we want to summarize the result by taking the mean distance between each pair or keep all
  
  # Calculate the NND distances
  nn_contacts <- lapply(pp_list, function(pp) {
    nndist(pp, by = factor(pp$marks), k=1)  
  })
  
  # Get the index of the first NND
  nn_contacts_which <- lapply(pp_list, function(pp) {
    nnwhich(pp, by = factor(pp$marks), k=1)  
  })
  
  nnd_per_pair <- lapply(seq_along(pp_list), function(i) {
    
    # Extract data for the current dataset
    ID <- samples[i]
    ppp_data <- pp_list[[i]]
    nnd_matrix <- nn_contacts[[i]]          # Matrix of nearest neighbor distances
    cell_types <- colnames(nnd_matrix)      # Extract cell type names
    
    # Use a list to collect data frames, then combine once
    result_list <- list()
    idx <- 1
    
    # Loop through all pairs of cell types
    for (cell_type_1 in cell_types) {
      for (cell_type_2 in cell_types) {
        # Get distances from cell_type_1 to cell_type_2
        distances <- nnd_matrix[, cell_type_2]
        
        # Remove NAs
        distances <- distances[!is.na(distances)]
        
        if (length(distances) > 0) {
          if (take_mean) {
            # Take the mean distance
            result_list[[idx]] <- data.frame(
              imageID = as.character(ID),
              cell_type_1 = cell_type_1,
              cell_type_2 = cell_type_2,
              distance = mean(distances)
            )
            idx <- idx + 1
          } else {
            # Keep individual distances
            result_list[[idx]] <- data.frame(
              imageID = as.character(ID),
              cell_type_1 = cell_type_1,
              cell_type_2 = cell_type_2,
              distance = distances
            )
            idx <- idx + 1
          }
        }
      }
    }
    
    # Combine all results for this image at once
    return(do.call(rbind, result_list))
  })
  
  # Combine results across all datasets
  NND.df <- do.call(rbind, nnd_per_pair)
  
  NND.df <- NND.df %>%
    left_join(colData(spe) %>% as.data.frame() %>% 
                dplyr::select(imageID, tumour_type) %>% 
                unique(), by = join_by("imageID")) %>%
    dplyr::filter(!is.na(tumour_type)) %>%
    mutate(cell_pair = paste0(cell_type_1, "_", cell_type_2))
  
  return(NND.df)
}


calculate_NNC <- function(pp_list, r) {
  # Calculate pairwise distances for each cell type
  nn_cell_types <- lapply(seq_along(pp_list), function(j) {
    ID <- samples[j]
    ppp_data <- pp_list[[j]]
    distances <- pairdist(ppp_data)
    
    # Set a distance threshold
    threshold <- r  # Define threshold
    
    pairwise_nnd <- data.frame()
    
    # Create a dataframe to store the counts
    cell_types <- unique(marks(ppp_data))  # Get unique cell types
    neighbor_counts <- data.frame(matrix(0, nrow = npoints(ppp_data), ncol = length(cell_types)))
    colnames(neighbor_counts) <- cell_types  # Name columns by cell types
    
    # Calculate counts of each cell type within the threshold for each point
    for (i in seq_len(npoints(ppp_data))) {
      within_threshold <- which(distances[i, ] <= threshold & distances[i, ] > 0)  # Exclude self
      neighbor_types <- marks(ppp_data)[within_threshold]  # Cell types within threshold
      
      # Count occurrences of each cell type
      cell_type_counts <- table(neighbor_types)
      
      # Update the counts in the dataframe
      neighbor_counts[i, names(cell_type_counts)] <- as.numeric(cell_type_counts)
    }
    
    # Combine results with the original cell types
    neighbor_summary <- rbind(pairwise_nnd, data.frame(
      imageID = ID,
      cell_type = marks(ppp_data),
      neighbor_counts
    )
    )
    
    neighbor_summary$total_neighbors <- rowSums(neighbor_counts)
    
    # Return the counts data frame for this dataset
    return(neighbor_summary)
  })
  
  
  #some images do not have all cell types so we need to fill them with 0s
  all_columns <- unique(unlist(lapply(nn_cell_types, names)))
  # Add missing columns to each data frame
  nn_cell_types <- lapply(nn_cell_types, function(df) {
    missing_columns <- setdiff(all_columns, names(df))
    for (col in missing_columns) {
      df[[col]] <- 0  # Add missing column with 0 values
    }
    return(df[, all_columns])  # Reorder columns to match the complete set
  })
  
  #Combine results across all datasets
  NNC.df <- do.call(rbind, nn_cell_types)

  return(NNC.df)
}


#Lest
calculate_LCross <- function(pp_list, r_test) {
  LCross_res <- lapply(seq_along(pp_list), function(i) {
    
    # Extract data for the current dataset
    ID <- samples[i]  # Image ID
    ppp_data <- pp_list[[i]]      # Current point pattern
    cell_types <- unique(marks(ppp_data))  # Extract cell type names
    
    # Initialize a data frame to store results
    LCross_df <- data.frame()
    LCross_all <- list()
    # Loop through all pairs of cell types
    for (cell_type_1 in cell_types) {
      for (cell_type_2 in cell_types) {
        
        # Skip self-pairs (otherwise extreme values will be present)
        if (cell_type_1 == cell_type_2) {
          next
        }
        
        # Calculate the pcfcross for cell_type_1 and cell_type_2
        LCross_res <- Lcross(ppp_data, i = cell_type_1, j = cell_type_2, r = r_test)
        
        # Check if the result contains valid data
        if (length(LCross_res$r) > 0 && length(LCross_res$iso) > 0) {
          # Extract distances (r) and g(r) values from LCross_res
          r <- LCross_res$r
          iso <- LCross_res$iso
          
          centered_iso <- iso - r
          
          # Combine the results into a data frame
          result_df <- data.frame(
            imageID = rep(as.character(ID), length(r)),    # Repeat ID for all r
            cell_type_1 = rep(as.character(cell_type_1), length(r)),  # Repeat cell_type_1
            cell_type_2 = rep(as.character(cell_type_2), length(r)),  # Repeat cell_type_2
            r = r,  # Distances (r values)
            iso = centered_iso     # isotropic corrected values
          )
          LCross_all[[i]] <- LCross_res
          
          # Append to the final data frame
          LCross_df <- rbind(LCross_df, result_df)
        } else {
          message("Skipping pair: ", cell_type_1, " vs ", cell_type_2, " in image ", ID, " due to empty result.")
        }
      }
    }
    
    # Return the pairwise LCross data frame for this dataset
    return(LCross_df)
  })
  
  LCross_res <- do.call(rbind, LCross_res)
  return(LCross_res)
}



#pair correlation function
calculate_pcf <- function(pp_list, r_test) {
  pcfcross_res <- lapply(seq_along(pp_list), function(i) {
    
    # Extract data for the current dataset
    ID <- samples[i]  # Image ID
    ppp_data <- pp_list[[i]]      # Current point pattern
    cell_types <- unique(marks(ppp_data))  # Extract cell type names
    
    # Initialize a data frame to store results
    pcfcross_df <- data.frame()
    
    # Loop through all pairs of cell types
    for (cell_type_1 in cell_types) {
      for (cell_type_2 in cell_types) {
        
        
        # Calculate the pcfcross for cell_type_1 and cell_type_2
        
        pcfcross_result <- pcfcross(ppp_data, i = cell_type_1, j = cell_type_2, r = r_test)
        
        
        # Extract distances (r) and g(r) values from pcfcross_result
        r <- pcfcross_result$r
        iso <- pcfcross_result$iso
        
        
        # Combine the results into a data frame
        result_df <- data.frame(
          imageID = rep(as.character(ID), length(r)),    # Repeat ID for all r
          cell_type_1 = rep(as.character(cell_type_1), length(r)),  # Repeat cell_type_1
          cell_type_2 = rep(as.character(cell_type_2), length(r)),  # Repeat cell_type_2
          r = r,  # Distances (r values)
          iso = iso     # isotropic corrected values
        )
        
        # Append to the final data frame
        pcfcross_df <- rbind(pcfcross_df, result_df)
      }
    }
    
    # Return the pairwise pcfcross data frame for this dataset
    return(pcfcross_df)
  })
  
  #Combine results across all datasets
  pcfcross_res <- do.call(rbind, pcfcross_res)
  return(pcfcross_res)
}



format_spatial_result <- function(result.df, radius) {
  #function for formatting cross L function and pcf
  result.df <- result.df %>% 
    filter(r == radius) %>%
    mutate(cell_pair = paste0(cell_type_1, "_", cell_type_2))
  
  
  result.mat <-  result.df %>%
    group_by(imageID, cell_pair) %>%
    summarize(iso = mean(iso, na.rm = TRUE)) %>%
    ungroup()
  
  
  result.mat <- result.mat %>%
    pivot_wider(names_from = cell_pair, values_from = iso) %>%
    mutate(across(everything(), ~ replace(., is.na(.), 0))) %>%
    column_to_rownames("imageID")
  
  return(result.mat)
}