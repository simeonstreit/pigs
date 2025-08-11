# utils.R — Project helper functions
# -----------------------------------
# Load this file with: source("code/utils.R")


# 1. metaphlan to phyloseq function

metaphlan_to_phyloseq <- function(mpa,
                                  metadata = NULL,
                                  version = 4,
                                  verbose = TRUE,
                                  tax_lvl = "Species") {
  
  if (version == 4) {
    if (is.character(mpa)) {
      # load raw metaphlan data
      mpa <- data.table::fread(mpa) %>%
        as.data.frame()
    }
  }
  
  if (version == 3) {
    if (is.character(mpa)) {
      # load raw metaphlan data
      mpa <- data.table::fread(mpa, skip = 1) %>%
        as.data.frame()
    }
  }
  
  # find for each row, to which depth of taxonomy it arrives (as integers)
  tax_lengths <- mpa$clade_name %>%
    strsplit("|", fixed = TRUE) %>%
    sapply(length)
  
  # get integer equivalent of the taxonomic level we want
  tax_lvl_int <- match(tax_lvl, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB"))
  
  # subset the otu table, keep UNASSIGNED and selected taxonomic level
  otu_cleaned <- mpa[tax_lengths == tax_lvl_int, 3:ncol(mpa)]
  
  if (!is.null(metadata)) {
    inters_names <- intersect(colnames(otu_cleaned), rownames(metadata))
    
    if (verbose) {
      if (length(colnames(otu_cleaned)[!(colnames(otu_cleaned) %in% inters_names)]) != 0) {
        cat("Metaphlan table samples lost: ")
        cat(colnames(otu_cleaned)[!(colnames(otu_cleaned) %in% inters_names)], sep = " ")
        cat("\n\n")
      }
      
      if (length(rownames(metadata)[!(rownames(metadata) %in% inters_names)]) != 0) {
        cat("Metadata table samples lost: ")
        cat(rownames(metadata)[!(rownames(metadata) %in% inters_names)], sep = " ")
        cat("\n\n")
      }
    }
    
    otu_cleaned <- otu_cleaned[, inters_names]
    metadata_cleaned <- metadata[match(inters_names, rownames(metadata)), ]
    
  } else {
    metadata_cleaned <- metadata
  }
  
  # --- updated taxonomy parsing ---
  all_tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB")
  
  taxonomy_tab <- mpa[tax_lengths == tax_lvl_int, 1] %>%
    as.data.frame() %>%
    set_names("clade_name") %>%
    tidyr::separate(
      col = "clade_name",
      into = all_tax_ranks,
      sep = "\\|",
      fill = "right",
      remove = FALSE
    ) %>%
    dplyr::select(all_of(all_tax_ranks[1:tax_lvl_int])) %>%
    as.matrix()
  
  # rename the UNKNOWN clade
  taxonomy_tab[1, ] <- rep("UNKNOWN", tax_lvl_int)
  
  # assign rownames
  rownames(taxonomy_tab) <- taxonomy_tab[, tax_lvl_int]
  rownames(otu_cleaned) <- rownames(taxonomy_tab)
  
  profiles <- phyloseq(
    otu_table(as.matrix(otu_cleaned), taxa_are_rows = TRUE),
    tax_table(taxonomy_tab),
    sample_data(metadata_cleaned, errorIfNULL = FALSE)
  )
  
  return(profiles)
}


# 2. Volcano plot function
exp_volcano <- function(de_table, comp){
  p <- EnhancedVolcano(de_table,
                       lab = de_table$SGB,
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = paste("Differential analysis of ", comp, sep=""),
                       pCutoff = 0.1,
                       FCcutoff = 2.0,
                       pointSize = 3.0,
                       labSize = 4.0,
                       col=c('black', 'black', 'black', 'red3'),
                       colAlpha = 0.8,
                       legendPosition = 'none',
                       drawConnectors = TRUE,
                       widthConnectors = 0.75,
                       subtitle = NULL,
                       axisLabSize=13.0)
  return(p)
}

# 3. filtering function (minimal 10 read counts and 10% prevalence per group)
  calculate_prevalence_per_group <- function(ps, min_count = 10, group_var = "animal")
  {
    
    
    # Get OTU table and sample data
    otu <- phyloseq::otu_table(ps) |> as.matrix()
    # Inverted logic here to have samples as rows (as in tidy data)
    if (phyloseq::taxa_are_rows(ps)) {
      otu <- t(otu)
    }
    
    sampdat <- data.frame(phyloseq::sample_data(ps))
    
    
    assertthat::assert_that(
      group_var %in% colnames(sampdat),
      msg = paste("Group variable", group_var, "not found in sample data.")
    )
    
    assertthat::assert_that(all(rownames(sampdat) == rownames(otu)),
                            msg = "Sample names in sample data and OTU table do not match."
    )
    
    
    group_counts <- as.vector(table(sampdat[[group_var]]))
    
    # Step 1: Set low-count values to NA
    presence_absence <- (otu>=min_count) * 1
    
    presence_absence |> rowsum(group = sampdat[[group_var]]) / group_counts -> frequency_per_group
    
    return(frequency_per_group)
    
    
  }
  
  
  
  
  
  
  
  
  filter_species_by_group_prevalence <- function(ps, min_count = 10, min_prevalence = 0.1, group_var = "animal") {
    
    
    frequency_per_group <- calculate_prevalence_per_group(ps, min_count, group_var)
    
    
    keep_taxa = colAnys(frequency_per_group> min_prevalence)
    
    
    
    
    message(sprintf("Keeping %d taxa out of %d based on group prevalence criteria.", 
                    sum(keep_taxa), length(taxa_names(ps)) 
    ))
    
    
    
    
    
    # Step 3: Filter phyloseq object to keep selected taxa
    ps_filtered <- phyloseq::prune_taxa(keep_taxa, ps)
    return(ps_filtered)
  }