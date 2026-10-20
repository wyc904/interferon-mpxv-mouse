# --------------------------------------------------------------------------
# Module E: Visualization of GO enrichment results for merging all comparison groups (optional)
# --------------------------------------------------------------------------

# If there are GO results from multiple comparison groups, comparative analysis can be performed

# Define the specific comparison groups to include in the merged plot
target_comparisons <- c("MPXV_IIb_vs_MOCK", "MPXV_Ib_vs_MOCK")

# --- 1. Combined Upregulated GO BP Dotplot (Targeting MOCK comparisons) ---
# Filter all_ego_up_list to only include the desired comparison groups
target_ego_up_list <- all_ego_up_list[names(all_ego_up_list) %in% target_comparisons]

if(length(target_ego_up_list) > 1) {
  cat("\n>>> Generate a combined dotplot for Upregulated GO BP results (MPXV_IIb_vs_MOCK & MPXV_Ib_vs_MOCK)...\n")
  
  # The GO results of upregulated genes were combined
  ego_up_combined_MOCK <- merge_result(target_ego_up_list)
  
  if(nrow(ego_up_combined_MOCK) > 0) {
    p_compare_up_MOCK <- dotplot(ego_up_combined_MOCK, showCategory=10) + 
      ggtitle("Comparison of GO BP Enrichment - Upregulated Genes (vs. MOCK)")
    
    # Save the new combined plot
    ggsave(file.path(output_dir, "GO_BP_Upregulated_Comparison_MOCK_Groups.pdf"), 
           p_compare_up_MOCK, width = 12, height = 8)
    
    cat("Combined dotplot for upregulated genes (vs. MOCK) saved as GO_BP_Upregulated_Comparison_MOCK_Groups.pdf\n")
    # This resulting diagram is a dotplot showing the enriched GO BP terms 
    # and how they compare between the two MPXV groups vs. MOCK.
    # 
  } else {
    cat("Warning: Merged upregulated GO BP results (vs. MOCK) are empty. Skipping combined plot.\n")
  }
} else {
  cat("\nWarning: Not enough upregulated GO BP results to create the specific combined plot (need both MPXV_IIb_vs_MOCK and MPXV_Ib_vs_MOCK).\n")
}

# --- 2. Original Combined Downregulated GO BP Dotplot (Retained for completeness) ---
if(length(all_ego_down_list) > 1) {
  cat("\n>>> Generate a comparison plot of GO enrichment results in the comparison group (Downregulated)...\n")
  # Combine the GO results of downregulated genes
  ego_down_combined <- merge_result(all_ego_down_list)
  if(nrow(ego_down_combined) > 0) {
    p_compare_down <- dotplot(ego_down_combined, showCategory=10) + 
      ggtitle("Comparison of GO BP Enrichment - Downregulated Genes (All Groups)")
    ggsave(file.path(output_dir, "GO_BP_Downregulated_Comparison.pdf"), p_compare_down, width = 12, height = 8)
  }
}