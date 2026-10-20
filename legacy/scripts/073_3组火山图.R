# 7.2 各比较组火山图（拆分版本）

# 7.2.1 MPXV_IIb_vs_MOCK 火山图
if("MPXV_IIb_vs_MOCK" %in% names(deg_list)) {
  res_df <- deg_list[["MPXV_IIb_vs_MOCK"]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(15) %>% 
      pull(gene)
    
    p_vol_IIb_vs_MOCK <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Down-regulated" = color_down, 
                                    "Not Significant" = color_ns, 
                                    "Up-regulated" = color_up)) +
      geom_vline(xintercept = c(-FC_cutoff, FC_cutoff), 
                 linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(P_cutoff), 
                 linetype = "dashed", color = "grey50") +
      geom_text_repel(data = subset(res_df, gene %in% top_genes),
                      aes(label = gene, color = change),
                      size = 2.5, show.legend = FALSE, max.overlaps = 15) +
      labs(title = "Volcano: MPXV IIb vs MOCK", 
           x = expression(Log[2]~Fold~Change), 
           y = expression(-Log[10]~Adjusted~P-value)) +
      xlim(-13, 13) +  # 添加x轴范围限制
      theme_bw() + 
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(p_vol_IIb_vs_MOCK)
    ggsave("Volcano_MPXV_IIb_vs_MOCK.pdf", p_vol_IIb_vs_MOCK, width = 7, height = 6)
    cat("已生成 MPXV IIb vs MOCK 火山图\n")
  }
}






# 7.2.2 MPXV_Ib_vs_MOCK 火山图
if("MPXV_Ib_vs_MOCK" %in% names(deg_list)) {
  res_df <- deg_list[["MPXV_Ib_vs_MOCK"]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(15) %>% 
      pull(gene)
    
    p_vol_Ib_vs_MOCK <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Down-regulated" = color_down, 
                                    "Not Significant" = color_ns, 
                                    "Up-regulated" = color_up)) +
      geom_vline(xintercept = c(-FC_cutoff, FC_cutoff), 
                 linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(P_cutoff), 
                 linetype = "dashed", color = "grey50") +
      geom_text_repel(data = subset(res_df, gene %in% top_genes),
                      aes(label = gene, color = change),
                      size = 2.5, show.legend = FALSE, max.overlaps = 15) +
      labs(title = "Volcano: MPXV Ib vs MOCK", 
           x = expression(Log[2]~Fold~Change), 
           y = expression(-Log[10]~Adjusted~P-value)) +
      xlim(-15, 15) +  # 设置x轴范围为-15到15
      ylim(0, 200) +   # 设置y轴范围为0到200
      theme_bw() + 
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(p_vol_Ib_vs_MOCK)
    ggsave("Volcano_MPXV_Ib_vs_MOCK.pdf", p_vol_Ib_vs_MOCK, width = 7, height = 6)
    cat("已生成 MPXV Ib vs MOCK 火山图\n")
  }
}







# 7.2.3 MPXV_IIb_vs_MPXV_Ib 火山图
if("MPXV_IIb_vs_MPXV_Ib" %in% names(deg_list)) {
  res_df <- deg_list[["MPXV_IIb_vs_MPXV_Ib"]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(15) %>% 
      pull(gene)
    
    p_vol_IIb_vs_Ib <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Down-regulated" = color_down, 
                                    "Not Significant" = color_ns, 
                                    "Up-regulated" = color_up)) +
      geom_vline(xintercept = c(-FC_cutoff, FC_cutoff), 
                 linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(P_cutoff), 
                 linetype = "dashed", color = "grey50") +
      geom_text_repel(data = subset(res_df, gene %in% top_genes),
                      aes(label = gene, color = change),
                      size = 2.5, show.legend = FALSE, max.overlaps = 15) +
      labs(title = "Volcano: MPXV IIb vs MPXV Ib", 
           x = expression(Log[2]~Fold~Change), 
           y = expression(-Log[10]~Adjusted~P-value))+
      xlim(-13, 13) +  # 设置x轴范围为-15到15 +
      theme_bw() + 
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(p_vol_IIb_vs_Ib)
    ggsave("Volcano_MPXV_IIb_vs_MPXV_Ib.pdf", p_vol_IIb_vs_Ib, width = 7, height = 6)
    cat("已生成 MPXV IIb vs MPXV Ib 火山图\n")
  }
}





# 检查是否所有比较组都成功生成了火山图
generated_plots <- c()
if(exists("p_vol_IIb_vs_MOCK")) generated_plots <- c(generated_plots, "MPXV_IIb_vs_MOCK")
if(exists("p_vol_Ib_vs_MOCK")) generated_plots <- c(generated_plots, "MPXV_Ib_vs_MOCK")
if(exists("p_vol_IIb_vs_Ib")) generated_plots <- c(generated_plots, "MPXV_IIb_vs_MPXV_Ib")

cat("成功生成的火山图:", paste(generated_plots, collapse = ", "), "\n")


















