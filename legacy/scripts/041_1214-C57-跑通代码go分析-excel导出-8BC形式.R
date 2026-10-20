



# 加载必要的包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)    # 用于富集分析
library(org.Mm.eg.db)       # 小鼠基因注释数据库
library(DOSE)               # 疾病本体分析
library(openxlsx)           # 新增：用于Excel输出
library(viridis)            # 新增：用于颜色方案

# 设置输出目录
output_dir <- "D:/Desktop/C57_MPXV_analysis_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}






# 1. 读取并整理计数数据 -------------------------------------------------
counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/C57-counts_gene.txt", 
                     header = TRUE, sep = "\t", comment.char = "#", row.names = 1)

# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 6:20]

# 重命名列名为简洁的样本名
colnames(count_matrix) <- c("C57_10d_1_P", "C57_10d_2_P", "C57_10d_3_P",
                            "C57_3d_1_P", "C57_3d_2_P", "C57_3d_3_P",
                            "C57_5d_1_P", "C57_5d_2_P", "C57_5d_3_P",
                            "C57_7d_1_P", "C57_7d_2_P", "C57_7d_3_P",
                            "c57_mock_1_P", "c57_mock_2_P", "c57_mock_3_P")

print("重命名后的列名:")
print(colnames(count_matrix))








# 2. 创建样本信息数据框 -------------------------------------------------
sample_info <- data.frame(
  sample = colnames(count_matrix),
  timepoint = factor(c(rep("10d", 3), rep("3d", 3), rep("5d", 3), 
                       rep("7d", 3), rep("mock", 3)),
                     levels = c("mock", "3d", "5d", "7d", "10d")),
  condition = factor(c(rep("MPXV_lb", 12), rep("mock", 3)),  # 修正：10d只有1个样本
                     levels = c("mock", "MPXV_lb")),
  group_label = c(rep("MPXV_lb", 12), rep("MOCK", 3))
)

rownames(sample_info) <- sample_info$sample

print("样本信息:")
print(sample_info)

# 3. 创建DESeq2对象 ----------------------------------------------------
# 确保计数矩阵只包含数值数据
count_matrix_numeric <- as.matrix(count_matrix)
mode(count_matrix_numeric) <- "numeric"

# 创建DESeq2数据集
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData = sample_info,
  design = ~ condition
)

print("DESeq2对象创建成功！")
print(paste("基因数:", nrow(dds)))
print(paste("样本数:", ncol(dds)))

# 4. 过滤低表达基因 ----------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3  # 至少在3个样本中表达量≥10
dds <- dds[keep, ]
print(paste("过滤后保留", nrow(dds), "个基因"))

# 5. 运行DESeq2分析 ----------------------------------------------------
dds <- DESeq(dds)

# 6. 按时间点进行差异分析 ----------------------------------------------
timepoints <- c("3d", "5d", "7d")  # 移除10d，因为只有一个重复
results_list <- list()
deg_list <- list()

# 定义差异基因阈值
FC_cutoff <- 1.0  # log2FoldChange > 1 (即2倍差异)
P_cutoff <- 0.05  # padj < 0.05




for(tp in timepoints) {
  print(paste("分析时间点:", tp))
  
  # 选择当前时间点和mock的样本
  tp_samples <- sample_info$sample[sample_info$timepoint %in% c(tp, "mock")]
  dds_sub <- dds[, tp_samples]
  
  # 确保condition因子水平正确（mock作为参考）
  dds_sub$condition <- relevel(dds_sub$condition, ref = "mock")
  
  # 运行DESeq2
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("condition", "MPXV_lb", "mock"), alpha = 0.05)
  
  # 保存完整结果
  results_list[[tp]] <- res
  
  # 提取显著差异基因并标记上下调
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      change = case_when(
        padj < P_cutoff & log2FoldChange > FC_cutoff ~ "Up-regulated",
        padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down-regulated",
        TRUE ~ "Not Significant"
      ),
      change = factor(change, levels = c("Down-regulated", "Not Significant", "Up-regulated")),
      timepoint = tp  # 添加时间点信息
    ) %>%
    arrange(padj)
  
  deg_list[[tp]] <- res_df
  print(paste(tp, "vs mock: 发现", sum(res_df$change != "Not Significant"), "个显著DEGs"))
}

# ===============================================
# 7. 保存差异表达分析结果到文件
# ===============================================
message("保存差异表达分析结果...")

# 创建Excel工作簿
wb_deg <- createWorkbook()

# 保存每个时间点的详细差异表达结果
for(tp in timepoints) {
  if(!is.null(deg_list[[tp]])) {
    sheet_name <- paste0("DEGs_", tp)
    addWorksheet(wb_deg, sheet_name)
    writeData(wb_deg, sheet = sheet_name, deg_list[[tp]], rowNames = FALSE)
  }
}

# 合并所有时间点的DEGs并保存
all_degs <- bind_rows(deg_list)
addWorksheet(wb_deg, "All_DEGs_Combined")
writeData(wb_deg, sheet = "All_DEGs_Combined", all_degs, rowNames = FALSE)

# 保存显著DEGs汇总
sig_degs_summary <- all_degs %>%
  filter(change != "Not Significant") %>%
  group_by(timepoint, change) %>%
  summarise(Count = n(), .groups = "drop")

addWorksheet(wb_deg, "Significant_DEGs_Summary")
writeData(wb_deg, sheet = "Significant_DEGs_Summary", sig_degs_summary, rowNames = FALSE)

# 保存DEGs Excel文件
deg_excel_file <- file.path(output_dir, "8bc_MPXV_DEGs_Analysis.xlsx")
saveWorkbook(wb_deg, deg_excel_file, overwrite = TRUE)
message(paste("差异表达结果已保存到:", deg_excel_file))

# 同时保存CSV格式
write.csv(all_degs, file.path(output_dir, "8bc_MPXV_All_DEGs.csv"), row.names = FALSE)
message(paste("所有DEGs已保存为CSV:", file.path(output_dir, "8bc_MPXV_All_DEGs.csv")))

# ===============================================
# 8. 提取差异表达基因（DEGs）列表
# ===============================================
up_gene_list <- list()
down_gene_list <- list()
logFC_cutoff <- 1  # 定义 logFC 阈值

for(tp in timepoints) {
  res_df <- deg_list[[tp]]
  
  if(nrow(res_df) > 0) {
    # 提取上调基因 (log2FoldChange > 1)
    up_genes <- res_df$gene[res_df$log2FoldChange > logFC_cutoff & res_df$padj < P_cutoff]
    if(length(up_genes) > 0) {
      up_gene_list[[tp]] <- up_genes
    }
    
    # 提取下调基因 (log2FoldChange < -1)
    down_genes <- res_df$gene[res_df$log2FoldChange < -logFC_cutoff & res_df$padj < P_cutoff]
    if(length(down_genes) > 0) {
      down_gene_list[[tp]] <- down_genes
    }
  }
}

# 检查基因列表是否为空
if (length(up_gene_list) == 0) {
  warning("在任何时间点均未找到显著上调基因。")
}
if (length(down_gene_list) == 0) {
  warning("在任何时间点均未找到显著下调基因。")
}

# ===============================================
# 9. 运行 GO 富集分析
# ===============================================
message("进行GO富集分析...")

# 创建GO结果工作簿
wb_go <- createWorkbook()
# ===============================================
# 9. 运行 GO 富集分析
# ===============================================
message("进行GO富集分析...")

# 创建GO结果工作簿
wb_go <- createWorkbook()

# 定义GO分析参数
pvalue_cutoff <- 0.1      # p值阈值，可调节（0.01-0.1）
qvalue_cutoff <- 0.2      # q值（调整后p值）阈值，可调节（0.05-0.2）
min_gs_size <- 10         # GO term的最小基因数，可调节（5-20）
max_gs_size <- 1000       # GO term的最大基因数，可调节（200-1000）
p_adjust_method <- "BH"   # p值校正方法，可调节（"BH", "bonferroni", "fdr"等）
ontology <- "BP"          # GO子类，可调节（"BP", "CC", "MF" 或 "ALL"）
show_category <- 20       # 显示富集条目的数量，可调节
min_enrichment <- 1.5     # 最小富集因子，可调节（1.5-3）

# 9.1 上调基因GO分析
if (length(up_gene_list) > 0) {
  message("正在对上调基因进行GO富集分析...")
  message(paste("上调基因列表数量:", length(up_gene_list)))
  
  for (tp in names(up_gene_list)) {
    message(paste("  分析时间点", tp, ":", length(up_gene_list[[tp]]), "个上调基因"))
  }
  
  # 使用compareCluster对多个时间点进行比较分析
  ck_up <- compareCluster(
    geneCluster = up_gene_list,
    fun = "enrichGO",
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = ontology,
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    pAdjustMethod = p_adjust_method,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    readable = TRUE
  )
  
  if (!is.null(ck_up) && nrow(as.data.frame(ck_up)) > 0) {
    up_go_df <- as.data.frame(ck_up)
    
    # 计算富集因子
    up_go_df <- up_go_df %>%
      mutate(
        GeneRatio_num = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
        BgRatio_num = as.numeric(sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
        EnrichmentRatio = GeneRatio_num / BgRatio_num
      ) %>%
      filter(EnrichmentRatio > min_enrichment)  # 过滤低富集因子
    
    if (nrow(up_go_df) > 0) {
      addWorksheet(wb_go, "上调基因_GO富集")
      writeData(wb_go, sheet = "上调基因_GO富集", up_go_df, rowNames = FALSE)
      
      # 绘制并保存富集条形图
      message("绘制上调基因GO富集条形图...")
      for (tp in names(up_gene_list)) {
        # 提取当前时间点的富集结果
        tp_go_df <- up_go_df %>% filter(Cluster == tp)
        
        if (nrow(tp_go_df) > 0) {
          # 选择前show_category个最显著的条目
          up_go_top <- tp_go_df %>%
            arrange(p.adjust) %>%
            head(min(show_category, nrow(tp_go_df)))
          
          # 绘制条形图
          bar_plot <- ggplot(up_go_top, 
                             aes(x = reorder(Description, -log10(p.adjust)), 
                                 y = -log10(p.adjust),
                                 fill = -log10(p.adjust))) +
            geom_bar(stat = "identity", alpha = 0.8) +
            scale_fill_viridis_c(option = "C") +
            coord_flip() +
            theme_minimal(base_size = 12) +
            labs(x = "GO Term", y = "-log10(adjusted p-value)",
                 title = paste("上调基因GO富集分析 (", tp, ")", sep = ""),
                 subtitle = paste0("Top ", nrow(up_go_top), " terms")) +
            theme(
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 10, hjust = 0.5),
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              legend.position = "none"
            )
          
          # 保存条形图
          bar_plot_file <- file.path(output_dir, paste0("上调基因_GO富集_条形图_", tp, ".pdf"))
          ggsave(bar_plot_file, plot = bar_plot, width = 12, height = 8)
          message(paste("  条形图已保存到:", bar_plot_file))
        }
      }
      
      # 绘制点图（跨时间点比较）
      message("绘制上调基因GO富集点图（跨时间点比较）...")
      if (length(unique(up_go_df$Cluster)) > 1) {
        dot_plot <- dotplot(ck_up, 
                            showCategory = min(5, length(unique(up_go_df$ID))), 
                            font.size = 10,
                            title = "上调基因GO富集分析（跨时间点比较）") +
          theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
        
        dot_plot_file <- file.path(output_dir, "上调基因_GO富集_点图_跨时间点比较.pdf")
        ggsave(dot_plot_file, plot = dot_plot, width = 12, height = 8)
        message(paste("  点图已保存到:", dot_plot_file))
      }
      
      message("上调基因GO富集分析完成")
    } else {
      message("上调基因GO富集分析：过滤后无显著富集结果")
    }
  } else {
    message("上调基因GO富集分析：无显著富集结果")
  }
} else {
  warning("上调基因列表为空，跳过上调基因GO富集分析")
}

# 9.2 下调基因GO分析
if (length(down_gene_list) > 0) {
  message("正在对下调基因进行GO富集分析...")
  message(paste("下调基因列表数量:", length(down_gene_list)))
  
  for (tp in names(down_gene_list)) {
    message(paste("  分析时间点", tp, ":", length(down_gene_list[[tp]]), "个下调基因"))
  }
  
  # 使用compareCluster对多个时间点进行比较分析
  ck_down <- compareCluster(
    geneCluster = down_gene_list,
    fun = "enrichGO",
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = ontology,
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    pAdjustMethod = p_adjust_method,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    readable = TRUE
  )
  
  if (!is.null(ck_down) && nrow(as.data.frame(ck_down)) > 0) {
    down_go_df <- as.data.frame(ck_down)
    
    # 计算富集因子
    down_go_df <- down_go_df %>%
      mutate(
        GeneRatio_num = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
        BgRatio_num = as.numeric(sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
        EnrichmentRatio = GeneRatio_num / BgRatio_num
      ) %>%
      filter(EnrichmentRatio > min_enrichment)  # 过滤低富集因子
    
    if (nrow(down_go_df) > 0) {
      addWorksheet(wb_go, "下调基因_GO富集")
      writeData(wb_go, sheet = "下调基因_GO富集", down_go_df, rowNames = FALSE)
      
      # 绘制并保存富集条形图
      message("绘制下调基因GO富集条形图...")
      for (tp in names(down_gene_list)) {
        # 提取当前时间点的富集结果
        tp_go_df <- down_go_df %>% filter(Cluster == tp)
        
        if (nrow(tp_go_df) > 0) {
          # 选择前show_category个最显著的条目
          down_go_top <- tp_go_df %>%
            arrange(p.adjust) %>%
            head(min(show_category, nrow(tp_go_df)))
          
          # 绘制条形图
          bar_plot <- ggplot(down_go_top, 
                             aes(x = reorder(Description, -log10(p.adjust)), 
                                 y = -log10(p.adjust),
                                 fill = -log10(p.adjust))) +
            geom_bar(stat = "identity", alpha = 0.8) +
            scale_fill_viridis_c(option = "D") +
            coord_flip() +
            theme_minimal(base_size = 12) +
            labs(x = "GO Term", y = "-log10(adjusted p-value)",
                 title = paste("下调基因GO富集分析 (", tp, ")", sep = ""),
                 subtitle = paste0("Top ", nrow(down_go_top), " terms")) +
            theme(
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 10, hjust = 0.5),
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              legend.position = "none"
            )
          
          # 保存条形图
          bar_plot_file <- file.path(output_dir, paste0("下调基因_GO富集_条形图_", tp, ".pdf"))
          ggsave(bar_plot_file, plot = bar_plot, width = 12, height = 8)
          message(paste("  条形图已保存到:", bar_plot_file))
        }
      }
      
      # 绘制点图（跨时间点比较）
      message("绘制下调基因GO富集点图（跨时间点比较）...")
      if (length(unique(down_go_df$Cluster)) > 1) {
        dot_plot <- dotplot(ck_down, 
                            showCategory = min(5, length(unique(down_go_df$ID))), 
                            font.size = 10,
                            title = "下调基因GO富集分析（跨时间点比较）") +
          theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
        
        dot_plot_file <- file.path(output_dir, "下调基因_GO富集_点图_跨时间点比较.pdf")
        ggsave(dot_plot_file, plot = dot_plot, width = 12, height = 8)
        message(paste("  点图已保存到:", dot_plot_file))
      }
      
      message("下调基因GO富集分析完成")
    } else {
      message("下调基因GO富集分析：过滤后无显著富集结果")
    }
  } else {
    message("下调基因GO富集分析：无显著富集结果")
  }
} else {
  warning("下调基因列表为空，跳过下调基因GO富集分析")
}

# 9.3 合并所有差异基因的GO分析
message("进行合并所有差异基因的GO分析...")
all_sig_genes <- list()

# 合并各时间点上调基因
all_up_genes <- unique(unlist(up_gene_list))
if (length(all_up_genes) > 0) {
  all_sig_genes[["All_Up"]] <- all_up_genes
  message(paste("  合并上调基因数:", length(all_up_genes)))
}

# 合并各时间点下调基因
all_down_genes <- unique(unlist(down_gene_list))
if (length(all_down_genes) > 0) {
  all_sig_genes[["All_Down"]] <- all_down_genes
  message(paste("  合并下调基因数:", length(all_down_genes)))
}

# 合并所有差异基因
all_deg_genes <- unique(c(all_up_genes, all_down_genes))
if (length(all_deg_genes) > 0) {
  all_sig_genes[["All_DEGs"]] <- all_deg_genes
  message(paste("  总差异基因数:", length(all_deg_genes)))
  
  # 对合并的差异基因进行GO富集分析
  if (length(all_sig_genes) > 0) {
    for (group_name in names(all_sig_genes)) {
      message(paste("  分析", group_name, "..."))
      genes <- all_sig_genes[[group_name]]
      
      if (length(genes) >= 5) {  # 至少5个基因才进行富集分析
        go_result <- enrichGO(
          gene = genes,
          OrgDb = org.Mm.eg.db,
          keyType = "SYMBOL",
          ont = ontology,
          pvalueCutoff = pvalue_cutoff,
          qvalueCutoff = qvalue_cutoff,
          pAdjustMethod = p_adjust_method,
          minGSSize = min_gs_size,
          maxGSSize = max_gs_size,
          readable = TRUE
        )
        
        if (!is.null(go_result) && nrow(as.data.frame(go_result)) > 0) {
          go_df <- as.data.frame(go_result)
          
          # 计算富集因子
          go_df <- go_df %>%
            mutate(
              GeneRatio_num = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
              BgRatio_num = as.numeric(sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
              EnrichmentRatio = GeneRatio_num / BgRatio_num
            ) %>%
            filter(EnrichmentRatio > min_enrichment)  # 过滤低富集因子
          
          if (nrow(go_df) > 0) {
            sheet_name <- paste0(group_name, "_GO富集")
            addWorksheet(wb_go, sheet_name)
            writeData(wb_go, sheet = sheet_name, go_df, rowNames = FALSE)
            
            # 绘制条形图
            go_top <- go_df %>%
              arrange(p.adjust) %>%
              head(min(show_category, nrow(go_df)))
            
            bar_plot <- ggplot(go_top, 
                               aes(x = reorder(Description, -log10(p.adjust)), 
                                   y = -log10(p.adjust),
                                   fill = -log10(p.adjust))) +
              geom_bar(stat = "identity", alpha = 0.8) +
              scale_fill_viridis_c(option = ifelse(group_name == "All_Up", "C", "D")) +
              coord_flip() +
              theme_minimal(base_size = 12) +
              labs(x = "GO Term", y = "-log10(adjusted p-value)",
                   title = paste(group_name, "GO富集分析"),
                   subtitle = paste0("Top ", nrow(go_top), " terms")) +
              theme(
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(size = 10, hjust = 0.5),
                axis.text = element_text(size = 10),
                axis.title = element_text(size = 12),
                legend.position = "none"
              )
            
            # 保存条形图
            bar_plot_file <- file.path(output_dir, paste0(group_name, "_GO富集_条形图.pdf"))
            ggsave(bar_plot_file, plot = bar_plot, width = 12, height = 8)
            message(paste("    条形图已保存到:", bar_plot_file))
            
            # 绘制点图
            dot_plot <- dotplot(go_result, 
                                showCategory = min(10, nrow(go_df)), 
                                font.size = 10,
                                title = paste(group_name, "GO富集分析")) +
              theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
            
            dot_plot_file <- file.path(output_dir, paste0(group_name, "_GO富集_点图.pdf"))
            ggsave(dot_plot_file, plot = dot_plot, width = 12, height = 8)
            message(paste("    点图已保存到:", dot_plot_file))
          }
        }
      } else {
        message(paste("  ", group_name, "基因数不足5个，跳过GO富集分析"))
      }
    }
  }
}

# 保存GO结果Excel文件
go_excel_file <- file.path(output_dir, "8bc_MPXV_GO_Analysis.xlsx")
saveWorkbook(wb_go, go_excel_file, overwrite = TRUE)
message(paste("GO富集分析结果已保存到:", go_excel_file))

# 保存参数设置
params_df <- data.frame(
  Parameter = c("P-value cutoff", "Q-value cutoff", "Min GSSize", "Max GSSize", 
                "pAdjustMethod", "Ontology", "Min Enrichment Ratio"),
  Value = c(pvalue_cutoff, qvalue_cutoff, min_gs_size, max_gs_size,
            p_adjust_method, ontology, min_enrichment)
)

write.csv(params_df, file.path(output_dir, "GO_Analysis_Parameters.csv"), row.names = FALSE)
message("GO分析参数已保存")

message("GO富集分析完成！")