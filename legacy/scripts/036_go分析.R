suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)   # 小鼠：如果不是小鼠，看下方说明替换
  library(ggplot2)
})

# ==============================================================================
# 0. 前置检查
# ==============================================================================
if(!exists("count_matrix")) stop("错误：请先加载 count_matrix！（行=gene，列=sample）")

count_matrix <- round(as.matrix(count_matrix))
storage.mode(count_matrix) <- "integer"
cat(">>> count_matrix 维度：", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# ==============================================================================
# 1. 样本信息：Tecovirimat / DMSO / mock
# ==============================================================================
tec_samples <- c(
  "JZ25070172.AGB6.Tecovirimat.1.P.AGB6.Tecovirimat.1.P_combined",
  "JZ25070173.AGB6.Tecovirimat.2.P.AGB6.Tecovirimat.2.P_combined",
  "JZ25070174.AGB6.Tecovirimat.3.P.AGB6.Tecovirimat.3.P_combined"
)

dmso_samples <- c(
  "JZ25070169.AGB6.MPXV_lb.1.P.AGB6.MPXV_lb.1.P_combined",
  "JZ25070170.AGB6.MPXV_lb.2.P.AGB6.MPXV_lb.2.P_combined",
  "JZ25070171.AGB6.MPXV_lb.3.P.AGB6.MPXV_lb.3.P_combined"
)

mock_samples <- c(
  "JZ25070166.AGB6.MOCK.1.P.AGB6.MOCK.1.P_combined",
  "JZ25070167.AGB6.MOCK.2.P.AGB6.MOCK.2.P_combined",
  "JZ25070168.AGB6.MOCK.3.P.AGB6.MOCK.3.P_combined"
)

sample_info <- data.frame(
  sample = c(tec_samples, dmso_samples, mock_samples),
  condition = c(rep("Tecovirimat", 3), rep("DMSO", 3), rep("mock", 3)),
  replicate = c(1:3, 1:3, 1:3),
  stringsAsFactors = FALSE
)

# 检查样本是否都在count_matrix中
missing <- setdiff(sample_info$sample, colnames(count_matrix))
if(length(missing) > 0){
  stop(paste0("错误：以下样本在 count_matrix 中未找到：\n", paste(missing, collapse = "\n")))
}

coldata <- sample_info %>%
  mutate(condition = factor(condition, levels = c("mock", "DMSO", "Tecovirimat"))) %>%
  column_to_rownames("sample")

counts_subset <- count_matrix[, rownames(coldata)]
counts_subset <- counts_subset[, rownames(coldata)]  # 保序
cat(">>> 样本数：", ncol(counts_subset), " | 组别：", paste(levels(coldata$condition), collapse = ", "), "\n")

# ==============================================================================
# 2. DESeq2：一次建模
# ==============================================================================
dds <- DESeqDataSetFromMatrix(countData = counts_subset, colData = coldata, design = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
cat(">>> 过滤后保留基因数:", nrow(dds), "\n")

dds <- DESeq(dds)

# ==============================================================================
# 3. 差异阈值
# ==============================================================================
FC_cutoff <- 1.0   # |log2FC| > 1
P_cutoff  <- 0.05  # padj < 0.05

# ==============================================================================
# 4. 工具函数：从某个对比提取 Up/Down 基因，并做 GO compareCluster
#    - 这里每个对比只有一个“基因集”，所以 compareCluster 不是必须；
#      但你明确要 compareCluster，我按“Up”和“Down”分别跑一次 enrichGO（更合适）。
#      如果你未来要跨多个对比一起 compareCluster，我也给你封装好了（见第5部分）。
# ==============================================================================

get_deg_lists <- function(dds, groupA, groupB, FC_cutoff=1, P_cutoff=0.05){
  res <- results(dds, contrast = c("condition", groupA, groupB))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(change = case_when(
      padj < P_cutoff & log2FoldChange >  FC_cutoff ~ "Up",
      padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down",
      TRUE ~ "NS"
    )) %>%
    arrange(padj)
  
  up_genes   <- res_df$gene[res_df$change == "Up"]
  down_genes <- res_df$gene[res_df$change == "Down"]
  
  list(res_df=res_df, up=unique(up_genes), down=unique(down_genes))
}

run_go_enrich <- function(gene_vec, title_prefix,
                          OrgDb=org.Mm.eg.db, keyType="SYMBOL",
                          ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.1,
                          showCategory=10,
                          out_prefix="GO"){
  
  if(length(gene_vec) < 5){
    message(">>> ", title_prefix, "：基因数过少（<5），跳过 GO。")
    return(NULL)
  }
  
  ego <- enrichGO(
    gene          = gene_vec,
    OrgDb         = OrgDb,
    keyType       = keyType,
    ont           = ont,
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff,
    readable      = TRUE
  )
  
  if(is.null(ego) || nrow(as.data.frame(ego)) == 0){
    message(">>> ", title_prefix, "：GO结果为空。")
    return(NULL)
  }
  
  # 保存结果表
  write.csv(as.data.frame(ego), paste0(out_prefix, "_", title_prefix, "_enrichGO.csv"), row.names = FALSE)
  
  # dotplot
  p <- dotplot(ego, showCategory = showCategory, title = paste0(title_prefix, " GO-", ont)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0(out_prefix, "_", title_prefix, "_dotplot.pdf"), p, width = 9, height = 6)
  print(p)
  
  return(ego)
}

# ==============================================================================
# 5. 你要的“分组之间分析 + GO”：三组两两比较
# ==============================================================================
comparisons <- list(
  c("Tecovirimat", "DMSO"),
  c("Tecovirimat", "mock"),
  c("DMSO", "mock")
)

all_deg <- list()
all_up  <- list()
all_down<- list()

for(cmp in comparisons){
  A <- cmp[1]; B <- cmp[2]
  tag <- paste0(A, "_vs_", B)
  message("\n==================== ", tag, " ====================")
  
  deg <- get_deg_lists(dds, A, B, FC_cutoff, P_cutoff)
  all_deg[[tag]] <- deg$res_df
  
  # 保存DESeq2结果
  write.csv(deg$res_df, paste0("DESeq2_", tag, "_all.csv"), row.names = FALSE)
  
  # 保存Up/Down gene list（文本）
  write.table(deg$up,   file = paste0("DEG_", tag, "_Up_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(deg$down, file = paste0("DEG_", tag, "_Down_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  message(">>> ", tag, "：Up=", length(deg$up), " | Down=", length(deg$down))
  
  # 分别做 GO（Up / Down）
  all_up[[tag]] <- deg$up
  all_down[[tag]] <- deg$down
  
  run_go_enrich(deg$up,   title_prefix = paste0(tag, "_Up"),
                out_prefix = "GO")
  run_go_enrich(deg$down, title_prefix = paste0(tag, "_Down"),
                out_prefix = "GO")
}

message("\n>>> 全部完成：DESeq2结果 + Up/Down基因列表 + GO富集dotplot（PDF）已输出。")

# ==============================================================================
# 6. （可选）如果你想“跨多个对比”一起 compareCluster（真正发挥 compareCluster 的意义）
#    例如：把每个对比的 Up 基因放进一个列表，一次性 compareCluster 比较
# ==============================================================================
# 上调 compareCluster（跨比较）
if(length(all_up) > 1 && any(sapply(all_up, length) >= 5)){
  message("\n>>> compareCluster：跨对比的 Up 基因 GO-BP 比较...")
  ck_up <- compareCluster(
    geneCluster   = all_up,
    fun           = "enrichGO",
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.1,
    readable      = TRUE
  )
  
  if(!is.null(ck_up) && nrow(as.data.frame(ck_up)) > 0){
    write.csv(as.data.frame(ck_up), "compareCluster_Up_GO_BP.csv", row.names = FALSE)
    p_ck_up <- dotplot(ck_up, showCategory = 10, title = "compareCluster Up: GO-BP") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("compareCluster_Up_GO_BP_dotplot.pdf", p_ck_up, width = 10, height = 6)
    print(p_ck_up)
  } else {
    message(">>> compareCluster Up：结果为空，跳过绘图。")
  }
}

# 下调 compareCluster（跨比较）
if(length(all_down) > 1 && any(sapply(all_down, length) >= 5)){
  message("\n>>> compareCluster：跨对比的 Down 基因 GO-BP 比较...")
  ck_down <- compareCluster(
    geneCluster   = all_down,
    fun           = "enrichGO",
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.1,
    readable      = TRUE
  )
  
  if(!is.null(ck_down) && nrow(as.data.frame(ck_down)) > 0){
    write.csv(as.data.frame(ck_down), "compareCluster_Down_GO_BP.csv", row.names = FALSE)
    p_ck_down <- dotplot(ck_down, showCategory = 10, title = "compareCluster Down: GO-BP") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("compareCluster_Down_GO_BP_dotplot.pdf", p_ck_down, width = 10, height = 6)
    print(p_ck_down)
  } else {
    message(">>> compareCluster Down：结果为空，跳过绘图。")
  }
}










