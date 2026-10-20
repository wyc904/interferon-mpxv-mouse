
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

setwd("D:\\project/hhy/mouse/second spleen/")
set.seed(123)
library(ggplot2)
#########################################################################
#读取xlsx
library(readxl)

# 读取Excel文件
# path参数指定你的.xlsx文件的路径
exp <- read_excel("D:\\project/hhy/second_spleen.genes.expression.annot.xlsx")
exp<-as.data.frame(exp)
ids<-exp[,1]
length(unique(ids))
#查看每个基因对应的探针数,tail查看后六个数，table生产频数统计表
tail(sort(table(ids)),30)

exp1<-exp[!duplicated(exp$Symbol),]
rownames(exp1)<-exp1[,1]


exp2<-exp1[,c(1,17:31)]
exp3<-exp1[,c(17:31)]
keep_genes <- rowSums(exp3) >= 5 &rowMeans(exp3) >= 1 &rowSums(exp3)<1000000

filtered_counts <- exp3[keep_genes, ]
exp3<-filtered_counts

rownames(exp2)<-exp2[,1]
rownames(exp3)<-exp2[keep_genes,1]
#获取分组信息
if (T) {
  
  library(stringr)
  group_list<-data.frame(
    sample<-colnames(exp2)[2:16] , # 第一列（字符型）
    group<-rep(c("AGB6-MOCK","AGB6-MPXVlb","AGB6-Tecovirimat","C57-MOCK","C57-NPXVlb"),each=3)          # 第二列（数值型）
  )
  colnames(group_list)=c('sample','group')
  rownames(group_list)<-group_list$sample
  boxplot(exp3,outline=FALSE, notch=T,col=group_list, las=2)
  
}
#基因标准化
if(T){
  library(DESeq2)
  exp3<-round(exp3)
  condition<-factor(group_list$group,levels=c("AGB6-MOCK","AGB6-MPXVlb","AGB6-Tecovirimat","C57-MOCK","C57-NPXVlb"))
  colData<-data.frame(row.names=colnames(exp3),condition)
  dds<-DESeqDataSetFromMatrix(exp3,colData = colData,design=~condition)
  #将上面的矩阵构建一个dds矩阵
  dds<-dds[rowSums(counts(dds))>1,]
  
  #将所有样本基因表达量之和小于1的基因过滤掉，这一步有待深究
  dds<-DESeq(dds)
  #results(dds,c("condition","MOCK","Tecovirimat"))
  
  normalized_count<-counts(dds,normalized=TRUE)
  exp3<-normalized_count
  rld <-  assay(rlog(dds, blind=T))
  #vsd<-assay(vst(dds,blind=FALSE))
  exp3_log<-rld
  boxplot(exp3_log,outline=FALSE, notch=T,col=group_list, las=2)
  dev.off()
}

#cluster聚类
if(F){
  dir.create("cluster_plot", recursive = TRUE)

  
  library(factoextra)
  library(ggplot2)
  dev.new()
  # 设置PDF输出参数（建议使用更合理的尺寸）
  pdf("cluster_plot/Cluster_Plots.pdf", width = 7, height = 10)  # 调整为更适合期刊排版的尺寸
  
  # 使用gridExtra包创建更灵活的布局
  library(gridExtra)
  
  # --- 绘制四个子图 ---
  plots <- list()
  
  # 图1: TE_MOCK_DMSO
  cor_dist1 <- as.dist(1 - cor(exp3_log[, c(1:9)], method = "pearson"))
  hc1 <- hclust(cor_dist1, method = "ward.D2")
  p1 <- fviz_dend(hc1,
                  k = 2,
                  rect = TRUE,
                  rect_border = "gray70",
                  lower_rect = -0.01,
                  rect_lty = 2,
                  rect_fill = FALSE,
                  lwd = 0.8,
                  cex = 0.7,  # 增大标签字体
                  color_labels_by_k = TRUE,
                  horiz = TRUE,
                  palette = "jco",  # 使用更清晰的JCO配色
                  main = "AGB6_MOCK_Ib_TE") +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  #plot(as.dendrogram(hc), nodePar = list(lab.cex=0.6, pch=c(NA,19), cex=0.7, col="blue"), horiz=TRUE, main="TE_MOCK_DMSO")
  
  # 图4: MOCK_Ib_IIb (直接创建第二个图)
  cor_dist4 <- as.dist(1 - cor(exp3_log[, c(10:15)], method = "pearson"))
  hc4 <- hclust(cor_dist4, method = "ward.D2")
  p4 <- fviz_dend(hc4,
                  k = 2,
                  rect = TRUE,
                  rect_border = "gray70",
                  lower_rect = -0.011,
                  rect_lty = 2,
                  rect_fill = FALSE,
                  lwd = 0.8,
                  cex = 0.7,
                  color_labels_by_k = TRUE,
                  horiz = TRUE,
                  palette = "jco",
                  main = "C57_MOCK_Ib") +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  
  # 组合图形（2行1列布局）
  grid.arrange(p1, p4, nrow = 2, heights = c(1, 1))
  
  # 关闭图形设备
  dev.off()
  ggsave("cluster_plot/Cluster_Plots.tiff", plot = grid.arrange(p1,p4,nrow=2), 
         width = 7, height = 10, units = "in", dpi = 300)
  
}




###############################################################################################################
########################################################################################
#model.matrix将分类变量变为哑变量矩阵，每个分类变量都用一个二进制向量表示，其中1表示某个变量的某个类别，0表示不是该类别
#design分组矩阵将表达矩阵的列分为几组，通过model.matrix得到
if(T){
  library(dplyr)
  
  contrasts <- list(
    c("condition", "C57-NPXVlb", "C57-MOCK"),
    c("condition", "AGB6-MPXVlb", "AGB6-MOCK"),
    c("condition", "AGB6-Tecovirimat","AGB6-MOCK" ),
    c("condition", "AGB6-MPXVlb","AGB6-Tecovirimat" )
  )
  
  # 遍历所有对比并保存结果
  for (contrast in contrasts) {
    # 提取结果（使用lfcShrink增强log2FC稳定性）
    res <- lfcShrink(
      dds,
      contrast = contrast,
      type = "ashr"  # 或"apeglm"
    )
    
    # 添加基因列（将行名转为数据框列）
    deg <- as.data.frame(res)
    deg$gene <- rownames(deg)
    log2FoldChange_t=1
    padj_t = 0.05
    genes_to_remove <- c("")
    
    
    deg <- deg %>%
      filter(!is.na(padj)) %>%  # 先去除NA行
      filter(!gene %in% genes_to_remove) %>%
      mutate(
        change = case_when(
          padj < padj_t & log2FoldChange > log2FoldChange_t ~ "up",
          padj < padj_t & log2FoldChange < -log2FoldChange_t ~ "down",
          TRUE ~ "stable"  # 其他情况标记为stable
        )
      )
    
    # 生成文件名（格式：对比组_vs_对照组.csv）
    filename <- paste0(
      "DEG_", 
      contrast[2],  # 实验组
      "_vs_", 
      contrast[3],  # 对照组
      ".csv"
    )
    
    # 保存为CSV（保留完整结果）
    dir.create("DEG", recursive = TRUE)
    setwd("cluster_plot/")
    write.csv(
      deg[, c("gene", "baseMean", "log2FoldChange", "pvalue", "padj","change")], 
      file = filename,
      row.names = FALSE
    )
    setwd("D:\\project/hhy/mouse/second spleen/")
  }
}




######################################################################################
#1.火山图----

if(T){
  library(gridExtra)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  # 生成所有火山图的函数（新增exclude_genes参数）
  generate_volcano <- function(dat, contrast_name,filtered_P,filtered_FC , exclude_genes = c("")) {
    # 移除不需要的基因
    dat_filtered <- dat %>% 
      filter(!gene %in% exclude_genes) %>% 
      mutate(log10FC = log2FoldChange / log2(10))  # 转换为log10刻度
    
    # 筛选满足条件的标签基因
    for_label <- dat_filtered %>%
      filter(
        -log10(padj) > filtered_P |       # p值显著性阈值
          abs(log10FC) > filtered_FC          # 表达变化阈值
      )
    
    # 绘制火山图
    p <- ggplot(dat_filtered, 
                aes(x = log10FC, 
                    y = -log10(padj))) +
      # 数据点
      geom_point(aes(color = change), alpha = 0.6, size = 2.5) +
      scale_color_manual(
        values = c("down" = "#377EB8", 
                   "stable" = "grey60", 
                   "up" = "#E41A1C"),
        name = "Expression"
      ) +
      
      # --- 原始参考线（保留不变）---
      geom_vline(xintercept = c(-1/log2(10), 1/log2(10)),        # 原始log2FC阈值线
                 linetype = "solid", color = "black", linewidth = 0.5) +
      geom_hline(yintercept = -log10(0.01),
                 linetype = "solid", color = "black", linewidth = 0.5) +
      
      # --- 新增参考线（突出显示新阈值）---
      geom_vline(xintercept = c(-filtered_FC, filtered_FC),  # 使用log10刻度
                 linetype = "dashed", color = "gray50", linewidth = 0.7) +
      geom_hline(yintercept = filtered_P, 
                 linetype = "dashed", color = "gray50", linewidth = 0.7) +
      
      # 坐标轴和标题
      ggtitle(contrast_name) +
      labs(x = expression(log[10](Fold~Change)), 
           y = expression(-log[10](P~adj))) +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major = element_blank(), # 移除主要网格线
        panel.grid.minor = element_blank()  # 移除次要网格线
      )
    
    # 添加符合条件的基因标签
    if (nrow(for_label) > 0) {
      p <- p +
        geom_label_repel(
          data = for_label,
          aes(label = gene),
          size = 4,
          box.padding = 0.5,
          segment.color = "grey50",
          max.overlaps = Inf,
          min.segment.length = 0.1,
          fill = alpha("white", 0.7)
        ) +
        geom_point(
          data = for_label,
          shape = 21, 
          color = "black", 
          size = 3, 
          stroke = 1
        )
    }
    
    return(p)
  }
  
  # 使用示例
  all_contrasts <- c("AGB6-MPXVlb_vs_AGB6-MOCK","C57-NPXVlb_vs_C57-MOCK","AGB6-Tecovirimat_vs_AGB6-MOCK","AGB6-MPXVlb_vs_AGB6-Tecovirimat")
  dev.new(width = 30, height = 15)
  dir.create("volcano", recursive = TRUE, showWarnings = FALSE)
  i=1
  for (coef in all_contrasts) {
    deg <- read.csv(paste0("DEG_", coef, ".csv"), stringsAsFactors = FALSE)
    deg$gene <- as.character(deg$gene)
    print(coef)
    
    filtered_P_list = c(30, 25, 20,30)
    filtered_FC_list = c(1.2, 0.9, 1.2,1.2)
    filtered_P<-filtered_P_list[i]
    filtered_FC<-filtered_FC_list[i]
    p1 <- generate_volcano(deg, paste0(coef, "\n(Excluded Genes)"), filtered_P,filtered_FC,exclude_genes = c("Rag1","Casp14"))
    p2 <- generate_volcano(deg, paste0(coef, "\n(All Genes)"),filtered_P,filtered_FC, exclude_genes = NULL)
    i=i+1
    # 并排组合图片
    
    combined_plot <- grid.arrange(
      p1, p2,
      ncol = 2,  # 并排显示
      widths = c(1, 1)  # 两列等宽
    )
    
    if(T){
      ggsave(
        paste0("volcano/volcano_", coef, ".png"),
        plot = p2,
        width = 15,
        height = 12,
        dpi = 600,
        bg = "white"
      )
    }
    # 保存组合图
    if(F){
      ggsave(
        paste0("volcano/volcano_paired_2_", coef, ".png"),
        plot = combined_plot,
        width = 24,  # 增加宽度以适应两图并排
        height = 12,   # 保持原高度
        dpi = 600,
        bg = "white"
      )
      
      
    }
    
  }
}

#############################################################################################
#############################################################################################
if(F){
  library(dplyr)
  
  # 创建输出目录
  if (!dir.exists("DEG")) dir.create("DEG")
  
  # 定义基因筛选函数
  select_genes <- function(deg_df, group_name) {
    # 提取前50变化倍数差异基因基因
    up_genes <- deg_df %>% filter(change != "stable") %>% 
      arrange(abs(log2FoldChange)) %>% 
      slice_head(n = 50) %>% pull(gene)
    
    #write.csv(combined, paste0("DEG/Selected_", group_name, ".csv"), row.names = FALSE)
    return(select_genes)
  }
  deg_files<-c("DEG_C57-NPXVlb_vs_C57-MOCK.csv","DEG_AGB6-MPXVlb_vs_AGB6-MOCK.csv",
               "DEG_AGB6-Tecovirimat_vs_AGB6-MOCK.csv")
  
  for (files in deg_files){
    deg<-read.csv(files)
    select_genes(deg,sub(".csv$", "", files))
  }
  
  
  
}
##########################################################################

if(T){
  library(pheatmap)
  # 读取差异分析结果
  
  deg1<-read.csv("DEG_C57-NPXVlb_vs_C57-MOCK.csv")
  deg2<-read.csv("DEG_AGB6-MPXVlb_vs_AGB6-MOCK.csv")
  deg3 <- read.csv("DEG_AGB6-Tecovirimat_vs_AGB6-MOCK.csv")
  deg4<-read.csv("DEG_AGB6-MPXVlb_vs_AGB6-Tecovirimat.csv")
  library(dplyr)
  # 定义函数：提取上下调基因（返回列表包含up和down）------------------------------
  get_top_genes <- function(deg, n = 50) {
    deg %>%
      filter(!is.na(log2FoldChange)) %>%       # 去除NA值
      arrange(desc(abs(log2FoldChange))) %>%   # 按绝对值降序排列
      slice_head(n = n) %>%                    # 取前n个基因
      pull(gene)                               # 提取基因名
  }
  if(T){
    # 提取每个对比组的上下调基因 -------------------------------------------------
    cg <- get_top_genes(deg1, 50)
    exp_sub <- exp3_log[cg,c(10:15) ]
    group_select<-group_list$group[c(10:15)]
    names(group_select)<-group_list$sample[c(10:15)]
    # 创建列注释（分组信息）
    annotation_col <- data.frame(Group = group_select)
    rownames(annotation_col) <- names(group_select)
    # 绘制热图
    p1 <- pheatmap(exp_sub,
                   scale = "row",          # 按行标准化
                   show_rownames = T,   # 隐藏行名（基因名过多时）
                   show_colnames = T,
                   annotation_col = annotation_col,
                   cluster_cols = FALSE,   # 不聚类列（保持样本顺序）
                   cluster_rows = T,  
                   main = "TOP50_C57-NPXVlb_vs_C57-MOCK",  # 标题为对比名称
                   silent = TRUE  ,         # 静默模式，避免直接绘图
                   border_color = NA  
    )
    
    ggsave("heatmap_TOP50_C57-NPXVlb_vs_C57-MOCK.pdf", p1, width = 8, height = 10,limitsize = FALSE)
    
  }
  
  if(T){
    # 提取每个对比组的上下调基因 -------------------------------------------------
    cg <- get_top_genes(deg2, 50)
    exp_sub <- exp3_log[cg,c(1:6) ]
    group_select<-group_list$group[c(1:6)]
    names(group_select)<-group_list$sample[c(1:6)]
    # 创建列注释（分组信息）
    annotation_col <- data.frame(Group = group_select)
    rownames(annotation_col) <- names(group_select)
    # 绘制热图
    p1 <- pheatmap(exp_sub,
                   scale = "row",          # 按行标准化
                   show_rownames = T,   # 隐藏行名（基因名过多时）
                   show_colnames = T,
                   annotation_col = annotation_col,
                   cluster_cols = FALSE,   # 不聚类列（保持样本顺序）
                   cluster_rows = T,  
                   main = "TOP50_AGB6-MPXVlb_vs_AGB6-MOCK",  # 标题为对比名称
                   silent = TRUE  ,         # 静默模式，避免直接绘图
                   border_color = NA  
    )
    
    ggsave("heatmap_TOP50_AGB6-MPXVlb_vs_AGB6-MOCK.pdf", p1, width = 8, height = 10,limitsize = FALSE)
    
  }
  if(T){
    # 提取每个对比组的上下调基因 -------------------------------------------------
    cg <- get_top_genes(deg3, 50)
    exp_sub <- exp3_log[cg,c(1:3,7:9) ]
    group_select<-group_list$group[c(1:3,7:9)]
    names(group_select)<-group_list$sample[c(1:3,7:9)]
    # 创建列注释（分组信息）
    annotation_col <- data.frame(Group = group_select)
    rownames(annotation_col) <- names(group_select)
    # 绘制热图
    p1 <- pheatmap(exp_sub,
                   scale = "row",          # 按行标准化
                   show_rownames = T,   # 隐藏行名（基因名过多时）
                   show_colnames = T,
                   annotation_col = annotation_col,
                   cluster_cols = FALSE,   # 不聚类列（保持样本顺序）
                   cluster_rows = T,  
                   main = "TOP50_AGB6-Tecovirimat_vs_AGB6-MOCK",  # 标题为对比名称
                   silent = TRUE  ,         # 静默模式，避免直接绘图
                   border_color = NA  
    )
    
    ggsave("heatmap_TOP50_AGB6-Tecovirimat_vs_AGB6-MOCK.pdf", p1, width = 8, height = 10,limitsize = FALSE)
    
  }
  if(T){
    # 提取每个对比组的上下调基因 -------------------------------------------------
    cg <- get_top_genes(deg4, 50)
    exp_sub <- exp3_log[cg,c(4:9) ]
    group_select<-group_list$group[c(4:9)]
    names(group_select)<-group_list$sample[c(4:9)]
    # 创建列注释（分组信息）
    annotation_col <- data.frame(Group = group_select)
    rownames(annotation_col) <- names(group_select)
    # 绘制热图
    p1 <- pheatmap(exp_sub,
                   scale = "row",          # 按行标准化
                   show_rownames = T,   # 隐藏行名（基因名过多时）
                   show_colnames = T,
                   annotation_col = annotation_col,
                   cluster_cols = FALSE,   # 不聚类列（保持样本顺序）
                   cluster_rows = T,  
                   main = "TOP50_AGB6-MPXVlb_vs_AGB6-Tecovirimat",  # 标题为对比名称
                   silent = TRUE  ,         # 静默模式，避免直接绘图
                   border_color = NA  
    )
    
    ggsave("heatmap_TOP50_AGB6-MPXVlb_vs_AGB6-Tecovirimat.pdf", p1, width = 8, height = 10,limitsize = FALSE)
    
  }
}
############################################################################################################
#GO分析
library(clusterProfiler)
library(org.Mm.eg.db)  # 小鼠示例，若为人类替换为 org.Mm.eg.db
library(ggplot2)
library(purrr)
if(F){
  generate_GO <- function(deg_file) {
    deg <- read.csv(deg_file)
    
    # 筛选差异基因 (|log2FoldChange| >1 & p <0.01)
    cg <- deg$gene[deg$change != "stable"]
    rownames(exp1)<-exp1$Symbol
    gene_id<-exp1[cg,32]
    entrez_ids <- bitr(gene_id, 
                       fromType = "ENSEMBL", 
                       toType = "ENTREZID", 
                       OrgDb = org.Mm.eg.db)
    if(T){
      
      #生物过程
      ego_BP <- enrichGO(gene = entrez_ids$ENTREZID,
                         OrgDb= org.Mm.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         minGSSize = 2,
                         keyType = "ENTREZID",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.01,
                         readable =TRUE)
      kegg_result <- enrichKEGG(
        gene = entrez_ids$ENTREZID,
        organism = "mmu",  # 小鼠：mmu，人类：hsa，大鼠：rno
        keyType = "kegg",
        pvalueCutoff = 0.01,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.01,
        minGSSize = 2
      )
      if (T) {
        library(ggplot2)
        library(gridExtra)
        library(enrichplot)  # 确保加载enrichplot包
        
        top_terms <- 20
        
        # 自定义主题函数统一调整字体和条宽
        custom_theme <- function() {
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # 标题字体
            axis.text.x = element_text(size =25),  # X轴标签字体
            axis.text.y = element_text(size = 25),  # Y轴标签字体
            axis.title.x = element_text(size = 30), # X轴标题字体
            axis.title.y = element_blank()          # 隐藏Y轴标题（通常富集分析图不需要）
          )
        }
        
        # 生成BP图（直接通过enrichplot::barplot调整参数）
        bp_plot <- barplot(
          ego_BP, 
          showCategory = top_terms, 
          title = "Biological Process",
          font.size = 25,      # 富集项名称字体大小
          label_format = 50,    # 标签换行宽度（字符数）
          geom = "bar",          # 指定使用条形图
          width = 0.4            # 调整条宽（值越小条越窄）
        ) + 
          custom_theme() +
          theme(plot.margin = margin(10, 20, 10, 10))  # 调整边距
        
        # 生成KEGG图
        kegg_plot <- barplot(
          kegg_result, 
          showCategory = top_terms, 
          title = "KEGG Pathways",
          font.size = 25,
          label_format = 50,
          geom = "bar",          # 指定使用条形图
          width = 0.4            # 调整条宽（值越小条越窄）
        ) + 
          custom_theme() +
          theme(plot.margin = margin(10, 20, 10, 10))
        
        # 合并图形并调整布局
        title_plot <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = sub(".csv$", "", file), size = 6) +
          theme_void()
        
        combined_plot <- grid.arrange(
          title_plot,
          arrangeGrob(bp_plot, kegg_plot, ncol = 2),
          nrow = 2,
          heights = c(0.1, 0.9)  # 标题占10%高度
        )
        
        # 保存高清图（调整宽高比优化显示）
        ggsave(
          paste0("GO_KEGG/GO_KEGG_", substr(file, 1, nchar(file)-4), ".jpg"), 
          combined_plot, 
          width = 30,        # 增加宽度避免文字重叠
          height = 16,       # 适当降低高度
          dpi = 300
        )}
      return(combined_plot)
    }
    
  }
  deg_files<-c("DEG_MPXV_Ib_vs_MPXV_IIb.csv","DEG_MPXV_IIb_vs_MOCK.csv",
               "DEG_MPXV_Ib_vs_MOCK.csv","DEG_Tecovirimat_vs_MOCK.csv",
               "DEG_DMSO_vs_MOCK.csv" ,"DEG_Tecovirimat_vs_DMSO.csv" )
  
  dir.create("GO_KEGG")
  for (file in deg_files) {
    generate_GO(file)
    print(file)
    
  }
}
