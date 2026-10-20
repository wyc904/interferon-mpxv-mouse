# 加载必要的包
library(dplyr)
library(tidyr)

# 创建数据框（使用您提供的数据）
sample_info <- data.frame(
  sample = c(
    "JZ25070166.AGB6.MOCK.1.P.AGB6.MOCK.1.P_combined",
    "JZ25070167.AGB6.MOCK.2.P.AGB6.MOCK.2.P_combined", 
    "JZ25070168.AGB6.MOCK.3.P.AGB6.MOCK.3.P_combined",
    "JZ25070169.AGB6.MPXV_lb.1.P.AGB6.MPXV_lb.1.P_combined",
    "JZ25070170.AGB6.MPXV_lb.2.P.AGB6.MPXV_lb.2.P_combined",
    "JZ25070171.AGB6.MPXV_lb.3.P.AGB6.MPXV_lb.3.P_combined",
    "JZ25070172.AGB6.Tecovirimat.1.P.AGB6.Tecovirimat.1.P_combined",
    "JZ25070173.AGB6.Tecovirimat.2.P.AGB6.Tecovirimat.2.P_combined",
    "JZ25070174.AGB6.Tecovirimat.3.P.AGB6.Tecovirimat.3.P_combined",
    "JZ25070175.C57.MOCK.1.P.C57.MOCK.1.P_combined",
    "JZ25070176.C57.MOCK.2.P.C57.MOCK.2.P_combined",
    "JZ25070177.C57.MOCK.3.P.C57.MOCK.3.P_combined",
    "JZ25070178.C57.MPXV_lb.1.P.C57.MPXV_lb.1.P_combined",
    "JZ25070179.C57.MPXV_lb.2.P.C57.MPXV_lb.2.P_combined"
  ),
  group = c(
    "combined_NA", "combined_NA", "combined_NA",
    "lb.1.P.AGB6.MPXV_lb.1.P", "lb.2.P.AGB6.MPXV_lb.2.P", "lb.3.P.AGB6.MPXV_lb.3.P",
    "combined_NA", "combined_NA", "combined_NA",
    "combined_NA", "combined_NA", "combined_NA", 
    "lb.1.P.C57.MPXV_lb.1.P", "lb.2.P.C57.MPXV_lb.2.P"
  ),
  replicate = c(NA, NA, NA, "combined", "combined", "combined", 
                NA, NA, NA, NA, NA, NA, "combined", "combined")
)

# 方法1：基于样本名称的模式匹配进行分组
sample_info_cleaned <- sample_info %>%
  mutate(
    # 从样本名称中提取关键信息
    strain = ifelse(grepl("AGB6", sample), "AGB6", 
                    ifelse(grepl("C57", sample), "C57", NA)),
    treatment = case_when(
      grepl("MOCK", sample) ~ "MOCK",
      grepl("MPXV_lb", sample) ~ "MPXV_lb", 
      grepl("Tecovirimat", sample) ~ "Tecovirimat",
      TRUE ~ NA_character_
    ),
    replicate_num = as.numeric(gsub(".*?(\\d)\\.P.*", "\\1", sample))
  ) %>%
  # 创建清晰的分组标签
  mutate(
    group_clear = paste(strain, treatment, sep = "."),
    sample_id = paste(strain, treatment, replicate_num, sep = ".")
  )


print(sample_info_cleaned %>% select(sample, strain, treatment, replicate_num, group_clear))

# 方法2：筛选出您需要的两个特定组（AGB6.MOCK 和 AGB6.MPXV_lb）
target_groups <- sample_info_cleaned %>%
  filter(group_clear %in% c("AGB6.MOCK", "AGB6.MPXV_lb"))

print("目标分析组（AGB6.MOCK 和 AGB6.MPXV_lb）：")
print(target_groups)

# 方法3：创建分组汇总统计
group_summary <- target_groups %>%
  group_by(group_clear) %>%
  summarise(
    n_replicates = n(),
    replicate_numbers = paste(sort(unique(replicate_num)), collapse = ", "),
    sample_ids = paste(sample_id, collapse = "; ")
  )

print("分组汇总信息：")
print(group_summary)

# 方法4：为后续分析创建样本列表
# 创建AGB6.MOCK组的样本列表
agb6_mock_samples <- target_groups %>%
  filter(group_clear == "AGB6.MOCK") %>%
  pull(sample)

# 创建AGB6.MPXV_lb组的样本列表  
agb6_mpxv_lb_samples <- target_groups %>%
  filter(group_clear == "AGB6.MPXV_lb") %>%
  pull(sample)

print("AGB6.MOCK组样本：")
print(agb6_mock_samples)

print("AGB6.MPXV_lb组样本：")
print(agb6_mpxv_lb_samples)

# 方法5：创建用于统计分析的数据框
analysis_design <- data.frame(
  sample = c(agb6_mock_samples, agb6_mpxv_lb_samples),
  group = c(rep("AGB6.MOCK", length(agb6_mock_samples)),
            rep("AGB6.MPXV_lb", length(agb6_mpxv_lb_samples))),
  replicate = c(1:3, 1:3)  # 每个组都有3个重复
)

print("最终分析设计：")
print(analysis_design)


# 保存分组信息到文件
write.csv(analysis_design, "sample_grouping_analysis.csv", row.names = FALSE)
write.csv(target_groups, "detailed_sample_info.csv", row.names = FALSE)

