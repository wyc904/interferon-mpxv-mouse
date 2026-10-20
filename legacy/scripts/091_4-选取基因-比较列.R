






library(openxlsx)
# 最简单的版本：筛选B列同时大于C、D、E、F列的数据
library(readxl)

# 读取数据
df <- read_excel("D:/Desktop/时间点平均表达.xlsx")

# # 筛选条件：B列 > C列 且 B列 > D列 且 B列 > E列 且 B列 > F列
# result <- df[df[[5]] > df[[2]] &
#                df[[5]] > df[[3]] &
#                df[[5]] > df[[4]] &
#                df[[5]] > df[[6]], ]


# result <- df[df[[3]] > df[[2]] &
#                df[[3]] > df[[5]] &
#                df[[3]] > df[[6]], ]
# 

# 或者使用10%浮动
result <- df[df[[5]] >= 0.9 * df[[2]] &
                         df[[5]] >= 0.9 * df[[3]] &
                         df[[5]] >= 0.9 * df[[4]] &
                         df[[5]] >= 0.9 * df[[6]], ]

# 输出结果
print(result)








write.xlsx(result, "7天.xlsx")



