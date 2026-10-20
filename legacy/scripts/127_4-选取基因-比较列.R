library(readxl)
library(openxlsx)

# 读取数据
df <- read_excel("D:/Desktop/绘制3体重热图.xlsx")

# 筛选条件：第2、3、4行均大于第5、6、7行
# 即：第2列 > 第5列 且 第3列 > 第6列 且 第4列 > 第7列
# result <- df[df[[2]] > df[[5]] &
#                df[[3]] > df[[6]] &
#                df[[4]] > df[[7]], ]


result <- df[df[[5]] > df[[2]] &
               df[[6]] > df[[3]] &
               df[[7]] > df[[4]], ]






# 输出结果
write.xlsx(result, "绘制3体重热图筛选结果上调.xlsx")














