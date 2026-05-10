# 设置路径一定要正确
setwd('C:\\Users\\13918\\Desktop\\IMcompre\\IMcompre1')
library(dplyr)
library(data.table)
#install.packages("MendelianRandomization")
#安装R包
library(MendelianRandomization)

# 载入函数

source('R\\function_pheno_out.R')
# 调用函数,根据自己表型输入,不区分大小写

outcome_keywords <- c("HannumAA")
###免疫细胞
get_phe_out("731immune_cells", outcome_keywords)
###炎症因子
get_phe_out("91inflammation_cytokines", outcome_keywords)
###代谢物
###486代谢物
get_phe_out("486metabolites", outcome_keywords)
###1400代谢物
get_phe_out("1400metabolites", outcome_keywords)
###kegg代谢物，从1400代谢物中来
get_phe_out("keggmetabolites", outcome_keywords)



ebi-a-GCST90014289

