setwd('C:/Users/13918/Desktop/IMcompre/IMcompre1/')
rm(list = ls())
gc()
####使用并行计算处理文件列表###########################
library(tidyverse)
library(data.table)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(parallel)

# 载入函数
source('R/function_uvmr.R') ## 如果要运行mrpresso多效性则用这个
source('R/function_uvmr_without_MRPRESSO.R')


# 使用示例，逐个运行！
# kegg metabolites为筛选kegg id后的代谢物，从1400代谢物中来
exposure_datasets <- c('91inflammation_cytokines')#"91inflammation_cytokines", "731immune_cells", "486metabolites","1400metabolites",'keggmetabolites'
outcome_dir <- "./outcome/" # 路径可能需要根据您的环境进行调整
get_uvmr(exposure_datasets, outcome_dir)
#

