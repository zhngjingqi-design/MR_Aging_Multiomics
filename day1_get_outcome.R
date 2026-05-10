setwd("C:/Users/13918/Desktop/IMcompre/IMcompre1/") 

# 安装 Bioconductor 管理的包（包括 VariantAnnotation）

# 安装 CRAN 上的 gwasglue 包

# 安装完成后加载包
library(VariantAnnotation)
library(gwasglue)
.libPaths()
source("gwasglue-master/相关文件.R")
library(devtools)
source('./R/function_get_vcf_local.R')

### if you want to read a single vcf file:
#get_local("eqtl-a-ENSG00000119917.vcf.gz", type="outcome", outfile="./outcome/IFIT2.rda")

### multiple，支持多个转换，但推荐一个个来
# 假设vcf_files是你需要处理的文件名列表
vcf_files <- list.files(path = '.',pattern = 'vcf.gz')
#vcf_files <- vcf_files[2]
# 对每个VCF文件进行操作
for (vcf_file in vcf_files) {
  # 生成输出文件名（使用gsub函数替换文件扩展名）
  outfile <- paste0("./outcome/",gsub("\\.vcf\\.gz$", ".rda", vcf_file))
  
  # 运行函数
  get_local(vcf_file, type="outcome", outfile=outfile)
}
 

