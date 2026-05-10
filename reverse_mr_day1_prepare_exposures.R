# ============================================================================
# 反向MR - Day1: 准备暴露数据（表观遗传年龄加速指标 → 从VCF提取工具变量）
# Bidirectional MR - Step 1: Prepare Exposure Data
#
# 说明：
#   将4个表观遗传年龄加速指标的VCF文件转换为暴露格式，
#   提取显著SNP (P < 5e-8)，进行LD clumping，
#   保存为与正向MR兼容的 meta_f_select 列表格式。
#
# 输出：
#   ./reverse_exposure/aging_phenotypes_exposure.rda
#     包含 meta_f_select (命名列表，每个元素是一个衰老表型的工具变量数据框)
#
# 使用前：
#   1. 确保已安装 gwasglue 包
#   2. 确保VCF文件存在于指定路径
# ============================================================================

setwd('C:/Users/13918/Desktop/IMcompre/IMcompre1')

# ---- 加载包 ----
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)
library(dplyr)

# ---- 加载VCF转换函数 ----
source('./R/function_get_vcf_local.R')

# ---- 创建输出目录 ----
if (!dir.exists("reverse_exposure")) {
  dir.create("reverse_exposure")
}

# ============================================================================
# 辅助函数: 本地LD clumping（避免OpenGWAS API的token需求）
# ============================================================================
run_local_clumping <- function(dat, clump_kb = 10000, clump_r2 = 0.001,
                                clump_p = 5e-8, pop = "EUR") {
  # 尝试多种方法进行LD clumping

  # ---- 方法1: 使用 ieugwasr 本地 plink ----
  if (requireNamespace("ieugwasr", quietly = TRUE)) {
    plink_exe <- tryCatch({
      ieugwasr::get_plink_exe()
    }, error = function(e) NULL)

    if (!is.null(plink_exe)) {
      ld_ref_dir <- "./ldref/"
      if (!dir.exists(ld_ref_dir)) dir.create(ld_ref_dir, showWarnings = FALSE)
      bfile_path <- paste0(ld_ref_dir, pop)

      if (!file.exists(paste0(bfile_path, ".bed"))) {
        cat("    Downloading 1000G", pop, "LD reference (approx. 600MB)...\n")
        cat("    This is a one-time download.\n")
        base_url <- "https://api.bioinfo-icts.ut.ee/cazy/1000G_phase3/"
        for (ext in c(".bed", ".bim", ".fam")) {
          tryCatch({
            download.file(paste0(base_url, pop, ext),
                          destfile = paste0(bfile_path, ext),
                          mode = "wb", quiet = TRUE)
          }, error = function(e) {})
        }
        if (!file.exists(paste0(bfile_path, ".bed"))) {
          cat("    Trying alternative LD reference source...\n")
          base_url2 <- "https://github.com/explodecomputer/plink-files/raw/master/"
          for (ext in c(".bed", ".bim", ".fam")) {
            tryCatch({
              download.file(paste0(base_url2, pop, ext),
                            destfile = paste0(bfile_path, ext),
                            mode = "wb", quiet = TRUE)
            }, error = function(e) {})
          }
        }
      }

      if (file.exists(paste0(bfile_path, ".bed"))) {
        cat("    Using local LD clumping with plink...\n")
        result <- tryCatch({
          ieugwasr::ld_clump(dat = dat, clump_kb = clump_kb,
                             clump_r2 = clump_r2, clump_p = clump_p,
                             pop = pop, plink_bin = plink_exe,
                             bfile = bfile_path)
        }, error = function(e) {
          cat("    Local clumping failed:", e$message, "\n")
          NULL
        })
        if (!is.null(result)) return(result)
      } else {
        cat("    LD reference download failed.\n")
      }
    }
  }

  # ---- 方法2: 基于距离的简单clumping ----
  cat("    Using distance-based clumping (genomic position)...\n")
  if (!"chr.exposure" %in% names(dat) || !"pos.exposure" %in% names(dat)) {
    cat("    No genomic position info. Returning all SNPs.\n")
    return(dat)
  }
  dat <- dat[order(dat$pval.exposure), ]
  keep_snps <- c()
  for (chr in unique(dat$chr.exposure)) {
    chr_dat <- dat[dat$chr.exposure == chr, ]
    chr_dat <- chr_dat[order(chr_dat$pos.exposure), ]
    while (nrow(chr_dat) > 0) {
      keep_snps <- c(keep_snps, chr_dat$SNP[1])
      pos <- chr_dat$pos.exposure[1]
      chr_dat <- chr_dat[chr_dat$pos.exposure < (pos - clump_kb) |
                           chr_dat$pos.exposure > (pos + clump_kb), ]
    }
  }
  result <- dat[dat$SNP %in% keep_snps, ]
  return(result)
}

# ============================================================================
# Step 1: 定义4个表观遗传年龄加速指标的VCF文件路径
# ============================================================================
# 根据您的表格：
#   HannumAA → ebi-a-GCST90014289
#   IEAA     → ebi-a-GCST90014290
#   PhenoAA  → ebi-a-GCST90014292
#   GrimAA   → ebi-a-GCST90014300
# （可选的: Frailty → ebi-a-GCST90020053, Telomere → ieu-b-4879）

aging_phenotypes <- list(
  HannumAA = list(
    vcf_path = "./OUTCOMErdata/ebi-a-GCST90014289.vcf.gz",
    trait_name = "HannumAA"
  ),
  IEAA = list(
    vcf_path = "./OUTCOMErdata/ebi-a-GCST90014290.vcf.gz",
    trait_name = "IEAA"
  ),
  PhenoAA = list(
    vcf_path = "./OUTCOMErdata/ebi-a-GCST90014292.vcf.gz",
    trait_name = "PhenoAA"
  ),
  GrimAA = list(
    vcf_path = "./OUTCOMErdata/ebi-a-GCST90014300.vcf.gz",
    trait_name = "GrimAA"
  ),
  Frailty = list(
    vcf_path = "./OUTCOMErdata/ebi-a-GCST90020053.vcf.gz",
    trait_name = "Frailty"
  ),
  Telomere = list(
    vcf_path = "./outcome/ieu-b-4879.vcf.gz",
    trait_name = "Telomere_length"
  )
)

# ============================================================================
# Step 2: 逐个转换VCF为暴露格式，提取工具变量
# ============================================================================

meta_f_select <- list()  # 存储所有衰老表型的工具变量

for (name in names(aging_phenotypes)) {
  cat("\n========================================\n")
  cat("Processing:", name, "\n")
  cat("========================================\n")

  info <- aging_phenotypes[[name]]

  # ---- 2a: 检查VCF文件是否存在 ----
  if (!file.exists(info$vcf_path)) {
    warning(paste("VCF file not found:", info$vcf_path))
    next
  }

  # ---- 2b: 读取VCF并转换为暴露格式 ----
  cat("  Reading VCF...\n")
  vcf <- VariantAnnotation::readVcf(info$vcf_path)

  cat("  Converting to exposure format...\n")
  exposure_raw <- gwasglue::gwasvcf_to_TwoSampleMR(vcf = vcf, type = "exposure")

  cat("  Total SNPs in VCF:", nrow(exposure_raw), "\n")

  # ---- 2c: 设置暴露名称 ----
  exposure_raw$exposure <- info$trait_name
  exposure_raw$id.exposure <- info$trait_name

  # ---- 2d: 筛选显著SNP (P < 5e-8) ----
  # GWAS标准显著性阈值
  exposure_sig <- subset(exposure_raw, pval.exposure < 5e-8)
  cat("  SNPs with P < 5e-8:", nrow(exposure_sig), "\n")

  # 如果P < 5e-8没有SNP，放宽到5e-7
  if (nrow(exposure_sig) == 0) {
    cat("  No SNP at P < 5e-8. Relaxing threshold to P < 5e-7...\n")
    exposure_sig <- subset(exposure_raw, pval.exposure < 5e-7)
    cat("  SNPs with P < 5e-7:", nrow(exposure_sig), "\n")
  }

  # 如果仍然没有，再放宽到5e-6
  if (nrow(exposure_sig) == 0) {
    cat("  No SNP at P < 5e-7. Relaxing threshold to P < 5e-6...\n")
    exposure_sig <- subset(exposure_raw, pval.exposure < 5e-6)
    cat("  SNPs with P < 5e-6:", nrow(exposure_sig), "\n")
  }

  # 如果还是没有，跳过该表型
  if (nrow(exposure_sig) == 0) {
    warning(paste("No significant SNPs found for", name, "at P < 5e-6. Skipping."))
    next
  }

  # ---- 2e: LD Clumping（使用本地方法，避免API token问题）----
  cat("  Running LD clumping...\n")
  if (nrow(exposure_sig) >= 2) {
    exposure_clumped <- run_local_clumping(
      exposure_sig,
      clump_kb = 10000,
      clump_r2 = 0.001,
      clump_p = 5e-8,
      pop = "EUR"
    )
  } else {
    exposure_clumped <- exposure_sig
  }

  cat("  SNPs after clumping:", nrow(exposure_clumped), "\n")

  # ---- 2f: 计算F统计量 ----
  exposure_clumped <- exposure_clumped %>%
    mutate(
      # F = (beta^2) / (beta^2 + N * se^2) 的变体
      # 标准F统计量: F = (beta/se)^2
      f_statistic = (beta.exposure / se.exposure)^2
    )

  cat("  F-statistics: min =", round(min(exposure_clumped$f_statistic), 2),
      ", mean =", round(mean(exposure_clumped$f_statistic), 2),
      ", max =", round(max(exposure_clumped$f_statistic), 2), "\n")

  # 剔除弱工具变量 (F < 10)
  n_before <- nrow(exposure_clumped)
  exposure_clumped <- subset(exposure_clumped, f_statistic > 10)
  n_removed <- n_before - nrow(exposure_clumped)
  if (n_removed > 0) {
    cat("  Removed", n_removed, "weak instruments (F < 10)\n")
  }

  cat("  Final instruments:", nrow(exposure_clumped), "SNPs\n")

  # ---- 2g: 存储到列表 ----
  meta_f_select[[name]] <- exposure_clumped
}

# ============================================================================
# Step 3: 保存暴露数据
# ============================================================================

cat("\n========================================\n")
cat("Saving exposure data...\n")
cat("Total phenotypes with valid instruments:", length(meta_f_select), "\n")
for (name in names(meta_f_select)) {
  cat("  -", name, ":", nrow(meta_f_select[[name]]), "SNPs\n")
}

save(meta_f_select, file = "./reverse_exposure/aging_phenotypes_exposure.rda")
cat("\nExposure data saved to: ./reverse_exposure/aging_phenotypes_exposure.rda\n")
cat("Done!\n")
