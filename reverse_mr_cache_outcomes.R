# ============================================================================
# 反向MR - 预处理：将TSV.GZ缓存为小体积.rda文件
# Reverse MR - Preprocess: Cache TSV.GZ to compact .rda files
#
# 用法: source("reverse_mr_cache_outcomes.R")
# 说明:
#   - 前提：先运行 reverse_mr_day1_prepare_exposures.R 生成暴露数据
#           + reverse_mr_download_outcomes.R 下载91个TSV.GZ文件
#   - 读入每个TSV.GZ（~306MB），通过chr:pos匹配暴露SNP
#   - 保存为小体积.rda缓存文件（仅含匹配的几十个SNP，几KB）
#   - 全部缓存完成后，reverse_mr_day2_analysis.R 可直接秒加载
#
# 输出: ./reverse_outcome_tsv/{cytokine}.rda （小体积）
# ============================================================================

setwd('C:/Users/13918/Desktop/IMcompre/IMcompre1')

library(data.table)
library(dplyr)

# ---- 加载暴露数据（获取chr:pos查找表）----
exposure_file <- "./reverse_exposure/aging_phenotypes_exposure.rda"
if (!file.exists(exposure_file)) {
  stop("Exposure file not found! Please run reverse_mr_day1_prepare_exposures.R first.")
}
load(exposure_file)

cat("Loaded exposure data:\n")
for (nm in names(meta_f_select)) {
  cat("  -", nm, ":", nrow(meta_f_select[[nm]]), "instrument SNPs\n")
}

# 构建chr:pos查找表（来自所有暴露表型的工具变量）
exposure_snp_lookup <- do.call(rbind, lapply(names(meta_f_select), function(nm) {
  df <- meta_f_select[[nm]]
  data.frame(
    SNP = df$SNP,
    chr = gsub("^chr", "", as.character(df$chr.exposure)),
    pos = as.integer(df$pos.exposure),
    stringsAsFactors = FALSE
  )
}))
exposure_snp_lookup <- unique(exposure_snp_lookup)
cat("\nUnique exposure SNPs for matching:", nrow(exposure_snp_lookup), "\n")
cat("  -> Will match by (chr, pos) coordinates\n\n")

# ---- 加载炎症因子信息表 ----
info91 <- read.csv("./data/info91_with_abbreviations.csv", stringsAsFactors = FALSE)
cat("Cytokine records:", nrow(info91), "\n\n")

# ---- 逐文件缓存 ----
tsv_dir <- "./reverse_outcome_tsv"
existing_tsv <- list.files(tsv_dir, pattern = "\\.tsv\\.gz$")
cat("Found", length(existing_tsv), "TSV.GZ files in", tsv_dir, "\n\n")

cached <- 0
skipped <- 0
failed_cyt <- c()

for (i in 1:nrow(info91)) {
  accession <- info91$ACCESSION[i]
  cyt_abbr  <- info91$TRAIT_Abbreviations[i]

  rda_file <- paste0(tsv_dir, "/", cyt_abbr, ".rda")
  tsv_file <- paste0(tsv_dir, "/", accession, ".tsv.gz")

  # 跳过已缓存的
  if (file.exists(rda_file)) {
    skipped <- skipped + 1
    next
  }

  # 检查TSV.GZ是否存在
  if (!file.exists(tsv_file)) {
    cat(sprintf("[%3d/%3d] %-12s - TSV.GZ not found, skipping\n", i, nrow(info91), cyt_abbr))
    failed_cyt <- c(failed_cyt, cyt_abbr)
    next
  }

  cat(sprintf("[%3d/%3d] %-12s reading...", i, nrow(info91), cyt_abbr))

  tryCatch({
    # 用fread快速读入所需列
    outcome_raw <- data.table::fread(
      tsv_file, sep = "\t", header = TRUE,
      select = c("variant_id", "beta", "standard_error",
                 "effect_allele", "other_allele",
                 "p_value", "effect_allele_frequency", "n"),
      showProgress = FALSE, nThread = 4
    )

    # 解析variant_id提取chr和pos
    if (grepl(":", outcome_raw$variant_id[1])) {
      parts <- strsplit(outcome_raw$variant_id, ":")
    } else {
      parts <- strsplit(outcome_raw$variant_id, "_")
    }
    outcome_raw$chr <- gsub("^chr", "", as.character(sapply(parts, `[`, 1)))
    outcome_raw$pos <- as.integer(sapply(parts, `[`, 2))

    # chr:pos匹配
    outcome_matched <- merge(
      outcome_raw, exposure_snp_lookup,
      by = c("chr", "pos"),
      all.x = FALSE, all.y = FALSE
    )

    if (nrow(outcome_matched) == 0) {
      cat(" 0 matched, skipping\n")
      failed_cyt <- c(failed_cyt, cyt_abbr)
      next
    }

    # 构建TwoSampleMR格式
    outcome <- data.frame(
      SNP = outcome_matched$SNP,
      chr.outcome = outcome_matched$chr,
      pos.outcome = outcome_matched$pos,
      beta.outcome = outcome_matched$beta,
      se.outcome = outcome_matched$standard_error,
      effect_allele.outcome = outcome_matched$effect_allele,
      other_allele.outcome = outcome_matched$other_allele,
      pval.outcome = outcome_matched$p_value,
      eaf.outcome = outcome_matched$effect_allele_frequency,
      samplesize.outcome = outcome_matched$n,
      outcome = cyt_abbr,
      id.outcome = cyt_abbr,
      stringsAsFactors = FALSE
    )

    save(outcome, file = rda_file)
    cached <- cached + 1
    cat(sprintf(" %d matched -> %s (%.1fKB)\n",
                nrow(outcome), basename(rda_file),
                file.info(rda_file)$size / 1024))

  }, error = function(e) {
    cat(" ERROR:", conditionMessage(e), "\n")
    failed_cyt <- c(failed_cyt, cyt_abbr)
  })
}

# ---- 汇总 ----
cat("\n========================================\n")
cat("Cache Complete!\n")
cat("========================================\n")
cat("  Total cytokines: ", nrow(info91), "\n")
cat("  Newly cached:    ", cached, "\n")
cat("  Previously done: ", skipped, "\n")
cat("  Failed/empty:    ", length(failed_cyt), "\n")
cat("  Ready for MR:    ", skipped + cached, "/", nrow(info91), "\n")

if (length(failed_cyt) > 0) {
  cat("\n  Failed cytokines:\n")
  for (f in failed_cyt) {
    acc <- info91$ACCESSION[info91$TRAIT_Abbreviations == f]
    cat("    ", f, "(", acc, ")\n")
  }
}

if (skipped + cached == nrow(info91)) {
  cat("\n  All cytokines cached! Now run reverse_mr_day2_analysis.R\n")
}

cat("\nDone!\n")
