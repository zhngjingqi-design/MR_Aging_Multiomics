# ============================================================================
# 反向MR - 下载所有91个炎症因子结局GWAS数据
# Reverse MR - Download All 91 Cytokine Outcome GWAS Files
#
# 用法: source("reverse_mr_download_outcomes.R")
# 说明:
#   - 从EBI FTP下载91个炎症因子的完整GWAS TSV.GZ文件（每个约306MB）
#   - 总计下载 ~28GB，建议在稳定网络下运行
#   - 已下载完成的文件会自动跳过（断点续传）
#   - 下载失败的文件会在最后汇总报告
#
# 输出: ./reverse_outcome_tsv/GCST*.tsv.gz
# ============================================================================

setwd('C:/Users/13918/Desktop/IMcompre/IMcompre1')

# ---- 加载信息表 ----
info91 <- read.csv("./data/info91_with_abbreviations.csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(info91), "cytokine records\n\n")

# ---- 创建下载目录 ----
if (!dir.exists("reverse_outcome_tsv")) dir.create("reverse_outcome_tsv")

# ---- 下载参数 ----
min_file_size <- 50 * 1024 * 1024  # TSV.GZ至少50MB才算成功（完整约306MB）
total_to_download <- sum(!sapply(info91$ACCESSION, function(acc) {
  f <- paste0("./reverse_outcome_tsv/", acc, ".tsv.gz")
  file.exists(f) && file.info(f)$size >= min_file_size
}))
cat("Files to download:", total_to_download, "/", nrow(info91), "\n")
cat("Estimated total size: ~", round(total_to_download * 0.3, 1), "GB\n\n")

if (total_to_download == 0) {
  cat("All files already downloaded. Nothing to do.\n")
} else {
  cat("Starting downloads...\n")
  cat("========================================\n")
}

# ---- 下载循环 ----
downloaded <- 0
skipped   <- 0
failed    <- 0
failed_list <- c()

for (i in 1:nrow(info91)) {
  accession <- info91$ACCESSION[i]
  cyt_abbr  <- info91$TRAIT_Abbreviations[i]
  full_name <- info91$TRAIT[i]

  tsv_file <- paste0("./reverse_outcome_tsv/", accession, ".tsv.gz")

  # 检查是否已完成
  if (file.exists(tsv_file) && file.info(tsv_file)$size >= min_file_size) {
    skipped <- skipped + 1
    next
  }

  # 清理不完整的旧文件
  if (file.exists(tsv_file)) file.remove(tsv_file)

  tsv_url <- sprintf(
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/%s/%s.tsv.gz",
    accession, accession
  )

  cat(sprintf("[%3d/%3d] %-12s %s", i, nrow(info91), cyt_abbr, full_name), "\n")
  cat("       ", accession, "\n")

  # 下载：wininet优先（EBI FTP上auto/libcurl常失败）
  tmp_file <- paste0(tsv_file, ".tmp")
  download_ok <- FALSE

  for (method in c("wininet", "libcurl", "curl", "auto")) {
    if (download_ok) break
    tryCatch({
      if (method == "auto") {
        download.file(tsv_url, destfile = tmp_file, mode = "wb",
                      quiet = FALSE, timeout = 300)
      } else {
        download.file(tsv_url, destfile = tmp_file, mode = "wb",
                      method = method, quiet = FALSE, timeout = 300)
      }
      if (file.exists(tmp_file) && file.info(tmp_file)$size >= min_file_size) {
        file.rename(tmp_file, tsv_file)
        download_ok <- TRUE
        downloaded <- downloaded + 1
        fsize_mb <- round(file.info(tsv_file)$size / 1024 / 1024, 1)
        cat("       OK (", method, ", ", fsize_mb, "MB) [",
            downloaded, " downloaded]\n\n", sep = "")
      } else {
        cat("       File too small (",
            round(file.info(tmp_file)$size / 1024, 1), "KB), retrying...\n", sep = "")
      }
    }, error = function(e) {
      cat("       Method", method, "failed:", conditionMessage(e), "\n")
    })
    if (file.exists(tmp_file)) file.remove(tmp_file)
  }

  if (!download_ok) {
    failed <- failed + 1
    failed_list <- c(failed_list, cyt_abbr)
    cat("       FAILED after all methods. Will report at end.\n\n")

    # 给出手动下载提示
    cat("       Manual download:\n")
    cat("        URL: ", tsv_url, "\n")
    cat("        Save to: ", tsv_file, "\n\n")
  }
}

# ---- 汇总报告 ----
cat("========================================\n")
cat("Download Summary\n")
cat("========================================\n")
cat("  Total cytokines:    ", nrow(info91), "\n")
cat("  Already downloaded: ", skipped, "\n")
cat("  Newly downloaded:   ", downloaded, "\n")
cat("  Failed:             ", failed, "\n")
cat("  Files ready:        ", skipped + downloaded, "/", nrow(info91), "\n")

if (length(failed_list) > 0) {
  cat("\n  Failed cytokines (", length(failed_list), "):\n", sep = "")
  for (f in failed_list) {
    acc <- info91$ACCESSION[info91$TRAIT_Abbreviations == f]
    cat("    ", f, " (", acc, ")\n", sep = "")
    cat("    URL: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/",
        acc, "/", acc, ".tsv.gz\n", sep = "")
  }
  cat("\n  Re-run this script to retry failed downloads.\n")
}

if (skipped + downloaded == nrow(info91)) {
  cat("\n  All 91 files ready! You can now run reverse_mr_day2_analysis.R\n")
}

cat("\nDone!\n")
