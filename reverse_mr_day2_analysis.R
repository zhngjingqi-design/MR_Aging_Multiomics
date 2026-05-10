# ============================================================================
# 反向MR - Day2: 运行双向MR分析（使用本地下载的TSV.GZ数据）
# Bidirectional MR - Step 2: Run Reverse MR Analysis
#
# 说明：
#   将4个表观遗传年龄加速指标作为暴露，炎症因子作为结局，
#   运行反向MR分析。
#   本脚本从EBI FTP下载炎症因子的完整GWAS数据（TSV.GZ格式），
#   避免了OpenGWAS API的JWT token认证问题。
#
# 前提：
#   先运行 reverse_mr_day1_prepare_exposures.R 生成暴露数据
#
# 输出：
#   ./reverse_result/*/ 目录下包含MR分析结果
# ============================================================================

setwd('C:/Users/13918/Desktop/IMcompre/IMcompre1')

# ---- 加载包 ----
library(TwoSampleMR)
library(dplyr)
library(parallel)
library(ggplot2)
library(data.table)

# ---- 创建输出目录 ----
if (!dir.exists("reverse_result")) dir.create("reverse_result")
if (!dir.exists("reverse_outcome_tsv")) dir.create("reverse_outcome_tsv")

# ============================================================================
# Step 1: 加载暴露数据（衰老表型工具变量）
# ============================================================================
exposure_file <- "./reverse_exposure/aging_phenotypes_exposure.rda"
if (!file.exists(exposure_file)) {
  stop("Exposure file not found! Please run reverse_mr_day1_prepare_exposures.R first.")
}
load(exposure_file)

cat("Loaded aging phenotype exposures:\n")
for (nm in names(meta_f_select)) {
  cat("  -", nm, ":", nrow(meta_f_select[[nm]]), "instrument SNPs\n")
}

# 所有暴露SNP的并集（用于提取结局数据）
all_exposure_snps <- unique(unlist(lapply(meta_f_select, function(df) df$SNP)))
cat("Total unique exposure SNPs:", length(all_exposure_snps), "\n")

# ============================================================================
# Step 2: 获取炎症因子结局数据
# ============================================================================
# 从 info91 文件读取OpenGWAS ID与因子名称映射
info91 <- read.csv("./data/info91_with_abbreviations.csv", stringsAsFactors = FALSE)
cat("\nLoaded", nrow(info91), "cytokine records\n")

# ---- 选择要分析的炎症因子 ----
# 分析全部91个炎症因子作为结局
cytokines_to_analyze <- info91$TRAIT_Abbreviations
cat("Analyzing", length(cytokines_to_analyze), "cytokines:\n  ",
    paste(cytokines_to_analyze, collapse = ", "), "\n")

# ---- 构建暴露SNP的chr:pos查找表 ----
# 用于后续与EBI GWAS文件的variant_id匹配
exposure_snp_lookup <- do.call(rbind, lapply(names(meta_f_select), function(nm) {
  df <- meta_f_select[[nm]]
  data.frame(
    SNP = df$SNP,
    chr = gsub("^chr", "", as.character(df$chr.exposure)),  # 标准化: 去掉chr前缀
    pos = as.integer(df$pos.exposure),
    stringsAsFactors = FALSE
  )
}))
exposure_snp_lookup <- unique(exposure_snp_lookup)
cat("\nBuilt chr:pos lookup table for", nrow(exposure_snp_lookup), "exposure SNPs\n")

# ---- 下载并读取结局数据 ----
outcome_list <- list()

for (cyt_abbr in cytokines_to_analyze) {
  cat("\n----------------------------------------\n")
  cat("Processing cytokine:", cyt_abbr, "\n")

  # 获取该因子的信息
  idx <- which(info91$TRAIT_Abbreviations == cyt_abbr)
  accession <- info91$ACCESSION[idx]
  full_name <- info91$TRAIT[idx]

  # 检查是否已有已下载的RDA缓存
  local_rda <- paste0("./reverse_outcome_tsv/", cyt_abbr, ".rda")
  if (file.exists(local_rda)) {
    cat("  Loading cached outcome data...\n")
    load(local_rda)  # 加载 outcome 变量
    if (exists("outcome") && nrow(outcome) > 0) {
      outcome_list[[cyt_abbr]] <- outcome
      cat("  Cached data has", nrow(outcome), "SNPs\n")
      next
    }
  }

  # 下载TSV.GZ文件
  tsv_url <- sprintf(
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/%s/%s.tsv.gz",
    accession, accession
  )
  tsv_file <- paste0("./reverse_outcome_tsv/", accession, ".tsv.gz")
  min_file_size <- 100000  # 至少100KB才算成功下载

  if (!file.exists(tsv_file) || file.info(tsv_file)$size < min_file_size) {
    # 清理不完整的旧文件
    if (file.exists(tsv_file) && file.info(tsv_file)$size < min_file_size) {
      cat("  Removing incomplete file from previous attempt...\n")
      file.remove(tsv_file)
    }

    cat("  Downloading:", accession, "(~306MB)...\n")
    cat("  URL:", tsv_url, "\n")

    # 下载到临时文件，成功后重命名，防止部分下载污染缓存
    tmp_file <- paste0(tsv_file, ".tmp")
    download_ok <- FALSE

    # 尝试多种下载方法（auto和libcurl对EBI FTP常失败，优先用wininet）
    for (method in c("wininet", "libcurl", "curl", "auto")) {
      if (download_ok) break
      tryCatch({
        cat("    Trying method:", method, "...\n")
        if (method == "auto") {
          download.file(tsv_url, destfile = tmp_file, mode = "wb", quiet = TRUE)
        } else {
          download.file(tsv_url, destfile = tmp_file, mode = "wb",
                        method = method, quiet = TRUE)
        }
        if (file.exists(tmp_file) && file.info(tmp_file)$size >= min_file_size) {
          file.rename(tmp_file, tsv_file)
          download_ok <- TRUE
          cat("    Download complete! (",
              round(file.info(tsv_file)$size / 1024 / 1024, 1), "MB)\n", sep = "")
        }
      }, error = function(e) {
        cat("    Method", method, "failed:", e$message, "\n")
      })
      # 每次尝试后清理临时文件
      if (file.exists(tmp_file)) file.remove(tmp_file)
    }

    if (!download_ok) {
      cat("  Skipping", cyt_abbr, "- all download methods failed.\n")
      cat("  You can manually download and place the file at:\n")
      cat("    ", tsv_url, "\n")
      cat("  -> ", tsv_file, "\n")
      next
    }
  } else {
    cat("  Using existing file:", basename(tsv_file),
        "(", round(file.info(tsv_file)$size / 1024 / 1024, 1), "MB)\n")
  }

  # ---- 读取TSV.GZ并通过chr:pos匹配SNP ----
  cat("  Reading outcome file (306MB, this may take 1-2 min)...\n")
  tryCatch({
    # 使用data.table::fread快速读取完整TSV
    suppressMessages({
      outcome_raw <- data.table::fread(
        tsv_file, sep = "\t", header = TRUE,
        select = c("variant_id", "beta", "standard_error",
                   "effect_allele", "other_allele",
                   "p_value", "effect_allele_frequency", "n"),
        showProgress = FALSE, nThread = 4
      )
    })
    cat("  Read", nrow(outcome_raw), "rows from outcome file\n")

    # 解析variant_id提取chr和pos
    # EBI格式: variant_id = "chr_pos_ref_alt" 或 "chr:pos:ref:alt"
    # 处理两种分隔符
    if (grepl(":", outcome_raw$variant_id[1])) {
      parts <- strsplit(outcome_raw$variant_id, ":")
    } else {
      parts <- strsplit(outcome_raw$variant_id, "_")
    }
    outcome_raw$chr  <- gsub("^chr", "", as.character(sapply(parts, `[`, 1)))  # 标准化
    outcome_raw$pos  <- as.integer(sapply(parts, `[`, 2))
    outcome_raw$ref  <- sapply(parts, `[`, 3)
    outcome_raw$alt  <- sapply(parts, `[`, 4)

    # 通过chr:pos与暴露SNP匹配
    outcome_matched <- merge(
      outcome_raw,
      exposure_snp_lookup,
      by = c("chr", "pos"),
      all.x = FALSE, all.y = FALSE
    )
    cat("  Matched", nrow(outcome_matched), "SNPs by chr:pos\n")

    if (nrow(outcome_matched) == 0) {
      cat("  No matching SNPs. Skipping this cytokine.\n")
      next
    }

    # 构建TwoSampleMR格式的outcome数据框
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

    # 保存缓存（仅保存匹配的SNP以节省空间）
    save(outcome, file = local_rda)

    cat("  Prepared outcome data:", nrow(outcome), "matching SNPs\n")
    outcome_list[[cyt_abbr]] <- outcome

  }, error = function(e) {
    cat("  Error reading file:", e$message, "\n")
  })
}

cat("\n========================================\n")
cat("Successfully prepared", length(outcome_list), "cytokine outcomes for MR\n")

if (length(outcome_list) == 0) {
  stop("No outcome data available. Please check internet connection.")
}

# ============================================================================
# Step 3: 运行MR分析
# ============================================================================
cat("\n========================================\n")
cat("Step 3: Running Reverse MR Analysis\n")
cat("========================================\n")

for (exposure_name in names(meta_f_select)) {
  cat("\n----------------------------------------\n")
  cat("Exposure:", exposure_name, "\n")
  cat("----------------------------------------\n")

  exposure_dat <- meta_f_select[[exposure_name]]
  result_dir <- paste0("./reverse_result/", exposure_name, "/")
  if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

  all_mr_results <- list()
  all_heterogeneity <- list()
  all_pleiotropy <- list()

  for (outcome_name in names(outcome_list)) {
    outcome_dat <- outcome_list[[outcome_name]]

    # 取共同SNP
    common_snps <- intersect(exposure_dat$SNP, outcome_dat$SNP)
    if (length(common_snps) == 0) {
      cat("  ", outcome_name, ": 0 overlapping SNPs, skipping\n")
      next
    }

    exposure_sub <- exposure_dat[exposure_dat$SNP %in% common_snps, ]
    outcome_sub <- outcome_dat[outcome_dat$SNP %in% common_snps, ]

    tryCatch({
      harm_dat <- harmonise_data(exposure_sub, outcome_sub, action = 2)
      harm_true <- harm_dat[harm_dat$mr_keep == TRUE, ]

      if (nrow(harm_true) < 2) {
        cat("  ", outcome_name, ":", nrow(harm_true), "valid SNPs (<2), skipping\n")
        next
      }

      # MR分析（全部方法）
      mr_res <- mr(harm_true, method_list = c(
        "mr_ivw", "mr_wald_ratio",
        "mr_weighted_median",
        "mr_egger_regression",
        "mr_simple_mode",
        "mr_weighted_mode"
      ))

      # 异质性检验
      heter_res <- mr_heterogeneity(harm_true)

      # 多效性检验
      pleio_res <- tryCatch(mr_pleiotropy_test(harm_true),
                            error = function(e) NULL)

      # 计算PVE
      pve <- harm_true %>%
        group_by(id.exposure) %>%
        summarise(pve = sum((beta.exposure^2) /
                             (beta.exposure^2 +
                                ifelse(is.na(samplesize.exposure), 10000,
                                       samplesize.exposure) * se.exposure^2)))

      mr_res <- mr_res %>% left_join(pve, by = "id.exposure")

      all_mr_results[[outcome_name]] <- mr_res
      all_heterogeneity[[outcome_name]] <- heter_res
      if (!is.null(pleio_res)) all_pleiotropy[[outcome_name]] <- pleio_res

      # 输出结果概要
      ivw_row <- mr_res[mr_res$method == "Inverse variance weighted", ]
      if (nrow(ivw_row) > 0) {
        cat("  ", outcome_name, ":", nrow(harm_true), "SNPs,",
            "IVW b =", round(ivw_row$b, 4),
            "p =", format(ivw_row$pval, digits = 3, scientific = TRUE), "\n")
      } else {
        cat("  ", outcome_name, ":", nrow(harm_true), "SNPs, no IVW result\n")
      }

    }, error = function(e) {
      cat("  ", outcome_name, ": Error -", e$message, "\n")
    })
  }

  # ---- 保存该暴露的所有结果 ----
  cat("\n  Saving results for", exposure_name, "...\n")

  if (length(all_mr_results) > 0) {
    res_df <- do.call(rbind, all_mr_results)
    res_df$FDR <- p.adjust(res_df$pval, method = "fdr")
    write.csv(res_df, paste0(result_dir, "effect.csv"), row.names = FALSE)
    save(all_mr_results, file = paste0(result_dir, "mr.rda"))

    sig_df <- res_df %>%
      filter(FDR < 0.05,
             method %in% c("Wald ratio", "Inverse variance weighted"))
    write.csv(sig_df, paste0(result_dir, "fdr.csv"), row.names = FALSE)

    cat("  Total associations:", nrow(res_df),
        ", FDR significant:", nrow(sig_df), "\n")
  }

  if (length(all_heterogeneity) > 0) {
    heter_df <- do.call(rbind, all_heterogeneity)
    write.csv(heter_df, paste0(result_dir, "heterogeneity.csv"), row.names = FALSE)
    save(all_heterogeneity, file = paste0(result_dir, "heterogeneity.rda"))
  }

  if (length(all_pleiotropy) > 0) {
    pleio_df <- do.call(rbind, all_pleiotropy)
    write.csv(pleio_df, paste0(result_dir, "pleiotropy.csv"), row.names = FALSE)
    save(all_pleiotropy, file = paste0(result_dir, "pleiotropy.rda"))
  }
}

# ============================================================================
# Step 4: 汇总报告
# ============================================================================
cat("\n========================================\n")
cat("Reverse MR Analysis Complete!\n")
cat("========================================\n\n")
cat("Results in ./reverse_result/:\n")

for (exposure_name in names(meta_f_select)) {
  effect_file <- paste0("./reverse_result/", exposure_name, "/effect.csv")
  if (file.exists(effect_file)) {
    res <- read.csv(effect_file)
    sig <- res[res$FDR < 0.05 &
               res$method %in% c("Wald ratio", "Inverse variance weighted"), ]
    cat("  ", exposure_name, ":")
    if (nrow(sig) > 0) {
      cat(" FDR significant ->", paste(unique(sig$outcome), collapse = ", "), "\n")
    } else {
      cat(" No FDR significant reverse associations\n")
    }
  }
}

cat("\nDone!\n")
