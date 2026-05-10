# ============================================================================
# 反向MR - 绘图：FDR显著结果 + 敏感性指标（对齐修复版）
# Usage: source("reverse_mr_plot_results.R")
# Output: ./reverse_plots/ 下6个PDF + 1个汇总火山图
# ============================================================================

setwd('C:/Users/13918/Desktop/IMcompre/IMcompre1')

for (pkg in c("dplyr","ggplot2","gridExtra","cowplot")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}

if (!dir.exists("reverse_plots")) dir.create("reverse_plots")

exposures <- c("HannumAA","IEAA","PhenoAA","GrimAA","Frailty","Telomere")
exposure_labels <- c(HannumAA="HannumAA", IEAA="IEAA", PhenoAA="PhenoAA",
                     GrimAA="GrimAA", Frailty="Frailty Index", Telomere="Telomere Length")

total_plots <- 0

for (exp_name in exposures) {
  cat("\n=== Processing:", exp_name, "===\n")
  rd <- paste0("./reverse_result/", exp_name, "/")

  if (!file.exists(paste0(rd,"effect.csv"))) { cat("  No data\n"); next }
  if (!file.exists(paste0(rd,"heterogeneity.csv"))) { cat("  No hetero\n"); next }
  if (!file.exists(paste0(rd,"pleiotropy.csv"))) { cat("  No pleio\n"); next }

  mr_df    <- read.csv(paste0(rd,"effect.csv"), stringsAsFactors=FALSE)
  heter_df <- read.csv(paste0(rd,"heterogeneity.csv"), stringsAsFactors=FALSE)
  pleio_df <- read.csv(paste0(rd,"pleiotropy.csv"), stringsAsFactors=FALSE)

  ivw <- mr_df[mr_df$method == "Inverse variance weighted", ]
  ivw_sig <- ivw[ivw$FDR < 0.05 & !is.na(ivw$FDR), ]
  if (nrow(ivw_sig) == 0) { cat("  No FDR significant\n"); next }

  ivw_sig$ci_lo <- ivw_sig$b - 1.96 * ivw_sig$se
  ivw_sig$ci_up <- ivw_sig$b + 1.96 * ivw_sig$se
  ivw_sig$neg_log10_FDR <- -log10(ivw_sig$FDR)

  het <- heter_df[heter_df$method == "Inverse variance weighted", c("outcome","Q","Q_pval")]
  ivw_sig <- merge(ivw_sig, het, by="outcome", all.x=TRUE)
  ple <- pleio_df[, c("outcome","egger_intercept","pval")]
  names(ple)[3] <- "pleiotropy_pval"
  ivw_sig <- merge(ivw_sig, ple, by="outcome", all.x=TRUE)

  ivw_sig <- ivw_sig[order(ivw_sig$b), ]
  ivw_sig$outcome <- factor(ivw_sig$outcome, levels=ivw_sig$outcome)
  cat("  Plotting", nrow(ivw_sig), "cytokines\n")

  # ---- 通用theme（保证三面板y轴对齐）----
  base_theme <- theme_bw(base_size=9) +
    theme(panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
          axis.text.y=element_text(size=8), plot.margin=margin(3,3,3,3))

  # ---- Panel 1: Dot plot ----
  x_lim <- max(abs(ivw_sig$ci_up), abs(ivw_sig$ci_lo), na.rm=TRUE) * 1.2
  p1 <- ggplot(ivw_sig, aes(x=b, y=outcome)) +
    geom_vline(xintercept=0, linetype="dashed", color="grey60", linewidth=0.3) +
    geom_errorbarh(aes(xmin=ci_lo, xmax=ci_up), height=0, linewidth=0.5, color="grey50") +
    geom_point(aes(color=neg_log10_FDR), size=2.5) +
    scale_color_gradient(low="#FDAE61", high="#D73027", name=expression(-log[10](FDR))) +
    scale_x_continuous(limits=c(-x_lim, x_lim)) +
    labs(title=NULL, x="IVW Beta (95% CI)", y=NULL) +
    base_theme + theme(legend.position="bottom")

  # ---- Panel 2: Heterogeneity ----
  ivw_sig$Qp <- sprintf("%.3f", ivw_sig$Q_pval)
  ivw_sig$Qs <- ifelse(ivw_sig$Q_pval>=0.05, "OK (p>=0.05)", "Het p<0.05")
  p2 <- ggplot(ivw_sig, aes(x=1, y=outcome)) +
    geom_tile(aes(fill=Qs), color="white", width=0.9) +
    geom_text(aes(label=Qp), size=2.5) +
    scale_fill_manual(values=c("OK (p>=0.05)"="#A1D99B","Het p<0.05"="#FC8D59"), name=NULL) +
    labs(title=expression(Hetero~italic(Q)[pval]), x=NULL, y=NULL) +
    base_theme + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                       axis.title.x=element_blank(), legend.position="none")

  # ---- Panel 3: Pleiotropy ----
  ivw_sig$Pp <- sprintf("%.3f", ivw_sig$pleiotropy_pval)
  ivw_sig$Ps <- ifelse(ivw_sig$pleiotropy_pval>=0.05, "OK (p>=0.05)", "Pleio p<0.05")
  p3 <- ggplot(ivw_sig, aes(x=1, y=outcome)) +
    geom_tile(aes(fill=Ps), color="white", width=0.9) +
    geom_text(aes(label=Pp), size=2.5) +
    scale_fill_manual(values=c("OK (p>=0.05)"="#A1D99B","Pleio p<0.05"="#FC8D59"), name=NULL) +
    labs(title=expression(Pleio~intercept[pval]), x=NULL, y=NULL) +
    base_theme + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                       axis.title.x=element_blank(), legend.position="none")

  # ---- Combine with alignment ----
  ph <- max(4, nrow(ivw_sig)*0.22)
  pall <- cowplot::plot_grid(p1, p2, p3, ncol=3, rel_widths=c(4,1.2,1.2),
                             align="h", axis="tb")
  pall <- cowplot::ggdraw() +
    cowplot::draw_label(paste0("Reverse MR: ", exposure_labels[exp_name],
                                " (FDR sig: ", nrow(ivw_sig), "/91)"),
                        x=0.5, y=0.98, size=12, fontface="bold") +
    cowplot::draw_plot(pall, x=0, y=0, width=1, height=0.96)
  out <- paste0("reverse_plots/MR_Plot_", exp_name, ".jpg")
  ggsave(out, plot=pall, width=11, height=ph, dpi=300, limitsize=FALSE)
  total_plots <- total_plots + 1
  cat("  Saved:", out, "\n")
}

# ---- Summary volcano with labels ----
cat("\n=== Summary volcano ===\n")
all_ivw <- data.frame()
for (exp_name in exposures) {
  f <- paste0("./reverse_result/", exp_name, "/effect.csv")
  if (!file.exists(f)) next
  d <- read.csv(f, stringsAsFactors=FALSE)
  d <- d[d$method == "Inverse variance weighted", ]
  d$Exposure <- exp_name
  all_ivw <- rbind(all_ivw, d)
}
all_ivw$neg_log10_pval <- -log10(all_ivw$pval)
all_ivw$sig <- ifelse(all_ivw$FDR<0.05 & !is.na(all_ivw$FDR), "FDR<0.05", "NS")

# 为每个facet选取top15最显著的标注
label_df <- all_ivw %>%
  group_by(Exposure) %>%
  filter(FDR<0.05 & !is.na(FDR)) %>%
  slice_max(order_by=neg_log10_pval, n=15)

pv <- ggplot(all_ivw, aes(x=b, y=neg_log10_pval)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey60", linewidth=0.3) +
  geom_hline(yintercept=-log10(0.05), linetype="dotted", color="grey40", linewidth=0.3) +
  geom_point(aes(color=sig), size=1, alpha=0.6) +
  scale_color_manual(values=c("FDR<0.05"="#D73027","NS"="#4575B4")) +
  geom_text(data=label_df, aes(label=outcome), size=2.3, hjust=-0.1, vjust=0,
            check_overlap=TRUE, color="grey20") +
  facet_wrap(~Exposure, ncol=3, scales="free_x") +
  scale_x_continuous(expand=expansion(mult=c(0.1, 0.35))) +
  labs(title="Reverse MR: All Exposures -> 91 Cytokines (labeled: top15 FDR)", x="IVW Beta",
       y=expression(-log[10](P))) +
  theme_bw(base_size=9) +
  theme(legend.position="bottom", strip.text=element_text(face="bold"),
        panel.grid.minor=element_blank(), plot.margin=margin(10,15,10,10))
ggsave("reverse_plots/MR_Summary_Volcano.jpg", pv, width=10, height=7, dpi=300)
cat("  Saved: MR_Summary_Volcano.pdf\n")

cat("\nDone! Generated", total_plots, "JPG plots + 1 volcano JPG\nLocation: ./reverse_plots/\n")
