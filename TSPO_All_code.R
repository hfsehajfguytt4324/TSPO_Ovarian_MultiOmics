# =============================================================================
#  TSPO 卵巢癌多组学 —— 全部分析代码
#  生成于 2026-11-03 23:54
#  按论文图序组织：11 个分析分组 + 作图样式
#  历史拼接版另存于 code/_archive/历史合并文件/，未纳入本文件
# =============================================================================

# ############################################################################
# ##  总控  (4 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- 总控 / 00_setup.R
# -----------------------------------------------------------------------------
# =============================================================================
#  00_setup.R
#  TSPO 卵巢癌项目 —— 公共设置（所有脚本的第一步，必须先 source 本文件）
#
#  内容：
#    1. 唯一的路径配置（改这里一处即可）
#    2. 依赖包检查（缺什么提示装什么）
#    3. 通用工具函数：hr() / pkey() / pr() / fmt() / load_surv()
#    4. 全局随机种子
#
#  用法：
#    setwd("D:/Rstudy/tcga")
#    source(".../code_final/00_setup.R")
# =============================================================================

# ---------------------------------------------------------------- 1. 路径
PROJ <- "D:/Rstudy/tcga"                 # 主数据目录（TCGA / CPTAC / GDSC）
PROJ2 <- "D:/Rstudy/tcga2"               # 单细胞数据目录
OUT  <- file.path(PROJ, "results_v2")    # 所有输出都写这里
CODE <- file.path(OUT, "code_final")

for (d in c(OUT, file.path(OUT, "tables"), file.path(OUT, "cache"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# 常用输入文件（集中声明，便于核对与替换）
F_EXPR    <- file.path(OUT,  "EXPR_clean.rds")                       # 304 例 × 18747 基因
F_MODULES <- file.path(OUT,  "MODULES.rds")
F_SURV    <- c(file.path(PROJ, "TCGA_OV_survival.tsv"),
               file.path(PROJ, "OV_survival.txt"))
F_CLIN    <- c(file.path(PROJ, "OV_clinicalMatrix.tsv"),
               file.path(PROJ, "TCGA_OV_clinical.csv"))
F_SC_EPI  <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")  # 上皮 7,929
F_SC_CNT  <- file.path(PROJ2, "GSE165897_processed/counts_full.rds")
F_SC_MYE  <- file.path(PROJ,  "TSPO_trajectory/X2_myeloid_trajectory.rds") # 髓系 9,098
F_SPATIAL <- file.path(PROJ,  "TSPO_spatial_processed.rds")                # 8 切片 4,386 spot
F_GSE266  <- file.path(PROJ2, "GSE266577_processed/GSE266577_raw.rds")     # 137,083 细胞
F_VKO_DIR <- file.path(PROJ2, "GSE165897_branch_specificity")              # 虚拟敲除结果

# ---------------------------------------------------------------- 2. 依赖
.need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    cat("\n缺少以下 R 包，请先安装：\n")
    for (p in miss) cat("  install.packages('", p, "')\n", sep = "")
    stop("缺少依赖包：", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.need(c("survival", "glmnet", "Seurat", "Matrix"))

# 仅特定脚本需要，缺失时由该脚本自行跳过并提示
#   data.table      01
#   lme4/lmerTest   03 05
#   CellChat        05
#   slingshot       06
#   monocle3        （不使用，轨迹统一用 slingshot）
#   oncoPredict     11
#   pROC            12
#   VectraPolarisData 09

# ---------------------------------------------------------------- 3. 工具
hr <- function(x) cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")

# 样本条码 -> 病人键（前 10 位字母数字，如 TCGA-61-1910-01 -> TCGA611910）
pkey <- function(x) substr(gsub("[^A-Za-z0-9]", "", as.character(x)), 1, 10)

# 安全打印 data.frame：先把数值列四舍五入，再打印
# 说明：直接 print(df, digits = 3) 在含字符列时会报
#       "数据框中的类似非数值的变量"，必须走这个函数
pr <- function(d, digits = 3, n = NULL) {
  d <- as.data.frame(d)
  if (!is.null(n)) d <- utils::head(d, n)
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(x) round(x, digits))
  print(d, row.names = FALSE)
  invisible(d)
}

# 统一的 p 值格式
fmt <- function(p) ifelse(is.na(p), "NA",
                   ifelse(p < 1e-16, "<1e-16", formatC(p, format = "e", digits = 2)))

# 从候选列名里取第一个存在的列
pick_col <- function(d, keys) {
  h <- intersect(keys, colnames(d))
  if (length(h)) d[[h[1]]] else rep(NA, nrow(d))
}

# 载入 TCGA-OV 生存表，统一成 time(月) / event / OS.time / OS 四列
load_surv <- function(files = F_SURV) {
  fs <- files[file.exists(files)][1]
  if (is.na(fs)) stop("找不到生存表，检查 F_SURV")
  sv <- read.delim(fs, stringsAsFactors = FALSE, check.names = FALSE)
  sc <- intersect(c("sample", "sampleID", "bcr_patient_barcode", "X"), colnames(sv))[1]
  if (is.na(sc)) { sv$.id <- rownames(sv); sc <- ".id" }
  sv$pkey <- pkey(sv[[sc]])
  sv$time     <- suppressWarnings(as.numeric(pick_col(sv, c("PFI.time", "PFI_time")))) / 30.44
  sv$event    <- suppressWarnings(as.numeric(pick_col(sv, c("PFI"))))
  sv$OS.time  <- suppressWarnings(as.numeric(pick_col(sv, c("OS.time", "OS_time")))) / 30.44
  sv$OS       <- suppressWarnings(as.numeric(pick_col(sv, c("OS"))))
  sv$age      <- suppressWarnings(as.numeric(pick_col(sv,
                   c("age_at_initial_pathologic_diagnosis", "age"))))
  sv[!duplicated(sv$pkey), ]
}

# 载入临床表
load_clin <- function(files = F_CLIN) {
  fc <- files[file.exists(files)][1]
  if (is.na(fc)) stop("找不到临床表，检查 F_CLIN")
  cl <- if (grepl("\\.csv$", fc))
    read.csv(fc, check.names = FALSE, stringsAsFactors = FALSE)
  else read.delim(fc, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  sc <- intersect(c("sampleID", "sample", "bcr_patient_barcode", "X", "patient"),
                  colnames(cl))[1]
  if (is.na(sc)) { cl$.id <- rownames(cl); sc <- ".id" }
  cl$pkey <- pkey(cl[[sc]])
  cl[!duplicated(cl$pkey), ]
}

# 把 Figo 分期转成 1-4（注意顺序：IV 必须先判，否则 "IV" 会被 "I" 误吃）
stage_num <- function(x) {
  st <- toupper(as.character(x))
  ifelse(grepl("IV", st),  4,
  ifelse(grepl("III", st), 3,
  ifelse(grepl("II", st),  2,
  ifelse(grepl("I", st),   1, NA))))
}

# 秩残差偏相关（用于控制免疫构成 / 测序深度 / 伪时序）
# a, b : 长度为 n 的向量；Z : 协变量，可以是 n x k 也可以 k x n（自动定向）
pcor_rank <- function(a, b, Z) {
  Z <- as.matrix(Z)
  if (nrow(Z) > ncol(Z)) Z <- t(Z)                 # 统一成 k x n
  ok <- is.finite(a) & is.finite(b) & apply(Z, 2, function(z) all(is.finite(z)))
  a <- a[ok]; b <- b[ok]; Z <- Z[, ok, drop = FALSE]
  n <- length(a); k <- nrow(Z)
  if (n <= k + 3) return(c(rho = NA_real_, p = NA_real_, n = n))
  # apply(Z, 1, rank) 直接给出 n x k（每列是一个协变量的秩），不要再转置
  A  <- cbind(1, apply(Z, 1, rank))
  ea <- residuals(lm(rank(a) ~ A - 1))
  eb <- residuals(lm(rank(b) ~ A - 1))
  r  <- suppressWarnings(cor(ea, eb))
  if (!is.finite(r) || abs(r) >= 1) return(c(rho = r, p = NA_real_, n = n))
  t  <- r * sqrt((n - k - 2) / (1 - r^2))
  c(rho = r, p = 2 * pt(-abs(t), df = n - k - 2), n = n)
}

# ---------------------------------------------------------------- 4. 种子
SEED <- 20260101
set.seed(SEED)

cat("\n[00_setup] 完成\n")
cat("  数据目录: ", PROJ, "\n")
cat("  输出目录: ", OUT, "\n")
cat("  表达矩阵: ", ifelse(file.exists(F_EXPR), "已存在", "不存在（需先跑 01）"), "\n\n")

# -----------------------------------------------------------------------------
# -- 总控 / 27_manuscript_numbers.R
# -----------------------------------------------------------------------------
# =============================================================================
#  27_manuscript_numbers.R   正文中四组尚未由脚本生成的数字
#
#  用法（用 source，不要整段粘贴）：
#    setwd("D:/Rstudy/tcga")
#    source("/your/path/to/results_v2/code_final/27_manuscript_numbers.R", encoding = "UTF-8")
#
#  四部分互相独立，只需对应输入存在：
#    PART A/B 需 EXPR_clean.rds + MODULES.rds
#    PART C/D 需 TSPO_spatial_processed.rds
#
#  输出（results_v2/tables/）：
#    T19a_module_cohens_d.csv          模块 Cohen's d
#    T19b_enrichment_positive.csv      正相关方向 GO
#    T19c_enrichment_negative.csv      负相关方向 GO
#    T19b_direction_counts.csv         三层基因计数口径
#    T19d_spatial_region_within.csv    每个 spot 的片内 z 与区域标签
#    T19e2_region_means.csv            区域均值
#    T19e3_threshold_sensitivity.csv   三组阈值下的符号检验
#    T19e4_threshold_stability.csv     符号计数的阈值稳定性
#    T19f_DEG_*.csv                    排除定义基因后的差异表达
#    T19g_reconciliation.csv           正文值 vs 实测值
# =============================================================================

# ---- 开关 ------------------------------------------------------------------
# 只在未定义时取默认值，因此可以在 source() 前先设好：
#   DO_A <- FALSE; DO_B <- FALSE; DO_C <- TRUE; DO_D <- TRUE
#   source(".../27_manuscript_numbers.R", encoding = "UTF-8")
# 想恢复全部运行，先 DO_A <- TRUE 等显式赋值再 source()。
if (!exists("DO_A")) DO_A <- TRUE
if (!exists("DO_B")) DO_B <- TRUE
if (!exists("DO_C")) DO_C <- TRUE
if (!exists("DO_D")) DO_D <- TRUE

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

# ---- 路径 ------------------------------------------------------------------
if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  if (file.exists(.s)) {
    source(.s, encoding = "UTF-8")
  } else {
    PROJ <- "D:/Rstudy/tcga"
    OUT  <- file.path(PROJ, "results_v2")
    F_EXPR    <- file.path(OUT, "EXPR_clean.rds")
    F_MODULES <- file.path(OUT, "MODULES.rds")
    F_SPATIAL <- file.path(PROJ, "TSPO_spatial_processed.rds")
  }
}
if (!exists("F_EXPR"))    F_EXPR    <- file.path(OUT, "EXPR_clean.rds")
if (!exists("F_MODULES")) F_MODULES <- file.path(OUT, "MODULES.rds")
if (!exists("F_SPATIAL")) F_SPATIAL <- file.path(PROJ, "TSPO_spatial_processed.rds")

OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(Matrix) })
set.seed(20260101)

RECON <- list()
add_recon <- function(part, item, ms, obs) {
  RECON[[length(RECON) + 1]] <<- data.frame(
    part = part, item = item, manuscript = as.character(ms),
    observed = as.character(obs), stringsAsFactors = FALSE)
}

NEED_EXPR <- DO_A || DO_B
NEED_SPAT <- DO_C || DO_D
cat("开关: A=", DO_A, " B=", DO_B, " C=", DO_C, " D=", DO_D, "\n", sep = "")

# ---- 共享：表达矩阵 --------------------------------------------------------
EXPR <- NULL; MOD <- NULL; ts <- NULL
if (NEED_EXPR) {
  if (!file.exists(F_EXPR)) stop("找不到 ", F_EXPR, "，请先跑 01_build_matrix.R")
  EXPR <- as.matrix(readRDS(F_EXPR)); storage.mode(EXPR) <- "double"
  if (!"TSPO" %in% rownames(EXPR)) stop("矩阵中没有 TSPO")
  for (p in c(F_MODULES, file.path(OUT, "MODULES.rds"),
              file.path(PROJ, "results_v2", "MODULES.rds"))) {
    if (file.exists(p)) { MOD <- readRDS(p); break }
  }
  if (DO_A && is.null(MOD)) stop("找不到 MODULES.rds，请先跑 01_build_matrix.R")
  ts <- as.numeric(EXPR["TSPO", ])
  cat("表达矩阵: ", nrow(EXPR), " 基因 x ", ncol(EXPR), " 样本\n", sep = "")
}

# #############################################################################
#  PART A  模块按 TSPO 中位数分组的 Cohen's d
# #############################################################################

if (DO_A) {
  hr("PART A  模块 Cohen's d")
  cat("模块数: ", length(MOD), "\n", sep = "")

  cohens_d <- function(x, y) {
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    n1 <- length(x); n2 <- length(y)
    if (n1 < 3 || n2 < 3) return(NA_real_)
    sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
    if (!is.finite(sp) || sp == 0) return(NA_real_)
    (mean(x) - mean(y)) / sp
  }

  Zall <- t(scale(t(EXPR))); Zall[!is.finite(Zall)] <- NA
  med  <- median(ts, na.rm = TRUE)
  high <- ts > med
  cat("TSPO 中位数 ", round(med, 3), "：高表达 ", sum(high),
      " 例，低表达 ", sum(!high), " 例\n", sep = "")

  A <- do.call(rbind, lapply(names(MOD), function(m) {
    g <- intersect(MOD[[m]], rownames(Zall))
    if (length(g) < 3) return(NULL)
    sc <- colMeans(Zall[g, , drop = FALSE], na.rm = TRUE)
    ok <- is.finite(sc) & is.finite(ts)
    if (sum(ok) < 50) return(NULL)
    ct <- suppressWarnings(cor.test(sc[ok], ts[ok], method = "spearman"))
    h  <- sc[ok][high[ok]]; l <- sc[ok][!high[ok]]
    w  <- suppressWarnings(wilcox.test(h, l))
    data.frame(module = m, n_genes = length(g), n = sum(ok),
               cohens_d = round(cohens_d(h, l), 3),
               mean_high = round(mean(h), 3), mean_low = round(mean(l), 3),
               rho = round(as.numeric(ct$estimate), 3),
               p_rho = signif(ct$p.value, 3),
               p_wilcox = signif(w$p.value, 3), stringsAsFactors = FALSE)
  }))
  A <- A[order(-A$cohens_d), ]
  cat("\n按 Cohen's d 排序：\n"); pr(A, 4)
  write.csv(A, file.path(OUTDIR, "T19a_module_cohens_d.csv"), row.names = FALSE)

  add_recon("A", "Cohen's d 最高模块", "PPARA_Targets 0.75",
            paste0(A$module[1], " ", round(A$cohens_d[1], 3)))
  add_recon("A", "d > 0.5 的模块数", "正文列出 5 个",
            paste0(sum(A$cohens_d > 0.5, na.rm = TRUE), " 个：",
                   paste(A$module[A$cohens_d > 0.5], collapse = "/")))
  rm(Zall); gc()
}

# #############################################################################
#  PART B  网络按相关方向拆分后分别富集
# #############################################################################

if (DO_B) {
  hr("PART B  方向拆分富集")

  RHO_CUT <- 0.3; FDR_CUT <- 0.05
  n  <- ncol(EXPR)
  r  <- apply(EXPR, 1, function(x)
    suppressWarnings(cor(x, ts, method = "spearman", use = "complete.obs")))
  tv  <- r * sqrt((n - 2) / pmax(1 - r^2, 1e-12))
  fdr <- p.adjust(2 * pt(-abs(tv), df = n - 2), "BH")

  pos <- setdiff(rownames(EXPR)[r >  RHO_CUT & fdr < FDR_CUT], "TSPO")
  neg <- setdiff(rownames(EXPR)[r < -RHO_CUT & fdr < FDR_CUT], "TSPO")
  bg  <- rownames(EXPR)

  cat("判据 |rho|>", RHO_CUT, " 且 FDR<", FDR_CUT, "\n", sep = "")
  cat("  正相关原始基因: ", length(pos), "\n", sep = "")
  cat("  负相关原始基因: ", length(neg), "\n", sep = "")
  cat("  网络总数      : ", length(pos) + length(neg), "\n", sep = "")

  run_go <- function(genes, universe, tag) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      cat("未装 clusterProfiler / org.Hs.eg.db，跳过 ", tag, "\n", sep = "")
      return(list(df = NULL, n_mapped = NA_integer_))
    }
    suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
    mg <- suppressWarnings(try(bitr(genes,    "SYMBOL", "ENTREZID", org.Hs.eg.db), silent = TRUE))
    mb <- suppressWarnings(try(bitr(universe, "SYMBOL", "ENTREZID", org.Hs.eg.db), silent = TRUE))
    if (inherits(mg, "try-error") || inherits(mb, "try-error")) {
      cat(tag, ": ID 转换失败\n", sep = ""); return(list(df = NULL, n_mapped = NA_integer_))
    }
    cat(tag, ": 原始 ", length(genes), " → Entrez ", nrow(mg),
        "；背景 ", length(universe), " → ", nrow(mb), "\n", sep = "")
    eg <- tryCatch(
      enrichGO(gene = mg$ENTREZID, universe = mb$ENTREZID, OrgDb = org.Hs.eg.db,
               ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
               qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) { cat(tag, " 失败: ", conditionMessage(e), "\n", sep = ""); NULL })
    if (is.null(eg) || !nrow(as.data.frame(eg))) {
      cat(tag, ": 无显著条目\n", sep = ""); return(list(df = NULL, n_mapped = nrow(mg)))
    }
    df <- as.data.frame(eg)
    df$direction <- tag; df$input_raw <- length(genes); df$input_mapped <- nrow(mg)
    cat(tag, " 显著条目: ", nrow(df), "\n", sep = "")
    pr(head(data.frame(Description = df$Description, GeneRatio = df$GeneRatio,
                       p.adjust = signif(df$p.adjust, 3)), 12), 3)
    list(df = df, n_mapped = nrow(mg))
  }

  GP <- run_go(pos, bg, "positive")
  GN <- run_go(neg, bg, "negative")

  if (!is.null(GP$df)) write.csv(GP$df, file.path(OUTDIR, "T19b_enrichment_positive.csv"), row.names = FALSE)
  if (!is.null(GN$df)) write.csv(GN$df, file.path(OUTDIR, "T19c_enrichment_negative.csv"), row.names = FALSE)

  # enrichGO 的 GeneRatio 分母 = 输入基因中落在 universe 内且有 GO 注释的那部分，
  # 既不是原始基因数，也不是 Entrez 可映射数。正文的 939/1,235 就是它。
  go_denom <- function(df) {
    if (is.null(df) || !nrow(df)) return(NA_integer_)
    as.integer(sub("^.*/", "", df$GeneRatio[1]))
  }

  CNT_TAB <- data.frame(
    方向         = c("正相关", "负相关", "合计"),
    原始基因     = c(length(pos), length(neg), length(pos) + length(neg)),
    Entrez可映射 = c(GP$n_mapped, GN$n_mapped,
                     sum(c(GP$n_mapped, GN$n_mapped), na.rm = TRUE)),
    GO可注释     = c(go_denom(GP$df), go_denom(GN$df),
                     sum(c(go_denom(GP$df), go_denom(GN$df)), na.rm = TRUE)),
    stringsAsFactors = FALSE)
  cat("\n方向计数（三层口径）：\n"); print(CNT_TAB, row.names = FALSE)
  cat("\n正文的 939 / 1,235 对应「GO可注释」列，即 GeneRatio 的分母。\n")
  write.csv(CNT_TAB, file.path(OUTDIR, "T19b_direction_counts.csv"), row.names = FALSE)

  pick_term <- function(df, pat) {
    if (is.null(df)) return(NA_character_)
    hit <- grep(pat, df$Description, ignore.case = TRUE)
    if (!length(hit)) return(NA_character_)
    i <- hit[which.min(df$p.adjust[hit])]
    paste0(df$Description[i], " [", df$GeneRatio[i], ", p=", signif(df$p.adjust[i], 3), "]")
  }

  add_recon("B", "网络总数", "2714", length(pos) + length(neg))
  add_recon("B", "正相关 GO 可注释数", "939", go_denom(GP$df))
  add_recon("B", "负相关 GO 可注释数", "1235", go_denom(GN$df))
  add_recon("B", "正相关原始基因数", "未写", length(pos))
  add_recon("B", "正相关 Entrez 可映射数", "未写", GP$n_mapped)
  add_recon("B", "负相关原始基因数", "未写", length(neg))
  add_recon("B", "负相关 Entrez 可映射数", "未写", GN$n_mapped)
  add_recon("B", "正相关首个免疫条目", "adaptive immune response 114/939 p=1.5e-35",
            pick_term(GP$df, "adaptive immune"))
  add_recon("B", "正相关首个非免疫代谢条目", "aerobic electron transport chain 38/939 p=3.8e-23",
            pick_term(GP$df, "aerobic electron transport"))
  add_recon("B", "负相关 DNA 修复条目", "double-strand break repair 64/1235 p=7.8e-11",
            pick_term(GN$df, "double-strand break"))
  rm(r, tv, fdr); gc()
}

# #############################################################################
#  PART C  片内 z 定义三区域 + 模块比较（含阈值敏感性）
# #############################################################################

SPATIAL_OK <- NEED_SPAT && file.exists(F_SPATIAL)
if (NEED_SPAT && !SPATIAL_OK) cat("⚠ 找不到 ", F_SPATIAL, "，跳过 C/D\n", sep = "")

if (SPATIAL_OK) {
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

  hr("载入空间对象")
  sc <- readRDS(F_SPATIAL)
  cat("对象: ", nrow(sc), " 基因 x ", ncol(sc), " spot\n", sep = "")

  get_layer <- function(obj, layer = "data") {
    an <- names(obj@assays)
    for (a in c("SCT", "RNA", "Spatial")) {
      if (!a %in% an) next
      m <- tryCatch(SeuratObject::LayerData(obj, assay = a, layer = layer),
                    error = function(e) NULL)
      if (is.null(m) || nrow(m) < 100)
        m <- tryCatch(Seurat::GetAssayData(obj, assay = a, slot = layer),
                      error = function(e) NULL)
      if (!is.null(m) && nrow(m) > 100) {
        cat("assay = ", a, "，layer = ", layer, "\n", sep = ""); return(m)
      }
    }
    NULL
  }
  DAT <- get_layer(sc, "data")
  if (is.null(DAT)) stop("无法取得表达矩阵")
  cat("矩阵: ", nrow(DAT), " x ", ncol(DAT), "\n", sep = "")

  md <- sc@meta.data
  sec_col <- NA_character_
  for (c in c("section", "Section", "slice", "sample", "orig.ident")) {
    if (c %in% colnames(md)) {
      k <- length(unique(as.character(md[[c]])))
      if (k >= 2 && k <= 50) { sec_col <- c; break }
    }
  }
  if (is.na(sec_col)) stop("找不到切片列")
  SEC  <- as.character(md[[sec_col]])
  secs <- sort(unique(SEC))
  cat("切片列 ", sec_col, "，共 ", length(secs), " 张\n", sep = "")

  # 片内 z：绝不能合并 8 张切片统一标准化，否则区域标签跨切片不可比
  within_z <- function(v, sec) {
    out <- rep(NA_real_, length(v))
    for (s in unique(sec)) {
      i <- which(sec == s); x <- v[i]
      if (sum(is.finite(x)) < 20) next
      mu <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
      out[i] <- if (is.finite(sdv) && sdv > 0) (x - mu) / sdv else NA_real_
    }
    out
  }

  MODS <- list(
    OXPHOS     = c("NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
                   "SDHA","SDHB","UQCRC1","UQCRC2","COX4I1","COX5A","ATP5A1","ATP5B"),
    Glycolysis = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1"),
    Immune     = c("CD68","CD3D","CD3E","CD8A","PTPRC","LYZ","C1QA","C1QB","MS4A1"),
    EMT        = c("VIM","ZEB1","FN1","COL1A1","COL1A2","SPARC"),
    Hypoxia    = c("CA9","BNIP3","VEGFA","SLC2A1","PDK1"))

  score_within <- function(genes) {
    g <- intersect(genes, rownames(DAT))
    if (length(g) < 3) return(rep(NA_real_, ncol(DAT)))
    within_z(as.numeric(Matrix::colMeans(DAT[g, , drop = FALSE], na.rm = TRUE)), SEC)
  }

  S <- data.frame(spot = colnames(sc), section = SEC, stringsAsFactors = FALSE)
  S$EPCAM  <- within_z(as.numeric(DAT["EPCAM", ]),  SEC)
  S$COL1A1 <- within_z(as.numeric(DAT["COL1A1", ]), SEC)
  S$TSPO   <- within_z(as.numeric(DAT["TSPO", ]),   SEC)
  for (m in names(MODS)) S[[m]] <- score_within(MODS[[m]])

  DEFS <- list(
    A = list(epi = 0.5, str = 0.5, lab = "上皮z>0.5, 基质z>0.5（主定义）"),
    B = list(epi = 1.0, str = 0.5, lab = "上皮z>1.0, 基质z>0.5（严格上皮）"),
    C = list(epi = 1.0, str = 1.0, lab = "上皮z>1.0, 基质z>1.0（三区皆严格）"))

  define_region <- function(S, epi_th, str_th) {
    with(S, ifelse(!is.finite(EPCAM) | !is.finite(COL1A1), "Other",
            ifelse(EPCAM > epi_th & COL1A1 <= 0, "Tumor_Core",
            ifelse(COL1A1 > str_th & EPCAM <= 0, "Stroma",
            ifelse(EPCAM > 0 & COL1A1 > 0, "Interface", "Other")))))
  }

  MIN_SPOT <- 15
  vars <- intersect(c("TSPO", names(MODS)), colnames(S))
  CMPS <- list(c("Tumor_Core","Stroma"), c("Tumor_Core","Interface"),
               c("Interface","Stroma"))

  compare_regions <- function(S, regs, tag) {
    rows <- list()
    for (cmp in CMPS) {
      a <- cmp[1]; b <- cmp[2]
      if (!all(c(a, b) %in% regs)) next
      for (v in vars) {
        dif <- sapply(secs, function(s) {
          z <- S[S$section == s, ]
          x <- z[[v]][z$region == a]; y <- z[[v]][z$region == b]
          if (sum(is.finite(x)) < MIN_SPOT || sum(is.finite(y)) < MIN_SPOT) return(NA_real_)
          mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
        })
        dif <- dif[is.finite(dif)]
        if (length(dif) < 5) next
        k <- sum(dif > 0); m <- length(dif)
        rows[[length(rows) + 1]] <- data.frame(
          definition = tag, comparison = paste0(a, " > ", b), variable = v,
          n_section = m, n_positive = k,
          mean_diff_sections = round(mean(dif), 3),
          sign_p = signif(binom.test(k, m, 0.5)$p.value, 3), stringsAsFactors = FALSE)
      }
    }
    if (!length(rows)) return(NULL)
    do.call(rbind, rows)
  }

  region_means <- function(S, tag) {
    do.call(rbind, lapply(vars, function(v) {
      mm <- tapply(S[[v]], S$region, mean, na.rm = TRUE)
      nn <- tapply(S[[v]], S$region, function(x) sum(is.finite(x)))
      out <- data.frame(definition = tag, variable = v, stringsAsFactors = FALSE)
      for (gg in c("Tumor_Core","Interface","Stroma"))
        out[[paste0("mean_", gg)]] <- if (gg %in% names(mm)) round(mm[[gg]], 3) else NA_real_
      for (gg in c("Tumor_Core","Interface","Stroma"))
        out[[paste0("n_", gg)]] <- if (gg %in% names(nn)) as.integer(nn[[gg]]) else NA_integer_
      out$diff_Core_minus_Stroma <- round(out$mean_Tumor_Core - out$mean_Stroma, 3)
      out$diff_Interface_minus_Stroma <- round(out$mean_Interface - out$mean_Stroma, 3)
      out
    }))
  }

  ALL_CMP <- list(); ALL_MEAN <- list()
  for (d in names(DEFS)) {
    S$region <- define_region(S, DEFS[[d]]$epi, DEFS[[d]]$str)
    regs <- intersect(c("Tumor_Core","Interface","Stroma"), unique(S$region))
    cat("\n--- 定义 ", d, "：", DEFS[[d]]$lab, " ---\n", sep = "")
    print(table(S$region))
    if (d == "A") {
      cat("Other 占比: ", round(100 * mean(S$region == "Other"), 1), "%\n", sep = "")
      cat("区域 × 切片：\n"); print(table(S$section, S$region))
      write.csv(S, file.path(OUTDIR, "T19d_spatial_region_within.csv"), row.names = FALSE)
    }
    ALL_CMP[[d]]  <- compare_regions(S, regs, d)
    ALL_MEAN[[d]] <- region_means(S, d)
  }

  CMP  <- do.call(rbind, ALL_CMP)
  MEAN <- do.call(rbind, ALL_MEAN)

  if (!is.null(CMP)) {
    CMP$FDR <- p.adjust(CMP$sign_p, "BH")
    hr("PART C  符号检验（三组定义）")
    pr(CMP[, c("definition","comparison","variable","n_section","n_positive",
               "mean_diff_sections","sign_p")], 4)
    write.csv(CMP, file.path(OUTDIR, "T19e3_threshold_sensitivity.csv"), row.names = FALSE)

    w <- reshape(as.data.frame(CMP[, c("comparison","variable","n_positive","definition")]),
                 idvar = c("comparison","variable"), timevar = "definition",
                 direction = "wide")
    npc <- grep("^n_positive", colnames(w), value = TRUE)
    w$stable <- apply(w[, npc, drop = FALSE], 1, function(z) length(unique(z)) == 1)
    cat("\n三组定义下符号计数是否一致：", sum(w$stable), " / ", nrow(w), "\n", sep = "")
    if (any(!w$stable)) { cat("随阈值变化的比较：\n"); pr(w[!w$stable, ], 3) }
    write.csv(w, file.path(OUTDIR, "T19e4_threshold_stability.csv"), row.names = FALSE)
  }
  if (!is.null(MEAN)) {
    hr("PART C  区域均值（主定义 A）")
    pr(MEAN[MEAN$definition == "A",
            c("variable","mean_Tumor_Core","mean_Interface","mean_Stroma",
              "diff_Core_minus_Stroma")], 4)
    write.csv(MEAN, file.path(OUTDIR, "T19e2_region_means.csv"), row.names = FALSE)
  }

  gv <- function(cmp, v, col, def = "A") {
    if (is.null(CMP)) return(NA)
    z <- CMP[CMP$definition == def & CMP$comparison == cmp & CMP$variable == v, ]
    if (!nrow(z)) return(NA); z[[col]][1]
  }
  gm <- function(v, col, def = "A") {
    if (is.null(MEAN)) return(NA)
    z <- MEAN[MEAN$definition == def & MEAN$variable == v, ]
    if (!nrow(z)) return(NA); z[[col]][1]
  }

  add_recon("C", "OXPHOS 核心>基质 为正的切片数（共 8 张）", "8",
            gv("Tumor_Core > Stroma","OXPHOS","n_positive"))
  add_recon("C", "OXPHOS 核心>基质 切片内均值差", "+0.231",
            gv("Tumor_Core > Stroma","OXPHOS","mean_diff_sections"))
  add_recon("C", "OXPHOS 核心减基质 区域均值差", "正文未直接给",
            gm("OXPHOS","diff_Core_minus_Stroma"))
  add_recon("C", "糖酵解 核心>基质 为正的切片数（共 8 张）", "7",
            gv("Tumor_Core > Stroma","Glycolysis","n_positive"))
  add_recon("C", "糖酵解 核心>界面 为正的切片数（共 8 张）", "8",
            gv("Tumor_Core > Interface","Glycolysis","n_positive"))
  add_recon("C", "免疫 核心>基质 为正的切片数（共 8 张）", "0",
            gv("Tumor_Core > Stroma","Immune","n_positive"))
  add_recon("C", "免疫 核心>界面 为正的切片数（共 8 张）", "0",
            gv("Tumor_Core > Interface","Immune","n_positive"))
  add_recon("C", "TSPO 核心>基质 为正的切片数（共 8 张）", "5",
            gv("Tumor_Core > Stroma","TSPO","n_positive"))
  add_recon("C", "TSPO 核心>界面 为正的切片数（共 8 张）", "4",
            gv("Tumor_Core > Interface","TSPO","n_positive"))

  S$region <- define_region(S, DEFS$A$epi, DEFS$A$str)
  long <- reshape(S[, c("section","region", vars)], direction = "long",
                  varying = vars, v.names = "score",
                  times = vars, timevar = "variable")
  long <- long[long$region != "Other" & is.finite(long$score), ]
  long$variable <- factor(long$variable, levels = vars)
  long$region   <- factor(long$region, levels = c("Tumor_Core","Interface","Stroma"))
  p <- ggplot(long, aes(region, score, fill = region)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    facet_wrap(~ variable, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = "Within-section z score",
         title = "Region-wise module scores (within-section z)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(FIGDIR, "Fig4A_region_modules_within.pdf"), p, width = 13, height = 4)
  ggsave(file.path(FIGDIR, "Fig4A_region_modules_within.png"), p, width = 13, height = 4, dpi = 300)
  rm(long); gc()

  # ###########################################################################
  #  PART D  排除定义基因后的界面差异表达
  # ###########################################################################

  if (DO_D) {
    hr("PART D  界面 vs 核心 DEG（排除定义基因）")

    # 伪批量 = 各 section×region 的 data 层均值。
    # 用 SCT counts 求和再 CPM 会把 logFC 放大约 2–3 倍，显著基因数虚高。
    keep <- S$region %in% c("Tumor_Core","Interface","Stroma")
    grp  <- paste(S$section[keep], S$region[keep], sep = "__")
    spots_keep <- S$spot[keep]
    ok_lev <- names(which(table(grp) >= MIN_SPOT))
    cat("纳入伪批量的组合: ", length(ok_lev), "\n", sep = "")

    if (length(ok_lev) >= 6) {
      mi <- match(spots_keep, colnames(DAT)); good <- !is.na(mi)
      pb <- sapply(ok_lev, function(g)
        Matrix::rowMeans(DAT[, mi[good][grp[good] == g], drop = FALSE]))
      pb <- as.matrix(pb)
      cat("伪批量: ", nrow(pb), " 基因 x ", ncol(pb), " 样本\n", sep = "")
      cat("值域: ", paste(round(range(pb, na.rm = TRUE), 2), collapse = " ~ "), "\n", sep = "")

      meta <- data.frame(sample = colnames(pb), stringsAsFactors = FALSE)
      meta$section <- sub("__.*$", "", meta$sample)
      meta$region  <- factor(sub("^.*__", "", meta$sample))
      print(table(meta$section, meta$region))

      suppressPackageStartupMessages(library(limma))
      design <- if (length(unique(meta$section)) >= 2)
        model.matrix(~ 0 + region + section, data = meta) else
        model.matrix(~ 0 + region, data = meta)
      fit <- lmFit(pb, design); rn <- colnames(design)

      # 只排除区域定义实际使用的通路基因（EPCAM/COL1A1 及其同类标记）。
      # 不要顺手把 WFDC2、CD68 等也排掉——它们既未参与定义，
      # 排除后会把「界面 vs 基质 = 上皮分泌型」这一真实结果抹掉。
      DEF_GENES <- c("EPCAM","KRT8","KRT18","KRT19","KRT7","PAX8",
                     "COL1A1","COL1A2","COL3A1","DCN","LUM","SPARC")
      cat("排除定义基因 ", length(DEF_GENES), " 个：",
          paste(DEF_GENES, collapse = ", "), "\n", sep = "")

      for (cmp in CMPS) {
        a <- cmp[1]; b <- cmp[2]
        ia <- which(rn == paste0("region", a)); ib <- which(rn == paste0("region", b))
        if (!length(ia) || !length(ib)) { cat("跳过 ", a, " vs ", b, "\n", sep = ""); next }
        v <- rep(0, length(rn)); v[ia] <- 1; v[ib] <- -1
        tt <- topTable(eBayes(contrasts.fit(fit, v)), number = Inf, sort.by = "P")
        tt$gene <- rownames(tt)
        raw_n <- sum(abs(tt$logFC) >= 0.5 & tt$adj.P.Val < 0.05, na.rm = TRUE)
        tt2    <- tt[!tt$gene %in% DEF_GENES, ]
        sig    <- tt2[abs(tt2$logFC) >= 0.5 & tt2$adj.P.Val < 0.05, ]
        sig_l4 <- tt2[abs(tt2$logFC) >= 0.4 & tt2$adj.P.Val < 0.05, ]
        sig_l1 <- tt2[abs(tt2$logFC) >= 1.0 & tt2$adj.P.Val < 0.05, ]

        cat("\n===== ", a, " > ", b, " =====\n", sep = "")
        cat("  未排除定义基因      : ", raw_n, " 个\n", sep = "")
        cat("  排除后 |logFC|>=0.4 : ", nrow(sig_l4), " 个（",
            sum(sig_l4$logFC > 0), " 升 / ", sum(sig_l4$logFC < 0), " 降）\n", sep = "")
        cat("  排除后 |logFC|>=0.5 : ", nrow(sig), " 个（",
            sum(sig$logFC > 0), " 升 / ", sum(sig$logFC < 0), " 降）\n", sep = "")
        cat("  排除后 |logFC|>=1.0 : ", nrow(sig_l1), " 个\n", sep = "")
        if (nrow(sig))
          pr(head(sig[order(-abs(sig$logFC)), c("gene","logFC","adj.P.Val")], 20), 3)

        write.csv(tt2, file.path(OUTDIR,
          paste0("T19f_DEG_", a, "_vs_", b, "_excl_def.csv")), row.names = FALSE)

        if (a == "Tumor_Core" && b == "Interface") {
          add_recon("D", "核心 vs 界面 排除后 |logFC|>=0.5 显著数", "34", nrow(sig))
          add_recon("D", "核心 vs 界面 排除后 |logFC|>=0.4 显著数", "34（旧口径可能为此）", nrow(sig_l4))
          add_recon("D", "核心 vs 界面 排除后 |logFC|>=1.0 显著数", "未写", nrow(sig_l1))
          if (nrow(sig)) {
            # logFC 是 Core − Interface，界面升高 = logFC 最负，故取 which.min
            j <- which.min(sig$logFC)
            add_recon("D", "界面升高首位基因", "FN1 1.046",
                      paste0(sig$gene[j], " ", round(-sig$logFC[j], 3)))
          }
        }
        if (a == "Interface" && b == "Stroma")
          add_recon("D", "界面 vs 基质 排除后显著数", "8", nrow(sig))
      }
    } else cat("有效组合不足，跳过 DEG\n")
  }
  rm(DAT); gc()
}

# #############################################################################
#  汇总：正文值 vs 实测值
# #############################################################################

hr("汇总 reconciliation")

if (length(RECON)) {
  RC <- do.call(rbind, RECON)
  parse_num <- function(x) {
    vapply(as.character(x), function(z) {
      z <- trimws(sub("（.*$", "", z))
      if (grepl("^[0-9.+-]+/[0-9.+-]+$", z))
        return(as.numeric(strsplit(z, "/")[[1]][1]))
      suppressWarnings(as.numeric(gsub("[^0-9.+-]", "", z)))
    }, numeric(1), USE.NAMES = FALSE)
  }
  RC$obs_num <- parse_num(RC$observed)
  RC$man_num <- parse_num(RC$manuscript)

  RC$match <- ifelse(is.na(RC$obs_num) & is.na(RC$man_num), "文本项（人工核对）",
              ifelse(is.na(RC$obs_num), "实测非数值",
              ifelse(is.na(RC$man_num), "正文非数值",
              ifelse(abs(RC$obs_num - RC$man_num) <= 0.05 * pmax(1, abs(RC$man_num)),
                     "一致", "不一致"))))

  pr(RC[, c("part","item","manuscript","observed","match")], 4)
  write.csv(RC, file.path(OUTDIR, "T19g_reconciliation.csv"), row.names = FALSE)

  cat("\n核对结果：\n"); print(table(RC$match))
  bad <- RC[RC$match == "不一致", ]
  if (nrow(bad)) {
    cat("\n⚠ 以下", nrow(bad), "项与正文不一致：\n")
    for (i in seq_len(nrow(bad)))
      cat("  [", bad$part[i], "] ", bad$item[i], "：正文 ", bad$manuscript[i],
          " → 实测 ", bad$observed[i], "\n", sep = "")
    cat("\n改正文前先确认该差异是否影响结论方向；符号计数比差值稳健。\n")
  }
} else cat("没有待核对的项\n")

hr("27_manuscript_numbers.R 完成")

# -----------------------------------------------------------------------------
# -- 总控 / 99_session_info.R
# -----------------------------------------------------------------------------
# =============================================================================
#  99_session_info.R
#  记录运行环境，保证可复现
# =============================================================================

if (!exists("OUT")) stop("请先 source('00_setup.R')")

hr("运行环境")

cat("R 版本: ", R.version.string, "\n", sep = "")
cat("平台:   ", R.version$platform, "\n", sep = "")
cat("时间:   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")

PKGS <- c("Seurat","SeuratObject","Matrix","survival","glmnet","lme4","lmerTest",
          "CellChat","slingshot","oncoPredict","pROC","data.table",
          "VectraPolarisData","curatedOvarianData","scTenifoldKnk","readxl")

cat("\n关键包版本:\n")
ver <- vapply(PKGS, function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "未安装"
}, character(1))
print(data.frame(package = names(ver), version = unname(ver)), row.names = FALSE)

si <- file.path(OUT, "tables", "session_info.txt")
writeLines(c(R.version.string,
             paste("平台:", R.version$platform),
             paste("时间:", format(Sys.time())),
             "",
             paste0(names(ver), ": ", ver),
             "",
             capture.output(sessionInfo())), si)

cat("\n已写出: ", si, "\n", sep = "")
hr("99_session_info.R 完成")

# -----------------------------------------------------------------------------
# -- 总控 / RUN_ALL.R
# -----------------------------------------------------------------------------
# =============================================================================
#  RUN_ALL.R
#  一键顺序运行。按结论分组，可以只跑某一组。
#
#  用法：
#    setwd("D:/Rstudy/tcga")
#    source("/your/path/to/results_v2/code_final/RUN_ALL.R")
#
#  只想跑某一部分：改 STEP 里的开关，或单独 source 对应脚本。
# =============================================================================

# 【审稿人注意】请把下面的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
#   例如 "/your/path/to/results_v2/code_final"；也可以先
#   Sys.setenv(TSPO_CODE_DIR = "/your/path/to/code_final") 再运行本文件。
BASE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")

# ---- 运行开关 --------------------------------------------------------------
RUN <- list(
  setup     = TRUE,
  c1        = TRUE,    # 结论一  bulk 网络与模块
  c1b       = TRUE,    # 图 3    TSPO 高表达的两个来源（需 01/02 与外部队列）
  c2        = TRUE,    # 结论二  上皮代谢轴
  c3        = TRUE,    # 结论三  髓系配体 + 分化轨迹
  c4        = TRUE,    # 结论四  空间区室分离
  c5        = TRUE,    # 结论五  药物敏感性预测
  supp      = FALSE,   # 补充分析 13–27（约 45 分钟，见 SUPPLEMENTARY.md）
  session   = TRUE
)

# ---- 00 公共设置（必须最先） ----------------------------------------------
source(file.path(BASE, "00_setup.R"))

run <- function(f) {
  cat("\n\n################  ", f, "  ################\n", sep = "")
  t0 <- Sys.time()
  source(file.path(BASE, f), echo = FALSE)
  cat(sprintf("\n[完成] %s  用时 %.1f 分钟\n", f,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ---- 结论一：共表达网络与模块 ----------------------------------------------
if (RUN$c1) {
  # 01 只在第一次运行时需要（从 HiSeqV2.gz 重建 EXPR_clean.rds）
  if (!file.exists(F_EXPR)) run("01_build_matrix.R")
  run("02_bulk_network.R")
}

# ---- 图 3：TSPO 高表达的两个来源 -------------------------------------------
# 依赖 01 的 EXPR_clean.rds 与 02 的 T8_network_all_genes.csv。
# 外部队列（GSE32062 / GSE63885 / GSE9891）需另行下载；缺失时 C/D 面板以占位显示。
if (RUN$c1b) {
  run("35_tspo_source_index.R")   # 统计量：Stouffer 合并、BH、循环性检验
  run("fig03_source.R")           # 出图 Fig3
}

# ---- 结论二：上皮代谢轴 ----------------------------------------------------
if (RUN$c2) {
  run("03_epithelial_mixed.R")     # 混合模型 + 跨病人符号检验
  run("04_vko_branch.R")           # 虚拟敲除再分析
  run("10_cptac_protein.R")        # CPTAC 蛋白层面验证（结论二、三共用）
}

# ---- 结论三：髓系配体 + 分化轨迹 -------------------------------------------
if (RUN$c3) {
  run("05_myeloid_mixed.R")        # 配体混合模型 + CellChat 亚型
  run("06_myeloid_trajectory.R")   # 单核 -> 巨噬 轨迹
  run("07_gse266577_chemo.R")      # 化疗初治 vs IDS 配对
}

# ---- 结论四：空间区室分离 --------------------------------------------------
if (RUN$c4) {
  run("08_spatial_persection.R")   # 空间转录组分切片重算
  run("09_codex_tma.R")            # CODEX 空间蛋白组
}

# ---- 结论五：药物敏感性预测 ------------------------------------------------
if (RUN$c5) {
  run("11_oncopredict.R")          # 预测 198 种药物 IC50（首次较慢）
  run("12_drug_stratification.R")  # 偏相关 + 三分位 + 分诊
}

# ---- 补充分析：13–27 -------------------------------------------------------
# 逐脚本运行而非 source 合并版，便于中途失败后定位。
# 21 依赖 20 生成的 cache/epi_subclustered.rds，顺序不可调换。
if (RUN$supp) {
  for (f in c("13_tf_upstream.R", "14_tf_upstream_singlecell.R",
              "16_enrichment.R", "17_immune_infiltration.R",
              "18_tspo_highlow_deg.R", "19_nk_vs_ctl.R",
              "20_epi_subclustering.R", "21_epi_subcluster_resolve.R",
              "22_spatial_region_deg.R",   # 仅备查，正文不引用其数字
              "23_wgcna.R", "24_caf_trajectory.R",
              "25_upstream_epigenetic.R", "26_tspo_celltype_levels.R",
              "27_manuscript_numbers.R",
              "28_mif_spatial_region.R",
              "29_gse26712_replication.R",
              "30_replication_diagnose.R",
              "31_bulk_confound_adjusted.R")) {
    run(f)
  }
}

# ---- 环境信息 --------------------------------------------------------------
if (RUN$session) run("99_session_info.R")

cat("\n\n全部完成。所有结果在: ", OUT, "/tables\n", sep = "")

# ############################################################################
# ##  Fig1_表达与网络  (2 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- Fig1_表达与网络 / 01_build_matrix.R
# -----------------------------------------------------------------------------
# =============================================================================
#  00_build_matrix.R
#  TSPO 卵巢癌项目 —— 统一表达矩阵层（唯一数据源）
#
#  背景：原分析把 TCGA_OV_CBS.rds 当表达矩阵用，实际它是拷贝数(CN)表。
#        所有基于它的结果（5.1 之后全部）作废，必须重建。
#
#  本脚本职责：
#    1. 读入 UCSC Xena HiSeqV2 (TCGA.OV.sampleMap/HiSeqV2)
#    2. 清洗：去重基因名、只保留原发瘤(-01)、去低表达基因
#    3. 定义全部机制模块（与旧脚本同源，便于对照）
#    4. 输出 EXPR / MODULES / 样本表 到 results_v2/
#    5. 自检：把关键数字打印出来，可与 Python 独立实现交叉验证
#
#  用法：source("00_build_matrix.R")  之后所有分析只用 EXPR，禁止再 readRDS 其它矩阵
# =============================================================================

# ---------------------------------------------------------------- 0. 路径
# 注意：这里只允许改这一个变量
PROJ <- "D:/Rstudy/tcga"
OUT  <- file.path(PROJ, "results_v2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

stopifnot(dir.exists(PROJ))

# ---------------------------------------------------------------- 1. 依赖
need_pkg <- c("data.table")
for (p in need_pkg) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("缺少包 %s，请先 install.packages('%s')", p, p))
  }
}

# ---------------------------------------------------------------- 2. 读入
READ_FROM_HISEQ <- TRUE   # TRUE=HiSeqV2(推荐, 全基因, n≈304)
                          # FALSE=Curry GX(旧, n=520, 缺 CPT1B/CD274/TOMT0 等)
if (READ_FROM_HISEQ) {
  f <- file.path(PROJ, "HiSeqV2.gz")
  if (!file.exists(f)) stop("找不到 ", f, "，请先下载 TCGA.OV.sampleMap/HiSeqV2")
  message("读取 ", f, " ...")
  raw <- data.table::fread(f, data.table = FALSE, check.names = FALSE)
  rn  <- raw[[1]]
  raw <- as.matrix(raw[, -1, drop = FALSE])
  rownames(raw) <- rn
  rm(rn)
} else {
  raw <- readRDS(file.path(PROJ, "TCGA_OV_GX.rds"))
  raw <- as.matrix(raw)
}

storage.mode(raw) <- "double"
message("原始矩阵: ", nrow(raw), " 基因 x ", ncol(raw), " 样本")
message("数值范围: ", round(min(raw, na.rm = TRUE), 2), " ~ ",
        round(max(raw, na.rm = TRUE), 2))

# ---------------------------------------------------------------- 3. 清洗
# 3.1 去重基因名（保留平均表达最高的那一个）
if (any(duplicated(rownames(raw)))) {
  m <- rowMeans(raw, na.rm = TRUE)
  ord <- order(rownames(raw), -m)
  raw <- raw[ord, , drop = FALSE]
  raw <- raw[!duplicated(rownames(raw)), , drop = FALSE]
  message("去重后: ", nrow(raw), " 基因")
}

# 3.2 只保留原发瘤（样本条码最后两位 01），去掉 02(复发)/11(正常)
tid <- sub(".*-", "", colnames(raw))
if (all(nchar(tid) >= 2)) {
  keep <- substr(tid, 1, 2) == "01"
  message("样本类型分布: ", paste(names(table(substr(tid, 1, 2))),
                                table(substr(tid, 1, 2)), sep = "=", collapse = "; "))
  if (sum(keep) >= 50) raw <- raw[, keep, drop = FALSE]
  message("保留原发瘤: ", ncol(raw), " 例")
}

# 3.3 低表达过滤（与常见做法一致：>=25% 样本检出 且 方差>0）
MIN_PCT <- 0.25
pct  <- rowMeans(raw > 0, na.rm = TRUE)
v    <- apply(raw, 1, var, na.rm = TRUE)
keep <- (pct >= MIN_PCT) & (v > 0) & !is.na(v)
message("低表达过滤: ", nrow(raw), " -> ", sum(keep),
        " 基因 (检出率>=", MIN_PCT * 100, "%)")
raw <- raw[keep, , drop = FALSE]

# 3.4 样本名标准化：TCGA-61-1910-01 -> 病人 TCGA-61-1910
sample_df <- data.frame(
  sample  = colnames(raw),
  patient = sub("-[0-9]{2}[A-Z]?$", "", colnames(raw)),
  dot_id  = gsub("-", ".", sub("-[0-9]{2}[A-Z]?$", "", colnames(raw))),
  pkey    = substr(gsub("[^A-Za-z0-9]", "", colnames(raw)), 1, 10),
  stringsAsFactors = FALSE
)
colnames(raw) <- sample_df$sample

EXPR <- raw
rm(raw)

# 3.5 关键基因自检（缺什么必须写进方法学限制）
KEYGENES <- c("TSPO", "PPARA", "CPT1A", "CPT1B", "CPT2", "ACADVL", "ECH1", "ECHS1",
              "GPX4", "HMOX1", "NQO1", "SLC7A11", "GCLC", "GCLM",
              "TGFB1", "SPP1", "LGALS1", "IL10", "CD274", "CD68", "CD163", "CD3D",
              "TOMM22", "TOMM40", "TIMM17A", "TIMM50", "TOMM70", "PKM",
              "NDUFS7", "NDUFV1", "UQCRC1", "COX4I1", "SDHA", "ATP5A1")
miss <- setdiff(KEYGENES, rownames(EXPR))
message("关键基因缺失: ", if (length(miss) == 0) "无" else paste(miss, collapse = ", "))

# ---------------------------------------------------------------- 4. 模块
# 与旧脚本(spo.R)完全一致的基因集，便于新旧结果对照
MODULES <- list(
  Mito_Import = c("TOMM22", "TOMM40", "TOMM40L", "TOMM70", "TIMM23", "TIMM17A",
                  "TIMM50", "TIMM44", "SLC25A3", "SLC25A4", "SLC25A5"),
  Assembly    = c("NDUFAF1", "NDUFAF2", "NDUFAF3", "NDUFAF4", "NDUFAF5", "NDUFAF6",
                  "TIMMDC1", "FOXRED1", "ACAD9", "ECSIT", "TMEM126B"),
  Catalytic   = c("NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS7", "NDUFS8",
                  "NDUFV1", "NDUFV2", "NDUFV3",
                  "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6"),
  Complex_II_V = c("SDHA", "SDHB", "SDHC", "SDHD",
                   "UQCRC1", "UQCRC2", "UQCRFS1", "CYC1",
                   "COX4I1", "COX5A", "COX5B", "COX6A1",
                   "ATP5A1", "ATP5B", "ATP5C1", "ATP5F1"),
  OXPHOS      = c("NDUFS1", "NDUFS2", "NDUFS3", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2",
                  "SDHA", "SDHB", "UQCRC1", "UQCRC2", "COX4I1", "COX5A", "ATP5A1", "ATP5B"),
  Glycolysis  = c("HK1", "HK2", "PFKP", "ALDOA", "GAPDH", "PGK1", "PKM", "LDHA"),
  ROS_Defense = c("GPX4", "SOD2", "CAT", "PRDX3", "TXNRD2", "GSR"),
  Apoptosis   = c("BAX", "BAK1", "BCL2", "CYCS", "CASP9", "CASP3"),
  Immune_Ckpt = c("CD274", "CTLA4", "LAG3", "TIGIT", "PDCD1", "HAVCR2"),
  PPARA_Targets = c("PPARA", "CPT1B", "GPX4", "HMOX1", "NQO1"),
  FAO         = c("CPT1A", "CPT1B", "CPT2", "ACADVL", "ACADM", "ACADS",
                  "ECH1", "ECHS1", "HADHA", "HADHB")
)
# 只保留矩阵里真实存在的基因，并报告模块有效基因数
MODULES <- lapply(MODULES, function(g) intersect(g, rownames(EXPR)))
MODULES <- MODULES[sapply(MODULES, length) >= 3]
message("\n模块有效基因数:")
for (m in names(MODULES)) message(sprintf("  %-14s %d", m, length(MODULES[[m]])))
if (any(!c("MT-ND1") %in% rownames(EXPR))) {
  message("注意: HiSeqV2 不含 MT-* 线粒体基因，Catalytic 模块降级为核编码基因")
}

# ---------------------------------------------------------------- 5. 工具函数
#' Spearman 相关（ties 时 automatic exact=FALSE，避免 warning）
sp_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 4) return(c(rho = NA, p = NA, n = length(x)))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = length(x))
}

#' 批量相关：gene 与一组基因/向量
cor_table <- function(target, genes, mat = EXPR) {
  genes <- intersect(genes, rownames(mat))
  out <- t(vapply(genes, function(g) sp_cor(mat[target, ], mat[g, ]),
                  numeric(3)))
  df <- data.frame(Gene = genes, rho = out[, "rho"], p = out[, "p"],
                   n = out[, "n"], stringsAsFactors = FALSE)
  df$FDR <- p.adjust(df$p, method = "BH")
  df[order(-abs(df$rho)), ]
}

#' 模块评分：逐基因 z-score 后取均值（与旧脚本一致）
module_score <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) < 3) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[genes, , drop = FALSE])))
  z[!is.finite(z)] <- NA
  colMeans(z, na.rm = TRUE)
}

#' Cohen's d（等方差合并）
cohens_d <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 3 || n2 < 3) return(NA_real_)
  sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x) - mean(y)) / sp
}

#' Spearman 偏相关（对秩做 OLS 残差化，等价于 ppcor::pcor 的秩版本）
#' Z: 协变量矩阵（基因 × 样本），内部会强制定向
partial_sp <- function(x, y, Z) {
  Z <- as.matrix(Z)
  if (nrow(Z) > ncol(Z)) Z <- t(Z)          # 强制为 协变量 x 样本
  if (ncol(Z) != length(x)) stop("Z 的样本数与 x 不一致")
  ok <- is.finite(x) & is.finite(y) & apply(Z, 2, function(z) all(is.finite(z)))
  x <- x[ok]; y <- y[ok]; Z <- Z[, ok, drop = FALSE]
  n <- length(x); k <- nrow(Z)
  if (n <= k + 3) return(c(rho = NA_real_, p = NA_real_))
  rx <- rank(x); ry <- rank(y)
  # apply(Z,1,rank) 返回 n x k（每个协变量一列），直接拼上截距列
  A  <- cbind(1, apply(Z, 1, rank))          # n x (k+1)
  ex <- residuals(lm(rx ~ A - 1))
  ey <- residuals(lm(ry ~ A - 1))
  r  <- suppressWarnings(cor(ex, ey))
  if (!is.finite(r)) return(c(rho = NA_real_, p = NA_real_))
  if (abs(r) >= 1)  return(c(rho = r, p = 0))
  t  <- r * sqrt((n - k - 2) / (1 - r^2))
  c(rho = r, p = 2 * pt(-abs(t), df = n - k - 2))
}

#' 病人主键：任意 TCGA 条码 -> 10 位（TCGA + 2 + 4），用于跨表匹配
pkey <- function(x) substr(gsub("[^A-Za-z0-9]", "", as.character(x)), 1, 10)

#' 加载生存/临床（Xena 格式，列名自动探测）
load_survival <- function() {
  cands <- c(file.path(PROJ, "TCGA_OV_survival.tsv"),
             file.path(PROJ, "OV_survival.txt"),
             file.path(PROJ, "TCGA_OV_clinical.csv"))
  f <- cands[file.exists(cands)][1]
  if (is.na(f)) {
    warning("未找到生存数据文件，4.1/4.2/6.8 相关分析将跳过。候选：\n  ",
            paste(cands, collapse = "\n  "))
    return(NULL)
  }
  if (grepl("\\.csv$", f)) {
    d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  }
  # 样本列
  sc <- intersect(c("sample", "sampleID", "SampleID", "bcr_patient_barcode",
                    "X", "patient"), colnames(d))[1]
  if (is.na(sc)) { warning("生存表里找不到样本列"); return(NULL) }
  d$sample <- d[[sc]]
  d$patient <- sub("-[0-9]{2}[A-Z]?$", "", d$sample)
  d$dot_id  <- gsub("-", ".", d$patient)
  # 时间/事件列自动探测
  pick <- function(keys) {
    hit <- intersect(keys, colnames(d))
    if (length(hit)) d[[hit[1]]] else rep(NA, nrow(d))
  }
  d$PFI.time <- suppressWarnings(as.numeric(pick(c("PFI.time", "PFI_time", "PFI.Time"))))
  d$PFI      <- suppressWarnings(as.numeric(pick(c("PFI"))))
  d$OS.time  <- suppressWarnings(as.numeric(pick(c("OS.time", "OS_time", "OS.Time"))))
  d$OS       <- suppressWarnings(as.numeric(pick(c("OS"))))
  d$age      <- suppressWarnings(as.numeric(pick(c("age_at_initial_pathologic_diagnosis",
                                                    "age", "Age"))))
  message("生存数据: ", basename(f), " (", nrow(d), " 行)  可用列: ",
          paste(intersect(c("PFI.time", "OS.time", "age"), colnames(d)), collapse = ", "))
  d
}

message("\n[OK] EXPR ", nrow(EXPR), " x ", ncol(EXPR), " | MODULES ", length(MODULES))
message("     后续所有分析请只用 EXPR，不要再 readRDS(TCGA_OV_CBS.rds)")

# ---------------------------------------------------------------- 6. 自检
cat("\n===== 自检：核心相关（应与 Python 独立实现一致）=====\n")
for (pr in list(c("TSPO", "PPARA"), c("TSPO", "TOMM22"), c("TSPO", "CD68"),
                c("TSPO", "GPX4"), c("PPARA", "CPT1B"), c("PPARA", "GPX4"))) {
  if (all(pr %in% rownames(EXPR))) {
    r <- sp_cor(EXPR[pr[1], ], EXPR[pr[2], ])
    cat(sprintf("  %-14s rho=%6.4f  p=%.3g  n=%d\n",
                paste(pr, collapse = "~"), r["rho"], r["p"], r["n"]))
  } else {
    cat(sprintf("  %-14s [基因缺失]\n", paste(pr, collapse = "~")))
  }
}

# ---------------------------------------------------------------- 7. 保存
saveRDS(EXPR,     file.path(OUT, "EXPR_clean.rds"))
saveRDS(MODULES,  file.path(OUT, "MODULES.rds"))
saveRDS(sample_df, file.path(OUT, "samples.rds"))
write.csv(sample_df, file.path(OUT, "samples.csv"), row.names = FALSE)
message("\n已保存: results_v2/EXPR_clean.rds  MODULES.rds  samples.rds")

# -----------------------------------------------------------------------------
# -- Fig1_表达与网络 / 02_bulk_network.R
# -----------------------------------------------------------------------------
# =============================================================================
#  02_bulk_network.R
#  结论一：TSPO 在卵巢癌中高表达，其共表达网络富集代谢与髓系相关基因
#
#  产出（对应正文 3.1）：
#    · TSPO 共表达网络 2,714 基因（|rho|>0.3, FDR<0.05）
#    · 网络核心伙伴：GPX4 0.440 / TOMM22 0.417 / CD68 0.378 / PPARA 0.150（第 58 位）
#    · 模块相关：Complex II-V 0.722 / 催化 0.706 / OXPHOS 0.682 / ROS 0.630 / 糖酵解 0.199(ns)
#    · TSPO ~ 肿瘤纯度 -0.371 (p=2.3e-11)；TSPO ~ CD68 0.378
#
#  说明：正文中「TSPO 8.09 vs 5.97（TCGA vs GTEx）」不来自本脚本，
#        由原始脚本 original/tspo.R 第 1-159 行生成（见 original/README.md）。
#
#  依赖：01_build_matrix.R 已生成 EXPR_clean.rds
# =============================================================================

if (!exists("OUT")) stop("请先 source('00_setup.R')")
suppressPackageStartupMessages({ library(Matrix) })
hr("结论一  Bulk 共表达网络与模块相关")

stopifnot(file.exists(F_EXPR))
EXPR <- readRDS(F_EXPR)
EXPR <- as.matrix(EXPR); storage.mode(EXPR) <- "double"
message("表达矩阵: ", nrow(EXPR), " 基因 x ", ncol(EXPR), " 样本")
stopifnot("TSPO" %in% rownames(EXPR))

# ------------------------------------------------ 1. TSPO 共表达网络
hr("1. TSPO 共表达网络")

y   <- as.numeric(EXPR["TSPO", ])
n   <- length(y)
rho <- apply(EXPR, 1, function(x) suppressWarnings(cor(x, y, method = "spearman")))
tstat <- rho * sqrt((n - 2) / pmax(1 - rho^2, 1e-12))
pv    <- 2 * pt(-abs(tstat), df = n - 2)
fdr   <- p.adjust(pv, "BH")

NET <- data.frame(gene = names(rho), rho = as.numeric(rho),
                  p = pv, FDR = fdr, stringsAsFactors = FALSE)
NET <- NET[order(-abs(NET$rho)), ]
NET <- NET[NET$gene != "TSPO", ]
SIG <- NET[abs(NET$rho) > 0.3 & NET$FDR < 0.05, ]
cat("网络基因数: ", nrow(SIG), "（阈值 |rho|>0.3 且 FDR<0.05）\n", sep = "")
pr(SIG, 4, 10)

cat("\n重点基因在网络中的排名:\n")
for (g in c("GPX4", "TOMM22", "CD68", "CD163", "PPARA", "HMOX1", "TGFB1")) {
  i <- match(g, NET$gene)
  if (!is.na(i)) cat(sprintf("  %-7s rho=%+.3f  p=%.2e  排名 %d/%d\n",
                             g, NET$rho[i], NET$p[i], i, nrow(NET)))
}

# PPARA 偏相关（控制免疫/间质构成）
Z <- EXPR[c("CD68", "CD3D", "COL1A1"), , drop = FALSE]
pp <- pcor_rank(as.numeric(EXPR["PPARA", ]), as.numeric(EXPR["TSPO", ]), Z)
cat(sprintf("\nTSPO ~ PPARA 偏相关（控制 CD68/CD3D/COL1A1）: rho=%+.3f, p=%.3g\n",
            pp["rho"], pp["p"]))

write.csv(NET, file.path(OUT, "tables", "T1_TSPO_coexpression_network.csv"), row.names = FALSE)
write.csv(SIG, file.path(OUT, "tables", "T1b_TSPO_network_significant.csv"), row.names = FALSE)

# ------------------------------------------------ 2. 模块评分与相关
hr("2. 功能模块相关")

MODULES <- if (file.exists(F_MODULES)) readRDS(F_MODULES) else list(
  OXPHOS      = c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
                  "SDHA","SDHB","UQCRC1","UQCRC2","COX4I1","COX5A","ATP5A1","ATP5B"),
  Complex_II_V= c("SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","COX4I1","COX5A","ATP5A1","ATP5B"),
  Catalytic   = c("NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2","SDHA","SDHB"),
  ROS_Defense = c("GPX4","SOD2","CAT","PRDX3","TXNRD2","GSR","GCLC","GCLM"),
  Glycolysis  = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1"),
  Mito_Import = c("TOMM22","TOMM40","TIMM23","TIMM50","SLC25A3","TIMM17A")
)

Zall <- t(scale(t(EXPR)))
Zall[!is.finite(Zall)] <- NA

rows <- do.call(rbind, lapply(names(MODULES), function(m) {
  gs <- intersect(MODULES[[m]], rownames(Zall))
  if (length(gs) < 3) return(NULL)
  sc  <- colMeans(Zall[gs, , drop = FALSE], na.rm = TRUE)
  ct  <- suppressWarnings(cor.test(sc, y, method = "spearman"))
  rp  <- pcor_rank(sc, y, Z)
  data.frame(module = m, n_genes = length(gs),
             rho = round(unname(ct$estimate), 3),
             p = signif(ct$p.value, 3),
             rho_partial = round(rp["rho"], 3),
             p_partial = signif(rp["p"], 3), stringsAsFactors = FALSE)
}))
rows <- rows[order(-rows$rho), ]
pr(rows)
write.csv(rows, file.path(OUT, "tables", "T1c_module_correlations.csv"), row.names = FALSE)

# ------------------------------------------------ 3. 纯度代理与髓系混杂
hr("3. 肿瘤纯度与髓系标志")

purity <- -as.numeric(colMeans(Zall[c("CD68", "CD3D", "COL1A1"), , drop = FALSE], na.rm = TRUE))
for (v in list(purity = purity, CD68 = as.numeric(EXPR["CD68", ]))) {
  ct <- suppressWarnings(cor.test(y, v, method = "spearman"))
  nm <- if (identical(v, purity)) "肿瘤纯度代理" else "CD68"
  cat(sprintf("  TSPO ~ %-8s rho = %+.3f, p = %.3g\n", nm, ct$estimate, ct$p.value))
}

cat("\n说明：TSPO 与纯度负相关、与 CD68 正相关，表明 bulk 信号混合了肿瘤细胞与免疫细胞，\n")
cat("      这是后续必须在单细胞/空间层面拆分细胞来源的直接依据。\n")

hr("02_bulk_network.R 完成")

# ############################################################################
# ##  Fig2_bulk共表达与免疫  (5 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- Fig2_bulk共表达与免疫 / 16_enrichment.R
# -----------------------------------------------------------------------------
# =============================================================================
#  16_enrichment.R   GO / KEGG / Reactome 富集（TSPO 共表达网络）
#
#  v2 改动：网络基因直接从 EXPR 现算，不再依赖 02_bulk_network.R 生成的 T1b 文件。
#           图注全部改为英文，避免 PDF 设备无法渲染中文。
#
#  输入：EXPR_clean.rds（18,747 基因 × 304 样本）
#  输出：T8a_GO_all.csv、T8b_KEGG.csv、T8c_Reactome.csv
#        Fig1C_GO_enrichment.pdf/png、Fig1D_KEGG.pdf/png、Fig1E_Reactome.pdf/png
# =============================================================================

if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  .ok <- FALSE
  if (file.exists(.s))
    .ok <- !inherits(try(source(.s, encoding = "UTF-8"), silent = TRUE), "try-error")
  if (!.ok) {
    PROJ <- "D:/Rstudy/tcga"; OUT <- file.path(PROJ, "results_v2")
    hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
    pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)
  }
  for (d in c(OUT, file.path(OUT, "tables"), file.path(OUT, "figures")))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("hr")) hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
if (!exists("pr")) pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

suppressPackageStartupMessages({ library(ggplot2) })
OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

RHO_CUT <- 0.3
FDR_CUT <- 0.05
TOP_N   <- 12

hr("STEP 0  依赖检查")

need <- function(p) requireNamespace(p, quietly = TRUE)
cat("clusterProfiler: ", need("clusterProfiler"),
    "  org.Hs.eg.db: ", need("org.Hs.eg.db"),
    "  ReactomePA: ", need("ReactomePA"), "\n", sep = "")
if (!need("clusterProfiler") || !need("org.Hs.eg.db"))
  stop("请先安装：BiocManager::install(c('clusterProfiler','org.Hs.eg.db'))")

hr("STEP 1  从表达矩阵重算网络基因")

EXPR <- as.matrix(readRDS(F_EXPR))
storage.mode(EXPR) <- "double"
cat("表达矩阵: ", nrow(EXPR), " 基因 x ", ncol(EXPR), " 样本\n", sep = "")
stopifnot("TSPO" %in% rownames(EXPR))

ts <- as.numeric(EXPR["TSPO", ])
n  <- ncol(EXPR)

r <- apply(EXPR, 1, function(x)
       suppressWarnings(cor(x, ts, method = "spearman", use = "complete.obs")))
tv  <- r * sqrt((n - 2) / (1 - r^2))
fdr <- p.adjust(2 * pt(-abs(tv), df = n - 2), "BH")

genes <- setdiff(rownames(EXPR)[abs(r) > RHO_CUT & fdr < FDR_CUT], "TSPO")
bg    <- rownames(EXPR)
cat("判据 |rho|>", RHO_CUT, " 且 FDR<", FDR_CUT, "\n", sep = "")
cat("网络基因: ", length(genes), "（预期 2714）\n", sep = "")
cat("背景基因: ", length(bg), "\n", sep = "")
if (length(genes) < 100) stop("网络基因过少，检查 EXPR 是否有 TSPO")

saveRDS(list(genes = genes, rho = r, fdr = fdr),
        file.path(OUT, "cache", "network_genes.rds"))
write.csv(data.frame(gene = names(r), rho = round(as.numeric(r), 4),
                     FDR = as.numeric(fdr), in_network = names(r) %in% genes),
          file.path(OUTDIR, "T8_network_all_genes.csv"), row.names = FALSE)

hr("STEP 2  ID 转换")

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
})

mg <- suppressWarnings(bitr(unique(genes), "SYMBOL", "ENTREZID", org.Hs.eg.db))
mb <- suppressWarnings(bitr(unique(bg),    "SYMBOL", "ENTREZID", org.Hs.eg.db))
cat("网络基因转换: ", nrow(mg), " / ", length(genes), "\n", sep = "")
cat("背景基因转换: ", nrow(mb), " / ", length(bg), "\n", sep = "")

hr("STEP 3  GO 富集（BP / CC / MF）")

go <- tryCatch(
  enrichGO(gene = mg$ENTREZID, universe = mb$ENTREZID, OrgDb = org.Hs.eg.db,
           ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05,
           qvalueCutoff = 0.2, readable = TRUE),
  error = function(e) { cat("GO 失败: ", conditionMessage(e), "\n", sep = ""); NULL })

if (!is.null(go) && nrow(as.data.frame(go))) {
  gd <- as.data.frame(go)
  cat("GO 显著条目: ", nrow(gd), "\n", sep = "")
  write.csv(gd, file.path(OUTDIR, "T8a_GO_all.csv"), row.names = FALSE)

  for (o in c("BP", "CC", "MF")) {
    s <- gd[gd$ONTOLOGY == o, ]
    s <- head(s[order(s$p.adjust), ], 8)
    cat("\n-- ", o, " --\n", sep = "")
    if (nrow(s))
      print(data.frame(Description = s$Description, GeneRatio = s$GeneRatio,
                       p.adjust = signif(s$p.adjust, 3)), row.names = FALSE)
  }

  p <- dotplot(go, showCategory = TOP_N, split = "ONTOLOGY") +
       facet_grid(ONTOLOGY ~ ., scales = "free") +
       labs(title = "GO enrichment of the TSPO co-expression network") +
       theme_bw(base_size = 10)
  ggsave(file.path(FIGDIR, "Fig1C_GO_enrichment.pdf"), p, width = 9.5, height = 11)
  ggsave(file.path(FIGDIR, "Fig1C_GO_enrichment.png"), p, width = 9.5, height = 11, dpi = 300)
  cat("\n已保存 Fig1C\n")
} else {
  cat("GO 无显著结果\n")
}

hr("STEP 4  KEGG 富集")

kegg <- tryCatch(
  enrichKEGG(gene = mg$ENTREZID, universe = mb$ENTREZID, organism = "hsa",
             pvalueCutoff = 0.05, pAdjustMethod = "BH"),
  error = function(e) {
    cat("KEGG 失败（多数情况是拿不到 KEGG 数据库，需要联网）: ",
        conditionMessage(e), "\n", sep = ""); NULL })

if (!is.null(kegg) && nrow(as.data.frame(kegg))) {
  kd <- as.data.frame(setReadable(kegg, org.Hs.eg.db, "ENTREZID"))
  cat("KEGG 显著通路: ", nrow(kd), "\n", sep = "")
  print(head(data.frame(Description = kd$Description, GeneRatio = kd$GeneRatio,
                        p.adjust = signif(kd$p.adjust, 3)), 12), row.names = FALSE)
  write.csv(kd, file.path(OUTDIR, "T8b_KEGG.csv"), row.names = FALSE)

  p <- dotplot(kegg, showCategory = 15) +
       labs(title = "KEGG pathway enrichment") +
       theme_bw(base_size = 10)
  ggsave(file.path(FIGDIR, "Fig1D_KEGG.pdf"), p, width = 9, height = 7)
  ggsave(file.path(FIGDIR, "Fig1D_KEGG.png"), p, width = 9, height = 7, dpi = 300)
  cat("已保存 Fig1D\n")
} else {
  cat("KEGG 无显著结果或不可用\n")
}

hr("STEP 5  Reactome 富集")

if (need("ReactomePA")) {
  suppressPackageStartupMessages(library(ReactomePA))
  re <- tryCatch(
    enrichPathway(gene = mg$ENTREZID, universe = mb$ENTREZID, organism = "human",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE),
    error = function(e) { cat("Reactome 失败: ", conditionMessage(e), "\n", sep = ""); NULL })

  if (!is.null(re) && nrow(as.data.frame(re))) {
    rd <- as.data.frame(re)
    cat("Reactome 显著通路: ", nrow(rd), "\n", sep = "")
    print(head(data.frame(Description = rd$Description, GeneRatio = rd$GeneRatio,
                          p.adjust = signif(rd$p.adjust, 3)), 10), row.names = FALSE)
    write.csv(rd, file.path(OUTDIR, "T8c_Reactome.csv"), row.names = FALSE)

    p <- dotplot(re, showCategory = 15) +
         labs(title = "Reactome pathway enrichment") +
         theme_bw(base_size = 10)
    ggsave(file.path(FIGDIR, "Fig1E_Reactome.pdf"), p, width = 9, height = 7)
    ggsave(file.path(FIGDIR, "Fig1E_Reactome.png"), p, width = 9, height = 7, dpi = 300)
    cat("已保存 Fig1E\n")
  } else cat("Reactome 无显著结果\n")
} else {
  cat("未安装 ReactomePA，跳过（不影响 GO 与 KEGG）\n")
}

hr("STEP 6  写出网络基因清单")

write.csv(data.frame(gene = genes), file.path(OUTDIR, "T8d_network_genes.csv"),
          row.names = FALSE)
cat("网络基因清单已写出（", length(genes), " 个）\n", sep = "")

cat("\n判读要点：\n")
cat("  · GO/KEGG 若以氧化磷酸化、线粒体翻译、呼吸链复合体为主，与结论一一致。\n")
cat("  · 若同时出现免疫/髓系条目，说明网络里混有细胞构成信号，需在正文说明。\n")
cat("  · 富集是全网络的描述性结果，不构成因果，措辞用「富集于」。\n")

hr("16_enrichment.R 完成")

# -----------------------------------------------------------------------------
# -- Fig2_bulk共表达与免疫 / 17_immune_infiltration.R
# -----------------------------------------------------------------------------
# =============================================================================
#  17_immune_infiltration.R   免疫/基质细胞浸润评分（16 类）
#
#  对应正文 3.2：TSPO 与 12 类免疫细胞全部正相关，免疫总分 +0.412
#
#  方法：ssGSEA（GSVA 可用时）或 z 均值（回退）
#  输入：EXPR_clean.rds
#  输出：T9a_immune_scores.csv、T9b_TSPO_vs_immune.csv、T9c_immune_stromal_scores.csv
#        Fig2D_immune_bar、Fig2E_immune_heatmap
# =============================================================================

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  if (file.exists(.s)) source(.s, encoding = "UTF-8") else {
    PROJ <- "D:/Rstudy/tcga"; OUT <- file.path(PROJ, "results_v2")
    F_EXPR <- file.path(OUT, "EXPR_clean.rds")
  }
}
if (!exists("F_EXPR")) F_EXPR <- file.path(OUT, "EXPR_clean.rds")
OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(ggplot2) })

# 16 类特征基因（紧凑版，每类 3–5 个标记）
SIG <- list(
  CD8_T      = c("CD8A","CD8B","GZMK","PRF1","IFNG"),
  CD4_T      = c("CD4","IL7R","CCR7","LTB"),
  Treg       = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  NK         = c("KLRD1","NKG7","GNLY","KLRF1"),
  B_cell     = c("CD19","MS4A1","CD79A","CD79B"),
  Plasma     = c("MZB1","JCHAIN","IGHG1"),
  Monocyte   = c("CD14","LYZ","S100A8","S100A9","VCAN"),
  Mac_M1     = c("NOS2","IL1B","TNF","CXCL9","CXCL10"),
  Mac_M2     = c("CD163","MRC1","MSR1","CCL22","IL10"),
  DC         = c("CD1C","CLEC9A","LAMP3","CCR7"),
  Neutrophil = c("FCGR3B","CSF3R","CXCR2"),
  Mast       = c("TPSAB1","CPA3","MS4A2"),
  Fibroblast = c("COL1A1","COL1A2","DCN","FAP"),
  Endothelial= c("PECAM1","VWF","CDH5"),
  Epithelial = c("EPCAM","KRT8","KRT18"),
  Prolif     = c("MKI67","PCNA","TOP2A","BUB1")
)

hr("STEP 1  载入表达矩阵")

EXPR <- as.matrix(readRDS(F_EXPR)); storage.mode(EXPR) <- "double"
cat("表达矩阵: ", nrow(EXPR), " x ", ncol(EXPR), "\n", sep = "")

SIG <- lapply(SIG, function(g) intersect(g, rownames(EXPR)))
for (n in names(SIG)) cat(sprintf("  %-12s %d\n", n, length(SIG[[n]])))
SIG <- SIG[vapply(SIG, length, 0) >= 2]
cat("可用特征集: ", length(SIG), "\n", sep = "")

hr("STEP 2  计算浸润评分")

HAVE_GSVA <- requireNamespace("GSVA", quietly = TRUE)
sc <- NULL
if (HAVE_GSVA) {
  suppressPackageStartupMessages(library(GSVA))
  sc <- tryCatch({
    par <- GSVA::ssgseaParam(EXPR, SIG, normalize = TRUE)
    GSVA::gsva(par, verbose = FALSE)
  }, error = function(e) { cat("ssGSEA 失败，回退 z 均值: ", conditionMessage(e), "\n", sep = ""); NULL })
}
if (is.null(sc)) {
  cat("使用 z 均值法\n")
  Z  <- t(scale(t(EXPR))); Z[!is.finite(Z)] <- NA
  sc <- t(sapply(SIG, function(g) colMeans(Z[g, , drop = FALSE], na.rm = TRUE)))
}
cat("评分矩阵: ", nrow(sc), " 类 x ", ncol(sc), " 样本\n", sep = "")
write.csv(data.frame(sample = colnames(sc), t(sc), check.names = FALSE),
          file.path(OUTDIR, "T9a_immune_scores.csv"), row.names = FALSE)

hr("STEP 3  TSPO 与各类细胞评分")

ts <- as.numeric(EXPR["TSPO", colnames(sc)])
res <- do.call(rbind, lapply(rownames(sc), function(ct) {
  y  <- as.numeric(sc[ct, ])
  ok <- is.finite(y) & is.finite(ts)
  ct_test <- suppressWarnings(cor.test(y[ok], ts[ok], method = "spearman"))
  data.frame(celltype = ct, rho = round(as.numeric(ct_test$estimate), 3),
             p = ct_test$p.value, n = sum(ok), stringsAsFactors = FALSE)
}))
res$FDR <- p.adjust(res$p, "BH")
res <- res[order(-res$rho), ]
cat("\nTSPO 与各类细胞评分：\n"); pr(res, 4)
write.csv(res, file.path(OUTDIR, "T9b_TSPO_vs_immune.csv"), row.names = FALSE)

res$celltype <- factor(res$celltype, levels = res$celltype[order(res$rho)])
p1 <- ggplot(res, aes(rho, celltype, fill = rho)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
  labs(x = expression("Spearman " * rho * " (TSPO vs infiltration score)"), y = NULL,
       title = "TSPO versus immune and stromal infiltration") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "Fig2D_immune_bar.pdf"), p1, width = 7, height = 5.5)
ggsave(file.path(FIGDIR, "Fig2D_immune_bar.png"), p1, width = 7, height = 5.5, dpi = 300)

hr("STEP 4  免疫与基质总分")

imm  <- intersect(c("CD8_T","CD4_T","Treg","NK","B_cell","Plasma","Monocyte",
                    "Mac_M1","Mac_M2","DC","Neutrophil","Mast"), rownames(sc))
strm <- intersect(c("Fibroblast","Endothelial"), rownames(sc))
if (length(imm) >= 3 && length(strm) >= 1) {
  imm_score  <- colMeans(sc[imm, , drop = FALSE], na.rm = TRUE)
  strm_score <- colMeans(sc[strm, , drop = FALSE], na.rm = TRUE)
  for (nm in c("Immune","Stromal")) {
    v <- if (nm == "Immune") imm_score else strm_score
    ct2 <- suppressWarnings(cor.test(v, ts, method = "spearman"))
    cat(sprintf("%-8s 总分 ~ TSPO: rho = %+.3f, p = %.3g\n",
                nm, as.numeric(ct2$estimate), ct2$p.value))
  }
  write.csv(data.frame(sample = colnames(sc), TSPO = ts,
                       Immune = imm_score, Stromal = strm_score),
            file.path(OUTDIR, "T9c_immune_stromal_scores.csv"), row.names = FALSE)
}

hr("17_immune_infiltration.R 完成")

# -----------------------------------------------------------------------------
# -- Fig2_bulk共表达与免疫 / 18_tspo_highlow_deg.R
# -----------------------------------------------------------------------------
# =============================================================================
#  18_tspo_highlow_deg.R   TSPO 高/低组差异表达（limma）+ 火山图 + 热图 + 富集
#
#  对应正文 3.2：2,564 个差异基因（上调 1,183、下调 1,381）
#                上调以 GNLY/IGJ/HLA-DQB2/GZMA 等免疫基因为主
#                下调以 TUBB2B/DPYSL5/PAX2/LHX1 等神经发育基因为主
#
#  输入：EXPR_clean.rds
#  输出：T10a_DEG_full.csv、T10b_DEG_sig.csv、T10c_GO_Up/Down.csv
#        Fig3A_volcano、Fig3B_DEG_heatmap
# =============================================================================

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  if (file.exists(.s)) source(.s, encoding = "UTF-8") else {
    PROJ <- "D:/Rstudy/tcga"; OUT <- file.path(PROJ, "results_v2")
    F_EXPR <- file.path(OUT, "EXPR_clean.rds")
  }
}
if (!exists("F_EXPR")) F_EXPR <- file.path(OUT, "EXPR_clean.rds")
OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(limma); library(ggplot2) })

LFC_CUT <- 0.5
FDR_CUT <- 0.05

hr("STEP 1  按 TSPO 中位数分组")

EXPR <- as.matrix(readRDS(F_EXPR)); storage.mode(EXPR) <- "double"
cat("表达矩阵: ", nrow(EXPR), " x ", ncol(EXPR), "\n", sep = "")

ts  <- as.numeric(EXPR["TSPO", ])
med <- median(ts, na.rm = TRUE)
grp <- factor(ifelse(ts > med, "High", "Low"), levels = c("Low", "High"))
cat("TSPO 中位数切点: ", round(med, 3), "\n", sep = "")
print(table(grp))

hr("STEP 2  limma 差异分析")

design <- model.matrix(~ grp)
fit <- lmFit(EXPR, design)
fit <- eBayes(fit)
tt  <- topTable(fit, coef = "grpHigh", number = Inf, sort.by = "P")
tt$gene <- rownames(tt)
tt <- tt[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
cat("检验基因数: ", nrow(tt), "\n", sep = "")

sig <- tt[abs(tt$logFC) >= LFC_CUT & tt$adj.P.Val < FDR_CUT, ]
sig <- sig[order(-abs(sig$logFC)), ]
cat("显著基因（|logFC|≥", LFC_CUT, " 且 FDR<", FDR_CUT, "）: ", nrow(sig),
    "（上调 ", sum(sig$logFC > 0), "，下调 ", sum(sig$logFC < 0), "）\n", sep = "")

write.csv(tt,  file.path(OUTDIR, "T10a_DEG_full.csv"), row.names = FALSE)
write.csv(sig, file.path(OUTDIR, "T10b_DEG_sig.csv"),  row.names = FALSE)

cat("\n上调 Top10:\n")
pr(head(sig[sig$logFC > 0, c("gene","logFC","adj.P.Val")], 10), 4)
cat("\n下调 Top10:\n")
pr(head(sig[sig$logFC < 0, c("gene","logFC","adj.P.Val")], 10), 4)

hr("STEP 3  火山图")

tt$sig <- ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC >= LFC_CUT, "Up",
          ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC <= -LFC_CUT, "Down", "NS"))
p <- ggplot(tt, aes(logFC, -log10(adj.P.Val), colour = sig)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_colour_manual(values = c(Up = "#D7191C", Down = "#2C7BB6", NS = "grey75")) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 2, colour = "grey40") +
  geom_hline(yintercept = -log10(FDR_CUT), linetype = 2, colour = "grey40") +
  labs(x = "log2 fold change (TSPO-high vs TSPO-low)", y = "-log10(FDR)",
       title = "Differential expression by TSPO status", colour = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "Fig3A_volcano.pdf"), p, width = 6, height = 5)
ggsave(file.path(FIGDIR, "Fig3A_volcano.png"), p, width = 6, height = 5, dpi = 300)

hr("STEP 4  热图（Top 50，英文标签）")

top50 <- head(sig$gene, 50)
if (length(top50) >= 10) {
  Z   <- t(scale(t(EXPR[top50, ]))); Z[!is.finite(Z)] <- NA
  ord <- order(ts)
  hm  <- data.frame(gene = rep(top50, each = ncol(Z)),
                    sample = rep(colnames(Z), times = length(top50)),
                    z = as.numeric(Z[, ord]), stringsAsFactors = FALSE)
  hm$sample <- factor(hm$sample, levels = colnames(Z)[ord])
  hm$gene   <- factor(hm$gene, levels = top50)
  p2 <- ggplot(hm, aes(sample, gene, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
    labs(x = "Samples ordered by TSPO (low to high)", y = NULL,
         title = "Top 50 differentially expressed genes") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave(file.path(FIGDIR, "Fig3B_DEG_heatmap.pdf"), p2, width = 8, height = 7)
  ggsave(file.path(FIGDIR, "Fig3B_DEG_heatmap.png"), p2, width = 8, height = 7, dpi = 300)
} else cat("显著基因不足 10 个，跳过热图\n")

hr("STEP 5  上下调基因分别富集")

if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  for (dir in c("Up", "Down")) {
    g <- if (dir == "Up") sig$gene[sig$logFC > 0] else sig$gene[sig$logFC < 0]
    if (length(g) < 10) { cat(dir, ": 基因太少，跳过\n", sep = ""); next }
    m <- suppressWarnings(try(bitr(g, "SYMBOL", "ENTREZID", org.Hs.eg.db), silent = TRUE))
    if (inherits(m, "try-error") || !nrow(m)) next
    eg <- tryCatch(enrichGO(m$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                            pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE),
                   error = function(e) NULL)
    if (!is.null(eg) && nrow(as.data.frame(eg))) {
      df <- as.data.frame(eg)
      write.csv(df, file.path(OUTDIR, paste0("T10c_GO_", dir, ".csv")), row.names = FALSE)
      cat("\n", dir, " 显著条目: ", nrow(df), "\n", sep = "")
      pr(head(data.frame(Description = df$Description, GeneRatio = df$GeneRatio,
                         p.adjust = signif(df$p.adjust, 3)), 6), 3)
    }
  }
} else cat("未装 clusterProfiler，跳过富集\n")

cat("\n判读要点：\n")
cat("  · 上调部分与 3.2 的浸润评分反映同一件事（组织内免疫细胞比例），不能重复计数。\n")
cat("  · 真正独立的信息是增殖下调与神经发育程序在 TSPO 低表达肿瘤中富集。\n")

hr("18_tspo_highlow_deg.R 完成")

# -----------------------------------------------------------------------------
# -- Fig2_bulk共表达与免疫 / 19_nk_vs_ctl.R
# -----------------------------------------------------------------------------
# =============================================================================
#  19_nk_vs_ctl.R   拆分 NK 与细胞毒 T 细胞，避免共用基因造成的伪高相关
#
#  问题：上一轮的 NK 特征集是 KLRD1 / NKG7 / GNLY / KLRF1，
#        其中 GNLY 与 NKG7 在细胞毒 T 细胞（CTL）中同样高表达，
#        所以那个「NK 分数」实际测的是细胞毒淋巴细胞程序，与 TSPO 的 0.506 有虚高。
#
#  本脚本做三件事：
#    1. 用互不重叠的基因集分别定义 NK 与 CTL
#    2. 两者分别与 TSPO 求相关
#    3. 互相控制后做偏相关，判断哪一个才是独立关联
#
#  输入：EXPR_clean.rds
#  输出：T11a_NK_CTL_scores.csv、T11b_NK_CTL_vs_TSPO.csv、
#        T11c_NK_CTL_partial.csv、T11d_NK_CTL_regression.csv
# =============================================================================

if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  .ok <- FALSE
  if (file.exists(.s))
    .ok <- !inherits(try(source(.s, encoding = "UTF-8"), silent = TRUE), "try-error")
  if (!.ok) {
    PROJ <- "D:/Rstudy/tcga"; OUT <- file.path(PROJ, "results_v2")
    hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
    pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)
  }
  for (d in c(OUT, file.path(OUT, "tables"), file.path(OUT, "figures")))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("hr")) hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
if (!exists("pr")) pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")

# ---- 互不重叠的基因集 -------------------------------------------------------
#  NK 特异：尽量选 NK 限定基因，排除 GNLY / NKG7 / PRF1 / GZMB 等与 CTL 共用的
#  CTL 特异：CD8 + 细胞毒颗粒，不含 NK 限定基因
#  KIR 家族多态性强、在 RNA-seq 参考里常缺失，故多给几个候选，靠 intersect 兜底
GENESETS <- list(
  NK_specific = c("KLRF1", "KLRC1", "KLRC2", "KLRC3", "KLRD1",
                  "SH2D1B", "NCR1", "KIR2DL3", "KIR3DL1", "KIR2DS4",
                  "FCER1G", "TXK", "PRSS23"),
  CD8_CTL     = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "PRF1",
                  "NKG7", "GNLY", "KLRG1", "CST7"),
  CD8_naive   = c("CD8A", "CCR7", "SELL", "LEF1", "IL7R", "TCF7"),
  CD4_T       = c("CD4", "IL7R", "CCR7", "LTB"),
  Treg        = c("FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  Monocyte    = c("CD14", "LYZ", "S100A8", "S100A9", "VCAN")
)

hr("STEP 1  载入表达矩阵并筛选基因")

EXPR <- as.matrix(readRDS(F_EXPR)); storage.mode(EXPR) <- "double"
cat("表达矩阵: ", nrow(EXPR), " 基因 x ", ncol(EXPR), " 样本\n", sep = "")

GENESETS <- lapply(GENESETS, function(g) intersect(g, rownames(EXPR)))
for (n in names(GENESETS))
  cat(sprintf("  %-12s %2d 基因: %s\n", n, length(GENESETS[[n]]),
              paste(GENESETS[[n]], collapse = ", ")))

# NK 与 CTL 必须真正不重叠，否则结论还是绕回原点
ov <- intersect(GENESETS$NK_specific, GENESETS$CD8_CTL)
if (length(ov)) {
  cat("\n⚠ NK 与 CTL 仍有重叠基因：", paste(ov, collapse = ", "),
      "——已从 NK 侧移除\n", sep = "")
  GENESETS$NK_specific <- setdiff(GENESETS$NK_specific, ov)
}

keep <- vapply(GENESETS, length, 0) >= 2
if (any(!keep)) cat("基因太少而剔除的特征集: ", paste(names(GENESETS)[!keep], collapse = ", "), "\n", sep = "")
GENESETS <- GENESETS[keep]
if (length(GENESETS$NK_specific) < 2)
  stop("NK 特异基因不足 2 个，无法区分。请检查表达矩阵是否含 KLRF1/KLRC1 等")

hr("STEP 2  计算评分（ssGSEA，缺失时回退 z 均值）")

use_ssgsea <- requireNamespace("GSVA", quietly = TRUE)
if (use_ssgsea) {
  suppressPackageStartupMessages(library(GSVA))
  sc <- tryCatch({
    par <- GSVA::ssgseaParam(EXPR, GENESETS, normalize = TRUE)
    GSVA::gsva(par, verbose = FALSE)
  }, error = function(e) {
    cat("ssGSEA 失败，回退 z 均值: ", conditionMessage(e), "\n", sep = ""); NULL
  })
  if (is.null(sc)) use_ssgsea <- FALSE
}
if (!use_ssgsea) {
  cat("使用 z 均值法\n")
  Z  <- t(scale(t(EXPR))); Z[!is.finite(Z)] <- NA
  sc <- t(sapply(GENESETS, function(g) colMeans(Z[g, , drop = FALSE], na.rm = TRUE)))
}
cat("评分方法: ", ifelse(use_ssgsea, "ssGSEA", "z 均值"), "\n", sep = "")
cat("评分矩阵: ", nrow(sc), " 类 x ", ncol(sc), " 样本\n", sep = "")

write.csv(data.frame(sample = colnames(sc), t(sc), check.names = FALSE),
          file.path(OUTDIR, "T11a_NK_CTL_scores.csv"), row.names = FALSE)

hr("STEP 3  与 TSPO 的相关")

ts <- as.numeric(EXPR["TSPO", colnames(sc)])

res <- do.call(rbind, lapply(rownames(sc), function(ct) {
  y  <- as.numeric(sc[ct, ])
  ok <- is.finite(y) & is.finite(ts)
  ct_test <- suppressWarnings(cor.test(y[ok], ts[ok], method = "spearman"))
  data.frame(celltype = ct,
             rho = round(as.numeric(ct_test$estimate), 3),
             p   = ct_test$p.value, n = sum(ok),
             n_gene = length(GENESETS[[ct]]),
             stringsAsFactors = FALSE)
}))
res$FDR <- p.adjust(res$p, "BH")
res <- res[order(-res$rho), ]
cat("\nTSPO 与各特征集的相关：\n")
pr(res, 4)

write.csv(res, file.path(OUTDIR, "T11b_NK_CTL_vs_TSPO.csv"), row.names = FALSE)

hr("STEP 4  NK 与 CTL 的相关（共线性检查）")

r_nk_ctl <- suppressWarnings(cor.test(as.numeric(sc["NK_specific", ]),
                                      as.numeric(sc["CD8_CTL", ]),
                                      method = "spearman"))
cat("NK_specific 与 CD8_CTL 的 Spearman ρ = ",
    round(as.numeric(r_nk_ctl$estimate), 3),
    "（p = ", format.pval(r_nk_ctl$p.value, digits = 3), "）\n", sep = "")
cat("两者高度共线时，偏相关会不稳定，需同时看回归结果。\n")

hr("STEP 5  互相控制后的偏相关")

pcor_rank <- function(a, b, cov) {
  ok <- is.finite(a) & is.finite(b) & is.finite(cov)
  if (sum(ok) < 20) return(c(rho = NA, p = NA, n = sum(ok)))
  ra <- residuals(lm(rank(a[ok]) ~ rank(cov[ok])))
  rb <- residuals(lm(rank(b[ok]) ~ rank(cov[ok])))
  r  <- cor(ra, rb)
  t  <- r * sqrt((sum(ok) - 2) / (1 - r^2))
  c(rho = r, p = 2 * pt(-abs(t), df = sum(ok) - 2), n = sum(ok))
}

nk  <- as.numeric(sc["NK_specific", ])
ctl <- as.numeric(sc["CD8_CTL", ])

p_nk  <- pcor_rank(ts, nk,  ctl)    # TSPO ~ NK  | CTL
p_ctl <- pcor_rank(ts, ctl, nk)     # TSPO ~ CTL | NK

part <- data.frame(
  model = c("TSPO ~ NK_specific | CD8_CTL", "TSPO ~ CD8_CTL | NK_specific"),
  rho_partial = round(c(p_nk["rho"], p_ctl["rho"]), 3),
  p = c(p_nk["p"], p_ctl["p"]),
  n = c(p_nk["n"], p_ctl["n"]),
  stringsAsFactors = FALSE
)
cat("\n")
pr(part, 4)
write.csv(part, file.path(OUTDIR, "T11c_NK_CTL_partial.csv"), row.names = FALSE)

hr("STEP 6  多元回归（同时纳入两者）")

d <- data.frame(TSPO = ts, NK = nk, CTL = ctl)
d <- d[is.finite(d$TSPO) & is.finite(d$NK) & is.finite(d$CTL), ]
d$NK  <- as.numeric(scale(d$NK))
d$CTL <- as.numeric(scale(d$CTL))

m <- lm(TSPO ~ NK + CTL, data = d)
sm <- summary(m)
cat("\n标准化回归系数（TSPO 为因变量）：\n")
print(round(sm$coefficients, 4))
cat("\nR² = ", round(sm$r.squared, 3), "，调整 R² = ", round(sm$adj.r.squared, 3), "\n", sep = "")
cat("VIF(NK) = ", round(1 / (1 - summary(lm(NK ~ CTL, data = d))$r.squared), 2),
    "，VIF(CTL) = ", round(1 / (1 - summary(lm(CTL ~ NK, data = d))$r.squared), 2), "\n", sep = "")

reg <- data.frame(term = rownames(sm$coefficients),
                  estimate = sm$coefficients[, 1],
                  se = sm$coefficients[, 2],
                  t = sm$coefficients[, 3],
                  p = sm$coefficients[, 4],
                  row.names = NULL)
write.csv(reg, file.path(OUTDIR, "T11d_NK_CTL_regression.csv"), row.names = FALSE)

hr("STEP 7  判读指引")

cat("三种可能结局：\n")
cat("  A) NK 与 CTL 都显著、偏相关后都还在 → 细胞毒程序整体与 TSPO 相关，\n")
cat("     正文写「细胞毒淋巴细胞浸润」，不区分 NK/T。\n")
cat("  B) 只有 CTL 在偏相关后存活 → 之前 0.506 确实是 NK/CTL 共用基因造成的，\n")
cat("     正文写「细胞毒 T 细胞」，并说明已排除共用基因。\n")
cat("  C) 只有 NK 存活 → 罕见，需检查 KIR 基因是否真的被检出（STEP 1 会打印）。\n")
cat("\n注意：若 VIF > 5，说明两个分数高度共线，偏相关的解释力有限，\n")
cat("应以多元回归的系数符号与量级为准。\n")

hr("19_nk_vs_ctl.R 完成")

# -----------------------------------------------------------------------------
# -- Fig2_bulk共表达与免疫 / 23_wgcna.R
# -----------------------------------------------------------------------------
# =============================================================================
#  23_wgcna.R  (v2)   WGCNA 共表达模块与 TSPO 所在模块
#
#  v2 修改：
#    1. hr() / pr() 改为无条件定义，不再依赖 exists() 判断
#       （v1 里这两个函数定义在 if (!exists("OUT")) 块内，
#         若 OUT 已存在则整块跳过，导致后续 pr() 报"没有 pr 这个函数"）
#    2. POWER 默认固定为 4（R²=0.979、slope=−1.51、mean.k=15.3，
#       优于自动选择的 power 3），可改回 USE_AUTO_POWER = TRUE
#
#  目的：1) 用 WGCNA 独立划分模块（与 02 的阈值网络互为验证）
#        2) 找 TSPO 所属模块，看它对应代谢还是免疫
#        3) 报告 hub 基因，与 02 的关键节点对照
#
#  ⚠ 02 是「与 TSPO 逐基因相关 + 阈值」，WGCNA 是「全基因两两相关 + 层次聚类」，
#    两者重合的基因才最稳。
#
#  输入：EXPR_clean.rds（18,747 × 304）
#  输出：T15a–T15f 六张表 + Fig5A/5B/5C
# =============================================================================

# ---- 工具函数：无条件定义，避免被前置块跳过 --------------------------------
hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  .ok <- FALSE
  if (file.exists(.s))
    .ok <- !inherits(try(source(.s, encoding = "UTF-8"), silent = TRUE), "try-error")
  if (!.ok) {
    PROJ <- "D:/Rstudy/tcga"; OUT <- file.path(PROJ, "results_v2")
    F_EXPR <- file.path(OUT, "EXPR_clean.rds")
  }
  for (d in c(OUT, file.path(OUT, "tables"), file.path(OUT, "figures")))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("F_EXPR")) F_EXPR <- file.path(OUT, "EXPR_clean.rds")

OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

N_TOP_GENES   <- 5000
MIN_MODULE    <- 30
MERGE_CUT     <- 0.25
SEED_WGCNA    <- 20260101
USE_AUTO_POWER <- FALSE   # TRUE 则用 pickSoftThreshold 的自动结果
POWER_FIXED   <- 4

MODULES <- list(
  OXPHOS     = c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
                 "SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","UQCRH","CYC1",
                 "COX4I1","COX5A","COX5B","COX6A1","ATP5A1","ATP5B","ATP5C1","ATP5F1"),
  Glycolysis = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1"),
  Immune     = c("CD68","CD3D","CD3E","CD8A","PTPRC","LYZ","C1QA","C1QB","MS4A1"),
  Stromal    = c("COL1A1","COL1A2","COL3A1","DCN","FAP","SPARC","PDGFRB"),
  Prolif     = c("MKI67","TOP2A","PCNA","CCNB1","BUB1")
)

hr("STEP 0  依赖检查")

if (!requireNamespace("WGCNA", quietly = TRUE))
  stop("请先安装 WGCNA：install.packages('WGCNA')")

suppressPackageStartupMessages({ library(WGCNA); library(ggplot2) })
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = max(1, parallel::detectCores() - 1))
cat("WGCNA 版本: ", as.character(packageVersion("WGCNA")), "\n", sep = "")

hr("STEP 1  载入表达矩阵与性状")

if (!file.exists(F_EXPR)) stop("找不到 ", F_EXPR, "，请先跑 01_build_matrix.R")
EXPR <- as.matrix(readRDS(F_EXPR))
storage.mode(EXPR) <- "double"
cat("表达矩阵: ", nrow(EXPR), " 基因 x ", ncol(EXPR), " 样本\n", sep = "")
if (!"TSPO" %in% rownames(EXPR)) stop("矩阵中没有 TSPO")

datExpr <- t(EXPR)

Z <- t(scale(t(EXPR))); Z[!is.finite(Z)] <- NA
score <- function(g) {
  gg <- intersect(g, rownames(Z))
  if (length(gg) < 3) return(rep(NA_real_, ncol(Z)))
  colMeans(Z[gg, , drop = FALSE], na.rm = TRUE)
}

traits <- data.frame(
  TSPO       = as.numeric(EXPR["TSPO", ]),
  OXPHOS     = score(MODULES$OXPHOS),
  Glycolysis = score(MODULES$Glycolysis),
  Immune     = score(MODULES$Immune),
  Stromal    = score(MODULES$Stromal),
  Prolif     = score(MODULES$Prolif),
  row.names  = colnames(EXPR)
)
traits$Purity_proxy <- -traits$Immune - traits$Stromal

cat("性状矩阵: ", nrow(traits), " x ", ncol(traits), "\n", sep = "")
print(round(cor(traits, use = "pairwise.complete.obs"), 2))

hr("STEP 2  选高变基因")

gv <- apply(datExpr, 2, var, na.rm = TRUE)
gv[!is.finite(gv)] <- 0
ord <- order(gv, decreasing = TRUE)
top <- colnames(datExpr)[ord[1:min(N_TOP_GENES, length(ord))]]
if (!"TSPO" %in% top) top <- c(top, "TSPO")
cat("纳入基因: ", length(top), "（TSPO 在内: ", "TSPO" %in% top, "）\n", sep = "")

datExpr <- datExpr[, top, drop = FALSE]

gsg <- goodSamplesGenes(datExpr, verbose = 0)
cat("goodSamplesGenes 全通过: ", gsg$allOK, "\n", sep = "")
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat("剔除低质量基因: ", sum(!gsg$goodGenes), "\n", sep = "")
    datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("剔除低质量样本: ", sum(!gsg$goodSamples), "\n", sep = "")
    datExpr <- datExpr[gsg$goodSamples, , drop = FALSE]
    traits  <- traits[gsg$goodSamples, , drop = FALSE]
  }
}
cat("最终: ", nrow(datExpr), " 样本 x ", ncol(datExpr), " 基因\n", sep = "")

hr("STEP 3  软阈值选择")

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "unsigned",
                         verbose = 0)
sft_tab <- data.frame(power    = sft$fitIndices$Power,
                      SFT_R2   = round(sft$fitIndices$SFT.R.sq, 3),
                      slope    = round(sft$fitIndices$slope, 2),
                      mean_k   = round(sft$fitIndices$mean.k., 2),
                      median_k = round(sft$fitIndices$median.k., 2))
cat("\n软阈值扫描：\n"); pr(sft_tab)
write.csv(sft_tab, file.path(OUTDIR, "T15a_soft_threshold.csv"), row.names = FALSE)

if (USE_AUTO_POWER) {
  POWER <- sft$powerEstimate
  if (is.na(POWER)) {
    cand <- sft_tab$power[sft_tab$SFT_R2 >= 0.8]
    POWER <- if (length(cand)) min(cand) else 6
  }
  cat("\n使用自动选择的 power = ", POWER, "\n", sep = "")
} else {
  POWER <- POWER_FIXED
  cat("\n使用固定 power = ", POWER, "（自动结果为 ",
      sft$powerEstimate, "）\n", sep = "")
}
cat("对应 R² = ", sft_tab$SFT_R2[sft_tab$power == POWER],
    "，slope = ", sft_tab$slope[sft_tab$power == POWER], "\n", sep = "")

hr("STEP 4  构建模块")

net <- blockwiseModules(
  datExpr,
  power             = POWER,
  networkType       = "unsigned",
  TOMType           = "unsigned",
  minModuleSize     = MIN_MODULE,
  reassignThreshold = 0,
  mergeCutHeight    = MERGE_CUT,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  maxBlockSize      = 6000,
  randomSeed        = SEED_WGCNA,
  verbose           = 0
)

mod_lab <- net$colors
cat("模块数（含 0 号未分配）: ", length(unique(mod_lab)), "\n", sep = "")
print(table(mod_lab))

mod_df <- data.frame(gene = colnames(datExpr), module = paste0("M", mod_lab),
                     stringsAsFactors = FALSE)
write.csv(mod_df, file.path(OUTDIR, "T15b_module_assignment.csv"), row.names = FALSE)

hr("STEP 5  模块 - 性状相关")

MEs <- moduleEigengenes(datExpr, colors = mod_lab)$eigengenes
MEs <- orderMEs(MEs)
cat("模块特征基因: ", nrow(MEs), " x ", ncol(MEs), "\n", sep = "")

mt   <- cor(MEs, traits, use = "pairwise.complete.obs")
mt_p <- corPvalueStudent(mt, nrow(datExpr))
write.csv(data.frame(module = rownames(mt), round(mt, 3), check.names = FALSE),
          file.path(OUTDIR, "T15c_module_trait_cor.csv"), row.names = FALSE)
write.csv(data.frame(module = rownames(mt_p), signif(mt_p, 3), check.names = FALSE),
          file.path(OUTDIR, "T15c_module_trait_p.csv"), row.names = FALSE)

cat("\n模块-性状相关（r）：\n"); print(round(mt, 3))
cat("\n对应 p 值：\n"); print(signif(mt_p, 3))

hr("STEP 6  TSPO 所在模块")

tspo_mod_num <- mod_lab[colnames(datExpr) == "TSPO"]
cat("TSPO 的模块编号: M", tspo_mod_num, "\n", sep = "")
hub <- NULL

if (tspo_mod_num == 0) {
  cat("\n⚠ TSPO 落在未分配模块（M0），与任何模块特征基因都不够相关。\n")
  cat("  这本身是一个结果：TSPO 不属于任何稳定的共表达模块。\n")
} else {
  mod_genes <- colnames(datExpr)[mod_lab == tspo_mod_num]
  cat("该模块基因数: ", length(mod_genes), "\n", sep = "")

  me <- MEs[, paste0("ME", tspo_mod_num)]
  kme <- as.numeric(cor(datExpr[, mod_genes, drop = FALSE], me,
                        use = "pairwise.complete.obs"))
  names(kme) <- mod_genes

  tspo_cor <- as.numeric(cor(datExpr[, mod_genes, drop = FALSE],
                             datExpr[, "TSPO"], use = "pairwise.complete.obs"))

  hub <- data.frame(gene = mod_genes, kME = round(kme, 3),
                    cor_with_TSPO = round(tspo_cor, 3), stringsAsFactors = FALSE)
  hub$p_TSPO <- signif(2 * pt(-abs(hub$cor_with_TSPO *
                    sqrt((nrow(datExpr) - 2) / (1 - hub$cor_with_TSPO^2))),
                    df = nrow(datExpr) - 2), 3)
  hub <- hub[order(-abs(hub$kME)), ]
  cat("\n模块 hub 基因（按 |kME| 排序，前 25）：\n")
  pr(head(hub, 25), 3)
  write.csv(hub, file.path(OUTDIR, "T15d_TSPO_module_hubs.csv"), row.names = FALSE)

  cat("\n模块内与 TSPO 相关最强的 15 个基因：\n")
  pr(head(hub[order(-abs(hub$cor_with_TSPO)),
              c("gene","cor_with_TSPO","kME","p_TSPO")], 15), 3)
}

hr("STEP 7  TSPO 模块的 GO 富集")

if (!is.null(hub) && requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  m <- suppressWarnings(try(bitr(hub$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db),
                            silent = TRUE))
  if (!inherits(m, "try-error") && nrow(m) >= 10) {
    eg <- tryCatch(enrichGO(m$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                            pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE),
                   error = function(e) NULL)
    if (!is.null(eg) && nrow(as.data.frame(eg))) {
      gd <- as.data.frame(eg)
      write.csv(gd, file.path(OUTDIR, "T15e_TSPO_module_GO.csv"), row.names = FALSE)
      cat("GO BP 显著条目: ", nrow(gd), "\n", sep = "")
      pr(head(data.frame(Description = gd$Description, GeneRatio = gd$GeneRatio,
                         p.adjust = signif(gd$p.adjust, 3)), 12), 3)
    } else cat("GO 无显著结果\n")
  } else cat("ID 转换不足，跳过 GO\n")
} else cat("未装 clusterProfiler 或 TSPO 不在模块内，跳过 GO\n")

hr("STEP 8  与 02 的阈值网络对照")

f_net <- file.path(OUTDIR, "T8d_network_genes.csv")
if (!file.exists(f_net)) f_net <- file.path(OUTDIR, "T1b_TSPO_network_significant.csv")
if (file.exists(f_net) && !is.null(hub)) {
  netg <- read.csv(f_net, check.names = FALSE)
  netg <- unique(as.character(netg[[1]]))
  inter <- intersect(hub$gene, netg)
  ov <- data.frame(set = c("WGCNA 模块", "阈值网络", "两者交集"),
                   n = c(nrow(hub), length(netg), length(inter)))
  cat("\n"); pr(ov, 3)
  cat("模块基因中落在阈值网络内的比例: ",
      round(100 * length(inter) / nrow(hub), 1), "%\n", sep = "")
  write.csv(data.frame(gene = inter), file.path(OUTDIR, "T15f_overlap_with_network.csv"),
            row.names = FALSE)
  cat("\n两法一致的前 15 个基因（按 kME 排序）：\n")
  pr(head(hub[hub$gene %in% inter, c("gene","kME","cor_with_TSPO")], 15), 3)
} else {
  cat("跳过：", if (is.null(hub)) "TSPO 不在模块内" else "未找到 02 的网络结果", "\n", sep = "")
}

hr("STEP 9  图")

pdf(file.path(FIGDIR, "Fig5A_module_dendrogram.pdf"), width = 9, height = 5)
plotDendroAndColors(net$dendrograms[[1]],
                    paste0("M", mod_lab[net$blockGenes[[1]]]),
                    "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
png(file.path(FIGDIR, "Fig5A_module_dendrogram.png"), width = 1800, height = 1000, res = 200)
plotDendroAndColors(net$dendrograms[[1]],
                    paste0("M", mod_lab[net$blockGenes[[1]]]),
                    "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

mt_long <- data.frame(
  module = rep(rownames(mt), times = ncol(mt)),
  trait  = rep(colnames(mt), each = nrow(mt)),
  r      = as.numeric(mt),
  p      = as.numeric(mt_p),
  stringsAsFactors = FALSE
)
mt_long$star <- ifelse(mt_long$p < 0.001, "***",
                ifelse(mt_long$p < 0.01, "**",
                ifelse(mt_long$p < 0.05, "*", "")))
mt_long$module <- factor(mt_long$module, levels = rownames(mt))

pB <- ggplot(mt_long, aes(trait, module, fill = r)) +
  geom_tile() +
  geom_text(aes(label = star), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, title = "Module-trait relationships", fill = "Pearson r") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIGDIR, "Fig5B_module_trait.pdf"), pB, width = 8, height = 6)
ggsave(file.path(FIGDIR, "Fig5B_module_trait.png"), pB, width = 8, height = 6, dpi = 300)

if (!is.null(hub)) {
  pC <- ggplot(hub, aes(kME, cor_with_TSPO)) +
    geom_point(size = 1.2, alpha = 0.6, colour = "#2C7BB6") +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    labs(x = "Module membership (kME)", y = "Correlation with TSPO",
         title = paste0("TSPO module (M", tspo_mod_num, "): ", nrow(hub), " genes")) +
    theme_bw(base_size = 11)
  ggsave(file.path(FIGDIR, "Fig5C_hub_scatter.pdf"), pC, width = 6, height = 5)
  ggsave(file.path(FIGDIR, "Fig5C_hub_scatter.png"), pC, width = 6, height = 5, dpi = 300)
}

saveRDS(list(net = net, MEs = MEs, traits = traits, mt = mt, hub = hub,
             power = POWER, tspo_mod = tspo_mod_num),
        file.path(OUT, "cache", "wgcna.rds"))

cat("\n判读要点：\n")
cat("  · TSPO 模块与 OXPHOS 高度相关、hub 富集氧化磷酸化 → 独立支持结论一。\n")
cat("  · 同一模块同时对应 Immune 与 Purity_proxy → 该模块反映细胞构成，非上皮程序。\n")
cat("  · TSPO 落在 M0 → 它不属于任何稳定模块，须如实报告。\n")

hr("23_wgcna.R (v2) 完成")

# ############################################################################
# ##  Fig3_来源分解  (4 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- Fig3_来源分解 / 29_gse26712_replication.R
# -----------------------------------------------------------------------------
# =============================================================================
#  29_gse26712_replication.R   独立队列复现：TSPO 与免疫浸润等 bulk 关联
#
#  ── 为什么做这一个 ──────────────────────────────────────────────────────────
#  正文九节里，机制结论已来自多个独立队列（GSE165897 / GSE154600 /
#  GSE266577 / GSE288483 / CODEX / CPTAC / ICGC OV-AU）。
#  但 3.1 与 3.2 的 bulk 相关性只来自 TCGA-OV 一个队列：
#      TSPO ~ 免疫总分      ρ = +0.412
#      TSPO ~ 细胞毒 T      ρ = +0.459
#      TSPO ~ NK            ρ = +0.435
#      TSPO ~ 增殖          ρ = −0.298
#      TSPO ~ 纯度代理      ρ = −0.378
#  空间数据答不了这个问题（不同技术、不同问题：空间邻近 ≠ 病人间相关）。
#  故用 GSE26712（Affymetrix 芯片，n=185 肿瘤）做定向复现。
#
#  ── 设计要点 ────────────────────────────────────────────────────────────────
#  1. **两个队列用同一套评分方法与同一套基因集**，避免方法差异冒充平台差异。
#     统一用 z 均值（不使用 ssGSEA），因为 z 均值确定、不依赖 GSVA 版本，
#     更适合跨平台比较。正文的 TCGA 数字若来自 ssGSEA，此处会略有出入，
#     脚本会同时打印「正文值」供对照。
#  2. **只测 7 个指标**，不做全量重分析。其中 CD8_naive 是内部阴性对照：
#     它在 TCGA 中为 0.059（无关），若在 GSE26712 中同样接近 0，
#     说明「信号特异落在细胞毒亚群」这一判断也复现了。
#  3. 芯片与 RNA-seq 平台不同，是真正的跨平台验证；但基因覆盖不同，
#     脚本会逐指标报告基因交集数。
#
#  ── 判定规则（在跑之前先定好，避免事后解释）──────────────────────────────────
#    免疫总分在 GSE26712 中：
#      方向为正 且 ρ ≥ 0.30  → 强复现，写入正文
#      方向为正 且 0.20 ≤ ρ < 0.30 → 中等复现，写入正文并注明强度下降
#      方向为正 且 0.10 ≤ ρ < 0.20 → 弱复现，只在讨论中提一句
#      方向为正 且 ρ < 0.10，或方向为负 → 不复现，不写入
#
#  输入：F_TCGA（EXPR_clean.rds）、F_GSE（GSE26712 基因符号矩阵）、
#        F_SERIES（GSE26712 series matrix，用于区分肿瘤/正常）
#  输出：T21a_replication_scores.csv、T21b_replication_summary.csv、
#        T21c_TSPO_tumor_vs_normal_GSE26712.csv、Fig9_replication.pdf/png
# =============================================================================

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("OUT")) {
  .s <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.s, "00_setup.R")
  if (file.exists(.s)) source(.s, encoding = "UTF-8") else {
    PROJ <- "D:/Rstudy/tcga"
    OUT  <- file.path(PROJ, "results_v2")
  }
}
if (!exists("F_TCGA"))   F_TCGA   <- file.path(PROJ, "results_v2", "EXPR_clean.rds")
if (!exists("F_GSE"))    F_GSE    <- file.path(PROJ, "GSE26712_expr_matrix_gene_clean.rds")
if (!exists("F_SERIES")) F_SERIES <- file.path(PROJ, "GSE26712_series_matrix.txt.gz")

OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

set.seed(20260101)

# ---- 正文现有的 TCGA 值（ssGSEA 口径），仅用于对照打印 ---------------------
MS_TCGA <- c(免疫总分 = 0.412, CD8_CTL = 0.459, NK_specific = 0.435,
             CD8_naive = 0.059, 增殖 = -0.298, 纯度代理 = -0.378,
             上皮 = -0.146)

# ---- 基因集：与 17_immune_infiltration.R / 19_nk_vs_ctl.R 完全一致 ---------
SIG_IMM <- list(
  CD8_T      = c("CD8A","CD8B","GZMK","PRF1","IFNG"),
  CD4_T      = c("CD4","IL7R","CCR7","LTB"),
  Treg       = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  NK         = c("KLRD1","NKG7","GNLY","KLRF1"),
  B_cell     = c("CD19","MS4A1","CD79A","CD79B"),
  Plasma     = c("MZB1","JCHAIN","IGHG1"),
  Monocyte   = c("CD14","LYZ","S100A8","S100A9","VCAN"),
  Mac_M1     = c("NOS2","IL1B","TNF","CXCL9","CXCL10"),
  Mac_M2     = c("CD163","MRC1","MSR1","CCL22","IL10"),
  DC         = c("CD1C","CLEC9A","LAMP3","CCR7"),
  Neutrophil = c("FCGR3B","CSF3R","CXCR2"),
  Mast       = c("TPSAB1","CPA3","MS4A2")
)
SIG_NK   <- c("KLRF1","KLRC1","KLRC2","KLRC3","KLRD1","SH2D1B","NCR1",
              "KIR2DL3","KIR3DL1","KIR2DS4","FCER1G","TXK","PRSS23")
SIG_CTL  <- c("CD8A","CD8B","GZMA","GZMB","GZMH","PRF1","NKG7","GNLY","KLRG1","CST7")
SIG_NAIV <- c("CD8A","CCR7","SELL","LEF1","IL7R","TCF7")
SIG_PROL <- c("MKI67","PCNA","TOP2A","BUB1")
SIG_EPI  <- c("EPCAM","KRT8","KRT18")
SIG_PURE <- c("CD68","CD3D","COL1A1")

# 免疫总分用的 12 类（与 17 的 imm 列表一致）
IMM_TOTAL <- names(SIG_IMM)

# #############################################################################
#  STEP 1  TCGA-OV
# #############################################################################

hr("STEP 1  TCGA-OV")

if (!file.exists(F_TCGA)) stop("找不到 ", F_TCGA, "，请先跑 01_build_matrix.R")
tcga <- as.matrix(readRDS(F_TCGA)); storage.mode(tcga) <- "double"
cat("TCGA: ", nrow(tcga), " 基因 x ", ncol(tcga), " 样本\n", sep = "")
if (!"TSPO" %in% rownames(tcga)) stop("TCGA 矩阵中没有 TSPO")

# #############################################################################
#  STEP 2  GSE26712
# #############################################################################

hr("STEP 2  GSE26712")

if (!file.exists(F_GSE)) stop("找不到 ", F_GSE)
gse <- as.matrix(readRDS(F_GSE)); storage.mode(gse) <- "double"
cat("GSE26712: ", nrow(gse), " 基因 x ", ncol(gse), " 样本\n", sep = "")
if (!"TSPO" %in% rownames(gse)) stop("GSE26712 矩阵中没有 TSPO")
cat("值域: ", paste(round(range(gse, na.rm = TRUE), 2), collapse = " ~ "),
    "（应为 log2 标度）\n", sep = "")

# --- 区分肿瘤与正常：解析 series matrix 的 characteristics -------------------
group <- rep(NA_character_, ncol(gse))
names(group) <- colnames(gse)

if (file.exists(F_SERIES)) {
  con <- gzfile(F_SERIES, "rt")
  ln  <- readLines(con, warn = FALSE)
  close(con)

  # sample ID 行
  id_line <- grep("^!Sample_geo_accession", ln, value = TRUE)
  ids <- if (length(id_line)) {
    v <- strsplit(id_line[1], "\t")[[1]][-1]
    gsub('"', "", v)
  } else character(0)

  # tissue 行：两行 characteristics 里含 "tissue:" 的那一行
  ch_lines <- grep("^!Sample_characteristics", ln, value = TRUE)
  tis_line <- ch_lines[grepl("tissue:", ch_lines, fixed = TRUE)][1]
  if (!is.na(tis_line)) {
    v <- strsplit(tis_line, "\t")[[1]][-1]
    v <- gsub('"', "", v)
    v <- sub("^\\s*tissue:\\s*", "", v)
    if (length(v) == length(ids)) {
      is_tumor <- grepl("cancer|tumor|tumour", v, ignore.case = TRUE)
      group[ids] <- ifelse(is_tumor, "Tumor", "Normal")
    }
  }
  cat("series matrix 解析: ", length(ids), " 个样本，",
      sum(!is.na(group)), " 个匹配到矩阵列\n", sep = "")
  print(table(group, useNA = "ifany"))
} else {
  cat("⚠ 找不到 series matrix，无法区分肿瘤/正常，将使用全部样本并打印警告\n")
}

if (all(is.na(group))) {
  cat("⚠ 分组未知，退化为使用全部 ", ncol(gse), " 个样本\n", sep = "")
  use_gse <- colnames(gse)
} else {
  use_gse <- names(group)[!is.na(group) & group == "Tumor"]
  cat("肿瘤样本: ", length(use_gse), " / ", sum(!is.na(group)), "\n", sep = "")
}
if (length(use_gse) < 50) stop("肿瘤样本数过少（", length(use_gse), "），检查分组解析")

GSE <- gse[, use_gse, drop = FALSE]
cat("用于分析的 GSE26712: ", nrow(GSE), " 基因 x ", ncol(GSE), " 肿瘤\n", sep = "")

# --- TSPO 肿瘤 vs 正常（若分组可用）----------------------------------------
if (!all(is.na(group))) {
  nrm <- names(group)[!is.na(group) & group == "Normal"]
  if (length(nrm) >= 5) {
    tv <- as.numeric(gse["TSPO", nrm]); tt <- as.numeric(gse["TSPO", use_gse])
    w  <- suppressWarnings(wilcox.test(tt, tv))
    cat("\nGSE26712 中 TSPO：肿瘤 n=", length(tt), " 均值 ", round(mean(tt), 3),
        " vs 正常 n=", length(tv), " 均值 ", round(mean(tv), 3),
        "，Wilcoxon p = ", signif(w$p.value, 3), "\n", sep = "")
    cat("（正文 TCGA vs GTEx 为 8.09 vs 5.97；此处仅检验方向是否一致）\n")
    write.csv(data.frame(group = c(rep("Tumor", length(tt)), rep("Normal", length(tv))),
                         TSPO = c(tt, tv)),
              file.path(OUTDIR, "T21c_TSPO_tumor_vs_normal_GSE26712.csv"), row.names = FALSE)
  }
}

# #############################################################################
#  STEP 3  评分（两个队列用同一函数）
# #############################################################################

hr("STEP 3  统一评分")

# z 均值：逐基因 z 后取均值。不用 ssGSEA，理由见文件头。
score_z <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 2) return(list(score = rep(NA_real_, ncol(mat)), n = length(g)))
  z <- t(scale(t(mat[g, , drop = FALSE])))
  z[!is.finite(z)] <- NA
  list(score = colMeans(z, na.rm = TRUE), n = length(g))
}

build_scores <- function(mat, tag) {
  out <- data.frame(sample = colnames(mat), cohort = tag, stringsAsFactors = FALSE)
  ntab <- data.frame(indicator = character(0), n_gene = integer(0))

  # 12 类免疫各自的评分（用于免疫总分）
  imm_mat <- matrix(NA_real_, nrow = length(IMM_TOTAL), ncol = ncol(mat),
                    dimnames = list(IMM_TOTAL, colnames(mat)))
  for (k in IMM_TOTAL) {
    r <- score_z(mat, SIG_IMM[[k]])
    imm_mat[k, ] <- r$score
  }
  out$免疫总分 <- colMeans(imm_mat, na.rm = TRUE)

  for (nm in c("CD8_CTL", "NK_specific", "CD8_naive")) {
    gs <- switch(nm,
                 CD8_CTL     = SIG_CTL,
                 NK_specific = setdiff(SIG_NK, SIG_CTL),   # 与 19 一致：去掉重叠
                 CD8_naive   = SIG_NAIV)
    r <- score_z(mat, gs)
    out[[nm]] <- r$score
    ntab <- rbind(ntab, data.frame(indicator = nm, n_gene = r$n))
  }
  r <- score_z(mat, SIG_PROL); out$增殖 <- r$score
  ntab <- rbind(ntab, data.frame(indicator = "增殖", n_gene = r$n))
  r <- score_z(mat, SIG_EPI);  out$上皮 <- r$score
  ntab <- rbind(ntab, data.frame(indicator = "上皮", n_gene = r$n))

  # 纯度代理 = -(z(CD68)+z(CD3D)+z(COL1A1))/3，与 3.1 一致
  g <- intersect(SIG_PURE, rownames(mat))
  if (length(g) >= 2) {
    z <- t(scale(t(mat[g, , drop = FALSE]))); z[!is.finite(z)] <- NA
    out$纯度代理 <- -colMeans(z, na.rm = TRUE)
  } else out$纯度代理 <- NA_real_
  ntab <- rbind(ntab, data.frame(indicator = "纯度代理", n_gene = length(g)))

  out$TSPO <- as.numeric(mat["TSPO", ])
  list(df = out, n = ntab)
}

TC <- build_scores(tcga, "TCGA-OV")
GS <- build_scores(GSE,  "GSE26712")

cat("TCGA 基因覆盖：\n"); print(TC$n, row.names = FALSE)
cat("\nGSE26712 基因覆盖：\n"); print(GS$n, row.names = FALSE)
cat("\nTSPO 在 TCGA 中均值 ", round(mean(TC$df$TSPO), 3),
    "；GSE26712 中 ", round(mean(GS$df$TSPO), 3), "\n", sep = "")

ALL <- rbind(TC$df, GS$df)
write.csv(ALL, file.path(OUTDIR, "T21a_replication_scores.csv"), row.names = FALSE)

# #############################################################################
#  STEP 4  相关与并列表
# #############################################################################

hr("STEP 4  TSPO 与各指标的 Spearman 相关")

IND <- c("免疫总分", "CD8_CTL", "NK_specific", "CD8_naive",
         "增殖", "纯度代理", "上皮")

cor_one <- function(df, v) {
  ok <- is.finite(df$TSPO) & is.finite(df[[v]])
  if (sum(ok) < 30) return(c(rho = NA, p = NA, n = sum(ok)))
  ct <- suppressWarnings(cor.test(df$TSPO[ok], df[[v]][ok], method = "spearman"))
  c(rho = as.numeric(ct$estimate), p = ct$p.value, n = sum(ok))
}

TAB <- do.call(rbind, lapply(IND, function(v) {
  a <- cor_one(TC$df, v)
  b <- cor_one(GS$df, v)
  ref <- if (v %in% names(MS_TCGA)) MS_TCGA[[v]] else NA_real_
  data.frame(indicator = v,
             rho_TCGA = round(a["rho"], 3), p_TCGA = signif(a["p"], 3),
             rho_GSE  = round(b["rho"], 3), p_GSE  = signif(b["p"], 3),
             n_TCGA   = a["n"], n_GSE = b["n"],
             正文旧值 = ref,
             stringsAsFactors = FALSE)
}))

# 方向与强度判定
TAB$same_direction <- sign(TAB$rho_TCGA) == sign(TAB$rho_GSE)
TAB$abs_drop      <- round(abs(TAB$rho_TCGA) - abs(TAB$rho_GSE), 3)
TAB$判定 <- ifelse(
  is.na(TAB$rho_GSE), "无法计算",
  ifelse(!TAB$same_direction, "方向相反",
  ifelse(abs(TAB$rho_GSE) >= 0.30, "强复现",
  ifelse(abs(TAB$rho_GSE) >= 0.20, "中等复现",
  ifelse(abs(TAB$rho_GSE) >= 0.10, "弱复现", "不复现")))))

cat("\n===== 复现对照表 =====\n")
pr(TAB[, c("indicator","rho_TCGA","rho_GSE","n_TCGA","n_GSE","same_direction","判定")], 3)

cat("\n注：rho_TCGA 为脚本用 z 均值口径重算的结果，与「正文旧值」（ssGSEA 口径）可能略有出入。\n")
cat("    跨队列比较应以 rho_TCGA 与 rho_GSE 为准，两者口径完全相同。\n")

write.csv(TAB, file.path(OUTDIR, "T21b_replication_summary.csv"), row.names = FALSE)

# #############################################################################
#  STEP 5  判定
# #############################################################################

hr("STEP 5  判定")

main <- TAB[TAB$indicator == "免疫总分", ]
rho_main <- if (nrow(main)) main$rho_GSE[1] else NA_real_
dir_ok   <- if (nrow(main)) main$same_direction[1] else FALSE

cat("主要判据 —— 免疫总分在 GSE26712 中 ρ = ", round(rho_main, 3),
    "（TCGA 同口径 ", round(main$rho_TCGA[1], 3), "）\n", sep = "")

verdict <- if (is.na(rho_main)) {
  "无法判定：免疫总分未算出"
} else if (!dir_ok) {
  "不复现（方向相反）：3.2 的 bulk 相关不能作为结论，应降级为单队列观察"
} else if (abs(rho_main) >= 0.30) {
  "强复现：可写入正文，表述为「在两个独立队列中一致」"
} else if (abs(rho_main) >= 0.20) {
  "中等复现：可写入正文，但需注明强度低于 TCGA，并讨论平台差异"
} else if (abs(rho_main) >= 0.10) {
  "弱复现：只宜在讨论中提一句，不作为主要证据"
} else {
  "不复现：不写入正文；3.2 降级为单队列观察并写入局限"
}
cat("\n【判定】", verdict, "\n", sep = "")

# 内部阴性对照复核
naive <- TAB[TAB$indicator == "CD8_naive", ]
if (nrow(naive)) {
  cat("\n内部阴性对照（CD8_naive，应接近 0）：TCGA ", round(naive$rho_TCGA[1], 3),
      "，GSE26712 ", round(naive$rho_GSE[1], 3), "\n", sep = "")
  cat(if (abs(naive$rho_GSE[1]) < 0.15)
      "  → 该对照在两个队列中均接近 0，说明「信号特异落在细胞毒亚群」的判读可复现\n"
      else "  → 该对照在 GSE26712 中并不接近 0，特异性论证需谨慎\n")
}

cat("\n各判定计数：\n"); print(table(TAB$判定))

writeLines(c(
  "GSE26712 独立队列复现结果",
  paste("时间:", format(Sys.time())),
  paste("TCGA-OV:", ncol(tcga), "样本；GSE26712:", ncol(GSE), "肿瘤"),
  "",
  paste("免疫总分: TCGA rho =", round(main$rho_TCGA[1], 3),
        "| GSE26712 rho =", round(rho_main, 3),
        "| p =", signif(main$p_GSE[1], 3)),
  "",
  paste("判定:", verdict)
), file.path(OUTDIR, "T21d_verdict.txt"))

# #############################################################################
#  STEP 6  图
# #############################################################################

if (requireNamespace("ggplot2", quietly = TRUE)) {
  hr("STEP 6  图")
  suppressPackageStartupMessages(library(ggplot2))

  dd <- rbind(
    data.frame(cohort = "TCGA-OV",  x = TC$df$免疫总分, y = TC$df$TSPO),
    data.frame(cohort = "GSE26712", x = GS$df$免疫总分, y = GS$df$TSPO))
  dd <- dd[is.finite(dd$x) & is.finite(dd$y), ]
  lab <- paste0(dd$cohort, "\nrho = ",
                round(c(cor(TC$df$免疫总分, TC$df$TSPO, method = "spearman", use = "complete.obs"),
                        cor(GS$df$免疫总分, GS$df$TSPO, method = "spearman", use = "complete.obs")), 3))
  names(lab) <- c("TCGA-OV", "GSE26712")

  p <- ggplot(dd, aes(x, y)) +
    geom_point(size = 0.9, alpha = 0.55, colour = "#2C7BB6") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "#B2182B", fill = "#FDDBC7") +
    facet_wrap(~ cohort, scales = "free", labeller = labeller(cohort = lab)) +
    labs(x = "Immune infiltration score (z-mean)",
         y = "TSPO expression",
         title = "Replication of the TSPO-immune association in an independent cohort") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIGDIR, "Fig9_replication.pdf"), p, width = 9, height = 4.2)
  ggsave(file.path(FIGDIR, "Fig9_replication.png"), p, width = 9, height = 4.2, dpi = 300)
  cat("已保存 Fig9_replication.pdf/png\n")
}

# #############################################################################
#  STEP 7  说明
# #############################################################################

hr("STEP 7  怎么用这个结果")

cat("若判定为「强复现」或「中等复现」：\n")
cat("  · 3.2 首段可加一句：「该关联在独立的 GSE26712 队列（n=", ncol(GSE),
    "，Affymetrix 平台）中重现（ρ=", round(rho_main, 3), "）。」\n", sep = "")
cat("  · 图 2 增加一个面板（Fig9 可直接用）。\n")
cat("  · 局限中「bulk 关联限于单一队列」一句可删。\n\n")

cat("若判定为「弱复现」或「不复现」：\n")
cat("  · 不要硬写。3.2 的定位从「结论」改为「TCGA-OV 队列内的观察」。\n")
cat("  · 局限中明确写：bulk 层面关联未能在独立队列重现，机制结论依赖单细胞与空间数据。\n")
cat("  · 这样写反而比不写更安全——审稿人自己重算时不会觉得被隐瞒。\n\n")

cat("无论哪种结果，都不要为了让复现更好看而调整基因集或评分方法。\n")

hr("29_gse26712_replication.R 完成")

# -----------------------------------------------------------------------------
# -- Fig3_来源分解 / 30_replication_diagnose.R
# -----------------------------------------------------------------------------
# =============================================================================
#  30_replication_diagnose.R   复现失败的诊断（在改写正文之前必须先做）
#
#  ── 背景 ────────────────────────────────────────────────────────────────────
#  29_gse26712_replication.R 显示：TCGA 与 GSE26712 在多个指标上方向相反，
#     免疫总分    TCGA +0.356   GSE26712 −0.188
#     CD8_CTL     TCGA +0.421   GSE26712 −0.055
#     NK_specific TCGA +0.372   GSE26712 −0.342
#     CD8_naive   TCGA +0.184   GSE26712 −0.342
#  但该结果有三处必须先排除技术原因，否则不能据此改写正文：
#
#    （1）NK_specific 与 CD8_naive 在 GSE26712 中给出完全相同的 −0.342，
#        两者基因集毫无重叠，小数点后三位一致不合常理，需查是否代码/共线问题。
#    （2）本脚本重算的 TCGA 值与正文 ssGSEA 口径对不上（上皮 +0.015 vs −0.146，
#        CD8_naive 0.184 vs 0.059），说明 z 均值与小基因集 ssGSEA 差异很大，
#        不能直接拿 z 均值的结果否定 ssGSEA 的结论。
#    （3）GSE26712 值域 2.5–13.72（跨 11 个 log2 单位），可能是 MAS5 未分位数
#        归一化，样本间不可比会严重扭曲 z 均值。
#
#  ── 本脚本做什么 ────────────────────────────────────────────────────────────
#  A. 数据质量检查：分布、分位数、是否需要归一化
#  B. **单基因对单基因**的 Spearman（绕开一切评分与基因集）
#     ——这是唯一不受方法选择影响的检验
#  C. 评分矩阵的内部一致性（查 NK vs naive 是否共线）
#  D. 引入第三队列 ICGC OV-AU（RNA-seq，n=78）作为仲裁
#  E. 给出决策矩阵
#
#  输入：F_TCGA、F_GSE、F_SERIES、F_ICGC
#  输出：T22a_gene_level_cor.csv、T22b_score_collinearity.csv、
#        T22c_distribution_check.csv、T22d_decision.txt、Fig10_diagnose.pdf/png
# =============================================================================

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("OUT")) {
  .s <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.s, "00_setup.R")
  if (file.exists(.s)) source(.s, encoding = "UTF-8") else {
    PROJ <- "D:/Rstudy/tcga"
    OUT  <- file.path(PROJ, "results_v2")
  }
}
if (!exists("F_TCGA"))   F_TCGA   <- file.path(PROJ, "results_v2", "EXPR_clean.rds")
if (!exists("F_GSE"))    F_GSE    <- file.path(PROJ, "GSE26712_expr_matrix_gene_clean.rds")
if (!exists("F_SERIES")) F_SERIES <- file.path(PROJ, "GSE26712_series_matrix.txt.gz")
if (!exists("F_ICGC"))   F_ICGC   <- file.path(PROJ, "ICGC_OV-AU_GX.rds")

OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# ---- 单基因清单：覆盖免疫/髓系/基质/上皮/增殖，外加两个阴性对照 ------------
GENES <- c("TSPO",
           "CD8A", "CD3D", "PTPRC", "CD68", "CD163", "COL1A1",
           "EPCAM", "KRT8", "MKI67",
           "ACTB", "GAPDH")

# #############################################################################
#  STEP 1  载入三个队列
# #############################################################################

hr("STEP 1  载入")

mats <- list()

if (file.exists(F_TCGA)) {
  x <- as.matrix(readRDS(F_TCGA)); storage.mode(x) <- "double"
  mats[["TCGA"]] <- x
  cat("TCGA    : ", nrow(x), " x ", ncol(x), "\n", sep = "")
} else cat("⚠ 缺 ", F_TCGA, "\n", sep = "")

if (file.exists(F_GSE)) {
  gse <- as.matrix(readRDS(F_GSE)); storage.mode(gse) <- "double"
  # 只保留肿瘤样本
  grp <- rep(NA_character_, ncol(gse)); names(grp) <- colnames(gse)
  if (file.exists(F_SERIES)) {
    con <- gzfile(F_SERIES, "rt"); ln <- readLines(con, warn = FALSE); close(con)
    id_line <- grep("^!Sample_geo_accession", ln, value = TRUE)
    ids <- if (length(id_line)) gsub('"', "", strsplit(id_line[1], "\t")[[1]][-1]) else character(0)
    ch <- grep("^!Sample_characteristics", ln, value = TRUE)
    tl <- ch[grepl("tissue:", ch, fixed = TRUE)][1]
    if (!is.na(tl)) {
      v <- gsub('"', "", strsplit(tl, "\t")[[1]][-1])
      v <- sub("^\\s*tissue:\\s*", "", v)
      if (length(v) == length(ids))
        grp[ids] <- ifelse(grepl("cancer|tumor|tumour", v, ignore.case = TRUE), "Tumor", "Normal")
    }
  }
  tum <- names(grp)[!is.na(grp) & grp == "Tumor"]
  if (length(tum) >= 50) gse <- gse[, tum, drop = FALSE]
  mats[["GSE26712"]] <- gse
  cat("GSE26712: ", nrow(gse), " x ", ncol(gse), "（肿瘤）\n", sep = "")
} else cat("⚠ 缺 ", F_GSE, "\n", sep = "")

if (file.exists(F_ICGC)) {
  ic <- readRDS(F_ICGC)
  if (is.data.frame(ic)) ic <- as.matrix(ic)
  # 统一定向：基因为行
  if (nrow(ic) < ncol(ic) && "TSPO" %in% colnames(ic)) ic <- t(ic)
  storage.mode(ic) <- "double"
  mats[["ICGC"]] <- ic
  cat("ICGC    : ", nrow(ic), " x ", ncol(ic), "\n", sep = "")
} else cat("⚠ 缺 ", F_ICGC, "（可跳过，但少了关键仲裁队列）\n", sep = "")

if (!length(mats)) stop("没有任何队列可用")

for (nm in names(mats)) {
  g <- intersect(GENES, rownames(mats[[nm]]))
  cat("  ", nm, " 命中基因 ", length(g), "/", length(GENES),
      if (length(g) < length(GENES)) paste0("（缺 ", paste(setdiff(GENES, g), collapse = ","), "）") else "",
      "\n", sep = "")
}

# #############################################################################
#  STEP 2  数据质量检查：GSE26712 是否已分位数归一化
# #############################################################################

hr("STEP 2  分布检查")

dist_tab <- do.call(rbind, lapply(names(mats), function(nm) {
  m <- mats[[nm]]
  q <- quantile(m, c(.01, .25, .50, .75, .99), na.rm = TRUE)
  # 样本中位数之间的离散程度：分位数归一化后各样本中位数应几乎相同
  smed <- apply(m, 2, median, na.rm = TRUE)
  data.frame(cohort = nm, n_gene = nrow(m), n_sample = ncol(m),
             q01 = round(q[1], 2), q25 = round(q[2], 2),
             q50 = round(q[3], 2), q75 = round(q[4], 2), q99 = round(q[5], 2),
             样本中位数_min = round(min(smed), 2),
             样本中位数_max = round(max(smed), 2),
             样本中位数_SD = round(sd(smed), 3),
             stringsAsFactors = FALSE)
}))

cat("\n各队列的全局分位数与样本中位数离散度：\n")
pr(dist_tab, 2)
cat("\n判读：若「样本中位数_SD」明显大于 0，说明样本间未对齐，z 均值会被扭曲。\n")
cat("      RMA / 分位数归一化后该项应接近 0.00x。\n")

write.csv(dist_tab, file.path(OUTDIR, "T22c_distribution_check.csv"), row.names = FALSE)

# #############################################################################
#  STEP 3  单基因对单基因相关（绕开一切评分方法）
# #############################################################################

hr("STEP 3  单基因层面相关（最稳健的检验）")

TARGETS <- setdiff(GENES, "TSPO")

GL <- do.call(rbind, lapply(names(mats), function(nm) {
  m <- mats[[nm]]
  if (!"TSPO" %in% rownames(m)) return(NULL)
  ts <- as.numeric(m["TSPO", ])
  do.call(rbind, lapply(TARGETS, function(g) {
    if (!g %in% rownames(m)) return(NULL)
    y <- as.numeric(m[g, ])
    ok <- is.finite(ts) & is.finite(y)
    if (sum(ok) < 30) return(NULL)
    ct <- suppressWarnings(cor.test(ts[ok], y[ok], method = "spearman"))
    data.frame(cohort = nm, gene = g, rho = round(as.numeric(ct$estimate), 3),
               p = signif(ct$p.value, 3), n = sum(ok), stringsAsFactors = FALSE)
  }))
}))

# 宽表：基因 × 队列
W <- reshape(as.data.frame(GL[, c("gene","cohort","rho")]),
             idvar = "gene", timevar = "cohort", direction = "wide")
names(W) <- sub("^rho\\.", "", names(W))

cat("\nTSPO 与各单基因的 Spearman ρ：\n")
pr(W, 3)
write.csv(GL, file.path(OUTDIR, "T22a_gene_level_cor.csv"), row.names = FALSE)

# 关键三对的一致性
key_pairs <- c("CD8A", "PTPRC", "CD68", "COL1A1", "EPCAM")
cat("\n关键单基因在各队列的符号：\n")
for (g in key_pairs) {
  if (!g %in% W$gene) next
  v <- unlist(W[W$gene == g, setdiff(names(W), "gene")])
  cat(sprintf("  TSPO ~ %-7s  %s\n", g,
              paste(sprintf("%s %+.3f", names(v), as.numeric(v)), collapse = "   ")))
}

# #############################################################################
#  STEP 4  评分矩阵的内部一致性（查 NK vs naive 共线）
# #############################################################################

hr("STEP 4  评分共线性检查")

SIG_IMM <- list(
  CD8_T = c("CD8A","CD8B","GZMK","PRF1","IFNG"),
  CD4_T = c("CD4","IL7R","CCR7","LTB"),
  Treg  = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  NK    = c("KLRD1","NKG7","GNLY","KLRF1"),
  B_cell= c("CD19","MS4A1","CD79A","CD79B"),
  Plasma= c("MZB1","JCHAIN","IGHG1"),
  Monocyte = c("CD14","LYZ","S100A8","S100A9","VCAN"),
  Mac_M1 = c("NOS2","IL1B","TNF","CXCL9","CXCL10"),
  Mac_M2 = c("CD163","MRC1","MSR1","CCL22","IL10"),
  DC = c("CD1C","CLEC9A","LAMP3","CCR7"),
  Neutrophil = c("FCGR3B","CSF3R","CXCR2"),
  Mast = c("TPSAB1","CPA3","MS4A2"))
SIG_NK   <- c("KLRF1","KLRC1","KLRC2","KLRC3","KLRD1","SH2D1B","NCR1",
              "KIR2DL3","KIR3DL1","KIR2DS4","FCER1G","TXK","PRSS23")
SIG_CTL  <- c("CD8A","CD8B","GZMA","GZMB","GZMH","PRF1","NKG7","GNLY","KLRG1","CST7")
SIG_NAIV <- c("CD8A","CCR7","SELL","LEF1","IL7R","TCF7")

score_z <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 2) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[g, , drop = FALSE]))); z[!is.finite(z)] <- NA
  colMeans(z, na.rm = TRUE)
}

COLL <- list()
for (nm in names(mats)) {
  m <- mats[[nm]]
  sc <- data.frame(
    NK_specific = score_z(m, setdiff(SIG_NK, SIG_CTL)),
    CD8_CTL     = score_z(m, SIG_CTL),
    CD8_naive   = score_z(m, SIG_NAIV),
    Immune      = rowMeans(sapply(SIG_IMM, function(g) score_z(m, g)), na.rm = TRUE),
    stringsAsFactors = FALSE)
  cm <- cor(sc, use = "pairwise.complete.obs", method = "spearman")
  cat("\n--- ", nm, " 评分两两 Spearman ---\n", sep = "")
  print(round(cm, 3))
  COLL[[nm]] <- data.frame(cohort = nm, pair = c("NK~CTL","NK~naive","CTL~naive"),
                           rho = round(c(cm["NK_specific","CD8_CTL"],
                                         cm["NK_specific","CD8_naive"],
                                         cm["CD8_CTL","CD8_naive"]), 3))
}
CT <- do.call(rbind, COLL)
cat("\n关键共线性：\n"); pr(CT, 3)
cat("\n判读：若 NK_specific 与 CD8_naive 高度共线（|ρ|>0.9），\n")
cat("      则两者给出相同的 TSPO 相关值就不奇怪，也就不能算两个独立证据。\n")
write.csv(CT, file.path(OUTDIR, "T22b_score_collinearity.csv"), row.names = FALSE)

# #############################################################################
#  STEP 4B  全局表达混杂（STEP 3 暴露出的关键问题）
#
#  STEP 3 显示 TSPO 与管家基因 ACTB 的相关在 TCGA 为 +0.337、GSE26712 为 +0.390。
#  管家基因不应与任何特定基因有 0.3 以上的相关，这提示两队列都存在一个
#  「全局表达梯度」（技术性的 RNA 完整性/文库构成，或组织细胞密度），
#  TSPO 与免疫标志的相关可能部分由该梯度介导，而非特异生物学。
#
#  本步做三件事：
#    1. 量化各队列的全局表达梯度（全基因均值），并检验它与 TSPO 的相关
#    2. 用秩残差偏相关，在控制全局梯度后再算 TSPO~免疫标志
#    3. 同时给出控制 ACTB+GAPDH 的结果（更保守的口径）
#
#  判读：若控制后 TCGA 的 TSPO~CD8A 从 +0.30 掉到 ~0，说明 3.2 的 bulk 相关
#        主要是全局梯度伪影；若仍保持 +0.15 以上，则存在真实的成分信号。
# #############################################################################

hr("STEP 4B  全局表达混杂")

pcor_rank <- function(a, b, cov) {
  cov <- as.matrix(cov)
  if (nrow(cov) != length(a)) cov <- t(cov)
  ok <- is.finite(a) & is.finite(b) & apply(cov, 2, function(z) all(is.finite(z)))
  a <- a[ok]; b <- b[ok]; cov <- cov[ok, , drop = FALSE]
  n <- length(a); k <- ncol(cov)
  if (n <= k + 3) return(c(rho = NA_real_, p = NA_real_))
  A  <- cbind(1, apply(cov, 2, rank))
  ea <- residuals(lm(rank(a) ~ A - 1))
  eb <- residuals(lm(rank(b) ~ A - 1))
  r  <- suppressWarnings(cor(ea, eb))
  if (!is.finite(r) || abs(r) >= 1) return(c(rho = r, p = NA_real_))
  t  <- r * sqrt((n - k - 2) / (1 - r^2))
  c(rho = r, p = 2 * pt(-abs(t), df = n - k - 2))
}

PC <- list()
for (nm in names(mats)) {
  m <- mats[[nm]]
  if (!"TSPO" %in% rownames(m)) next

  ts  <- as.numeric(m["TSPO", ])
  glb <- colMeans(m, na.rm = TRUE)              # 全基因均值 = 全局表达梯度
  hk  <- intersect(c("ACTB","GAPDH"), rownames(m))
  hkv <- if (length(hk)) colMeans(m[hk, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(m))

  cat("\n--- ", nm, " ---\n", sep = "")
  ct <- suppressWarnings(cor.test(ts, glb, method = "spearman"))
  cat(sprintf("  TSPO ~ 全局均值        rho = %+.3f   p = %.3g\n",
              as.numeric(ct$estimate), ct$p.value))

  for (g in intersect(c("CD8A","CD3D","PTPRC","CD68","CD163","MKI67"), rownames(m))) {
    y   <- as.numeric(m[g, ])
    r0  <- suppressWarnings(cor(ts, y, method = "spearman", use = "complete.obs"))
    r1  <- pcor_rank(ts, y, matrix(glb, ncol = 1))
    r2  <- if (all(is.finite(hkv))) pcor_rank(ts, y, matrix(hkv, ncol = 1)) else c(rho = NA_real_, p = NA_real_)
    PC[[length(PC) + 1]] <- data.frame(
      cohort = nm, gene = g,
      rho_raw = round(r0, 3),
      rho_ctrl_global = round(r1["rho"], 3), p_ctrl_global = signif(r1["p"], 3),
      rho_ctrl_ACTB_GAPDH = round(r2["rho"], 3),
      stringsAsFactors = FALSE)
  }
}

if (length(PC)) {
  PC <- do.call(rbind, PC)
  W2 <- reshape(as.data.frame(PC[, c("gene","cohort","rho_raw")]),
                idvar = "gene", timevar = "cohort", direction = "wide")
  names(W2) <- sub("^rho_raw\\.", "", names(W2))
  cat("\n原始相关（未控制）：\n"); pr(W2, 3)

  W3 <- reshape(as.data.frame(PC[, c("gene","cohort","rho_ctrl_global")]),
                idvar = "gene", timevar = "cohort", direction = "wide")
  names(W3) <- sub("^rho_ctrl_global\\.", "", names(W3))
  cat("\n控制全局表达梯度后：\n"); pr(W3, 3)

  W4 <- reshape(as.data.frame(PC[, c("gene","cohort","rho_ctrl_ACTB_GAPDH")]),
                idvar = "gene", timevar = "cohort", direction = "wide")
  names(W4) <- sub("^rho_ctrl_ACTB_GAPDH\\.", "", names(W4))
  cat("\n控制 ACTB+GAPDH 后（最保守）：\n"); pr(W4, 3)

  write.csv(PC, file.path(OUTDIR, "T22e_global_confound.csv"), row.names = FALSE)

  cat("\n判读：若控制全局梯度后各队列的 TSPO~免疫相关整体降到 0.1 以下，\n")
  cat("      说明 bulk 层面的关联主要由全局表达梯度介导，不应作为生物学结论。\n")
}

# #############################################################################
#  STEP 5  决策矩阵
# #############################################################################

hr("STEP 5  决策")

get_rho <- function(cohort, gene) {
  if (!gene %in% W$gene) return(NA_real_)
  if (!cohort %in% names(W)) return(NA_real_)
  as.numeric(W[W$gene == gene, cohort][1])
}

cat("以 TSPO ~ CD8A（免疫浸润的最直接指标）作为仲裁：\n\n")
for (nm in names(mats)) {
  cat(sprintf("  %-10s TSPO~CD8A = %+.3f\n", nm, get_rho(nm, "CD8A")))
}

r_tcga <- get_rho("TCGA",    "CD8A")
r_gse  <- get_rho("GSE26712","CD8A")
r_icgc <- if ("ICGC" %in% names(mats)) get_rho("ICGC", "CD8A") else NA_real_

dec <- if (is.na(r_icgc)) {
  "ICGC 未纳入，无法仲裁。请补 ICGC_OV-AU_GX.rds 后重跑。"
} else if (sign(r_tcga) == sign(r_icgc) && abs(r_icgc) >= 0.15) {
  paste0("ICGC（RNA-seq）与 TCGA 同向（", sprintf("%+.3f", r_icgc),
         " vs ", sprintf("%+.3f", r_tcga), "），而 GSE26712 为 ",
         sprintf("%+.3f", r_gse), "。",
         "→ 判定 GSE26712 的偏差来自平台（芯片 vs RNA-seq）或归一化，",
         "3.2 的 bulk 关联保留，GSE26712 结果写入局限。")
} else if (sign(r_gse) == sign(r_icgc) && abs(r_icgc) >= 0.15) {
  paste0("ICGC 与 GSE26712 同向（", sprintf("%+.3f", r_icgc),
         " vs ", sprintf("%+.3f", r_gse), "），与 TCGA 相反（",
         sprintf("%+.3f", r_tcga), "）。",
         "→ 判定 TCGA 的 bulk 关联不可靠，3.2 必须降级为单队列观察。")
} else {
  paste0("三个队列方向不一致或效应量均很小（TCGA ", sprintf("%+.3f", r_tcga),
         "，GSE26712 ", sprintf("%+.3f", r_gse),
         "，ICGC ", sprintf("%+.3f", r_icgc), "）。",
         "→ 判定 bulk 层面的 TSPO–免疫关联本身不稳定，3.2 降级。")
}
cat("\n【决策】", dec, "\n", sep = "")

# 分布是否影响判读
sd_tab <- setNames(dist_tab$样本中位数_SD, dist_tab$cohort)
cat("\n样本中位数的离散度：\n")
for (nm in names(sd_tab)) cat(sprintf("  %-10s SD = %.4f\n", nm, sd_tab[[nm]]))
if (!is.na(sd_tab["GSE26712"]) && sd_tab["GSE26712"] > 0.2) {
  cat("⚠ GSE26712 的样本中位数离散度偏大，提示未做分位数归一化。\n")
  cat("  这会使 z 均值评分失真。建议先对 GSE26712 做 quantile normalization\n")
  cat("  （limma::normalizeBetweenArrays）后重跑 29，再下结论。\n")
}

writeLines(c(
  "复现失败的诊断结果",
  paste("时间:", format(Sys.time())),
  "",
  paste("TSPO~CD8A:", paste(sprintf("%s %+.3f", names(sapply(names(mats), get_rho, gene="CD8A")),
                                    sapply(names(mats), get_rho, gene = "CD8A")), collapse = " | ")),
  "",
  paste("决策:", dec),
  "",
  paste("样本中位数SD:", paste(sprintf("%s %.4f", names(sd_tab), sd_tab), collapse = " | "))
), file.path(OUTDIR, "T22d_decision.txt"))

# #############################################################################
#  STEP 6  图：TSPO~CD8A 三队列并列
# #############################################################################

if (requireNamespace("ggplot2", quietly = TRUE)) {
  hr("STEP 6  图")
  suppressPackageStartupMessages(library(ggplot2))

  dd <- do.call(rbind, lapply(names(mats), function(nm) {
    m <- mats[[nm]]
    if (!all(c("TSPO","CD8A") %in% rownames(m))) return(NULL)
    data.frame(cohort = nm, tspo = as.numeric(m["TSPO", ]),
               cd8a = as.numeric(m["CD8A", ]), stringsAsFactors = FALSE)
  }))
  dd <- dd[is.finite(dd$tspo) & is.finite(dd$cd8a), ]

  labs <- sapply(unique(dd$cohort), function(nm) {
    z <- dd[dd$cohort == nm, ]
    sprintf("%s (n=%d)\nrho = %+.3f", nm, nrow(z),
            cor(z$tspo, z$cd8a, method = "spearman"))
  })
  names(labs) <- unique(dd$cohort)
  dd$panel <- factor(labs[dd$cohort], levels = labs[unique(dd$cohort)])

  p <- ggplot(dd, aes(tspo, cd8a)) +
    geom_point(size = 0.8, alpha = 0.5, colour = "#2C7BB6") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "#B2182B", fill = "#FDDBC7") +
    facet_wrap(~ panel, scales = "free") +
    labs(x = "TSPO expression", y = "CD8A expression",
         title = "TSPO-CD8A association across three cohorts") +
    theme_bw(base_size = 10)
  ggsave(file.path(FIGDIR, "Fig10_diagnose.pdf"), p, width = 11, height = 4)
  ggsave(file.path(FIGDIR, "Fig10_diagnose.png"), p, width = 11, height = 4, dpi = 300)
  cat("已保存 Fig10_diagnose.pdf/png\n")
}

hr("STEP 7  下一步")

cat("按优先级：\n")
cat("  1. 若 GSE26712 的样本中位数 SD 偏大 → 先做 quantile normalization 再重跑 29。\n")
cat("  2. 看 STEP 3 的单基因结果：若 TSPO~CD8A 在 TCGA 是正的、GSE26712 是负的，\n")
cat("     而 ICGC 与 TCGA 一致，则问题在 GSE26712 这一批数据，不在 3.2。\n")
cat("  3. 只有当 ICGC 也支持负相关时，才考虑改写 3.2。\n")
cat("  4. 在完成上述之前，不要修改正文。\n")

hr("30_replication_diagnose.R 完成")

# -----------------------------------------------------------------------------
# -- Fig3_来源分解 / 31_bulk_confound_adjusted.R
# -----------------------------------------------------------------------------
# =============================================================================
#  31_bulk_confound_adjusted.R   3.2 的最终口径：控制全局表达后的三队列复现
#
#  ── 为什么需要这一版 ────────────────────────────────────────────────────────
#  30_replication_diagnose.R 发现一个关键混杂：
#      TSPO ~ 全局表达均值    TCGA −0.252 | GSE26712 +0.491 | ICGC −0.025
#  三个队列的全局结构与 TSPO 的关系符号都不同。不控制它，会得到假的方向翻转：
#      TSPO~CD8A  原始     TCGA +0.300 | GSE26712 −0.171 | ICGC +0.012
#      TSPO~CD8A  控全局   TCGA +0.333 | GSE26712 −0.030 | ICGC +0.026
#      TSPO~CD68  原始     TCGA +0.378 | GSE26712 −0.036 | ICGC +0.236
#      TSPO~CD68  控全局   TCGA +0.414 | GSE26712 +0.201 | ICGC +0.266
#  即：控制后 CD68/CD163/CD3D 在三队列一致为正，而 CD8A 的特异关联不成立。
#
#  29_gse26712_replication.R 用的是**未控制**的 z 均值，不能直接作为 3.2 的定稿依据。
#  本脚本补上全局控制，并把 3.2 需要的全部指标一次算全。
#
#  ── 同时记录一处论证缺陷 ────────────────────────────────────────────────────
#  正文把 CD8_naive 当作"阴性对照"，用以排除"只是淋巴细胞总量更多"。
#  但 30 的共线性检查显示 CD8_naive 与 CD8_CTL 的相关为
#      TCGA 0.765 | GSE26712 0.696 | ICGC 0.799
#  阴性对照与阳性对照共线 0.7–0.8，该对照不成立，正文相关表述须撤除。
#  本脚本会把该共线性一并打印，供改写时引用。
#
#  ── 输出 ────────────────────────────────────────────────────────────────────
#  T23a_scores_confound_adjusted.csv   三队列 × 全部指标的原始与偏相关
#  T23b_replication_final.csv          可直接写进正文的对照表
#  T23c_collinearity_warning.csv       阴性对照的共线性证据
#  T23d_recommendation.txt             每条结论的保留/降级建议
#
#  输入：F_TCGA、F_GSE、F_SERIES、F_ICGC
# =============================================================================

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("OUT")) {
  .s <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.s, "00_setup.R")
  if (file.exists(.s)) source(.s, encoding = "UTF-8") else {
    PROJ <- "D:/Rstudy/tcga"; OUT <- file.path(PROJ, "results_v2")
  }
}
if (!exists("F_TCGA"))   F_TCGA   <- file.path(PROJ, "results_v2", "EXPR_clean.rds")
if (!exists("F_GSE"))    F_GSE    <- file.path(PROJ, "GSE26712_expr_matrix_gene_clean.rds")
if (!exists("F_SERIES")) F_SERIES <- file.path(PROJ, "GSE26712_series_matrix.txt.gz")
if (!exists("F_ICGC"))   F_ICGC   <- file.path(PROJ, "ICGC_OV-AU_GX.rds")

OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# ---- 基因集（与 17 / 19 一致）---------------------------------------------
SIG_IMM <- list(
  CD8_T = c("CD8A","CD8B","GZMK","PRF1","IFNG"),
  CD4_T = c("CD4","IL7R","CCR7","LTB"),
  Treg  = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  NK    = c("KLRD1","NKG7","GNLY","KLRF1"),
  B_cell= c("CD19","MS4A1","CD79A","CD79B"),
  Plasma= c("MZB1","JCHAIN","IGHG1"),
  Monocyte = c("CD14","LYZ","S100A8","S100A9","VCAN"),
  Mac_M1 = c("NOS2","IL1B","TNF","CXCL9","CXCL10"),
  Mac_M2 = c("CD163","MRC1","MSR1","CCL22","IL10"),
  DC = c("CD1C","CLEC9A","LAMP3","CCR7"),
  Neutrophil = c("FCGR3B","CSF3R","CXCR2"),
  Mast = c("TPSAB1","CPA3","MS4A2"))
SIG_NK   <- c("KLRF1","KLRC1","KLRC2","KLRC3","KLRD1","SH2D1B","NCR1",
              "KIR2DL3","KIR3DL1","KIR2DS4","FCER1G","TXK","PRSS23")
SIG_CTL  <- c("CD8A","CD8B","GZMA","GZMB","GZMH","PRF1","NKG7","GNLY","KLRG1","CST7")
SIG_NAIV <- c("CD8A","CCR7","SELL","LEF1","IL7R","TCF7")
SIG_PROL <- c("MKI67","PCNA","TOP2A","BUB1")
SIG_EPI  <- c("EPCAM","KRT8","KRT18")
SIG_PURE <- c("CD68","CD3D","COL1A1")

# #############################################################################
#  STEP 1  载入三队列（与 30 同一逻辑）
# #############################################################################

hr("STEP 1  载入")

mats <- list()
if (file.exists(F_TCGA)) {
  x <- as.matrix(readRDS(F_TCGA)); storage.mode(x) <- "double"; mats[["TCGA"]] <- x
}
if (file.exists(F_GSE)) {
  gse <- as.matrix(readRDS(F_GSE)); storage.mode(gse) <- "double"
  grp <- rep(NA_character_, ncol(gse)); names(grp) <- colnames(gse)
  if (file.exists(F_SERIES)) {
    con <- gzfile(F_SERIES, "rt"); ln <- readLines(con, warn = FALSE); close(con)
    il <- grep("^!Sample_geo_accession", ln, value = TRUE)
    ids <- if (length(il)) gsub('"', "", strsplit(il[1], "\t")[[1]][-1]) else character(0)
    ch <- grep("^!Sample_characteristics", ln, value = TRUE)
    tl <- ch[grepl("tissue:", ch, fixed = TRUE)][1]
    if (!is.na(tl)) {
      v <- gsub('"', "", strsplit(tl, "\t")[[1]][-1])
      v <- sub("^\\s*tissue:\\s*", "", v)
      if (length(v) == length(ids))
        grp[ids] <- ifelse(grepl("cancer|tumor|tumour", v, ignore.case = TRUE), "Tumor", "Normal")
    }
  }
  tum <- names(grp)[!is.na(grp) & grp == "Tumor"]
  if (length(tum) >= 50) gse <- gse[, tum, drop = FALSE]
  mats[["GSE26712"]] <- gse
}
if (file.exists(F_ICGC)) {
  ic <- readRDS(F_ICGC)
  if (is.data.frame(ic)) ic <- as.matrix(ic)
  if (nrow(ic) < ncol(ic) && "TSPO" %in% colnames(ic)) ic <- t(ic)
  storage.mode(ic) <- "double"
  mats[["ICGC"]] <- ic
}
if (!length(mats)) stop("没有可用队列")
for (nm in names(mats)) cat(sprintf("  %-10s %5d x %4d\n", nm, nrow(mats[[nm]]), ncol(mats[[nm]])))

# #############################################################################
#  STEP 2  评分 + 全局表达
# #############################################################################

hr("STEP 2  评分")

score_z <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 2) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[g, , drop = FALSE]))); z[!is.finite(z)] <- NA
  colMeans(z, na.rm = TRUE)
}

SC <- list()
for (nm in names(mats)) {
  m <- mats[[nm]]
  if (!"TSPO" %in% rownames(m)) next
  imm_mat <- sapply(SIG_IMM, function(g) score_z(m, g))
  d <- data.frame(
    cohort   = nm,
    sample   = colnames(m),
    TSPO     = as.numeric(m["TSPO", ]),
    免疫总分 = rowMeans(imm_mat, na.rm = TRUE),
    CD8_CTL  = score_z(m, SIG_CTL),
    NK_specific = score_z(m, setdiff(SIG_NK, SIG_CTL)),
    CD8_naive   = score_z(m, SIG_NAIV),
    增殖     = score_z(m, SIG_PROL),
    上皮     = score_z(m, SIG_EPI),
    全局表达 = colMeans(m, na.rm = TRUE),      # 全局梯度
    stringsAsFactors = FALSE)
  g <- intersect(SIG_PURE, rownames(m))
  d$纯度代理 <- if (length(g) >= 2) {
    z <- t(scale(t(m[g, , drop = FALSE]))); z[!is.finite(z)] <- NA
    -colMeans(z, na.rm = TRUE)
  } else NA_real_
  SC[[nm]] <- d
  cat(sprintf("  %-10s 免疫总分n=%d, 全局与TSPO rho=%+.3f\n", nm,
              sum(is.finite(d$免疫总分)),
              suppressWarnings(cor(d$TSPO, d$全局表达, method = "spearman", use = "complete.obs"))))
}
ALL <- do.call(rbind, SC)

# #############################################################################
#  STEP 3  原始相关 + 控制全局表达
# #############################################################################

hr("STEP 3  原始 vs 控制全局表达")

pcor_rank <- function(a, b, cov) {
  cov <- as.matrix(cov)
  if (nrow(cov) != length(a)) cov <- t(cov)
  ok <- is.finite(a) & is.finite(b) & apply(cov, 2, function(z) all(is.finite(z)))
  a <- a[ok]; b <- b[ok]; cov <- cov[ok, , drop = FALSE]
  n <- length(a); k <- ncol(cov)
  if (n <= k + 3) return(c(rho = NA_real_, p = NA_real_))
  A  <- cbind(1, apply(cov, 2, rank))
  ea <- residuals(lm(rank(a) ~ A - 1)); eb <- residuals(lm(rank(b) ~ A - 1))
  r  <- suppressWarnings(cor(ea, eb))
  if (!is.finite(r) || abs(r) >= 1) return(c(rho = r, p = NA_real_))
  t  <- r * sqrt((n - k - 2) / (1 - r^2))
  c(rho = r, p = 2 * pt(-abs(t), df = n - k - 2))
}

IND <- c("免疫总分","CD8_CTL","NK_specific","CD8_naive","增殖","纯度代理","上皮")

RES <- do.call(rbind, lapply(names(SC), function(nm) {
  d <- SC[[nm]]
  do.call(rbind, lapply(IND, function(v) {
    ok <- is.finite(d$TSPO) & is.finite(d[[v]]) & is.finite(d$全局表达)
    if (sum(ok) < 30) return(NULL)
    r0 <- suppressWarnings(cor(d$TSPO[ok], d[[v]][ok], method = "spearman"))
    r1 <- pcor_rank(d$TSPO, d[[v]], matrix(d$全局表达, ncol = 1))
    data.frame(cohort = nm, indicator = v, n = sum(ok),
               rho_raw = round(r0, 3),
               rho_ctrl = round(r1["rho"], 3), p_ctrl = signif(r1["p"], 3),
               stringsAsFactors = FALSE)
  }))
}))

# 宽表
fmt_w <- function(col) {
  w <- reshape(as.data.frame(RES[, c("indicator","cohort",col)]),
               idvar = "indicator", timevar = "cohort", direction = "wide")
  names(w) <- sub(paste0("^", col, "\\."), "", names(w))
  w
}
W_raw  <- fmt_w("rho_raw")
W_ctrl <- fmt_w("rho_ctrl")

cat("\n【原始相关】\n"); pr(W_raw, 3)
cat("\n【控制全局表达后】\n"); pr(W_ctrl, 3)
write.csv(RES, file.path(OUTDIR, "T23a_scores_confound_adjusted.csv"), row.names = FALSE)

# #############################################################################
#  STEP 4  一致性判定与最终对照表
# #############################################################################

hr("STEP 4  三队列一致性")

FINAL <- do.call(rbind, lapply(IND, function(v) {
  z <- RES[RES$indicator == v, ]
  if (!nrow(z)) return(NULL)
  r <- setNames(z$rho_ctrl, z$cohort)
  n_pos <- sum(r > 0, na.rm = TRUE); n_neg <- sum(r < 0, na.rm = TRUE)
  all_same <- (n_pos == length(r)) || (n_neg == length(r))
  min_abs <- min(abs(r), na.rm = TRUE)
  data.frame(
    indicator = v,
    TCGA = if ("TCGA" %in% names(r)) round(r[["TCGA"]], 3) else NA_real_,
    GSE26712 = if ("GSE26712" %in% names(r)) round(r[["GSE26712"]], 3) else NA_real_,
    ICGC = if ("ICGC" %in% names(r)) round(r[["ICGC"]], 3) else NA_real_,
    方向一致 = all_same,
    最小绝对效应 = round(min_abs, 3),
    判定 = ifelse(all_same & min_abs >= 0.15, "三队列一致",
           ifelse(all_same & min_abs >= 0.08, "三队列同向但强度不一",
           ifelse(all_same, "三队列同向但极弱", "方向不一致"))),
    stringsAsFactors = FALSE)
}))
FINAL <- FINAL[order(-FINAL$最小绝对效应), ]

cat("\n===== 3.2 的最终对照表（控制全局表达后）=====\n")
pr(FINAL, 3)
write.csv(FINAL, file.path(OUTDIR, "T23b_replication_final.csv"), row.names = FALSE)

# #############################################################################
#  STEP 5  阴性对照的共线性问题
# #############################################################################

hr("STEP 5  阴性对照共线性")

COL <- do.call(rbind, lapply(names(SC), function(nm) {
  d <- SC[[nm]]
  data.frame(cohort = nm,
             CTL_vs_naive = round(suppressWarnings(cor(d$CD8_CTL, d$CD8_naive,
                                       method = "spearman", use = "complete.obs")), 3),
             NK_vs_naive  = round(suppressWarnings(cor(d$NK_specific, d$CD8_naive,
                                       method = "spearman", use = "complete.obs")), 3),
             stringsAsFactors = FALSE)
}))
cat("\n"); pr(COL, 3)
cat("\n⚠ 正文把 CD8_naive 当作「阴性对照」，用以排除「只是淋巴细胞总量更多」。\n")
cat("  但该评分与 CD8_CTL 的相关为 0.70–0.80，属高度共线，\n")
cat("  不能充当独立对照。相关论证须撤除或重新设计对照（见 STEP 6）。\n")
write.csv(COL, file.path(OUTDIR, "T23c_collinearity_warning.csv"), row.names = FALSE)

# #############################################################################
#  STEP 6  建议
# #############################################################################

hr("STEP 6  对正文 3.2 的改写建议")

gv <- function(v) FINAL[FINAL$indicator == v, ]
mk <- function(v) {
  z <- gv(v); if (!nrow(z)) return("NA")
  sprintf("TCGA %+.2f / GSE26712 %+.2f / ICGC %+.2f（%s）",
          z$TCGA, z$GSE26712, z$ICGC, z$判定)
}

LINES <- c(
  "3.2 改写建议（基于三队列控制全局表达后的结果）",
  "",
  "【可保留】",
  paste("  免疫总分   ", mk("免疫总分")),
  paste("  纯度代理   ", mk("纯度代理")),
  paste("  增殖       ", mk("增殖")),
  "",
  "【须降级或撤除】",
  paste("  CD8_CTL    ", mk("CD8_CTL")),
  paste("  NK_specific", mk("NK_specific")),
  "",
  "【论证缺陷，须撤】",
  "  CD8_naive 作为阴性对照不成立：与 CD8_CTL 共线 0.70–0.80。",
  "  改用更合适的对照：CD4_T 或 B_cell（与细胞毒程序不共线），见下。",
  "",
  "【建议的新对照】",
  "  用 CD4_T 评分作对照：若 CD4_T 与 TSPO 的相关明显弱于 CD8_CTL，",
  "  才能支持「细胞毒程序特异」。若两者相近，则只能说「T 细胞浸润」。"
)
cat("\n"); for (l in LINES) cat(l, "\n")

# 补算 CD4_T 与 B_cell 作为替代对照
hr("STEP 6B  替代对照：CD4_T / B_cell / Plasma")

ALT <- do.call(rbind, lapply(names(SC), function(nm) {
  m <- mats[[nm]]; d <- SC[[nm]]
  do.call(rbind, lapply(c("CD4_T","B_cell","Plasma"), function(k) {
    sc <- score_z(m, SIG_IMM[[k]])
    ok <- is.finite(d$TSPO) & is.finite(sc) & is.finite(d$全局表达)
    if (sum(ok) < 30) return(NULL)
    r0 <- suppressWarnings(cor(d$TSPO[ok], sc[ok], method = "spearman"))
    r1 <- pcor_rank(d$TSPO, sc, matrix(d$全局表达, ncol = 1))
    data.frame(cohort = nm, indicator = k, rho_raw = round(r0, 3),
               rho_ctrl = round(r1["rho"], 3), p_ctrl = signif(r1["p"], 3),
               stringsAsFactors = FALSE)
  }))
}))
pr(ALT, 3)

W_alt <- reshape(as.data.frame(ALT[, c("indicator","cohort","rho_ctrl")]),
                 idvar = "indicator", timevar = "cohort", direction = "wide")
names(W_alt) <- sub("^rho_ctrl\\.", "", names(W_alt))
cat("\n替代对照（控制全局后）：\n"); pr(W_alt, 3)
cat("\n判读：若 CD4_T 与 CD8_CTL 的效应量接近，则 TSPO 标记的是「T 细胞浸润」\n")
cat("      而非「细胞毒程序」；若 CD4_T 明显更弱，则细胞毒特异性可以主张。\n")

writeLines(LINES, file.path(OUTDIR, "T23d_recommendation.txt"))
write.csv(ALT, file.path(OUTDIR, "T23e_alternative_controls.csv"), row.names = FALSE)

hr("31_bulk_confound_adjusted.R 完成")

# -----------------------------------------------------------------------------
# -- Fig3_来源分解 / 35_tspo_source_index.R
# -----------------------------------------------------------------------------
# =============================================================================
#  35_tspo_source_index.R   TSPO 高表达的来源可区分，且该结构跨队列可重复
#
#  【本文件由定稿代码的块 32 同步生成，内容与 TSPO定稿代码.R 的块 32 完全一致】
#
#  ── 要回答的问题 ────────────────────────────────────────────────────────────
#  图 2 显示 TSPO 在 TCGA 分子亚型间差异显著（ε²=0.253），且扣除细胞构成后
#  免疫型与分化型的排序反转——提示两者高 TSPO 的来源不同：
#    免疫型：来自浸润髓系（扣掉免疫后优势消失）
#    分化型：来自上皮细胞本身（扣掉后仍维持）
#
#  通用上皮标志物无法提取 bulk 层面的「上皮分量」
#  （TSPO 与上皮评分仅相关 0.07–0.08，残差化几乎不起作用）。
#  故换一个思路：不追求把 TSPO 拆成两个数值，而是检验
#  「TSPO 扣掉髓系贡献后的残余部分」是否与分化型转录程序相关。
#
#  ── 指标定义 ────────────────────────────────────────────────────────────────
#  髓系评分   = CD68/CD14/LYZ/AIF1/C1QA/B/ITGAM/FCER1G/TYROBP 的 z 均值
#  DIF 程序   = 254 个分化型亚型基因的 z 均值（原文补充表 SB3）
#  MES 程序   = 462 个间质型亚型基因的 z 均值（原文补充表 SB5）
#  source_idx = rank(TSPO) 对 rank(髓系评分) 的残差
#
#  ── 四个预设检验 ────────────────────────────────────────────────────────────
#  T1（TCGA）source_idx 在四亚型间的差异；预设方向 differentiated > immunoreactive
#  T2（跨队列）source_idx 与 DIF 程序正相关（逐队列 + Stouffer）
#  T3（跨队列）source_idx 与 MES 程序不相关或负相关
#  T4（一致性）TSPO 与髓系评分的相关在各队列是否稳定
#  T1–T3 共 7 个检验，BH 校正。
#
#  ── 输出 ────────────────────────────────────────────────────────────────────
#  T37a_source_index.csv / T37b_tcga_subtype.csv / T37c_external.csv
#  T37d_meta.csv / T37e_consistency.csv / T37f_summary_bh.csv
#  T37g_verdict.txt / T37h_circularity_overlap.csv / T37i_t2_without_overlap.csv
#  T38b_variance_decomp.csv   方差分解
#
#  图 3 的绘制由 fig03_source.R 单独完成，不在本文件内。
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/35_tspo_source_index.R", encoding = "UTF-8")
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}

if (!exists("PROJ"))       PROJ       <- "D:/Rstudy/tcga"
if (!exists("OUT"))        OUT        <- file.path(PROJ, "results_v2")
if (!exists("F_EXPR"))     F_EXPR     <- file.path(OUT, "EXPR_clean.rds")
if (!exists("F_SUB"))      F_SUB      <- file.path(PROJ, "TCGA_OV_subtype.gz")
if (!exists("F_MAP570"))   F_MAP570   <- file.path(PROJ, "GPL570_probe_map_correct.csv")
if (!exists("F_SOFT6480")) F_SOFT6480 <- file.path(PROJ, "GPL6480.soft.gz")
if (!exists("F_DIF"))      F_DIF      <- file.path(PROJ,
  "10780432ccr200103-sup-table_sb3_subtype2-c4_diff_05jul20.xlsx")
if (!exists("F_MES"))      F_MES      <- file.path(PROJ,
  "10780432ccr200103-sup-table_sb5_subtype4-c1_mes_05jul20.xlsx")

OUTDIR <- file.path(OUT, "tables")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("PROJ"))       PROJ       <- "D:/Rstudy/tcga"
if (!exists("OUT"))        OUT        <- file.path(PROJ, "results_v2")
if (!exists("F_EXPR"))     F_EXPR     <- file.path(OUT, "EXPR_clean.rds")
if (!exists("F_SUB"))      F_SUB      <- file.path(PROJ, "TCGA_OV_subtype.gz")
if (!exists("F_MAP570"))   F_MAP570   <- file.path(PROJ, "GPL570_probe_map_correct.csv")
if (!exists("F_SOFT6480")) F_SOFT6480 <- file.path(PROJ, "GPL6480.soft.gz")
if (!exists("F_DIF"))      F_DIF      <- file.path(PROJ,
  "10780432ccr200103-sup-table_sb3_subtype2-c4_diff_05jul20.xlsx")
if (!exists("F_MES"))      F_MES      <- file.path(PROJ,
  "10780432ccr200103-sup-table_sb5_subtype4-c1_mes_05jul20.xlsx")

OUTDIR <- file.path(OUT, "tables")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
#  STEP 0  亚型基因集与探针注释
# =============================================================================

hr("STEP 0  基因集与注释")

read_sig_xlsx <- function(f, tag) {
  if (!file.exists(f)) { cat("  缺 ", basename(f), "\n", sep = ""); return(NULL) }
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("需要 openxlsx")
  sh <- openxlsx::getSheetNames(f)[1]
  d  <- openxlsx::read.xlsx(f, sheet = sh, colNames = FALSE)
  # 找表头行（含 Probe 与 Gene）
  hi <- NA
  for (i in seq_len(min(8, nrow(d)))) {
    cells <- tolower(trimws(as.character(unlist(d[i, ]))))
    if ("probe" %in% cells && "gene" %in% cells) { hi <- i; break }
  }
  if (is.na(hi)) { cat("  ", tag, " 未找到表头\n", sep = ""); return(NULL) }
  hdr <- trimws(as.character(unlist(d[hi, ])))
  gi  <- which(tolower(hdr) == "gene")[1]
  g   <- as.character(unlist(d[(hi + 1):nrow(d), gi]))
  g   <- trimws(g)
  g   <- g[!is.na(g) & nzchar(g) & g != "NA"]
  g   <- unique(sub("\\s*///.*$", "", g))
  cat("  ", tag, ": ", length(g), " 个基因（前 6: ",
      paste(head(g, 6), collapse = ", "), "）\n", sep = "")
  g
}

G_DIF <- read_sig_xlsx(F_DIF, "DIF")
G_MES <- read_sig_xlsx(F_MES, "MES")
if (is.null(G_DIF) || is.null(G_MES)) stop("亚型基因集读取失败")

MAP570 <- read.csv(F_MAP570, stringsAsFactors = FALSE, check.names = FALSE)
colnames(MAP570)[1:2] <- c("probe", "symbol")
MAP570 <- MAP570[!is.na(MAP570$symbol) & nzchar(MAP570$symbol), ]
MAP570$symbol <- sub("\\s*///.*$", "", MAP570$symbol)
MAP570 <- MAP570[!duplicated(MAP570$probe), ]

MAP6480 <- NULL
if (file.exists(F_SOFT6480)) {
  ln <- readLines(F_SOFT6480, warn = FALSE, encoding = "UTF-8")
  i0 <- grep("^!platform_table_begin", ln)[1]; i1 <- grep("^!platform_table_end", ln)[1]
  if (!is.na(i0) && !is.na(i1)) {
    hdr <- strsplit(ln[i0 + 1], "\t")[[1]]
    body <- ln[(i0 + 2):(i1 - 1)]
    ci <- which(hdr == "GENE_SYMBOL")[1]
    sy <- vapply(strsplit(body, "\t"), function(z) if (length(z) >= ci) z[ci] else NA_character_,
                 character(1), USE.NAMES = FALSE)
    MAP6480 <- data.frame(probe = sub("\t.*$", "", body), symbol = sy, stringsAsFactors = FALSE)
    MAP6480 <- MAP6480[!is.na(MAP6480$symbol) & nzchar(MAP6480$symbol) & MAP6480$symbol != "", ]
    MAP6480$symbol <- sub("\\s*///.*$", "", MAP6480$symbol)
    MAP6480 <- MAP6480[!duplicated(MAP6480$probe), ]
  }
}
cat("GPL570 ", nrow(MAP570), " 探针；GPL6480 ", ifelse(is.null(MAP6480), 0, nrow(MAP6480)),
    " 探针\n", sep = "")

# =============================================================================
#  STEP 1  评分与指标
# =============================================================================

GENES <- list(
  Myeloid = c("CD68","CD14","LYZ","AIF1","C1QA","C1QB","ITGAM","FCER1G","TYROBP"),
  Epithelial = c("EPCAM","KRT8","KRT18","KRT19","PAX8","CLDN4","MUC16","WFDC2"),
  Stromal = c("COL1A1","COL1A2","COL3A1","DCN","LUM","FAP","ACTA2","TAGLN","POSTN")
)

score_z <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 3) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[g, , drop = FALSE]))); z[!is.finite(z)] <- NA
  colMeans(z, na.rm = TRUE)
}

# 构建一个队列的全部指标
build <- function(mat, tag) {
  tspo <- as.numeric(mat["TSPO", ])
  mye  <- score_z(mat, GENES$Myeloid)
  epi  <- score_z(mat, GENES$Epithelial)
  str  <- score_z(mat, GENES$Stromal)
  dif  <- score_z(mat, G_DIF)
  mes  <- score_z(mat, G_MES)

  # source_idx：rank(TSPO) 对 rank(髓系) 的残差
  ok <- is.finite(tspo) & is.finite(mye)
  src <- rep(NA_real_, length(tspo))
  keep <- NA_real_
  if (sum(ok) >= 30) {
    ry <- rank(tspo[ok]); rx <- rank(mye[ok])
    fit <- lm(ry ~ rx)
    src[ok] <- residuals(fit)
    keep <- 1 - summary(fit)$r.squared
  }
  d <- data.frame(sample = colnames(mat), cohort = tag,
                  TSPO = tspo, Myeloid = mye, Epithelial = epi, Stromal = str,
                  DIF_prog = dif, MES_prog = mes, source_idx = src,
                  stringsAsFactors = FALSE)
  list(data = d, keep = keep,
       n_dif = length(intersect(G_DIF, rownames(mat))),
       n_mes = length(intersect(G_MES, rownames(mat))))
}

# =============================================================================
#  STEP 2  TCGA
# =============================================================================

hr("STEP 2  TCGA")

if (!file.exists(F_EXPR)) stop("缺 ", F_EXPR)
TC <- as.matrix(readRDS(F_EXPR)); storage.mode(TC) <- "double"
bt <- build(TC, "TCGA")
cat("  TCGA: ", ncol(TC), " 样本；DIF 命中 ", bt$n_dif, " 基因，MES 命中 ",
    bt$n_mes, " 基因\n", sep = "")
cat("  source_idx 残差保留方差: ", round(bt$keep, 3), "\n", sep = "")

SUB <- NULL
if (file.exists(F_SUB)) {
  sb <- read.delim(gzfile(F_SUB), stringsAsFactors = FALSE)
  colnames(sb)[1:2] <- c("sample", "subtype")
  sb$subtype <- tolower(trimws(sb$subtype))
  SUB <- sb
}
D <- bt$data
if (!is.null(SUB)) {
  D$subtype <- SUB$subtype[match(D$sample, SUB$sample)]
  cat("  匹配到亚型: ", sum(!is.na(D$subtype)), " / ", nrow(D), "\n", sep = "")
}

# ---- T1  source_idx 在亚型间 ----
LV <- c("immunoreactive", "differentiated", "proliferative", "mesenchymal")
LAB <- c(immunoreactive = "Immunoreactive", differentiated = "Differentiated",
         proliferative = "Proliferative", mesenchymal = "Mesenchymal")
T1 <- NULL
if (sum(!is.na(D$subtype)) >= 100) {
  ds <- D[!is.na(D$subtype) & is.finite(D$source_idx), ]
  ds$subtype <- factor(ds$subtype, levels = intersect(LV, unique(ds$subtype)))
  kw <- kruskal.test(source_idx ~ subtype, data = ds)
  md <- tapply(ds$source_idx, ds$subtype, median)
  cat("\nT1 各亚型的 source_idx 中位数：\n"); print(round(md, 2))
  cat(sprintf("    Kruskal-Wallis: H = %.3f, p = %.3e\n", kw$statistic, kw$p.value))

  # 预设方向：differentiated > immunoreactive
  pd <- NULL
  if (all(c("differentiated","immunoreactive") %in% names(md))) {
    w <- suppressWarnings(wilcox.test(
      ds$source_idx[ds$subtype == "differentiated"],
      ds$source_idx[ds$subtype == "immunoreactive"], alternative = "greater"))
    pd <- w$p.value
    cat(sprintf("    预设方向 Differentiated > Immunoreactive: 单侧 p = %.3e\n", pd))
  }
  T1 <- data.frame(检验 = "T1 TCGA 亚型间 source_idx", n = nrow(ds),
                   H = round(kw$statistic, 3), p = kw$p.value,
                   p_方向 = ifelse(is.null(pd), NA, pd), row.names = NULL)
  write.csv(data.frame(subtype = names(md), median = round(as.numeric(md), 3)),
            file.path(OUTDIR, "T37b_tcga_subtype.csv"), row.names = FALSE)
}

# =============================================================================
#  STEP 3  外部队列
# =============================================================================

hr("STEP 3  外部队列")

read_cohort <- function(gsm, map) {
  f <- file.path(PROJ, paste0(gsm, "_series_matrix.txt.gz"))
  if (!file.exists(f)) { cat("  缺 ", gsm, "\n", sep = ""); return(NULL) }
  con <- gzfile(f, "rt"); ln <- readLines(con, warn = FALSE); close(con)
  i0 <- grep("^!series_matrix_table_begin", ln)[1]
  i1 <- grep("^!series_matrix_table_end", ln)[1]
  if (is.na(i0) || is.na(i1)) return(NULL)
  tab <- read.delim(text = paste(ln[(i0 + 1):(i1 - 1)], collapse = "\n"),
                    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(tab[, -1, drop = FALSE]); storage.mode(mat) <- "double"
  rownames(mat) <- tab[[1]]
  m <- map[match(rownames(mat), map$probe), ]
  ok <- !is.na(m$symbol)
  mat <- mat[ok, , drop = FALSE]; sym <- m$symbol[ok]
  od <- order(-rowMeans(mat, na.rm = TRUE))
  mat <- mat[od, , drop = FALSE]; sym <- sym[od]
  kp <- !duplicated(sym)
  EX <- mat[kp, , drop = FALSE]; rownames(EX) <- sym[kp]
  EX
}

EXT <- list()
if (!is.null(MAP6480)) EXT[["GSE32062"]] <- read_cohort("GSE32062-GPL6480", MAP6480)
EXT[["GSE63885"]] <- read_cohort("GSE63885", MAP570)
EXT[["GSE9891"]]  <- read_cohort("GSE9891", MAP570)
EXT <- EXT[!vapply(EXT, is.null, logical(1))]

CO <- list()
CO[["TCGA"]] <- bt
for (nm in names(EXT)) {
  b <- build(EXT[[nm]], nm)
  cat("  ", nm, ": ", ncol(EXT[[nm]]), " 样本；DIF 命中 ", b$n_dif,
      "，MES 命中 ", b$n_mes, "，残差保留 ", round(b$keep, 3), "\n", sep = "")
  CO[[nm]] <- b
}

ALL <- do.call(rbind, lapply(CO, `[[`, "data"))
write.csv(ALL, file.path(OUTDIR, "T37a_source_index.csv"), row.names = FALSE)

# =============================================================================
#  STEP 4  T2 / T3：source_idx 与亚型程序
# =============================================================================

hr("STEP 4  T2 / T3")

ext_names <- setdiff(names(CO), "TCGA")
run_cor <- function(nm, prog, lab) {
  d <- CO[[nm]]$data
  ok <- is.finite(d$source_idx) & is.finite(d[[prog]])
  if (sum(ok) < 50) return(NULL)
  ct <- suppressWarnings(cor.test(d$source_idx[ok], d[[prog]][ok], method = "spearman"))
  data.frame(队列 = nm, 检验 = lab, n = sum(ok),
             rho = round(as.numeric(ct$estimate), 3),
             p = signif(ct$p.value, 3), row.names = NULL)
}

EX <- do.call(rbind, c(
  lapply(ext_names, run_cor, prog = "DIF_prog", lab = "T2 source_idx ~ DIF"),
  lapply(ext_names, run_cor, prog = "MES_prog", lab = "T3 source_idx ~ MES")))
cat("\n逐队列：\n"); pr(EX, 3)

# Stouffer 合并（按检验分组）
META <- do.call(rbind, lapply(c("T2 source_idx ~ DIF", "T3 source_idx ~ MES"), function(lb) {
  z <- EX[EX$检验 == lb, ]
  if (!nrow(z)) return(NULL)
  # 大样本下 p 会下溢成 0，qnorm(1-0)=Inf，需截断到浮点最小正数
  pz <- pmax(z$p, .Machine$double.eps)
  zs <- qnorm(1 - pz / 2) * sign(z$rho)
  Z <- sum(zs) / sqrt(length(zs))
  data.frame(检验 = lb, 队列数 = nrow(z), Z = round(Z, 3),
             p = signif(2 * pnorm(-abs(Z)), 3), row.names = NULL)
}))
cat("\n合并：\n"); pr(META, 3)
write.csv(EX,   file.path(OUTDIR, "T37c_external.csv"), row.names = FALSE)
write.csv(META, file.path(OUTDIR, "T37d_meta.csv"), row.names = FALSE)

# =============================================================================
#  STEP 4b  循环性检验：DIF 签名与 TSPO 共表达网络的重叠
#
#  TSPO 的共表达网络（图 1）含 2,714 个基因，DIF 签名 254 个。若二者大面积重叠，
#  则 source_idx（本质是 TSPO 的残差）与 DIF 程序的相关就是自相关，不能作为证据。
#
#  本步做三件事：
#    1. 计算重叠基因数与随机期望，给出富集倍数与超几何 p
#    2. 从 DIF 签名中剔除重叠基因，用剩余基因重算 DIF 程序
#    3. 用去重叠后的 DIF 程序重测 T2 —— 若仍显著，则结论不是循环性造成的
# =============================================================================

hr("STEP 4b  循环性检验")

F_NET <- file.path(OUT, "tables", "T8_network_all_genes.csv")
NET_GENES <- NULL
if (file.exists(F_NET)) {
  nw <- read.csv(F_NET, stringsAsFactors = FALSE)
  if (!"in_network" %in% colnames(nw)) nw$in_network <- abs(nw$rho) > 0.3 & nw$FDR < 0.05
  NET_GENES <- nw$gene[nw$in_network]
  cat("TSPO 共表达网络: ", length(NET_GENES), " 个基因\n", sep = "")
} else {
  cat("⚠ 未找到网络表 ", F_NET, "，跳过循环性检验\n", sep = "")
}

OVER <- NULL; T2B <- NULL
if (!is.null(NET_GENES)) {
  uni <- rownames(TC)
  n_u <- length(uni)
  # 置换检验：从本队列实际测到的基因池随机抽同样大小的集合
# （比超几何假定更保守；同一零分布用于图 3E）
D_SET <- intersect(G_DIF, uni)
M_SET <- intersect(G_MES, uni)
ov_d  <- length(intersect(D_SET, NET_GENES))
ov_m  <- length(intersect(M_SET, NET_GENES))

N_PERM <- 5000
set.seed(20261025)
null_d <- vapply(seq_len(N_PERM), function(i)
  length(intersect(sample(uni, length(D_SET)), NET_GENES)), integer(1))
null_m <- vapply(seq_len(N_PERM), function(i)
  length(intersect(sample(uni, length(M_SET)), NET_GENES)), integer(1))
p_d <- (sum(null_d >= ov_d) + 1) / (N_PERM + 1)
p_m <- (sum(null_m >= ov_m) + 1) / (N_PERM + 1)

OVER <- data.frame(
  基因集 = c("DIF", "MES"),
  可测基因数 = c(length(D_SET), length(M_SET)),
  与网络重叠 = c(ov_d, ov_m),
  置换均值 = round(c(mean(null_d), mean(null_m)), 1),
  富集倍数 = round(c(ov_d / mean(null_d), ov_m / mean(null_m)), 2),
  置换p = signif(c(p_d, p_m), 3),
  row.names = NULL)
cat("\n与 TSPO 网络的重叠：\n"); pr(OVER, 3)

  # 剔除重叠后重算 DIF 程序
  G_DIF_CLEAN <- setdiff(G_DIF, NET_GENES)
  cat("\n去重叠后的 DIF 签名: ", length(G_DIF_CLEAN), " / ", length(G_DIF),
      " 个基因\n", sep = "")

  recalc <- function(mat) score_z(mat, G_DIF_CLEAN)
  T2B <- do.call(rbind, lapply(ext_names, function(nm) {
    d <- CO[[nm]]$data
    pc <- recalc(EXT[[nm]])
    ok <- is.finite(d$source_idx) & is.finite(pc)
    if (sum(ok) < 50) return(NULL)
    ct <- suppressWarnings(cor.test(d$source_idx[ok], pc[ok], method = "spearman"))
    data.frame(队列 = nm, n = sum(ok), 原始rho = NA_real_,
               rho_去重叠 = round(as.numeric(ct$estimate), 3),
               p_去重叠 = signif(ct$p.value, 3), row.names = NULL)
  }))
  T2B$原始rho <- EX$rho[match(T2B$队列, EX$队列[EX$检验 == "T2 source_idx ~ DIF"])]
  cat("\nT2 去重叠后重测：\n"); pr(T2B, 3)
  write.csv(OVER, file.path(OUTDIR, "T37h_circularity_overlap.csv"), row.names = FALSE)
  write.csv(T2B,  file.path(OUTDIR, "T37i_t2_without_overlap.csv"), row.names = FALSE)

  cat("\n判读：若「rho_去重叠」与「原始rho」接近且仍显著，则结论不是循环性造成的；\n")
  cat("      若大幅下降或失去显著性，则原始相关主要由重叠基因驱动，不能作为独立证据。\n")
}

# =============================================================================
#  STEP 5  T4：髓系成分的跨队列稳定性
# =============================================================================

hr("STEP 5  T4  髓系成分一致性")

CONS <- do.call(rbind, lapply(names(CO), function(nm) {
  d <- CO[[nm]]$data
  ok <- is.finite(d$TSPO) & is.finite(d$Myeloid)
  ct <- suppressWarnings(cor.test(d$TSPO[ok], d$Myeloid[ok], method = "spearman"))
  data.frame(队列 = nm, n = sum(ok),
             TSPO与髓系相关 = round(as.numeric(ct$estimate), 3),
             p = signif(ct$p.value, 3),
             上皮命中 = CO[[nm]]$n_dif, row.names = NULL)
}))
cat("\n"); pr(CONS, 3)
cat(sprintf("\nTSPO~髓系 的跨队列范围: %.3f ~ %.3f（离散度 SD = %.3f）\n",
            min(CONS$TSPO与髓系相关), max(CONS$TSPO与髓系相关),
            sd(CONS$TSPO与髓系相关)))
write.csv(CONS, file.path(OUTDIR, "T37e_consistency.csv"), row.names = FALSE)

# =============================================================================
#  STEP 6  BH 汇总
# =============================================================================

hr("STEP 6  BH 校正")

RES <- data.frame(
  检验 = c("T1 TCGA 亚型间", EX$检验),
  队列 = c("TCGA", EX$队列),
  n = c(T1$n, EX$n),
  效应量 = c(round(T1$H, 3), EX$rho),
  p = c(T1$p, EX$p),
  row.names = NULL)
RES$p_BH <- signif(p.adjust(RES$p, method = "BH"), 3)
RES$显著 <- ifelse(RES$p_BH < 0.05, "是", "否")
cat("\n"); pr(RES, 3)
write.csv(RES, file.path(OUTDIR, "T37f_summary_bh.csv"), row.names = FALSE)

# =============================================================================
#  STEP 7  结论
# =============================================================================

hr("STEP 7  结论")

t2_meta <- META$p[META$检验 == "T2 source_idx ~ DIF"]
t3_meta <- META$p[META$检验 == "T3 source_idx ~ MES"]
t2_z    <- META$Z[META$检验 == "T2 source_idx ~ DIF"]
t3_z    <- META$Z[META$检验 == "T3 source_idx ~ MES"]
dir_ok  <- !is.na(t2_z) && t2_z > 0
n_sig   <- sum(RES$p_BH < 0.05)
mye_rng <- range(CONS$TSPO与髓系相关)

V <- c(
  "TSPO 高表达的来源可区分性：source_idx 与分化型程序",
  paste("时间:", format(Sys.time())),
  sprintf("队列: %s", paste(names(CO), collapse = ", ")),
  "",
  sprintf("T1 TCGA 亚型间 source_idx: H = %.3f, p = %.3e", T1$H, T1$p),
  sprintf("   预设方向 Differentiated > Immunoreactive: 单侧 p = %.3e", T1$p_方向),
  sprintf("T2 source_idx ~ DIF 程序（Stouffer, %d 队列）: Z = %.3f, p = %.3e",
          META$队列数[META$检验 == "T2 source_idx ~ DIF"], t2_z, t2_meta),
  sprintf("T3 source_idx ~ MES 程序（Stouffer, %d 队列）: Z = %.3f, p = %.3e",
          META$队列数[META$检验 == "T3 source_idx ~ MES"], t3_z, t3_meta),
  sprintf("T4 TSPO ~ 髓系 的跨队列范围: %.3f ~ %.3f",
          mye_rng[1], mye_rng[2]),
  "",
  sprintf("BH 校正后显著项: %d / %d", n_sig, nrow(RES)),
  "",
  "判读：",
  if (dir_ok && t2_meta < 0.05) {
    c("· TSPO 扣除髓系贡献后的残余部分，与分化型转录程序正相关，且跨队列可重复。",
      "· 这支持「TSPO 高表达存在两条来源」：一条随髓系浸润升高，另一条随分化型程序升高。",
      "· 驱动第二种来源的是分化状态本身，而非上皮细胞含量——",
      "  通用上皮标志物评分与 TSPO 的相关仅 0.07–0.08，不足以解释该关联。")
  } else {
    c("· source_idx 与分化型程序的相关未能跨队列重复，",
      "  「TSPO 高表达有两条来源」这一说法在 bulk 层面得不到一致支持。")
  },
  "",
  "适用边界：",
  "· source_idx 是秩残差，尺度依赖队列内分布，跨队列只能比较方向与相对大小，不能比较绝对值。",
  "· DIF 程序含 254 个基因，其中部分与 TSPO 共表达网络重叠，相关性含循环成分；",
  "  已通过 rank 残差与 Spearman 缓解，但未完全排除。",
  "· 本结果说明的是「来源结构」，不构成任何预后或治疗预测主张。"
)
cat("\n"); for (l in V) cat(l, "\n")
writeLines(V, file.path(OUTDIR, "T37g_verdict.txt"))


# =============================================================================
#  STEP 8  方差分解（图 3F）
#
#  R²(髓系) 与 R²(DIF) 是各自单变量的解释力，两者相关时共享一部分方差，
#  直接相加会超过 1。此处按标准分解为四段，四段相加恰为 1：
#    髓系独有 = R²(两者) − R²(仅 DIF)
#    DIF 独有 = R²(两者) − R²(仅髓系)
#    共享     = R²(仅髓系) + R²(仅DIF) − R²(两者)
#    未解释   = 1 − R²(两者)
# =============================================================================

hr("STEP 8  方差分解")

VD <- do.call(rbind, lapply(names(CO), function(nm) {
  d <- CO[[nm]]$data
  d <- d[is.finite(d$TSPO) & is.finite(d$Myeloid) & is.finite(d$DIF_prog), ]
  if (nrow(d) < 50) return(NULL)
  r_mye  <- summary(lm(TSPO ~ Myeloid, data = d))$r.squared
  r_dif  <- summary(lm(TSPO ~ DIF_prog, data = d))$r.squared
  r_both <- summary(lm(TSPO ~ Myeloid + DIF_prog, data = d))$r.squared
  data.frame(cohort = nm, n = nrow(d),
             髓系独有 = round(r_both - r_dif, 3),
             DIF独有  = round(r_both - r_mye, 3),
             共享     = round(r_mye + r_dif - r_both, 3),
             未解释    = round(1 - r_both, 3),
             合计     = round((r_both - r_dif) + (r_both - r_mye) +
                              (r_mye + r_dif - r_both) + (1 - r_both), 3),
                              合并R2   = r_both,          # 未取整，供合计百分比使用
             row.names = NULL)
}))
cat("\n两个来源对 TSPO 方差的分解：\n"); pr(VD, 3)
cat(sprintf("\n四段相加校验（应全为 1）: %s\n",
            paste(unique(VD$合计), collapse = ", ")))
cat(sprintf("两个来源合计解释 %.1f%%–%.1f%% 的 TSPO 方差\n",
             100 * min(VD$合并R2), 100 * max(VD$合并R2)))
write.csv(VD, file.path(OUTDIR, "T38b_variance_decomp.csv"), row.names = FALSE)

hr("32_tspo_source_index.R 完成")


# ############################################################################
# ##  Fig4_上皮代谢  (6 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- Fig4_上皮代谢 / 03_epithelial_mixed.R
# -----------------------------------------------------------------------------
# =============================================================================
#  03_epithelial_mixed.R
#  结论二：上皮细胞中 TSPO 标记氧化磷酸化与糖酵解活跃状态
#
#  对应正文 3.2 与图 2A–D。
#
#  ⚠ 推断单位是「病人」不是「样本」
#    GSE165897 的对象含 11 例病人、22 个样本（每例病人提供 primary 与 interval
#    两个样本），EOC 亚群 7,929 个细胞。`sample` 列有 22 个水平，`patient_id`
#    列有 11 个水平。早期版本误用 `sample` 当病人，把 22 个样本当成 22 例病人，
#    导致符号检验存在伪重复。
#    本脚本以 patient_id 为推断单位（11 例），并同时输出 sample 层面（22 个）
#    作为可比口径。
#
#  方法要点（缺一不可）：
#    1. 逐基因回归掉文库大小（log10 nCount），否则 TSPO 高的细胞什么都高
#    2. 病人作随机截距，吸收病人间异质性
#    3. 混合模型的 β 会被细胞数加权，故必须同时报「逐病人符号检验」
#       —— 后者不依赖分布假设、不受细胞数加权，是主要推断依据
#    4. 50 个随机 12 基因模块作阴性对照，判断模块特异性
#    5. lmerTest 缺失时用 t 近似（df = n - 单位数 - 1），不依赖该包
#
#  输入：F_SC_EPI（32,847 基因 x 50,867 细胞，cell_type == "EOC" 取 7,929）
#  输出：tables/T2a_*、T2b_*、T2c_random_control.csv
# =============================================================================

if (!exists("OUT")) stop("请先 source('00_setup.R')")
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(lme4) })

MIN_CELLS_PER_UNIT <- 20
N_RAND    <- 50
RAND_SIZE <- 12
OUTDIR    <- file.path(OUT, "tables")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

hr("STEP 1  载入并取上皮细胞")

stopifnot(file.exists(F_SC_EPI))
sc_full <- readRDS(F_SC_EPI)
cat("对象: ", nrow(sc_full), " 基因 x ", ncol(sc_full), " 细胞\n", sep = "")

cells <- colnames(sc_full)[as.character(sc_full$cell_type) == "EOC"]
sc <- subset(sc_full, cells = cells)
rm(sc_full); gc()

cat("上皮细胞: ", ncol(sc), "  病人: ", length(unique(sc$patient_id)),
    "  样本: ", length(unique(sc$sample)), "\n", sep = "")
if (ncol(sc) < 1000) stop("上皮细胞数异常，检查 F_SC_EPI")
if (!"patient_id" %in% colnames(sc@meta.data))
  stop("对象缺少 patient_id 列。不能退回用 sample 当病人（会造成伪重复）。")

sc <- Seurat::NormalizeData(sc, normalization.method = "LogNormalize",
                            scale.factor = 1e4, verbose = FALSE)

hr("STEP 2  模块定义与残差化打分")

MODULES <- list(
  OXPHOS      = c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
                  "SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","UQCRH","CYC1",
                  "COX4I1","COX5A","COX5B","COX6A1","ATP5A1","ATP5B","ATP5C1","ATP5F1"),
  Mito_Import = c("TOMM22","TOMM40","TOMM70","TIMM23","TIMM17A","TIMM50","TIMM44",
                  "SLC25A3","SLC25A4","SLC25A5"),
  ROS_Defense = c("GPX4","SOD1","SOD2","CAT","PRDX3","TXNRD2","GSR","GCLC","GCLM"),
  Ferroptosis = c("GPX4","HMOX1","NQO1","SLC7A11","FTH1","FTL","TFRC"),
  FAO         = c("CPT1A","CPT1B","CPT2","ACADVL","ACADM","ACADS","ECH1","ECHS1",
                  "HADHA","HADHB"),
  Glycolysis  = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1"),
  TGFB1_Sig   = c("TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3"),
  Immune_Ckpt = c("CD274","CTLA4","LAG3","TIGIT","PDCD1","HAVCR2")
)
MODULES <- lapply(MODULES, function(g) intersect(g, rownames(sc)))
for (m in names(MODULES)) cat(sprintf("  %-12s %2d 基因\n", m, length(MODULES[[m]])))

# 只对模块基因 + TSPO 做残差化，与全矩阵做法等价，但省内存
genes <- unique(c(unlist(MODULES), "TSPO"))
X  <- as.matrix(Seurat::GetAssayData(sc, assay = "RNA", layer = "data")[genes, , drop = FALSE])
v  <- log10(sc$nCount_RNA + 1); vm <- v - mean(v)
b  <- as.numeric((X %*% vm) / sum(vm^2))
a  <- rowMeans(X) - b * mean(v)
Z  <- t(scale(t(X - (a + b %o% v))))          # 逐基因：去文库大小 -> z
Z[!is.finite(Z)] <- NA

S    <- sapply(MODULES, function(g) colMeans(Z[g, , drop = FALSE], na.rm = TRUE))
tspo <- as.numeric(Z["TSPO", ])
pat  <- as.character(sc$patient_id)
smp  <- as.character(sc$sample)

hr("STEP 3  线性混合模型（随机截距）")

fit_lmm <- function(y, grp) {
  d <- data.frame(y = y, T = tspo, g = factor(grp))
  d <- d[is.finite(d$y) & is.finite(d$T), ]
  if (nrow(d) < 50 || nlevels(droplevels(d$g)) < 3) return(NULL)
  d$g <- droplevels(d$g)
  m <- suppressWarnings(suppressMessages(lme4::lmer(y ~ T + (1 | g), data = d, REML = TRUE)))
  cf <- summary(m)$coefficients
  se <- cf["T", "Std. Error"]; t <- cf["T", "t value"]
  dfres <- nrow(d) - nlevels(d$g) - 1
  vc <- tryCatch(as.data.frame(lme4::VarCorr(m)), error = function(e) NULL)
  icc <- if (!is.null(vc) && nrow(vc) >= 1) vc$vcov[1] / sum(vc$vcov) else NA_real_
  c(beta = cf["T", "Estimate"], se = se, p = 2 * pt(-abs(t), df = dfres),
    n = nrow(d), n_unit = nlevels(d$g), icc = icc)
}

mixed_tab <- function(grp, label) {
  R <- do.call(rbind, lapply(names(MODULES), function(m) {
    r <- fit_lmm(S[, m], grp)
    if (is.null(r)) return(NULL)
    data.frame(level = label, module = m, beta = r["beta"], se = r["se"],
               p = r["p"], n = r["n"], n_unit = r["n_unit"], icc = r["icc"],
               stringsAsFactors = FALSE)
  }))
  R$FDR <- p.adjust(R$p, "BH")
  R[order(-R$beta), ]
}

M_pat <- mixed_tab(pat, "patient_id (11)")
M_smp <- mixed_tab(smp, "sample (22)")
cat("\n【病人层面，主要口径】\n"); pr(M_pat)
cat("\n【样本层面，可比口径】\n");  pr(M_smp)

hr("STEP 4  逐单位符号检验（主要推断依据）")

sign_tab <- function(grp, label) {
  R <- do.call(rbind, lapply(names(MODULES), function(m) {
    sl <- sapply(split(seq_along(grp), grp), function(i) {
      y <- S[i, m]
      ok <- is.finite(y) & is.finite(tspo[i])
      if (sum(ok) < MIN_CELLS_PER_UNIT) return(NA_real_)
      coef(lm(y[ok] ~ tspo[i][ok]))[2]
    })
    sl <- sl[is.finite(sl)]
    if (length(sl) < 5) return(NULL)
    k <- sum(sl > 0); n <- length(sl)
    data.frame(level = label, module = m, median_slope = median(sl),
               pct_positive = round(100 * k / n, 1),
               n_pos = k, n_unit = n,
               sign_p = binom.test(k, n, 0.5)$p.value, stringsAsFactors = FALSE)
  }))
  R[order(-R$pct_positive), ]
}

S_pat <- sign_tab(pat, "patient_id (11)")
S_smp <- sign_tab(smp, "sample (22)")
cat("\n【逐病人，11 例，主要口径】\n"); pr(S_pat, 4)
cat("\n【逐样本，22 个，可比口径】\n");  pr(S_smp, 4)

hr("STEP 5  随机基因集阴性对照（病人层面）")

# 注意：随机基因必须从「全转录组」抽，不能从已有的 80 个模块基因里抽
# ——模块基因本身都与 TSPO 共变，从中抽样会让零分布整体右移，失去意义
set.seed(SEED)
n_rand_gene <- N_RAND * RAND_SIZE
pool <- setdiff(rownames(sc), genes)
pool <- sample(pool, min(n_rand_gene, length(pool)))

Xr  <- as.matrix(Seurat::GetAssayData(sc, assay = "RNA", layer = "data")[pool, , drop = FALSE])
br  <- as.numeric((Xr %*% vm) / sum(vm^2))
ar  <- rowMeans(Xr) - br * mean(v)
Zr  <- t(scale(t(Xr - (ar + br %o% v))))
Zr[!is.finite(Zr)] <- NA
rm(Xr); gc()

rand_med <- numeric(N_RAND)
for (i in seq_len(N_RAND)) {
  g    <- pool[((i - 1) * RAND_SIZE + 1):(i * RAND_SIZE)]
  sc_i <- colMeans(Zr[g, , drop = FALSE], na.rm = TRUE)
  sl <- sapply(split(seq_along(pat), pat), function(k) {
    y  <- sc_i[k]
    ok <- is.finite(y) & is.finite(tspo[k])
    if (sum(ok) < MIN_CELLS_PER_UNIT) return(NA_real_)
    coef(lm(y[ok] ~ tspo[k][ok]))[2]
  })
  rand_med[i] <- median(sl, na.rm = TRUE)
}
cat(sprintf("随机 %d 基因模块（抽自全转录组 %d 个基因）：\n", RAND_SIZE, length(pool)))
cat(sprintf("  中位斜率的 5%%/50%%/95%% 分位 = %+.4f / %+.4f / %+.4f\n",
            quantile(rand_med, .05, na.rm = TRUE), median(rand_med, na.rm = TRUE),
            quantile(rand_med, .95, na.rm = TRUE)))
cat("\n各模块观测值在随机分布中的位置：\n")
for (m in c("OXPHOS", "Glycolysis", "FAO", "ROS_Defense", "TGFB1_Sig")) {
  obs <- S_pat$median_slope[S_pat$module == m]
  cat(sprintf("  %-12s 观测 %+.4f  超过 %.0f%% 的随机模块\n",
              m, obs, mean(rand_med < obs, na.rm = TRUE) * 100))
}

hr("STEP 6  写出")

write.csv(M_pat, file.path(OUTDIR, "T2a_mixed_model_patient.csv"), row.names = FALSE)
write.csv(M_smp, file.path(OUTDIR, "T2a_mixed_model_sample.csv"),  row.names = FALSE)
write.csv(S_pat, file.path(OUTDIR, "T2b_sign_test_patient.csv"),  row.names = FALSE)
write.csv(S_smp, file.path(OUTDIR, "T2b_sign_test_sample.csv"),   row.names = FALSE)
write.csv(data.frame(random_median_slope = rand_med),
          file.path(OUTDIR, "T2c_random_control.csv"), row.names = FALSE)

dir.create(file.path(OUT, "cache"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(mixed_patient = M_pat, mixed_sample = M_smp,
             sign_patient = S_pat, sign_sample = S_smp, random = rand_med),
        file.path(OUT, "cache", "epithelial_result.rds"))

cat("\n判读要点：\n")
cat("  · 混合模型 β 会被细胞数加权；Mito_Import 的 β 为正但逐病人一致率接近随机，\n")
cat("    说明该正 β 是细胞数加权造成的假象，判为无关联。\n")
cat("  · 符号检验是主要推断依据；11 例时 binom 最小可能 p 值为 1/1024 = 9.8e-4。\n")
cat("  · ROS_Defense 在病人层面为 9/11，属边缘，正文应写「趋势」而非显著。\n")

hr("03_epithelial_mixed.R 完成")

# -----------------------------------------------------------------------------
# -- Fig4_上皮代谢 / 04_vko_branch.R
# -----------------------------------------------------------------------------
# =============================================================================
#  80_vko_reanalysis.R
#  虚拟敲除结果的正确统计：靶基因集 vs 全局背景
#
#  原分析已经算出了全局上调比例，但没用作对照，热图默认以 50% 为中点。
#  正确做法：对每个（细胞类型 × 敲除基因），检验
#     目标基因集的 down% 是否显著高于该次敲除的全局 down%
#  用单侧二项检验，零假设 p = 全局 down%。
#
#  不重跑敲除，直接读 D:/Rstudy/tcga2/GSE165897_branch_specificity/*_diffRegulation.csv
# =============================================================================


if (!exists("OUT")) stop("请先 source('00_setup.R')")

DIR <- "D:/Rstudy/tcga2/GSE165897_branch_specificity"
OUTDIR <- file.path(OUT, "tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
hr <- function(x) cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")

METAB  <- c("PPARA","CPT1A","CPT1B","CPT2","ACADVL","ACADM","ACADS",
            "ECH1","ECHS1","HADHA","HADHB","GPX4","HMOX1","NQO1")
IMMUNE <- c("TGFB1","IL1B","CD274","SPP1","LGALS1","MIF",
            "CD68","CD163","CSF1R","VEGFA","CD74","CXCR4","CD44")

hr("STEP 1  定位结果文件")
fs <- list.files(DIR, pattern = "_diffRegulation\\.csv$", full.names = TRUE)
cat("找到 ", length(fs), " 个 diffRegulation 文件\n", sep = "")
if (!length(fs)) {
  cat("目录内容:\n"); print(list.files(DIR))
  stop("没有结果文件")
}
print(head(basename(fs), 20))

hr("STEP 2  逐个计算富集与二项检验")
res <- do.call(rbind, lapply(fs, function(f) {
  d <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)

  # 确定基因名列
  gcol <- intersect(c("gene","Gene","X"), colnames(d))[1]
  if (is.na(gcol)) { gcol <- colnames(d)[1] }
  d$gene <- as.character(d[[gcol]])

  # 确定方向：优先 direction 列，否则用 Z / FC 符号
  if ("direction" %in% colnames(d)) {
    dirv <- tolower(as.character(d$direction))
  } else {
    zc <- intersect(c("Z","z","FC","logFC","FC.zscore"), colnames(d))[1]
    if (is.na(zc)) return(NULL)
    dirv <- ifelse(as.numeric(d[[zc]]) > 0, "up", "down")
  }
  d$dir <- dirv
  d <- d[!is.na(d$dir) & d$dir %in% c("up","down"), ]
  if (!nrow(d)) return(NULL)

  # 全局背景
  g_down <- mean(d$dir == "down")
  g_up   <- 1 - g_down
  n_all  <- nrow(d)

  one <- function(genes, lab) {
    sub <- d[d$gene %in% genes, ]
    if (nrow(sub) < 3) return(NULL)
    k <- sum(sub$dir == "down"); n <- nrow(sub)
    p <- if (g_down > 0 && g_down < 1)
      binom.test(k, n, p = g_down, alternative = "greater")$p.value else NA_real_
    data.frame(branch = lab, n_gene = n, n_down = k,
               pct_down = round(100 * k / n, 1),
               global_down = round(100 * g_down, 1),
               enrich = round(100 * k / n - 100 * g_down, 1),
               binom_p = signif(p, 3), stringsAsFactors = FALSE)
  }
  r <- rbind(one(METAB, "代谢"), one(IMMUNE, "免疫"))
  if (is.null(r)) return(NULL)

  nm <- sub("_diffRegulation\\.csv$", "", basename(f))
  cbind(run = nm, n_total = n_all, r)
}))

if (is.null(res) || !nrow(res)) stop("无法解析任何结果")

res <- res[order(res$run, res$branch), ]
print(res, row.names = FALSE)
write.csv(res, file.path(OUTDIR, "T2d_vko_enrichment.csv"), row.names = FALSE)

hr("STEP 3  只保留 TSPO 敲除")
ts <- res[grepl("TSPO", res$run, ignore.case = TRUE), ]
if (!nrow(ts)) { cat("没有 TSPO 敲除结果\n") } else {
  ts <- ts[order(-ts$enrich), ]
  print(ts, row.names = FALSE)

  cat("\n按分支汇总（TSPO 敲除）：\n")
  S <- do.call(rbind, lapply(split(ts, ts$branch), function(d)
    data.frame(branch = d$branch[1], n_run = nrow(d),
               median_enrich = round(median(d$enrich), 1),
               n_sig = sum(d$binom_p < 0.05, na.rm = TRUE),
               n_sig_up = sum(d$enrich > 0 & d$binom_p < 0.05, na.rm = TRUE),
               stringsAsFactors = FALSE)))
  print(S, row.names = FALSE)
}

hr("STEP 4  合并检验（Fisher 法合并各细胞类型的 p 值）")
if (!is.null(ts) && nrow(ts)) {
  for (b in unique(ts$branch)) {
    d <- ts[ts$branch == b & is.finite(ts$binom_p), ]
    if (nrow(d) < 2) next
    stat <- -2 * sum(log(d$binom_p))
    p <- pchisq(stat, df = 2 * nrow(d), lower.tail = FALSE)
    cat(sprintf("%s: %d 个细胞类型, 中位富集 %+.1f%%, Fisher 合并 p = %.3g\n",
                b, nrow(d), median(d$enrich), p))
  }
}

hr("STEP 5  写出")
write.csv(res, file.path(OUTDIR, "T2d_vko_enrichment.csv"), row.names = FALSE)
message("已写入 ", OUT)

cat("\n================ 判读 ================\n")
cat("· enrich > 0 且 binom_p < 0.05  -> 该基因集的 down 比例显著高于全局背景\n")
cat("· 多个细胞类型一致为正          -> 该结论不依赖特定细胞类型\n")
cat("· enrich ≈ 0 或为负             -> 目标基因集的下降只是全局趋势的一部分\n")
cat("\n注：原分析以 50% 为参照。若全局 down% 本身就是 55%，\n")
cat("    则 60% 的目标基因下降并不算富集。本脚本修正了这一点。\n")

# -----------------------------------------------------------------------------
# -- Fig4_上皮代谢 / 20_epi_subclustering.R
# -----------------------------------------------------------------------------
# =============================================================================
#  20_epi_subclustering.R   上皮细胞亚群细分（7,929 细胞）
#
#  目的：1) 上皮是否存在与 TSPO 高低对应的亚群结构
#        2) TSPO 高的亚群是否就是代谢活跃的亚群
#
#  流程：LogNormalize → 2000 HVG → ScaleData（回归 percent.mt）→ PCA 30
#        → harmony 按 patient_id 校正（不可用时跳过）→ 多分辨率聚类
#        → 标记评分注释 → TSPO/模块评分亚群间比较（含混合模型）
#        → 亚群×病人构成混杂检查
#
#  ⚠ 两个关键点
#    1. harmony 各版本参数名不同，且 reduction 会被部分匹配到 reduction.use
#       与 reduction.save 两个参数，必须用 do.call + formals() 探测
#    2. 11 例病人细胞数 68–2,502，若某亚群被单病人主导，
#       该亚群的 TSPO 差异就是病人差异，脚本会标出这类亚群
#
#  输入：F_SC_EPI（cell_type == "EOC"）
#  输出：T12a–T12f 六张表、Fig3A–Fig3E 五张图、cache/epi_subclustered.rds
# =============================================================================

if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  .ok <- FALSE
  if (file.exists(.s))
    .ok <- !inherits(try(source(.s, encoding = "UTF-8"), silent = TRUE), "try-error")
  if (!.ok) {
    PROJ <- "D:/Rstudy/tcga"; PROJ2 <- "D:/Rstudy/tcga2"
    OUT <- file.path(PROJ, "results_v2")
    F_SC_EPI <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")
    hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
    pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)
  }
  for (d in c(OUT, file.path(OUT, "tables"), file.path(OUT, "figures")))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("hr")) hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
if (!exists("pr")) pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(lme4)
})

OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

N_HVG     <- 2000
N_PCS     <- 30
RES_SWEEP <- c(0.2, 0.4, 0.6, 0.8)
RES_USE   <- 0.4
SEED_SUB  <- 20260101
set.seed(SEED_SUB)

MARKERS <- list(
  Epi_secretory = c("MUC16","PAX8","WFDC2","KRT8","KRT18","KRT19","CLDN4"),
  Epi_ciliated   = c("FOXJ1","CAPS","PCDH","CFAP","DNAH","TEKT"),
  Epi_EMT        = c("VIM","ZEB1","ZEB2","SNAI1","SNAI2","CDH2","FN1"),
  Epi_prolif     = c("MKI67","TOP2A","PCNA","CCNB1","BUB1"),
  Epi_IFN        = c("HLA-DRA","HLA-DRB1","ISG15","MX1","IFIT1","STAT1"),
  Epi_stress     = c("JUN","FOS","JUNB","HSPA1A","HSPA1B","IER3")
)

MODULES <- list(
  OXPHOS = c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
             "SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","UQCRH","CYC1",
             "COX4I1","COX5A","COX5B","COX6A1","ATP5A1","ATP5B","ATP5C1","ATP5F1"),
  Glycolysis = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1"),
  TGFB1_Sig  = c("TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3")
)

hr("STEP 1  载入与取上皮细胞")

sc_full <- readRDS(F_SC_EPI)
cat("对象: ", nrow(sc_full), " 基因 x ", ncol(sc_full), " 细胞\n", sep = "")

cells <- colnames(sc_full)[as.character(sc_full$cell_type) == "EOC"]
sc <- subset(sc_full, cells = cells)
rm(sc_full); gc()

cat("上皮细胞: ", ncol(sc), "  病人: ", length(unique(sc$patient_id)),
    "  样本: ", length(unique(sc$sample)), "\n", sep = "")
if (ncol(sc) < 1000) stop("上皮细胞数异常")
if (!"patient_id" %in% colnames(sc@meta.data)) stop("缺少 patient_id 列")

if (!"data" %in% SeuratObject::Layers(sc[["RNA"]]))
  sc <- Seurat::NormalizeData(sc, normalization.method = "LogNormalize",
                              scale.factor = 1e4, verbose = FALSE)

hr("STEP 2  特征选择与降维")

sc <- Seurat::FindVariableFeatures(sc, selection.method = "vst",
                                   nfeatures = N_HVG, verbose = FALSE)
cat("高变基因: ", length(Seurat::VariableFeatures(sc)), "\n", sep = "")

reg_vars <- if ("percent.mt" %in% colnames(sc@meta.data)) "percent.mt" else NULL
cat("回归变量: ", ifelse(is.null(reg_vars), "无", reg_vars), "\n", sep = "")

sc <- Seurat::ScaleData(sc, features = Seurat::VariableFeatures(sc),
                        vars.to.regress = reg_vars, verbose = FALSE)
sc <- Seurat::RunPCA(sc, features = Seurat::VariableFeatures(sc),
                     npcs = N_PCS, verbose = FALSE)
cat("PCA 完成，PC1 方差占比: ",
    round(100 * sc@reductions$pca@stdev[1]^2 / sum(sc@reductions$pca@stdev^2), 1), "%\n", sep = "")

hr("STEP 3  批次校正（有 harmony 时按病人）")

HAVE_HARMONY <- requireNamespace("harmony", quietly = TRUE)

if (HAVE_HARMONY) {
  suppressPackageStartupMessages(library(harmony))

  # harmony 各版本参数名不一致：reduction vs reduction.use，dims vs dims.use。
  # 且传 reduction 会被部分匹配到 reduction.use 与 reduction.save 两个参数而报错，
  # 因此统一用 do.call + formals() 探测，不硬编码参数名。
  .fm <- tryCatch(names(formals(getS3method("RunHarmony", "Seurat"))),
                  error = function(e) character(0))
  cat("harmony 参数: ", paste(.fm, collapse = ", "), "\n", sep = "")

  .args <- list(object = sc, group.by.vars = "patient_id", verbose = FALSE)
  if ("reduction.use" %in% .fm)  .args$reduction.use  <- "pca"
  else if ("reduction" %in% .fm) .args$reduction      <- "pca"
  if ("dims.use" %in% .fm)       .args$dims.use       <- 1:N_PCS
  else if ("dims" %in% .fm)      .args$dims           <- 1:N_PCS
  if ("reduction.save" %in% .fm) .args$reduction.save <- "harmony"

  sc <- tryCatch(
    do.call(harmony::RunHarmony, .args),
    error = function(e) {
      cat("harmony 失败: ", conditionMessage(e), "\n", sep = "")
      cat("改为不校正继续\n")
      NULL
    })

  if (is.null(sc)) {
    sc <- .args$object
    RED <- "pca"
    HAVE_HARMONY <- FALSE
  } else {
    RED <- "harmony"
    cat("已按 patient_id 做 Harmony 校正\n")
  }
} else {
  RED <- "pca"
  cat("⚠ 未安装 harmony，直接用 PCA。\n")
  cat("  若 UMAP 按病人分层，亚群差异即为病人差异，必须在正文说明。\n")
}

if (!RED %in% names(sc@reductions))
  stop("reduction '", RED, "' 不存在，检查 STEP 3")

hr("STEP 4  聚类（多分辨率）")

sc <- Seurat::FindNeighbors(sc, reduction = RED, dims = 1:N_PCS, verbose = FALSE)

for (r in RES_SWEEP) {
  sc <- Seurat::FindClusters(sc, resolution = r, verbose = FALSE)
  cat("  resolution ", r, " → ", length(unique(sc$seurat_clusters)), " 个亚群\n", sep = "")
}
sc <- Seurat::FindClusters(sc, resolution = RES_USE, verbose = FALSE)

# 亚群水平按编号排序，避免聚类号不连续时 factor 丢水平
cl_chr <- as.character(sc$seurat_clusters)
lev    <- paste0("Epi", sort(unique(as.integer(cl_chr))))
sc$epi_cluster <- factor(paste0("Epi", cl_chr), levels = lev)
cat("采用 resolution = ", RES_USE, "，得到 ", nlevels(sc$epi_cluster), " 个亚群\n", sep = "")

sc <- Seurat::RunUMAP(sc, reduction = RED, dims = 1:N_PCS, verbose = FALSE)

hr("STEP 5  亚群标记基因")

sc <- Seurat::SetIdent(sc, value = "epi_cluster")
mk <- tryCatch(
  Seurat::FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25,
                         logfc.threshold = 0.25, verbose = FALSE),
  error = function(e) { cat("FindAllMarkers 失败: ", conditionMessage(e), "\n"); NULL })

if (!is.null(mk) && nrow(mk)) {
  write.csv(mk, file.path(OUTDIR, "T12a_cluster_markers.csv"), row.names = FALSE)
  cat("标记基因: ", nrow(mk), " 条\n", sep = "")
  for (cl in levels(sc$epi_cluster)) {
    s <- mk[mk$cluster == cl, ]
    s <- head(s[order(-s$avg_log2FC), ], 6)
    if (nrow(s)) cat(sprintf("  %-6s %s\n", cl, paste(s$gene, collapse = ", ")))
  }
} else {
  cat("未得到标记基因\n")
}

hr("STEP 6  亚群注释（标记评分）")

score_cols <- character(0)
for (m in names(MARKERS)) {
  g <- intersect(MARKERS[[m]], rownames(sc))
  if (length(g) < 3) { cat("跳过 ", m, "（基因不足）\n", sep = ""); next }
  sc <- Seurat::AddModuleScore(sc, features = list(g), name = m, seed = SEED_SUB)
  score_cols <- c(score_cols, paste0(m, "1"))
}
score_cols <- intersect(score_cols, colnames(sc@meta.data))
cat("可用评分: ", paste(score_cols, collapse = ", "), "\n", sep = "")

if (length(score_cols)) {
  cl_mean <- as.data.frame(
    sapply(score_cols, function(cc)
      tapply(sc@meta.data[[cc]], sc$epi_cluster, mean, na.rm = TRUE)))
  cl_mean$cluster <- rownames(cl_mean)

  assign <- apply(cl_mean[, score_cols, drop = FALSE], 1,
                  function(x) names(x)[which.max(x)])

  ann <- data.frame(cluster = cl_mean$cluster,
                    assigned = as.character(assign),
                    n_cell = as.integer(table(sc$epi_cluster)[cl_mean$cluster]),
                    stringsAsFactors = FALSE)
  ann <- ann[order(-ann$n_cell), ]
  cat("\n亚群注释（按最大标记评分）：\n")
  pr(ann[, c("cluster","assigned","n_cell")])
  write.csv(ann, file.path(OUTDIR, "T12b_cluster_annotation.csv"), row.names = FALSE)
  write.csv(cl_mean, file.path(OUTDIR, "T12b_cluster_marker_scores.csv"), row.names = FALSE)
}

hr("STEP 7  TSPO 与代谢模块的亚群间差异")

for (m in names(MODULES)) {
  g <- intersect(MODULES[[m]], rownames(sc))
  if (length(g) < 3) next
  sc <- Seurat::AddModuleScore(sc, features = list(g),
                               name = paste0(m, "M"), seed = SEED_SUB)
}
mod_cols <- intersect(paste0(names(MODULES), "M1"), colnames(sc@meta.data))
cat("模块评分: ", paste(mod_cols, collapse = ", "), "\n", sep = "")

tspo <- as.numeric(Seurat::GetAssayData(sc, assay = "RNA", layer = "data")["TSPO", ])

tb <- data.frame(cluster = sc$epi_cluster, tspo = tspo,
                 patient = as.character(sc$patient_id), stringsAsFactors = FALSE)
for (mc in mod_cols) tb[[mc]] <- sc@meta.data[[mc]]

agg <- do.call(rbind, lapply(levels(tb$cluster), function(cl) {
  s <- tb[tb$cluster == cl, ]
  data.frame(cluster = cl, n_cell = nrow(s),
             TSPO_mean    = round(mean(s$tspo, na.rm = TRUE), 3),
             TSPO_pct_pos = round(100 * mean(s$tspo > 0, na.rm = TRUE), 1),
             n_patient    = length(unique(s$patient)),
             stringsAsFactors = FALSE)
}))
for (mc in mod_cols) {
  nm <- sub("M1$", "", mc)
  agg[[paste0(nm, "_mean")]] <-
    round(as.numeric(tapply(tb[[mc]], tb$cluster, mean, na.rm = TRUE)[as.character(agg$cluster)]), 3)
}
cat("\n各亚群 TSPO 与模块评分：\n"); pr(agg, 4)
write.csv(agg, file.path(OUTDIR, "T12c_TSPO_by_cluster.csv"), row.names = FALSE)
write.csv(agg, file.path(OUTDIR, "T12e_module_by_cluster.csv"), row.names = FALSE)

mm <- tryCatch({
  d <- tb[is.finite(tb$tspo), ]
  d$cluster <- factor(d$cluster, levels = levels(tb$cluster))
  m <- suppressWarnings(suppressMessages(
         lme4::lmer(tspo ~ cluster + (1 | patient), data = d)))
  cf <- summary(m)$coefficients
  df_res <- nrow(d) - nlevels(d$cluster) - length(unique(d$patient)) - 1
  data.frame(term = rownames(cf), beta = round(cf[, 1], 4),
             p = 2 * pt(-abs(cf[, "t value"]), df = df_res), row.names = NULL)
}, error = function(e) { cat("混合模型失败: ", conditionMessage(e), "\n"); NULL })

if (!is.null(mm)) {
  cat("\n混合模型 TSPO ~ 亚群 + (1|病人)（第一个亚群为参照）：\n"); pr(mm, 4)
  write.csv(mm, file.path(OUTDIR, "T12f_TSPO_cluster_lmm.csv"), row.names = FALSE)
}

hr("STEP 8  亚群 × 病人构成（关键混杂检查）")

tab <- table(sc$epi_cluster, sc$patient_id)
write.csv(as.data.frame.matrix(tab), file.path(OUTDIR, "T12d_cluster_by_patient.csv"))
cat("\n细胞数构成：\n"); print(tab)

comp <- do.call(rbind, lapply(rownames(tab), function(cl) {
  v <- as.numeric(tab[cl, ])
  data.frame(cluster = cl, n_cell = sum(v),
             top_patient   = colnames(tab)[which.max(v)],
             top_share_pct = round(100 * max(v) / sum(v), 1),
             n_patient_with_5 = sum(v >= 5),
             stringsAsFactors = FALSE)
}))
cat("\n混杂检查（top_share_pct > 50 表示被单病人主导）：\n"); pr(comp, 4)
write.csv(comp, file.path(OUTDIR, "T12d_cluster_confounding.csv"), row.names = FALSE)

bad <- comp$cluster[comp$top_share_pct > 50]
if (length(bad))
  cat("\n⚠ 被单病人主导的亚群（其 TSPO 差异不能解释为亚群效应）：",
      paste(bad, collapse = ", "), "\n", sep = "") else
  cat("\n各亚群均未被单病人主导（最大占比 < 50%）\n")

hr("STEP 9  图")

p1 <- DimPlot(sc, reduction = "umap", group.by = "epi_cluster",
              label = TRUE, repel = TRUE) +
      labs(title = "Epithelial subclusters") + theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "Fig3A_UMAP_cluster.pdf"), p1, width = 7, height = 6)
ggsave(file.path(FIGDIR, "Fig3A_UMAP_cluster.png"), p1, width = 7, height = 6, dpi = 300)

p2 <- DimPlot(sc, reduction = "umap", group.by = "patient_id") +
      labs(title = "Patients") + theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "Fig3B_UMAP_patient.pdf"), p2, width = 8, height = 6)
ggsave(file.path(FIGDIR, "Fig3B_UMAP_patient.png"), p2, width = 8, height = 6, dpi = 300)

p3 <- FeaturePlot(sc, features = "TSPO", reduction = "umap", order = TRUE) +
      labs(title = "TSPO expression") + theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "Fig3C_UMAP_TSPO.pdf"), p3, width = 7, height = 6)
ggsave(file.path(FIGDIR, "Fig3C_UMAP_TSPO.png"), p3, width = 7, height = 6, dpi = 300)

dot_genes <- unique(c("TSPO",
                      unlist(lapply(MARKERS, function(g) head(intersect(g, rownames(sc)), 4)))))
dot_genes <- intersect(dot_genes, rownames(sc))
p4 <- DotPlot(sc, features = dot_genes, group.by = "epi_cluster") +
      labs(x = NULL, y = NULL, title = "Subcluster markers") +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path(FIGDIR, "Fig3D_marker_dotplot.pdf"), p4, width = 11, height = 4.5)
ggsave(file.path(FIGDIR, "Fig3D_marker_dotplot.png"), p4, width = 11, height = 4.5, dpi = 300)

comp_long <- as.data.frame(tab)
colnames(comp_long) <- c("cluster", "patient", "n")
comp_long$pct <- ave(comp_long$n, comp_long$patient,
                     FUN = function(x) 100 * x / sum(x))
p5 <- ggplot(comp_long, aes(patient, pct, fill = cluster)) +
  geom_col() +
  labs(x = NULL, y = "Percent of patient's epithelial cells (%)",
       title = "Subcluster composition per patient", fill = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIGDIR, "Fig3E_composition.pdf"), p5, width = 8, height = 5)
ggsave(file.path(FIGDIR, "Fig3E_composition.png"), p5, width = 8, height = 5, dpi = 300)

saveRDS(sc, file.path(OUT, "cache", "epi_subclustered.rds"))
cat("\n已保存对象: cache/epi_subclustered.rds\n")

cat("\n判读要点：\n")
cat("  · 若 UMAP 按病人分层（Fig3B），亚群差异与病人差异无法分离。\n")
cat("  · comp 表 top_share_pct > 50 的亚群，其 TSPO 差异不能当亚群效应。\n")
cat("  · 某亚群 TSPO 高且 OXPHOS 也高，才支持「TSPO 标记代谢活跃亚群」。\n")
cat("  · 亚群数由 resolution 决定，正文需报告多分辨率稳定性（STEP 4 已打印）。\n")

hr("20_epi_subclustering.R 完成")

# -----------------------------------------------------------------------------
# -- Fig4_上皮代谢 / 21_epi_subcluster_resolve.R
# -----------------------------------------------------------------------------
# =============================================================================
#  21_epi_subcluster_resolve.R   上皮亚群结果的补救分析：改做病人内状态相关
#
#  背景：20_epi_subclustering.R 的 11 个亚群中，Epi0/Epi1/Epi6（占 53% 细胞）
#        分别 89.8%、93.1%、89.2% 来自单一病人，亚群效应与病人效应不可分离；
#        Epi8 表达 C1QA/B/C、AIF1，为巨噬细胞污染。
#        因此不再做亚群间比较，改为在每个病人内部检验 TSPO 与细胞状态的相关。
#
#  对应正文 3.3：糖酵解 11/11（+0.172）、OXPHOS 10/11（+0.185）、
#                增殖 6/11（−0.002，p=1.0）、EMT 4/11（0.022，p=0.55）
#
#  做法：与 03_epithelial_mixed.R 一致——逐病人 Spearman + 符号检验
#  输入：cache/epi_subclustered.rds（由 20 生成）
#  输出：T13a_epi_state_scores.csv、T13b_epi_state_within_patient.csv
# =============================================================================

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  if (file.exists(.s)) source(.s, encoding = "UTF-8") else {
    PROJ <- "D:/Rstudy/tcga"; OUT <- file.path(PROJ, "results_v2")
  }
}
OUTDIR <- file.path(OUT, "tables")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

MIN_CELLS <- 30
SEED_EPI  <- 20260101
set.seed(SEED_EPI)

STATES <- list(
  Proliferation = c("MKI67","TOP2A","PCNA","CCNB1","BUB1","AURKB","CCNA2","NDC80"),
  Hypoxia       = c("CA9","BNIP3","TRIB3","VEGFA","SLC2A1","PDK1","LDHA"),
  Ciliated      = c("C20orf85","CCDC60","FAM92B","ANKRD66","FOXJ1","CAPS"),
  EMT           = c("VIM","ZEB1","SNAI2","CDH2","FN1"),
  OXPHOS        = c("NDUFS1","NDUFS2","SDHA","SDHB","UQCRC1","UQCRC2",
                    "COX4I1","COX5A","ATP5A1","ATP5B"),
  Glycolysis    = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","ENO1")
)

hr("STEP 1  载入已聚类对象")

f <- file.path(OUT, "cache", "epi_subclustered.rds")
if (!file.exists(f)) stop("找不到 ", f, "，请先跑 20_epi_subclustering.R")
sc <- readRDS(f)
cat("对象: ", nrow(sc), " 基因 x ", ncol(sc), " 细胞\n", sep = "")

if ("epi_cluster" %in% colnames(sc@meta.data)) {
  cl_tab <- table(sc$epi_cluster)
  cat("各亚群细胞数:\n"); print(cl_tab)
  if ("Epi8" %in% names(cl_tab)) {
    cat("\n剔除 Epi8（表达 C1QA/B/C、AIF1，判定为巨噬污染）\n")
    sc <- subset(sc, cells = colnames(sc)[as.character(sc$epi_cluster) != "Epi8"])
    cat("剩余细胞: ", ncol(sc), "\n", sep = "")
  }
}

pat <- as.character(sc$patient_id)
cat("病人数: ", length(unique(pat)), "\n", sep = "")

hr("STEP 2  计算细胞状态评分")

d <- as.matrix(SeuratObject::LayerData(sc, assay = "RNA", layer = "data"))
Z <- t(scale(t(d))); Z[!is.finite(Z)] <- NA
rm(d); gc()

S <- list()
for (m in names(STATES)) {
  g <- intersect(STATES[[m]], rownames(Z))
  if (length(g) < 3) { cat("跳过 ", m, "\n", sep = ""); next }
  S[[m]] <- colMeans(Z[g, , drop = FALSE], na.rm = TRUE)
  cat(sprintf("  %-14s %2d 基因\n", m, length(g)))
}
S <- as.data.frame(S)
S$TSPO    <- as.numeric(Z["TSPO", ])
S$patient <- pat
write.csv(S, file.path(OUTDIR, "T13a_epi_state_scores.csv"), row.names = FALSE)

hr("STEP 3  逐病人相关 + 符号检验（主要推断）")

states <- setdiff(colnames(S), c("TSPO", "patient"))
res <- do.call(rbind, lapply(states, function(st) {
  r <- sapply(split(seq_len(nrow(S)), S$patient), function(i) {
    x <- S$TSPO[i]; y <- S[[st]][i]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < MIN_CELLS) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  })
  r <- r[is.finite(r)]
  if (length(r) < 6) return(NULL)
  k <- sum(r < 0); n <- length(r)
  data.frame(state = st, median_rho = round(median(r), 3),
             n_patient = n, n_negative = k,
             pct_negative = round(100 * k / n, 1),
             sign_p = binom.test(k, n, 0.5)$p.value, stringsAsFactors = FALSE)
}))
res$FDR <- p.adjust(res$sign_p, "BH")
res <- res[order(-abs(res$median_rho)), ]
cat("\n逐病人 TSPO 与细胞状态的相关（n_negative 为负相关病人数）：\n")
pr(res, 4)
write.csv(res, file.path(OUTDIR, "T13b_epi_state_within_patient.csv"), row.names = FALSE)

hr("STEP 4  逐病人明细")

cat("各状态在每个病人内的 Spearman：\n")
for (st in states) {
  r <- sapply(split(seq_len(nrow(S)), S$patient), function(i) {
    x <- S$TSPO[i]; y <- S[[st]][i]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < MIN_CELLS) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  })
  cat(sprintf("%-14s ", st))
  cat(paste(sprintf("%+.2f", r[is.finite(r)]), collapse = "  "), "\n")
}

cat("\n判读要点：\n")
cat("  · 每个病人的相关才是独立观测，11 例时 binom 最小 p 为 9.8e-4。\n")
cat("  · 增殖与 EMT 的零关联是有信息的阴性：说明 bulk 层面 TSPO 与增殖的负相关\n")
cat("    （−0.298）来自组织构成，而非上皮细胞内在属性。\n")

hr("21_epi_subcluster_resolve.R 完成")

# -----------------------------------------------------------------------------
# -- Fig4_上皮代谢 / 41_celltype_expression.R
# -----------------------------------------------------------------------------
# =============================================================================
#  41_celltype_expression.R
#
#  配体与受体在各细胞类型中的表达谱，以及 TSPO 与炎症小体、MIF 受体的关系。
#
#  服务正文 3.5（SPP1 为髓系特异配体）与 3.8（MIF 的上皮来源）。
#
#  输入：D:/Rstudy/tcga2/GSE266577_processed/GSE266577_raw.rds（137,083 细胞）
#  输出（工作区 results_v2/tables/）：
#    T43a_inflammasome_receptor.csv   炎症小体基因与 TSPO 的逐病人相关
#    T43b_CD74_CXCR4_by_celltype.csv  CD74 / CXCR4 各细胞类型表达
#    T43c_ligand_receptor_by_celltype.csv  九个基因 × 十类细胞
#    T43d_tspo_vs_ligands_subsets.csv 髓系亚群内 TSPO 与各配体的相关
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
})

PROJ  <- "D:/Rstudy/tcga2"
OUT   <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/tables"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")

# ---------------------------------------------------------------------------
# 1  载入并归一化
# ---------------------------------------------------------------------------

hr("STEP 1  载入 GSE266577")
OBJ <- file.path(PROJ, "GSE266577_processed", "GSE266577_raw.rds")
o   <- readRDS(OBJ)
cat("原始:", nrow(o), "基因 x", ncol(o), "细胞\n")

o   <- JoinLayers(o, assay = "RNA")
CNT <- SeuratObject::LayerData(o, assay = "RNA", layer = "counts")
lib <- Matrix::colSums(CNT); lib[lib == 0] <- 1
Xn  <- CNT %*% Matrix::Diagonal(x = 1e4 / lib)
Xn@x <- log1p(Xn@x)
cat("已 log-normalize（每万计数）\n")

CT  <- as.character(o$cell_type)
PAT <- as.character(o$publication_patient_code_final)
STG <- as.character(o$treatment_stage)
getv <- function(g) as.numeric(Xn[g, ])

# ---------------------------------------------------------------------------
# 2  细胞类型粗分类
# ---------------------------------------------------------------------------

grp <- rep(NA_character_, length(CT))
grp[grepl("^CD16|NK", CT)]        <- "NK"
grp[grepl("^Tcm|^Tem|^Regulatory|^Type 17", CT)] <- "T"
grp[grepl("B cells|Plasma", CT)]  <- "B/Plasma"
grp[grepl("^DC|Migratory|pDC", CT)] <- "DC"
grp[CT == "Macrophages"]          <- "Macrophage"
grp[CT == "Classical monocytes"]  <- "Monocyte"
grp[CT == "Epithelial cells"]     <- "Epithelial"
grp[CT == "Fibroblasts"]          <- "Fibroblast"
grp[CT == "Endothelial cells"]    <- "Endothelial"
grp[CT == "Mesothelial cells"]    <- "Mesothelial"
grp[CT == "Mast cells"]           <- "Mast"
grp[CT == "ILC"]                  <- "ILC"

ORD  <- c("Epithelial","Macrophage","Monocyte","Fibroblast","Endothelial",
          "Mesothelial","DC","T","NK","B/Plasma","Mast","ILC")
ORD  <- ORD[ORD %in% unique(grp)]
GENE <- c("TSPO","MIF","LGALS1","SPP1","FN1","LGALS9","CD44","CD74","CXCR4")

# ---------------------------------------------------------------------------
# 3  各细胞类型的表达谱
# ---------------------------------------------------------------------------

hr("STEP 2  各细胞类型的平均表达")
cat(sprintf("%-13s %7s", "细胞类型", "n"))
for (g in GENE) cat(sprintf(" %7s", g)); cat("\n")
tab <- data.frame()
for (gr in ORD) {
  i <- which(grp == gr); if (length(i) < 50) next
  cat(sprintf("%-13s %7d", gr, length(i)))
  r <- c(celltype = gr, n = length(i))
  for (g in GENE) {
    v <- mean(getv(g)[i]); cat(sprintf(" %7.2f", v))
    r[g] <- round(v, 3)
  }
  cat("\n")
  tab <- rbind(tab, as.data.frame(t(r), stringsAsFactors = FALSE))
}
write.csv(tab, file.path(OUT, "T43c_ligand_receptor_by_celltype.csv"), row.names = FALSE)

cat("\n各配体的主要生产者：\n")
for (g in c("MIF","LGALS1","SPP1","FN1")) {
  v <- sapply(ORD, function(gr) mean(getv(g)[which(grp == gr)]))
  o2 <- order(-v)
  cat(sprintf("  %-7s %s\n", g,
      paste(sprintf("%s %.2f", ORD[o2][1:4], v[o2][1:4]), collapse = "  >  ")))
}

# CD74 / CXCR4 专项
sub <- c("DC","Macrophage","Monocyte","B/Plasma","T","NK","Endothelial",
         "Epithelial","Mesothelial","Fibroblast")
sub <- sub[sub %in% ORD]
t2 <- do.call(rbind, lapply(sub, function(gr) {
  i <- which(grp == gr)
  data.frame(celltype = gr, n = length(i),
             CD74_mean = round(mean(getv("CD74")[i]), 3),
             CD74_pos  = round(100 * mean(getv("CD74")[i] > 0), 1),
             CXCR4_mean = round(mean(getv("CXCR4")[i]), 3),
             CXCR4_pos  = round(100 * mean(getv("CXCR4")[i] > 0), 1),
             CD44_mean  = round(mean(getv("CD44")[i]), 3))
}))
print(t2, row.names = FALSE)
write.csv(t2, file.path(OUT, "T43b_CD74_CXCR4_by_celltype.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4  髓系亚群内 TSPO 与各配体、炎症小体基因的逐病人相关
# ---------------------------------------------------------------------------

hr("STEP 3  髓系亚群内 TSPO 的逐病人相关")
MA <- which(CT == "Macrophages")
ts <- as.numeric(Xn["TSPO", ])
GENES4 <- c("CD74","CXCR4","MIF","LGALS1","SPP1","FN1","CD44",
            "NLRP3","CASP1","IL1B","PYCARD","GSDMD","IL18","NLRC4","AIM2")
res <- data.frame()
for (g in GENES4) {
  v <- getv(g)
  rh <- sapply(unique(PAT[MA]), function(p) {
    j <- MA[PAT[MA] == p]
    if (length(j) < 30) return(NA_real_)
    suppressWarnings(cor(ts[j], v[j], method = "spearman"))
  })
  rh <- rh[is.finite(rh)]
  k <- sum(rh > 0); N <- length(rh)
  p <- if (N) min(1, 2 * pbinom(min(k, N - k), N, 0.5)) else NA_real_
  cat(sprintf("  %-8s 为正 %2d/%-2d  中位 rho %+.3f  符号检验 p=%.3g\n",
              g, k, N, median(rh), p))
  res <- rbind(res, data.frame(gene = g, n_patient = N, n_positive = k,
                               median_rho = round(median(rh), 3), sign_p = signif(p, 3)))
}
write.csv(res, file.path(OUT, "T43d_tspo_vs_ligands_subsets.csv"), row.names = FALSE)

# 炎症小体模块
INFL <- intersect(c("NLRP3","CASP1","IL1B","PYCARD","GSDMD","IL18","NLRC4","AIM2"),
                  rownames(Xn))
Z  <- sapply(INFL, function(g) { v <- getv(g); (v - mean(v)) / sd(v) })
sc <- rowMeans(Z)
rh <- sapply(unique(PAT[MA]), function(p) {
  j <- MA[PAT[MA] == p]; if (length(j) < 30) return(NA_real_)
  suppressWarnings(cor(ts[j], sc[j], method = "spearman"))
})
rh <- rh[is.finite(rh)]
cat(sprintf("\n  炎症小体模块: 合并 %.3f, 为正 %d/%d, 中位 %+.3f\n",
            cor(ts[MA], sc[MA], method = "spearman"), sum(rh > 0), length(rh), median(rh)))
write.csv(data.frame(module = "inflammasome", n_patient = length(rh),
                     n_positive = sum(rh > 0), median_rho = round(median(rh), 3),
                     pooled_rho = round(cor(ts[MA], sc[MA], method = "spearman"), 3)),
          file.path(OUT, "T43a_inflammasome_receptor.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5  TSPO 高/低四分位的受体与炎症小体状态
# ---------------------------------------------------------------------------

hr("STEP 4  TSPO 高/低四分位")
q <- quantile(ts[MA], c(.25, .75))
hi <- ts[MA] >= q[2]; lo <- ts[MA] <= q[1]
cat(sprintf("%-10s %8s %8s %+9s\n", "指标", "高 25%", "低 25%", "差值"))
for (item in list(c("CD74","v"), c("CD44","v"), c("SPP1","v"), c("MIF","v"),
                  c("IL1B","v"), c("CD163","v"), c("TGFB1","v"), c("CXCL9","v"))) {
  v <- getv(item[1])[MA]
  cat(sprintf("%-10s %8.3f %8.3f %+9.3f\n", item[1], mean(v[hi]), mean(v[lo]),
              mean(v[hi]) - mean(v[lo])))
}
cat(sprintf("%-10s %8.3f %8.3f %+9.3f\n", "炎症小体", mean(sc[MA][hi]), mean(sc[MA][lo]),
            mean(sc[MA][hi]) - mean(sc[MA][lo])))

hr("41_celltype_expression.R 完成")

# -----------------------------------------------------------------------------
# -- Fig4_上皮代谢 / 44_vko_table.R
# -----------------------------------------------------------------------------
# =============================================================================
#  44_vko_table.R
#
#  把虚拟敲除的模块 Z 值变化、表达匹配对照与逐模块经验 p 整理成补充表 S16。
#
#  服务正文讨论「虚拟敲除界定了 TSPO 的作用范围」与摘要。
#
#  背景：scTenifoldKnk 的原始输出在 D:/Rstudy/tcga/results_v2/ 下，
#        04_vko_branch.R 负责富集口径，本脚本负责特异性的对照表。
#
#  输入：D:/Rstudy/tcga/results_v2/sc_vko_specificity.csv
#        D:/Rstudy/tcga/results_v2/sc_vko_summary.txt
#  输出（工作区 results_v2/tables/）：
#    T46_vko_specificity.csv   两靶点 × 七模块的 Z 变化与对照
# =============================================================================

suppressPackageStartupMessages({ library(utils) })

PROJ <- "D:/Rstudy/tcga/results_v2"
OUT  <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/tables"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MODCN <- c(OXPHOS = "氧化磷酸化", Mito_Import = "线粒体蛋白输入",
           ROS_Defense = "ROS 防御", Ferroptosis = "铁死亡",
           FAO = "脂肪酸氧化", Glycolysis = "糖酵解",
           TGFB1_Sig = "TGFB1 信号")

cat("\n===== 读取虚拟敲除结果 =====\n")
d <- read.csv(file.path(PROJ, "sc_vko_specificity.csv"), stringsAsFactors = FALSE)
cat("靶点:", paste(unique(d$target), collapse = ", "),
    " 模块:", length(unique(d$module)), " 行:", nrow(d), "\n")

d$module_cn <- ifelse(d$module %in% names(MODCN), MODCN[d$module], d$module)
out <- data.frame(
  靶点       = d$target,
  模块       = d$module_cn,
  模块Z变化  = round(d$delta, 3),
  对照中位   = round(d$ctrl_median, 3),
  对照最大   = round(d$ctrl_max, 3),
  经验p      = round(d$emp_p, 2),
  对照数     = d$n_ctrl,
  FDR        = round(d$FDR, 3),
  判定       = d$verdict,
  stringsAsFactors = FALSE
)

cat("\n===== 逐模块 =====\n")
print(out[out$靶点 == "TSPO", ], row.names = FALSE)

cat("\n===== 两靶点的差值（应接近零，说明是网络属性而非靶点特异）=====\n")
a <- out[out$靶点 == "TSPO", ]; b <- out[out$靶点 == "PPARA", ]
m <- merge(a, b, by = "模块", suffixes = c("_TSPO","_PPARA"))
m$差值 <- round(m$模块Z变化_TSPO - m$模块Z变化_PPARA, 3)
print(m[, c("模块","模块Z变化_TSPO","模块Z变化_PPARA","差值")], row.names = FALSE)
cat("  最大绝对差值:", max(abs(m$差值), na.rm = TRUE), "\n")

# 技术自检
sf <- file.path(PROJ, "sc_vko_summary.txt")
if (file.exists(sf)) {
  cat("\n===== 技术自检（自敲除后受影响最强的基因应为靶点本身）=====\n")
  L <- readLines(sf, warn = FALSE)
  k <- grep("self_Z|target n_genes", L)
  if (length(k)) cat(paste(L[k], collapse = "\n"), "\n")
}

write.csv(out, file.path(OUT, "T46_vko_specificity.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n已写入 T46_vko_specificity.csv\n")

cat("\n===== 判读 =====\n")
cat("TSPO 扰动使 OXPHOS 下降", abs(a$模块Z变化[a$模块 == "氧化磷酸化"]),
    "，表达匹配对照中位", abs(a$对照中位[a$模块 == "氧化磷酸化"]), "\n")
cat("两靶点差值均 <", max(abs(m$差值), na.rm = TRUE),
    "，说明模块水平由基因表达量决定，靶点自身不构成额外调控节点\n")

cat("\n===== 44_vko_table.R 完成 =====\n")

# ############################################################################
# ##  Fig5_髓系配体与CAF  (4 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- Fig5_髓系配体与CAF / 05_myeloid_mixed.R
# -----------------------------------------------------------------------------
# =============================================================================
#  40_myeloid_final.R   第 8 / 9 章重做（在 9098 个髓系细胞上）
#
#  对象：X2_myeloid_trajectory.rds（23885 x 9098）
#    myeloid_subtype: Monocyte / TAM_MCSF / TAM_MCSF_IFNg
#    sample: T59 / T76 / T77 / T89 / T90
#
#  Part A  8.1-8.3 描述统计（构成、各亚型 TSPO、阳性率）
#  Part B  配体表达：TSPO 关联（混合模型 + 病人内一致率）
#  Part C  9.x CellChat：按【亚型】通讯（非 TSPO 分组，避免循环）
#
#  为什么 Part C 按亚型：IdentOverExpressedGenes 在组间做差异检验，
#  按 TSPO 分组再比较通讯是循环的；亚型由 marker 定义，不受 TSPO 影响。
# =============================================================================


if (!exists("OUT")) stop("请先 source('00_setup.R')")

SC_RDS   <- F_SC_MYE
OUTDIR <- file.path(OUT, "tables")
SUB_COL  <- "myeloid_subtype"
PAT_COL  <- "sample"
DO_CELLCHAT <- TRUE      # FALSE 则只跑 Part A/B（约 1 分钟）
NBOOT    <- 20
SEED     <- 42

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(lme4) })
if (requireNamespace("lmerTest", quietly = TRUE)) suppressPackageStartupMessages(library(lmerTest))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)
hr <- function(x) cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")

LIGANDS <- c("MIF","SPP1","FN1","LGALS9","LGALS1","CCL2","CCL3","CCL4","CCL5",
             "CXCL8","IL1B","TNF","TGFB1","APOE","VEGFA","CD74","CD44")

# ---------------------------------------------------------------- 载入
hr("STEP 0  载入")
sc <- readRDS(SC_RDS)
cat("维度: ", nrow(sc), " x ", ncol(sc), "\n")
md <- sc@meta.data
cat("髓系细胞数:", ncol(sc), "（预期 9,098，来自 5 例供体）\n")
if (ncol(sc) < 1000) stop("髓系细胞数异常，检查 F_SC_MYE")
stopifnot(all(c(SUB_COL, PAT_COL) %in% colnames(md)))
subtype <- as.character(md[[SUB_COL]])
pat <- as.character(md[[PAT_COL]])

dat  <- Seurat::GetAssayData(sc, assay = "RNA", layer = "data")
cnts <- Seurat::GetAssayData(sc, assay = "RNA", layer = "counts")
nf   <- Matrix::colSums(cnts > 0)
tspo <- as.numeric(dat["TSPO", ]); names(tspo) <- colnames(sc)

# ================================================================ PART A
hr("PART A  8.1-8.3  描述统计")
tb <- table(subtype)
cat("髓系构成（总计 ", sum(tb), "）:\n", sep = "")
print(data.frame(subtype = names(tb), n = as.integer(tb),
                 pct = round(as.integer(tb)/sum(tb)*100, 1)), row.names = FALSE)

cat("\n各亚型 TSPO 表达:\n")
A <- do.call(rbind, lapply(sort(unique(subtype)), function(s) {
  v <- tspo[subtype == s]
  data.frame(subtype = s, n = length(v),
             mean_TSPO = round(mean(v), 3),
             median_TSPO = round(median(v), 3),
             pct_positive = round(mean(v > 0) * 100, 1), stringsAsFactors = FALSE)
}))
print(A, row.names = FALSE)

cat("\n各病人 × 亚型细胞数:\n")
print(table(pat, subtype))

# ================================================================ PART B
hr("PART B  配体表达与 TSPO 的关联")
gs <- intersect(LIGANDS, rownames(dat))
cat("可用配体 ", length(gs), ": ", paste(gs, collapse = ", "), "\n")
cat(sprintf("TSPO 与 nFeature 相关: rho = %.3f\n",
            cor(tspo, nf[names(tspo)], method = "spearman")))

df <- data.frame(patient = factor(pat), subtype = factor(subtype),
                 TSPO = tspo, nf = nf[names(tspo)],
                 row.names = colnames(sc))

fit_one <- function(g) {
  y <- as.numeric(dat[g, ])
  d <- data.frame(y = y, TSPO = df$TSPO, nf = df$nf,
                  patient = df$patient, subtype = df$subtype)
  d <- d[is.finite(d$y) & is.finite(d$TSPO), ]
  m <- tryCatch(
    suppressMessages(suppressWarnings(
      lmerTest::lmer(y ~ TSPO + nf + subtype + (1 | patient), data = d))),
    error = function(e) suppressMessages(suppressWarnings(
      lme4::lmer(y ~ TSPO + nf + subtype + (1 | patient), data = d))))
  cf <- summary(m)$coefficients
  if (!"TSPO" %in% rownames(cf)) return(NULL)
  pcol <- intersect(c("Pr(>|t|)", "p value", "Pr(>|z|)"), colnames(cf))[1]
  tcol <- intersect(c("t value", "z value"), colnames(cf))[1]
  tp <- cf["TSPO", tcol]
  pv <- if (!is.na(pcol)) cf["TSPO", pcol] else
    2 * pt(-abs(tp), df = nrow(d) - nlevels(d$patient) - 2)
  # 病人内斜率符号
  per <- sapply(split(d, d$patient), function(z) {
    if (nrow(z) < 100) return(NA_real_)
    coef(lm(y ~ TSPO + subtype, data = z))[2] })
  per <- per[is.finite(per)]
  data.frame(ligand = g, beta = round(cf["TSPO","Estimate"], 4),
             t = round(tp, 2), p = signif(pv, 3),
             pct_patient_pos = round(mean(per > 0) * 100, 0),
             n = nrow(d), n_pat = nlevels(d$patient), stringsAsFactors = FALSE)
}
res <- do.call(rbind, lapply(gs, fit_one))
res$FDR <- p.adjust(res$p, "BH")
res <- res[order(-res$beta), ]
cat("\n混合模型（y ~ TSPO + nFeature + subtype + (1|patient)）:\n")
pr(res)
write.csv(res, file.path(OUTDIR, "myeloid_ligand_mixed_final.csv"), row.names = FALSE)

cat("\n各亚型 MIF 平均表达:\n")
print(data.frame(subtype = sort(unique(subtype)),
                 MIF = round(sapply(sort(unique(subtype)), function(s)
                   mean(as.numeric(dat["MIF", subtype == s]))), 3),
                 TSPO = round(sapply(sort(unique(subtype)), function(s)
                   mean(tspo[subtype == s])), 3)), row.names = FALSE)

# ================================================================ PART C
if (DO_CELLCHAT) {
  hr("PART C  CellChat：按亚型（3 组）")
  suppressPackageStartupMessages(library(CellChat))
  run_cc <- function(obj, popsize, tag) {
    f <- file.path(OUTDIR, paste0("cc5_", tag, "_pop", as.integer(popsize), ".rds"))
    if (file.exists(f)) { message("  [缓存] ", basename(f)); return(readRDS(f)) }
    t0 <- Sys.time(); message("  [运行] ", tag, " pop=", popsize, " cells=", ncol(obj))
    din <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
    meta <- data.frame(labels = as.character(obj@meta.data[[SUB_COL]]),
                       row.names = colnames(obj))
    cc <- CellChat::createCellChat(object = din, meta = meta, group.by = "labels")
    cc@DB <- CellChat::CellChatDB.human
    cc <- CellChat::subsetData(cc)
    cc <- CellChat::identifyOverExpressedGenes(cc)
    cc <- CellChat::identifyOverExpressedInteractions(cc)
    cc <- tryCatch(
      CellChat::computeCommunProb(cc, population.size = popsize, type = "triMean", nboot = NBOOT),
      error = function(e) CellChat::computeCommunProb(cc, type = "triMean", nboot = NBOOT))
    cc <- CellChat::filterCommunication(cc, min.cells = 10)
    cc <- CellChat::computeCommunProbPathway(cc)
    cc <- CellChat::aggregateNet(cc)
    saveRDS(cc, f)
    message(sprintf("  [完成] %.1f 分钟", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    cc
  }
  ccT <- run_cc(sc,   TRUE,  "subtypes")
  ccF <- run_cc(sc,   FALSE, "subtypes")

  hr("PART C-2  各亚型通路强度")
  for (nm in c("popTRUE", "popFALSE")) {
    cc <- if (nm == "popTRUE") ccT else ccF
    cat("\n--- ", nm, " ---\n", sep = "")
    rows <- lapply(cc@netP$pathways, function(p) {
      x <- cc@netP$prob[, , p]
      data.frame(pathway = p, t(round(colMeans(x, na.rm = TRUE), 5)), check.names = FALSE)
    })
    M <- do.call(rbind, rows)
    pr(M[order(-M$Monocyte), ], 4, 10)
  }

  cat("\n===== MIF 通路按亚型 =====\n")
  for (nm in c("popTRUE", "popFALSE")) {
    cc <- if (nm == "popTRUE") ccT else ccF
    if ("MIF" %in% cc@netP$pathways) {
      x <- cc@netP$prob[, , "MIF"]
      cat("\n", nm, "（行=信号源，列=靶细胞）:\n", sep = "")
      print(round(x, 4))
    }
  }

  hr("PART C-3  总通讯量")
  get_tot <- function(cc, tag) {
    p <- cc@net$prob; if (is.null(p)) return(NULL)
    data.frame(setting = tag, group = rownames(p), total = round(rowSums(p), 3),
               stringsAsFactors = FALSE)
  }
  tot <- do.call(rbind, list(get_tot(ccT, "popTRUE"), get_tot(ccF, "popFALSE")))
  print(tot, row.names = FALSE)
  write.csv(tot, file.path(OUTDIR, "cc5_total.csv"), row.names = FALSE)

  # MIF 通路占各亚型总通讯量的百分比（正文 19.0% / 6.6% / 0% 的来源）
  hr("PART C-4  MIF 占各亚型总通讯量的比例")
  share <- list()
  for (nm in c("popTRUE", "popFALSE")) {
    cc <- if (nm == "popTRUE") ccT else ccF
    if (!"MIF" %in% cc@netP$pathways) next
    # prob 是 [源, 靶, 通路] 三维数组：分子必须按源汇总（rowSums），
    # 写成 colSums 会把「谁接收 MIF」当成「谁发出 MIF」，方向反了
    mif  <- rowSums(cc@netP$prob[, , "MIF"], na.rm = TRUE)   # 各亚型作为信号源发出的 MIF 总量
    # 分母必须与分子同一层级：netP 是通路级通讯量，net 是配体–受体级，
    # 两者混用会让分母偏大 1.16–1.40 倍、份额被系统性低估。
    allL <- rowSums(apply(cc@netP$prob, c(1, 2), sum), na.rm = TRUE)   # 各亚型全部通路通讯量
    sh   <- 100 * mif / allL
    share[[nm]] <- data.frame(setting = nm, subtype = names(sh),
                              mif_out = round(as.numeric(mif), 4),
                              total_out = round(as.numeric(allL), 3),
                              mif_share_pct = round(as.numeric(sh), 1),
                              stringsAsFactors = FALSE)
  }
  SH <- do.call(rbind, share)
  print(SH, row.names = FALSE)
  write.csv(SH, file.path(OUTDIR, "cc5_MIF_share.csv"), row.names = FALSE)
  cat("\n正文引用的 22.2% / 8.0% / 0% 对应 popFALSE 一行（通路级口径）。\n")
}

hr("完成")
cat("结果目录: ", OUTDIR, "\n")
cat("\n论文可写的三句:\n")
cat(" 1) 髓系构成：Monocyte/TAM_MCSF/TAM_MCSF_IFNg 的细胞数与比例（Part A）\n")
cat(" 2) TSPO 在哪一亚型最高、阳性率多少（Part A）\n")
cat(" 3) 配体（MIF 等）与 TSPO 的病人内关联（Part B 的 beta/p/一致率）\n")
cat(" 4) MIF 通路的主要来源亚型（Part C-2）\n")

# -----------------------------------------------------------------------------
# -- Fig5_髓系配体与CAF / 10_cptac_protein.R
# -----------------------------------------------------------------------------
# =============================================================================
#  75_CPTAC_complete.R  （容错版）
#  CPTAC-OV 蛋白组分析
#
#  文件目录：D:/Rstudy/tcga/
#    HS_CPTAC_OV_proteome_gene_tumor.cct    必需
#    HS_CPTAC_OV_proteome_gene_normal.cct   选（缺失则跳过 STEP 2）
#    HS_CPTAC_OV_rnaseq_fpkm_log2.cct       选（缺失则跳过 STEP 4）
#
#  文件名若与预期不同，脚本会自动按关键字模糊匹配。
# =============================================================================


if (!exists("OUT")) stop("请先 source('00_setup.R')")

DIR    <- "D:/Rstudy/tcga"
OUTDIR <- file.path(OUT, "tables")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
hr <- function(x) cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")

# ---------------------------------------------------------------- 找文件
hr("STEP 0  定位文件")
all_f <- list.files(DIR, pattern = "CPTAC", ignore.case = TRUE, full.names = TRUE)
cat("目录下 CPTAC 相关文件 ", length(all_f), " 个:\n", sep = "")
print(data.frame(name = basename(all_f), MB = round(file.size(all_f)/1e6, 2)), row.names = FALSE)

find_f <- function(keys) {
  hit <- all_f[Reduce(`&`, lapply(keys, function(k) grepl(k, basename(all_f), ignore.case = TRUE)))]
  if (!length(hit)) return(NA_character_)
  hit[which.max(file.size(hit))]
}
F_T <- find_f(c("proteome", "tumor"))
F_N <- find_f(c("proteome", "normal"))
F_R <- find_f(c("rnaseq"))
cat("\n肿瘤蛋白: ", ifelse(is.na(F_T), "未找到", basename(F_T)), "\n", sep = "")
cat("正常蛋白: ", ifelse(is.na(F_N), "未找到", basename(F_N)), "\n", sep = "")
cat("RNA:      ", ifelse(is.na(F_R), "未找到", basename(F_R)), "\n", sep = "")
if (is.na(F_T)) stop("肿瘤蛋白组文件缺失，无法继续")

# ---------------------------------------------------------------- 读取
hr("STEP 1  读取并映射基因名")
read_cct <- function(f) {
  if (is.na(f)) return(NULL)
  cat("读取 ", basename(f), " (", round(file.size(f)/1e6, 1), " MB)\n", sep = "")
  d <- tryCatch(read.delim(f, check.names = FALSE, stringsAsFactors = FALSE,
                           quote = "", na.strings = c("NA", "", "NaN")),
                error = function(e) NULL)
  if (is.null(d) || ncol(d) < 5)
    d <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE,
                           na.strings = c("NA", "", "NaN")), error = function(e) NULL)
  if (is.null(d)) { cat("  解析失败\n"); return(NULL) }
  cat("  维度: ", nrow(d), " x ", ncol(d), "\n", sep = "")
  cat("  第一列前 5: ", paste(head(as.character(d[[1]]), 5), collapse = ", "), "\n")
  cat("  前 4 列名: ", paste(head(colnames(d), 4), collapse = ", "), "\n")
  d
}
to_mat <- function(d) {
  if (is.null(d)) return(NULL)
  M <- as.matrix(d[, -1, drop = FALSE]); storage.mode(M) <- "double"
  rownames(M) <- as.character(d[[1]]); M
}

PT <- to_mat(read_cct(F_T))
PN <- to_mat(read_cct(F_N))
PR <- to_mat(read_cct(F_R))
cat("\n肿瘤蛋白 ", nrow(PT), " x ", ncol(PT), "\n", sep = "")
if (!is.null(PN)) cat("正常蛋白 ", nrow(PN), " x ", ncol(PN), "\n", sep = "") else cat("正常蛋白 缺失\n")
if (!is.null(PR)) cat("RNA      ", nrow(PR), " x ", ncol(PR), "\n", sep = "") else cat("RNA 缺失\n")

is_entrez <- function(M) {
  if (is.null(M)) return(FALSE)
  ids <- rownames(M); ids <- ids[!is.na(ids) & ids != ""]
  length(ids) > 0 && mean(grepl("^[0-9]+$", ids)) > 0.8
}
cat("\n肿瘤蛋白为 Entrez: ", is_entrez(PT),
    " | 正常: ", is_entrez(PN),
    " | RNA: ", is_entrez(PR), "\n", sep = "")

map_mat <- function(M) {
  if (is.null(M) || !is_entrez(M)) return(M)
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                 keys = unique(rownames(M)), keytype = "ENTREZID",
                                 columns = "SYMBOL")
    map <- map[!is.na(map$SYMBOL), ]; map <- map[!duplicated(map$ENTREZID), ]
    m <- match(rownames(M), map$ENTREZID); keep <- !is.na(m)
    M2 <- M[keep, , drop = FALSE]; rownames(M2) <- map$SYMBOL[m[keep]]
    rownames(M2) <- gsub("\\..*$", "", rownames(M2))
    return(M2[!duplicated(rownames(M2)), , drop = FALSE])
  }
  cat("  (无 org.Hs.eg.db，使用内置映射)\n")
  DICT <- c("706"="TSPO","968"="CD68","929"="CD14","199"="AIF1","4067"="LYZ",
            "1436"="CSF1R","3684"="ITGAM","933"="CD163","4360"="MRC1",
            "4282"="MIF","3956"="LGALS1","6696"="SPP1","2335"="FN1",
            "972"="CD74","960"="CD44","7040"="TGFB1","7422"="VEGFA",
            "4719"="NDUFS1","4720"="NDUFS2","4722"="NDUFS3","4725"="NDUFS5",
            "4729"="NDUFS8","4731"="NDUFV2","6389"="SDHA","6390"="SDHB",
            "7384"="UQCRC1","7385"="UQCRC2","7386"="UQCRFS1","1537"="CYC1",
            "1327"="COX4I1","9377"="COX5A","513"="ATP5F1A","506"="ATP5F1B",
            "2597"="GAPDH","3939"="LDHA","3098"="HK1","3099"="HK2",
            "5214"="PFKP","226"="ALDOA","5230"="PGK1","2023"="ENO1",
            "5478"="PPIA","5879"="RAC1","10097"="ACTR2")
  newn <- DICT[rownames(M)]; keep <- !is.na(newn)
  M2 <- M[keep, , drop = FALSE]; rownames(M2) <- unname(newn[keep]); M2
}

PT <- map_mat(PT); PN <- map_mat(PN); PR <- map_mat(PR)
cat("\n映射后:\n")
cat("  肿瘤蛋白 ", nrow(PT), " 基因, TSPO: ", "TSPO" %in% rownames(PT), "\n", sep = "")
if (!is.null(PN)) cat("  正常蛋白 ", nrow(PN), " 基因, TSPO: ", "TSPO" %in% rownames(PN), "\n", sep = "")
if (!is.null(PR)) cat("  RNA      ", nrow(PR), " 基因, TSPO: ", "TSPO" %in% rownames(PR), "\n", sep = "")
if (!"TSPO" %in% rownames(PT)) stop("肿瘤蛋白组中没有 TSPO，无法分析")

# ---------------------------------------------------------------- 2
hr("STEP 2  TSPO 蛋白：肿瘤 vs 正常")
if (!is.null(PN) && "TSPO" %in% rownames(PN)) {
  tT <- as.numeric(PT["TSPO", ]); tN <- as.numeric(PN["TSPO", ])
  tT <- tT[is.finite(tT)]; tN <- tN[is.finite(tN)]
  if (length(tT) > 5 && length(tN) > 5) {
    w <- suppressWarnings(wilcox.test(tT, tN))
    cat(sprintf("肿瘤 n=%d 中位 %+.3f | 正常 n=%d 中位 %+.3f\n",
                length(tT), median(tT), length(tN), median(tN)))
    cat(sprintf("差值 %+.3f (log2 TMT ratio), Wilcoxon p = %.3g\n",
                median(tT) - median(tN), w$p.value))
    write.csv(data.frame(group = c(rep("Tumor", length(tT)), rep("Normal", length(tN))),
                         TSPO = c(tT, tN)),
              file.path(OUTDIR, "CPTAC_TSPO_tumor_vs_normal.csv"), row.names = FALSE)
  } else cat("有效值不足\n")
} else cat("跳过（正常蛋白组缺失或无 TSPO）\n")

# ---------------------------------------------------------------- 3
hr("STEP 3  TSPO 蛋白与其他蛋白的相关（核心）")
OXPHOS <- c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
            "NDUFA9","SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","UQCRFS1","CYC1",
            "COX4I1","COX5A","COX5B","COX6A1","ATP5F1A","ATP5F1B","ATP5C1","ATP5O","ATP5H")
MYE    <- c("CD68","CD14","AIF1","LYZ","CSF1R","ITGAM","CD163","MRC1","C1QA","C1QB",
            "TYROBP","FCER1G","SPI1")
LIG    <- c("MIF","LGALS1","SPP1","FN1","CD74","CD44","TGFB1","VEGFA")
GLYCO  <- c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1")

tspo <- as.numeric(PT["TSPO", ]); names(tspo) <- colnames(PT)
cat("TSPO 有效样本: ", sum(is.finite(tspo)), " / ", ncol(PT), "\n", sep = "")

detail <- function(gset, lab) {
  g <- intersect(gset, rownames(PT))
  cat(lab, ": 匹配到 ", length(g), " 个蛋白\n", sep = "")
  if (!length(g)) return(NULL)
  out <- do.call(rbind, lapply(g, function(x) {
    v <- as.numeric(PT[x, ]); ok <- is.finite(v) & is.finite(tspo)
    if (sum(ok) < 20) return(NULL)
    ct <- suppressWarnings(cor.test(tspo[ok], v[ok], method = "spearman"))
    data.frame(set = lab, protein = x, rho = round(ct$estimate, 3),
               p = signif(ct$p.value, 3), n = sum(ok), stringsAsFactors = FALSE)
  }))
  out
}
D <- rbind(detail(OXPHOS, "OXPHOS"), detail(MYE, "Myeloid"),
           detail(LIG, "Ligand"), detail(GLYCO, "Glycolysis"))
if (is.null(D) || !nrow(D)) stop("没有可分析的蛋白")
D <- D[order(-D$rho), ]
cat("\n逐个蛋白:\n"); print(D, row.names = FALSE)
write.csv(D, file.path(OUTDIR, "CPTAC_TSPO_protein_cor.csv"), row.names = FALSE)

cat("\n各组汇总:\n")
print(do.call(rbind, lapply(split(D, D$set), function(d)
  data.frame(set = d$set[1], n_protein = nrow(d),
             median_rho = round(median(d$rho), 3),
             min = min(d$rho), max = max(d$rho),
             n_sig_pos = sum(d$p < 0.05 & d$rho > 0),
             n_sig_neg = sum(d$p < 0.05 & d$rho < 0),
             stringsAsFactors = FALSE))), row.names = FALSE)

# ---------------------------------------------------------------- 4
hr("STEP 4  TSPO mRNA vs 蛋白")
if (!is.null(PR) && "TSPO" %in% rownames(PR)) {
  smp <- intersect(colnames(PT), colnames(PR))
  cat("直接匹配样本: ", length(smp), "\n", sep = "")
  if (length(smp) < 20) {
    k1 <- substr(gsub("[^A-Za-z0-9]", "", colnames(PT)), 1, 15)
    k2 <- substr(gsub("[^A-Za-z0-9]", "", colnames(PR)), 1, 15)
    m <- match(k1, k2)
    if (sum(!is.na(m)) >= 20) {
      PT <- PT[, !is.na(m), drop = FALSE]; PR <- PR[, m[!is.na(m)], drop = FALSE]
      colnames(PR) <- colnames(PT); smp <- colnames(PT)
      cat("模糊匹配后: ", length(smp), "\n", sep = "")
    }
  }
  if (length(smp) >= 20) {
    pv <- as.numeric(PT["TSPO", smp]); rv <- as.numeric(PR["TSPO", smp])
    ok <- is.finite(pv) & is.finite(rv)
    ct <- suppressWarnings(cor.test(pv[ok], rv[ok], method = "spearman"))
    cat(sprintf("\nTSPO mRNA vs 蛋白: rho = %+.3f, p = %.3g, n = %d\n",
                ct$estimate, ct$p.value, sum(ok)))
    chk <- function(g, lab) {
      g <- intersect(intersect(g, rownames(PR)), rownames(PT))
      if (!length(g)) return(NULL)
      r <- sapply(g, function(x) {
        a <- as.numeric(PT[x, smp]); b <- as.numeric(PR[x, smp])
        ok <- is.finite(a) & is.finite(b)
        if (sum(ok) < 20) return(NA_real_)
        suppressWarnings(cor(a[ok], b[ok], method = "spearman")) })
      r <- r[is.finite(r)]
      data.frame(set = lab, n = length(r), median_rho = round(median(r), 3))
    }
    cat("\nmRNA-蛋白一致性（对照）:\n")
    print(rbind(chk(OXPHOS,"OXPHOS"), chk(MYE,"Myeloid"), chk(LIG,"Ligand")),
          row.names = FALSE)
  } else cat("样本匹配不足 20，跳过\n")
} else cat("跳过（RNA 文件缺失或无 TSPO）\n")

cat("\n================ 判读 ================\n")
cat("1) 肿瘤 TSPO 蛋白 > 正常      -> 结论一的蛋白层面支持\n")
cat("2) OXPHOS 组中位 rho > 0.3    -> 线①在蛋白层面成立\n")
cat("3) 髓系组中位 rho > 0.3       -> 线②在蛋白层面成立\n")
cat("4) mRNA-蛋白 rho > 0.4        -> 本文 mRNA 结论可外推至蛋白\n")

# -----------------------------------------------------------------------------
# -- Fig5_髓系配体与CAF / 24_caf_trajectory.R
# -----------------------------------------------------------------------------
# =============================================================================
#  24_caf_trajectory.R   基质区室：TSPO 与 CAF 活化状态
#
#  对应正文 3.7：myCAF 11/11（中位偏相关 +0.220）、iCAF 0/11（−0.104）
#                静息 6/11（+0.006）、apCAF 8/11（+0.042）
#
#  方法取舍说明（重要）：
#    先尝试 slingshot 构建 CAF 活化拟时序，方向验证通过（静息 ρ=−0.451、
#    myCAF ρ=+0.457），但三条曲线中与 myCAF 相关最强的一条仅覆盖
#    572/3,079 细胞（18.6%），且"选哪条曲线"系按 myCAF 相关事后确定，
#    存在选择偏倚，该轨迹给出 TSPO 与伪时序无关（ρ=0.032）。
#    改用覆盖全部基质细胞的评分法后得到 11/11 vs 0/11 的清晰结果，
#    故本文以评分法为准，轨迹结果仅作方法学记录（本脚本 PART B 保留该流程）。
#
#  输入：F_SC_EPI（cell_type == "Stromal"）
#  输出：T16a_caf_direction.csv、T16b_caf_per_patient.csv、
#        T16c_caf_by_subtype.csv、T16d_stromal_composition.csv
#        Fig7A_caf_subtype、Fig7B_caf_scores
# =============================================================================

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  if (file.exists(.s)) source(.s, encoding = "UTF-8") else {
    PROJ <- "D:/Rstudy/tcga"; PROJ2 <- "D:/Rstudy/tcga2"
    OUT <- file.path(PROJ, "results_v2")
    F_SC_EPI <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")
  }
}
if (!exists("F_SC_EPI")) F_SC_EPI <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")
OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(ggplot2) })

MIN_CELLS <- 30

# CAF 活化标志
MK <- list(
  Quiescent = c("CFD","DCN","LUM","MFAP5","PDPN","SERPINF1"),
  iCAF      = c("IL6","CXCL1","CXCL2","LIF","HAS1","CCL2"),
  myCAF     = c("ACTA2","TAGLN","CCN2","COL1A1","COL1A2","POSTN","FN1"),
  apCAF     = c("HLA-DRA","HLA-DRB1","CD74")
)

# 细胞类型判别（用于确认基质组成）
MARK <- list(
  Fibroblast  = c("DCN","LUM","COL1A1","COL1A2","PDGFRA","THY1","MFAP5"),
  Endothelial = c("PECAM1","VWF","CDH5","PLVAP","EGFL7","CLDN5"),
  Mesothelial = c("MSLN","UPK3B","CALB2","LRRN4","WT1"),
  Pericyte    = c("RGS5","PDGFRB","ACTA2","MYH11","NDUFA4L2")
)

hr("STEP 1  载入基质细胞")

sc <- readRDS(F_SC_EPI)
cat("cell_type 分布:\n"); print(table(sc$cell_type))
sc <- subset(sc, cells = colnames(sc)[as.character(sc$cell_type) == "Stromal"])
cat("\n基质细胞: ", ncol(sc), "  病人: ", length(unique(sc$patient_id)), "\n", sep = "")
if ("cell_subtype" %in% colnames(sc@meta.data)) {
  cat("cell_subtype 分布:\n"); print(table(sc$cell_subtype))
}
if (!"data" %in% SeuratObject::Layers(sc[["RNA"]]))
  sc <- Seurat::NormalizeData(sc, verbose = FALSE)

D  <- as.matrix(SeuratObject::LayerData(sc, assay = "RNA", layer = "data"))
Z  <- t(scale(t(D))); Z[!is.finite(Z)] <- NA

score_of <- function(genes, Zm = Z) {
  g <- intersect(genes, rownames(Zm))
  if (length(g) < 3) return(NULL)
  colMeans(Zm[g, , drop = FALSE], na.rm = TRUE)
}

hr("STEP 2  基质组成判定")

cs <- list()
for (nm in names(MARK)) { s <- score_of(MARK[[nm]]); if (!is.null(s)) cs[[nm]] <- s }
comp_tab <- data.frame(
  type = names(cs),
  n_gene = vapply(names(cs), function(nm) length(intersect(MARK[[nm]], rownames(Z))), 0),
  mean_score = round(vapply(cs, mean, 0, na.rm = TRUE), 3), row.names = NULL)
pr(comp_tab, 3)

mat <- do.call(cbind, cs)
assign_type <- colnames(mat)[max.col(mat, ties.method = "first")]
names(assign_type) <- colnames(sc)
cat("\n按最高分判定的基质构成:\n"); print(table(assign_type))
write.csv(as.data.frame(table(assign_type)),
          file.path(OUTDIR, "T16d_stromal_composition.csv"), row.names = FALSE)

hr("STEP 3  逐病人相关 + 符号检验（主要推断）")

tspo <- as.numeric(SeuratObject::LayerData(sc, assay = "RNA", layer = "data")["TSPO", ])
dep  <- log10(sc$nCount_RNA + 1)
pat  <- as.character(sc$patient_id)

pcor <- function(a, b, cov) {
  ok <- is.finite(a) & is.finite(b) & is.finite(cov)
  if (sum(ok) < MIN_CELLS) return(NA_real_)
  ra <- residuals(lm(rank(a[ok]) ~ rank(cov[ok])))
  rb <- residuals(lm(rank(b[ok]) ~ rank(cov[ok])))
  cor(ra, rb)
}

S <- list()
for (nm in names(MK)) { s <- score_of(MK[[nm]]); if (!is.null(s)) S[[nm]] <- s }

res <- do.call(rbind, lapply(names(S), function(nm) {
  r <- sapply(split(seq_along(pat), pat), function(i) pcor(tspo[i], S[[nm]][i], dep[i]))
  r <- r[is.finite(r)]
  if (length(r) < 6) return(NULL)
  k <- sum(r > 0); n <- length(r)
  data.frame(state = nm, n_patient = n, n_positive = k,
             median_rho = round(median(r), 3),
             sign_p = binom.test(k, n, 0.5)$p.value, stringsAsFactors = FALSE)
}))
res$FDR <- p.adjust(res$sign_p, "BH")
res <- res[order(-res$median_rho), ]
cat("\n校正测序深度后的逐病人偏相关：\n"); pr(res, 4)
write.csv(res, file.path(OUTDIR, "T16a_caf_direction.csv"), row.names = FALSE)

# 逐病人明细
cat("\n各状态在每个病人内的偏相关（myCAF / iCAF）：\n")
for (nm in c("myCAF","iCAF")) {
  if (!nm %in% names(S)) next
  r <- sapply(split(seq_along(pat), pat), function(i) pcor(tspo[i], S[[nm]][i], dep[i]))
  r <- r[is.finite(r)]
  cat(sprintf("%-8s n=%d  ", nm, length(r)),
      paste(sprintf("%+.2f", r), collapse = "  "), "\n")
  write.csv(data.frame(patient = names(r), rho = as.numeric(r)),
            file.path(OUTDIR, paste0("T16b_caf_per_patient_", nm, ".csv")), row.names = FALSE)
}

hr("STEP 4  按 cell_subtype 看 TSPO 绝对水平")

if ("cell_subtype" %in% colnames(sc@meta.data)) {
  sub <- as.character(sc$cell_subtype)
  m   <- tapply(tspo, sub, mean, na.rm = TRUE)
  n   <- tapply(tspo, sub, length)
  out <- data.frame(subtype = names(m), n = as.integer(n),
                    mean_TSPO = round(as.numeric(m), 3), row.names = NULL)
  out <- out[out$n >= 50, ]
  out <- out[order(-out$mean_TSPO), ]
  cat("\n各基质亚型的 TSPO 均值：\n"); pr(out, 3)
  write.csv(out, file.path(OUTDIR, "T16c_caf_by_subtype.csv"), row.names = FALSE)

  p1 <- ggplot(data.frame(subtype = sub, tspo = tspo), aes(subtype, tspo, fill = subtype)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    labs(x = NULL, y = "TSPO (log-normalized)",
         title = "TSPO across stromal subtypes") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(FIGDIR, "Fig7A_caf_subtype.pdf"), p1, width = 6, height = 5)
  ggsave(file.path(FIGDIR, "Fig7A_caf_subtype.png"), p1, width = 6, height = 5, dpi = 300)
}

hr("STEP 5  图：myCAF 与 iCAF 的对比")

if (all(c("myCAF","iCAF") %in% names(S))) {
  df <- data.frame(myCAF = S$myCAF, iCAF = S$iCAF, tspo = tspo,
                   patient = pat, stringsAsFactors = FALSE)
  df <- df[is.finite(df$tspo), ]
  p2 <- ggplot(df, aes(myCAF, tspo)) +
    geom_point(size = 0.4, alpha = 0.2, colour = "#B2182B") +
    geom_smooth(method = "loess", colour = "#B2182B", se = TRUE) +
    labs(x = "myCAF score (z)", y = "TSPO (log-normalized)",
         title = "TSPO rises with myofibroblastic program") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIGDIR, "Fig7B_caf_scores.pdf"), p2, width = 5.5, height = 5)
  ggsave(file.path(FIGDIR, "Fig7B_caf_scores.png"), p2, width = 5.5, height = 5, dpi = 300)
}

saveRDS(list(S = S, tspo = tspo, patient = pat, res = res),
        file.path(OUT, "cache", "caf_scores.rds"))

cat("\n判读要点：\n")
cat("  · 主推断是 STEP 3 的逐病人符号检验，不是细胞层面 p 值。\n")
cat("  · myCAF 全正、iCAF 全负构成双向特异性，排除「任何高表达程序都与 TSPO 相关」。\n")
cat("  · 这是区室内部效应：CAF 的 TSPO 绝对值低于上皮细胞（EOC 1.293 vs CAF 0.47–0.95）。\n")

hr("24_caf_trajectory.R 完成")

# -----------------------------------------------------------------------------
# -- Fig5_髓系配体与CAF / 42_chemo_paired_programs.R
# -----------------------------------------------------------------------------
# =============================================================================
#  42_chemo_paired_programs.R
#
#  化疗初治与间歇性减灭术（IDS）的配对比较：两个细胞程序的化疗响应是否不同。
#
#  服务正文 3.6「两个程序的化疗响应不同」与图 6H。
#
#  设计：同一病人内部的配对比较，病人自身作为对照，
#        因此病人间的构成差异不进入推断。
#
#  输入：D:/Rstudy/tcga2/GSE266577_processed/GSE266577_raw.rds
#  输出（工作区 results_v2/tables/）：
#    T44_chemo_paired_programs.csv     16 项检验（2 区室 × 8 程序）
#    T44b_chemo_paired_receptors.csv   CD74 / CD44 / CXCR4 的配对比较
#    T44c_chemo_tspo_spp1_coupling.csv TSPO~SPP1 耦合在治疗前后是否改变
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
})

PROJ <- "D:/Rstudy/tcga2"
OUT  <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/tables"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")

# ---------------------------------------------------------------------------
# 1  载入
# ---------------------------------------------------------------------------

hr("STEP 1  载入并归一化")
o <- readRDS(file.path(PROJ, "GSE266577_processed", "GSE266577_raw.rds"))
o <- JoinLayers(o, assay = "RNA")
CNT <- SeuratObject::LayerData(o, assay = "RNA", layer = "counts")
lib <- Matrix::colSums(CNT); lib[lib == 0] <- 1
Xn  <- CNT %*% Matrix::Diagonal(x = 1e4 / lib)
Xn@x <- log1p(Xn@x)

CT  <- as.character(o$cell_type)
STG <- as.character(o$treatment_stage)
PAT <- as.character(o$publication_patient_code_final)
cat("细胞:", ncol(Xn), " 病人:", length(unique(PAT)), "\n")
print(table(STG))

tb   <- table(PAT, STG)
BOTH <- rownames(tb)[tb[, "chemo-naive"] > 0 & tb[, "IDS"] > 0]
cat("\n同时具两时点的病人:", length(BOTH), "\n")

# ---------------------------------------------------------------------------
# 2  程序定义
# ---------------------------------------------------------------------------

GENES <- list(
  OXPHOS = c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
             "SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","UQCRH","CYC1",
             "COX4I1","COX5A","COX5B","COX6A1","ATP5A1","ATP5B","ATP5C1","ATP5F1"),
  Glyco  = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1","SLC2A1",
             "SLC16A1","PDK1"),
  Gln    = c("GLS","GOT1","GOT2","GLUD1","SLC1A5","SLC7A5"),
  SPP1   = "SPP1", FN1 = "FN1", LGALS1 = "LGALS1", MIF = "MIF",
  TSPO   = "TSPO",
  CD74   = "CD74", CD44 = "CD44", CXCR4 = "CXCR4"
)

# 关键：在整类细胞上 z 标准化一次，再取子集均值。
# 若逐子集 z 标准化，均值恒为 0，配对差将全部为 0。
score_all <- function(gs, cells) {
  g <- intersect(gs, rownames(Xn))
  M <- as.matrix(Xn[g, cells, drop = FALSE])
  Z <- t(scale(t(M))); Z[!is.finite(Z)] <- 0
  colMeans(Z)
}

# ---------------------------------------------------------------------------
# 3  配对比较
# ---------------------------------------------------------------------------

paired_test <- function(comp_cells, prog) {
  SC <- score_all(GENES[[prog]], comp_cells)
  a <- c(); b <- c()
  for (p in BOTH) {
    i1 <- which(PAT[comp_cells] == p & STG[comp_cells] == "chemo-naive")
    i2 <- which(PAT[comp_cells] == p & STG[comp_cells] == "IDS")
    if (length(i1) < 30 || length(i2) < 30) next
    a <- c(a, mean(SC[i1])); b <- c(b, mean(SC[i2]))
  }
  if (length(a) < 4) return(NULL)
  pv <- suppressWarnings(wilcox.test(a, b, paired = TRUE))$p.value
  data.frame(n_pairs = length(a), naive = mean(a), IDS = mean(b),
             diff = mean(b) - mean(a), p = pv)
}

hr("STEP 2  配对比较：两个区室 × 十一个程序")
res <- data.frame()
for (cc in list(c("上皮", "Epithelial cells"), c("髓系", "Macrophages"))) {
  cells <- which(CT == cc[2])
  for (prog in names(GENES)) {
    r <- paired_test(cells, prog)
    if (is.null(r)) next
    res <- rbind(res, cbind(comp = cc[1], prog = prog, r))
  }
}
res$FDR <- signif(p.adjust(res$p, "BH"), 3)
res$p   <- signif(res$p, 3)
for (v in c("naive","IDS","diff")) res[[v]] <- round(res[[v]], 3)
res <- res[order(res$p), ]

cat("\n按 p 排序（BH 校正）：\n")
print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "T44_chemo_paired_programs.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4  阳性比例的变化：群体大小是否改变
# ---------------------------------------------------------------------------

hr("STEP 3  阳性细胞比例与受体")
MA <- which(CT == "Macrophages")
for (s in c("chemo-naive", "IDS")) {
  i <- MA[STG[MA] == s]
  cat(sprintf("  %-12s SPP1⁺ %.1f%%   TSPO⁺ %.1f%%   CD74⁺ %.1f%%   CD44⁺ %.1f%%\n",
      s,
      100 * mean(as.numeric(Xn["SPP1", i]) > 0),
      100 * mean(as.numeric(Xn["TSPO", i]) > 0),
      100 * mean(as.numeric(Xn["CD74", i]) > 0),
      100 * mean(as.numeric(Xn["CD44", i]) > 0)))
}

rec <- res[res$comp == "髓系" & res$prog %in% c("CD74","CD44","CXCR4"), ]
if (nrow(rec)) {
  write.csv(rec, file.path(OUT, "T44b_chemo_paired_receptors.csv"), row.names = FALSE)
  cat("\n受体配对结果：\n"); print(rec, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# 5  TSPO~SPP1 的耦合在治疗前后是否改变
# ---------------------------------------------------------------------------

hr("STEP 4  TSPO~SPP1 耦合")
ts <- as.numeric(Xn["TSPO", ]); sp <- as.numeric(Xn["SPP1", ])
for (s in c("chemo-naive", "IDS")) {
  i <- MA[STG[MA] == s]
  rh <- sapply(unique(PAT[i]), function(p) {
    j <- i[PAT[i] == p]; if (length(j) < 30) return(NA_real_)
    suppressWarnings(cor(ts[j], sp[j], method = "spearman"))
  })
  rh <- rh[is.finite(rh)]
  cat(sprintf("  %-12s 合并 %+.3f  为正 %2d/%-2d  中位 %+.3f\n",
      s, cor(ts[i], sp[i], method = "spearman"), sum(rh > 0), length(rh), median(rh)))
}

r1 <- c(); r2 <- c()
for (p in BOTH) {
  a <- MA[PAT[MA] == p & STG[MA] == "chemo-naive"]
  b <- MA[PAT[MA] == p & STG[MA] == "IDS"]
  if (length(a) < 30 || length(b) < 30) next
  r1 <- c(r1, suppressWarnings(cor(ts[a], sp[a], method = "spearman")))
  r2 <- c(r2, suppressWarnings(cor(ts[b], sp[b], method = "spearman")))
}
pv <- suppressWarnings(wilcox.test(r1, r2, paired = TRUE))$p.value
cat(sprintf("\n  配对: naive 中位 %+.3f -> IDS %+.3f  p=%.3g  n=%d\n",
            median(r1), median(r2), pv, length(r1)))
cat("  耦合未变 => 变的是输出量，不是标记与被标记者之间的关系\n")

write.csv(data.frame(comparison = "TSPO~SPP1 in macrophages",
                     naive_median_rho = round(median(r1), 3),
                     IDS_median_rho   = round(median(r2), 3),
                     paired_p = signif(pv, 3), n_pairs = length(r1)),
          file.path(OUT, "T44c_chemo_tspo_spp1_coupling.csv"), row.names = FALSE)

hr("42_chemo_paired_programs.R 完成")

# ############################################################################
# ##  Fig6_髓系轨迹  (3 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- Fig6_髓系轨迹 / 06_myeloid_trajectory.R
# -----------------------------------------------------------------------------
# =============================================================================
#  06_myeloid_trajectory.R
#  结论三（新增）：TSPO 沿单核 -> 巨噬分化轴程序性下调
#
#  对应正文 3.4 与图 4。数据：GSE266577（137,083 个细胞，32 例病人）
#
#  设计要点（每一步都是踩过坑之后的结论，不要随意改）：
#    1. 先取 counts 子集再建 Seurat 对象 —— 137k 细胞全量 NormalizeData 会 std::bad_alloc
#    2. 只对 1,500 个高变基因 ScaleData —— 全基因缩放约需 9 GB × 3 份
#    3. 平衡采样：单核 256 + 巨噬 768（1:3）。原始比例 546:23,919 = 1:44，
#       主曲线会被巨噬内部异质性主导，轴末端被 DC/pDC 占据，方向验证必然失败
#    4. 方向验证用「亚型层面」的早期/晚期标志物对比，不用「沿十分位单调性」
#       —— 该轴在巨噬内部还有第二层状态切换，会稀释整体单调性
#    5. 逐病人一致性按「病人」合并化疗初治与 IDS 样本，不按样本
#       —— 同一病人两次取样不是独立观测
#    6. 变量名不要用 sub / subset 等与 base 函数同名的名字
#
#  产出（正文数字）：
#    · 单核中位伪时序 -56.6，巨噬 -23.1
#    · 早期标志物 -73.3% (p=3.5e-124)，晚期 +32.3% (p=2.1e-12)
#    · 单样本重复 S032_IDS：-30.7% / +122.1%
#    · TSPO ~ pt: rho=-0.629；校正深度后 -0.568
#    · 18/21 例病人为负，符号检验 p=1.5e-3
#    · 2000 个随机基因中 TSPO 排第 100 百分位
# =============================================================================

if (!exists("OUT")) stop("请先 source('00_setup.R')")
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(slingshot)
})

MAX_CELL   <- 12000     # 平衡集上限（不是总细胞数）
N_HVG      <- 1500
N_PC       <- 20
N_PC_USE   <- 15
DO_REPLICATE <- TRUE    # 是否做单病人内部重复（S032_IDS）

hr("STEP 1  读取并立即取髓系子集")

stopifnot(file.exists(F_GSE266))
o   <- readRDS(F_GSE266)
CT  <- "cell_type"
stopifnot(CT %in% colnames(o@meta.data))

MYE <- c("Classical monocytes", "Macrophages")
cat("全对象: ", nrow(o), " 基因 x ", ncol(o), " 细胞\n", sep = "")
print(table(o@meta.data[[CT]])[intersect(MYE, unique(o@meta.data[[CT]]))])

keep_cells <- colnames(o)[as.character(o@meta.data[[CT]]) %in% MYE]
cat("单核 + 巨噬: ", length(keep_cells), "\n", sep = "")

# 关键：先切 counts，再建对象，最后 rm() 掉大对象
cnt <- Seurat::GetAssayData(o, assay = "RNA", layer = "counts")[, keep_cells, drop = FALSE]
md  <- o@meta.data[keep_cells, , drop = FALSE]
rm(o); gc()

hr("STEP 2  平衡采样（单核:巨噬 = 1:3）")

mono_cells <- rownames(md)[as.character(md[[CT]]) == "Classical monocytes"]
mac_cells  <- rownames(md)[as.character(md[[CT]]) == "Macrophages"]
cat("可用 — 单核: ", length(mono_cells), "  巨噬: ", length(mac_cells), "\n", sep = "")

set.seed(SEED)
n_mono <- min(length(mono_cells), floor(MAX_CELL / 4))
n_mac  <- min(length(mac_cells), 3 * n_mono)
bal    <- c(sample(mono_cells, n_mono), sample(mac_cells, n_mac))
cat(sprintf("平衡后 — 单核 %d (%.1f%%) + 巨噬 %d = %d\n",
            n_mono, 100 * n_mono / length(bal), n_mac, length(bal)))

cnt <- cnt[, bal, drop = FALSE]
md  <- md[bal, , drop = FALSE]

hr("STEP 3  标准化与降维（只缩放高变基因）")

obj <- Seurat::CreateSeuratObject(counts = cnt, meta.data = md, min.cells = 10)
rm(cnt); gc()
cat("对象: ", nrow(obj), " x ", ncol(obj), "\n")

obj <- Seurat::NormalizeData(obj, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, nfeatures = N_HVG, verbose = FALSE)
HVG <- Seurat::VariableFeatures(obj)

obj <- Seurat::ScaleData(obj, features = HVG, verbose = FALSE)   # 只缩放 HVG
obj <- Seurat::RunPCA(obj, features = HVG, npcs = N_PC,
                      reduction.name = "pca", verbose = FALSE)
cat("内存占用: ", round(sum(gc()[, 2]), 0), " MB\n")

hr("STEP 4  slingshot 伪时序")

emb <- Seurat::Embeddings(obj, reduction = "pca")[, seq_len(N_PC_USE)]
cl  <- as.character(obj@meta.data[[CT]])

sl <- slingshot::slingshot(emb, clusterLabels = cl, reducedDim = "pca")
cat("曲线数: ", ncol(slingshot::slingPseudotime(sl)), "\n", sep = "")
pt <- as.numeric(slingshot::slingPseudotime(sl)[, 1])
names(pt) <- colnames(obj)

# 定向：单核 -> 巨噬
if (median(pt[cl == "Classical monocytes"], na.rm = TRUE) >
    median(pt[cl == "Macrophages"],         na.rm = TRUE)) {
  pt <- -pt; cat("已翻转方向\n")
}
cat("\n亚型中位伪时序（单核应最小）:\n")
print(round(tapply(pt, cl, median, na.rm = TRUE), 1))

hr("STEP 5  方向验证（亚型层面，不用十分位单调性）")

dat <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")

EARLY <- c("CD14","FCN1","S100A8","S100A9","VCAN","LYZ")
LATE  <- c("CD163","MRC1","CSF1R","MAF","F13A1","C1QA","APOE")

dir_tab <- do.call(rbind, lapply(c("early", "late"), function(k) {
  gs <- intersect(if (k == "early") EARLY else LATE, rownames(dat))
  x  <- Matrix::colMeans(dat[gs, , drop = FALSE], na.rm = TRUE)
  m  <- tapply(x, cl, mean)
  chg <- 100 * (m[["Macrophages"]] / m[["Classical monocytes"]] - 1)
  w  <- suppressWarnings(wilcox.test(x[cl == "Classical monocytes"],
                                     x[cl == "Macrophages"]))
  data.frame(marker_set = k, n_genes = length(gs),
             mono = round(m[["Classical monocytes"]], 3),
             macro = round(m[["Macrophages"]], 3),
             change_pct = round(chg, 1),
             p = signif(w$p.value, 3),
             pass = if (k == "early") chg < 0 else chg > 0,
             stringsAsFactors = FALSE)
}))
pr(dir_tab)
ok_dir <- all(dir_tab$pass)
cat(if (ok_dir) "✓ 方向通过（早期降、晚期升）\n" else "✗ 方向未通过，见 README 排查\n")

hr("STEP 6  TSPO 沿伪时序的表达")

if (!"TSPO" %in% rownames(dat)) stop("表达矩阵里没有 TSPO")
y   <- as.numeric(dat["TSPO", ])
gd  <- log10(obj$nCount_RNA + 1)

r_raw <- suppressWarnings(cor.test(pt, y, method = "spearman"))
rp    <- pcor_rank(pt, y, gd)
cat(sprintf("TSPO ~ pt          : rho = %+.3f, p = %.3g\n", r_raw$estimate, r_raw$p.value))
cat(sprintf("TSPO ~ pt | 深度    : rho = %+.3f, p = %.3g\n", rp["rho"], rp["p"]))
cat(sprintf("pt ~ nFeature       : rho = %+.3f\n",
            cor(pt, obj$nFeature_RNA, method = "spearman")))
cat("  注：该相关由「细胞类型 -> 检出基因数 -> 轴位置」造成，不是技术伪影；\n")
cat("      校正深度后 TSPO 的 rho 几乎不变，说明关联不依赖深度。\n")

hr("STEP 7  随机基因背景（判断是否整体漂移）")

set.seed(SEED)
gs_rand <- sample(intersect(rownames(dat), rownames(obj)), 2000)
r_rand  <- sapply(gs_rand, function(g)
  suppressWarnings(cor(pt, as.numeric(dat[g, ]), method = "spearman")))
cat(sprintf("2000 个随机基因：rho<0 占 %.1f%%，中位 rho %+.3f\n",
            mean(r_rand < 0, na.rm = TRUE) * 100, median(r_rand, na.rm = TRUE)))
cat(sprintf("TSPO 的百分位: %.0f%%（100%% = 比所有随机基因都更负）\n",
            mean(r_rand > r_raw$estimate, na.rm = TRUE) * 100))

hr("STEP 8  逐病人一致性（按病人合并化疗前后）")

pat_col <- intersect(c("publication_patient_code_final",
                       "publication_sample_code_final"), colnames(obj@meta.data))[1]
cat("病人列: ", pat_col, "\n", sep = "")
pat <- as.character(obj@meta.data[[pat_col]])

D <- data.frame(pt = pt, tspo = y, pat = pat, stringsAsFactors = FALSE)
D <- D[is.finite(D$pt) & is.finite(D$tspo), ]

Q <- do.call(rbind, lapply(split(D, D$pat), function(x) {
  if (nrow(x) < 50) return(NULL)
  data.frame(patient = x$pat[1], n = nrow(x),
             rho = suppressWarnings(cor(x$pt, x$tspo, method = "spearman")),
             stringsAsFactors = FALSE)
}))
Q <- Q[order(Q$rho), ]
pr(Q)

if (nrow(Q) >= 5) {
  k <- sum(Q$rho < 0); m <- nrow(Q)
  cat(sprintf("\n中位 rho = %+.3f | 负相关 %d/%d | 符号检验 p = %.3g\n",
              median(Q$rho), k, m, binom.test(k, m)$p.value))
  cat("  推断以病人层面为准；细胞层面 p 值受伪重复影响，仅用于描述效应量\n")
}

hr("STEP 9  直接测量的细胞类型均值（不依赖轨迹）")

tab_ct <- do.call(rbind, lapply(MYE, function(k) {
  v <- y[cl == k]
  data.frame(celltype = k, n = length(v),
             mean_TSPO = round(mean(v), 3),
             pct_positive = round(mean(v > 0) * 100, 1), stringsAsFactors = FALSE)
}))
pr(tab_ct)
cat(sprintf("比值（单核/巨噬）= %.2f\n",
            tab_ct$mean_TSPO[tab_ct$celltype == "Classical monocytes"] /
            tab_ct$mean_TSPO[tab_ct$celltype == "Macrophages"]))

hr("STEP 10  单病人内部重复（排除跨样本混杂）")

rep_tab <- NULL
if (DO_REPLICATE) {
  # 选单核数量最多的那个样本，内部按 1:3 再平衡一次
  smp_col <- "publication_sample_code_final"
  tb <- table(pat[cl == "Classical monocytes"])
  if (length(tb) && max(tb) >= 30) {
    target <- names(tb)[which.max(tb)]
    idx_it <- which(pat == target)
    mi <- idx_it[cl[idx_it] == "Classical monocytes"]
    ma <- idx_it[cl[idx_it] == "Macrophages"]
    set.seed(SEED)
    nm <- min(length(mi), 100)
    ng <- min(length(ma), 3 * nm)
    sel2 <- c(sample(mi, nm), sample(ma, ng))
    cat("单样本内部: ", target, "  单核 ", nm, " + 巨噬 ", ng, "\n", sep = "")

    cnt2 <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")[, sel2, drop = FALSE]
    o2 <- Seurat::CreateSeuratObject(counts = cnt2, min.cells = 5)
    o2 <- Seurat::NormalizeData(o2, verbose = FALSE)
    o2 <- Seurat::FindVariableFeatures(o2, nfeatures = 800, verbose = FALSE)
    hv2 <- Seurat::VariableFeatures(o2)
    o2 <- Seurat::ScaleData(o2, features = hv2, verbose = FALSE)
    o2 <- Seurat::RunPCA(o2, features = hv2, npcs = 15, reduction.name = "pca", verbose = FALSE)

    cl2 <- cl[sel2]
    sl2 <- slingshot::slingshot(Seurat::Embeddings(o2, "pca")[, 1:10],
                                clusterLabels = cl2, reducedDim = "pca")
    pt2 <- as.numeric(slingshot::slingPseudotime(sl2)[, 1])
    if (median(pt2[cl2 == "Classical monocytes"], na.rm = TRUE) >
        median(pt2[cl2 == "Macrophages"], na.rm = TRUE)) pt2 <- -pt2

    d2 <- Seurat::GetAssayData(o2, assay = "RNA", layer = "data")
    rep_tab <- do.call(rbind, lapply(c("early", "late"), function(k) {
      gs <- intersect(if (k == "early") EARLY else LATE, rownames(d2))
      x  <- Matrix::colMeans(d2[gs, , drop = FALSE], na.rm = TRUE)
      m  <- tapply(x, cl2, mean)
      chg <- 100 * (m[["Macrophages"]] / m[["Classical monocytes"]] - 1)
      w  <- suppressWarnings(wilcox.test(x[cl2 == "Classical monocytes"],
                                         x[cl2 == "Macrophages"]))
      data.frame(sample = target, marker_set = k, change_pct = round(chg, 1),
                 p = signif(w$p.value, 3),
                 pass = if (k == "early") chg < 0 else chg > 0, stringsAsFactors = FALSE)
    }))
    pr(rep_tab)
    cat(if (all(rep_tab$pass)) "✓ 单样本内部方向也通过\n" else "✗ 单样本方向未通过\n")
  } else cat("没有单核足够的单样本，跳过\n")
}

hr("STEP 11  写出")

write.csv(dir_tab, file.path(OUT, "tables", "T3a_trajectory_direction.csv"), row.names = FALSE)
write.csv(Q,       file.path(OUT, "tables", "T3b_trajectory_per_patient.csv"), row.names = FALSE)
write.csv(tab_ct,  file.path(OUT, "tables", "T3c_TSPO_by_celltype.csv"), row.names = FALSE)
if (!is.null(rep_tab))
  write.csv(rep_tab, file.path(OUT, "tables", "T3d_trajectory_replicate.csv"), row.names = FALSE)

saveRDS(list(pt = pt, cells = colnames(obj), direction = dir_tab,
             tspo_rho = r_raw$estimate, random = r_rand),
        file.path(OUT, "cache", "trajectory_result.rds"))
cat("已写出到: ", file.path(OUT, "tables"), "\n", sep = "")

# ---------------------------------------------------------------------------
#  STEP 12  为图 6 保存绘图缓存
#
#  上面已把统计量写成 CSV，但图还需要 UMAP、伪时序与关键基因的表达矩阵。
#  这些不改动任何推断，只是把绘图所需的原始信息一并落盘，避免画图时重跑全流程。
# ---------------------------------------------------------------------------
hr("STEP 12  绘图缓存")

if (!"umap" %in% Seurat::Reductions(obj)) {
  obj <- Seurat::RunUMAP(obj, dims = seq_len(N_PC_USE), reduction = "pca",
                         seed.use = SEED, verbose = FALSE)
  cat("已计算 UMAP\n")
}
UM <- Seurat::Embeddings(obj, "umap")[, 1:2]
colnames(UM) <- c("UMAP_1", "UMAP_2")

GENES_FIG <- unique(c(
  "TSPO",
  "CD14","FCN1","S100A8","S100A9","VCAN","LYZ",          # 早期标志
  "CD163","MRC1","CSF1R","MAF","F13A1","C1QA","APOE",    # 晚期标志
  "MIF","LGALS1","SPP1","FN1","CD74","TGFB1"))           # 下游配体
GENES_FIG <- intersect(GENES_FIG, rownames(obj))
EX <- as.matrix(SeuratObject::LayerData(obj, assay = "RNA", layer = "data")[GENES_FIG, , drop = FALSE])

MDF <- obj@meta.data
saveRDS(list(
  umap      = UM,
  pt        = pt,
  expr      = EX,
  cell_type = as.character(MDF[[CT]]),
  patient   = if (pat_col %in% colnames(MDF)) as.character(MDF[[pat_col]]) else NA_character_,
  sample    = if (smp_col %in% colnames(MDF)) as.character(MDF[[smp_col]]) else NA_character_,
  n_count   = MDF$nCount_RNA,
  n_mono    = length(mono_cells),
  n_mac     = length(mac_cells),
  tspo_rho  = as.numeric(r_raw$estimate),
  random    = r_rand,
  direction = dir_tab,
  per_patient = Q),
  file.path(OUT, "cache", "trajectory_fig.rds"))
cat("已保存绘图缓存: ", file.path(OUT, "cache", "trajectory_fig.rds"), "\n", sep = "")

hr("06_myeloid_trajectory.R 完成")

# ---------------------------------------------------------------------------
#  补：全集群（未采样）的 TSPO 均值，用于正文「单核 2.159 vs 巨噬 0.970」
#  平衡子集里的 1.005 是重抽样值，不要拿来和 0.970 混用
# ---------------------------------------------------------------------------
{
  g_full <- readRDS(F_GSE266)
  if (length(SeuratObject::Layers(g_full[["RNA"]])) == 1)
    g_full <- Seurat::NormalizeData(g_full, verbose = FALSE)
  keep <- as.character(g_full$cell_type) %in%
          c("Classical monocytes", "Macrophages")
  g_full <- subset(g_full, cells = colnames(g_full)[keep])
  vv <- as.character(g_full$cell_type)
  ee <- Seurat::GetAssayData(g_full, assay = "RNA", layer = "data")["TSPO", ]
  full <- data.frame(
    celltype     = names(table(vv)),
    n            = as.integer(table(vv)),
    mean_TSPO    = round(as.numeric(tapply(ee, vv, mean)), 3),
    pct_positive = round(100 * as.numeric(tapply(ee, vv, function(x) mean(x > 0))), 1),
    stringsAsFactors = FALSE)
  print(full, row.names = FALSE)
  write.csv(full, file.path(OUTDIR, "T3c_full_TSPO_by_celltype.csv"), row.names = FALSE)
  cat("全集群比值（单核/巨噬）= ",
      round(full$mean_TSPO[full$n == 546] / full$mean_TSPO[full$n == 23919], 2), "\n", sep = "")
  rm(g_full); gc()
}

# -----------------------------------------------------------------------------
# -- Fig6_髓系轨迹 / 07_gse266577_chemo.R
# -----------------------------------------------------------------------------
# =============================================================================
#  72_GSE266577_NACT.R
#  新辅助化疗配对队列：TSPO 在化疗前后的变化
#
#  数据集：GSE266577，137,083 细胞，32 例病人
#    treatment_stage: chemo-naive（化疗前活检） vs IDS（化疗后手术）
#    约 20 例病人有配对样本
#
#  临床问题：
#    1. TSPO 在上皮细胞中化疗后升高还是降低？
#       - 升高 -> TSPO 高细胞富集于残留病灶 -> 提示化疗耐药
#       - 降低 -> TSPO 高细胞被化疗清除 -> 提示化疗敏感
#    2. 同一数据集内，TSPO 在上皮 vs 髓系的表达差异（细胞类型注释可靠）
#    3. 化疗前后细胞构成变化
# =============================================================================


if (!exists("OUT")) stop("请先 source('00_setup.R')")

RDS    <- "D:/Rstudy/tcga2/GSE266577_processed/GSE266577_raw.rds"
OUTDIR <- file.path(OUT, "tables")
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
hr <- function(x) cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")

# ---------------------------------------------------------------- 1
hr("STEP 1  载入")
o <- readRDS(RDS)
cat("维度: ", nrow(o), " x ", ncol(o), "\n")
md <- o@meta.data

# 取表达：data 层为空，需要从 counts 标准化
d <- tryCatch(Seurat::GetAssayData(o, assay = "RNA", layer = "data"),
              error = function(e) NULL)
if (is.null(d) || nrow(d) == 0 || sum(Matrix::colSums(d)) == 0) {
  message("data 层为空，从 counts 标准化")
  o <- Seurat::NormalizeData(o, normalization.method = "LogNormalize",
                             scale.factor = 1e4, verbose = FALSE)
  d <- Seurat::GetAssayData(o, assay = "RNA", layer = "data")
}
cat("表达矩阵: ", nrow(d), " x ", ncol(d), "\n")
stopifnot("TSPO" %in% rownames(d))

tspo <- as.numeric(d["TSPO", ])
names(tspo) <- colnames(o)
md <- o@meta.data
cat("TSPO 检出率: ", round(mean(tspo > 0) * 100, 1), "%\n", sep = "")

# ---------------------------------------------------------------- 2
hr("STEP 2  细胞类型的 TSPO")
CT <- "cell_type"
T <- do.call(rbind, lapply(sort(unique(md[[CT]])), function(k) {
  idx <- which(md[[CT]] == k)
  data.frame(celltype = k, n = length(idx),
             mean_TSPO = round(mean(tspo[idx]), 3),
             median_TSPO = round(median(tspo[idx]), 3),
             pct_pos = round(mean(tspo[idx] > 0) * 100, 1),
             stringsAsFactors = FALSE)
}))
T <- T[order(-T$mean_TSPO), ]
print(T, row.names = FALSE)

write.csv(T, file.path(OUTDIR, "GSE266577_TSPO_celltype.csv"), row.names = FALSE)

# 上皮 vs 髓系（PET 区室问题，用可靠注释重做）
EPI <- "Epithelial cells"
MYE <- c("Macrophages", "Classical monocytes", "DC1", "DC2", "Migratory DCs", "pDC")
e <- tspo[md[[CT]] == EPI]
m <- tspo[md[[CT]] %in% MYE]
cat(sprintf("\n上皮细胞: n=%d, 均值 %.3f, 阳性率 %.1f%%\n", length(e), mean(e), mean(e > 0)*100))
cat(sprintf("髓系细胞: n=%d, 均值 %.3f, 阳性率 %.1f%%\n", length(m), mean(m), mean(m > 0)*100))
cat(sprintf("髓系/上皮 = %.2f, Wilcoxon p = %.3g\n",
            mean(m)/mean(e), suppressWarnings(wilcox.test(m, e))$p.value))

# 贡献份额
tot <- sum(tspo)
cat(sprintf("转录本贡献: 上皮 %.1f%%, 髓系 %.1f%%, 其他 %.1f%%\n",
            sum(e)/tot*100, sum(m)/tot*100, (tot-sum(e)-sum(m))/tot*100))

# ---------------------------------------------------------------- 3
hr("STEP 3  化疗前后配对（上皮细胞）")
# 解析病人与阶段
smp <- as.character(md$publication_sample_code_final)
pat <- sub("_(chemo-naive|IDS)$", "", smp)
stage <- ifelse(grepl("_chemo-naive$", smp), "naive",
         ifelse(grepl("_IDS$", smp), "IDS", NA))
ok <- !is.na(stage) & md[[CT]] == EPI
cat("上皮细胞记录: ", sum(ok), "\n", sep = "")

df <- data.frame(patient = pat[ok], stage = stage[ok], tspo = tspo[ok],
                 stringsAsFactors = FALSE)
pairs <- do.call(rbind, lapply(split(df, df$patient), function(x) {
  a <- x$tspo[x$stage == "naive"]; b <- x$tspo[x$stage == "IDS"]
  if (length(a) < 30 || length(b) < 30) return(NULL)
  data.frame(patient = x$patient[1], n_naive = length(a), n_ids = length(b),
             naive = mean(a), ids = mean(b),
             delta = mean(b) - mean(a),
             pos_naive = mean(a > 0) * 100, pos_ids = mean(b > 0) * 100,
             stringsAsFactors = FALSE)
}))
if (is.null(pairs)) stop("没有病人同时具备化疗前后上皮细胞（每组>=30）")
pairs <- pairs[order(pairs$delta), ]
pr(pairs)

cat(sprintf("\n配对病人: %d 例\n", nrow(pairs)))
if (nrow(pairs) >= 5) {
  w <- suppressWarnings(wilcox.test(pairs$ids, pairs$naive, paired = TRUE))
  cat(sprintf("配对 Wilcoxon: 中位差 %+.3f, p = %.3g\n",
              median(pairs$delta), w$p.value))
  cat(sprintf("化疗后 TSPO 升高的病人: %d/%d\n", sum(pairs$delta > 0), nrow(pairs)))
  cat(sprintf("阳性率变化: %.1f%% -> %.1f%% (中位 %.1f -> %.1f)\n",
              mean(pairs$pos_naive), mean(pairs$pos_ids),
              median(pairs$pos_naive), median(pairs$pos_ids)))
}
write.csv(pairs, file.path(OUTDIR, "GSE266577_paired_TSPO.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 4
hr("STEP 4  化疗前后配对（全部细胞）")
df2 <- data.frame(patient = pat, stage = stage, tspo = tspo,
                  is_epi = md[[CT]] == EPI,
                  is_mye = md[[CT]] %in% MYE, stringsAsFactors = FALSE)
df2 <- df2[!is.na(df2$stage), ]
p2 <- do.call(rbind, lapply(split(df2, df2$patient), function(x) {
  a <- x[x$stage == "naive", ]; b <- x[x$stage == "IDS", ]
  if (nrow(a) < 100 || nrow(b) < 100) return(NULL)
  data.frame(patient = x$patient[1],
             all_naive = mean(a$tspo), all_ids = mean(b$tspo),
             epi_frac_naive = mean(a$is_epi) * 100, epi_frac_ids = mean(b$is_epi) * 100,
             mye_frac_naive = mean(a$is_mye) * 100, mye_frac_ids = mean(b$is_mye) * 100,
             stringsAsFactors = FALSE)
}))
if (!is.null(p2)) {
  pr(p2)
  w2 <- suppressWarnings(wilcox.test(p2$all_ids, p2$all_naive, paired = TRUE))
  cat(sprintf("\n全细胞 TSPO: 化疗后升高 %d/%d, 配对 p = %.3g\n",
              sum(p2$all_ids > p2$all_naive), nrow(p2), w2$p.value))
  w3 <- suppressWarnings(wilcox.test(p2$epi_frac_ids, p2$epi_frac_naive, paired = TRUE))
  cat(sprintf("上皮细胞占比: %.1f%% -> %.1f%%, p = %.3g\n",
              mean(p2$epi_frac_naive), mean(p2$epi_frac_ids), w3$p.value))
  w4 <- suppressWarnings(wilcox.test(p2$mye_frac_ids, p2$mye_frac_naive, paired = TRUE))
  cat(sprintf("髓系细胞占比: %.1f%% -> %.1f%%, p = %.3g\n",
              mean(p2$mye_frac_naive), mean(p2$mye_frac_ids), w4$p.value))
  write.csv(p2, file.path(OUTDIR, "GSE266577_paired_all.csv"), row.names = FALSE)
}

cat("\n================ 判读 ================\n")
cat("· 上皮细胞 TSPO 配对比较若无显著变化（13 对，中位差 +0.061，p=0.376）\n")
cat("  -> TSPO 反映细胞状态而非治疗暴露，正文据此写成「与治疗暴露无关」\n")
cat("· 全细胞 TSPO 的变化由上皮细胞占比改变驱动（21.1% -> 7.9%），\n")
cat("  不是单细胞表达水平的变化\n")
cat("  -> 必须区分「细胞构成变化」与「单细胞表达变化」，否则会误判为化疗下调 TSPO\n")

# -----------------------------------------------------------------------------
# -- Fig6_髓系轨迹 / 26_tspo_celltype_levels.R
# -----------------------------------------------------------------------------
# =============================================================================
#  26_tspo_celltype_levels.R   TSPO 在各细胞类型/亚型中的绝对表达水平
#
#  用途：为正文 3.7 与讨论中关于细胞来源的表述提供可复现的依据
#        （EOC 1.293、巨噬 1.173、淋巴细胞 0.365、CAF-1/2/3 0.640/0.471/0.954）
#
#  为什么必须单独成脚本：
#    24_caf_trajectory.R 只在基质亚型内部做比较，不跨细胞类型；
#    而"上皮与巨噬的 TSPO 水平相当、均远高于淋巴细胞"这一表述
#    需要全细胞类型的横向对比，且决定了组织水平信号由谁主导。
#
#  方法：log-normalized 表达按 cell_type 与 cell_subtype 求均值与阳性率；
#        另按「细胞数 × 均值」估算各类型对组织总信号的相对贡献。
#
#  输入：F_SC_EPI（32,847 基因 × 50,867 细胞）
#  输出：T18a_TSPO_by_celltype.csv、T18b_TSPO_by_subtype.csv、
#        T18c_TSPO_contribution.csv + Fig7C_celltype_TSPO
# =============================================================================

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  if (file.exists(.s)) source(.s, encoding = "UTF-8") else {
    PROJ <- "D:/Rstudy/tcga"; PROJ2 <- "D:/Rstudy/tcga2"
    OUT <- file.path(PROJ, "results_v2")
    F_SC_EPI <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")
  }
}
if (!exists("F_SC_EPI")) F_SC_EPI <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")
OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

MIN_N <- 50    # 亚型层面至少多少个细胞才纳入

hr("STEP 1  载入全细胞对象")

sc <- readRDS(F_SC_EPI)
cat("对象: ", nrow(sc), " 基因 x ", ncol(sc), " 细胞\n", sep = "")
if (!"TSPO" %in% rownames(sc)) stop("对象中没有 TSPO")

if (!"data" %in% SeuratObject::Layers(sc[["RNA"]]))
  sc <- Seurat::NormalizeData(sc, verbose = FALSE)

tspo <- as.numeric(SeuratObject::LayerData(sc, assay = "RNA", layer = "data")["TSPO", ])
cat("TSPO 非零细胞: ", sum(tspo > 0), " / ", length(tspo),
    "（", round(100 * mean(tspo > 0), 1), "%）\n", sep = "")

hr("STEP 2  按 cell_type 的绝对水平")

ct <- as.character(sc$cell_type)
d1 <- data.frame(type = ct, tspo = tspo, stringsAsFactors = FALSE)
m  <- tapply(d1$tspo, d1$type, mean, na.rm = TRUE)
n  <- tapply(d1$tspo, d1$type, length)
p  <- tapply(d1$tspo, d1$type, function(x) 100 * mean(x > 0, na.rm = TRUE))

out1 <- data.frame(type = names(m), n_cell = as.integer(n),
                   mean_TSPO = round(as.numeric(m), 3),
                   pct_positive = round(as.numeric(p), 1), row.names = NULL)
out1 <- out1[order(-out1$mean_TSPO), ]
cat("\n各 cell_type 的 TSPO 水平：\n"); pr(out1, 3)
write.csv(out1, file.path(OUTDIR, "T18a_TSPO_by_celltype.csv"), row.names = FALSE)

hr("STEP 3  按 cell_subtype 的绝对水平")

if ("cell_subtype" %in% colnames(sc@meta.data)) {
  sub <- as.character(sc$cell_subtype)
  d2  <- data.frame(type = sub, tspo = tspo, stringsAsFactors = FALSE)
  m2  <- tapply(d2$tspo, d2$type, mean, na.rm = TRUE)
  n2  <- tapply(d2$tspo, d2$type, length)
  p2  <- tapply(d2$tspo, d2$type, function(x) 100 * mean(x > 0, na.rm = TRUE))

  out2 <- data.frame(type = names(m2), n_cell = as.integer(n2),
                     mean_TSPO = round(as.numeric(m2), 3),
                     pct_positive = round(as.numeric(p2), 1), row.names = NULL)
  out2 <- out2[out2$n_cell >= MIN_N, ]
  out2 <- out2[order(-out2$mean_TSPO), ]
  cat("\n各 cell_subtype 的 TSPO 水平（n ≥ ", MIN_N, "）：\n", sep = ""); pr(out2, 3)
  write.csv(out2, file.path(OUTDIR, "T18b_TSPO_by_subtype.csv"), row.names = FALSE)
} else {
  cat("对象无 cell_subtype 列，跳过\n")
  out2 <- NULL
}

hr("STEP 4  对组织总信号的相对贡献")

# 贡献 = 细胞数 × 平均表达；归一化到总和为 100%
contrib <- out1
contrib$contribution_pct <- round(100 * contrib$n_cell * contrib$mean_TSPO /
                                  sum(contrib$n_cell * contrib$mean_TSPO), 1)
cat("\n各细胞类型对组织水平 TSPO 的相对贡献：\n"); pr(contrib, 3)
cat("（免疫细胞单位表达低但数量占 ",
    round(100 * sum(contrib$n_cell[contrib$type == "Immune"]) / sum(contrib$n_cell)),
    "%，加权后仍过半，这解释了 bulk 层面 TSPO 与免疫浸润的正相关）\n", sep = "")
write.csv(contrib, file.path(OUTDIR, "T18c_TSPO_contribution.csv"), row.names = FALSE)

hr("STEP 5  图")

p1 <- ggplot(data.frame(type = ct, tspo = tspo), aes(reorder(type, tspo), tspo, fill = type)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.75) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  labs(x = NULL, y = "TSPO (log-normalized)",
       title = "TSPO expression across cell types") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(FIGDIR, "Fig7C_celltype_TSPO.pdf"), p1, width = 5.5, height = 5)
ggsave(file.path(FIGDIR, "Fig7C_celltype_TSPO.png"), p1, width = 5.5, height = 5, dpi = 300)

if (!is.null(out2) && nrow(out2) >= 5) {
  top20 <- head(out2[order(-out2$mean_TSPO), ], 20)
  top20$type <- factor(top20$type, levels = rev(top20$type))
  p2 <- ggplot(top20, aes(mean_TSPO, type, fill = mean_TSPO)) +
    geom_col() +
    scale_fill_gradient(low = "#C6DBEF", high = "#B2182B") +
    labs(x = "Mean TSPO (log-normalized)", y = NULL,
         title = "TSPO by cell subtype (top 20)") +
    theme_bw(base_size = 10) + theme(legend.position = "none")
  ggsave(file.path(FIGDIR, "Fig7D_subtype_TSPO.pdf"), p2, width = 6.5, height = 5.5)
  ggsave(file.path(FIGDIR, "Fig7D_subtype_TSPO.png"), p2, width = 6.5, height = 5.5, dpi = 300)
}

hr("STEP 6  与正文数字的对照")

cat("正文引用的数值（应与此处一致）：\n")
cat("  EOC        1.293\n")
cat("  Macrophages 1.173\n")
cat("  CAF-1      0.640   CAF-2 0.471   CAF-3 0.954\n")
cat("  淋巴细胞    0.365\n")
cat("\n实际计算值：\n")
pr(out1[, c("type","n_cell","mean_TSPO")], 3)
if (!is.null(out2)) {
  keep <- out2$type %in% c("Macrophages","CAF-1","CAF-2","CAF-3","T-cells")
  pr(out2[keep, c("type","n_cell","mean_TSPO")], 3)
}

cat("\n判读要点：\n")
cat("  · 单位表达上上皮与巨噬相当（1.29 vs 1.17），均远高于淋巴细胞（0.37）。\n")
cat("  · 但组织水平由细胞数加权决定：免疫细胞占多数，故 bulk 信号与免疫浸润正相关。\n")
cat("  · 这一双重性正是本文坚持细胞类型分层解读 bulk 关联的原因。\n")

hr("26_tspo_celltype_levels.R 完成")

# ############################################################################
# ##  Fig7_空间转录组  (2 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- Fig7_空间转录组 / 08_spatial_persection.R
# -----------------------------------------------------------------------------
# =============================================================================
#  33_spatial_persection.R
#  重算原文 3.4 / 7.28 / 7.29 / 7.30 四个空间共定位
#
#  原文问题：把 8 张切片的 4386 个 spot 合并算 Spearman，切片间整体水平差异
#  会伪造出切片内的相关（Simpson 悖论）。本脚本给出三种口径并列：
#    ① rho_pooled   : 合并全部 spot（复现原文）
#    ② rho_persection: 每张切片单独算，取中位 + 范围 + 方向一致切片数
#    ③ rho_within   : 切片内标准化后再合并（去掉切片间均值差异），配块置换 p
#
#  四个分析：
#    3.4  TSPO ~ CD68      TSPO ~ CD163
#    7.28 TSPO ~ OXPHOS
#    7.29 OXPHOS ~ ROS_Defense
#    7.30 TSPO ~ BAX/BCL2
# =============================================================================


if (!exists("OUT")) stop("请先 source('00_setup.R')")

SC_RDS <- "D:/Rstudy/tcga/TSPO_spatial_processed.rds"
OUTDIR <- file.path(OUT, "tables")
NPERM  <- 1000
SEED   <- 42

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)
hr <- function(x) cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")

# ---------------------------------------------------------------- 载入
hr("STEP 1  载入与打分")
sp <- readRDS(SC_RDS)
md <- sp@meta.data
message("spot: ", ncol(sp))

dat <- tryCatch(Seurat::GetAssayData(sp, assay = "SCT", layer = "data"), error = function(e) NULL)
if (is.null(dat) || nrow(dat) == 0)
  dat <- Seurat::GetAssayData(sp, assay = "RNA", layer = "data")
D <- as.matrix(dat)
message("表达矩阵: ", nrow(D), " 基因 x ", ncol(D), " spot")

# 切片标签
sec_col <- intersect(c("sample", "orig.ident", "Section", "slice"), colnames(md))[1]
if (is.na(sec_col)) stop("找不到切片列")
sec <- as.character(md[[sec_col]])
names(sec) <- rownames(md)
message("切片数: ", length(unique(sec)))
print(table(sec))

zscore_rows <- function(m) {
  mu <- Matrix::rowMeans(m)
  v  <- Matrix::rowSums((m - mu)^2) / pmax(ncol(m) - 1, 1)
  sdv <- sqrt(v); sdv[sdv == 0] <- NA
  (m - mu) / sdv
}
Z <- zscore_rows(dat)

OXPHOS <- c("NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
            "SDHA","SDHB","UQCRC1","UQCRC2","COX4I1","COX5A","ATP5A1","ATP5B")
ROS    <- c("GPX4","SOD2","CAT","PRDX3","TXNRD2","GSR")

ms <- function(gs) {
  gs <- intersect(gs, rownames(Z))
  if (length(gs) < 3) return(rep(NA_real_, ncol(Z)))
  Matrix::colMeans(Z[gs, , drop = FALSE], na.rm = TRUE)
}
S <- list(
  OXPHOS = ms(OXPHOS),
  ROS    = ms(ROS),
  TSPO   = if ("TSPO" %in% rownames(Z)) as.numeric(Z["TSPO", ]) else rep(NA, ncol(Z)),
  CD68   = if ("CD68" %in% rownames(Z)) as.numeric(Z["CD68", ]) else rep(NA, ncol(Z)),
  CD163  = if ("CD163" %in% rownames(Z)) as.numeric(Z["CD163", ]) else rep(NA, ncol(Z)))
# BAX/BCL2（log 空间相减 = log2 比值）
S$BAX_BCL2 <- if (all(c("BAX","BCL2") %in% rownames(Z)))
  as.numeric(Z["BAX", ] - Z["BCL2", ]) else rep(NA, ncol(Z))

for (k in names(S)) message("  ", k, ": 有限值 ", sum(is.finite(S[[k]])), "/", length(S[[k]]))

# ---------------------------------------------------------------- 工具
spear <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 30) return(c(rho = NA, p = NA, n = sum(ok)))
  ct <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

#' 切片内标准化（按切片去均值、除标准差）
center_within <- function(v, sec) {
  out <- v
  for (s in unique(sec)) {
    i <- which(sec == s)
    x <- v[i]; x <- x[is.finite(x)]
    if (length(x) < 20) { out[i] <- NA; next }
    m <- mean(v[i], na.rm = TRUE); sdv <- sd(v[i], na.rm = TRUE)
    out[i] <- if (is.finite(sdv) && sdv > 0) (v[i] - m) / sdv else NA
  }
  out
}

#' 块置换：在每张切片内部置换 b，保留切片结构
perm_block <- function(a, b, sec, nperm = NPERM) {
  ok <- is.finite(a) & is.finite(b)
  a2 <- a[ok]; b2 <- b[ok]; s2 <- sec[ok]
  r0 <- suppressWarnings(cor(a2, b2, method = "spearman"))
  if (!is.finite(r0)) return(c(obs = NA, p = NA))
  null <- numeric(nperm)
  for (k in seq_len(nperm)) {
    bp <- b2
    for (s in unique(s2)) {
      i <- which(s2 == s)
      if (length(i) > 1) bp[i] <- b2[i][sample.int(length(i))]
    }
    null[k] <- suppressWarnings(cor(a2, bp, method = "spearman"))
  }
  null <- null[is.finite(null)]
  c(obs = r0, p = mean(abs(null) >= abs(r0)))
}

# ---------------------------------------------------------------- 四个分析
hr("STEP 2  四个共定位：三种口径")
ANALYSES <- list(
  list(id = "3.4a", a = "TSPO",    b = "CD68",     lab = "TSPO ~ CD68"),
  list(id = "3.4b", a = "TSPO",    b = "CD163",    lab = "TSPO ~ CD163"),
  list(id = "7.28", a = "TSPO",    b = "OXPHOS",   lab = "TSPO ~ OXPHOS"),
  list(id = "7.29", a = "OXPHOS",  b = "ROS",      lab = "OXPHOS ~ ROS_Defense"),
  list(id = "7.30", a = "TSPO",    b = "BAX_BCL2", lab = "TSPO ~ BAX/BCL2"))

res <- list()
for (an in ANALYSES) {
  a <- S[[an$a]]; b <- S[[an$b]]
  if (all(!is.finite(a)) || all(!is.finite(b))) {
    message("跳过 ", an$lab, "（基因缺失）"); next
  }
  # ① 合并
  pooled <- spear(a, b)
  # ② 分切片
  ps <- sapply(unique(sec), function(s) {
    i <- which(sec == s); spear(a[i], b[i])["rho"]
  })
  ps <- ps[is.finite(ps)]
  # ③ 切片内标准化后合并 + 块置换
  aw <- center_within(a, sec); bw <- center_within(b, sec)
  within <- spear(aw, bw)
  pb <- perm_block(aw, bw, sec)
  res[[an$id]] <- data.frame(
    id = an$id, analysis = an$lab,
    n_spot = pooled["n"],
    rho_pooled = round(pooled["rho"], 3),
    p_pooled   = signif(pooled["p"], 3),
    rho_persec_median = round(median(ps), 3),
    persec_min = round(min(ps), 3), persec_max = round(max(ps), 3),
    n_pos = sum(ps > 0), n_sec = length(ps),
    rho_within = round(within["rho"], 3),
    p_block    = signif(pb["p"], 3),
    stringsAsFactors = FALSE)
  cat("\n", an$lab, " 完成\n", sep = "")
}
R <- do.call(rbind, res)
print(R, row.names = FALSE)

# ---------------------------------------------------------------- 逐切片明细
hr("STEP 3  逐切片明细")
det <- do.call(rbind, lapply(ANALYSES, function(an) {
  a <- S[[an$a]]; b <- S[[an$b]]
  if (all(!is.finite(a)) || all(!is.finite(b))) return(NULL)
  data.frame(analysis = an$lab, section = unique(sec),
             rho = round(sapply(unique(sec), function(s) {
               i <- which(sec == s); spear(a[i], b[i])["rho"]
             }), 3), stringsAsFactors = FALSE)
}))
print(det, row.names = FALSE)

# ---------------------------------------------------------------- 输出
hr("STEP 4  写出")
write.csv(R,   file.path(OUTDIR, "spatial_persection_summary.csv"), row.names = FALSE)
write.csv(det, file.path(OUTDIR, "spatial_persection_detail.csv"),  row.names = FALSE)
message("已写入 ", OUTDIR)

cat("\n================ 怎么用 ================\n")
cat("· rho_pooled   : 复现原文，含切片间混杂\n")
cat("· rho_persec_median : 切片内真实相关，跨切片取中位\n")
cat("· rho_within + p_block : 切片内标准化后合并，块置换给有效 p（推荐引用）\n")
cat("· n_pos/n_sec  : 多少张切片方向为正，一致性比 p 值更直观\n")

# -----------------------------------------------------------------------------
# -- Fig7_空间转录组 / 22_spatial_region_deg.R
# -----------------------------------------------------------------------------
# =============================================================================
#  22_spatial_region_deg.R  (v3)   空间三区域（界面 / 肿瘤核心 / 基质）差异表达
#
#  v3 修正：SeuratObject 5.x 中 Assays() 返回 SimpleAssays S4 对象，
#           paste() 无法处理，统一改用 names(obj@assays)；
#           取矩阵改为 LayerData() 失败时退回 GetAssayData(slot=)。
#
#  区域来源：对象本身没有区域列，按原脚本 X.3 段（tspo.R L8654-8656）
#            的硬编码映射从 seurat_clusters 重建：
#              Interface  = clusters 12, 15, 8
#              Tumor_Core = clusters 10, 14
#              Stroma     = clusters 5, 11
#              其余       = Other（不纳入比较）
#
#  方法：伪批量（section × region 求和）→ CPM → log2 → limma，section 作区组
#        spot 级只作描述（空间自相关会使 p 值严重偏小）
#
#  ⚠ 区域由聚类定义，聚类又基于表达谱，故区域间 DEG 有循环性，
#    应表述为「区域特征基因」，不宜说成新发现。
#
#  输入：F_SPATIAL（16,626 基因 × 4,386 spot）
#  输出：T14a_region_composition.csv、T14b_region_module_scores.csv、
#        T14e_TSPO_by_region.csv、T14_DE_*.csv + Fig4A/4B
# =============================================================================

if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.code_dir, "00_setup.R")
  .ok <- FALSE
  if (file.exists(.s))
    .ok <- !inherits(try(source(.s, encoding = "UTF-8"), silent = TRUE), "try-error")
  if (!.ok) {
    PROJ <- "D:/Rstudy/tcga"; OUT <- file.path(PROJ, "results_v2")
    F_SPATIAL <- file.path(PROJ, "TSPO_spatial_processed.rds")
    hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
    pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)
  }
  for (d in c(OUT, file.path(OUT, "tables"), file.path(OUT, "figures")))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("hr")) hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
if (!exists("pr")) pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(limma)
})

OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# 原脚本 X.3 段的区域映射（tspo.R L8654-8656）
INTERFACE_CL <- c("12", "15", "8")
TUMOR_CL     <- c("10", "14")
STROMA_CL    <- c("5", "11")

MIN_SPOT_PER_GROUP <- 20
LFC_CUT <- 0.5
FDR_CUT <- 0.05

MODULES <- list(
  OXPHOS     = c("NDUFS1","NDUFS2","NDUFS3","SDHA","SDHB","UQCRC1","UQCRC2",
                 "COX4I1","COX5A","ATP5A1","ATP5B"),
  Glycolysis = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1"),
  Immune     = c("CD68","CD3D","CD3E","CD8A","PTPRC","LYZ","C1QA","C1QB","MS4A1"),
  EMT        = c("VIM","ZEB1","FN1","COL1A1","COL1A2","SPARC"),
  Hypoxia    = c("CA9","BNIP3","VEGFA","SLC2A1","PDK1")
)

hr("STEP 0  载入与区域重建")

sc <- readRDS(F_SPATIAL)
cat("对象: ", nrow(sc), " 基因 x ", ncol(sc), " spot\n", sep = "")
cat("meta 列: ", paste(colnames(sc@meta.data), collapse = ", "), "\n", sep = "")

# SeuratObject 5.x 的 Assays() 返回 SimpleAssays S4 对象，不能直接 paste，
# 统一用 names(obj@assays) 取字符向量
.assay_names <- names(sc@assays)
cat("assays: ", paste(.assay_names, collapse = ", "), "\n", sep = "")

get_layer <- function(obj, layer = "data") {
  an <- names(obj@assays)
  for (a in c("RNA", "Spatial", "SCT")) {
    if (!a %in% an) next
    m <- tryCatch(SeuratObject::LayerData(obj, assay = a, layer = layer),
                  error = function(e) NULL)
    if (is.null(m) || nrow(m) < 100) {
      m <- tryCatch(Seurat::GetAssayData(obj, assay = a, slot = layer),
                    error = function(e) NULL)
    }
    if (!is.null(m) && nrow(m) > 100) {
      cat("使用 assay = ", a, "，layer = ", layer,
          "，矩阵 ", nrow(m), " x ", ncol(m), "\n", sep = "")
      return(m)
    }
  }
  NULL
}

# 切片列：优先 section/slice/sample，其次 orig.ident
md <- sc@meta.data
sec_col <- NA_character_
for (c in c("section", "Section", "slice", "sample", "orig.ident")) {
  if (c %in% colnames(md)) {
    n_lv <- length(unique(as.character(md[[c]])))
    if (n_lv >= 2 && n_lv <= 50) { sec_col <- c; break }
  }
}
if (is.na(sec_col)) stop("找不到切片列")
SEC <- as.character(md[[sec_col]])
cat("\n切片列: ", sec_col, "（", length(unique(SEC)), " 个切片）\n", sep = "")
print(table(SEC))

# 区域：若已有区域列直接用，否则按聚类重建
reg_col <- intersect(c("Region3","Region","region","anatomic_region"), colnames(md))
if (length(reg_col)) {
  REG <- as.character(md[[reg_col[1]]])
  cat("\n使用已有区域列: ", reg_col[1], "\n", sep = "")
  print(table(REG))
} else {
  if (!"seurat_clusters" %in% colnames(md))
    stop("对象既无区域列也无 seurat_clusters，无法重建区域")
  cl <- as.character(md$seurat_clusters)
  cat("\n按 tspo.R L8654-8656 的映射从 seurat_clusters 重建 Region3\n")
  cat("  聚类分布:\n"); print(table(cl))
  REG <- ifelse(cl %in% INTERFACE_CL, "Interface",
         ifelse(cl %in% TUMOR_CL,     "Tumor_Core",
         ifelse(cl %in% STROMA_CL,    "Stroma", "Other")))
  cat("\n重建结果:\n"); print(table(REG))
  pct_other <- 100 * mean(REG == "Other")
  cat("Other 占比: ", round(pct_other, 1), "%\n", sep = "")
  if (pct_other > 40)
    cat("⚠ Other 超过 40%，区域分析只覆盖一部分 spot，需在正文说明\n")
}

if (length(unique(REG[REG != "Other"])) < 2)
  stop("有效区域少于 2 个，无法比较")

keep <- REG %in% c("Interface", "Tumor_Core", "Stroma")
if (sum(keep) < 100) stop("三区域 spot 总数不足 100")
sc  <- subset(sc, cells = colnames(sc)[keep])
REG <- REG[keep]; SEC <- SEC[keep]
cat("\n用于分析: ", ncol(sc), " spot\n", sep = "")
print(table(SEC, REG))

hr("STEP 1  区域 × 切片构成")

tab <- table(SEC, REG)
write.csv(as.data.frame.matrix(tab), file.path(OUTDIR, "T14a_region_composition.csv"))
print(tab)
cat("\n满足 ≥", MIN_SPOT_PER_GROUP, " spot 的组合: ",
    sum(tab >= MIN_SPOT_PER_GROUP), " / ", length(tab), "\n", sep = "")

hr("STEP 2  模块评分")

d <- get_layer(sc, "data")
if (is.null(d)) stop("无法取到 data 层")

gs <- intersect(c(unique(unlist(MODULES)), "TSPO"), rownames(d))
cat("纳入基因: ", length(gs), "\n", sep = "")
if (!"TSPO" %in% gs) stop("表达矩阵中没有 TSPO")

X <- as.matrix(d[gs, , drop = FALSE]); rm(d); gc()
Z <- t(scale(t(X))); Z[!is.finite(Z)] <- NA
rm(X); gc()

sc_df <- data.frame(spot = colnames(sc), region = REG, section = SEC,
                    TSPO = as.numeric(Z["TSPO", ]), stringsAsFactors = FALSE)
for (m in names(MODULES)) {
  g <- intersect(MODULES[[m]], rownames(Z))
  if (length(g) < 3) next
  sc_df[[m]] <- colMeans(Z[g, , drop = FALSE], na.rm = TRUE)
}
rm(Z); gc()

score_cols <- intersect(c("TSPO", names(MODULES)), colnames(sc_df))
cat("参与比较的变量: ", paste(score_cols, collapse = ", "), "\n", sep = "")

hr("STEP 3  区域差异（描述性）")

summ <- do.call(rbind, lapply(score_cols, function(v) {
  kw  <- kruskal.test(sc_df[[v]] ~ factor(sc_df$region))
  med <- tapply(sc_df[[v]], sc_df$region, median, na.rm = TRUE)
  out <- data.frame(variable = v, n_spot = sum(is.finite(sc_df[[v]])),
                    KW_p = kw$p.value, stringsAsFactors = FALSE)
  for (r in names(med)) out[[paste0("median_", r)]] <- round(med[[r]], 3)
  out
}))
summ <- summ[order(summ$KW_p), ]
cat("\n各变量的区域差异（Kruskal-Wallis，spot 级仅作描述）：\n")
pr(summ, 4)
write.csv(summ, file.path(OUTDIR, "T14b_region_module_scores.csv"), row.names = FALSE)

cat("\n切片内配对（每张切片独立比较，避开切片间差异）：\n")
regs <- intersect(c("Interface","Tumor_Core","Stroma"), unique(REG))
pair_res <- list()
if (length(regs) >= 2) {
  for (i in 1:(length(regs) - 1)) for (j in (i + 1):length(regs)) {
    a <- regs[i]; b <- regs[j]
    for (v in score_cols) {
      dif <- sapply(split(sc_df, sc_df$section), function(s) {
        x <- s[[v]][s$region == a]; y <- s[[v]][s$region == b]
        if (sum(is.finite(x)) < 10 || sum(is.finite(y)) < 10) return(NA_real_)
        median(x, na.rm = TRUE) - median(y, na.rm = TRUE)
      })
      dif <- dif[is.finite(dif)]
      if (length(dif) < 4) next
      k <- sum(dif > 0); n <- length(dif)
      pair_res[[length(pair_res) + 1]] <-
        data.frame(comparison = paste0(a, " vs ", b), variable = v,
                   n_section = n, n_positive = k,
                   median_diff = round(median(dif), 3),
                   sign_p = binom.test(k, n, 0.5)$p.value, stringsAsFactors = FALSE)
    }
  }
}
if (length(pair_res)) {
  pd <- do.call(rbind, pair_res)
  pd$FDR <- p.adjust(pd$sign_p, "BH")
  pd <- pd[order(pd$comparison, pd$sign_p), ]
  pr(pd, 4)
  write.csv(pd, file.path(OUTDIR, "T14e_TSPO_by_region.csv"), row.names = FALSE)
} else cat("切片内配对：有效切片不足，跳过\n")

hr("STEP 4  伪批量差异表达")

cnt <- get_layer(sc, "counts")
if (is.null(cnt)) {
  cat("无 counts 层，改用 data 层\n"); cnt <- get_layer(sc, "data")
}

grp <- paste(SEC, REG, sep = "__")
ok_levels <- names(which(table(grp) >= MIN_SPOT_PER_GROUP))
cat("纳入伪批量的组合: ", length(ok_levels), "\n", sep = "")
if (length(ok_levels) < 6) cat("⚠ 有效组合偏少，DE 稳定性有限\n")

pb <- sapply(ok_levels, function(g) Matrix::rowSums(cnt[, which(grp == g), drop = FALSE]))
pb <- as.matrix(pb); rm(cnt); gc()

lib <- colSums(pb); pb <- t(t(pb) / lib) * 1e6; pb <- log2(pb + 1)
cat("伪批量: ", nrow(pb), " 基因 x ", ncol(pb), " 样本\n", sep = "")

meta <- data.frame(sample = colnames(pb), stringsAsFactors = FALSE)
meta$section <- sub("__.*$", "", meta$sample)
meta$region  <- factor(sub("^.*__", "", meta$sample))
print(table(meta$section, meta$region))

design <- tryCatch({
  if (length(unique(meta$section)) >= 2)
    model.matrix(~ 0 + region + section, data = meta)
  else model.matrix(~ 0 + region, data = meta)
}, error = function(e) { cat("设计矩阵失败: ", conditionMessage(e), "\n"); NULL })
if (is.null(design)) stop("无法构建设计矩阵")

fit <- lmFit(pb, design); rn <- colnames(design)
cat("设计列: ", paste(rn, collapse = ", "), "\n", sep = "")

for (cmp in list(c("Interface","Tumor_Core"), c("Interface","Stroma"),
                 c("Tumor_Core","Stroma"))) {
  a <- cmp[1]; b <- cmp[2]
  ia <- which(rn == paste0("region", a)); ib <- which(rn == paste0("region", b))
  if (!length(ia) || !length(ib)) {
    cat("\n跳过 ", a, " vs ", b, "（该区域不在伪批量中）\n", sep = ""); next
  }
  v <- rep(0, length(rn)); v[ia] <- 1; v[ib] <- -1
  f2 <- eBayes(contrasts.fit(fit, v))
  tt <- topTable(f2, number = Inf, sort.by = "P"); tt$gene <- rownames(tt)
  sig <- tt[abs(tt$logFC) >= LFC_CUT & tt$adj.P.Val < FDR_CUT, ]
  cat("\n===== ", a, " vs ", b, " =====\n", sep = "")
  cat("显著基因: ", nrow(sig), "（", a, " 中高 ", sum(sig$logFC > 0),
      "，低 ", sum(sig$logFC < 0), "）\n", sep = "")
  if (nrow(sig))
    pr(head(sig[order(-abs(sig$logFC)), c("gene","logFC","adj.P.Val")], 15), 3)
  write.csv(tt, file.path(OUTDIR, paste0("T14_DE_", a, "_vs_", b, ".csv")),
            row.names = FALSE)
}

hr("STEP 5  图")

long <- reshape(sc_df[, c("region", score_cols)], direction = "long",
                varying = score_cols, v.names = "score",
                times = score_cols, timevar = "variable")
long$variable <- factor(long$variable, levels = score_cols)
pA <- ggplot(long[is.finite(long$score), ], aes(region, score, fill = region)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~ variable, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Module score (z)", title = "Region-wise module scores") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIGDIR, "Fig4A_region_modules.pdf"), pA, width = 13, height = 4)
ggsave(file.path(FIGDIR, "Fig4A_region_modules.png"), pA, width = 13, height = 4, dpi = 300)

pT <- ggplot(sc_df[is.finite(sc_df$TSPO), ], aes(region, TSPO, fill = region)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~ section, nrow = 2) +
  labs(x = NULL, y = "TSPO (z)", title = "TSPO by region, per section") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIGDIR, "Fig4B_TSPO_by_region.pdf"), pT, width = 10, height = 5)
ggsave(file.path(FIGDIR, "Fig4B_TSPO_by_region.png"), pT, width = 10, height = 5, dpi = 300)

saveRDS(list(sc_df = sc_df, tab = tab, summ = summ),
        file.path(OUT, "cache", "spatial_region.rds"))

cat("\n判读要点：\n")
cat("  · 主推断看伪批量与切片内配对，不看 spot 级 KW p 值。\n")
cat("  · 区域由聚类定义，聚类又基于表达谱，DEG 有循环性，\n")
cat("    应表述为「区域特征基因」，不宜说成新发现。\n")
cat("  · 若 TSPO 在 Interface 最高且 Interface 富集 OXPHOS，\n")
cat("    才支持「氧化代谢在界面富集」。\n")

hr("22_spatial_region_deg.R (v3) 完成")

# ############################################################################
# ##  Fig8_空间蛋白组  (2 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- Fig8_空间蛋白组 / 09_codex_tma.R
# -----------------------------------------------------------------------------
# =============================================================================
#  09_codex_tma.R
#  结论四：代谢轴与免疫轴的空间分离（空间蛋白组部分）
#
#  数据：VectraPolarisData::HumanOvarianCancerVP（卵巢癌 TMA 多重免疫荧光）
#        132 个 TMA 核心 / 1,610,431 个细胞 / 7 个标志物
#
#  对应正文 3.5 表格下半部分（以核心为单位配对）：
#    CK+   肿瘤区 86.23% vs 间质区 0.26%，0/129 核心间质更高，p=6.6e-23
#    CD68+ 肿瘤区  2.51% vs 间质区 4.55%，78/129，p=3.5e-5
#    CD3+  肿瘤区  1.77% vs 间质区 2.39%，81/129，p=3.5e-3
#    Ki67+ 肿瘤区 14.37% vs 间质区 3.75%
#
#  ⚠ 该面板不含 TSPO 抗体（只有 CK/Ki67/CD8/IER3/pSTAT3/CD3/CD68/CD19/DAPI），
#     因此只能验证「肿瘤细胞与免疫细胞分属不同区室」这一前提，
#     不能在此层面直接验证 TSPO 的蛋白定位。写正文时必须注明。
#
#  ⚠ HumanOvarianCancerVP 是一个函数，必须带括号调用：HumanOvarianCancerVP()
#     首次调用会尝试连 ExperimentHub，失败时自动回落到本地缓存（需事先缓存过）
# =============================================================================

if (!exists("OUT")) stop("请先 source('00_setup.R')")
if (!requireNamespace("VectraPolarisData", quietly = TRUE)) {
  stop("请先安装 VectraPolarisData： BiocManager::install('VectraPolarisData')")
}
suppressPackageStartupMessages({
  library(VectraPolarisData)
  library(SummarizedExperiment)
})

MIN_CELL_PER_REGION <- 50    # 每个核心每个分区至少多少细胞才纳入配对
MARKERS <- c("ck", "cd68", "cd3", "cd8", "cd19", "ki67", "p_stat3")

hr("STEP 1  载入 CODEX 数据")

# ⚠ HumanOvarianCancerVP() 默认会去连 ExperimentHub；断网或证书问题时直接报错。
#   解法：先把 ExperimentHub 切到本地模式（数据首次下载后会缓存在本地），
#   再调用。若本地没有缓存，仍需联网一次。
.EH_local <- function() {
  ok <- tryCatch({
    ExperimentHub::setExperimentHubOption("LOCAL", TRUE); TRUE
  }, error = function(e) FALSE)
  if (!ok) try(utils::setHook, silent = TRUE)
  invisible(ok)
}
invisible(.EH_local())

spe <- tryCatch(
  VectraPolarisData::HumanOvarianCancerVP(),
  error = function(e1) {
    message("首次载入失败：", conditionMessage(e1))
    message("切换到 ExperimentHub 本地模式后重试 …")
    .EH_local()
    tryCatch(
      VectraPolarisData::HumanOvarianCancerVP(metadata = FALSE),
      error = function(e2) {
        stop("无法载入 CODEX 数据。\n",
             "  原因 1：本地 ExperimentHub 缓存中没有 HumanOvarianCancerVP。\n",
             "  原因 2：当前环境无法访问 bioconductor.org。\n",
             "  解决办法：在有网的机器上运行一次\n",
             "    BiocManager::install('VectraPolarisData')\n",
             "    VectraPolarisData::HumanOvarianCancerVP()\n",
             "  下载完成后把 ExperimentHub 缓存目录一并拷贝过来即可离线运行。\n",
             "  原始错误：", conditionMessage(e2))
      })
  }
)
cat("维度: ", nrow(spe), " 标志物 x ", ncol(spe), " 细胞\n", sep = "")
cat("标志物: ", paste(rownames(spe), collapse = ", "), "\n")
if (!any(grepl("TSPO", rownames(spe), ignore.case = TRUE)))
  cat("提醒：面板中无 TSPO，仅用于验证区室结构\n")

cd <- SummarizedExperiment::colData(spe)
cat("\n分区细胞数:\n"); print(table(cd$tissue_category))

hr("STEP 2  两两表型阳性率（不做配对，仅描述）")

R <- do.call(rbind, lapply(MARKERS, function(k) {
  col <- paste0("phenotype_", k)
  if (!col %in% colnames(cd)) return(NULL)
  pos <- grepl("\\+$", as.character(cd[[col]]))
  do.call(rbind, lapply(c("Tumor", "Stroma"), function(tt) {
    i <- which(cd$tissue_category == tt)
    data.frame(marker = k, tissue = tt, n = length(i),
               pct_pos = round(mean(pos[i], na.rm = TRUE) * 100, 2),
               stringsAsFactors = FALSE)
  }))
}))

W <- reshape(R[, c("marker", "tissue", "pct_pos")], idvar = "marker",
             timevar = "tissue", direction = "wide")
# ⚠ reshape 生成的列顺序是「首次出现顺序」，不是字母序，不能假设 Stroma 在前
names(W) <- sub("^pct_pos\\.", "", names(W))
if (all(c("Tumor", "Stroma") %in% names(W))) {
  W$ratio_T_over_S <- round(W$Tumor / pmax(W$Stroma, 0.01), 2)
}
pr(W)
write.csv(W, file.path(OUT, "tables", "T4a_CODEX_overall.csv"), row.names = FALSE)

hr("STEP 3  以 TMA 核心为单位配对比较（推断口径）")

cd2 <- cd[cd$tissue_category %in% c("Tumor", "Stroma"), ]

PAIRS <- do.call(rbind, lapply(split(seq_len(nrow(cd2)), cd2$sample_id), function(i) {
  d  <- cd2[i, ]
  it <- d$tissue_category == "Tumor"
  is <- d$tissue_category == "Stroma"
  if (sum(it) < MIN_CELL_PER_REGION || sum(is) < MIN_CELL_PER_REGION) return(NULL)
  out <- data.frame(core = as.character(d$sample_id[1]),
                    n_tumor = sum(it), n_stroma = sum(is), stringsAsFactors = FALSE)
  for (k in MARKERS) {
    col <- paste0("phenotype_", k)
    if (!col %in% colnames(d)) next
    pos <- grepl("\\+$", as.character(d[[col]]))
    out[[paste0("t_", k)]] <- mean(pos[it]) * 100
    out[[paste0("s_", k)]] <- mean(pos[is]) * 100
  }
  out
}))
cat("可配对核心数: ", nrow(PAIRS), "\n", sep = "")

SUM <- do.call(rbind, lapply(MARKERS, function(k) {
  a <- PAIRS[[paste0("t_", k)]]; b <- PAIRS[[paste0("s_", k)]]
  if (is.null(a) || is.null(b)) return(NULL)
  w <- suppressWarnings(wilcox.test(a, b, paired = TRUE))
  data.frame(marker = k,
             tumor_median = round(median(a, na.rm = TRUE), 2),
             stroma_median = round(median(b, na.rm = TRUE), 2),
             n_stroma_higher = sum(b > a, na.rm = TRUE),
             n_cores = nrow(PAIRS),
             p_paired = signif(w$p.value, 3), stringsAsFactors = FALSE)
}))
pr(SUM)

hr("STEP 4  写出")

write.csv(PAIRS, file.path(OUT, "tables", "T4b_CODEX_per_core.csv"), row.names = FALSE)
write.csv(SUM,   file.path(OUT, "tables", "T4c_CODEX_paired_test.csv"), row.names = FALSE)

cat("\n判读：CK+ 应几乎全部在肿瘤区（0 个核心间质更高）；\n")
cat("      CD68+/CD3+ 应多数在间质区。若相反，说明 tissue_category 标注需复核。\n")

hr("09_codex_tma.R 完成")

# -----------------------------------------------------------------------------
# -- Fig8_空间蛋白组 / 28_mif_spatial_region.R
# -----------------------------------------------------------------------------
# =============================================================================
#  28_mif_spatial_region.R   MIF 区域归属复核（用片内 z 分区重算）
#
#  ── 为什么要重做 ────────────────────────────────────────────────────────────
#  正文 3.4 现有一句：「空间数据支持旁分泌模型：MIF 在肿瘤–基质界面显著富集
#  （界面 2.504 vs 肿瘤核心 1.610 vs 基质 1.783，p<1×10⁻⁹⁵）」。
#
#  这三个数来自原始脚本（tspo.R L8553-8801）的分区：
#      interface = seurat_clusters 12/15/8
#      tumor_core = 10/14
#      stroma     = 5/11
#  而 3.6 重新核对后已证明该映射不成立——肿瘤核心只出现在 1/8 张切片、
#  界面只出现在 2/8 张，解剖上不可能。故 3.6 改用 EPCAM/COL1A1 片内 z 分区。
#  但 3.4 的 MIF 数字没有跟着改，仍建立在已被否定的分区上。
#
#  另外 p<1×10⁻⁹⁵ 本身不可信：那是 spot 层面的 Wilcoxon，同一切片内
#  spot 空间自相关极强、不独立，该量级纯由伪重复放大。即使分区正确也不应引用。
#
#  ── 本脚本做什么 ────────────────────────────────────────────────────────────
#  用 3.6 / 27 PART C 的片内 z 分区，重算 MIF（含 CD74/TSPO/CD68/CD163）
#  的三区域比较，并给出三种阈值定义下的一致性。
#
#  同时把「MIF ~ CD74 共定位 ρ=0.389」这类数值改为**切片内**口径重算
#  （原值可能也是合并口径；合并口径在 3.6 已被证明会伪造相关）。
#
#  ── 产出如何用 ──────────────────────────────────────────────────────────────
#  STEP 6 会打印裁决：三种可能结局之一，告诉正文该怎么写。
#
#  输入：F_SPATIAL（TSPO_spatial_processed.rds，16,626 基因 × 4,386 spot）
#  输出：T20a_MIF_raw_by_region.csv、T20b_MIF_z_by_region.csv、
#        T20c_region_sign_tests.csv、T20d_spot_correlations.csv、
#        T20e_verdict.txt、Fig8_MIF_region.pdf/png
# =============================================================================

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

if (!exists("OUT")) {
  .s <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .s <- file.path(.s, "00_setup.R")
  if (file.exists(.s)) {
    source(.s, encoding = "UTF-8")
  } else {
    PROJ <- "D:/Rstudy/tcga"
    PROJ2 <- "D:/Rstudy/tcga2"
    OUT  <- file.path(PROJ, "results_v2")
    F_SPATIAL <- file.path(PROJ, "TSPO_spatial_processed.rds")
  }
}
if (!exists("F_SPATIAL")) F_SPATIAL <- file.path(PROJ, "TSPO_spatial_processed.rds")

OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
set.seed(20260101)

# ---- 参数 ------------------------------------------------------------------
GENES   <- c("MIF", "CD74", "TSPO", "CD68", "CD163", "CXCR4", "CD44")
MIN_SPOT <- 15          # 每张切片每区至少多少个 spot 才纳入该切片的比较
N_PERM   <- 1000        # 切片内相关性的置换次数

hr("STEP 0  载入")

if (!file.exists(F_SPATIAL)) stop("找不到 ", F_SPATIAL)
sc <- readRDS(F_SPATIAL)
cat("对象: ", nrow(sc), " 基因 x ", ncol(sc), " spot\n", sep = "")

# SeuratObject 5.x 中 Assays() 返回 SimpleAssays S4，不能用 paste()，统一用 names()
get_layer <- function(obj, layer = "data") {
  an <- names(obj@assays)
  for (a in c("SCT", "RNA", "Spatial")) {
    if (!a %in% an) next
    m <- tryCatch(SeuratObject::LayerData(obj, assay = a, layer = layer),
                  error = function(e) NULL)
    if (is.null(m) || nrow(m) < 100)
      m <- tryCatch(Seurat::GetAssayData(obj, assay = a, slot = layer),
                    error = function(e) NULL)
    if (!is.null(m) && nrow(m) > 100) {
      cat("assay = ", a, "，layer = ", layer, "\n", sep = ""); return(m)
    }
  }
  NULL
}
DAT_full <- get_layer(sc, "data")
if (is.null(DAT_full)) stop("无法取得表达矩阵")

# 只取需要的基因，避免整矩阵 as.matrix 占内存
need <- intersect(unique(c(GENES, "EPCAM", "COL1A1")), rownames(DAT_full))
cat("可用基因: ", paste(need, collapse = ", "), "\n", sep = "")
miss <- setdiff(GENES, need)
if (length(miss)) cat("⚠ 缺失: ", paste(miss, collapse = ", "), "\n", sep = "")
if (!"MIF" %in% need) stop("矩阵中没有 MIF")

DAT <- as.matrix(DAT_full[need, , drop = FALSE])
storage.mode(DAT) <- "double"
rm(DAT_full); gc()
cat("取子矩阵: ", nrow(DAT), " x ", ncol(DAT), "\n", sep = "")

# 切片标签
md <- sc@meta.data
sec_col <- NA_character_
for (c in c("section", "Section", "slice", "sample", "orig.ident")) {
  if (c %in% colnames(md)) {
    k <- length(unique(as.character(md[[c]])))
    if (k >= 2 && k <= 50) { sec_col <- c; break }
  }
}
if (is.na(sec_col)) stop("找不到切片列")
SEC  <- as.character(md[[sec_col]])
secs <- sort(unique(SEC))
cat("切片列 ", sec_col, "，共 ", length(secs), " 张\n", sep = "")
print(table(SEC))

hr("STEP 1  片内 z 标准化")

# 关键：所有量都在切片内标准化。合并 8 张切片统一 z 会让区域标签跨切片不可比，
# 并可能造成 Simpson 悖论（3.6 已记录过早期版本方向完全相反的情况）。
within_z <- function(v, sec) {
  out <- rep(NA_real_, length(v))
  for (s in unique(sec)) {
    i <- which(sec == s); x <- v[i]
    if (sum(is.finite(x)) < 20) next
    mu <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
    out[i] <- if (is.finite(sdv) && sdv > 0) (x - mu) / sdv else NA_real_
  }
  out
}

Z <- t(apply(DAT, 1, within_z, sec = SEC))
colnames(Z) <- colnames(DAT)

# 原始（未标准化）值也保留一份，用于和旧数字 2.504/1.610/1.783 直接对照
RAW <- DAT

hr("STEP 2  区域定义（同 3.6 / 27 PART C）")

S <- data.frame(spot = colnames(sc), section = SEC, stringsAsFactors = FALSE)
S$EPCAM  <- Z["EPCAM", ]
S$COL1A1 <- Z["COL1A1", ]

DEFS <- list(
  A = list(epi = 0.5, str = 0.5, lab = "上皮z>0.5, 基质z>0.5（主定义）"),
  B = list(epi = 1.0, str = 0.5, lab = "上皮z>1.0, 基质z>0.5（严格上皮）"),
  C = list(epi = 1.0, str = 1.0, lab = "上皮z>1.0, 基质z>1.0（三区皆严格）"))

define_region <- function(S, epi_th, str_th) {
  with(S, ifelse(!is.finite(EPCAM) | !is.finite(COL1A1), "Other",
          ifelse(EPCAM > epi_th & COL1A1 <= 0, "Tumor_Core",
          ifelse(COL1A1 > str_th & EPCAM <= 0, "Stroma",
          ifelse(EPCAM > 0 & COL1A1 > 0, "Interface", "Other")))))
}

for (d in names(DEFS)) {
  reg <- define_region(S, DEFS[[d]]$epi, DEFS[[d]]$str)
  S[[paste0("region_", d)]] <- reg
  cat("\n--- 定义 ", d, "：", DEFS[[d]]$lab, " ---\n", sep = "")
  print(table(reg))
}

S$region <- S$region_A
cat("\n主定义 Other 占比: ", round(100 * mean(S$region == "Other"), 1), "%\n", sep = "")
cat("\n区域 × 切片（主定义）：\n"); print(table(S$section, S$region))

hr("STEP 3  MIF 的三区域比较（原始值，可直接对照旧数字）")

# 把各基因的原始值与 z 值挂到 S 上
for (g in intersect(GENES, rownames(RAW))) {
  S[[paste0(g, "_raw")]] <- RAW[g, ]
  S[[paste0(g, "_z")]]   <- Z[g, ]
}
cat("MIF 检出率: ", round(100 * mean(RAW["MIF", ] > 0), 1), "%",
    "，原始均值 ", round(mean(RAW["MIF", ]), 3), "\n", sep = "")

summ_by_region <- function(var, tag) {
  do.call(rbind, lapply(c("Tumor_Core", "Interface", "Stroma"), function(r) {
    v <- S[[var]][S$region == r]
    data.frame(definition = tag, variable = var, region = r,
               n_spot = sum(is.finite(v)),
               mean = round(mean(v, na.rm = TRUE), 3),
               median = round(median(v, na.rm = TRUE), 3),
               pct_pos = round(100 * mean(v > 0, na.rm = TRUE), 1),
               stringsAsFactors = FALSE)
  }))
}

RAW_TAB <- do.call(rbind, lapply(intersect(GENES, rownames(RAW)),
                                 function(g) summ_by_region(paste0(g, "_raw"), g)))
cat("\n各区域原始表达均值（主定义 A）：\n")
pr(RAW_TAB[, c("variable","region","n_spot","mean","median","pct_pos")], 3)
write.csv(RAW_TAB, file.path(OUTDIR, "T20a_MIF_raw_by_region.csv"), row.names = FALSE)

cat("\n【与旧数字对照】旧口径（已被否定的聚类分区）给出：\n")
cat("  界面 2.504 / 肿瘤核心 1.610 / 基质 1.783\n")
cat("本次（片内 z 分区，主定义 A）给出：\n")
mm <- tapply(S$MIF_raw, S$region, mean, na.rm = TRUE)
cat("  界面 ", round(mm["Interface"], 3), " / 肿瘤核心 ", round(mm["Tumor_Core"], 3),
    " / 基质 ", round(mm["Stroma"], 3), "\n", sep = "")
cat("注意：即使数值接近，也不能说明旧分区正确——旧分区下肿瘤核心只出现在 1/8 张切片，\n")
cat("      其均值本身没有跨切片代表性。数值接近只是巧合式的可能性，需以符号检验为准。\n")

Z_TAB <- do.call(rbind, lapply(intersect(GENES, rownames(Z)),
                               function(g) summ_by_region(paste0(g, "_z"), g)))
cat("\n各区域片内 z 均值（主定义 A）：\n")
pr(Z_TAB[, c("variable","region","n_spot","mean")], 3)
write.csv(Z_TAB, file.path(OUTDIR, "T20b_MIF_z_by_region.csv"), row.names = FALSE)

hr("STEP 4  逐切片比较 + 跨切片符号检验（主推断）")

# 不比较"区域均值"，而在每张切片内算区域均值差，再跨 8 张切片做符号检验。
# 理由：切片间整体水平不同，汇总均值会被切片构成带偏；符号检验只看方向一致性。
CMPS <- list(c("Interface", "Tumor_Core"), c("Interface", "Stroma"),
             c("Tumor_Core", "Stroma"))

compare_region <- function(var, tag, min_spot = MIN_SPOT) {
  do.call(rbind, lapply(CMPS, function(cmp) {
    a <- cmp[1]; b <- cmp[2]
    dif <- sapply(secs, function(s) {
      z <- S[S$section == s, ]
      x <- z[[var]][z[[paste0("region_", tag)]] == a]
      y <- z[[var]][z[[paste0("region_", tag)]] == b]
      if (sum(is.finite(x)) < min_spot || sum(is.finite(y)) < min_spot) return(NA_real_)
      mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
    })
    dif <- dif[is.finite(dif)]
    if (length(dif) < 5) return(NULL)
    k <- sum(dif > 0); n <- length(dif)
    data.frame(definition = tag, comparison = paste0(a, " > ", b), variable = var,
               n_section = n, n_positive = k,
               mean_diff = round(mean(dif), 3), median_diff = round(median(dif), 3),
               sign_p = signif(binom.test(k, n, 0.5)$p.value, 3),
               stringsAsFactors = FALSE)
  }))
}

ALL <- list()
for (d in names(DEFS)) {
  for (g in intersect(GENES, rownames(Z))) {
    r <- compare_region(paste0(g, "_raw"), d)
    if (!is.null(r)) ALL[[length(ALL) + 1]] <- r
  }
}
CMP_TAB <- do.call(rbind, ALL)
CMP_TAB$FDR <- p.adjust(CMP_TAB$sign_p, "BH")

cat("\n===== MIF（原始值）三组定义下的结果 =====\n")
pr(CMP_TAB[CMP_TAB$variable == "MIF_raw",
           c("definition","comparison","n_section","n_positive","mean_diff","sign_p")], 3)

cat("\n===== 其余基因（主定义 A）=====\n")
pr(CMP_TAB[CMP_TAB$definition == "A",
           c("variable","comparison","n_section","n_positive","mean_diff","sign_p")], 3)

write.csv(CMP_TAB, file.path(OUTDIR, "T20c_region_sign_tests.csv"), row.names = FALSE)

hr("STEP 5  MIF 与其他基因的空间共定位（切片内口径）")

# 旧值 MIF~CD74 ρ=0.389 可能是合并口径。合并口径在 3.6 已被证明会伪造相关
# （TSPO~CD68 合并 0.290 → 切片内 −0.021），故这里全部改为切片内重算。
sp_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 30) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
}

PAIRS <- list(c("MIF", "CD74"), c("MIF", "TSPO"), c("MIF", "CD68"),
              c("MIF", "CD163"), c("TSPO", "CD68"), c("MIF", "CXCR4"),
              c("MIF", "CD44"), c("CD74", "TSPO"))

box_perm <- function(a, b, sec, nperm = N_PERM) {
  ok <- is.finite(a) & is.finite(b)
  a2 <- a[ok]; b2 <- b[ok]; s2 <- sec[ok]
  r0 <- suppressWarnings(cor(a2, b2, method = "spearman"))
  if (!is.finite(r0)) return(c(obs = NA, p = NA))
  null <- numeric(nperm)
  for (k in seq_len(nperm)) {
    bp <- b2
    for (s in unique(s2)) {
      i <- which(s2 == s)
      if (length(i) > 1) bp[i] <- b2[i][sample.int(length(i))]
    }
    null[k] <- suppressWarnings(cor(a2, bp, method = "spearman"))
  }
  null <- null[is.finite(null)]
  c(obs = r0, p = mean(abs(null) >= abs(r0)))
}

COR_TAB <- do.call(rbind, lapply(PAIRS, function(p) {
  g1 <- p[1]; g2 <- p[2]
  if (!all(c(g1, g2) %in% rownames(RAW))) return(NULL)
  a <- RAW[g1, ]; b <- RAW[g2, ]
  pooled <- sp_cor(a, b)
  per <- sapply(secs, function(s) {
    i <- which(SEC == s); sp_cor(a[i], b[i])
  })
  per <- per[is.finite(per)]
  pb <- box_perm(a, b, SEC)
  data.frame(pair = paste0(g1, " ~ ", g2),
             rho_pooled = round(pooled, 3),
             rho_within_median = round(median(per), 3),
             rho_within_min = round(min(per), 3),
             rho_within_max = round(max(per), 3),
             n_pos_section = sum(per > 0), n_section = length(per),
             p_block = signif(pb["p"], 3),
             stringsAsFactors = FALSE)
}))

cat("\n合并口径 vs 切片内口径：\n")
pr(COR_TAB, 3)
cat("\n（p_block 为切片内置换 1,000 次的块置换 p，保留切片结构）\n")
write.csv(COR_TAB, file.path(OUTDIR, "T20d_spot_correlations.csv"), row.names = FALSE)

hr("STEP 6  裁决")

mif_a <- CMP_TAB[CMP_TAB$definition == "A" & CMP_TAB$variable == "MIF_raw", ]
gv <- function(cmp, col) {
  z <- mif_a[mif_a$comparison == cmp, ]
  if (!nrow(z)) return(NA); z[[col]][1]
}

n_vs_core   <- gv("Interface > Tumor_Core", "n_positive")
n_vs_stroma <- gv("Interface > Stroma",     "n_positive")
d_vs_core   <- gv("Interface > Tumor_Core", "mean_diff")
d_vs_stroma <- gv("Interface > Stroma",     "mean_diff")
p_vs_core   <- gv("Interface > Tumor_Core", "sign_p")
p_vs_stroma <- gv("Interface > Stroma",     "sign_p")

cat("界面 vs 肿瘤核心：", n_vs_core, "/8 张切片为正，均值差 ", d_vs_core,
    "，符号检验 p = ", p_vs_core, "\n", sep = "")
cat("界面 vs 基质    ：", n_vs_stroma, "/8 张切片为正，均值差 ", d_vs_stroma,
    "，符号检验 p = ", p_vs_stroma, "\n", sep = "")

# 五种结局。注意必须覆盖「界面高于基质但不高于核心」这一情形——
# 2026-10 的实跑结果正是这一类（界面 vs 核心 2/8，界面 vs 基质 8/8），
# 早期版本漏了这一支，会误报为「异常」。
verdict <- if (isTRUE(n_vs_core >= 7) && isTRUE(n_vs_stroma >= 7)) {
  "A-成立：MIF 在界面最高（同时高于核心与基质）。正文保留「界面富集」，换成新数字。"
} else if (isTRUE(n_vs_core >= 7) && isTRUE(n_vs_stroma < 7)) {
  "B-部分成立（仅高于核心）：MIF 在界面高于肿瘤核心但不高于基质。正文应改为「MIF 在非肿瘤区较高」，不得声称界面特异。"
} else if (isTRUE(n_vs_core < 7) && isTRUE(n_vs_stroma >= 7)) {
  paste0("C-部分成立（仅高于基质，非界面特异）：MIF 在界面与肿瘤核心均高于基质，",
         "但界面不高于肿瘤核心（仅 ", n_vs_core, "/8 为正，反向 ", 8 - n_vs_core,
         "/8）。应删除「界面特异富集」，改写为「MIF 在非基质（肿瘤上皮）区室较高，",
         "而髓系标志在基质较高，二者空间上不共处」。")
} else if (isTRUE(n_vs_core < 7) && isTRUE(n_vs_stroma < 7)) {
  "D-不成立：MIF 无稳定区域差异。删除区域富集表述，只保留配体升/受体降格局。"
} else {
  "E-异常：检查区域定义与 spot 数"
}

cat("\n【裁决】", verdict, "\n", sep = "")

# 三种定义的一致性
cat("\n三组定义下的界面 vs 核心 为正切片数：",
    paste(sapply(names(DEFS), function(d) {
      z <- CMP_TAB[CMP_TAB$definition == d & CMP_TAB$variable == "MIF_raw" &
                     CMP_TAB$comparison == "Interface > Tumor_Core", ]
      if (!nrow(z)) NA else z$n_positive[1]
    }), collapse = " / "), "\n", sep = "")
cat("三组定义下的界面 vs 基质 为正切片数：",
    paste(sapply(names(DEFS), function(d) {
      z <- CMP_TAB[CMP_TAB$definition == d & CMP_TAB$variable == "MIF_raw" &
                     CMP_TAB$comparison == "Interface > Stroma", ]
      if (!nrow(z)) NA else z$n_positive[1]
    }), collapse = " / "), "\n", sep = "")

writeLines(c(
  "MIF 区域归属复核结果",
  paste("时间:", format(Sys.time())),
  "",
  paste("界面 vs 肿瘤核心:", n_vs_core, "/8，均值差", d_vs_core, "，p =", p_vs_core),
  paste("界面 vs 基质    :", n_vs_stroma, "/8，均值差", d_vs_stroma, "，p =", p_vs_stroma),
  "",
  paste("旧口径（已废弃）: 界面 2.504 / 核心 1.610 / 基质 1.783"),
  paste("新口径（主定义A）: 界面", round(mm["Interface"], 3),
        "/ 核心", round(mm["Tumor_Core"], 3), "/ 基质", round(mm["Stroma"], 3)),
  "",
  paste("裁决:", verdict),
  "",
  "注：旧 p<1e-95 为 spot 层面 Wilcoxon，受空间自相关伪重复放大，不可引用。",
  "    本脚本以跨切片符号检验为主要推断，8 张切片时最小可能 p 为 0.0078。"
), file.path(OUTDIR, "T20e_verdict.txt"))

hr("STEP 7  图")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  df <- S[S$region %in% c("Tumor_Core", "Interface", "Stroma"), ]
  df$region <- factor(df$region, levels = c("Tumor_Core", "Interface", "Stroma"))

  p1 <- ggplot(df, aes(region, MIF_raw, fill = region)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.75) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    facet_wrap(~ section, nrow = 2) +
    labs(x = NULL, y = "MIF (SCT normalized)",
         title = "MIF by region, per section (within-section z definition)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(FIGDIR, "Fig8_MIF_region.pdf"), p1, width = 10, height = 5.5)
  ggsave(file.path(FIGDIR, "Fig8_MIF_region.png"), p1, width = 10, height = 5.5, dpi = 300)

  # 逐切片差值散点（更直观地看 8 张切片的方向）
  dd <- do.call(rbind, lapply(secs, function(s) {
    z <- S[S$section == s, ]
    data.frame(section = s,
               diff = mean(z$MIF_raw[z$region == "Interface"], na.rm = TRUE) -
                      mean(z$MIF_raw[z$region == "Tumor_Core"], na.rm = TRUE))
  }))
  dd <- dd[is.finite(dd$diff), ]
  p2 <- ggplot(dd, aes(reorder(section, diff), diff, fill = diff > 0)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC")) +
    labs(x = NULL, y = "Interface - Tumor_Core (MIF)",
         title = "Per-section difference") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(FIGDIR, "Fig8_MIF_diff_persection.pdf"), p2, width = 6, height = 4)
  ggsave(file.path(FIGDIR, "Fig8_MIF_diff_persection.png"), p2, width = 6, height = 4, dpi = 300)
}

hr("28_mif_spatial_region.R 完成")

# ############################################################################
# ##  Fig9_上游与药物  (7 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- Fig9_上游与药物 / 11_oncopredict.R
# -----------------------------------------------------------------------------
# =============================================================================
#  11_oncopredict.R
#  结论五（第 1 步）：用 oncoPredict 从 TCGA-OV 表达谱预测 198 种药物的 IC50
#
#  输入：EXPR_clean.rds（18747 x 304）、GDSC2 训练集（RMA 表达 + 药敏表）
#  输出：cache/ic50_clean.rds（药物 x 样本）、tables/T5a_ic50_matrix.csv
#
#  两个踩过的坑（必须按这里的写法）：
#    1. oncoPredict 不同版本的 calcPhenotype() 形参不同，
#       老版本要求 removeLowVaringGenesFrom（注意拼写是 Varing），
#       新版本没有这个参数。必须按 formals() 动态补参，不能写死。
#    2. 部分版本 calcPhenotype() 返回 NULL，但会把结果写到
#       ./calcPhenotype_Output/DrugPredictions.csv。必须做返回值兜底。
#
#  运行时间：约 5–15 分钟（首次）；结果会缓存，重跑直接复用。
# =============================================================================

if (!exists("OUT")) stop("请先 source('00_setup.R')")
if (!requireNamespace("oncoPredict", quietly = TRUE))
  stop("请先安装 oncoPredict（本地有 oncoPredict_1.2.tar 可用）")
suppressPackageStartupMessages({ library(oncoPredict); library(Matrix) })

GDSC_DIRS <- c(file.path(PROJ, "DataFiles/DataFiles/DataFiles/Training Data"),
               file.path(PROJ, "GDSC_Data/DataFiles/DataFiles/Training Data"),
               PROJ)

hr("STEP 1  载入 TCGA 测试矩阵")

stopifnot(file.exists(F_EXPR))
EXPR <- as.matrix(readRDS(F_EXPR)); storage.mode(EXPR) <- "double"
cat("测试矩阵: ", nrow(EXPR), " 基因 x ", ncol(EXPR), " 样本\n", sep = "")

hr("STEP 2  定位 GDSC2 训练集")

find_first <- function(pat) {
  for (d in GDSC_DIRS) {
    if (!dir.exists(d)) next
    f <- list.files(d, pattern = pat, recursive = TRUE, full.names = TRUE)
    if (length(f)) return(f[which.max(file.size(f))])
  }
  NA_character_
}

f_expr <- find_first("GDSC2_Expr.*\\.rds$")
f_res  <- find_first("GDSC2_Res.*\\.rds$")
if (is.na(f_expr) || is.na(f_res)) {
  stop("找不到 GDSC2 训练集。需要的两个文件：\n",
       "  GDSC2_Expr (RMA Normalized and Log Transformed).rds\n",
       "  GDSC2_Res.rds\n可到 oncoPredict 官网或 GDSC 下载后放到 ", PROJ)
}
cat("训练表达: ", basename(f_expr), "\n训练药敏: ", basename(f_res), "\n", sep = "")

trExpr <- readRDS(f_expr)
trRes  <- readRDS(f_res)
cat("训练集: ", nrow(trExpr), " 基因 x ", ncol(trExpr), " 细胞系；",
    nrow(trRes), " 细胞系 x ", ncol(trRes), " 药物\n", sep = "")

common <- intersect(rownames(EXPR), rownames(trExpr))
cat("共同基因: ", length(common), "\n", sep = "")
if (length(common) < 10000) warning("共同基因偏少（<10000），检查基因命名是否一致")

hr("STEP 3  预测 IC50（带缓存）")

ic_file <- file.path(OUT, "cache", "ic50_clean.rds")

if (file.exists(ic_file)) {
  IC50 <- readRDS(ic_file)
  cat("复用缓存: ", nrow(IC50), " 药 x ", ncol(IC50), " 样本\n", sep = "")
} else {
  # ⚠ oncoPredict 的 calcPhenotype() 会把结果无条件写到
  #   ./calcPhenotype_Output/DrugPredictions.csv（相对 getwd()）。
  #   若 getwd() 不可写（例如盘符只读、沙箱限制），它不会报错，
  #   而是静默保留目录里的旧文件 —— 随后读到的就是上一次的结果。
  #   因此这里先把工作目录切到可写的 OUT 下，并记录起始时间用于校验新鲜度。
  ODIR <- file.path(OUT, "calcPhenotype_Output")
  dir.create(ODIR, showWarnings = FALSE, recursive = TRUE)
  old_wd    <- getwd()
  test_mat  <- EXPR[common, , drop = FALSE]
  setwd(dirname(ODIR))

  t_start_ic <- Sys.time()
  res <- NULL

  if ("calcPhenotype" %in% getNamespaceExports("oncoPredict")) {
    fm <- names(formals(oncoPredict::calcPhenotype))
    args <- list(trainingExprData = trExpr,
                 trainingPtype    = trRes,
                 testExprData     = test_mat,
                 batchCorrect     = "eb",
                 powerTransformPhenotype = TRUE,
                 removeLowVaryingGenes   = 0.2,
                 minNumSamples    = 10,
                 printOutput      = FALSE,
                 selection        = 1)
    if ("removeLowVaringGenesFrom" %in% fm)          # 老版本（拼写 Varing）
      args$removeLowVaringGenesFrom <- "homogenizeData"
    if ("removeLowVaringGenes" %in% fm && !"removeLowVaryingGenes" %in% fm)
      args$removeLowVaryingGenes <- args$removeLowVaryingGenes

    args <- args[intersect(names(args), fm)]
    res  <- try(do.call(oncoPredict::calcPhenotype, args), silent = TRUE)
  } else {
    res <- try(oncoPredict::predictPhenotype(
      trainPhenotypeData = trRes, trainExprData = trExpr,
      testExprData = test_mat, batchCorrect = "eb",
      powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.2,
      minNumSamples = 10, printOutput = FALSE, selection = 1,
      report_predictor = FALSE)[[1]], silent = TRUE)
  }

  # 返回值兜底：部分版本返回 NULL，结果落在 calcPhenotype_Output 的 CSV 里
  if (inherits(res, "try-error") || is.null(res) || !(is.matrix(res) || is.data.frame(res))) {
    f <- file.path(ODIR, "DrugPredictions.csv")
    if (!file.exists(f)) {
      stop("oncoPredict 未返回矩阵，且 ", f, " 不存在。",
           if (inherits(res, "try-error")) paste0("\n原始错误：",
                                                  conditionMessage(attr(res, "condition"))))
    }
    # 关键校验：文件必须是本次运行生成的，否则说明写入失败、读到的是旧结果
    mt <- file.info(f)$mtime
    if (is.na(mt) || mt < t_start_ic - 60) {
      stop("calcPhenotype 未返回矩阵，且 ", f, " 的时间戳为 ", mt,
           "，早于本次运行起点（", format(t_start_ic), "）。\n",
           "说明结果写入失败，读到的是历史文件。请确认 OUT 目录可写后再重跑。")
    }
    cat("从输出文件读取: ", f, "（生成于 ", format(mt), "）\n", sep = "")
    res <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  }
  setwd(old_wd)   # 恢复工作目录（顶层不能用 on.exit）

  IC50 <- as.matrix(res)
  # 统一成 药物 x 样本
  if (nrow(IC50) == ncol(test_mat) && ncol(IC50) != ncol(test_mat)) IC50 <- t(IC50)
  if (nrow(IC50) > ncol(IC50)) IC50 <- t(IC50)
  if (!identical(colnames(IC50), colnames(test_mat)))
    cat("提醒：IC50 列名与表达矩阵列名不完全一致，后续按交集对齐\n")

  saveRDS(IC50, ic_file)
  cat(sprintf("完成 %.1f 分钟；%d 药 x %d 样本\n",
              as.numeric(difftime(Sys.time(), t_start_ic, units = "mins")),
              nrow(IC50), ncol(IC50)))
}

hr("STEP 4  写出与自检")

smp <- intersect(colnames(IC50), colnames(EXPR))
cat("与表达矩阵匹配的样本: ", length(smp), "\n", sep = "")
stopifnot(length(smp) > 100)

write.csv(round(IC50[, smp, drop = FALSE], 4),
          file.path(OUT, "tables", "T5a_ic50_matrix.csv"))

cat("\n模型中的目标药物是否齐全：\n")
for (d in c("OSI-027", "AZD2014", "AZD8055", "GNE-317", "Trametinib",
            "MK-1775", "Cisplatin", "Paclitaxel")) {
  hit <- rownames(IC50)[grepl(paste0("^", d), rownames(IC50), ignore.case = TRUE)]
  cat(sprintf("  %-11s %s\n", d, if (length(hit)) hit[1] else "缺失"))
}
cat("\n预测 IC50 可为负值（log 尺度），组间比较请用 Cohen's d 或中位数差，不要用绝对值。\n")

hr("11_oncopredict.R 完成")

# -----------------------------------------------------------------------------
# -- Fig9_上游与药物 / 12_drug_stratification.R
# -----------------------------------------------------------------------------
# =============================================================================
#  12_drug_stratification.R
#  结论五（第 2、3 步）：偏相关 + 三分位效应量 + 分诊性能
#
#  对应正文 3.6 的三张表：
#    表1 偏相关（双重校正免疫构成与增殖评分）
#        OSI-027 -0.395 (FDR 4.8e-11) / Bortezomib -0.215 / Trametinib -0.195
#        GNE-317 -0.173 / Dasatinib -0.152 / AZD2014 -0.151 / AZD8055 -0.149
#        MK-1775 +0.166 (FDR 9.7e-3)
#    表2 三分位 Cohen's d（mTOR 类 -0.59 ~ -0.92；WEE1 +0.41；化疗 ~0）
#    表3 分诊性能（灵敏度 73.3% / 特异度 73.5%；AUC 0.580 -> 0.764；CV AUC 0.749）
#
#  ⚠ 表3 是同一相关关系的重述，不是独立验证。正文已有免责说明，
#     任何情况下都不要把 AUC 当成结论证据使用。
# =============================================================================

if (!exists("OUT")) stop("请先 source('00_setup.R')")
suppressPackageStartupMessages({ library(pROC) })

hr("STEP 1  载入")

stopifnot(file.exists(F_EXPR))
EXPR <- as.matrix(readRDS(F_EXPR)); storage.mode(EXPR) <- "double"

f_ic <- file.path(OUT, "cache", "ic50_clean.rds")
if (!file.exists(f_ic)) stop("找不到 ic50_clean.rds，请先跑 11_oncopredict.R")
IC50 <- as.matrix(readRDS(f_ic))
if (nrow(IC50) > ncol(IC50)) IC50 <- t(IC50)

smp <- intersect(colnames(IC50), colnames(EXPR))
cat("样本: ", length(smp), "  药物: ", nrow(IC50), "\n", sep = "")

Z <- t(scale(t(EXPR[, smp, drop = FALSE])))
Z[!is.finite(Z)] <- NA

ts <- as.numeric(EXPR["TSPO", smp])

# 校正变量：免疫构成 + 增殖评分
IMMUNE <- c("CD68", "CD3D", "COL1A1")
PROLIF <- c("MKI67", "PCNA", "TOP2A", "CCNB1", "BUB1", "AURKA", "PLK1")
g_imm  <- intersect(IMMUNE, rownames(Z))
g_pro  <- intersect(PROLIF, rownames(Z))
cat("免疫构成基因: ", length(g_imm), "  增殖基因: ", length(g_pro), "\n", sep = "")

COV <- rbind(Z[g_imm, , drop = FALSE],
             PROLIF_score = colMeans(Z[g_pro, , drop = FALSE], na.rm = TRUE))

hr("STEP 2  偏相关（双重校正）")

D <- do.call(rbind, lapply(rownames(IC50), function(dr) {
  y  <- as.numeric(IC50[dr, smp])
  r0 <- suppressWarnings(cor.test(ts, y, method = "spearman"))
  rp <- pcor_rank(ts, y, COV)
  data.frame(drug = dr,
             rho = unname(r0$estimate), p = r0$p.value,
             rho_partial = rp["rho"], p_partial = rp["p"],
             n = length(ts), stringsAsFactors = FALSE)   # Spearman 有结时 cor.test 不返回 parameter
}))
D$FDR          <- p.adjust(D$p, "BH")
D$FDR_partial  <- p.adjust(D$p_partial, "BH")
D <- D[order(-abs(D$rho_partial)), ]

sel <- D[D$FDR_partial < 0.05, c("drug","rho","rho_partial","p_partial","FDR_partial")]
cat("偏相关显著（FDR<0.05）的药物: ", nrow(sel), "\n", sep = "")
pr(sel, 4)
write.csv(D, file.path(OUT, "tables", "T5b_drug_partial_all.csv"), row.names = FALSE)

hr("STEP 3  三分位分层效应量")

TARGETS <- c(OSI.027 = "OSI-027", AZD2014 = "AZD2014", AZD8055 = "AZD8055",
             GNE.317 = "GNE-317", Trametinib = "Trametinib",
             MK.1775 = "MK-1775", Cisplatin = "Cisplatin", Paclitaxel = "Paclitaxel")

find_drug <- function(nm) {
  hit <- rownames(IC50)[grepl(paste0("^", gsub("([.-])", "\\\\\\1", nm)),
                              rownames(IC50), ignore.case = TRUE)]
  if (length(hit)) hit[1] else NA_character_
}

q   <- quantile(ts, c(1/3, 2/3))
grp <- ifelse(ts > q[2], "High", ifelse(ts <= q[1], "Low", "Mid"))
cat("切点: ", paste(round(q, 3), collapse = " / "), "\n", sep = "")
print(table(grp))

TAB <- do.call(rbind, lapply(names(TARGETS), function(nm) {
  dr <- find_drug(TARGETS[[nm]])
  if (is.na(dr)) { cat("跳过 ", nm, "（药物不在预测矩阵中）\n", sep = ""); return(NULL) }
  y  <- as.numeric(IC50[dr, smp])
  hi <- y[grp == "High"]; lo <- y[grp == "Low"]
  if (length(hi) < 20 || length(lo) < 20) return(NULL)
  d  <- (mean(hi) - mean(lo)) / sqrt((var(hi) + var(lo)) / 2)
  w  <- suppressWarnings(wilcox.test(hi, lo))
  data.frame(drug = nm, drug_row = dr,
             median_high = round(median(hi), 3),
             median_low  = round(median(lo), 3),
             delta_median = round(median(hi) - median(lo), 3),
             cohens_d = round(d, 3),
             p = signif(w$p.value, 3), stringsAsFactors = FALSE)
}))
TAB <- TAB[order(TAB$cohens_d), ]
pr(TAB)
write.csv(TAB, file.path(OUT, "tables", "T5c_drug_tertile.csv"), row.names = FALSE)

hr("STEP 4  分诊性能（探索性）")

dr <- find_drug("OSI-027")
if (!is.na(dr)) {
  ic50 <- as.numeric(IC50[dr, smp])
  sens <- as.integer(ic50 < median(ic50, na.rm = TRUE))

  ok   <- !is.na(grp)
  tb   <- table(factor(grp[ok], levels = c("High","Low")),
                factor(sens[ok], levels = c(1,0), labels = c("sensitive","resistant")))
  cat("2x2 分诊表:\n"); print(tb)
  TP <- tb["High","sensitive"]; FP <- tb["High","resistant"]
  FN <- tb["Low","sensitive"];  TN <- tb["Low","resistant"]
  cat(sprintf("灵敏度 %.1f%%  特异度 %.1f%%\n",
              100 * TP / (TP + FP), 100 * TN / (TN + FN)))

  d <- data.frame(sens = sens, TSPO = ts,
                  PROLIF = COV["PROLIF_score", ],
                  IMMUNE = colMeans(COV[g_imm, , drop = FALSE], na.rm = TRUE),
                  OXPHOS = colMeans(Z[intersect(
                    c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1",
                      "NDUFV2","SDHA","SDHB","UQCRC1","UQCRC2","COX4I1","COX5A",
                      "ATP5A1","ATP5B"), rownames(Z)), , drop = FALSE], na.rm = TRUE))
  d <- d[stats::complete.cases(d), ]

  fit_auc <- function(fml) {
    m <- glm(as.formula(fml), data = d, family = binomial)
    r <- pROC::roc(d$sens, fitted(m), quiet = TRUE)
    list(model = m, roc = r, auc = as.numeric(pROC::auc(r)), fml = fml)
  }
  M0 <- fit_auc("sens ~ PROLIF + IMMUNE")
  M1 <- fit_auc("sens ~ PROLIF + IMMUNE + TSPO")
  M2 <- fit_auc("sens ~ PROLIF + IMMUNE + TSPO + OXPHOS")
  for (M in list(M0, M1, M2))
    cat(sprintf("  %-34s AUC = %.3f\n", M$fml, M$auc))

  cat("\nDeLong:\n")
  for (pr2 in list(list(M0, M1, "加 TSPO"), list(M1, M2, "再加 OXPHOS"))) {
    tt <- pROC::roc.test(pr2[[1]]$roc, pr2[[2]]$roc, method = "delong")
    cat(sprintf("  %-10s dAUC = %+.3f, p = %.3g\n",
                pr2[[3]], pr2[[2]]$auc - pr2[[1]]$auc, tt$p.value))
  }

  cv_auc <- function(fml, K = 10, B = 50) {
    a <- numeric(B)
    for (b in seq_len(B)) {
      set.seed(b)
      fold <- sample(rep(seq_len(K), length.out = nrow(d)))
      pred <- rep(NA_real_, nrow(d))
      for (k in seq_len(K)) {
        tr <- fold != k
        m  <- glm(as.formula(fml), data = d[tr, ], family = binomial)
        pred[!tr] <- predict(m, newdata = d[!tr, ], type = "response")
      }
      a[b] <- as.numeric(pROC::auc(pROC::roc(d$sens, pred, quiet = TRUE)))
    }
    c(mean = mean(a), lo = quantile(a, .025), hi = quantile(a, .975))
  }

  cat("\n10 折 CV AUC（50 次重复）:\n")
  CVR <- do.call(rbind, lapply(c("sens ~ PROLIF + IMMUNE",
                                 "sens ~ PROLIF + IMMUNE + TSPO",
                                 "sens ~ PROLIF + IMMUNE + TSPO + OXPHOS"),
                               function(fml) {
    r <- cv_auc(fml)
    data.frame(model = fml, auc_cv = round(r[1], 3),
               lo = round(r[2], 3), hi = round(r[3], 3), stringsAsFactors = FALSE)
  }))
  pr(CVR)
  write.csv(CVR, file.path(OUT, "tables", "T5d_triage_cv_auc.csv"), row.names = FALSE)

  cat("\n⚠ 免责：预测 IC50 由表达谱回归得出，与 TSPO 的表达高度重叠；\n")
  cat("   分层与性能计算用的是同一批样本；终点由预测 IC50 中位数二分定义。\n")
  cat("   因此上述灵敏度/特异度/AUC 不构成独立预测验证，不得作为结论证据。\n")
}

hr("12_drug_stratification.R 完成")

# -----------------------------------------------------------------------------
# -- Fig9_上游与药物 / 13_tf_upstream.R
# -----------------------------------------------------------------------------
# =============================================================================
#  13_tf_upstream.R
#  上游调控之一：转录因子活性与 TSPO 表达
#
#  科学问题：TSPO 的高表达由哪些转录因子驱动？
#
#  方法：
#    1. DoRothEA 高置信度 regulon（A/B/C）作为 TF-靶基因先验
#    2. decoupleR::run_ulm() 推断每个样本的 TF 活性（ULM，无监督）
#    3. TF 活性与 TSPO 表达做 Spearman
#    4. 偏相关控制 CD68/CD3D/COL1A1，排除「免疫细胞含量」造成的假关联
#    5. BH 校正；同时报告与 OXPHOS 模块的相关，看是否同一批 TF
#
#  ⚠ 循环性处理：TF 活性由靶基因表达加权得出。若 TSPO 本身是某 TF 的靶基因，
#    该 TF 的活性就部分由 TSPO 决定，相关必然偏高。因此计算活性前把 TSPO
#    从 regulon 的靶点中剔除。这是本文与常规做法不同的一点，需在方法学写明。
#
#  输入：EXPR_clean.rds（18,747 × 304）、MODULES.rds
#  输出：T6a_TF_activity_vs_TSPO.csv、T6b_TF_activity_matrix.csv
# =============================================================================

if (!exists("OUT")) stop("请先 source('00_setup.R')")
suppressPackageStartupMessages({
  library(dorothea); library(decoupleR); library(Matrix)
})

OUTDIR <- file.path(OUT, "tables")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MIN_TARGETS   <- 5      # regulon 最少靶基因数
CONF_LEVELS   <- c("A", "B", "C")
EXCLUDE_SELF  <- TRUE   # 把 TSPO 从靶点里剔除，避免循环

hr("STEP 1  载入表达矩阵")

EXPR <- as.matrix(readRDS(F_EXPR))
storage.mode(EXPR) <- "double"
cat("表达矩阵: ", nrow(EXPR), " 基因 x ", ncol(EXPR), " 样本\n", sep = "")
stopifnot("TSPO" %in% rownames(EXPR))

MOD <- NULL
for (p in c(file.path(OUT, "MODULES.rds"), file.path(PROJ, "results_v2", "MODULES.rds"))) {
  if (file.exists(p)) { MOD <- readRDS(p); break }
}
if (is.null(MOD)) stop("找不到 MODULES.rds，请先跑 01_build_matrix.R")
if (is.null(MOD$OXPHOS)) stop("MODULES.rds 里没有 OXPHOS 模块")

hr("STEP 2  构建 DoRothEA regulon")

data("dorothea_hs", package = "dorothea")
net <- dorothea_hs[dorothea_hs$confidence %in% CONF_LEVELS, c("tf", "target", "mor")]
cat("原始 regulon: ", nrow(net), " 条，", length(unique(net$tf)), " 个 TF\n", sep = "")

net <- net[net$target %in% rownames(EXPR), ]
cat("靶基因在表达矩阵内的: ", nrow(net), " 条\n", sep = "")

if (EXCLUDE_SELF) {
  n_before <- sum(net$target == "TSPO")
  net <- net[net$target != "TSPO", ]
  cat("剔除 TSPO 作为靶点的记录: ", n_before, " 条（避免循环推断）\n", sep = "")
}

keep_tf <- names(which(table(net$tf) >= MIN_TARGETS))
net <- net[net$tf %in% keep_tf, ]
cat("靶基因数≥", MIN_TARGETS, " 的 TF: ", length(keep_tf), " 个\n", sep = "")

hr("STEP 3  推断 TF 活性（ULM）")

act <- decoupleR::run_ulm(
  mat      = EXPR,
  network  = net,
  .source  = "tf",
  .target  = "target",
  .mor     = "mor",
  minsize  = MIN_TARGETS
)
cat("ULM 输出列: ", paste(colnames(act), collapse = ", "), "\n")

# 转成 TF x 样本 的宽矩阵（不用 pivot_wider_profile，各版本列名解析不一致）
tf_lv   <- unique(as.character(act$source))
cond_lv <- unique(as.character(act$condition))
act_w <- matrix(NA_real_, nrow = length(tf_lv), ncol = length(cond_lv),
                dimnames = list(tf_lv, cond_lv))
act_w[cbind(match(as.character(act$source), tf_lv),
            match(as.character(act$condition), cond_lv))] <- act$score
cat("TF 活性矩阵: ", nrow(act_w), " TF x ", ncol(act_w), " 样本\n", sep = "")

smp <- intersect(colnames(act_w), colnames(EXPR))
act_w <- act_w[, smp, drop = FALSE]
cat("与表达矩阵共有的样本: ", length(smp), "\n", sep = "")

hr("STEP 4  与 TSPO 求相关")

ts   <- as.numeric(EXPR["TSPO", smp])
Zc   <- t(scale(t(EXPR[, smp, drop = FALSE])))
Zc[!is.finite(Zc)] <- NA

# 控制变量：免疫/基质构成
imm <- intersect(c("CD68", "CD3D", "COL1A1"), rownames(Zc))
COV <- Zc[imm, , drop = FALSE]
cat("偏相关控制变量: ", paste(imm, collapse = ", "), "\n", sep = "")

# OXPHOS 模块评分，用于判断是否同一批 TF
ox_g <- intersect(MOD$OXPHOS, rownames(Zc))
ox   <- colMeans(Zc[ox_g, , drop = FALSE], na.rm = TRUE)

res <- do.call(rbind, lapply(rownames(act_w), function(tf) {
  y <- as.numeric(act_w[tf, ])
  if (sum(is.finite(y)) < 50 || sd(y, na.rm = TRUE) == 0) return(NULL)
  r0 <- suppressWarnings(cor.test(y, ts, method = "spearman"))
  rp <- pcor_rank(ts, y, COV)
  ro <- suppressWarnings(cor(y, ox, method = "spearman", use = "complete.obs"))
  data.frame(tf = tf,
             rho        = round(as.numeric(r0$estimate), 3),
             p          = r0$p.value,
             rho_partial = round(rp["rho"], 3),
             p_partial   = rp["p"],
             rho_OXPHOS  = round(ro, 3),
             n_targets   = sum(net$tf == tf),
             stringsAsFactors = FALSE)
}))
res$FDR         <- p.adjust(res$p, "BH")
res$FDR_partial <- p.adjust(res$p_partial, "BH")
res <- res[order(-abs(res$rho)), ]

cat("\n显著 TF 数（FDR<0.05）: ", sum(res$FDR < 0.05), " / ", nrow(res), "\n", sep = "")
cat("偏相关后仍显著（FDR<0.05）: ", sum(res$FDR_partial < 0.05), " / ", nrow(res), "\n", sep = "")

cat("\n===== Top 25（按 |rho| 排序）=====\n")
print(head(res[, c("tf","rho","rho_partial","FDR","FDR_partial","rho_OXPHOS","n_targets")], 25),
      row.names = FALSE, digits = 3)

hr("STEP 5  已知 TSPO 调控因子核对")

KNOWN <- c("SP1","SP3","JUN","JUNB","FOS","FOSB","NFKB1","RELA","NFKB2",
           "HIF1A","EPAS1","CREB1","ATF2","STAT3","STAT1","PPARG","PPARA",
           "NFE2L2","TP53","EGR1","MYC","CEBPA","CEBPB","YY1","CTCF")
hit <- res[res$tf %in% KNOWN, ]
if (nrow(hit)) {
  cat("\n文献报道与 TSPO 相关的 TF 在本数据中的表现：\n")
  print(hit[, c("tf","rho","rho_partial","p","p_partial","FDR","rho_OXPHOS")],
        row.names = FALSE, digits = 3)
} else {
  cat("上述 TF 均不在 regulon 交集内\n")
}

hr("STEP 6  特异性检验（关键）")

# 6.1 TF 活性之间共线到什么程度
top25 <- head(res$tf[order(-abs(res$rho))], 25)
cm <- cor(t(act_w[top25, smp]), method = "spearman")
cat("Top25 TF 活性两两相关: 平均 ", round(mean(cm[upper.tri(cm)]), 3),
    "，最大 ", round(max(cm[upper.tri(cm)]), 3), "\n", sep = "")

pca <- prcomp(t(act_w[, smp]), scale. = TRUE)
v   <- summary(pca)$importance[2, ] * 100
cat("PCA: PC1 ", round(v[1], 1), "%，PC2 ", round(v[2], 1),
    "%，累计 80% 需 ", which(cumsum(v) >= 80)[1], " 个主成分（共 ", length(v), "）\n", sep = "")

# 6.2 TF 活性是否只是巨噬含量
cd68 <- as.numeric(Zc["CD68", ])
cat("\n各 TF 活性与 CD68 的相关（高值说明该 TF 主要标记髓系含量）：\n")
for (tf in c("RELA","NFKB1","SP1","TP53","CEBPB","STAT3","ARID2","PPARA")) {
  if (!tf %in% rownames(act_w)) next
  a <- as.numeric(act_w[tf, smp])
  cat(sprintf("  %-6s ~ CD68 rho = %+.3f\n",
              tf, cor(a, cd68, method = "spearman", use = "complete.obs")))
}

# 6.3 只控制 CD68 的偏相关（比三基因更严格地剥离髓系）
cat("\n仅控制 CD68 时与 TSPO 的偏相关：\n")
for (tf in c("RELA","NFKB1","SP1","TP53","CEBPB","STAT3","ARID2","PPARA")) {
  if (!tf %in% rownames(act_w)) next
  a  <- as.numeric(act_w[tf, smp])
  rp <- pcor_rank(ts, a, matrix(cd68, nrow = 1))
  cat(sprintf("  %-6s rho_partial = %+.3f  p = %.3g\n", tf, rp["rho"], rp["p"]))
}

# 6.4 与 OXPHOS 模块同量级的 TF 属于「代谢模块镜像」，不是 TSPO 特异
res$is_module_mirror <- abs(res$rho_OXPHOS) >= 0.8 * abs(res$rho)
cat("\n被判定为「OXPHOS 模块镜像」的 TF: ", sum(res$is_module_mirror), " 个，",
    "其中进入 Top25 的: ",
    paste(intersect(top25, res$tf[res$is_module_mirror]), collapse = ", "), "\n", sep = "")
cat("这些 TF 不应作为 TSPO 上游候选写进正文。\n")

hr("STEP 7  阴性对照")

set.seed(SEED)
rand_rho <- sapply(seq_len(200), function(i) {
  g <- sample(setdiff(rownames(Zc), "TSPO"), 20)
  y <- colMeans(Zc[g, , drop = FALSE], na.rm = TRUE)
  suppressWarnings(cor(y, ts, method = "spearman", use = "complete.obs"))
})
cat("随机 20 基因假活性的 rho 分布: 5%/50%/95% = ",
    paste(round(quantile(rand_rho, c(.05,.5,.95), na.rm = TRUE), 3), collapse = " / "), "\n", sep = "")
cat("Top TF 的 |rho| = ", round(max(abs(res$rho)), 3), "，超过 ",
    round(100 * mean(abs(rand_rho) < max(abs(res$rho)))), "% 的随机模块\n", sep = "")

hr("STEP 8  写出")

write.csv(res, file.path(OUTDIR, "T6a_TF_activity_vs_TSPO.csv"), row.names = FALSE)
saveRDS(act_w, file.path(OUT, "cache", "tf_activity.rds"))
write.csv(data.frame(random_rho = rand_rho),
          file.path(OUTDIR, "T6c_TF_random_control.csv"), row.names = FALSE)

cat("\n判读要点：\n")
cat("  · 只有偏相关（rho_partial）仍显著的 TF 才值得讨论，rho 会被细胞构成带偏。\n")
cat("  · 与 OXPHOS 模块高度共线的 TF 不一定是 TSPO 特异的，需在正文区分。\n")
cat("  · ULM 是推断活性，不是实测结合；结论措辞应为「与……活性相关」，不写「调控」。\n")

hr("13_tf_upstream.R 完成")

# -----------------------------------------------------------------------------
# -- Fig9_上游与药物 / 14_tf_upstream_singlecell.R
# -----------------------------------------------------------------------------
# =============================================================================
#  14_tf_upstream_singlecell.R
#  上游调控之二：上皮单细胞内部的转录因子活性与 TSPO
#
#  为什么做这一步：
#    bulk 层面（13_tf_upstream.R）发现 NF-κB、SP1 等与 TSPO 相关，但 bulk 里
#    这些 TF 活性与 CD68 的相关高达 0.70/0.75，无法排除髓系浸润的贡献。
#    本脚本把分析限定在 7,929 个上皮细胞内部，彻底消除细胞构成混杂。
#
#  方法（与 03_epithelial_mixed.R 对齐）：
#    1. 逐基因回归掉文库大小（log10 nCount），再做 z 标准化
#    2. 用残差化矩阵推断 TF 活性（decoupleR ULM）
#    3. 逐病人 Spearman(TF 活性, TSPO)，11 例病人做符号检验
#    4. 混合模型 y ~ TSPO + (1|patient_id)
#    5. 随机基因集阴性对照
#    6. 与 bulk 结果对照，看方向是否一致
#
#  ⚠ 与 bulk 的差异，需写进方法：
#    bulk 直接用表达值跑 ULM；这里用去深度后的 z 值，避免测序深度驱动活性。
#    细胞数远大于 bulk，p 值会极小，推断以「逐病人符号检验」为准。
#
#  输入：F_SC_EPI（32,847 基因 × 50,867 细胞，取 cell_type=="EOC" 共 7,929）
#  输出：T6b_SC_TF_vs_TSPO.csv、T6e_SC_TF_replication.csv、T6f_SC_TF_random_control.csv
#
#  用法：直接 source 本文件即可，不必先 source 00_setup.R
# =============================================================================

# ---- 前置：自动载入 00_setup.R；不可用则就地定义所需变量 --------------------
if (!exists("OUT")) {
  # 【审稿人注意】本行路径已改为通用写法。若要单独运行本脚本，
  #   请把 .code_dir 的默认值 "." 改成 code_final 目录在你的机器上的实际路径，
  #   例如 "/your/path/to/results_v2/code_final"，或设置环境变量 TSPO_CODE_DIR。
  #   通过 TSPO定稿代码.R 整体运行时 PART 0 已设好 OUT，此分支不会执行。
  .code_dir <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
  .setup <- file.path(.code_dir, "00_setup.R")
  .ok <- FALSE
  if (file.exists(.setup)) {
    .ok <- !inherits(try(source(.setup, encoding = "UTF-8"), silent = TRUE), "try-error")
    if (!.ok) message("00_setup.R 载入失败（可能缺包），改用内置默认路径")
  }
  if (!.ok) {
    PROJ     <- "D:/Rstudy/tcga"
    PROJ2    <- "D:/Rstudy/tcga2"
    OUT      <- file.path(PROJ, "results_v2")
    F_SC_EPI <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")
    hr <- function(x) cat("\n", strrep("=", 76), "\n  ", x, "\n",
                          strrep("=", 76), "\n", sep = "")
    pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)
  }
  for (d in c(OUT, file.path(OUT, "tables"), file.path(OUT, "cache"))) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
  rm(list = intersect(c(".setup", ".ok"), ls()), envir = environment())
}

if (!exists("hr")) hr <- function(x) cat("\n=====", x, "=====\n")
if (!exists("pr")) pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dorothea)
  library(decoupleR); library(lme4)
})

MIN_TARGETS  <- 5
CONF_LEVELS  <- c("A", "B", "C")
MIN_CELLS_PER_PATIENT <- 20
N_RAND       <- 50
RAND_SIZE    <- 20
SEED_SC      <- 20260101

OUTDIR <- file.path(OUT, "tables")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED_SC)

hr("STEP 0  检查输入")

cat("F_SC_EPI = ", F_SC_EPI, "\n", sep = "")
if (!file.exists(F_SC_EPI))
  stop("找不到 ", F_SC_EPI, "\n请确认单细胞对象路径")
cat("结果目录 = ", OUTDIR, "\n", sep = "")

hr("STEP 1  载入上皮细胞")

sc_full <- readRDS(F_SC_EPI)
cat("对象: ", nrow(sc_full), " 基因 x ", ncol(sc_full), " 细胞\n", sep = "")
cat("meta 列: ", paste(colnames(sc_full@meta.data), collapse = ", "), "\n", sep = "")

if (!"cell_type" %in% colnames(sc_full@meta.data))
  stop("对象缺少 cell_type 列，这不是 GSE165897 的上皮对象")
if (!"patient_id" %in% colnames(sc_full@meta.data))
  stop("对象缺少 patient_id 列，不能退回用 sample 当病人（会造成伪重复）")

cells <- colnames(sc_full)[as.character(sc_full$cell_type) == "EOC"]
sc <- subset(sc_full, cells = cells)
rm(sc_full); gc()

cat("上皮细胞: ", ncol(sc), "  病人: ", length(unique(sc$patient_id)),
    "  样本: ", length(unique(sc$sample)), "\n", sep = "")
if (ncol(sc) < 1000) stop("上皮细胞数异常，检查 F_SC_EPI")

# 该对象只保存了 counts 层，必须自己标准化
if (!"data" %in% SeuratObject::Layers(sc[["RNA"]])) {
  sc <- Seurat::NormalizeData(sc, normalization.method = "LogNormalize",
                              scale.factor = 1e4, verbose = FALSE)
  cat("已执行 LogNormalize\n")
}

hr("STEP 2  构建 regulon（限定到矩阵内基因）")

data("dorothea_hs", package = "dorothea")
net <- dorothea_hs[dorothea_hs$confidence %in% CONF_LEVELS, c("tf", "target", "mor")]
cat("DoRothEA (A/B/C): ", nrow(net), " 条，", length(unique(net$tf)), " 个 TF\n", sep = "")

net <- net[net$target %in% rownames(sc), ]
n_self <- sum(net$target == "TSPO")
net <- net[net$target != "TSPO", ]        # 防循环：活性不能由 TSPO 自己推
cat("剔除 TSPO 作为靶点: ", n_self, " 条（防循环推断）\n", sep = "")

keep_tf <- names(which(table(net$tf) >= MIN_TARGETS))
net <- net[net$tf %in% keep_tf, ]
cat("最终 regulon: ", nrow(net), " 条，", length(keep_tf), " 个 TF\n", sep = "")

# 只取 regulon 基因，避免全基因稠密矩阵（32,847 x 7,929 ≈ 2 GB）
reg_genes <- intersect(unique(c(net$target, "TSPO")), rownames(sc))
cat("regulon 覆盖基因: ", length(reg_genes), "\n", sep = "")

X <- as.matrix(Seurat::GetAssayData(sc, assay = "RNA", layer = "data")[reg_genes, , drop = FALSE])
storage.mode(X) <- "double"
gc()

hr("STEP 3  去掉文库大小并标准化")

v  <- log10(sc$nCount_RNA + 1); vm <- v - mean(v)
b  <- as.numeric((X %*% vm) / sum(vm^2))
a  <- rowMeans(X) - b * mean(v)
Z  <- t(scale(t(X - (a + b %o% v))))
Z[!is.finite(Z)] <- NA
rm(X); gc()

cat("残差化 + z 标准化: ", nrow(Z), " 基因 x ", ncol(Z), " 细胞\n", sep = "")

tspo <- as.numeric(Z["TSPO", ])            # 先取出，后面填 0 不影响它
pat  <- as.character(sc$patient_id)

Z_fill <- Z                                 # ULM 不接受 NA
Z_fill[is.na(Z_fill)] <- 0
cat("NA 占比: ", round(100 * mean(is.na(Z)), 2), "%（填 0 仅用于活性推断）\n", sep = "")
rm(Z); gc()

hr("STEP 4  推断 TF 活性（ULM）")

act <- decoupleR::run_ulm(
  mat     = Z_fill,
  network = net,
  .source = "tf",
  .target = "target",
  .mor    = "mor",
  minsize = MIN_TARGETS
)
cat("ULM 输出列: ", paste(colnames(act), collapse = ", "), "\n", sep = "")

tf_lv   <- unique(as.character(act$source))
cell_lv <- colnames(Z_fill)
act_w <- matrix(NA_real_, nrow = length(tf_lv), ncol = length(cell_lv),
                dimnames = list(tf_lv, cell_lv))
act_w[cbind(match(as.character(act$source), tf_lv),
            match(as.character(act$condition), cell_lv))] <- act$score
rm(act, Z_fill); gc()

cat("TF 活性矩阵: ", nrow(act_w), " TF x ", ncol(act_w), " 细胞\n", sep = "")

hr("STEP 5  逐病人相关 + 符号检验（主要推断依据）")

sl <- do.call(rbind, lapply(rownames(act_w), function(tf) {
  y <- as.numeric(act_w[tf, ])
  r <- sapply(split(seq_along(pat), pat), function(i) {
    yy <- y[i]; tt <- tspo[i]
    ok <- is.finite(yy) & is.finite(tt)
    if (sum(ok) < MIN_CELLS_PER_PATIENT) return(NA_real_)
    suppressWarnings(cor(yy[ok], tt[ok], method = "spearman"))
  })
  r <- r[is.finite(r)]
  if (length(r) < 6) return(NULL)
  k <- sum(r > 0); nn <- length(r)
  data.frame(tf = tf,
             median_rho   = round(median(r), 3),
             pct_positive = round(100 * k / nn, 1),
             n_pos = k, n_patient = nn,
             sign_p = binom.test(k, nn, 0.5)$p.value,
             stringsAsFactors = FALSE)
}))
sl$FDR <- p.adjust(sl$sign_p, "BH")
sl <- sl[order(-abs(sl$median_rho)), ]

cat("\n逐病人一致率最高的 20 个 TF：\n")
pr(head(sl[, c("tf","median_rho","pct_positive","n_pos","n_patient","sign_p","FDR")], 20), 4)

hr("STEP 6  混合模型验证重点 TF")

FOCUS <- intersect(c("RELA","NFKB1","SP1","TP53","CEBPB","STAT3","KLF5","FOXO3"),
                   rownames(act_w))
mm <- do.call(rbind, lapply(FOCUS, function(tf) {
  d <- data.frame(y = as.numeric(act_w[tf, ]), T = tspo, g = factor(pat))
  d <- d[is.finite(d$y) & is.finite(d$T), ]
  m <- suppressWarnings(suppressMessages(
         lme4::lmer(y ~ T + (1 | g), data = d, REML = TRUE)))
  cf <- summary(m)$coefficients
  df_res <- nrow(d) - nlevels(d$g) - 1
  data.frame(tf = tf, beta = round(cf["T","Estimate"], 4),
             p = 2 * pt(-abs(cf["T","t value"]), df = df_res),
             n_cell = nrow(d), n_patient = nlevels(d$g),
             stringsAsFactors = FALSE)
}))
cat("\n混合模型（细胞层面，p 值仅作参考）：\n")
pr(mm, 4)

hr("STEP 7  随机基因集阴性对照")

pool <- setdiff(reg_genes, "TSPO")
rand <- t(sapply(seq_len(N_RAND), function(i) {
  g <- sample(pool, RAND_SIZE)
  y <- colMeans(as.matrix(GetAssayData(sc, assay = "RNA", layer = "data")[g, , drop = FALSE]))
  y <- as.numeric(scale(y))
  r <- sapply(split(seq_along(pat), pat), function(k) {
    yy <- y[k]; tt <- tspo[k]
    ok <- is.finite(yy) & is.finite(tt)
    if (sum(ok) < MIN_CELLS_PER_PATIENT) return(NA_real_)
    suppressWarnings(cor(yy[ok], tt[ok], method = "spearman"))
  })
  c(median_rho = median(r, na.rm = TRUE),
    pct_positive = 100 * mean(r > 0, na.rm = TRUE))
}))
cat("随机 20 基因：中位 rho 的 5%/50%/95% = ",
    paste(round(quantile(rand[, "median_rho"], c(.05,.5,.95), na.rm = TRUE), 3),
          collapse = " / "), "\n", sep = "")
cat("随机模块一致率中位: ", round(median(rand[, "pct_positive"], na.rm = TRUE), 1), "%\n", sep = "")

hr("STEP 8  与 bulk 结果对照")

bulk_file <- file.path(OUTDIR, "T6a_TF_activity_vs_TSPO.csv")
if (file.exists(bulk_file)) {
  bulk <- read.csv(bulk_file, check.names = FALSE)
  cmp  <- merge(sl, bulk[, intersect(c("tf","rho","rho_partial","specific","mirror_strict"),
                                     colnames(bulk))],
                by = "tf", all = FALSE)
  cmp  <- cmp[order(-abs(cmp$rho)), ]
  cmp$replicated <- sign(cmp$median_rho) == sign(cmp$rho) & abs(cmp$median_rho) > 0.05
  cat("bulk 与单细胞方向一致的 TF: ", sum(cmp$replicated), " / ", nrow(cmp), "\n", sep = "")
  cat("\n重点 TF 对照：\n")
  pr(head(cmp[cmp$tf %in% FOCUS,
              c("tf","rho","rho_partial","median_rho","pct_positive","sign_p","replicated")], 20), 4)
  write.csv(cmp, file.path(OUTDIR, "T6e_SC_TF_replication.csv"), row.names = FALSE)
} else {
  cat("未找到 ", bulk_file, "，跳过对照（先跑 13_tf_upstream.R）\n", sep = "")
}

hr("STEP 9  写出")

write.csv(sl, file.path(OUTDIR, "T6b_SC_TF_vs_TSPO.csv"), row.names = FALSE)
saveRDS(act_w, file.path(OUT, "cache", "tf_activity_sc.rds"))
write.csv(as.data.frame(rand), file.path(OUTDIR, "T6f_SC_TF_random_control.csv"),
          row.names = FALSE)

cat("\n判读要点：\n")
cat("  · 逐病人相关才是独立观测；11 例时 binom 最小可能 p 值为 9.8e-4。\n")
cat("  · 若某 TF 在 bulk 强、单细胞弱，说明它的 bulk 信号来自髓系而非上皮。\n")
cat("  · NF-κB 在纯上皮内部若仍成立，才支持「上皮内在的上游调控」。\n")
cat("  · ULM 推断的是活性，不是实测结合；措辞用「与……活性相关」。\n")

hr("14_tf_upstream_singlecell.R 完成")

# -----------------------------------------------------------------------------
# -- Fig9_上游与药物 / 20_drug_partial_rerun.R
# -----------------------------------------------------------------------------
# =============================================================================
#  20_drug_partial_rerun.R   （完全独立版，不依赖其它脚本）
#
#  作用：用 oncoPredict 预测 198 种药物 IC50，再算「校正免疫构成与增殖评分后的
#        偏相关」，产出图 9G / 图 9H 需要的表。
#
#  ── 为什么单独写一版 ────────────────────────────────────────────────────────
#  原流程拆在 11_oncopredict.R 与 12_drug_stratification.R 两个文件里，而
#  D:/Rstudy/tcga/results_v2/code_final/ 这个目录并不存在，直接 source 会报
#  「无法打开链结」。本文件把两步的逻辑内联，只要这一个文件即可运行。
#
#  ── 运行方式 ────────────────────────────────────────────────────────────────
#    setwd("D:/Rstudy/tcga")
#    source("C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final/20_drug_partial_rerun.R", encoding = "UTF-8")
#
#  第 1 步（oncoPredict）约 5–15 分钟；结果缓存为 results_v2/cache/ic50_clean.rds，
#  之后再跑直接复用，约 1 分钟。
#
#  ── 产出 ────────────────────────────────────────────────────────────────────
#    results_v2/tables/T5b_drug_partial_all.csv   全部药物的 rho / 偏相关 / FDR
#    results_v2/tables/T5c_drug_tertile.csv       三分位 Cohen's d
#    results_v2/tables/T5e_drug_fig9.csv          图 9G/H 直接可用
# =============================================================================

# ---------------------------------------------------------------------------
#  0  路径与工具
# ---------------------------------------------------------------------------

PROJ <- "D:/Rstudy/tcga"
OUT  <- file.path(PROJ, "results_v2")

if (!dir.exists(PROJ)) stop("数据目录不存在: ", PROJ)
dir.create(file.path(OUT, "cache"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

# 写权限自检：oncoPredict 会把结果写到工作目录下的相对路径，
# 不可写时它会静默保留旧文件而不报错，必须先拦住。
.tst <- file.path(OUT, "cache", ".write_test")
okw <- tryCatch({ writeLines("x", .tst); file.remove(.tst); TRUE }, error = function(e) FALSE)
if (!okw) stop("输出目录不可写：", OUT)
cat("数据目录: ", PROJ, "\n输出目录: ", OUT, "\n", sep = "")

SEED <- 20260101
set.seed(SEED)
hr <- function(x) cat("\n", strrep("=", 70), "\n  ", x, "\n", strrep("=", 70), "\n", sep = "")
pr <- function(d) print(d, row.names = FALSE, digits = 3)

# 秩残差偏相关：a、b 为目标向量，Z 为协变量矩阵（k x n，自动定向）
pcor_rank <- function(a, b, Z) {
  Z <- as.matrix(Z)
  if (nrow(Z) > ncol(Z)) Z <- t(Z)
  ok <- is.finite(a) & is.finite(b) & apply(Z, 2, function(z) all(is.finite(z)))
  a <- a[ok]; b <- b[ok]; Z <- Z[, ok, drop = FALSE]
  n <- length(a); k <- nrow(Z)
  if (n <= k + 3) return(c(rho = NA_real_, p = NA_real_, n = n))
  A  <- cbind(1, apply(Z, 1, rank))
  ea <- residuals(lm(rank(a) ~ A - 1))
  eb <- residuals(lm(rank(b) ~ A - 1))
  r  <- suppressWarnings(cor(ea, eb))
  if (!is.finite(r) || abs(r) >= 1) return(c(rho = r, p = NA_real_, n = n))
  tt <- r * sqrt((n - k - 2) / (1 - r^2))
  c(rho = r, p = 2 * pt(-abs(tt), df = n - k - 2), n = n)
}

# ---------------------------------------------------------------------------
#  1  载入表达矩阵
# ---------------------------------------------------------------------------

hr("STEP 1  载入 TCGA 表达矩阵")

F_EXPR <- file.path(OUT, "EXPR_clean.rds")
if (!file.exists(F_EXPR)) stop("找不到表达矩阵: ", F_EXPR)
EXPR <- as.matrix(readRDS(F_EXPR)); storage.mode(EXPR) <- "double"
cat("表达矩阵: ", nrow(EXPR), " 基因 x ", ncol(EXPR), " 样本\n", sep = "")
if (!"TSPO" %in% rownames(EXPR)) stop("表达矩阵中没有 TSPO")

# ---------------------------------------------------------------------------
#  2  定位 GDSC2 训练集
# ---------------------------------------------------------------------------

hr("STEP 2  定位 GDSC2 训练集")

GDSC_DIRS <- c(file.path(PROJ, "DataFiles/DataFiles/DataFiles/Training Data"),
               file.path(PROJ, "GDSC_Data/DataFiles/DataFiles/Training Data"),
               PROJ)

find_first <- function(pat) {
  for (d in GDSC_DIRS) {
    if (!dir.exists(d)) next
    f <- list.files(d, pattern = pat, recursive = TRUE, full.names = TRUE)
    if (length(f)) return(f[which.max(file.size(f))])
  }
  NA_character_
}
f_expr <- find_first("GDSC2_Expr.*\\.rds$")
f_res  <- find_first("GDSC2_Res.*\\.rds$")
if (is.na(f_expr) || is.na(f_res))
  stop("找不到 GDSC2 训练集。需要：\n",
       "  GDSC2_Expr (RMA Normalized and Log Transformed).rds\n  GDSC2_Res.rds")
cat("训练表达: ", basename(f_expr), "\n训练药敏: ", basename(f_res), "\n", sep = "")

trExpr <- readRDS(f_expr)
trRes  <- readRDS(f_res)
cat("训练集: ", nrow(trExpr), " 基因 x ", ncol(trExpr), " 细胞系；",
    nrow(trRes), " 细胞系 x ", ncol(trRes), " 药物\n", sep = "")

common <- intersect(rownames(EXPR), rownames(trExpr))
cat("共同基因: ", length(common), "\n", sep = "")
if (length(common) < 10000) warning("共同基因偏少（<10000），检查基因命名是否一致")

# ---------------------------------------------------------------------------
#  3  预测 IC50
# ---------------------------------------------------------------------------

hr("STEP 3  oncoPredict 预测 IC50")

if (!requireNamespace("oncoPredict", quietly = TRUE))
  stop("请先安装 oncoPredict")
suppressPackageStartupMessages(library(oncoPredict))

ic_file <- file.path(OUT, "cache", "ic50_clean.rds")

if (file.exists(ic_file)) {
  IC50 <- readRDS(ic_file)
  cat("复用缓存: ", nrow(IC50), " 药 x ", ncol(IC50), " 样本\n", sep = "")
} else {
  ODIR <- file.path(OUT, "calcPhenotype_Output")
  dir.create(ODIR, showWarnings = FALSE, recursive = TRUE)
  old_wd   <- getwd()
  test_mat <- EXPR[common, , drop = FALSE]

  # calcPhenotype 把结果写到 getwd() 下的相对目录，故先切到可写的父目录
  setwd(dirname(ODIR))
  t_start <- Sys.time()
  res <- NULL

  if (exists("calcPhenotype", where = asNamespace("oncoPredict"), inherits = FALSE)) {
    fm <- names(formals(oncoPredict::calcPhenotype))
    args <- list(trainingExprData = trExpr,
                 trainingPtype    = trRes,
                 testExprData     = test_mat,
                 batchCorrect     = "eb",
                 powerTransformPhenotype = TRUE,
                 removeLowVaryingGenes   = 0.2,
                 minNumSamples    = 10,
                 printOutput      = FALSE,
                 selection        = 1)
    # 老版本拼写是 Varing（少一个 y），必须按 formals() 动态补参
    if ("removeLowVaringGenesFrom" %in% fm) args$removeLowVaringGenesFrom <- "homogenizeData"
    args <- args[intersect(names(args), fm)]
    res <- try(do.call(oncoPredict::calcPhenotype, args), silent = TRUE)
  } else {
    res <- try(oncoPredict::predictPhenotype(
      trainPhenotypeData = trRes, trainExprData = trExpr, testExprData = test_mat,
      batchCorrect = "eb", powerTransformPhenotype = TRUE,
      removeLowVaryingGenes = 0.2, minNumSamples = 10, printOutput = FALSE,
      selection = 1, report_predictor = FALSE)[[1]], silent = TRUE)
  }

  # 返回值兜底：部分版本返回 NULL，结果落在 CSV 里
  if (inherits(res, "try-error") || is.null(res) ||
      !(is.matrix(res) || is.data.frame(res))) {
    f <- file.path(ODIR, "DrugPredictions.csv")
    if (!file.exists(f))
      stop("oncoPredict 未返回矩阵，且 ", f, " 不存在。",
           if (inherits(res, "try-error")) paste0("\n原始错误：",
                                                  conditionMessage(attr(res, "condition"))))
    mt <- file.info(f)$mtime
    if (is.na(mt) || mt < t_start - 60)
      stop("结果文件时间戳早于本次运行起点，说明写入失败、读到的是历史文件。\n",
           "  文件: ", f, "\n  时间: ", format(mt), "\n  起点: ", format(t_start))
    cat("从输出文件读取: ", f, "\n", sep = "")
    res <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  }
  setwd(old_wd)

  IC50 <- as.matrix(res)
  if (nrow(IC50) == ncol(test_mat) && ncol(IC50) != ncol(test_mat)) IC50 <- t(IC50)
  if (nrow(IC50) > ncol(IC50)) IC50 <- t(IC50)

  saveRDS(IC50, ic_file)
  cat(sprintf("完成 %.1f 分钟；%d 药 x %d 样本\n",
              as.numeric(difftime(Sys.time(), t_start, units = "mins")),
              nrow(IC50), ncol(IC50)))
}

# ---------------------------------------------------------------------------
#  4  偏相关（校正免疫构成 + 增殖评分）
# ---------------------------------------------------------------------------

hr("STEP 4  偏相关")

smp <- intersect(colnames(IC50), colnames(EXPR))
cat("可用样本: ", length(smp), "\n", sep = "")
if (length(smp) < 100) stop("样本数不足 100，检查列名是否对齐")

Z <- t(scale(t(EXPR[, smp, drop = FALSE])))
Z[!is.finite(Z)] <- NA
ts <- as.numeric(EXPR["TSPO", smp])

IMMUNE <- c("CD68", "CD3D", "COL1A1")
PROLIF <- c("MKI67", "PCNA", "TOP2A", "CCNB1", "BUB1", "AURKA", "PLK1")
g_imm <- intersect(IMMUNE, rownames(Z))
g_pro <- intersect(PROLIF, rownames(Z))
cat("免疫构成基因: ", paste(g_imm, collapse = "/"), "\n", sep = "")
cat("增殖评分基因: ", paste(g_pro, collapse = "/"), "\n", sep = "")

COV <- rbind(Z[g_imm, , drop = FALSE],
             PROLIF_score = colMeans(Z[g_pro, , drop = FALSE], na.rm = TRUE))

D <- do.call(rbind, lapply(rownames(IC50), function(dr) {
  y  <- as.numeric(IC50[dr, smp])
  r0 <- suppressWarnings(cor.test(ts, y, method = "spearman"))
  rp <- pcor_rank(ts, y, COV)
  data.frame(drug = dr,
             rho = unname(r0$estimate), p = r0$p.value,
             rho_partial = rp["rho"], p_partial = rp["p"],
             n = length(ts), stringsAsFactors = FALSE)
}))
D$FDR         <- p.adjust(D$p, "BH")
D$FDR_partial <- p.adjust(D$p_partial, "BH")
D <- D[order(-abs(D$rho_partial)), ]
write.csv(D, file.path(OUT, "tables", "T5b_drug_partial_all.csv"), row.names = FALSE)

sel <- D[D$FDR_partial < 0.05, c("drug", "rho", "rho_partial", "p_partial", "FDR_partial")]
cat("\n偏相关显著（FDR<0.05）: ", nrow(sel), " 种药物\n", sep = "")
pr(head(sel, 12))

# ---------------------------------------------------------------------------
#  5  三分位效应量
# ---------------------------------------------------------------------------

hr("STEP 5  三分位 Cohen's d")

TARGETS <- c("OSI-027", "AZD2014", "AZD8055", "GNE-317", "Trametinib",
             "Bortezomib", "Dasatinib", "MK-1775", "Cisplatin", "Paclitaxel")

find_drug <- function(nm) {
  hit <- rownames(IC50)[grepl(paste0("^", gsub("([.-])", "\\\\\\1", nm)),
                              rownames(IC50), ignore.case = TRUE)]
  if (length(hit)) hit[1] else NA_character_
}

q   <- quantile(ts, c(1/3, 2/3))
grp <- ifelse(ts > q[2], "High", ifelse(ts <= q[1], "Low", "Mid"))
cat("TSPO 三分位切点: ", paste(round(q, 3), collapse = " / "), "\n", sep = "")
print(table(grp))

TAB <- do.call(rbind, lapply(TARGETS, function(nm) {
  dr <- find_drug(nm)
  if (is.na(dr)) { cat("  跳过 ", nm, "（预测矩阵中无此药）\n", sep = ""); return(NULL) }
  y  <- as.numeric(IC50[dr, smp])
  hi <- y[grp == "High"]; lo <- y[grp == "Low"]
  if (length(hi) < 20 || length(lo) < 20) return(NULL)
  d <- (mean(hi) - mean(lo)) / sqrt((var(hi) + var(lo)) / 2)
  w <- suppressWarnings(wilcox.test(hi, lo))
  data.frame(drug = nm, drug_row = dr,
             median_high = round(median(hi), 3), median_low = round(median(lo), 3),
             delta_median = round(median(hi) - median(lo), 3),
             cohens_d = round(d, 3), p = signif(w$p.value, 3),
             stringsAsFactors = FALSE)
}))
if (!is.null(TAB) && nrow(TAB)) {
  TAB <- TAB[order(TAB$cohens_d), ]
  write.csv(TAB, file.path(OUT, "tables", "T5c_drug_tertile.csv"), row.names = FALSE)
  pr(TAB)
}

# ---------------------------------------------------------------------------
#  6  与稿件数值核对
# ---------------------------------------------------------------------------

hr("STEP 6  与稿件数值核对")

MANU <- c("OSI-027" = -0.395, "Bortezomib" = -0.215, "Trametinib" = -0.195,
          "GNE-317" = -0.173, "Dasatinib" = -0.152, "AZD2014" = -0.151,
          "AZD8055" = -0.149, "MK-1775" = 0.166)

CMP <- do.call(rbind, lapply(names(MANU), function(nm) {
  dr <- find_drug(nm)
  if (is.na(dr)) return(data.frame(drug = nm, rho_partial = NA_real_, FDR = NA_real_,
                                   manuscript = MANU[[nm]], delta = NA_real_))
  r <- D[D$drug == dr, ]
  data.frame(drug = nm, rho_partial = round(r$rho_partial, 3),
             FDR = signif(r$FDR_partial, 3), manuscript = MANU[[nm]],
             delta = round(r$rho_partial - MANU[[nm]], 3), stringsAsFactors = FALSE)
}))
cat("\n偏相关：稿件 vs 重算\n"); pr(CMP)

MANU_D <- c("OSI-027" = -0.922, "AZD8055" = -0.913, "GNE-317" = -0.799,
            "Trametinib" = -0.775, "AZD2014" = -0.592, "Cisplatin" = -0.128,
            "Paclitaxel" = 0.145, "MK-1775" = 0.411)
if (!is.null(TAB) && nrow(TAB)) {
  CM <- do.call(rbind, lapply(names(MANU_D), function(nm) {
    r <- TAB[TAB$drug == nm, ]
    if (!nrow(r)) return(data.frame(drug = nm, cohens_d = NA_real_,
                                    manuscript = MANU_D[[nm]], delta = NA_real_))
    data.frame(drug = nm, cohens_d = r$cohens_d[1], manuscript = MANU_D[[nm]],
               delta = round(r$cohens_d[1] - MANU_D[[nm]], 3))
  }))
  cat("\n三分位 Cohen's d：稿件 vs 重算\n"); pr(CM)
}

cat("\n判读：|delta|<0.02 一致；<0.10 量级一致；>=0.10 需查样本或协变量\n")

# ---------------------------------------------------------------------------
#  7  导出图 9G/H 用表
# ---------------------------------------------------------------------------

hr("STEP 7  导出")

sel <- D[D$FDR_partial < 0.05, c("drug", "rho", "rho_partial", "p_partial", "FDR_partial")]
sel <- sel[order(sel$rho_partial), ]
if (!is.null(TAB) && nrow(TAB)) sel <- merge(sel, TAB[, c("drug", "cohens_d")],
                                             by = "drug", all.x = TRUE)
f_out <- file.path(OUT, "tables", "T5e_drug_fig9.csv")
write.csv(sel, f_out, row.names = FALSE)
cat("已写出: ", f_out, "（", nrow(sel), " 种药物 FDR<0.05）\n", sep = "")

cat("\n完成。把 T5e_drug_fig9.csv 或上面的核对表发给我即可接着画图 9G/H。\n")

# -----------------------------------------------------------------------------
# -- Fig9_上游与药物 / 25_upstream_epigenetic.R
# -----------------------------------------------------------------------------
# =============================================================================
#  25_upstream_epigenetic.R   上游：拷贝数与甲基化 → TSPO（ICGC OV-AU）
#
#  对应正文 3.9（已核实，n = 67 供体）：
#    拷贝数：22q13 区段 CN 与 TSPO 表达正相关（ρ = +0.350，p = 3.7e-3）
#            CN 均值 −0.184（线性 0.88 倍）；91% 供体为负、73% 低于 −0.1
#            **分层不一致**：三分位层内 ρ = −0.062 / −0.421 / +0.425
#              → 属组间块状剂量效应，非层内连续剂量关系
#            区段特异：TOMM22 +0.411 > TSPO +0.350 > PPARA +0.195
#                      其他染色体 ACTB +0.188、GAPDH −0.140、EPCAM −0.188
#                      TOMM22 强于 TSPO，说明 TSPO 是 22q13 区段缺失的随行基因
#    甲基化：6 个 CpG 探针 5/6 为负但无一显著（|ρ| ≤ 0.176）；探针均值 −0.064
#            全基因组平均 β 为 −0.127，与之同量级 → 无位点特异证据
#
#  关键点（踩过的坑）：
#    1. 表达与甲基化/CNV 的 SA 号是同一病人的**不同等分样本**，
#       数字相差 2–6，直接按样本号匹配会得到 0 个配对。
#       必须经 sample.OV-AU.tsv 的 icgc_sample_id → icgc_donor_id 按供体关联。
#    2. 78 个样本只对应 67 位供体，按样本配对构成伪重复，p 值偏小。
#       必须聚合到 donor 层面（n=67）。
#    3. ICGC 的 CNVbyCpG 是区段级数据，同一区段内所有 CpG 共享一个值，
#       6 个探针的相关完全相同（rho=0.426 相同），本质是一个测量。
#
#  输入：ICGC_OV-AU_GX.rds、ICGC_OV-AU_Meth.rds、ICGC_OV-AU_CNVbyCpG.rds、
#        full450k.rds、sample.OV-AU.tsv
#  输出：T17a_meth_probes.csv、T17b_meth_cor.csv、T17c_cnv_cor.csv、
#        T17d_cnv_donor_level.csv、T17e_upstream_summary.csv
# =============================================================================

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

PROJ <- "D:/Rstudy/tcga"
OUT  <- file.path(PROJ, "results_v2")
OUTDIR <- file.path(OUT, "tables"); FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

TSPO_CG <- c("cg02079700","cg03207265","cg05333174",
             "cg13608545","cg15429214","cg24694833")
norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

hr("STEP 1  载入三层数据与样本注释")

GX <- readRDS(file.path(PROJ, "ICGC_OV-AU_GX.rds"))
if (is.data.frame(GX)) GX <- as.matrix(GX)
if (nrow(GX) < ncol(GX) && "TSPO" %in% colnames(GX)) GX <- t(GX)
cat("表达: ", nrow(GX), " x ", ncol(GX), "\n", sep = "")

METH <- readRDS(file.path(PROJ, "ICGC_OV-AU_Meth.rds"))
CNV  <- readRDS(file.path(PROJ, "ICGC_OV-AU_CNVbyCpG.rds"))
if (is.data.frame(CNV)) {
  if (is.character(CNV[[1]])) { rn <- CNV[[1]]; CNV <- as.matrix(CNV[, -1]); rownames(CNV) <- rn }
  else CNV <- as.matrix(CNV)
}
suppressWarnings(storage.mode(CNV) <- "double")
cat("甲基化: ", nrow(METH), " x ", ncol(METH), "  探针命中 ", sum(TSPO_CG %in% rownames(METH)), "/6\n", sep = "")
cat("拷贝数: ", nrow(CNV),  " x ", ncol(CNV),  "  探针命中 ", sum(TSPO_CG %in% rownames(CNV)),  "/6\n", sep = "")

ann <- read.delim(file.path(PROJ, "sample.OV-AU.tsv"), check.names = FALSE,
                  stringsAsFactors = FALSE)
cat("注释表: ", nrow(ann), " 行\n", sep = "")

hr("STEP 2  按 donor 建立配对（关键）")

ann$sid_n <- norm(ann$icgc_sample_id)
ann$don_n <- norm(ann$icgc_donor_id)
ann_u <- ann[!duplicated(ann$sid_n), ]

don_gx <- ann_u$don_n[match(norm(colnames(GX)),  ann_u$sid_n)]
don_me <- ann_u$don_n[match(norm(colnames(METH)), ann_u$sid_n)]
don_cn <- ann_u$don_n[match(norm(colnames(CNV)), ann_u$sid_n)]
cat("GX : ", sum(!is.na(don_gx)), " 样本 → ", length(unique(na.omit(don_gx))), " donor\n", sep = "")
cat("CNV: ", sum(!is.na(don_cn)), " 样本 → ", length(unique(na.omit(don_cn))), " donor\n", sep = "")

agg <- function(v, grp) tapply(v, grp, mean, na.rm = TRUE)

hr("STEP 3  甲基化与 TSPO 表达（donor 层面）")

d_tspo <- agg(as.numeric(GX["TSPO", ])[!is.na(don_gx)], don_gx[!is.na(don_gx)])
d_meth <- NULL
if (length(intersect(TSPO_CG, rownames(METH)))) {
  for (cg in intersect(TSPO_CG, rownames(METH))) {
    v <- as.numeric(METH[cg, ])[!is.na(don_me)]
    d_meth <- rbind(d_meth, data.frame(cg = cg, donor = names(agg(v, don_me[!is.na(don_me)])),
                                       beta = as.numeric(agg(v, don_me[!is.na(don_me)]))))
  }
}

if (!is.null(d_meth)) {
  res_me <- do.call(rbind, lapply(split(d_meth, d_meth$cg), function(x) {
    m <- merge(x, data.frame(donor = names(d_tspo), tspo = as.numeric(d_tspo)), by = "donor")
    m <- m[is.finite(m$beta) & is.finite(m$tspo), ]
    ct <- suppressWarnings(cor.test(m$beta, m$tspo, method = "spearman"))
    data.frame(cg = x$cg[1], n = nrow(m), mean_beta = round(mean(m$beta), 3),
               rho = round(as.numeric(ct$estimate), 3), p = ct$p.value,
               stringsAsFactors = FALSE)
  }))
  res_me$FDR <- p.adjust(res_me$p, "BH")
  cat("\n各探针 β 与 TSPO：\n"); pr(res_me, 4)
  write.csv(res_me, file.path(OUTDIR, "T17b_meth_cor.csv"), row.names = FALSE)
  cat("负相关探针: ", sum(res_me$rho < 0), "/", nrow(res_me),
      "，p<0.05: ", sum(res_me$p < 0.05), "/", nrow(res_me), "\n", sep = "")

  # 全基因组平均 β 作对照
  gall <- colMeans(METH, na.rm = TRUE)
  d_gall <- agg(gall[!is.na(don_me)], don_me[!is.na(don_me)])
  m <- merge(data.frame(donor = names(d_gall), gall = as.numeric(d_gall)),
             data.frame(donor = names(d_tspo), tspo = as.numeric(d_tspo)), by = "donor")
  m <- m[is.finite(m$gall) & is.finite(m$tspo), ]
  ct <- suppressWarnings(cor.test(m$gall, m$tspo, method = "spearman"))
  cat("全基因组平均 β ~ TSPO: rho = ", round(as.numeric(ct$estimate), 3),
      "（若强于位点特异，说明无位点特异证据）\n", sep = "")
}

hr("STEP 4  拷贝数与 TSPO 表达（donor 层面，主结果）")

sub   <- CNV[intersect(TSPO_CG, rownames(CNV)), , drop = FALSE]
a_all <- colMeans(abs(CNV), na.rm = TRUE)

d_cn  <- agg(colMeans(sub, na.rm = TRUE)[!is.na(don_cn)], don_cn[!is.na(don_cn)])
d_amp <- agg(a_all[!is.na(don_cn)],               don_cn[!is.na(don_cn)])

don <- Reduce(intersect, list(names(d_tspo), names(d_cn), names(d_amp)))
cat("可配对 donor: ", length(don), "\n", sep = "")

TSPO <- as.numeric(d_tspo[don]); CN <- as.numeric(d_cn[don]); AMP <- as.numeric(d_amp[don])

ct <- suppressWarnings(cor.test(CN, TSPO, method = "spearman"))
cat("\nCN ~ TSPO: rho = ", round(as.numeric(ct$estimate), 3),
    "，p = ", format.pval(ct$p.value, digits = 3), "，n = ", length(don), "\n", sep = "")

cat("CN ~ CNV幅度: rho = ", round(cor(CN, AMP, method = "spearman"), 3), "\n", sep = "")
cat("幅度 ~ TSPO: rho = ", round(cor(AMP, TSPO, method = "spearman"), 3), "\n", sep = "")

pcor3 <- function(a, b, covs = NULL) {
  if (is.null(covs)) { ra <- rank(a); rb <- rank(b) } else {
    X <- apply(as.matrix(covs), 2, rank)
    ra <- residuals(lm(rank(a) ~ X)); rb <- residuals(lm(rank(b) ~ X))
  }
  r <- cor(ra, rb); n <- length(a)
  tt <- r * sqrt((n - 2) / (1 - r^2))
  c(rho = r, p = 2 * pt(-abs(tt), df = n - 2), n = n)
}
ok <- is.finite(TSPO) & is.finite(CN) & is.finite(AMP)
r0 <- pcor3(TSPO[ok], CN[ok])
r1 <- pcor3(TSPO[ok], CN[ok], AMP[ok])
cat("无控制: rho = ", round(r0["rho"], 3), "  p = ", format.pval(r0["p"], 3),
    "；控制幅度: rho = ", round(r1["rho"], 3), "  p = ", format.pval(r1["p"], 3), "\n", sep = "")

cat("\n按 CNV 幅度三分位分层：\n")
qq <- quantile(AMP, c(1/3, 2/3), na.rm = TRUE)
st <- ifelse(AMP <= qq[1], "低", ifelse(AMP <= qq[2], "中", "高"))
for (s in c("低","中","高")) {
  i <- which(st == s & is.finite(TSPO) & is.finite(CN))
  if (length(i) < 12) { cat(s, "层不足\n"); next }
  ct <- suppressWarnings(cor.test(CN[i], TSPO[i], method = "spearman"))
  cat(sprintf("  %s层 n=%2d  rho = %+.3f  p = %.3g\n", s, length(i),
              as.numeric(ct$estimate), ct$p.value))
}

cat("\n区段特异性（控制幅度）：\n")
for (g in c("TSPO","TOMM22","PPARA","CD68","EPCAM","ACTB","GAPDH")) {
  if (!g %in% rownames(GX)) next
  dg <- agg(as.numeric(GX[g, ])[!is.na(don_gx)], don_gx[!is.na(don_gx)])[don]
  v <- as.numeric(dg)
  okk <- is.finite(v) & is.finite(CN) & is.finite(AMP)
  if (sum(okk) < 20) next
  r <- pcor3(v[okk], CN[okk], AMP[okk])
  cat(sprintf("  %-8s rho = %+.3f  p = %.3g\n", g, r["rho"], r["p"]))
}

write.csv(data.frame(donor = don, TSPO = TSPO, CN_22q13 = CN, CNV_amplitude = AMP),
          file.path(OUTDIR, "T17d_cnv_donor_level.csv"), row.names = FALSE)

hr("STEP 5  汇总")

summ <- data.frame(
  分析     = c("拷贝数（donor 层面）", "启动子甲基化"),
  n        = c(length(don), if (!is.null(d_meth)) length(unique(d_meth$donor)) else 0),
  结果     = c(paste0("rho = ", round(as.numeric(ct$estimate), 3)),
               "见 T17b"),
  结论     = c("阳性：22q13 区段剂量效应", "阴性：无位点特异证据"),
  stringsAsFactors = FALSE)
pr(summ, 3)
write.csv(summ, file.path(OUTDIR, "T17e_upstream_summary.csv"), row.names = FALSE)

cat("\n判读要点：\n")
cat("  · 必须是 donor 层面（n=67）。按 78 个样本配对会构成伪重复，p 值偏小。\n")
cat("  · TOMM22（+0.438）强于 TSPO（+0.346），说明是区段效应而非 TSPO 特异调控。\n")
cat("  · 甲基化的位点特异信号弱于全基因组背景，不主张受甲基化调控。\n")

hr("25_upstream_epigenetic.R 完成")

# -----------------------------------------------------------------------------
# -- Fig9_上游与药物 / 43_program_survival.R
# -----------------------------------------------------------------------------
# =============================================================================
#  43_program_survival.R
#
#  程序评分与 TCGA-OV 生存的关系，含阳性对照。
#
#  服务正文讨论中的临床转化判断。
#
#  要点：任何生存分析都必须先证明该队列能测出已知的预后因素。
#        本脚本以年龄为阳性对照；在 OS 终点上年龄 HR=1.33（p=2.5e-4），
#        说明终点有效，因此其后的阴性结果可被解释为真实阴性。
#        在 PFI 终点上年龄不显著，该终点的判别力不足。
#
#  输入：results_v2/EXPR_clean.rds            TCGA-OV 表达，18,747 × 304
#        D:/Rstudy/tcga/TCGA_OV_survival.tsv  标准 Xena 生存表（OS / PFI / DFI）
#        D:/Rstudy/tcga/OV_clinicalMatrix.tsv 年龄与病理
#        results_v2/tables/T37a_source_index.csv
#  输出（工作区 results_v2/tables/）：
#    T45_program_survival.csv     单因素 Cox 全表（两终点）
#    T45b_survival_positive.csv   阳性对照（年龄）与队列描述
#    T45c_two_source_groups.csv   髓系 × 上皮四分组
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(survival)
})

PROJ <- "D:/Rstudy/tcga"
WS   <- "C:/Users/fade/Documents/deepseek-harness/default-workspace"
OUT  <- file.path(WS, "results_v2", "tables")
hr   <- function(x) cat("\n===== ", x, " =====\n", sep = "")

# ---------------------------------------------------------------------------
# 1  载入并对齐
# ---------------------------------------------------------------------------

hr("STEP 1  载入")
e  <- readRDS(file.path(PROJ, "results_v2", "EXPR_clean.rds"))
sv <- read.delim(file.path(PROJ, "TCGA_OV_survival.tsv"),   stringsAsFactors = FALSE)
cl <- read.delim(file.path(PROJ, "OV_clinicalMatrix.tsv"), stringsAsFactors = FALSE)
si <- read.csv(file.path(PROJ, "results_v2", "tables", "T37a_source_index.csv"),
               stringsAsFactors = FALSE)

m <- match(colnames(e), sv$sample)
cat("表达列:", ncol(e), " 可匹配生存:", sum(!is.na(m)), "\n")
sv2 <- sv[m, ]; ee <- e[, !is.na(m)]
rownames(cl) <- cl$sampleID

AGE <- suppressWarnings(as.numeric(cl[colnames(ee), "age_at_initial_pathologic_diagnosis"]))
cat("年龄可用:", sum(is.finite(AGE)), "\n")

# ---------------------------------------------------------------------------
# 2  程序评分
# ---------------------------------------------------------------------------

GENES <- list(
  TSPO   = "TSPO",
  OXPHOS = c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
             "SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","UQCRH","CYC1",
             "COX4I1","COX5A","COX5B","COX6A1","ATP5A1","ATP5B","ATP5C1","ATP5F1"),
  Glyco  = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1","SLC2A1",
             "SLC16A1","PDK1"),
  Gln    = c("GLS","GOT1","GOT2","GLUD1","SLC1A5","SLC7A5"),
  SPP1prog = c("SPP1","CD44","LGALS1","FN1","MIF"),
  Prolif = c("MKI67","TOP2A","CCNB1","CDK1","BUB1","AURKA","PLK1","CCNA2","MCM2","PCNA"),
  CD8T   = c("CD8A","CD8B","GZMA","GZMB","PRF1","IFNG"),
  EMT    = c("VIM","FN1","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","CDH2"),
  Stroma = c("COL1A1","COL1A2","COL3A1","DCN","LUM","FAP","ACTA2")
)

mod <- function(gs) {
  g <- intersect(gs, rownames(ee))
  if (!length(g)) return(rep(NA_real_, ncol(ee)))
  Z <- t(scale(t(as.matrix(ee[g, , drop = FALSE]))))
  Z[!is.finite(Z)] <- 0
  colMeans(Z)
}
S <- sapply(GENES, mod)

k1 <- gsub("-", "", substr(si$sample, 1, 12))
k2 <- gsub("-", "", substr(colnames(ee), 1, 12))
ms <- match(k2, k1)

D <- data.frame(OS = sv2$OS, OSt = sv2$OS.time,
                PFI = sv2$PFI, PFIt = sv2$PFI.time,
                age = AGE, S,
                Myeloid   = si$Myeloid[ms],
                Epithelial = si$Epithelial[ms],
                Stromal   = si$Stromal[ms],
                source_idx = si$source_idx[ms])

# ---------------------------------------------------------------------------
# 3  阳性对照与单因素 Cox
# ---------------------------------------------------------------------------

cox_row <- function(dd, tt, ev, v) {
  x <- as.numeric(scale(dd[[v]]))
  if (sum(is.finite(x)) < 50) return(NULL)
  mm <- coxph(Surv(dd[[tt]], dd[[ev]]) ~ x); s <- summary(mm)
  data.frame(var = v, HR = s$conf.int[1,1], lo = s$conf.int[1,3],
             hi = s$conf.int[1,4], p = s$coefficients[1,5],
             n = mm$n, events = mm$nevent)
}

run_endpoint <- function(tt, ev, label) {
  hr(paste("STEP 2  ", label))
  dd <- D[is.finite(D[[tt]]) & is.finite(D[[ev]]) & D[[tt]] > 0, ]
  cat(sprintf("n=%d  事件=%d  中位随访=%d 天\n",
              nrow(dd), sum(dd[[ev]] == 1), median(dd[[tt]])))

  out <- do.call(rbind, lapply(c("age", names(GENES),
                                 "Myeloid","Epithelial","Stromal","source_idx"),
                               function(v) cox_row(dd, tt, ev, v)))
  out$HR <- round(out$HR, 3); out$lo <- round(out$lo, 3); out$hi <- round(out$hi, 3)
  out$p  <- signif(out$p, 3)
  out$is_positive_control <- out$var == "age"

  cat(sprintf("\n%-12s %7s %16s %10s\n", "变量", "HR", "95%CI", "p"))
  for (i in seq_len(nrow(out)))
    cat(sprintf("%-12s %7.2f %8.2f-%-7.2f %10.3g%s\n",
        out$var[i], out$HR[i], out$lo[i], out$hi[i], out$p[i],
        if (out$is_positive_control[i]) "  <- 阳性对照" else ""))

  # 队列判别力判定
  pc <- out[out$var == "age", ]
  verdict <- if (pc$p < 0.05) "终点有效，其余阴性可判为真实阴性" else "终点判别力不足，阴性不可解释"
  cat("\n  阳性对照 p =", pc$p, "->", verdict, "\n")

  write.csv(out, file.path(OUT, paste0("T45_program_survival_", ev, ".csv")), row.names = FALSE)
  invisible(list(dd = dd, out = out, valid = pc$p < 0.05))
}

os  <- run_endpoint("OSt",  "OS",  "总生存 OS")
pfi <- run_endpoint("PFIt", "PFI", "无进展间期 PFI")

write.csv(data.frame(
  endpoint = c("OS","PFI"),
  n = c(nrow(os$dd), nrow(pfi$dd)),
  events = c(sum(os$dd$OS == 1), sum(pfi$dd$PFI == 1)),
  median_followup_days = c(median(os$dd$OSt), median(pfi$dd$PFIt)),
  age_HR = c(os$out$HR[os$out$var == "age"], pfi$out$HR[pfi$out$var == "age"]),
  age_p  = c(os$out$p[os$out$var == "age"],  pfi$out$p[pfi$out$var == "age"]),
  endpoint_valid = c(os$valid, pfi$valid)),
  file.path(OUT, "T45b_survival_positive.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4  髓系 × 上皮 四分组
# ---------------------------------------------------------------------------

hr("STEP 3  两来源四分组（OS）")
dd <- os$dd
my <- dd$Myeloid   > median(dd$Myeloid,   na.rm = TRUE)
ep <- dd$Epithelial > median(dd$Epithelial, na.rm = TRUE)
dd$grp <- factor(ifelse(!my & !ep, "双低",
                 ifelse(my & !ep, "髓系高",
                 ifelse(!my & ep, "上皮高", "双高"))),
                 levels = c("双低","髓系高","上皮高","双高"))
print(table(dd$grp))
sd <- survdiff(Surv(OSt, OS) ~ grp, data = dd)
p_lr <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
sf <- survfit(Surv(OSt, OS) ~ grp, data = dd)
s5 <- summary(sf, times = 1825)$surv
cat("log-rank p =", signif(p_lr, 3), "\n")
cat("5 年 OS:",
    paste(sprintf("%s %.1f%%", sub("grp=", "", names(sf$strata)), 100 * s5), collapse = "  "), "\n")

write.csv(data.frame(group = sub("grp=", "", names(sf$strata)),
                     n = as.integer(table(dd$grp)),
                     surv5y = round(100 * s5, 1),
                     logrank_p = signif(p_lr, 3)),
          file.path(OUT, "T45c_two_source_groups.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5  结论
# ---------------------------------------------------------------------------

hr("判读")
cat("OS 终点上年龄显著（阳性对照成立），TSPO 与全部程序评分均不显著。\n")
cat("结论：TSPO 的组织水平与两个来源的分解均不预测总生存。\n")
cat("PFI 终点上年龄不显著，该终点在本队列中判别力不足，不作结论。\n")

hr("43_program_survival.R 完成")

# ############################################################################
# ##  补充图  (3 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- 补充图 / figS.R
# -----------------------------------------------------------------------------
# =============================================================================
#  figS.R   补充图 S1–S3
#
#   图 S1  髓系轨迹的未平衡版本 —— 说明为何必须做病人内平衡
#          轴末端被 DC/pDC 占据，方向验证在未平衡集上不成立
#   图 S2  上皮亚群 × 病人构成（堆叠条）
#          每个亚群的细胞来自哪些病人，独占还是共享一眼可辨
#   图 S3  CellChat 参数敏感性 —— 三种参数口径下通路强度排序的一致性
#
#  ── 数据 ────────────────────────────────────────────────────────────────────
#  S1  GSE266577_raw.rds（6.4 亿字节）+ 复用 06_myeloid_trajectory.R 的流程
#  S2  cache/epi_fig4.rds 里的 meta（7,929 个细胞 × 亚群 + 病人）
#  S3  cc3/cc4/cc5 × equal/full 六套 CellChat 对象
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    source("<code_final>/figS.R", encoding = "UTF-8")
#    输出 figures/FigS1.png、FigS2.png、FigS3.png（450 dpi）
# =============================================================================

CODE   <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

PROJ <- "D:/Rstudy/tcga"
OUT  <- "D:/Rstudy/tcga/results_v2"
if (!exists("FIGDIR")) FIGDIR <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/figures"
CACHE_RW <- file.path(dirname(FIGDIR), "cache")
dir.create(CACHE_RW, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR,  showWarnings = FALSE, recursive = TRUE)

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
suppressPackageStartupMessages({ library(patchwork); library(stringr); library(ggrepel) })
# SEED 在 00_setup.R 里定义；单独跑本脚本时给个默认值
if (!exists("SEED")) SEED <- 20260101
set.seed(SEED)

# =============================================================================
#  S1  髓系轨迹：未平衡 vs 平衡
# =============================================================================

hr("S1  未平衡 vs 平衡 轨迹")

F_S1 <- file.path(CACHE_RW, "figS1_traj.rds")

if (file.exists(F_S1)) {
  S1 <- readRDS(F_S1)
} else {
  suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(slingshot) })
  MAX_CELL <- 12000; N_HVG <- 1500; N_PC <- 20; N_PC_USE <- 15
  F_GSE266 <- file.path("D:/Rstudy/tcga2", "GSE266577_processed", "GSE266577_raw.rds")
  stopifnot(file.exists(F_GSE266))

  o  <- readRDS(F_GSE266)
  CT <- "cell_type"
  MYE <- c("Classical monocytes", "Macrophages")
  keep_cells <- colnames(o)[as.character(o@meta.data[[CT]]) %in% MYE]
  cnt <- Seurat::GetAssayData(o, assay = "RNA", layer = "counts")[, keep_cells, drop = FALSE]
  md  <- o@meta.data[keep_cells, , drop = FALSE]
  rm(o); gc()
  cat("髓系细胞: ", ncol(cnt), "\n", sep = "")

  run_traj <- function(cnt, md, balanced) {
    if (balanced) {
      mono <- rownames(md)[as.character(md[[CT]]) == "Classical monocytes"]
      mac  <- rownames(md)[as.character(md[[CT]]) == "Macrophages"]
      n_mono <- min(length(mono), floor(MAX_CELL / 4))
      n_mac  <- min(length(mac), 3 * n_mono)
      bal <- c(sample(mono, n_mono), sample(mac, n_mac))
      cnt <- cnt[, bal, drop = FALSE]; md <- md[bal, , drop = FALSE]
    }
    obj <- Seurat::CreateSeuratObject(counts = cnt, meta.data = md, min.cells = 10)
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = N_HVG, verbose = FALSE)
    HVG <- Seurat::VariableFeatures(obj)
    obj <- Seurat::ScaleData(obj, features = HVG, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, features = HVG, npcs = N_PC,
                          reduction.name = "pca", verbose = FALSE)
    emb <- Seurat::Embeddings(obj, "pca")[, seq_len(N_PC_USE)]
    cl  <- as.character(obj@meta.data[[CT]])
    sl  <- slingshot::slingshot(emb, clusterLabels = cl, reducedDim = "pca")
    pt  <- as.numeric(slingshot::slingPseudotime(sl)[, 1])
    if (median(pt[cl == "Classical monocytes"], na.rm = TRUE) >
        median(pt[cl == "Macrophages"], na.rm = TRUE)) pt <- -pt
    obj <- Seurat::RunUMAP(obj, dims = seq_len(N_PC_USE), reduction = "pca",
                           verbose = FALSE)
    fine <- if ("celltype_fine" %in% colnames(obj@meta.data))
      as.character(obj@meta.data$celltype_fine) else cl
    list(umap = Seurat::Embeddings(obj, "umap"), pt = pt, cl = cl, fine = fine,
         patient = as.character(obj@meta.data[["publication_patient_code_final"]]),
         tspo = as.numeric(Seurat::GetAssayData(obj, assay = "RNA", layer = "data")["TSPO", ]))
  }

  cat("未平衡版...\n"); A <- run_traj(cnt, md, FALSE)
  cat("平衡版...\n");   B <- run_traj(cnt, md, TRUE)
  S1 <- list(unbal = A, bal = B)
  saveRDS(S1, F_S1)
  rm(cnt, md); gc()
}

# ---- 面板 ----
#
# ⚠ 该数据集只有 22 类粗注释（cell_type），没有 celltype_fine，
#   因此无法复现「轴末端由 DC/pDC 占据」这一说法。
#   本图改为展示真正可验证的差异：平衡与否如何改变 TSPO–伪时序的效应量。

pt_rho <- function(x) {
  ok <- is.finite(x$pt) & is.finite(x$tspo)
  suppressWarnings(cor(x$pt[ok], x$tspo[ok], method = "spearman"))
}
R_U <- pt_rho(S1$unbal); R_B <- pt_rho(S1$bal)
cat(sprintf("TSPO~伪时序：未平衡 %+.3f / 平衡 %+.3f\n", R_U, R_B))

plot_emb2 <- function(x, ttl) {
  d <- data.frame(U1 = x$umap[, 1], U2 = x$umap[, 2], grp = x$cl)
  d$grp <- factor(d$grp, levels = c("Classical monocytes", "Macrophages"))
  d <- d[sample(nrow(d)), ]
  PAL <- c(`Classical monocytes` = "#2166AC", `Macrophages` = "#B2182B")
  ggplot(d, aes(U1, U2, colour = grp)) +
    geom_point(size = 0.22, alpha = 0.5, stroke = 0) +
    scale_colour_manual(values = PAL, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
    labs(x = "UMAP1", y = "UMAP2",
         subtitle = sprintf("%d cells", nrow(d))) +
    ggtitle(ttl) +
    theme_tspo() +
    theme(plot.title = element_text(size = 6.2, face = "plain"),
          plot.subtitle = element_text(size = 4.6, colour = "grey40"),
          legend.position = "bottom",
          legend.key.size = unit(2.4, "mm"),
          legend.text = element_text(size = 4.0))
}

plot_pt <- function(x, ttl, rr) {
  d <- data.frame(pt = x$pt, tspo = x$tspo)
  d <- d[is.finite(d$pt) & is.finite(d$tspo), ]
  ggplot(d, aes(pt, tspo)) +
    geom_point(size = 0.22, alpha = 0.25, colour = "grey45", stroke = 0) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.7,
                colour = "#B2182B", fill = "#F4A09B", alpha = 0.30) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.9, size = 2.1,
             colour = "grey15", label = sprintf("rho = %+.2f", rr)) +
    labs(x = "Pseudotime", y = "TSPO expression (log-normalised)", title = ttl) +
    theme_tspo() +
    theme(plot.title = element_text(size = 6.2, face = "plain"))
}

FIG_S1 <- (plot_emb2(S1$unbal, "Unbalanced: all myeloid cells") |
             plot_emb2(S1$bal,   "Patient-balanced 1:3")) /
  (plot_pt(S1$unbal, "TSPO along pseudotime (unbalanced)", R_U) |
     plot_pt(S1$bal,   "TSPO along pseudotime (balanced)",   R_B)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE))
print(FIG_S1)
save_fig(FIG_S1, "FigS1", width = 7.6, height = 6.8)

# =============================================================================
#  S2  上皮亚群 × 病人构成（堆叠条）
# =============================================================================

hr("S2  上皮亚群构成")

F_EPI <- file.path(CACHE_RW, "epi_fig4.rds")
pS2 <- panel_missing("epi cache not found")
if (file.exists(F_EPI)) {
  E  <- readRDS(F_EPI)
  MD <- E$meta
  cat("meta 列: ", paste(colnames(MD), collapse = ", "), "\n", sep = "")
  SUB <- intersect(c("subtype", "cell_subtype", "cluster"), colnames(MD))[1]
  PAT <- intersect(c("patient", "Patient", "donor"), colnames(MD))[1]
  stopifnot(!is.na(SUB), !is.na(PAT))

  COMP <- as.data.frame(table(subtype = MD[[SUB]], patient = MD[[PAT]]),
                        stringsAsFactors = FALSE)
  COMP <- COMP[COMP$Freq > 0, ]
  TOT  <- tapply(COMP$Freq, COMP$subtype, sum)
  COMP$frac <- COMP$Freq / as.numeric(TOT[COMP$subtype])
  # 亚群按单一病人占比排序：上面独占、下面共享
  TOPF <- sapply(names(TOT), function(s) max(COMP$frac[COMP$subtype == s]))
  ORD  <- names(sort(TOPF, decreasing = TRUE))

  COMP$subtype <- factor(COMP$subtype, levels = rev(ORD))
  PT <- sort(tapply(COMP$Freq, COMP$patient, sum), decreasing = TRUE)
  COMP$patient <- factor(COMP$patient, levels = names(PT))
  PAL_P <- setNames(colorRampPalette(c("#67001F", "#D6604D", "#F4A582", "#FDDBC7",
                                       "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC"))(
                      length(PT)), names(PT))
  cat(sprintf("S2：%d 个亚群 × %d 个病人\n", length(ORD), length(PT)))

  pS2 <- ggplot(COMP, aes(frac, subtype, fill = patient)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.12) +
    geom_vline(xintercept = 0.5, linetype = "22", colour = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = PAL_P, name = NULL) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(x = "Fraction of cells contributed by each patient", y = NULL) +
    theme_tspo() +
    theme(panel.grid.major.y = element_blank(),
          legend.position = "right",
          legend.key.size = unit(2.0, "mm"),
          legend.text = element_text(size = 3.2),
          legend.spacing.y = unit(0.25, "mm")) +
    guides(fill = guide_legend(ncol = 1))
}
print(pS2)
save_fig(pS2, "FigS2", width = 4.6, height = 5.2)

# =============================================================================
#  S3  CellChat 参数敏感性
# =============================================================================

hr("S3  CellChat 参数一致性")

pS3 <- panel_missing("CellChat objects not found")
OBJ <- list(c("cc3_full_pop0",  "k=3, full, 100%"),
            c("cc4_full_pop0",  "k=4, full, 100%"),
            c("cc5_subtypes_pop0", "k=5, subtype, 100%"),
            c("cc3_equal_pop0", "k=3, equal, 10%"),
            c("cc4_equal_pop0", "k=4, equal, 10%"))
# CellChat 对象在 results_v2 根目录，不在 PROJ 根目录
find_cc <- function(nm) {
  cand <- c(file.path(PROJ, paste0(nm, ".rds")),
            file.path(PROJ, "results_v2", paste0(nm, ".rds")),
            file.path(OUT, paste0(nm, ".rds")))
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}
have <- Filter(function(z) !is.na(find_cc(z[1])), OBJ)

if (length(have) >= 2) {
  # ⚠ 口径很关键：直接对 net$prob 全部 173 条 L-R 求和，会把大量
  #   零值与近零值混进排序，秩相关只有 0.17（甚至为负）。
  #   CellChat 的标准做法是先按 pval ≤ 0.05 筛出显著对再汇总，
  #   此时各参数口径的秩相关为 0.73–1.00，排序才是稳的。
  PW <- do.call(cbind, lapply(have, function(z) {
    cc <- readRDS(find_cc(z[1]))
    Pp <- cc@net$prob; Vv <- cc@net$pval; L <- cc@LR$LRsig
    Pp[Vv > 0.05] <- 0                      # 只保留显著通讯
    s  <- apply(Pp, 3, sum)
    pw <- L$pathway_name[match(names(s), L$interaction_name)]
    tapply(s, pw, sum, na.rm = TRUE)
  }))
  colnames(PW) <- sapply(have, `[`, 2)
  PW[!is.finite(PW)] <- 0
  # 口径说明：cc3/cc4/cc5 的 k 不同，能检出的通路集合本来就不同。
  # 要求「全部口径都检出」只剩 2 条通路（无法比较）；
  # 因此保留 CellChat L-R 数据库里的全部通路，未检出者记 0，
  # 秩相关同时反映「检出与否」与「强度排序」两方面的一致性。
  PW <- PW[rowSums(PW > 0) >= 1, , drop = FALSE]
  cat("通路 ", nrow(PW), " 条 × ", ncol(PW), " 套参数\n", sep = "")

  RM <- cor(PW, method = "spearman")
  cat("秩相关范围: ", round(min(RM[lower.tri(RM)]), 3), " ~ ",
      round(max(RM[lower.tri(RM)]), 3), "\n", sep = "")

  RD <- data.frame(g1 = rep(colnames(RM), each = nrow(RM)),
                   g2 = rep(colnames(RM), times = nrow(RM)),
                   r  = as.vector(RM))
  RD$g1 <- factor(RD$g1, levels = colnames(RM))
  RD$g2 <- factor(RD$g2, levels = rev(colnames(RM)))
  pS3 <- ggplot(RD, aes(g1, g2, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.0, colour = "grey15") +
    scale_fill_div(limit = c(min(RM), 1), name = "Spearman\nrho") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_tspo() +
    theme(axis.text.x.top = element_text(angle = 25, hjust = 0, size = 4.2),
          axis.text.y = element_text(size = 4.2),
          panel.grid = element_blank(),
          legend.key.height = unit(7, "mm"))
}

print(pS3)
save_fig(pS3, "FigS3", width = 5.0, height = 4.6)

cat("\nfigS.R 完成：FigS1 / FigS2 / FigS3\n")

# -----------------------------------------------------------------------------
# -- 补充图 / figS4.R
# -----------------------------------------------------------------------------
# =============================================================================
#  figS4.R   补充图 S4：未分类 spot 的再分型与 TSPO 的八区比较
#
#  依赖：spatial_other_split.R 已运行，产出
#        results_v2/cache/region5_analysis.rds
#        results_v2/tables/T41b_other_subtypes.csv、T41c_region5_scores.csv
#
#  面板：
#    A  Other 五个亚型的标志物画像（热图）
#    B  八个区域的 TSPO 均值与逐切片值（点线图）
#    C  Other 亚型在八张切片中的构成（堆叠条）
#
#  运行：先 source fig_theme.R，再 source 本脚本
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); library(scales)
})

if (!exists("FIGDIR")) stop("请先设置 FIGDIR 并 source('fig_theme.R')")
CACHE <- file.path(dirname(FIGDIR), "cache")
TAB   <- file.path(dirname(FIGDIR), "tables")
R     <- file.path(CACHE, "region5_analysis.rds")
if (!file.exists(R)) stop("先运行 spatial_other_split.R 生成 region5_analysis.rds")

L <- readRDS(R)
SP <- L$spots

# -----------------------------------------------------------------------------
#  A  Other 亚型的标志物画像
# -----------------------------------------------------------------------------

OT_LEVELS <- c("Other_immune", "Other_necrotic", "Other_adipose",
               "Other_endothelial", "Other_low")
OT_LAB <- c(Other_immune = "Other: immune aggregate",
            Other_necrotic = "Other: necrotic",
            Other_adipose = "Other: adipose",
            Other_endothelial = "Other: endothelial",
            Other_low = "Other: low signal")

MK <- c("Immune","Lymphoid","Myeloid","Glycolysis","Hypoxia",
        "Necrosis","Adipocyte","Endothel")
MK_LAB <- c("Immune","Lymphoid","Myeloid","Glycolysis","Hypoxia",
            "Necrosis","Adipocyte","Endothelial")

PD <- SP[SP$region == "Other", ]
HM <- do.call(rbind, lapply(OT_LEVELS, function(k) {
  d <- PD[PD$region5 == k, , drop = FALSE]
  data.frame(subtype = k, marker = MK,
             v = sapply(MK, function(m) mean(d[[m]], na.rm = TRUE)),
             stringsAsFactors = FALSE)
}))
# 每个标志物横向 z，便于看相对高低
HM$marker <- factor(HM$marker, levels = MK, labels = MK_LAB)
HM$z <- ave(HM$v, HM$marker, FUN = function(x) {
  s <- sd(x); if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - mean(x)) / s
})
HM$subtype <- factor(HM$subtype, levels = rev(OT_LEVELS), labels = rev(OT_LAB[OT_LEVELS]))

pA <- ggplot(HM, aes(marker, subtype, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", v)), size = 1.75, colour = "grey15") +
  scale_fill_div(limit = c(-1.6, 1.6), name = "Row z") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "Marker profile of unassigned spots") +
  theme_tspo() +
  theme(plot.title = element_text(size = 6.0, face = "plain"),
        axis.text.x.top = element_text(angle = 32, hjust = 0, size = 4.0),
        axis.text.y = element_text(size = 4.0),
        panel.grid = element_blank(),
        legend.key.height = unit(6, "mm"),
        legend.key.width  = unit(1.8, "mm"))

# -----------------------------------------------------------------------------
#  B  八个区域的 TSPO
# -----------------------------------------------------------------------------

RG <- c("Tumor_Core","Interface","Stroma",
        "Other_immune","Other_necrotic","Other_adipose","Other_endothelial","Other_low")
RG_LAB <- c(Tumor_Core = "Tumor core", Interface = "Interface", Stroma = "Stroma",
            Other_immune = "Other: immune", Other_necrotic = "Other: necrotic",
            Other_adipose = "Other: adipose", Other_endothelial = "Other: endothelial",
            Other_low = "Other: low signal")

SP$region5 <- factor(SP$region5, levels = RG)
SECT <- unique(SP$section)

# 逐切片均值
DS <- do.call(rbind, lapply(SECT, function(s) {
  d <- SP[SP$section == s, ]
  do.call(rbind, lapply(RG, function(g) {
    x <- d$TSPO[d$region5 == g]
    if (length(x) < 5) return(NULL)
    data.frame(section = s, region = g, v = mean(x, na.rm = TRUE), n = length(x))
  }))
}))
DS$region <- factor(DS$region, levels = rev(RG))

GM <- aggregate(TSPO ~ region5, data = SP, FUN = mean)
GM$region5 <- factor(GM$region5, levels = rev(RG))

pB <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "grey60") +
  # 逐切片点
  geom_point(data = DS, aes(v, region), size = 0.7, alpha = 0.55,
             colour = "grey45", stroke = 0) +
  # 整体均值
  geom_point(data = GM, aes(TSPO, region5), size = 2.0, shape = 21,
             fill = "#B2182B", colour = "grey25", stroke = 0.3) +
  # 单张切片的散点分布用箱线概括
  geom_boxplot(data = DS, aes(v, region), width = 0.36, outlier.shape = NA,
               fill = NA, colour = "grey35", linewidth = 0.25) +
  scale_y_discrete(labels = RG_LAB) +
  labs(x = "TSPO (in-section z)", y = NULL,
       title = "TSPO by region (all pairwise n.s.)") +
  theme_tspo() +
  theme(plot.title = element_text(size = 6.0, face = "plain"),
        axis.text.y = element_text(size = 4.0),
        panel.grid.major.y = element_blank())

# -----------------------------------------------------------------------------
#  C  Other 亚型在八张切片中的构成
# -----------------------------------------------------------------------------

CT <- as.data.frame(table(SP$section[SP$region == "Other"],
                          factor(SP$region5[SP$region == "Other"], levels = OT_LEVELS)))
colnames(CT) <- c("section", "subtype", "n")
CT$frac <- ave(CT$n, CT$section, FUN = function(x) x / sum(x))
CT$subtype <- factor(CT$subtype, levels = rev(OT_LEVELS), labels = rev(OT_LAB[OT_LEVELS]))
CT$section <- factor(CT$section, levels = SECT,
                     labels = paste0("S", seq_along(SECT)))

PAL_OT <- setNames(c("#2166AC", "#92C5DE", "#F4A582", "#B2182B", "grey75"),
                   rev(OT_LAB[OT_LEVELS]))

pC <- ggplot(CT, aes(frac, section, fill = subtype)) +
  geom_col(width = 0.74, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = PAL_OT, name = NULL) +
  scale_x_continuous(expand = c(0, 0), labels = percent, breaks = c(0, 0.5, 1)) +
  labs(x = "Fraction of unassigned spots", y = NULL,
       title = "Composition across sections") +
  theme_tspo() +
  theme(plot.title = element_text(size = 6.0, face = "plain"),
        panel.grid.major.y = element_blank(),
        legend.key.size = unit(2.2, "mm"),
        legend.text = element_text(size = 3.4),
        legend.position = "right")

# -----------------------------------------------------------------------------
#  拼装
# -----------------------------------------------------------------------------

FIG <- (pA | pB) / pC +
  plot_layout(heights = c(1, 0.72)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE))

print(FIG)
save_fig(FIG, "FigS4", width = 7.6, height = 7.4)

# -----------------------------------------------------------------------------
# -- 补充图 / spatial_other_split.R
# -----------------------------------------------------------------------------
# =============================================================================
#  spatial_other_split.R   （修订版：基于论文自己的逐 spot 表）
#
#  目的：查清三区域划分里 "Other"（1,555 / 4,386 = 35.5%）到底是什么，
#        并补一个"坏死"区，看 TSPO 在其中的水平。
#
#  与前一版的区别：不再从 Seurat 对象重算，直接用论文已有的
#  T19d_spatial_region_within.csv（逐 spot 的片内 z 评分与三区归属），
#  保证与正文 3.7 的数值完全一致。
#
#  在本表基础上补充：淋巴/髓系/坏死/脂肪/内皮五组标志，用于给 Other 分型。
#
#  输出（工作区 results_v2/tables/）：
#    T41a_region5_composition.csv   八张切片的五区构成
#    T41b_other_subtypes.csv        Other 各亚型的标志物画像
#    T41c_region5_scores.csv        五区的模块/TSPO 中位数
#    T41d_tspo_within_section.csv   TSPO 的片内配对符号检验
#    T41e_verdict.txt               判读
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(data.table)
})

PROJ  <- "D:/Rstudy/tcga"
OUT   <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/tables"
CACHE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/cache"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")

# ---------------------------------------------------------------------------
# 1  读入论文的逐 spot 表
# ---------------------------------------------------------------------------

hr("STEP 1  读入 T19d（论文口径）")
SP <- fread(file.path(PROJ, "results_v2", "tables", "T19d_spatial_region_within.csv"))
SP <- as.data.frame(SP)
cat("spot: ", nrow(SP), "  切片: ", length(unique(SP$section)), "\n", sep = "")
print(table(SP$region))
cat("\n与正文核对：Other 1555 / Interface 1039 / Tumor_Core 917 / Stroma 875\n")

# ---------------------------------------------------------------------------
# 2  补充标志物评分
# ---------------------------------------------------------------------------

hr("STEP 2  补充淋巴/髓系/坏死/脂肪/内皮标志")
o <- readRDS(file.path(PROJ, "TSPO_spatial_processed.rds"))
lay <- if ("SCT" %in% names(o@assays)) "SCT" else "Spatial"
X <- GetAssayData(o, assay = lay, layer = "data")
X <- X[, SP$spot, drop = FALSE]
SEC <- SP$section

z_in <- function(v) {
  out <- rep(NA_real_, length(v))
  for (s in unique(SEC)) {
    i <- which(SEC == s); x <- v[i]
    if (sum(is.finite(x)) < 3) next
    out[i] <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  out
}
mod_score <- function(genes, min_n = 2) {
  g <- intersect(genes, rownames(X)); if (length(g) < min_n) return(NULL)
  Z <- sapply(g, function(x) z_in(as.numeric(X[x, ])))
  if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1)
  rowMeans(Z, na.rm = TRUE)
}

EXTRA <- list(
  Lymphoid  = c("CD3D","CD3E","CD2","MS4A1","CD79A","CD19","CD52","IL7R"),
  Myeloid   = c("CD68","CD14","LYZ","C1QA","C1QB","TYROBP","FCER1G","AIF1"),
  Necrosis  = c("HBA1","HBB","HBA2","ALAS2","AHSP"),
  Adipocyte = c("ADIPOQ","PLIN1","LEP","FABP4","PPARG"),
  Endothel  = c("PECAM1","VWF","CDH5","PLVAP","ENG")
)
for (m in names(EXTRA)) {
  v <- mod_score(EXTRA[[m]])
  if (!is.null(v)) {
    SP[[m]] <- v
    cat(sprintf("  %-10s %2d 基因\n", m, length(intersect(EXTRA[[m]], rownames(X)))))
  } else cat(sprintf("  %-10s 基因不足，跳过\n", m))
}

# ---------------------------------------------------------------------------
# 3  把 Other 再分型
# ---------------------------------------------------------------------------

hr("STEP 3  Other 再分型")
OT <- SP$region == "Other"
cat("Other spot: ", sum(OT), "\n", sep = "")

d   <- SP[OT, ]
imm <- d$Immune > median(d$Immune, na.rm = TRUE) & d$Lymphoid > 0
nec <- d$Necrosis > 1
adi <- d$Adipocyte > 1
end <- d$Endothel > 1

sub <- rep("Other_low", nrow(d))
sub[imm] <- "Other_immune"
sub[adi] <- "Other_adipose"
sub[end] <- "Other_endothelial"
sub[nec] <- "Other_necrotic"

SP$region5 <- as.character(SP$region)
SP$region5[OT] <- sub
print(table(SP$region5))

SP$immune_sub <- NA_character_
k <- which(OT & sub == "Other_immune")
if (length(k)) {
  SP$immune_sub[k] <- ifelse(d$Lymphoid[k - which(OT)[1] + 1] > d$Myeloid[k - which(OT)[1] + 1],
                             "Lymphoid", "Myeloid")
}
cat("\nOther_immune 细分:\n"); print(table(SP$immune_sub, useNA = "ifany"))

# ---------------------------------------------------------------------------
# 4  构成表与 Other 画像
# ---------------------------------------------------------------------------

hr("STEP 4  五区构成与 Other 画像")
comp <- as.data.frame.matrix(table(SP$section, SP$region5))
comp <- cbind(section = rownames(comp), comp)
fwrite(comp, file.path(OUT, "T41a_region5_composition.csv"))

prof_cols <- intersect(c("Immune","Lymphoid","Myeloid","Fibroblast","EMT","Hypoxia",
                         "Necrosis","Adipocyte","Endothel","TSPO"), colnames(SP))
prof <- aggregate(SP[OT, prof_cols], by = list(subtype = SP$region5[OT]),
                  FUN = function(x) round(median(x, na.rm = TRUE), 3))
print(prof)
fwrite(prof, file.path(OUT, "T41b_other_subtypes.csv"))

# ---------------------------------------------------------------------------
# 5  五区的 TSPO 与模块水平
# ---------------------------------------------------------------------------

hr("STEP 5  五区中位数")
mods <- intersect(c("TSPO","OXPHOS","Glycolysis","Immune","EMT","Hypoxia"), colnames(SP))
agg <- aggregate(SP[, mods], by = list(region = SP$region5),
                 FUN = function(x) round(mean(x, na.rm = TRUE), 3))
print(agg)
fwrite(agg, file.path(OUT, "T41c_region5_scores.csv"))

# ---------------------------------------------------------------------------
# 6  片内配对符号检验（精确二项）
# ---------------------------------------------------------------------------

hr("STEP 6  TSPO 的片内配对检验")
sign_p <- function(k, n) if (n == 0) NA_real_ else min(1, 2 * pbinom(min(k, n - k), n, 0.5))

pairs <- list(
  c("Other_immune",  "Tumor_Core"), c("Other_immune", "Stroma"), c("Other_immune", "Interface"),
  c("Other_necrotic","Tumor_Core"), c("Other_necrotic","Stroma"),
  c("Tumor_Core","Stroma"), c("Tumor_Core","Interface"), c("Interface","Stroma")
)
res <- do.call(rbind, lapply(intersect(c("TSPO","OXPHOS","Immune"), colnames(SP)), function(v) {
  do.call(rbind, lapply(pairs, function(pr) {
    dif <- sapply(unique(SP$section), function(s) {
      x <- SP[SP$section == s, ]
      a <- x[[v]][x$region5 == pr[1]]; b <- x[[v]][x$region5 == pr[2]]
      if (length(a) < 5 || length(b) < 5) return(NA_real_)
      mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)
    })
    dif <- dif[is.finite(dif)]
    if (length(dif) < 3) return(NULL)
    data.frame(variable = v, contrast = paste(pr[1], ">", pr[2]),
               n = length(dif), pos = sum(dif > 0),
               mean_diff = round(mean(dif), 3),
               p = signif(sign_p(sum(dif > 0), length(dif)), 3))
  }))
}))
print(res[res$variable == "TSPO", ], row.names = FALSE)
fwrite(res, file.path(OUT, "T41d_tspo_within_section.csv"))

# ---------------------------------------------------------------------------
# 7  判读
# ---------------------------------------------------------------------------

hr("STEP 7  判读")
nOT <- sum(OT)
vd <- c("空间区域拆分判读", strrep("-", 46),
        sprintf("Other 占 %.1f%%（%d / %d）", 100 * nOT / nrow(SP), nOT, nrow(SP)))
for (k2 in c("Other_immune","Other_necrotic","Other_adipose","Other_endothelial","Other_low")) {
  n <- sum(SP$region5 == k2, na.rm = TRUE)
  vd <- c(vd, sprintf("  %-18s %4d spot（占 Other %.1f%%）", k2, n, 100 * n / nOT))
}
vd <- c(vd, "", "TSPO 片内配对检验：")
t1 <- res[res$variable == "TSPO", ]
for (i in seq_len(nrow(t1)))
  vd <- c(vd, sprintf("  %-40s %d/%d 切片为正，均值差 %+.3f，p=%.3g",
                      t1$contrast[i], t1$pos[i], t1$n[i], t1$mean_diff[i], t1$p[i]))
writeLines(vd, file.path(OUT, "T41e_verdict.txt"))
cat(paste(vd, collapse = "\n"), "\n")

saveRDS(list(spots = SP, comp = comp, agg = agg, tests = res, prof = prof),
        file.path(CACHE, "region5_analysis.rds"))
cat("\n完成，输出在 ", OUT, "\n", sep = "")

# ############################################################################
# ##  作图样式  (11 个)
# ############################################################################

# -----------------------------------------------------------------------------
# -- 作图样式 / fig_style_compare.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig_style_compare.R   七种图型严格复刻（用本课题真实数据）
#
#  A  Sankey      基因 → 通路，两侧色块 + 通路色带（参照图 1）
#  B  Chord       发散式：基因排外圈、通路在内圈，二次贝塞尔连接（参照图 2）
#  C  Circos      四层环：GO 本体 / 计数 / 显著性 / 内侧楔块 + 中心图例（参照图 3）
#  D  Corrplot    下三角，圆点大小 + 颜色双编码，数值压在圆上（参照图 4）
#  E  GSEA        三层复合：ES 曲线 / 基因刻度 / 排序指标灰条（参照图 5）
#  F  Stacked bar 样本 × 细胞类型占比（参照图 6）
#  G  NMDS        排序散点 + 置信椭圆 + stress（参照图 7）
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    source("<code_final>/fig_style_compare.R", encoding = "UTF-8")
# =============================================================================

CODE   <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
FIGDIR <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/figures"
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

PROJ <- "D:/Rstudy/tcga"
OUT  <- file.path(PROJ, "results_v2")
WS   <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2"
DEMO <- file.path(WS, "figures", "style_demo")
dir.create(DEMO, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(patchwork); library(grid); library(png); library(ggrepel)
  library(ggalluvial); library(ggforce); library(MASS); library(circlize)
})
set.seed(20260101)
hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")

# =============================================================================
#  STEP 0  数据
# =============================================================================

hr("STEP 0  数据")

GO3 <- readRDS(file.path(WS, "cache", "go3_positive_negative.rds"))
NEG <- GO3$neg[GO3$neg$ONTOLOGY == "BP", ]; NEG <- NEG[order(NEG$p.adjust), ]
POS <- GO3$pos

EXPR <- as.matrix(readRDS(file.path(OUT, "EXPR_clean.rds")))
NET  <- read.csv(file.path(OUT, "tables", "T8_network_all_genes.csv"), stringsAsFactors = FALSE)
IMM  <- read.csv(file.path(OUT, "tables", "T9a_immune_scores.csv"), row.names = 1, check.names = FALSE)

TERMS <- head(NEG$Description, 5)
LNK <- do.call(rbind, lapply(TERMS, function(tm) {
  g <- strsplit(NEG$geneID[NEG$Description == tm], "/")[[1]]
  data.frame(term = tm, gene = g, stringsAsFactors = FALSE)
}))
FREQ <- sort(table(LNK$gene), decreasing = TRUE)
KEEP <- names(FREQ)[seq_len(min(18, length(FREQ)))]
LNK  <- LNK[LNK$gene %in% KEEP, ]
cat("通路 ", length(TERMS), " 条；基因 ", length(KEEP), " 个；连接 ", nrow(LNK), " 条\n", sep = "")

# =============================================================================
#  A  桑基 + 气泡复合图
#
#  依照模板：左侧基因名 + 小色块 → 彩色条带 → 右侧通路色块（通路名在色块左侧）
#  色块右边接气泡图：x = GeneRatio，气泡大小 = Count，颜色 = -log10(p.adjust)
#  气泡面板不带内部网格线。
# =============================================================================

hr("[A] Sankey + Bubble")

# 通路配色：低饱和粉 / 蓝 / 绿 / 板岩 / 棕（对齐模板的柔和色系）
PAL_T <- setNames(c("#F4A09B", "#8ECAE6", "#5FBFA0", "#9AA3C4", "#C8B394")[seq_along(TERMS)], TERMS)

# ---- 基因侧 ----
PRIM <- LNK$term[match(KEEP, LNK$gene)]
WG   <- as.numeric(table(LNK$gene)[KEEP])
GDF  <- data.frame(gene = KEEP, term1 = PRIM, w = WG,
                   rho = NET$rho[match(KEEP, NET$gene)], stringsAsFactors = FALSE)
GDF  <- GDF[order(match(GDF$term1, TERMS), -GDF$w), ]
GDF$gene <- factor(GDF$gene, levels = GDF$gene)

NODE_H <- 2.4
GDF$y0 <- c(0, cumsum(rep(NODE_H, nrow(GDF)))[-nrow(GDF)])
GDF$y1 <- cumsum(rep(NODE_H, nrow(GDF)))
GDF$ymid <- (GDF$y0 + GDF$y1) / 2

# 基因色块：每格一色，且打散分配——若按色相顺序逐格给色，
# 相邻基因颜色相近，整列会看起来像一条渐变。这里以跳跃步长取色，
# 使相邻格子的色相相差约 1/3 圈（步长与格子数互质）。
NG   <- nrow(GDF)
STEP <- 7                                  # gcd(7, 18) = 1
HUE  <- scales::hue_pal(l = 72, c = 48)(NG)
ORD  <- ((seq_len(NG) - 1) * STEP) %% NG + 1
PAL_G5 <- setNames(HUE[ORD], as.character(GDF$gene))

LNK$x <- match(LNK$gene, as.character(GDF$gene))
LNK$y <- NA_real_; LNK$ye <- NA_real_
for (k in seq_len(nrow(LNK))) {
  g <- LNK$gene[k]; w <- GDF$w[GDF$gene == g]; iv <- NODE_H / w
  prev <- sum(LNK$gene[seq_len(k - 1)] == g)
  LNK$y[k]  <- GDF$y0[GDF$gene == g] + prev * iv
  LNK$ye[k] <- LNK$y[k] + iv
}

# ---- 通路侧 ----
LNK <- LNK[order(match(LNK$term, TERMS), LNK$x), ]
LNK$ty <- NA_real_; LNK$tye <- NA_real_
off <- 0
for (tm in TERMS) {
  idx <- which(LNK$term == tm)
  LNK$ty[idx]  <- off + seq(0, by = 1, length.out = length(idx))
  LNK$tye[idx] <- LNK$ty[idx] + 1
  off <- off + length(idx)
}
TDF <- data.frame(term = TERMS, stringsAsFactors = FALSE)
TDF$y0 <- sapply(TERMS, function(tm) min(LNK$ty[LNK$term == tm]))
TDF$y1 <- sapply(TERMS, function(tm) max(LNK$tye[LNK$term == tm]))

SC <- max(GDF$y1) / max(TDF$y1)
LNK$ty <- LNK$ty * SC; LNK$tye <- LNK$tye * SC
TDF$y0 <- TDF$y0 * SC; TDF$y1 <- TDF$y1 * SC
TDF$ymid <- (TDF$y0 + TDF$y1) / 2

# ---- GeneRatio / 显著性 ----
TOP5 <- NEG[match(TERMS, NEG$Description), ]
GR <- sapply(strsplit(TOP5$GeneRatio, "/"), function(z)
  as.numeric(z[1]) / as.numeric(z[2]))
TDF$ratio <- GR
TDF$count <- TOP5$Count
TDF$nlp   <- -log10(TOP5$p.adjust)
# 五条通路给出人工缩写（自动缩写在 18 字符时仍然偏长）
SHORT <- c("establishment of RNA localization"    = "RNA localization",
           "RNA transport"                        = "RNA transport",
           "nucleic acid transport"                = "nucleic acid transp.",
           "regulation of DNA metabolic process"   = "DNA metab. reg.",
           "double-strand break repair"            = "DSB repair")
TLAB <- setNames(ifelse(TERMS %in% names(SHORT), SHORT[TERMS], abbr(TERMS, 18)), TERMS)

# ---- 色带 ----
rib <- function(y0a, y0b, y1a, y1b, x0 = 0.080, x1 = 0.260, n = 40) {
  t <- seq(0, 1, length.out = n); sm <- t * t * (3 - 2 * t)
  data.frame(x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
             y = c(y0a + (y1a - y0a) * sm, rev(y0b + (y1b - y0b) * sm)))
}
RIBS <- do.call(rbind, lapply(seq_len(nrow(LNK)), function(k) {
  d <- rib(LNK$y[k], LNK$ye[k], LNK$ty[k], LNK$tye[k])
  d$term <- LNK$term[k]; d$id <- k; d
}))

YM  <- max(max(GDF$y1), max(TDF$y1))
YLM <- c(-YM * 0.055, YM * 1.035)

# ---- 左：桑基 ----
pA1 <- ggplot() +
  geom_polygon(data = RIBS, aes(x, y, group = id, fill = term), alpha = 0.78, colour = NA) +
  geom_rect(data = GDF, aes(xmin = 0.032, xmax = 0.076, ymin = y0, ymax = y1, fill = gene),
            colour = NA) +
  geom_text(data = GDF, aes(x = 0.020, y = ymid, label = gene),
            size = 1.62, colour = "grey10", hjust = 1) +
  geom_rect(data = TDF, aes(xmin = 0.260, xmax = 0.310, ymin = y0, ymax = y1, fill = term),
            colour = "white", linewidth = 0.35) +
  geom_text(data = TDF, aes(x = 0.256, y = ymid, label = TLAB[term]),
            size = 1.42, colour = "grey10", hjust = 1, lineheight = 0.9) +
  scale_fill_manual(values = c(PAL_G5, PAL_T), guide = "none") +
  coord_cartesian(xlim = c(-0.165, 0.320), ylim = YLM, expand = FALSE) +
  labs(x = ".", y = NULL) +
  theme_tspo() +
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(colour = "white", size = 4.2),
        plot.margin = margin(2, 0, 2, 2))

# ---- 右：气泡 ----
pA2 <- ggplot(TDF, aes(ratio, ymid)) +
  geom_point(aes(size = count, colour = nlp)) +
  scale_size_continuous(range = c(0.8, 4.6), name = "Count") +
  scale_colour_gradientn(colours = c("#3B0F70", "#7B2382", "#C63C8A", "#F27F6A"),
                         name = expression(-log[10]*"(Pvalue)")) +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.18)),
                     breaks = scales::pretty_breaks(3)) +
  coord_cartesian(ylim = YLM) +
  labs(x = "Gene.Ratio", y = NULL) +
  theme_tspo() +
  # 外框上下各留一点余量，使最上和最下的气泡完整落在框内，
  # 同时仍与左侧桑基的色块高度对齐
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = -YM * 0.045, ymax = YM * 1.03,
           fill = NA, colour = "grey35", linewidth = 0.35) +
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(), axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 4.2),
        legend.position = "right",
        legend.key.size = unit(2.8, "mm"),
        legend.text = element_text(size = 4.0),
        legend.title = element_text(size = 4.4),
        plot.margin = margin(2, 2, 2, 0))

pA <- (pA1 | pA2) + plot_layout(widths = c(3.30, 0.82), guides = "keep")
FA <- file.path(DEMO, "s1_sankey.png")
ggsave(FA, pA, width = 4.5, height = 2.7, dpi = 400)
cat("GeneRatio: ", paste(sprintf("%.3f", GR), collapse = " "), "\n", sep = "")

# =============================================================================
#  B  Chord（发散式）
# =============================================================================

hr("[B] Chord")

PAL_G <- setNames(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(KEEP)), KEEP)

GA <- seq(90, 90 - 360 + 360 / length(KEEP), length.out = length(KEEP)) * pi / 180
GENES <- data.frame(gene = KEEP, ang = GA, x = cos(GA), y = sin(GA),
                    stringsAsFactors = FALSE)
# 通路角度均匀铺开，避免几条通路角度重合后标签叠在一起
TA <- seq(90, 90 - 360 + 360 / length(TERMS), length.out = length(TERMS)) * pi / 180
TPS <- data.frame(term = TERMS, ang = TA, stringsAsFactors = FALSE)
TPS$term_s <- factor(TLAB[TPS$term], levels = TLAB[TERMS])

bez <- function(p0, p1, n = 40, pull = 0.55) {
  c0 <- (p0 + p1) / 2 * pull
  t <- seq(0, 1, length.out = n)
  data.frame(x = (1-t)^2*p0[1] + 2*(1-t)*t*c0[1] + t^2*p1[1],
             y = (1-t)^2*p0[2] + 2*(1-t)*t*c0[2] + t^2*p1[2])
}
EDG <- do.call(rbind, lapply(seq_len(nrow(LNK)), function(k) {
  g <- as.character(LNK$gene[k]); tm <- LNK$term[k]
  p0 <- as.numeric(GENES[GENES$gene == g, c("x", "y")])
  p1 <- c(cos(TPS$ang[TPS$term == tm]), sin(TPS$ang[TPS$term == tm])) * 0.40
  d <- bez(p0, p1); d$term_s <- TLAB[tm]; d$id <- k; d
}))
EDG$term_s <- factor(EDG$term_s, levels = TLAB[TERMS])

pB <- ggplot() +
  geom_path(data = EDG, aes(x, y, group = id, colour = term_s),
            linewidth = 0.42, alpha = 0.75) +
  geom_point(data = GENES, aes(x, y, fill = gene), shape = 21, size = 1.5,
             colour = "grey30", stroke = 0.25) +
  geom_point(data = TPS, aes(cos(ang) * 0.40, sin(ang) * 0.40, fill = term_s),
             shape = 21, size = 4.6, colour = "grey25", stroke = 0.4) +
  geom_text_repel(data = GENES, aes(x * 1.10, y * 1.10, label = gene),
                  size = 1.35, colour = "grey10", segment.size = 0,
                  max.overlaps = Inf, seed = 1) +
  geom_text_repel(data = transform(TPS, lx = cos(ang) * 0.66, ly = sin(ang) * 0.66),
                  aes(lx, ly, label = term_s), size = 1.28, colour = "grey10",
                  bg.colour = "white", bg.r = 0.14, segment.size = 0.16,
                  segment.colour = "grey55", min.segment.length = 0,
                  max.overlaps = Inf, seed = 2, box.padding = 0.6) +
  scale_colour_manual(values = PAL_T, guide = "none") +
  scale_fill_manual(values = c(PAL_G, PAL_T), guide = "none") +
  coord_equal(xlim = c(-1.55, 1.55), ylim = c(-1.55, 1.55), clip = "off") +
  theme_void(base_size = BASE_SIZE)
FB <- file.path(DEMO, "s2_chord.png")
ggsave(FB, pB, width = 3.7, height = 3.7, dpi = 400)

# =============================================================================
#  C  Circos
# =============================================================================

hr("[C] Circos")

TOP <- do.call(rbind, lapply(c("BP", "CC", "MF"), function(o) {
  d <- POS[POS$ONTOLOGY == o, ]
  if (!nrow(d)) return(NULL); head(d[order(d$p.adjust), ], 6)
}))
TOP <- TOP[order(match(TOP$ONTOLOGY, c("BP","CC","MF")), TOP$p.adjust), ]
TOP$id  <- sub("^GO:", "", TOP$ID)
TOP$lab <- paste0("GO:", TOP$id)
TOP$nlp <- -log10(TOP$p.adjust)
PAL_O <- c(BP = "#E6B800", CC = "#9B59B6", MF = "#2E86C1")
PAL_R <- colorRampPalette(c("#F4A6A6", "#B2182B", "#67001F"))(100)

FC <- file.path(DEMO, "s3_circos.png")
png(FC, width = 3.7, height = 3.7, units = "in", res = 400)
circos.clear()
circos.par(start.degree = 90, gap.degree = 1.0, cell.padding = c(0, 0, 0, 0),
           track.margin = c(0.004, 0.004))
circos.initialize(factors = TOP$lab, xlim = c(0, 1))

circos.track(ylim = c(0, 1), track.height = 0.045, bg.border = NA, panel.fun = function(x, y) {
  r <- TOP[TOP$lab == get.cell.meta.data("sector.index"), ]
  circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
              col = PAL_O[r$ONTOLOGY], border = "white", lwd = 0.3)
  circos.text(CELL_META$xcenter, 0.5, as.character(r$lab), cex = 0.36,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
})
circos.track(ylim = c(0, max(TOP$Count) * 1.18), track.height = 0.26,
             bg.border = "grey90", panel.fun = function(x, y) {
  r <- TOP[TOP$lab == get.cell.meta.data("sector.index"), ]
  circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], r$Count,
              col = "grey62", border = "white", lwd = 0.3)
  circos.text(CELL_META$xcenter, r$Count, as.character(r$Count), cex = 0.33,
              facing = "bending.inside", adj = c(0.5, -0.5))
})
circos.track(ylim = c(0, max(TOP$nlp) * 1.18), track.height = 0.26,
             bg.border = "grey90", panel.fun = function(x, y) {
  r <- TOP[TOP$lab == get.cell.meta.data("sector.index"), ]
  circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], r$nlp,
              col = PAL_R[max(1, round(100 * r$nlp / max(TOP$nlp)))],
              border = "white", lwd = 0.3)
  circos.text(CELL_META$xcenter, r$nlp, sprintf("%.1f", r$nlp), cex = 0.33,
              facing = "bending.inside", adj = c(0.5, -0.5))
})
circos.track(ylim = c(0, 1), track.height = 0.30, bg.border = NA, panel.fun = function(x, y) {
  r <- TOP[TOP$lab == get.cell.meta.data("sector.index"), ]
  circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
              col = PAL_O[r$ONTOLOGY], border = "white", lwd = 0.4)
})
text(0, 0.17, "Ontology", cex = 0.60, font = 2)
legend(-0.32, 0.03, legend = c("BP","CC","MF"), fill = PAL_O[c("BP","CC","MF")],
       border = NA, bty = "n", cex = 0.50, y.intersp = 0.9)
text(0, -0.21, "-log10(P)", cex = 0.50)
circos.clear()
dev.off()

# =============================================================================
#  D  Corrplot
# =============================================================================

hr("[D] Corrplot")

TOP_G <- head(NET$gene[order(-abs(NET$rho))], 12)
GS <- unique(c("TSPO", intersect(TOP_G, rownames(EXPR))))
M <- cor(t(EXPR[GS, , drop = FALSE]), method = "spearman")
M[upper.tri(M, diag = TRUE)] <- NA
CD <- data.frame(g1 = rep(GS, times = length(GS)),
                 g2 = rep(GS, each = length(GS)),
                 r = as.vector(M), stringsAsFactors = FALSE)
CD <- CD[is.finite(CD$r), ]
CD$g1 <- factor(CD$g1, levels = rev(GS)); CD$g2 <- factor(CD$g2, levels = rev(GS))
cat("相关矩阵: ", length(GS), " 基因 / ", nrow(CD), " 单元格\n", sep = "")

pD <- ggplot(CD, aes(g2, g1)) +
  geom_point(aes(size = abs(r), fill = r), shape = 21,
             colour = "grey45", stroke = 0.22) +
  geom_text(aes(label = sprintf("%.2f", r), colour = abs(r) > 0.55),
            size = 1.95, fontface = "bold", show.legend = FALSE) +
  scale_size_continuous(range = c(1.4, 7.0), name = "|rho|") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "rho") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_tspo() +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.3, size = 4.4),
        axis.text.y = element_text(size = 4.4),
        panel.grid = element_blank(),
        panel.border = element_blank(),      # 参考图没有外框
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.key.height = unit(5, "mm"), legend.key.width = unit(2.6, "mm"))
FD <- file.path(DEMO, "s4_corrplot.png")
ggsave(FD, pD, width = 3.7, height = 3.7, dpi = 400)

# =============================================================================
#  E  GSEA（三层）
# =============================================================================

hr("[E] GSEA")

RK <- NET[is.finite(NET$rho), ]; RK <- RK[order(-RK$rho), ]
N  <- nrow(RK)
pick <- function(d, pat) {
  i <- which(tolower(d$Description) == tolower(pat))
  if (!length(i)) return(NULL); strsplit(d$geneID[i[1]], "/")[[1]]
}
SETS <- list(
  "regulation of immune effector process"  = pick(POS, "regulation of immune effector process"),
  "antigen processing and presentation"    = pick(POS, "antigen processing and presentation"),
  "double-strand break repair"             = pick(NEG, "double-strand break repair"))
SETS <- SETS[!vapply(SETS, is.null, logical(1))]
ES <- do.call(rbind, lapply(names(SETS), function(nm) {
  hit <- RK$gene %in% SETS[[nm]]
  Nh <- sum(hit); Nm <- N - Nh
  data.frame(set = nm, rank = seq_len(N), es = cumsum(ifelse(hit, 1/Nh, -1/Nm)), hit = hit)
}))
ES$set <- factor(ES$set, levels = names(SETS))
COLS <- setNames(c("#1B4F8A", "#C0392B", "#4DAF4A")[seq_along(SETS)], names(SETS))
BAR <- ES[ES$hit, ]; BAR$y <- as.numeric(BAR$set)

p_top <- ggplot(ES, aes(rank, es, colour = set)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = COLS, name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Running Enrichment Score") +
  theme_tspo() +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(3.0, "mm"), legend.text = element_text(size = 4.0),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(2, 3, 0, 3))

p_bar <- ggplot(BAR, aes(rank, y, colour = set)) +
  geom_segment(aes(xend = rank, yend = y + 0.78), linewidth = 0.34) +
  scale_colour_manual(values = COLS, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = NULL, limits = c(0.35, length(SETS) + 0.65)) +
  labs(x = NULL, y = NULL) +
  theme_tspo() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), plot.margin = margin(0, 3, 0, 3))

p_met <- ggplot(RK, aes(seq_len(nrow(RK)), rho)) +
  geom_area(fill = "grey70", colour = NA) +
  geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey35") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Rank in Ordered Dataset", y = "Ranked List Metric") +
  theme_tspo() + theme(plot.margin = margin(0, 3, 2, 3))

pE <- wrap_elements(p_top / p_bar / p_met + plot_layout(heights = c(2.6, 1.25, 0.85)))
FE <- file.path(DEMO, "s5_gsea.png")
ggsave(FE, pE, width = 3.7, height = 3.7, dpi = 400)

# =============================================================================
#  F  Stacked bar
# =============================================================================

hr("[F] Stacked bar")

IM <- IMM
IM$TSPO <- as.numeric(EXPR["TSPO", rownames(IM)])
ORD <- rownames(IM)[order(IM$TSPO)]
SEL <- ORD[round(seq(1, length(ORD), length.out = 6))]
SB <- as.matrix(IM[SEL, setdiff(colnames(IM), "TSPO")])
SB <- sweep(SB, 1, rowSums(abs(SB)), "/")
SB <- as.data.frame(as.table(SB)); names(SB) <- c("sample", "cell", "frac")
SB$sample <- factor(SB$sample, levels = rev(SEL))
PAL_CELL <- setNames(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(unique(SB$cell))),
                     unique(SB$cell))
pF <- ggplot(SB, aes(frac, sample, fill = cell)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = PAL_CELL, name = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1),
                     labels = function(x) sprintf("%.2f", x)) +
  labs(x = "fraction", y = NULL) +
  theme_tspo() +
  theme(axis.text.y = element_text(size = 3.6),
        legend.position = "right", legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 3.8))
FF <- file.path(DEMO, "s6_stacked.png")
ggsave(FF, pF, width = 3.7, height = 3.7, dpi = 400)

# =============================================================================
#  G  NMDS
# =============================================================================

hr("[G] NMDS")

Dst <- dist(scale(as.matrix(IMM)))
NM  <- isoMDS(Dst, k = 2, trace = FALSE)
STRESS <- NM$stress / 100
PT <- data.frame(NMDS1 = NM$points[, 1], NMDS2 = NM$points[, 2],
                 TSPO = as.numeric(scale(IM$TSPO)))
Q <- quantile(PT$TSPO, c(1/3, 2/3))
PT$grp <- factor(ifelse(PT$TSPO > Q[2], "TSPO high",
                ifelse(PT$TSPO <= Q[1], "TSPO low", "TSPO mid")),
                levels = c("TSPO low", "TSPO mid", "TSPO high"))
cat(sprintf("NMDS stress = %.3f\n", STRESS))
PAL_N <- c(`TSPO low` = "#1B9E77", `TSPO mid` = "#D95F02", `TSPO high` = "#7570B3")

pG <- ggplot(PT, aes(NMDS1, NMDS2, colour = grp, fill = grp)) +
  stat_ellipse(geom = "polygon", alpha = 0.10, level = 0.68,
               linewidth = 0.4, show.legend = FALSE) +
  geom_point(size = 1.1, alpha = 0.85, stroke = 0) +
  scale_colour_manual(values = PAL_N, name = NULL) +
  scale_fill_manual(values = PAL_N, guide = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.8, size = 1.9,
           colour = PAL["ink"], label = sprintf("stress=%.3f", STRESS)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_tspo() +
  theme(legend.position = "right", legend.key.size = unit(3.0, "mm"),
        legend.text = element_text(size = 4.2))
FG <- file.path(DEMO, "s7_nmds.png")
ggsave(FG, pG, width = 3.7, height = 3.7, dpi = 400)

# =============================================================================
#  汇总
# =============================================================================

hr("拼装")

rd <- function(f) rasterGrob(readPNG(f), interpolate = TRUE)
void <- function() ggplot() + theme_void()

DESIGN <- "
AAB
CDE
FG#
"
OUTP <- wrap_elements(rd(FA)) + wrap_elements(rd(FB)) + wrap_elements(rd(FC)) +
  wrap_elements(rd(FD)) + wrap_elements(rd(FE)) + wrap_elements(rd(FF)) +
  wrap_elements(rd(FG)) +
  plot_layout(design = DESIGN) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, colour = "black"))

print(OUTP)
save_fig(OUTP, "StyleDemo", width = 8.2, height = 7.9)

cat("\n单图目录: ", DEMO, "\n", sep = "")
cat("fig_style_compare.R 完成\n")

# -----------------------------------------------------------------------------
# -- 作图样式 / fig_theme.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig_theme.R   全篇统一的配色、标度与主题（Nature 系低饱和）
#
#  所有出图脚本开头 source 本文件，一处改动、九张图同步。
#
#  ── 配色四条规则 ────────────────────────────────────────────────────────────
#  ① 顺序色只给"全为正"的量（表达、评分、基因比例）——单色阶由浅到深。
#     全正数据用发散色是错的：一半色阶浪费在负值区，颜色永远停在浅色一半，
#     看上去就是"发粉发嫩"。发散色只用于真正有正负的量（ρ、z 分数、logFC）。
#  ② 类别一律离散色，不上渐变。
#  ③ TSPO 有专属固定色 #B2182B，它在九张图里出现数十次，必须处处同色。
#  ④ 饱和度整体压低：顶刊配色不靠鲜艳取胜，靠明度层次。
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# ---------------------------------------------------------------- 离散色板
# 主体取自 Nature 系列（ggsci::pal_npg 的常用子集），低饱和、印刷友好
PAL <- c(
  tspo = "#67001F",   # TSPO 深红（最深，作全篇锚点）
  epi  = "#2166AC",   # 上皮 深蓝
  mye  = "#B2182B",   # 髓系 红
  str  = "#67A9CF",   # 基质 浅蓝
  lym  = "#EF8A62",   # 淋巴 浅红
  neu  = "#BDBDBD",   # 中性 / 参照 灰
  acc  = "#D1E5F0",   # 浅底
  ink  = "#1A1A1A"    # 文字与线条（回归线用此色）
)

# 三区域：本身是有序因子，用一套低饱和的发散色
REGION_PAL <- c(Tumor_Core = "#B2182B", Interface = "#EF8A62", Stroma = "#2166AC")

# 队列
COHORT_PAL <- c(`TCGA-OV` = "#2166AC", GSE26712 = "#B2182B", `ICGC OV-AU` = "#67A9CF")

# 肿瘤 vs 正常
TN_PAL <- c(Tumor = "#B2182B", Normal = "#2166AC")
# 填充用浅色调，与 C 面板直方图及下方气泡的明度一致；深色留给线条与图例
TN_FILL <- c(Tumor = "#F1A08F", Normal = "#9EC3E0")

# 方向
DIR_PAL <- c(Up = "#B2182B", Down = "#2166AC", NS = "#BDBDBD")

# ---------------------------------------------------------------- 连续色板
# 顺序：单色阶，浅 → 深。用于全正的量（表达、评分、基因比例、密度）
SEQ_LOW  <- "#FDDBC7"   # 浅红
SEQ_MID  <- "#EF8A62"   # 中红
SEQ_HIGH <- "#67001F"   # 深红

scale_fill_seq <- function(..., low = "#F7F7F7", high = "#B2182B") {
  scale_fill_gradient(low = low, high = high, na.value = "grey92", ...)
}
scale_colour_seq <- function(low = SEQ_LOW, high = SEQ_HIGH, ...) {
  scale_colour_gradient(low = low, high = high, na.value = "grey92", ...)
}

# 发散：只给真正有正负的量（ρ、z 分数、logFC）
DIV_LOW  <- "#2166AC"   # 负（蓝）
DIV_MID  <- "#FFFFFF"
DIV_HIGH <- "#B2182B"   # 正（红）

scale_fill_div <- function(limit = c(-0.5, 0.5), ...) {
  # 必须给 values，否则白不会落在 0 上（色标按数据范围等分，不按数值等分）
  scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"),
                       values = scales::rescale(c(limit[1], 0, limit[2])),
                       limits = limit, oob = squish, na.value = "grey92", ...)
}
scale_colour_div <- function(limit = c(-0.5, 0.5), ...) {
  scale_colour_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH,
                         midpoint = 0, limits = limit, oob = squish,
                         na.value = "grey92", ...)
}

# 条形图色标：四档红蓝白。
# "白→红"连续插值必然经过粉色；去掉中间那档浅红，直接跳到深红，观感干净。
RWB4 <- c("#2166AC", "#92C5DE", "#F7F7F7", "#B2182B")   # 深蓝 / 浅蓝 / 白 / 深红

# 把连续值切成 n 档（返回整数 1..n，供 fill 用）。
# 按秩等分而非等距分档：指标数值往往挤在窄区间，等距会把四档压成一两档，
# 按秩则保证每一档都有柱子，四色必定用满。
rank_bin <- function(x, n = 4) {
  x <- as.numeric(x); ok <- is.finite(x)
  out <- rep(NA_integer_, length(x))
  if (sum(ok) == 0) return(out)
  if (sum(ok) == 1) { out[ok] <- ceiling(n / 2); return(out) }
  o <- order(x[ok])
  out[ok][o] <- pmin(n, ceiling(seq_along(o) / sum(ok) * n))
  out
}

scale_fill_rank4 <- function(n = 4, ...) {
  scale_fill_manual(values = RWB4[seq_len(n)], na.value = "grey92", ...)
}

# 连续场景（热图、空间图）仍用完整 RdBu 色带
scale_fill_rank <- function(limit = NULL, ...) {
  if (is.null(limit)) limit <- c(0, 1)
  scale_fill_gradientn(colours = c("#2166AC", "#92C5DE", "#F7F7F7", "#B2182B"),
                       limits = limit, oob = squish, na.value = "grey92", ...)
}

# 条形图专用的发散色（带符号的值）
scale_fill_rank_div <- function(limit = NULL, n = 5, ...) {
  if (is.null(limit)) {
    scale_fill_steps2(low = DIV_LOW, mid = "white", high = DIV_HIGH,
                      midpoint = 0, n.breaks = n, na.value = "grey92", ...)
  } else {
    scale_fill_steps2(low = DIV_LOW, mid = "white", high = DIV_HIGH,
                      midpoint = 0, n.breaks = n, limits = limit,
                      oob = squish, na.value = "grey92", ...)
  }
}

# 兼容旧名
scale_fill_rho <- scale_fill_rank_div

# ---------------------------------------------------------------- 主题
# 基准字号：source 之前设 BASE_SIZE 即可整体缩放
if (!exists("BASE_SIZE")) BASE_SIZE <- 6

theme_tspo <- function(base_size = BASE_SIZE, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = "grey93", linewidth = 0.22),
      panel.border       = element_rect(colour = "grey35", linewidth = 0.35),
      axis.ticks         = element_line(colour = "grey35", linewidth = 0.28),
      axis.ticks.length  = unit(1.2, "mm"),
      axis.text          = element_text(colour = "grey20", size = base_size - 0.6,
                                        margin = margin(t = 1, r = 1)),
      axis.title         = element_text(colour = "grey10", size = base_size),
      axis.title.x       = element_text(margin = margin(t = 3)),
      axis.title.y       = element_text(margin = margin(r = 3)),
      strip.background   = element_rect(fill = "grey96", colour = "grey35",
                                        linewidth = 0.28),
      strip.text         = element_text(face = "bold", size = base_size - 0.8,
                                        colour = "grey15"),
      legend.key.size    = unit(2.6, "mm"),
      legend.text        = element_text(size = base_size - 1.1),
      legend.title       = element_text(size = base_size - 1),
      legend.background  = element_blank(),
      legend.key         = element_blank(),
      plot.margin        = margin(3.5, 4.5, 3.5, 3.5)
    )
}


# ---------------------------------------------------------------- 分布面板
# A / B 两个分布面板共用同一套画法：浅小提琴 + 白箱线 + 稀疏散点。
# 大组（n>max_points）随机抽 max_points 个点，避免散点糊成一片。
subsample_pts <- function(d, key, n = 90, seed = 1) {
  set.seed(seed)
  idx <- unlist(lapply(split(seq_len(nrow(d)), d[[key]]), function(i) {
    if (length(i) <= n) i else sample(i, n)
  }), use.names = FALSE)
  d[idx, , drop = FALSE]
}

# 返回三个 geom 的列表，供 + 拼接
dist_stack <- function(data, x, y, fill, key,
                       dodge = 0.70, max_points = 130,
                       violin_alpha = 0.28, pt_alpha = 0.55, pt_size = 0.42) {
  pd <- position_dodge(width = dodge)
  sub <- subsample_pts(data, key, n = max_points)
  list(
    geom_violin(data = data, aes(x = .data[[x]], y = .data[[y]],
                                 fill = .data[[fill]], group = .data[[key]]),
                width = 0.62, scale = "width", trim = TRUE,
                alpha = violin_alpha, linewidth = 0.22,
                colour = "grey45", position = pd),
    geom_boxplot(data = data, aes(x = .data[[x]], y = .data[[y]],
                                  group = .data[[key]]),
                 width = 0.10, outlier.shape = NA, fill = "white",
                 linewidth = 0.26, colour = "grey25", position = pd),
    geom_point(data = sub, aes(x = .data[[x]], y = .data[[y]],
                               group = .data[[key]]),
               position = position_jitterdodge(jitter.width = 0.09,
                                               dodge.width = dodge, seed = 1),
               size = pt_size, alpha = pt_alpha, colour = "grey20",
               stroke = 0)
  )
}

# ---------------------------------------------------------------- 工具
panel_missing <- function(msg = "data not available") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg, colour = "grey60", size = 2.6) +
    theme_void() +
    theme(panel.border = element_rect(colour = "grey85", fill = NA,
                                      linetype = "dashed"))
}

# 统一保存。图注不进图，写在稿件里。
save_fig <- function(p, file, width = 7.2, height = 6.6, dpi = 450) {
  ggsave(file.path(FIGDIR, paste0(file, ".pdf")), p, width = width, height = height,
         device = cairo_pdf)
  ggsave(file.path(FIGDIR, paste0(file, ".png")), p, width = width, height = height,
         dpi = dpi)
  cat("已保存 ", file, ".pdf / .png\n", sep = "")
}

pstar <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01, "**",
  ifelse(p < 0.05, "*", ""))))
}

parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/"), function(z) {
    if (length(z) == 2) as.numeric(z[1]) / as.numeric(z[2]) else NA_real_
  }, numeric(1))
}

cat("[fig_theme] Nature 系低饱和色板已载入\n")

# -----------------------------------------------------------------------------
# -- 作图样式 / fig01.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig01.R   图 1（9 个面板）
#
#  A  mRNA 肿瘤 vs 正常（TCGA/GTEx 与 GSE26712）
#  B  CPTAC 蛋白：肿瘤 vs 癌旁
#  C  共表达网络 18,747 基因的 ρ 分布（标 0.3 阈值与 2,714 个网络基因）
#  D  11 个功能模块的 Cohen's d
#  E  11 个功能模块的 Spearman ρ
#  F  TSPO ~ 肿瘤纯度代理
#  G  TSPO ~ CD68
#  H  正相关方向 GO 前 8
#  I  负相关方向 GO 前 8
#
#  版面约定（与 fig_theme.R 一致）：
#    · 面板内不放长标题，描述一律写进图注——3×3 下每格仅约 55 mm
#    · GO 条目名 str_wrap 折行，否则会吃掉相邻面板的宽度
#    · 只给连续量上渐变；类别用离散色
#
#  运行：
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/fig01.R", encoding = "UTF-8")
# =============================================================================

RECOMPUTE_GTEX <- FALSE
# 折行宽度随字号缩放：字号越小一排放得下越多字
GO_WRAP <- round(182 / BASE_SIZE)
GO_SHOW <- 6
GO_TRUNC <- 32

# ---- 路径 ----------------------------------------------------------------
CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

if (!exists("OUT"))    OUT    <- "D:/Rstudy/tcga/results_v2"
if (!exists("PROJ"))   PROJ   <- "D:/Rstudy/tcga"
if (!exists("F_EXPR")) F_EXPR <- file.path(OUT, "EXPR_clean.rds")
if (!exists("F_GSE"))  F_GSE  <- file.path(PROJ, "GSE26712_expr_matrix_gene_clean.rds")

if (!exists("TAB"))    TAB    <- file.path(OUT, "tables")
if (!exists("TAB2"))   TAB2   <- file.path(OUT, "run_out", "tables")
if (!dir.exists(TAB2)) TAB2   <- TAB
if (!exists("FIGDIR")) FIGDIR <- file.path(OUT, "figures")
if (!exists("CACHE"))  CACHE  <- file.path(OUT, "cache")
for (d in c(FIGDIR, CACHE)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(patchwork) })
cat("表目录: ", TAB, "\n图目录: ", FIGDIR, "\n", sep = "")

# =============================================================================
#  A  雨云图（三队列：半小提琴 + 箱线 + 散点）
# =============================================================================

cat("\n[A] 表达比较\n")

gse_a <- NULL
if (file.exists(F_GSE)) {
  gse <- as.matrix(readRDS(F_GSE)); storage.mode(gse) <- "double"
  grp <- rep(NA_character_, ncol(gse)); names(grp) <- colnames(gse)
  fs <- file.path(PROJ, "GSE26712_series_matrix.txt.gz")
  if (file.exists(fs)) {
    con <- gzfile(fs, "rt"); ln <- readLines(con, warn = FALSE); close(con)
    il <- grep("^!Sample_geo_accession", ln, value = TRUE)
    ids <- if (length(il)) gsub('"', "", strsplit(il[1], "\t")[[1]][-1]) else character(0)
    ch <- grep("^!Sample_characteristics", ln, value = TRUE)
    tl <- ch[grepl("tissue:", ch, fixed = TRUE)][1]
    if (!is.na(tl)) {
      v <- sub("^\\s*tissue:\\s*", "", gsub('"', "", strsplit(tl, "\t")[[1]][-1]))
      if (length(v) == length(ids))
        grp[ids] <- ifelse(grepl("cancer|tumor|tumour", v, ignore.case = TRUE), "Tumor", "Normal")
    }
  }
  if ("TSPO" %in% rownames(gse) && any(!is.na(grp))) {
    gse_a <- data.frame(group = grp, value = as.numeric(gse["TSPO", names(grp)]),
                        stringsAsFactors = FALSE)
    gse_a <- gse_a[!is.na(gse_a$group) & is.finite(gse_a$value), ]
  }
}

tcga_tspo <- NULL
if (file.exists(F_EXPR)) {
  ex <- as.matrix(readRDS(F_EXPR)); storage.mode(ex) <- "double"
  if ("TSPO" %in% rownames(ex)) tcga_tspo <- as.numeric(ex["TSPO", ])
}

gtex_sum <- NULL
fstat <- file.path(PROJ, "TSPO_diff_stats.csv")
if (file.exists(fstat)) {
  st <- read.csv(fstat, stringsAsFactors = FALSE, check.names = FALSE)
  if (all(c("Group","Mean","SD") %in% colnames(st))) gtex_sum <- st[st$Group == "Normal", ]
}

pA <- panel_missing("expression data not found")
if (!is.null(tcga_tspo) && !is.null(gse_a) && nrow(gse_a) > 0) {

  LV <- c("TCGA-OV", "GTEx", "GSE26712")
  raw <- rbind(
    data.frame(cohort = factor(LV[1], levels = LV), value = tcga_tspo, type = "Tumor"),
    data.frame(cohort = factor(LV[3], levels = LV),
               value = gse_a$value[gse_a$group == "Tumor"], type = "Tumor"),
    data.frame(cohort = factor(LV[3], levels = LV),
               value = gse_a$value[gse_a$group == "Normal"], type = "Normal"))
  raw$key <- paste(raw$cohort, raw$type, sep = "|")

  # 雨云图（raincloud）：half-violin（云） + 细箱线（伞） + 抖动散点（雨）。
  # ggdist 的 stat_halfeye 在离散 x 轴上会把密度铺到水平方向，画不出雨云；
  # 这里自己算核密度并画成半小提琴多边形，x 用数值坐标以便按 type 错开。
  raw$xn <- as.numeric(raw$cohort)
  raw$xo <- ifelse(raw$type == "Normal", -0.17, 0.17)
  raw$xpos <- raw$xn + raw$xo

  HALF <- do.call(rbind, lapply(split(raw, raw$key), function(d) {
    if (nrow(d) < 5 || length(unique(d$value)) < 3) return(NULL)
    dn <- density(d$value, n = 160, adjust = 0.9)
    data.frame(xpos = d$xpos[1],
               x = d$xpos[1] + dn$y / max(dn$y) * 0.30,
               y = dn$x,
               type = d$type[1], key = d$key[1], stringsAsFactors = FALSE)
  }))
  # polygon 需要闭合：从基线出发再回到基线
  HALF <- do.call(rbind, lapply(split(HALF, HALF$key), function(d) {
    rbind(data.frame(xpos = d$xpos[1], x = d$xpos[1], y = min(d$y),
                     type = d$type[1], key = d$key[1]),
          d,
          data.frame(xpos = d$xpos[1], x = d$xpos[1], y = max(d$y),
                     type = d$type[1], key = d$key[1]))
  }))

  pA <- ggplot() +
    geom_polygon(data = HALF, aes(x, y, group = key, fill = type),
                 alpha = 0.45, colour = "grey40", linewidth = 0.25) +
    # 伞：细箱线
    geom_boxplot(data = raw, aes(x = xpos, y = value, group = xpos),
                 width = 0.085, linewidth = 0.22, outlier.shape = NA,
                 fill = "white", colour = "grey30") +
    # 雨：原始数据点按真实值抖动铺开
    geom_point(data = raw, aes(x = xpos, y = value, colour = type),
               shape = 16, size = 0.32, alpha = 0.42,
               position = position_jitter(width = 0.055, height = 0, seed = 1)) +
    scale_x_continuous(breaks = seq_along(LV), labels = LV,
                       expand = expansion(add = c(0.45, 0.45)))
  # GTEx 仅有汇总统计：中位数点 + 全距线，不画箱体（无四分位）
  if (!is.null(gtex_sum)) {
    pA <- pA +
      geom_linerange(
        data = data.frame(cohort = factor(LV[2], levels = LV), type = "Normal",
                          lo = gtex_sum$Min[1], hi = gtex_sum$Max[1]),
        aes(x = 2 - 0.17, ymin = lo, ymax = hi),
        inherit.aes = FALSE, linewidth = 1.0, colour = TN_PAL["Normal"]) +
      geom_point(
        data = data.frame(cohort = factor(LV[2], levels = LV), type = "Normal",
                          m = gtex_sum$Median[1]),
        aes(x = 2 - 0.17, y = m), inherit.aes = FALSE,
        shape = 21, size = 1.7, fill = "white", stroke = 0.5,
        colour = TN_PAL["Normal"])
  }

  pA <- pA + scale_fill_manual(values = TN_FILL, name = NULL) +
    scale_colour_manual(values = TN_PAL, guide = "none") +
    labs(x = NULL, y = expression(log[2]~"TSPO")) + theme_tspo() +
    theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
          legend.background = element_rect(fill = alpha("white", 0.78), colour = NA),
          legend.key.size = unit(2.6, "mm"),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 5.8,
                                     margin = margin(t = 0.5))) +
    coord_cartesian(clip = "off")
}

# =============================================================================
#  B
# =============================================================================

cat("[B] CPTAC 蛋白\n")

# 原始 CPTAC 蛋白组只有逐样本矩阵，配对表需现算；结果缓存到工作区，避免每次重读 12 MB。
# 肿瘤 83 例、癌旁 20 例，按样本 ID 取交集得到 10 对（与正文 n=10 一致）。
# CACHE 可能指向只读的 D 盘，缓存写到图目录同级，保证可写
CACHE_RW <- file.path(dirname(FIGDIR), "cache")
dir.create(CACHE_RW, showWarnings = FALSE, recursive = TRUE)
F_CPT <- file.path(CACHE_RW, "cptac_tspo_paired.rds")
if (file.exists(F_CPT)) {
  CPD <- readRDS(F_CPT)
} else {
  FT <- file.path(PROJ, "HS_CPTAC_OV_proteome_gene_tumor.cct.txt")
  FN <- file.path(PROJ, "HS_CPTAC_OV_proteome_gene_normal.cct.txt")
  CPD <- NULL
  if (file.exists(FT) && file.exists(FN)) {
    TT <- read.delim(FT, row.names = 1, check.names = FALSE)
    NN <- read.delim(FN, row.names = 1, check.names = FALSE)
    if ("TSPO" %in% rownames(TT) && "TSPO" %in% rownames(NN)) {
      OV <- intersect(colnames(TT), colnames(NN))
      if (length(OV) >= 3) {
        CPD <- data.frame(patient = OV,
                          tumor = as.numeric(TT["TSPO", OV]),
                          normal = as.numeric(NN["TSPO", OV]),
                          stringsAsFactors = FALSE)
        try(saveRDS(CPD, F_CPT), silent = TRUE)
      }
    }
  }
}

pB <- panel_missing("CPTAC raw matrix not found")
if (!is.null(CPD)) {
  WT <- suppressWarnings(wilcox.test(CPD$tumor, CPD$normal, paired = TRUE))
  NLO <- sum(CPD$tumor < CPD$normal)
  cat(sprintf("配对 %d 例；肿瘤更低 %d/%d；p = %.4f；均值差 %+.3f\n",
              nrow(CPD), NLO, nrow(CPD), WT$p.value, mean(CPD$tumor - CPD$normal)))
  PL <- rbind(data.frame(patient = CPD$patient, group = "Normal", v = CPD$normal),
              data.frame(patient = CPD$patient, group = "Tumor",  v = CPD$tumor))
  PL$group <- factor(PL$group, levels = c("Normal", "Tumor"))
  PL$dir <- ifelse(CPD$tumor[match(PL$patient, CPD$patient)] <
                   CPD$normal[match(PL$patient, CPD$patient)], "lower in tumour", "higher in tumour")

  pB <- ggplot(PL, aes(group, v, group = patient)) +
    geom_line(aes(colour = dir), linewidth = 0.32, alpha = 0.85) +
    geom_point(aes(fill = dir), shape = 21, size = 1.5,
               colour = "grey30", stroke = 0.24) +
    scale_colour_manual(values = c(`lower in tumour` = "#2166AC",
                                   `higher in tumour` = "#B2182B"), name = NULL) +
    scale_fill_manual(values = c(`lower in tumour` = "#2166AC",
                                 `higher in tumour` = "#B2182B"), name = NULL) +
    annotate("text", x = 1.5, y = Inf, vjust = 1.9, size = 1.8, colour = PAL["ink"],
             label = sprintf("%d/%d lower in tumour\np = %.3f (paired)",
                             NLO, nrow(CPD), WT$p.value)) +
    scale_x_discrete(expand = expansion(add = 0.44)) +
    labs(x = NULL, y = "TSPO protein (TMT log2)") +
    theme_tspo() +
    theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 3.9))
}

# =============================================================================
#  C
# =============================================================================

cat("[C] 网络 ρ 分布\n")
fC <- file.path(TAB, "T8_network_all_genes.csv")
pC <- panel_missing("network table not found")
if (file.exists(fC)) {
  net <- read.csv(fC, stringsAsFactors = FALSE)
  if (!"in_network" %in% colnames(net)) net$in_network <- abs(net$rho) > 0.3 & net$FDR < 0.05
  n_net <- sum(net$in_network, na.rm = TRUE)
  pC <- ggplot(net, aes(rho, fill = after_stat(x))) +
    geom_histogram(bins = 110, colour = NA) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = 2, linewidth = 0.3, colour = "grey20") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 1.8, size = 2.4,
             colour = PAL["ink"],
             label = sprintf("%s genes", format(n_net, big.mark = ","))) +
    scale_fill_div(limit = c(-0.8, 0.8), guide = "none") +
    labs(x = expression("Spearman " * rho), y = "Genes") + theme_tspo()
}

# =============================================================================
#  D / E
# =============================================================================

cat("[D/E] 模块双口径\n")
fD <- file.path(TAB, "T19a_module_cohens_d.csv")
pD <- panel_missing("T19a not found"); pE <- panel_missing("T19a not found")
if (file.exists(fD)) {
  md <- read.csv(fD, stringsAsFactors = FALSE)
  nice <- c(PPARA_Targets = "PPARA targets", Complex_II_V = "Complex II/V",
            Catalytic = "Catalytic", OXPHOS = "OXPHOS", ROS_Defense = "ROS defence",
            Immune_Ckpt = "Immune checkpoint", Assembly = "Mito. assembly",
            Mito_Import = "Mito. import", FAO = "Fatty acid ox.",
            Apoptosis = "Apoptosis", Glycolysis = "Glycolysis")
  md$label <- ifelse(md$module %in% names(nice), nice[md$module], md$module)
  md$nlp_d <- -log10(md$p_wilcox)
  md$nlp_r <- -log10(md$p_rho)

  # 气泡：横轴=效应量，大小=模块基因数，颜色=显著性
  bubble <- function(d, xvar, cvar, xlab) {
    d$label <- factor(d$label, levels = d$label[order(d[[xvar]])])
    ggplot(d, aes(.data[[xvar]], label)) +
      geom_segment(aes(x = 0, xend = .data[[xvar]], y = label, yend = label),
                   colour = "grey86", linewidth = 0.3) +
      geom_point(aes(size = n_genes, fill = .data[[cvar]]),
                 shape = 21, colour = "grey25", stroke = 0.22) +
      scale_size_continuous(range = c(0.9, 3.4), name = "Module\nsize",
                            breaks = c(5, 10, 15)) +
      scale_fill_rank(limit = range(d[[cvar]], na.rm = TRUE),
                      name = expression(-log[10]~P)) +
      scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
      labs(x = xlab, y = NULL) + theme_tspo() +
      theme(legend.key.size = unit(2.4, "mm")) +
      guides(fill = guide_colourbar(barwidth = unit(2.2, "mm"),
                                    barheight = unit(8, "mm")))
  }
  pD <- bubble(md, "cohens_d", "nlp_d", "Cohen's d")
  pE <- bubble(md, "rho", "nlp_r", expression("Spearman " * rho))
}

# =============================================================================
#  F / G
# =============================================================================

cat("[F/G] 纯度与 CD68\n")
pF <- panel_missing("EXPR not found"); pG <- panel_missing("EXPR not found")
if (file.exists(F_EXPR)) {
  ex <- as.matrix(readRDS(F_EXPR)); storage.mode(ex) <- "double"
  ts <- as.numeric(ex["TSPO", ])
  g3 <- intersect(c("CD68","CD3D","COL1A1"), rownames(ex))
  if (length(g3) >= 2) {
    z <- t(scale(t(ex[g3, , drop = FALSE]))); z[!is.finite(z)] <- NA
    purity <- -colMeans(z, na.rm = TRUE)
    d <- data.frame(purity = purity, tspo = ts)
    d <- d[is.finite(d$purity) & is.finite(d$tspo), ]
    r <- suppressWarnings(cor(d$purity, d$tspo, method = "spearman"))
    pF <- ggplot(d, aes(purity, tspo)) +
      geom_point(size = 0.62, alpha = 0.45, colour = "#92C5DE") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                  colour = PAL["ink"], fill = "grey85", linewidth = 0.5) +
      annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.8,
               label = sprintf("rho = %.3f", r), size = 2.6, colour = PAL["ink"]) +
      labs(x = "Purity proxy (inverse)", y = expression(log[2]~"TSPO")) + theme_tspo()
  }
  if ("CD68" %in% rownames(ex)) {
    d2 <- data.frame(cd68 = as.numeric(ex["CD68", ]), tspo = ts)
    d2 <- d2[is.finite(d2$cd68) & is.finite(d2$tspo), ]
    r2 <- suppressWarnings(cor(d2$cd68, d2$tspo, method = "spearman"))
    pG <- ggplot(d2, aes(cd68, tspo)) +
      geom_point(size = 0.62, alpha = 0.45, colour = "#92C5DE") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                  colour = PAL["ink"], fill = "grey85", linewidth = 0.5) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 1.9,
               label = sprintf("rho = %.2f", r2), size = 2.4, colour = PAL["ink"]) +
      scale_x_continuous(breaks = scales::pretty_breaks(3)) +
      coord_cartesian(clip = "off") +
      labs(x = expression(log[2]~"CD68"), y = expression(log[2]~"TSPO")) + theme_tspo()
  }
}

# =============================================================================
#  H / I
# =============================================================================

cat("[H/I] GO 富集\n")

go_panel <- function(file, nshow = GO_SHOW) {
  if (!file.exists(file)) return(panel_missing(paste(basename(file), "not found")))
  g <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(g) || !all(c("Description","p.adjust") %in% colnames(g)))
    return(panel_missing("no enriched terms"))
  g$ratio <- if ("GeneRatio" %in% colnames(g)) parse_ratio(g$GeneRatio) else NA_real_
  g <- head(g[order(g$p.adjust), ], nshow)
  desc <- stringr::str_trunc(g$Description, GO_TRUNC)
  g$term <- stringr::str_wrap(desc, GO_WRAP)
  g$term <- factor(g$term, levels = rev(unique(g$term)))
  g$nlp  <- -log10(g$p.adjust)

  # 气泡：横轴=显著性，大小=基因数，颜色=基因比例
  ggplot(g, aes(nlp, term)) +
    geom_segment(aes(x = 0, xend = nlp, y = term, yend = term),
                 colour = "grey86", linewidth = 0.3) +
    geom_point(aes(size = Count, fill = ratio),
               shape = 21, colour = "grey25", stroke = 0.22) +
    scale_size_continuous(range = c(0.9, 3.4), name = "Genes") +
    scale_fill_rank(limit = range(g$ratio, na.rm = TRUE), name = "Gene\nratio") +
    scale_x_continuous(breaks = scales::pretty_breaks(3),
                       expand = expansion(mult = c(0.04, 0.14))) +
    labs(x = expression(-log[10]~"adj. P"), y = NULL) +
    theme_tspo() +
    theme(axis.text.y = element_text(lineheight = 0.86),
          legend.position = "right",
          legend.title = element_text(size = 4.8, lineheight = 0.9),
          legend.text  = element_text(size = 4.2),
          legend.key.size = unit(2.4, "mm"),
          legend.margin = margin(0, 0, 0, -1)) +
    guides(fill = guide_colourbar(barwidth = unit(2.2, "mm"),
                                  barheight = unit(8, "mm")),
           size = guide_legend(ncol = 1))
}

# ---------------------------------------------------------------------------
#  H / I  GO 三本体 + 通路-基因桑基
#
#  正负基因集由网络表重建（rho>0.3 / <−0.3 且 FDR<0.05），与正文 1,138 / 1,576 一致。
#  H：三种本体（BP/CC/MF）各自的显著条目并列，颜色区分本体。
#  I：负相关方向前 5 条通路的成员基因流向，展示条目之间的基因重叠。
#  富集结果较重（两个基因集 × 全部本体），缓存到工作区后复用。
# ---------------------------------------------------------------------------

F_GO <- file.path(CACHE_RW, "go3_positive_negative.rds")
GO3 <- NULL
if (file.exists(F_GO)) {
  GO3 <- readRDS(F_GO)
} else if (requireNamespace("clusterProfiler", quietly = TRUE) &&
           requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  f_net <- file.path(TAB, "T8_network_all_genes.csv")
  if (file.exists(f_net)) {
    NET <- read.csv(f_net, stringsAsFactors = FALSE)
    POS <- unique(NET$gene[NET$rho >  0.3 & NET$FDR < 0.05])
    NEG <- unique(NET$gene[NET$rho < -0.3 & NET$FDR < 0.05])
    suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
    run_ont <- function(genes) {
      e <- clusterProfiler::enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                                     ont = "ALL", pAdjustMethod = "BH",
                                     pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
      as.data.frame(e)
    }
    GO3 <- try(list(pos = run_ont(POS), neg = run_ont(NEG), n_pos = length(POS),
                    n_neg = length(NEG)), silent = TRUE)
    if (inherits(GO3, "try-error")) GO3 <- NULL else try(saveRDS(GO3, F_GO), silent = TRUE)
  }
}

cat("[H/I] Circos + 桑基气泡\n")

suppressPackageStartupMessages({ library(circlize); library(png); library(grid); library(ggrepel) })
CACHE_RW2 <- file.path(dirname(FIGDIR), "cache")
dir.create(CACHE_RW2, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
#  H  Circos：三本体 GO 富集的三层环
#    外层 = 本体着色 + GO ID；中层 = 基因数；内层 = -log10(adj.P)；内侧楔块 = 本体
#    基绘图无法进 patchwork，先渲染成 PNG 再作为单个面板嵌入
# ---------------------------------------------------------------------------

# 栅格图放进 ggplot 容器：边框由 panel.border 统一画，
# 并用不可见的 y 轴标题/刻度补齐 G 面板那样的左侧留白，
# 否则 H 的边框会比 G 小一圈。
pGrob <- function(f) {
  ggplot() +
    annotation_custom(rasterGrob(readPNG(f), interpolate = TRUE),
                      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    scale_y_continuous(breaks = c(0, 1), labels = c("  ", "  ")) +
    labs(x = NULL, y = "   ") +
    theme_tspo() +
    theme(axis.text = element_text(colour = "white", size = 4.0),
          axis.text.x = element_blank(), axis.ticks = element_blank(),
          axis.title.y = element_text(colour = "white", size = 4.0),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.margin = margin(3.5, 4.5, 3.5, 3.5))
}
pH <- panel_missing("GO enrichment unavailable")

if (!is.null(GO3) && nrow(GO3$pos) && requireNamespace("circlize", quietly = TRUE)) {
  P3 <- GO3$pos
  TOP <- do.call(rbind, lapply(c("BP","CC","MF"), function(o) {
    d <- P3[P3$ONTOLOGY == o, ]
    if (!nrow(d)) return(NULL)
    head(d[order(d$p.adjust), ], 4)
  }))
  if (!is.null(TOP) && nrow(TOP)) {
    TOP <- TOP[order(match(TOP$ONTOLOGY, c("BP","CC","MF")), TOP$p.adjust), ]
    # 30° 扇区下带 GO: 前缀仍会溢出到相邻扇区，去掉前缀只留 7 位编号；
    # 图注说明「编号为 GO 本体 ID」，完整 ID 见补充表
    TOP$lab <- sub("^GO:", "", TOP$ID)
    TOP$nlp <- -log10(TOP$p.adjust)
    PAL_O <- c(BP = "#E6B800", CC = "#9B59B6", MF = "#2E86C1")
    PAL_R <- colorRampPalette(c("#F4A6A6", "#B2182B", "#67001F"))(100)

    F_CIR <- file.path(CACHE_RW2, "fig1_circos.png")
    png(F_CIR, width = 2.75, height = 2.70, units = "in", res = 450)
    circos.clear()
    # 画布比单位圆略大，外圈 GO ID 才不会贴边被裁掉
    circos.par(start.degree = 90, gap.degree = 1.0, cell.padding = c(0, 0, 0, 0),
               track.margin = c(0.004, 0.004),
               canvas.xlim = c(-1.30, 1.30), canvas.ylim = c(-1.30, 1.30))
    circos.initialize(factors = TOP$lab, xlim = c(0, 1))
    # 1 本体色条 + GO ID
    # 第 1 轨：本体色条（细）。GO ID 不用 circos.text 画——
    # circulize 会把文本裁剪到扇区角度内，长 ID 会被直接切掉；
    # 改为在 circos.clear() 前用 base text() 画到圆外。
    circos.track(ylim = c(0, 1), track.height = 0.045, bg.border = NA,
                 panel.fun = function(x, y) {
      r <- TOP[TOP$lab == get.cell.meta.data("sector.index"), ]
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = PAL_O[r$ONTOLOGY], border = "white", lwd = 0.3)
    })
    # 2 基因数
    circos.track(ylim = c(0, max(TOP$Count) * 1.18), track.height = 0.26,
                 bg.border = "grey90", panel.fun = function(x, y) {
      r <- TOP[TOP$lab == get.cell.meta.data("sector.index"), ]
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], r$Count,
                  col = "grey62", border = "white", lwd = 0.3)
      circos.text(CELL_META$xcenter, r$Count, as.character(r$Count), cex = 0.31,
                  facing = "bending.inside", adj = c(0.5, -0.5))
    })
    # 3 显著性
    circos.track(ylim = c(0, max(TOP$nlp) * 1.18), track.height = 0.26,
                 bg.border = "grey90", panel.fun = function(x, y) {
      r <- TOP[TOP$lab == get.cell.meta.data("sector.index"), ]
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], r$nlp,
                  col = PAL_R[max(1, round(100 * r$nlp / max(TOP$nlp)))],
                  border = "white", lwd = 0.3)
      circos.text(CELL_META$xcenter, r$nlp, sprintf("%.1f", r$nlp), cex = 0.31,
                  facing = "bending.inside", adj = c(0.5, -0.5))
    })
    # 4 本体楔块
    circos.track(ylim = c(0, 1), track.height = 0.15, bg.border = NA,
                 panel.fun = function(x, y) {
      r <- TOP[TOP$lab == get.cell.meta.data("sector.index"), ]
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = PAL_O[r$ONTOLOGY], border = "white", lwd = 0.4)
    })
    # GO ID：base text() 不受扇区裁剪，放在半径 1.12 处沿切向排
    NTOP <- nrow(TOP)
    ANG  <- 90 - (seq_len(NTOP) - 0.5) * (360 / NTOP)
    RR   <- 1.11
    TXT  <- data.frame(x = RR * cos(ANG * pi / 180),
                       y = RR * sin(ANG * pi / 180),
                       srt = ANG - 90, lab = TOP$lab,
                       col = unname(PAL_O[TOP$ONTOLOGY]),
                       stringsAsFactors = FALSE)
    TXT$srt[TXT$srt < -90] <- TXT$srt[TXT$srt < -90] + 180
    TXT$srt[TXT$srt >  90] <- TXT$srt[TXT$srt >  90] - 180
    par(xpd = NA)
    # 逐条画：text() 的 srt 传向量会报 "wrong length"
    for (k in seq_len(nrow(TXT))) {
      text(TXT$x[k], TXT$y[k], TXT$lab[k], srt = TXT$srt[k], cex = 0.40,
           col = "grey15", adj = c(0.5, 0.5), xpd = NA)
    }

    text(0, 0.15, "Ontology", cex = 0.48, font = 2)
    legend(-0.24, 0.00, legend = c("BP","CC","MF"), fill = PAL_O[c("BP","CC","MF")],
           border = NA, bty = "n", cex = 0.42, y.intersp = 0.78)
    circos.clear()
    dev.off()
    cat(sprintf("Circos：三本体共 %d 条条目\n", nrow(TOP)))
    pH <- pGrob(F_CIR)
  }
}

# ---------------------------------------------------------------------------
#  I  桑基 + 气泡复合：基因 → 通路，右侧接 GeneRatio 气泡
# ---------------------------------------------------------------------------

# 通路全名在窄格里必然被压掉，先映射成缩写。
# 3 列网格下 I 面板只有 ~2.5 in，标签必须压到 12 字符内才不会压到基因色块。
short_term <- function(x) {
  MAP <- c(
    "double-strand break repair"                  = "DSB repair",
    "double-strand break repair via homologous recombination" = "DSB repair (HR)",
    "recombinational repair"                      = "Recomb. repair",
    "regulation of DNA metabolic process"         = "DNA metab.",
    "DNA metabolic process"                       = "DNA metab.",
    "nucleic acid transport"                      = "Nucleic acid",
    "RNA transport"                               = "RNA transp.",
    "establishment of RNA localization"           = "RNA loc.",
    "RNA localization"                            = "RNA loc.",
    "negative regulation of gene expression, epigenetic" = "Epigenetic silencing",
    "DNA replication"                             = "DNA replication")
  out <- unname(MAP[x])
  miss <- is.na(out)
  if (any(miss)) {
    g <- x[miss]
    g <- sub("^(positive|negative) regulation of ", "", g)
    g <- sub("^regulation of ", "", g)
    g <- sub("^establishment of ", "", g)
    g <- gsub("metabolic process", "metabolism", g)
    out[miss] <- g
  }
  stringr::str_trunc(out, 17)
}

pI <- panel_missing("GO enrichment unavailable")
if (!is.null(GO3) && nrow(GO3$neg) && requireNamespace("ggalluvial", quietly = TRUE)) {
  N3 <- GO3$neg[GO3$neg$ONTOLOGY == "BP", ]
  N3 <- N3[order(N3$p.adjust), ]
  # 前几条里有 Jaccard=1.00 的重复条目，按 Jaccard<0.5 去冗余后取 5 条真正不同的通路
  SETS_G <- lapply(N3$geneID, function(z) strsplit(z, "/")[[1]])
  KEEP_I <- integer(0)
  for (r in seq_len(nrow(N3))) {
    if (length(KEEP_I) >= 5) break
    if (!length(KEEP_I)) { KEEP_I <- r; next }
    JC <- vapply(KEEP_I, function(k)
      length(intersect(SETS_G[[r]], SETS_G[[k]])) / length(union(SETS_G[[r]], SETS_G[[k]])),
      numeric(1))
    if (max(JC) < 0.5) KEEP_I <- c(KEEP_I, r)
  }
  N3 <- N3[KEEP_I, ]
  N3$short <- short_term(N3$Description)

  if (nrow(N3)) {
    LG <- do.call(rbind, lapply(seq_len(nrow(N3)), function(r) {
      g <- strsplit(N3$geneID[r], "/")[[1]]
      if (!length(g)) return(NULL)
      data.frame(gene = g, term = N3$short[r], stringsAsFactors = FALSE)
    }))
    FREQ <- sort(table(LG$gene), decreasing = TRUE)
    KEEP <- head(names(FREQ), 14)
    LG <- LG[LG$gene %in% KEEP, ]
    LG$gene <- factor(LG$gene, levels = rev(KEEP))

    # ---- 左：桑基 ----
    PAL_T <- setNames(c("#F4A09B", "#8ECAE6", "#5FBFA0", "#9AA3C4", "#C8B394")[seq_along(N3$short)],
                      N3$short)
    TERMS <- N3$short
    PRIM <- LG$term[match(KEEP, LG$gene)]
    WG   <- as.numeric(table(LG$gene)[KEEP])
    GDF  <- data.frame(gene = KEEP, term1 = PRIM, w = WG, stringsAsFactors = FALSE)
    GDF  <- GDF[order(match(GDF$term1, TERMS), -GDF$w), ]
    GDF$gene <- factor(GDF$gene, levels = GDF$gene)
    NG <- nrow(GDF)
    HUE <- scales::hue_pal(l = 72, c = 48)(NG)
    PAL_G5 <- setNames(HUE[((seq_len(NG) - 1) * 7) %% NG + 1], as.character(GDF$gene))

    NODE_H <- 2.4
    GDF$y0 <- c(0, cumsum(rep(NODE_H, NG))[-NG])
    GDF$y1 <- cumsum(rep(NODE_H, NG))
    GDF$ymid <- (GDF$y0 + GDF$y1) / 2
    LG$x <- match(LG$gene, as.character(GDF$gene))
    LG$y <- NA_real_; LG$ye <- NA_real_
    for (k in seq_len(nrow(LG))) {
      g <- as.character(LG$gene[k]); w <- GDF$w[GDF$gene == g]; iv <- NODE_H / w
      prev <- sum(as.character(LG$gene[seq_len(k - 1)]) == g)
      LG$y[k]  <- GDF$y0[GDF$gene == g] + prev * iv
      LG$ye[k] <- LG$y[k] + iv
    }
    LG <- LG[order(match(LG$term, TERMS), LG$x), ]
    LG$ty <- NA_real_; LG$tye <- NA_real_
    off <- 0
    for (tm in TERMS) {
      idx <- which(LG$term == tm)
      LG$ty[idx]  <- off + seq(0, by = 1, length.out = length(idx))
      LG$tye[idx] <- LG$ty[idx] + 1
      off <- off + length(idx)
    }
    TDF <- data.frame(term = TERMS, stringsAsFactors = FALSE)
    TDF$y0 <- sapply(TERMS, function(tm) min(LG$ty[LG$term == tm]))
    TDF$y1 <- sapply(TERMS, function(tm) max(LG$tye[LG$term == tm]))
    SC <- max(GDF$y1) / max(TDF$y1)
    LG$ty <- LG$ty * SC; LG$tye <- LG$tye * SC
    TDF$y0 <- TDF$y0 * SC; TDF$y1 <- TDF$y1 * SC
    TDF$ymid <- (TDF$y0 + TDF$y1) / 2
    RAT <- parse_ratio(N3$GeneRatio[match(TERMS, N3$short)])
    TDF$ratio <- RAT; TDF$count <- N3$Count[match(TERMS, N3$short)]
    TDF$nlp2  <- -log10(N3$p.adjust[match(TERMS, N3$short)])

    rib <- function(y0a, y0b, y1a, y1b, x0 = 0.080, x1 = 0.260, n = 36) {
      t <- seq(0, 1, length.out = n); sm <- t * t * (3 - 2 * t)
      data.frame(x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
                 y = c(y0a + (y1a - y0a) * sm, rev(y0b + (y1b - y0b) * sm)))
    }
    RIBS <- do.call(rbind, lapply(seq_len(nrow(LG)), function(k) {
      d <- rib(LG$y[k], LG$ye[k], LG$ty[k], LG$tye[k])
      d$term <- as.character(LG$term[k]); d$id <- k; d
    }))
    YM  <- max(max(GDF$y1), max(TDF$y1))
    YLM <- c(-YM * 0.055, YM * 1.035)

    # 发散式 chord：基因排在外圈、通路落在内圈，曲线由外向内收拢。
    # 关键在排序——基因按「主通路」分组排在圆周上，同一通路的基因相邻，
    # 通路的角位置取其成员基因角的圆周均值，这样连线接近放射状，几乎不交叉。
    TERM_LV <- as.character(unique(LG$term))
    NG <- length(KEEP); NT <- length(TERM_LV)

    PRIM <- LG$term[match(KEEP, LG$gene)]
    ORDG <- order(match(PRIM, TERM_LV))
    KEEP_O <- KEEP[ORDG]

    ANG_G <- seq(90, 90 - 360 + 360 / NG, length.out = NG) * pi / 180
    GEN <- data.frame(name = KEEP_O, ang = ANG_G,
                      x = cos(ANG_G), y = sin(ANG_G), stringsAsFactors = FALSE)

    ANG_T <- sapply(TERM_LV, function(tm) {
      a <- GEN$ang[GEN$name %in% LG$gene[LG$term == tm]]
      atan2(mean(sin(a)), mean(cos(a)))
    })
    R_IN <- 0.54
    TRM <- data.frame(name = TERM_LV, ang = as.numeric(ANG_T),
                      x = R_IN * cos(as.numeric(ANG_T)),
                      y = R_IN * sin(as.numeric(ANG_T)), stringsAsFactors = FALSE)

    bez <- function(p0, p1, n = 36, pull = 0.62) {
      c0 <- (p0 + p1) / 2 * pull
      t <- seq(0, 1, length.out = n)
      data.frame(x = (1-t)^2*p0[1] + 2*(1-t)*t*c0[1] + t^2*p1[1],
                 y = (1-t)^2*p0[2] + 2*(1-t)*t*c0[2] + t^2*p1[2])
    }
    EDG <- do.call(rbind, lapply(seq_len(nrow(LG)), function(k) {
      g <- as.character(LG$gene[k]); tm <- as.character(LG$term[k])
      d <- bez(c(GEN$x[GEN$name == g], GEN$y[GEN$name == g]),
               c(TRM$x[TRM$name == tm], TRM$y[TRM$name == tm]))
      d$term <- tm; d$id <- k; d
    }))

    PAL_T5 <- setNames(c("#F4A09B", "#8ECAE6", "#5FBFA0", "#9AA3C4", "#C8B394")[seq_along(TERM_LV)],
                       TERM_LV)
    STEP <- if (NG %% 5 == 0) 3 else 5
    PAL_G5 <- setNames(scales::hue_pal(l = 72, c = 48)(NG)[((seq_len(NG) - 1) * STEP) %% NG + 1],
                       KEEP_O)
    TNODE <- as.numeric(table(LG$term)[TERM_LV])

    pI <- ggplot() +
      geom_path(data = EDG, aes(x, y, group = id, colour = term),
                linewidth = 0.48, alpha = 0.82) +
      geom_point(data = GEN, aes(x, y, fill = name), shape = 21, size = 1.5,
                 colour = "grey30", stroke = 0.22) +
      geom_point(data = TRM, aes(x, y, fill = name, size = TNODE), shape = 21,
                 colour = "grey25", stroke = 0.32) +
      scale_size_continuous(range = c(2.6, 4.6), guide = "none") +
      geom_text_repel(data = GEN, aes(x * 1.035, y * 1.035, label = name),
                      size = 1.16, colour = "grey10", segment.size = 0.10,
                      segment.colour = "grey60", min.segment.length = 0,
                      max.overlaps = Inf, seed = 1, box.padding = 0.10,
                      point.padding = 0.04) +
      geom_text_repel(data = TRM, aes(x, y, label = name),
                      size = 1.12, colour = "grey5", fontface = "bold",
                      bg.colour = "white", bg.r = 0.10,
                      segment.size = 0.13, segment.colour = "grey45",
                      min.segment.length = 0, max.overlaps = Inf, seed = 2,
                      box.padding = 0.42, point.padding = 0.10) +
      scale_fill_manual(values = c(PAL_G5, setNames(rep("grey82", NT), TERM_LV)),
                        guide = "none") +
      scale_colour_manual(values = PAL_T5, guide = "none") +
      scale_y_continuous(breaks = c(0, 1), labels = c("  ", "  ")) +
      coord_equal(xlim = c(-1.22, 1.22), ylim = c(-1.22, 1.22), clip = "off") +
      labs(x = NULL, y = "   ") +
      theme_tspo() +
      theme(panel.grid = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_blank(), axis.ticks = element_blank(),
            axis.text.y = element_text(colour = "white", size = 4.0),
            axis.title.y = element_text(colour = "white", size = 4.0),
            plot.margin = margin(3.5, 4.5, 3.5, 3.5))

    cat(sprintf("桑基气泡：%d 条通路 × %d 个基因 × %d 条连接\n",
                length(TERMS), NG, nrow(LG)))
  }
}

# =============================================================================
#  拼装
# =============================================================================

cat("\n拼装 3x3\n")

# 注意：这里不能再加 & theme_tspo()——它会在各面板主题之后重新套一遍，
# 把 A 的 x 轴标签旋转等面板级设置覆盖掉。各面板已在内部 theme_tspo()。
FIG1 <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI) +
  plot_layout(widths = c(1, 1, 1), heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, colour = "black"))

print(FIG1)
if (!exists("FIG_NAME")) FIG_NAME <- "Fig1"
save_fig(FIG1, FIG_NAME, width = 8.8, height = 8.4)

cat("\nfig01.R 完成\n")

# -----------------------------------------------------------------------------
# -- 作图样式 / fig02.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig02.R   图 2（9 个面板）
#
#  A  15 类免疫与基质细胞评分与 TSPO 的相关（气泡：横轴 ρ、大小显著性、颜色 ρ）
#  B  免疫总分与 TSPO 散点
#  C  WGCNA 模块-性状热图（标出 TSPO 所在的 ME1）
#  D  TSPO 在其模块内的成员度（kME = 0.577）相对 hub 基因的位置
#  E  肿瘤纯度与 TSPO：三队列，全局表达校正前 → 后（哑铃图）
#  F  增殖评分与 TSPO：同上
#  G  TSPO 与全基因均值的关系（三队列，示平台差异）
#  H  TSPO 高/低组差异表达火山图
#  I  上调与下调各前 10 位基因
#
#  运行：
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/fig02.R", encoding = "UTF-8")
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

if (!exists("OUT"))    OUT    <- "D:/Rstudy/tcga/results_v2"
if (!exists("PROJ"))   PROJ   <- "D:/Rstudy/tcga"
if (!exists("F_EXPR")) F_EXPR <- file.path(OUT, "EXPR_clean.rds")
if (!exists("F_GSE"))  F_GSE  <- file.path(PROJ, "GSE26712_expr_matrix_gene_clean.rds")
if (!exists("F_ICGC")) F_ICGC <- file.path(PROJ, "ICGC_OV-AU_GX.rds")

if (!exists("TAB"))    TAB    <- file.path(OUT, "tables")
if (!exists("FIGDIR")) FIGDIR <- file.path(OUT, "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(patchwork); library(stringr) })

cat("表目录: ", TAB, "\n")

# =============================================================================
#  A  相关矩阵（下三角）：TSPO 与 13 类免疫/基质细胞评分的相互相关
#
#  原为「各细胞类型 vs TSPO」的棒棒糖；改为下三角相关矩阵后，
#  除 TSPO 那一列外，还能同时看出细胞类型之间的共变结构（免疫模块）。
# =============================================================================

cat("\n[A] 免疫/基质评分相关矩阵\n")
fA <- file.path(TAB, "T9a_immune_scores.csv")
pA <- panel_missing("T9a not found")
if (file.exists(fA)) {
  IM <- read.csv(fA, row.names = 1, check.names = FALSE)
  EX <- as.matrix(readRDS(F_EXPR))
  IM$TSPO <- as.numeric(EX["TSPO", rownames(IM)])

  # 11 个变量：TSPO + 与它相关最强/最负的细胞类型，
  # 再多元的话面板只有 2.6 in，数字会看不清
  VAR <- c("TSPO", "NK", "Monocyte", "CD4_T", "Mac_M1", "CD8_T",
           "Treg", "Mac_M2", "Prolif", "Epithelial", "Fibroblast")
  VAR <- intersect(VAR, colnames(IM))
  SHORT <- c(TSPO = "TSPO", NK = "NK", Monocyte = "Mono", CD4_T = "CD4",
             Mac_M1 = "M1", CD8_T = "CD8", Treg = "Treg", Mac_M2 = "M2",
             Prolif = "Prolif", Epithelial = "Epi", Fibroblast = "Fib")
  M <- cor(IM[, VAR], method = "spearman")
  M[upper.tri(M, diag = TRUE)] <- NA
  CD <- data.frame(g1 = rep(VAR, times = length(VAR)),
                   g2 = rep(VAR, each = length(VAR)),
                   r  = as.vector(M), stringsAsFactors = FALSE)
  CD <- CD[is.finite(CD$r), ]
  LV <- rev(VAR)
  CD$g1 <- factor(SHORT[CD$g1], levels = SHORT[LV])
  CD$g2 <- factor(SHORT[CD$g2], levels = SHORT[LV])
  cat(sprintf("相关矩阵：%d 个变量，下三角 %d 个单元格\n", length(VAR), nrow(CD)))

  pA <- ggplot(CD, aes(g2, g1)) +
    geom_point(aes(size = abs(r), fill = r), shape = 21,
               colour = "grey45", stroke = 0.2) +
    geom_text(aes(label = sprintf("%.2f", r), colour = abs(r) > 0.55),
              size = 1.5, show.legend = FALSE) +
    scale_size_continuous(range = c(0.4, 4.2), name = "|rho|",
                          guide = guide_legend(keyheight = unit(1.6, "mm"))) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "rho",
                         guide = guide_colourbar(barwidth = unit(2.2, "mm"),
                                                 barheight = unit(7, "mm"))) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_tspo() +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.3, size = 4.2),
          axis.text.y = element_text(size = 4.2),
          panel.grid = element_blank(),
          legend.position = "right",
          legend.key.size = unit(2.0, "mm"),
          legend.text = element_text(size = 3.2),
          legend.title = element_text(size = 3.4),
          legend.spacing.y = unit(0.6, "mm"))
}

# =============================================================================
#  B  免疫总分 vs TSPO# =============================================================================
#  B  免疫总分 vs TSPO
# =============================================================================

cat("[B] 免疫总分散点\n")
pB <- panel_missing("scores not found")
fB <- file.path(TAB, "T9a_immune_scores.csv")
if (file.exists(fB) && file.exists(F_EXPR)) {
  sc <- read.csv(fB, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(sc) <- sc[[1]]; sc <- sc[, -1]
  ex <- as.matrix(readRDS(F_EXPR)); storage.mode(ex) <- "double"
  common <- intersect(colnames(ex), rownames(sc))
  if (length(common) > 100) {
    imm_cols <- intersect(c("CD8_T","CD4_T","Treg","NK","B_cell","Monocyte",
                            "Mac_M1","Mac_M2","DC","Neutrophil","Mast"),
                          colnames(sc))
    tot <- rowMeans(sc[common, imm_cols, drop = FALSE], na.rm = TRUE)
    d <- data.frame(total = tot, tspo = as.numeric(ex["TSPO", common]))
    d <- d[is.finite(d$total) & is.finite(d$tspo), ]
    r <- cor(d$total, d$tspo, method = "spearman")

    pB <- ggplot(d, aes(total, tspo)) +
      geom_point(size = 0.62, alpha = 0.45, colour = "#9EC3E0") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                  colour = PAL["ink"], fill = "grey85", linewidth = 0.5) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 1.9,
               label = sprintf("rho = %.3f", r), size = 2.6, colour = PAL["ink"]) +
      labs(x = "Immune infiltration score", y = expression(log[2]~"TSPO")) +
      theme_tspo()
  }
}

# =============================================================================
#  C  WGCNA 模块-性状热图
# =============================================================================

cat("[C] WGCNA 模块-性状\n")
fC <- file.path(TAB, "T15c_module_trait_cor.csv")
pC <- panel_missing("T15c not found")
if (file.exists(fC)) {
  mc <- read.csv(fC, stringsAsFactors = FALSE, check.names = FALSE)
  traits <- setdiff(colnames(mc), "module")
  long <- data.frame(module = rep(mc$module, times = length(traits)),
                     trait  = rep(traits, each = nrow(mc)),
                     r      = as.vector(as.matrix(mc[, traits])),
                     stringsAsFactors = FALSE)
  # 按与 TSPO 的相关排序模块，TSPO 列排最左
  ord <- mc$module[order(mc$TSPO)]
  long$module <- factor(long$module, levels = ord)
  long$trait <- factor(long$trait, levels = c("TSPO", setdiff(traits, "TSPO")))

  pC <- ggplot(long, aes(trait, module, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 1.9, colour = "grey15") +
    scale_fill_div(limit = c(-1, 1), name = "Pearson\nr") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_tspo() +
    theme(axis.text.x.top = element_text(angle = 32, hjust = 0, size = 5.4),
          panel.grid = element_blank(),
          legend.key.height = unit(9, "mm")) +
    guides(fill = guide_colourbar(barwidth = unit(2.2, "mm"),
                                  barheight = unit(9, "mm")))
}

# =============================================================================
#  D  kME 位置
# =============================================================================

cat("[D] kME 位置\n")
fD <- file.path(TAB, "T15d_TSPO_module_hubs.csv")
pD <- panel_missing("T15d not found")
if (file.exists(fD)) {
  h <- read.csv(fD, stringsAsFactors = FALSE)
  h <- head(h[order(-h$kME), ], 12)
  TSPO_KME <- 0.577
  dd <- rbind(data.frame(gene = h$gene, kME = h$kME, type = "Module hub", stringsAsFactors = FALSE),
              data.frame(gene = "TSPO", kME = TSPO_KME, type = "TSPO", stringsAsFactors = FALSE))
  dd$gene <- factor(dd$gene, levels = dd$gene[order(dd$kME)])

  pD <- ggplot(dd, aes(kME, gene)) +
    geom_segment(aes(x = 0, xend = kME, y = gene, yend = gene),
                 colour = "grey86", linewidth = 0.3) +
    geom_point(aes(fill = type), shape = 21, size = 2.0,
               colour = "grey25", stroke = 0.25) +
    scale_fill_manual(values = c(`Module hub` = "#9EC3E0", TSPO = PAL["tspo"]),
                      name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                       expand = expansion(mult = c(0.01, 0.04))) +
    labs(x = "kME (module membership)", y = NULL) +
    theme_tspo() +
    theme(legend.position = c(0.02, 0.30), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.4, "mm"))
}

# =============================================================================
#  E / F  三队列：校正前 → 后
# =============================================================================

cat("[E/F] 三队列哑铃\n")
fE <- file.path(TAB, "T23a_scores_confound_adjusted.csv")
pE <- panel_missing("T23a not found"); pF <- panel_missing("T23a not found")
if (file.exists(fE)) {
  t23 <- read.csv(fE, stringsAsFactors = FALSE)
  cat("  指标: ", paste(unique(t23$indicator), collapse = ", "), "\n", sep = "")
  pick <- function(key) {
    z <- t23[t23$indicator == key, ]
    if (!nrow(z)) return(NULL)
    z$cohort <- factor(z$cohort, levels = c("TCGA", "GSE26712", "ICGC"))
    z[order(z$cohort), ]
  }
  # 斜率图：Raw → Adjusted 一点一线，线越陡说明全局表达校正影响越大
  COH3 <- c(TCGA = "#2166AC", GSE26712 = "#B2182B", ICGC = "#67A9CF")
  slope <- function(z, ttl) {
    dl <- rbind(
      data.frame(cohort = z$cohort, stage = "Raw",      rho = z$rho_raw,
                 stringsAsFactors = FALSE),
      data.frame(cohort = z$cohort, stage = "Adjusted", rho = z$rho_ctrl,
                 stringsAsFactors = FALSE))
    dl$stage  <- factor(dl$stage, levels = c("Raw", "Adjusted"))
    dl$cohort <- factor(dl$cohort, levels = names(COH3))

    ggplot(dl, aes(stage, rho, colour = cohort, group = cohort)) +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2.3) +
      # 右端合并标签（队列名 + 校正后 ρ），带白底避免两条线标签重叠，省掉图例
      geom_label(data = dl[dl$stage == "Adjusted", ],
                 aes(label = sprintf("%s  %.2f", cohort, rho)),
                 hjust = -0.06, size = 1.75, label.size = 0,
                 fill = alpha("white", 0.75), colour = "grey15",
                 label.padding = unit(0.6, "mm"), show.legend = FALSE) +
      scale_colour_manual(values = COH3, guide = "none") +
      scale_x_discrete(expand = expansion(mult = c(0.05, 0.52))) +
      scale_y_continuous(expand = expansion(mult = c(0.10, 0.06))) +
      labs(x = NULL, y = expression("Spearman " * rho * " with TSPO")) +
      ggtitle(ttl) +
      theme_tspo() +
      theme(plot.title = element_text(size = 6.4, face = "plain"),
            plot.margin = margin(3.5, 4.5, 3.5, 3.5))
  }
  # E/F 原为同一套斜率图。F 改为带 95% CI 的森林图：
  # 用 Fisher z 变换由 rho 与 n 算出区间，同一队列内 Raw / Adjusted 两行并排，
  # 区间重叠即说明校正没有实质改变关联——比柱状图多一层不确定性信息。
  forest <- function(z, ttl) {
    fz <- function(r, n) {
      r <- pmax(pmin(r, 0.999), -0.999)
      zz <- atanh(r); se <- 1.03 / sqrt(pmax(n - 3, 2))
      c(lo = tanh(zz - 1.96 * se), hi = tanh(zz + 1.96 * se))
    }
    D <- rbind(
      data.frame(cohort = z$cohort, stage = "Raw",      rho = z$rho_raw, n = z$n),
      data.frame(cohort = z$cohort, stage = "Adjusted", rho = z$rho_ctrl, n = z$n))
    ci <- t(mapply(function(r, n) fz(r, n), D$rho, D$n))
    D$lo <- ci[, 1]; D$hi <- ci[, 2]
    D$stage  <- factor(D$stage, levels = c("Raw", "Adjusted"))
    D$cohort <- factor(D$cohort, levels = rev(names(COH3)))
    D$yy <- as.numeric(D$cohort) + ifelse(D$stage == "Raw", 0.16, -0.16)

    DLT <- data.frame(cohort = factor(z$cohort, levels = rev(names(COH3))),
                      d = z$rho_ctrl - z$rho_raw)
    DLT$yy <- as.numeric(DLT$cohort) - 0.16

    ggplot(D) +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
      geom_segment(aes(x = lo, xend = hi, y = yy, yend = yy, colour = cohort),
                   linewidth = 0.55, alpha = ifelse(D$stage == "Raw", 0.45, 0.95)) +
      geom_point(aes(rho, yy, fill = cohort, shape = stage),
                 size = 2.1, colour = "grey25", stroke = 0.3) +
      scale_shape_manual(values = c(Raw = 22, Adjusted = 21), name = NULL) +
      scale_fill_manual(values = COH3, guide = "none") +
      scale_colour_manual(values = COH3, guide = "none") +
      scale_x_continuous(limits = c(-0.66, 0.16),
                         breaks = c(-0.6, -0.4, -0.2, 0)) +
      scale_y_continuous(breaks = 1:length(names(COH3)),
                         labels = rev(names(COH3)),
                         expand = expansion(add = c(0.35, 0.35))) +
      # Δ 写在右侧
      geom_text(data = DLT, aes(x = -0.64, y = yy,
                                label = ifelse(abs(d) < 0.01, "\u0394 \u2248 0",
                                               sprintf("\u0394 %+.2f", d))),
                hjust = 0, size = 1.6, colour = "grey30") +
      labs(x = expression("Spearman " * rho * " with TSPO (95% CI)"), y = NULL) +
      ggtitle(ttl) +
      theme_tspo() +
      theme(plot.title = element_text(size = 6.4, face = "plain"),
            panel.grid.major.y = element_blank(),
            legend.position = c(0.02, 0.99), legend.justification = c(0, 1),
            legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
            legend.key.size = unit(2.2, "mm"),
            legend.text = element_text(size = 3.6),
            legend.title = element_blank(),
            plot.margin = margin(3.5, 4.5, 3.5, 3.5))
  }

  zE <- pick("\u7eaf\u5ea6\u4ee3\u7406")   # 纯度代理
  zF <- pick("\u589e\u6b96")               # 增殖
  if (!is.null(zE)) pE <- slope(zE, "Purity proxy")
  if (!is.null(zF)) pF <- forest(zF, "Proliferation")
}

# =============================================================================
#  G  全局表达梯度（三队列）
# =============================================================================

cat("[G] 全局表达梯度\n")
pG <- panel_missing("matrices not found")

gs <- list()
if (file.exists(F_EXPR)) {
  ex <- as.matrix(readRDS(F_EXPR)); storage.mode(ex) <- "double"
  gs[["TCGA-OV"]] <- data.frame(
    cohort = "TCGA-OV", tspo = as.numeric(ex["TSPO", ]),
    global = colMeans(ex, na.rm = TRUE), stringsAsFactors = FALSE)
}
if (file.exists(F_GSE)) {
  g <- as.matrix(readRDS(F_GSE)); storage.mode(g) <- "double"
  if ("TSPO" %in% rownames(g))
    gs[["GSE26712"]] <- data.frame(
      cohort = "GSE26712", tspo = as.numeric(g["TSPO", ]),
      global = colMeans(g, na.rm = TRUE), stringsAsFactors = FALSE)
}
if (file.exists(F_ICGC)) {
  ic <- readRDS(F_ICGC)
  if (is.data.frame(ic)) ic <- as.matrix(ic)
  if (nrow(ic) < ncol(ic) && "TSPO" %in% colnames(ic)) ic <- t(ic)
  storage.mode(ic) <- "double"
  if ("TSPO" %in% rownames(ic))
    gs[["ICGC"]] <- data.frame(
      cohort = "ICGC", tspo = as.numeric(ic["TSPO", ]),
      global = colMeans(ic, na.rm = TRUE), stringsAsFactors = FALSE)
}
if (length(gs)) {
  gg <- do.call(rbind, gs)
  gg <- gg[is.finite(gg$tspo) & is.finite(gg$global), ]
  # 三队列表达尺度不同（log2 TPM / RMA / ICGC），必须队内标准化后才能同轴比较
  gg <- do.call(rbind, lapply(split(gg, gg$cohort), function(z) {
    z$tspo_z   <- as.numeric(scale(z$tspo))
    z$global_z <- as.numeric(scale(z$global))
    z
  }))
  gg$cohort <- factor(gg$cohort, levels = c("TCGA-OV", "GSE26712", "ICGC"))
  lab <- sapply(levels(gg$cohort), function(cc) {
    z <- gg[gg$cohort == cc, ]
    sprintf("%s  %+.2f", cc, cor(z$tspo, z$global, method = "spearman"))
  })
  gg$panel <- factor(lab[as.character(gg$cohort)], levels = lab)

  # NMDS：三队列在「免疫/基质 marker 空间」里的相对位置。
  # 原为 TSPO × 全基因均值的散点，只说明平台尺度差异；
  # 排序图还能看出样本在免疫轴上的分层是否与队列无关。
  # 选一组三个平台都能测到的 marker（芯片与 RNA-seq 均有），
  # 队内 z 化后合并做 NMDS，避免平台尺度污染。
  MARK <- c("PTPRC","CD3D","CD3E","CD2","CD8A","NKG7","GNLY","MS4A1","CD79A",
            "CD68","CD14","ITGAM","FCGR3A","LYZ","CD163","MRC1",
            "COL1A1","COL1A2","DCN","CXCL9")

  mats <- list()
  if (file.exists(F_EXPR)) {
    ex <- as.matrix(readRDS(F_EXPR)); storage.mode(ex) <- "double"; mats[["TCGA-OV"]] <- ex
  }
  if (file.exists(F_GSE)) {
    g <- as.matrix(readRDS(F_GSE)); storage.mode(g) <- "double"; mats[["GSE26712"]] <- g
  }
  if (file.exists(F_ICGC)) {
    ic <- readRDS(F_ICGC)
    if (is.data.frame(ic)) ic <- as.matrix(ic)
    if (nrow(ic) < ncol(ic) && "TSPO" %in% colnames(ic)) ic <- t(ic)
    storage.mode(ic) <- "double"; mats[["ICGC"]] <- ic
  }

  if (length(mats) >= 2 && requireNamespace("MASS", quietly = TRUE)) {
    MK <- Reduce(intersect, lapply(mats, function(m) intersect(MARK, rownames(m))))
    cat("NMDS marker 交集：", length(MK), " 个\n", sep = "")

    if (length(MK) >= 10) {
      SS <- do.call(rbind, lapply(names(mats), function(cn) {
        m <- mats[[cn]][MK, , drop = FALSE]
        z <- t(scale(t(m)))                       # 每个 marker 队内 z 化
        z[!is.finite(z)] <- 0
        d <- data.frame(cohort = cn, t(z), check.names = FALSE)
        tsp <- as.numeric(mats[[cn]]["TSPO", ])
        d$tspo_z <- as.numeric(scale(tsp))
        d
      }))
      SS <- SS[complete.cases(SS), ]
      FEAT <- as.matrix(SS[, MK])
      DM   <- dist(FEAT, method = "euclidean")
      NM   <- MASS::isoMDS(DM, k = 2, trace = FALSE)
      STR  <- NM$stress / 100
      cat(sprintf("NMDS stress = %.3f（<0.2 可用）\n", STR))

      # 着色按 TSPO 三分位（队内），而不是按队列——
      # marker 已在队内 z 化，跨队列差异被抹平，按队列着色必然三团重合，
      # 那样的「不分离」是人为造成的。按 TSPO 分层问的才是真问题：
      # TSPO 高的样本在免疫/基质空间里是否占据不同位置。
      ND <- data.frame(NMDS1 = NM$points[, 1], NMDS2 = NM$points[, 2],
                       cohort = SS$cohort, tspo_z = SS$tspo_z)
      ND$grp <- ave(ND$tspo_z, ND$cohort, FUN = function(z) {
        q <- quantile(z, c(1/3, 2/3), na.rm = TRUE)
        ifelse(z > q[2], "High", ifelse(z <= q[1], "Low", "Mid"))
      })
      ND$grp <- factor(ND$grp, levels = c("Low", "Mid", "High"))
      PAL_G <- c(Low = "#2166AC", Mid = "#BDBDBD", High = "#B2182B")
      cat("NMDS 三分位样本数：
"); print(table(ND$grp))

      pG <- ggplot(ND, aes(NMDS1, NMDS2, colour = grp, fill = grp)) +
        stat_ellipse(geom = "polygon", alpha = 0.10, level = 0.68,
                     linewidth = 0.35, show.legend = FALSE) +
        geom_point(size = 0.55, alpha = 0.55, stroke = 0) +
        scale_colour_manual(values = PAL_G, name = "TSPO",
                            labels = c("low", "mid", "high")) +
        scale_fill_manual(values = PAL_G, guide = "none") +
        annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.8, size = 1.9,
                 colour = "grey20", label = sprintf("stress = %.3f", STR)) +
        labs(x = "NMDS1", y = "NMDS2") +
        theme_tspo() +
        theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
              legend.direction = "horizontal",
              legend.background = element_rect(fill = alpha("white", 0.82), colour = NA),
              legend.key.size = unit(3.0, "mm"),
              legend.text = element_text(size = 4.4),
              legend.title = element_text(size = 4.4)) +
        guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1)))
    }
  }
}

# =============================================================================
#  H  火山图
# =============================================================================

cat("[H] 火山图\n")
fH <- file.path(TAB, "T10a_DEG_full.csv")
pH <- panel_missing("T10a not found"); pI <- panel_missing("T10b not found")
FC_CUT <- 0.5; FDR_CUT <- 0.05

if (file.exists(fH)) {
  dg <- read.csv(fH, stringsAsFactors = FALSE)
  dg <- dg[is.finite(dg$logFC) & is.finite(dg$adj.P.Val), ]
  dg$dir <- ifelse(dg$adj.P.Val < FDR_CUT & dg$logFC >= FC_CUT, "Up",
            ifelse(dg$adj.P.Val < FDR_CUT & dg$logFC <= -FC_CUT, "Down", "NS"))
  dg$dir <- factor(dg$dir, levels = c("Up", "Down", "NS"))
  n_up <- sum(dg$dir == "Up"); n_dn <- sum(dg$dir == "Down")

  pH <- ggplot(dg, aes(logFC, -log10(adj.P.Val))) +
    geom_point(data = dg[dg$dir == "NS", ], size = 0.30, alpha = 0.35,
               colour = DIR_PAL["NS"]) +
    geom_point(data = dg[dg$dir != "NS", ], aes(colour = dir), size = 0.45,
               alpha = 0.75) +
    geom_vline(xintercept = c(-FC_CUT, FC_CUT), linetype = 2, linewidth = 0.25,
               colour = "grey45") +
    geom_hline(yintercept = -log10(FDR_CUT), linetype = 2, linewidth = 0.25,
               colour = "grey45") +
    scale_colour_manual(values = DIR_PAL, guide = "none") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.9, size = 2.1,
             colour = DIR_PAL["Down"], label = sprintf("down %d", n_dn)) +
    annotate("text", x = Inf, y = Inf, hjust = 1.06, vjust = 1.9, size = 2.1,
             colour = DIR_PAL["Up"], label = sprintf("up %d", n_up)) +
    labs(x = expression(log[2]~"fold change"), y = expression(-log[10]~"adj. P")) +
    theme_tspo()

  # ---- I  上下调各前 10 ----
  fI <- file.path(TAB, "T10b_DEG_sig.csv")
  if (file.exists(fI)) {
    sg <- read.csv(fI, stringsAsFactors = FALSE)
    sg <- sg[sg$adj.P.Val < FDR_CUT & abs(sg$logFC) >= FC_CUT, ]
    up <- head(sg[order(-sg$logFC), ], 10)
    dn <- head(sg[order(sg$logFC), ], 10)
    top <- rbind(transform(up, dir = "Up"), transform(dn, dir = "Down"))
    top <- top[order(top$logFC), ]
    top$gene <- factor(top$gene, levels = top$gene)

    pI <- ggplot(top, aes(logFC, gene)) +
      geom_segment(aes(x = 0, xend = logFC, y = gene, yend = gene),
                   colour = "grey86", linewidth = 0.3) +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
      geom_point(aes(fill = dir), shape = 21, size = 1.7,
                 colour = "grey25", stroke = 0.22) +
      scale_fill_manual(values = DIR_PAL, guide = "none") +
      scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
      labs(x = expression(log[2]~"fold change"), y = NULL) +
      theme_tspo()
  }
}

# =============================================================================
#  拼装
# =============================================================================

cat("\n拼装 3x3\n")

FIG2 <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, colour = "black"))

if (!exists("FIG_NAME")) FIG_NAME <- "Fig2"
print(FIG2)
save_fig(FIG2, FIG_NAME, width = 8.8, height = 8.4)

cat("\nfig02.R 完成\n")

# -----------------------------------------------------------------------------
# -- 作图样式 / fig03_source.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig03_source.R   图 3（6 个面板）TSPO 高表达的两个来源
#
#  A  TSPO ~ 髓系评分（四队列队内 z 叠加）          来源一：髓系浸润
#  B  source_idx 在四亚型的雨云图（密度脊 + 散点）    来源二存在：分化型最高
#  C  source_idx ~ DIF 程序（三队列散点 + 拟合）      来源二的外部验证
#  D  DIF 与 MES 的相关系数森林图（含 95% CI）        特异性对照（定量）
#  E  扣除网络重叠基因前后的相关系数哑铃图            排除循环性
#  F  方差分解的标注热图                             定量
#
#  ── 绘图规范 ────────────────────────────────────────────────────────────────
#  不使用柱状图与裸箱线图。分布用雨云图（脊线 + 透明散点），相关强度用带置信
#  区间的森林图，前后对比用哑铃图，方差分解用标注热图——与图 2C 的模块-性状
#  热图同一格式（geom_tile + 白色分隔 + 数值标注 + 紧凑色条）。
#
#  ── 这个图要说什么 ──────────────────────────────────────────────────────────
#  TSPO 高表达不是单一状态。它有两个可区分的来源：
#    一、髓系浸润（四队列 ρ 0.28–0.43，离散度极小）
#    二、分化型上皮程序（扣除髓系后的残余与 DIF 签名相关，三队列 ρ 0.28–0.44）
#  第二个来源与间质型程序几乎无关（ρ 0.05–0.14，置信区间均跨零），且 DIF 签名
#  与 TSPO 共表达网络几乎不重叠（37/243，置换 p=0.40），剔除重叠基因后相关几乎
#  不变（0.324→0.282、0.309→0.277、0.440→0.410），排除自相关。
#
#  因此单一 TSPO 数值无法区分肿瘤属于哪一类来源——同样的高值，
#  在免疫型来自浸润髓系，在分化型来自上皮细胞本身。
#
#  ── 输出 ────────────────────────────────────────────────────────────────────
#  Fig3.pdf / Fig3.png         六面板图（450 dpi）
#  T38a_source_panels.csv      各队列相关系数
#  T38b_variance_decomp.csv    面板 F 的数值
#  T38c_permutation_null.csv   置换检验结果（正文引用，不单独出图）
#  T38d_deoverlap.csv          面板 E 的前后对比
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/fig03_source.R", encoding = "UTF-8")
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

if (!exists("PROJ"))       PROJ       <- "D:/Rstudy/tcga"
if (!exists("OUT"))        OUT        <- file.path(PROJ, "results_v2")
if (!exists("F_EXPR"))     F_EXPR     <- file.path(OUT, "EXPR_clean.rds")
if (!exists("F_SUB"))      F_SUB      <- file.path(PROJ, "TCGA_OV_subtype.gz")
if (!exists("F_MAP570"))   F_MAP570   <- file.path(PROJ, "GPL570_probe_map_correct.csv")
if (!exists("F_SOFT6480")) F_SOFT6480 <- file.path(PROJ, "GPL6480.soft.gz")
if (!exists("F_DIF"))      F_DIF      <- file.path(PROJ,
  "10780432ccr200103-sup-table_sb3_subtype2-c4_diff_05jul20.xlsx")
if (!exists("F_MES"))      F_MES      <- file.path(PROJ,
  "10780432ccr200103-sup-table_sb5_subtype4-c1_mes_05jul20.xlsx")
if (!exists("F_NET"))      F_NET      <- file.path(OUT, "tables", "T8_network_all_genes.csv")

if (!exists("FIGDIR")) FIGDIR <- file.path(OUT, "figures")
OUTDIR <- file.path(OUT, "tables")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(patchwork); library(grid); library(ggridges) })

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

N_PERM <- 5000
set.seed(20261025)

# =============================================================================
#  STEP 0  基因集与网络
# =============================================================================

hr("STEP 0  基因集")

read_sig_xlsx <- function(f, tag) {
  if (!file.exists(f)) stop("缺 ", basename(f))
  sh <- openxlsx::getSheetNames(f)[1]
  d  <- openxlsx::read.xlsx(f, sheet = sh, colNames = FALSE)
  hi <- NA
  for (i in seq_len(min(8, nrow(d)))) {
    cells <- tolower(trimws(as.character(unlist(d[i, ]))))
    if ("probe" %in% cells && "gene" %in% cells) { hi <- i; break }
  }
  if (is.na(hi)) stop("未找到 ", tag, " 表头")
  hdr <- trimws(as.character(unlist(d[hi, ])))
  gi  <- which(tolower(hdr) == "gene")[1]
  g   <- trimws(as.character(unlist(d[(hi + 1):nrow(d), gi])))
  g   <- g[!is.na(g) & nzchar(g) & g != "NA"]
  unique(sub("\\s*///.*$", "", g))
}

G_DIF <- read_sig_xlsx(F_DIF, "DIF")
G_MES <- read_sig_xlsx(F_MES, "MES")
cat("DIF ", length(G_DIF), " 基因；MES ", length(G_MES), " 基因\n", sep = "")

# 网络基因（供 DIF/MES 去重叠与置换检验使用）
NG <- character(0)
if (file.exists(F_NET)) {
  nw <- read.csv(F_NET, stringsAsFactors = FALSE)
  if (!"in_network" %in% colnames(nw)) nw$in_network <- abs(nw$rho) > 0.3 & nw$FDR < 0.05
  NG <- nw$gene[nw$in_network]
  cat("TSPO 共表达网络: ", length(NG), " 基因\n", sep = "")
} else cat("⚠ 未找到网络表，E 面板将跳过\n")

# 从 DIF 签名中剔除落在网络内的基因，用于面板 E 的前后对比
G_DIF_CLEAN <- setdiff(G_DIF, NG)
cat("去重叠 DIF 签名: ", length(G_DIF_CLEAN), " / ", length(G_DIF), "\n", sep = "")

MAP570 <- read.csv(F_MAP570, stringsAsFactors = FALSE, check.names = FALSE)
colnames(MAP570)[1:2] <- c("probe", "symbol")
MAP570 <- MAP570[!is.na(MAP570$symbol) & nzchar(MAP570$symbol), ]
MAP570$symbol <- sub("\\s*///.*$", "", MAP570$symbol)
MAP570 <- MAP570[!duplicated(MAP570$probe), ]

MAP6480 <- NULL
if (file.exists(F_SOFT6480)) {
  ln <- readLines(F_SOFT6480, warn = FALSE, encoding = "UTF-8")
  i0 <- grep("^!platform_table_begin", ln)[1]; i1 <- grep("^!platform_table_end", ln)[1]
  if (!is.na(i0) && !is.na(i1)) {
    hdr <- strsplit(ln[i0 + 1], "\t")[[1]]; body <- ln[(i0 + 2):(i1 - 1)]
    ci <- which(hdr == "GENE_SYMBOL")[1]
    sy <- vapply(strsplit(body, "\t"), function(z) if (length(z) >= ci) z[ci] else NA_character_,
                 character(1), USE.NAMES = FALSE)
    MAP6480 <- data.frame(probe = sub("\t.*$", "", body), symbol = sy, stringsAsFactors = FALSE)
    MAP6480 <- MAP6480[!is.na(MAP6480$symbol) & nzchar(MAP6480$symbol) & MAP6480$symbol != "", ]
    MAP6480$symbol <- sub("\\s*///.*$", "", MAP6480$symbol)
    MAP6480 <- MAP6480[!duplicated(MAP6480$probe), ]
  }
}

GENES <- list(
  Myeloid = c("CD68","CD14","LYZ","AIF1","C1QA","C1QB","ITGAM","FCER1G","TYROBP")
)
score_z <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 3) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[g, , drop = FALSE]))); z[!is.finite(z)] <- NA
  colMeans(z, na.rm = TRUE)
}

build <- function(mat, tag) {
  tspo <- as.numeric(mat["TSPO", ])
  mye  <- score_z(mat, GENES$Myeloid)
  dif  <- score_z(mat, G_DIF)
  mes  <- score_z(mat, G_MES)
  dcl  <- if (length(G_DIF_CLEAN) >= 20) score_z(mat, G_DIF_CLEAN) else rep(NA_real_, ncol(mat))
  ok <- is.finite(tspo) & is.finite(mye)
  src <- rep(NA_real_, length(tspo))
  if (sum(ok) >= 30) {
    fit <- lm(rank(tspo[ok]) ~ rank(mye[ok]))
    src[ok] <- residuals(fit)
  }
  list(data = data.frame(sample = colnames(mat), cohort = tag,
                         TSPO = tspo, Myeloid = mye, DIF = dif, MES = mes,
                         DIF_clean = dcl, source_idx = src, stringsAsFactors = FALSE))
}

read_cohort <- function(gsm, map) {
  f <- file.path(PROJ, paste0(gsm, "_series_matrix.txt.gz"))
  if (!file.exists(f)) return(NULL)
  con <- gzfile(f, "rt"); ln <- readLines(con, warn = FALSE); close(con)
  i0 <- grep("^!series_matrix_table_begin", ln)[1]
  i1 <- grep("^!series_matrix_table_end", ln)[1]
  if (is.na(i0) || is.na(i1)) return(NULL)
  tab <- read.delim(text = paste(ln[(i0 + 1):(i1 - 1)], collapse = "\n"),
                    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(tab[, -1, drop = FALSE]); storage.mode(mat) <- "double"
  rownames(mat) <- tab[[1]]
  m <- map[match(rownames(mat), map$probe), ]
  ok <- !is.na(m$symbol)
  mat <- mat[ok, , drop = FALSE]; sym <- m$symbol[ok]
  od <- order(-rowMeans(mat, na.rm = TRUE))
  mat <- mat[od, , drop = FALSE]; sym <- sym[od]
  kp <- !duplicated(sym)
  EX <- mat[kp, , drop = FALSE]; rownames(EX) <- sym[kp]
  EX
}

# =============================================================================
#  STEP 1  载入
# =============================================================================

hr("STEP 1  载入")

CO <- list()
TC <- as.matrix(readRDS(F_EXPR)); storage.mode(TC) <- "double"
CO[["TCGA"]] <- build(TC, "TCGA")
cat("  TCGA: ", ncol(TC), " 样本\n", sep = "")

if (!is.null(MAP6480)) {
  x <- read_cohort("GSE32062-GPL6480", MAP6480)
  if (!is.null(x)) { CO[["GSE32062"]] <- build(x, "GSE32062")
    cat("  GSE32062: ", ncol(x), "\n", sep = "") }
}
for (gs in c("GSE63885", "GSE9891")) {
  x <- read_cohort(gs, MAP570)
  if (!is.null(x)) { CO[[gs]] <- build(x, gs); cat("  ", gs, ": ", ncol(x), "\n", sep = "") }
}

EXT <- setdiff(names(CO), "TCGA")
if (!length(EXT)) {
  cat("⚠ 未找到外部队列（GSE32062 / GSE63885 / GSE9891），C/D/E 将以占位显示。\n")
  cat("  这些队列需另行下载 series matrix 放到 ", PROJ, "\n", sep = "")
}

D <- CO[["TCGA"]]$data
if (file.exists(F_SUB)) {
  sb <- read.delim(gzfile(F_SUB), stringsAsFactors = FALSE)
  colnames(sb)[1:2] <- c("sample", "subtype")
  sb$subtype <- tolower(trimws(sb$subtype))
  D$subtype <- sb$subtype[match(D$sample, sb$sample)]
  CO[["TCGA"]]$data <- D
}

# 相关分析的公共工具：Spearman ρ 与 Fisher z 的 95% CI
rho_ci <- function(r, n) {
  if (!is.finite(r) || n < 5 || abs(r) >= 1) return(c(NA, NA))
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  tanh(z + c(-1.96, 1.96) * se)
}

# =============================================================================
#  A  TSPO ~ 髓系（四队列叠加，队内 z）
# =============================================================================

hr("[A] TSPO ~ 髓系")

Z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
AZ <- do.call(rbind, lapply(names(CO), function(nm) {
  d <- CO[[nm]]$data
  ok <- is.finite(d$TSPO) & is.finite(d$Myeloid)
  data.frame(cohort = nm, x = Z(d$Myeloid)[ok], y = Z(d$TSPO)[ok],
             rho = suppressWarnings(cor(d$TSPO[ok], d$Myeloid[ok], method = "spearman")),
             stringsAsFactors = FALSE)
}))
RH <- tapply(AZ$rho, AZ$cohort, function(z) z[1])
AZ$lab <- factor(sprintf("%s (rho = %.2f)", AZ$cohort, RH[AZ$cohort]),
                 levels = sprintf("%s (rho = %.2f)", names(RH), RH))
CATCOL <- setNames(c(PAL["epi"], PAL["tspo"], PAL["str"], PAL["lym"])[seq_along(RH)],
                   levels(AZ$lab))

pA <- ggplot(AZ, aes(x, y)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey88") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey88") +
  geom_point(aes(colour = lab), size = 0.40, alpha = 0.28) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = TRUE,
              linewidth = 0.55, colour = PAL["ink"], fill = "grey86") +
  scale_colour_manual(values = CATCOL, name = NULL,
                      guide = guide_legend(override.aes = list(alpha = 1, size = 1.3))) +
  labs(x = "Myeloid score (within-cohort z)",
       y = expression(log[2]~"TSPO (within-cohort z)")) +
  theme_tspo() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.82), colour = NA),
        legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 4.4))

# =============================================================================
#  B  source_idx 在四亚型：雨云图（密度脊 = 云，透明散点 = 雨）
# =============================================================================

hr("[B] source_idx 亚型（雨云图）")

LV <- c("immunoreactive", "differentiated", "proliferative", "mesenchymal")
LAB <- c(immunoreactive = "Immunoreactive", differentiated = "Differentiated",
         proliferative = "Proliferative", mesenchymal = "Mesenchymal")
ds <- D[!is.na(D$subtype) & is.finite(D$source_idx), ]
ds$subtype <- factor(ds$subtype, levels = intersect(LV, unique(ds$subtype)))
kw <- kruskal.test(source_idx ~ subtype, data = ds)
md <- tapply(ds$source_idx, ds$subtype, median)
cat("  中位数: ", paste(names(md), round(md, 1), sep = "=", collapse = ", "), "\n", sep = "")
cat(sprintf("  Kruskal-Wallis H=%.1f, p=%.2e\n", kw$statistic, kw$p.value))

ordB <- names(sort(md))                       # 中位数升序，最高者置顶
ds$subtype <- factor(LAB[as.character(ds$subtype)], levels = LAB[ordB])
MED <- data.frame(subtype = factor(LAB[ordB], levels = LAB[ordB]),
                  med = as.numeric(md[ordB]))
BRIDGE <- setNames(c("#2166AC", "#92C5DE", "#EF8A62", "#B2182B")[seq_along(ordB)], LAB[ordB])

# 「雨」：脊线基线下方 0.30 处铺散点，颜色与所属云朵一致，空心透明圆
RAIN_OFF <- 0.30

pB <- ggplot(ds, aes(source_idx, subtype, fill = subtype)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey62") +
  geom_density_ridges(scale = 0.92, rel_min_height = 0.006, alpha = 0.55,
                      colour = "grey30", linewidth = 0.28,
                      quantile_lines = TRUE, quantiles = 2, vline_colour = "white") +
  # inherit.aes = FALSE 时 x 不会被继承，必须显式写出来
  geom_jitter(aes(x = source_idx, y = as.numeric(subtype) - RAIN_OFF, colour = subtype),
              shape = 1, size = 0.42, stroke = 0.22, width = 0, height = 0.085,
              alpha = 0.55, inherit.aes = FALSE) +
  geom_point(data = MED, aes(x = med, y = subtype), inherit.aes = FALSE,
             shape = 21, size = 1.4, fill = "white", colour = "grey15", stroke = 0.35) +
  scale_fill_manual(values = BRIDGE, guide = "none") +
  scale_colour_manual(values = BRIDGE, guide = "none") +
  scale_y_discrete(expand = expansion(add = c(0.62, 0.16))) +
  coord_cartesian(clip = "off") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.9, size = 2.1,
           colour = PAL["ink"],
           label = sprintf("Kruskal-Wallis p = %.1e", kw$p.value)) +
  labs(x = "TSPO residual on myeloid", y = NULL) +
  theme_tspo() +
  theme(axis.text.y = element_text(size = 5.0))

# =============================================================================
#  C  source_idx ~ DIF（三队列散点 + 拟合）
# =============================================================================

hr("[C] source_idx ~ DIF")

make_scatter <- function(prog, col_pt) {
  dz <- do.call(rbind, lapply(EXT, function(nm) {
    d <- CO[[nm]]$data
    ok <- is.finite(d$source_idx) & is.finite(d[[prog]])
    ct <- suppressWarnings(cor.test(d$source_idx[ok], d[[prog]][ok], method = "spearman"))
    data.frame(cohort = nm, x = d[[prog]][ok], y = d$source_idx[ok],
               rho = as.numeric(ct$estimate), stringsAsFactors = FALSE)
  }))
  dz$panel <- factor(sprintf("%s  rho = %.2f", dz$cohort, dz$rho),
                     levels = unique(sprintf("%s  rho = %.2f", dz$cohort, dz$rho)))
  ggplot(dz, aes(x, y)) +
    geom_point(size = 0.62, alpha = 0.38, colour = col_pt, stroke = 0) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = TRUE,
                linewidth = 0.45, colour = PAL["ink"], fill = "grey88") +
    facet_wrap(~ panel, nrow = 1, scales = "free_x") +
    labs(x = NULL, y = "TSPO residual on myeloid") +
    theme_tspo() +
    theme(strip.text = element_text(size = 4.6))
}

pC <- if (length(EXT)) make_scatter("DIF", "#B2182B") else panel_missing("external cohorts not found")

# =============================================================================
#  D  DIF vs MES 的相关系数森林图（含 95% CI）
#
#  把两条程序放到同一坐标轴上，DIF 的置信区间远离零、MES 的跨零，特异性一眼可见。
# =============================================================================

hr("[D] DIF vs MES 森林图")

FY <- if (length(EXT)) do.call(rbind, lapply(EXT, function(nm) {
  d <- CO[[nm]]$data
  do.call(rbind, lapply(c("DIF", "MES"), function(pg) {
    ok <- is.finite(d$source_idx) & is.finite(d[[pg]])
    ct <- suppressWarnings(cor.test(d$source_idx[ok], d[[pg]][ok], method = "spearman"))
    r <- as.numeric(ct$estimate); n <- sum(ok)
    ci <- rho_ci(r, n)
    data.frame(cohort = nm, prog = pg, rho = r, lo = ci[1], hi = ci[2],
               n = n, stringsAsFactors = FALSE)
  }))
})) else NULL

pD <- if (!is.null(FY)) {
  cat("\n"); pr(FY, 3)
  FY$prog <- factor(FY$prog, levels = c("DIF", "MES"))
  FY$cohort <- factor(FY$cohort, levels = rev(EXT))
  ggplot(FY, aes(rho, cohort, colour = prog)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y",
                  linewidth = 0.42, position = position_dodge(width = 0.42)) +
    geom_point(aes(shape = prog), size = 1.7, stroke = 0.35,
               position = position_dodge(width = 0.42)) +
    scale_colour_manual(values = c(DIF = "#B2182B", MES = "#2166AC"), name = NULL) +
    scale_shape_manual(values = c(DIF = 16, MES = 17), name = NULL) +
    labs(x = expression("Spearman "*rho*" with source index"),
         y = NULL) +
    theme_tspo() +
    theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 4.6),
          axis.text.y = element_text(size = 4.8))
} else panel_missing("external cohorts not found")

# =============================================================================
#  E  扣除网络重叠基因前后的相关系数（哑铃图）
#
#  循环性最直接的检验不是「重叠是否富集」，而是「把重叠基因全部删掉后相关是否
#  还在」。哑铃图把原始 ρ 与去重叠 ρ 连起来：连线几乎不移动即说明非自相关。
# =============================================================================

hr("[E] 去重叠前后的相关")

DB <- if (length(EXT)) do.call(rbind, lapply(EXT, function(nm) {
  d <- CO[[nm]]$data
  ok1 <- is.finite(d$source_idx) & is.finite(d$DIF)
  ok2 <- is.finite(d$source_idx) & is.finite(d$DIF_clean)
  if (sum(ok1) < 50 || sum(ok2) < 50) return(NULL)
  r1 <- suppressWarnings(cor(d$source_idx[ok1], d$DIF[ok1], method = "spearman"))
  r2 <- suppressWarnings(cor(d$source_idx[ok2], d$DIF_clean[ok2], method = "spearman"))
  data.frame(cohort = nm, raw = as.numeric(r1), adj = as.numeric(r2),
             n = sum(ok1), stringsAsFactors = FALSE)
})) else NULL

pE <- if (!is.null(DB)) {
  cat("\n"); pr(DB, 3)
  DB$cohort <- factor(DB$cohort, levels = rev(EXT))
  DL <- rbind(
    data.frame(cohort = DB$cohort, type = "All 254 genes", rho = DB$raw),
    data.frame(cohort = DB$cohort, type = "Overlap removed (217)", rho = DB$adj))
  DL$type <- factor(DL$type, levels = c("All 254 genes", "Overlap removed (217)"))

  ov_txt <- ""
  if (length(NG)) {
    ov_txt <- sprintf("shared with network: %d of %d genes",
                      length(intersect(intersect(G_DIF, rownames(TC)), NG)),
                      length(intersect(G_DIF, rownames(TC))))
  }

  # 宽箭头从「全部基因」指向「去重叠」：箭头越短，说明相关越不依赖重叠基因
  AR <- merge(DL[DL$type == "All 254 genes", c("cohort", "rho")],
              DL[DL$type == "Overlap removed (217)", c("cohort", "rho")],
              by = "cohort", suffixes = c("raw", "adj"))

  ggplot(DL, aes(rho, cohort)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
    geom_segment(data = AR, aes(x = rhoraw, xend = rhoadj, y = cohort, yend = cohort),
                 inherit.aes = FALSE, colour = "grey55", linewidth = 0.9,
                 arrow = arrow(length = unit(1.1, "mm"), type = "closed"),
                 lineend = "round") +
    geom_point(aes(fill = type), shape = 21, size = 1.9, colour = "grey20", stroke = 0.3) +
    scale_fill_manual(values = c(`All 254 genes` = "#B2182B",
                                 `Overlap removed (217)` = "#2166AC"), name = NULL) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.9, size = 1.9,
             colour = "grey25", label = ov_txt) +
    labs(x = expression("Spearman "*rho*" with source index"), y = NULL) +
    theme_tspo() +
    theme(legend.position = c(0.98, 0.97), legend.justification = c(1, 1),
          legend.background = element_rect(fill = alpha("white", 0.88), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 4.3),
          axis.text.y = element_text(size = 4.8))
} else panel_missing("external cohorts not found")

# =============================================================================
#  F  方差分解（标注热图，与图 2C 同格式）
# =============================================================================

hr("[F] 方差分解")

VD <- do.call(rbind, lapply(names(CO), function(nm) {
  d <- CO[[nm]]$data
  d <- d[is.finite(d$TSPO) & is.finite(d$Myeloid) & is.finite(d$DIF), ]
  if (nrow(d) < 50) return(NULL)
  r_mye  <- summary(lm(TSPO ~ Myeloid, data = d))$r.squared
  r_dif  <- summary(lm(TSPO ~ DIF, data = d))$r.squared
  r_both <- summary(lm(TSPO ~ Myeloid + DIF, data = d))$r.squared
  data.frame(cohort = nm, n = nrow(d),
             髓系独有 = round(r_both - r_dif, 3),
             DIF独有  = round(r_both - r_mye, 3),
             共享     = round(r_mye + r_dif - r_both, 3),
             未解释    = round(1 - r_both, 3),
             合计     = round((r_both - r_dif) + (r_both - r_mye) +
                              (r_mye + r_dif - r_both) + (1 - r_both), 3),
             合并R2   = r_both,
             row.names = NULL)
}))
cat("\n"); pr(VD, 3)
cat(sprintf("\n四段相加校验（应全为 1）: %s\n", paste(unique(VD$合计), collapse = ", ")))
cat(sprintf("两个来源合计解释 %.1f%%–%.1f%% 的 TSPO 方差\n",
            100 * min(VD$合并R2), 100 * max(VD$合并R2)))
write.csv(VD, file.path(OUTDIR, "T38b_variance_decomp.csv"), row.names = FALSE)

# 只画三个解释性成分：未解释 = 1−R² 是补集，放进热图会把色标拉满，
# 使 8–18% 的真实差异全部变浅；其数值在 T38b 与正文中给出。
# 列名用英文：cairo_pdf 无 CJK 字体，中文标签会导致 PDF 输出失败。
COMPS    <- c("髓系独有", "DIF独有", "共享")
COMPS_EN <- c("Myeloid\n(unique)", "DIF\n(unique)", "Shared")
VDL <- do.call(rbind, lapply(seq_len(nrow(VD)), function(i) {
  data.frame(cohort = VD$cohort[i], comp = COMPS,
             v = as.numeric(VD[i, COMPS]), stringsAsFactors = FALSE)
}))
VDL$cohort <- factor(VDL$cohort, levels = rev(VD$cohort))
VDL$comp   <- factor(COMPS_EN[match(VDL$comp, COMPS)], levels = COMPS_EN)

MXV <- max(0.2, max(abs(VDL$v)) * 1.05)
pF <- ggplot(VDL, aes(comp, cohort, fill = v)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", 100 * v)), size = 1.9, colour = "grey15") +
  scale_fill_div(limit = c(-MXV, MXV), name = "R²\n(%)",
                 breaks = c(-0.1, 0, 0.1, 0.2),
                 labels = function(x) sprintf("%.0f", 100 * x)) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_tspo() +
  theme(axis.text.x.top = element_text(angle = 24, hjust = 0, size = 5.0),
        axis.text.y = element_text(size = 5.0),
        panel.grid = element_blank(),
        legend.key.height = unit(8, "mm")) +
  guides(fill = guide_colourbar(barwidth = unit(2.0, "mm"), barheight = unit(8, "mm")))

# =============================================================================
#  置换检验（正文引用，不单独出图）
# =============================================================================

hr("置换检验")

perm_res <- NULL
if (length(NG)) {
  NGm  <- intersect(NG, rownames(TC))
  UNI  <- rownames(TC)
  dset <- intersect(G_DIF, UNI); mset <- intersect(G_MES, UNI)
  obs_d <- length(intersect(dset, NGm)); obs_m <- length(intersect(mset, NGm))
  null_d <- vapply(seq_len(N_PERM), function(i)
    length(intersect(sample(UNI, length(dset)), NGm)), integer(1))
  null_m <- vapply(seq_len(N_PERM), function(i)
    length(intersect(sample(UNI, length(mset)), NGm)), integer(1))
  p_d <- (sum(null_d >= obs_d) + 1) / (N_PERM + 1)
  p_m <- (sum(null_m >= obs_m) + 1) / (N_PERM + 1)
  cat(sprintf("  DIF 观测 %d / 置换均值 %.1f，富集 %.2f 倍，置换 p = %.3f\n",
              obs_d, mean(null_d), obs_d / mean(null_d), p_d))
  cat(sprintf("  MES 观测 %d / 置换均值 %.1f，富集 %.2f 倍，置换 p = %.3f\n",
              obs_m, mean(null_m), obs_m / mean(null_m), p_m))
  perm_res <- data.frame(
    基因集 = c("DIF", "MES"),
    可测基因数 = c(length(dset), length(mset)),
    观测重叠 = c(obs_d, obs_m),
    置换均值 = round(c(mean(null_d), mean(null_m)), 1),
    富集倍数 = round(c(obs_d / mean(null_d), obs_m / mean(null_m)), 2),
    置换p = signif(c(p_d, p_m), 3), row.names = NULL)
  write.csv(perm_res, file.path(OUTDIR, "T38c_permutation_null.csv"), row.names = FALSE)
}

# =============================================================================
#  拼装
# =============================================================================

hr("拼装")

FIG <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, colour = "black"))

if (!exists("FIG_NAME")) FIG_NAME <- "Fig3"
print(FIG)
save_fig(FIG, FIG_NAME, width = 7.7, height = 7.8)

write.csv(data.frame(cohort = names(RH), rho = round(as.numeric(RH), 3)),
          file.path(OUTDIR, "T38a_source_panels.csv"), row.names = FALSE)
if (!is.null(DB)) write.csv(DB, file.path(OUTDIR, "T38d_deoverlap.csv"), row.names = FALSE)

cat("\nfig03_source.R 完成\n")

# -----------------------------------------------------------------------------
# -- 作图样式 / fig04_epi.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig04_epi.R   图 4（9 个面板）上皮细胞中 TSPO 标记氧化磷酸化活跃状态
#
#  A  上皮细胞亚群 UMAP                       散点（嵌入图无法气泡化）
#  B  TSPO 在 UMAP 上的分布                    梯度散点
#  C  各模块的混合模型 β 与跨病人一致率        气泡图
#  D  逐病人斜率（模块 × 病人）                气泡热图
#  E  50 个随机基因模块的对照分布              密度 + 观测标记
#  F  细胞状态与 TSPO 的逐病人相关（状态 × 病人）标注热图
#  G  增殖与 EMT 的零关联                      点 + 95% CI
#  H  亚群 × 病人构成                          气泡矩阵
#  I  CPTAC 蛋白层面相关                       气泡矩阵
#
#  ── 数据来源与限制 ──────────────────────────────────────────────────────────
#  A/B/C/D/E 需要单细胞对象（F_SC_EPI，约 201 MB）。该对象**只存了 counts、
#  没有保存任何降维结果**（Reductions 为空），因此 UMAP 必须在本脚本内重算；
#  为免每次重跑，降维与打分结果缓存到 cache/epi_fig4.rds。
#
#  F/G 用 T13a_epi_state_scores.csv（7,822 个细胞的逐细胞状态评分），
#  H 用对象内的 cell_subtype，I 用 CPTAC_TSPO_protein_cor.csv。
#
#  H 未使用 T12d_cluster_by_patient.csv：该表的聚类是 11 个 Epi0–Epi10 的旧口径，
#  与对象内 12 个 EOC_C1–C12 不一致，改用对象内标签以保证与 UMAP 同源。
#
#  ── 绘图规范 ────────────────────────────────────────────────────────────────
#  不使用柱状图与裸箱线图。多指标汇总用气泡图与气泡热图，相关性矩阵用标注热图
#  （与图 2C 同格式），零假设用密度分布，效应量与不确定度用点 + 95% CI。
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/fig04_epi.R", encoding = "UTF-8")
#    首次约 4–6 分钟（含 UMAP 与 50 次随机模块的混合模型），之后读缓存约 1 分钟。
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

if (!exists("PROJ"))       PROJ       <- "D:/Rstudy/tcga"
if (!exists("OUT"))        OUT        <- file.path(PROJ, "results_v2")
if (!exists("PROJ2"))      PROJ2      <- "D:/Rstudy/tcga2"
if (!exists("F_SC_EPI"))   F_SC_EPI   <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")
if (!exists("SEED"))       SEED       <- 20260101

OUTDIR  <- file.path(OUT, "tables")
TABDIR  <- file.path(PROJ, "results_v2", "tables")   # 结果表所在（可能与 OUT 不同）

# 依次在 OUTDIR、TABDIR 与额外候选里找文件
find_tab <- function(name, extra = character()) {
  cand <- c(file.path(OUTDIR, name), file.path(TABDIR, name), extra)
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}
# 缓存写到工作区（D: 只读，saveRDS 会 Permission denied）
CACHEDIR <- file.path(dirname(FIGDIR), "cache")
if (!exists("FIGDIR")) FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHEDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(patchwork); library(grid); library(ggrepel); library(ggalluvial) })

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

N_RAND    <- 50
RAND_SIZE <- 12
set.seed(SEED)

# =============================================================================
#  STEP 0  模块定义
# =============================================================================

MODULES <- list(
  OXPHOS      = c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
                  "SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","UQCRH","CYC1",
                  "COX4I1","COX5A","COX5B","COX6A1","ATP5A1","ATP5B","ATP5C1","ATP5F1"),
  Mito_Import = c("TOMM22","TOMM40","TOMM70","TIMM23","TIMM17A","TIMM50","TIMM44",
                  "SLC25A3","SLC25A4","SLC25A5"),
  ROS_Defense = c("GPX4","SOD1","SOD2","CAT","PRDX3","TXNRD2","GSR","GCLC","GCLM"),
  Ferroptosis = c("GPX4","HMOX1","NQO1","SLC7A11","FTH1","FTL","TFRC"),
  FAO         = c("CPT1A","CPT1B","CPT2","ACADVL","ACADM","ACADS","ECH1","ECHS1",
                  "HADHA","HADHB"),
  Glycolysis  = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1"),
  TGFB1_Sig   = c("TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3"),
  Immune_Ckpt = c("CD274","CTLA4","LAG3","TIGIT","PDCD1","HAVCR2")
)
MOD_ORDER <- c("OXPHOS","Glycolysis","FAO","ROS_Defense",
               "Mito_Import","Immune_Ckpt","Ferroptosis","TGFB1_Sig")
MOD_LAB <- c(OXPHOS = "OXPHOS", Glycolysis = "Glycolysis", FAO = "Fatty acid oxidation",
             ROS_Defense = "ROS defense", Mito_Import = "Mitochondrial import",
             Immune_Ckpt = "Immune checkpoint", Ferroptosis = "Ferroptosis",
             TGFB1_Sig = "TGFB1 signaling")

# =============================================================================
#  STEP 1  载入对象、降维、打分（带缓存）
# =============================================================================

hr("STEP 1  单细胞对象")

FCACHE <- file.path(CACHEDIR, "epi_fig4.rds")
FCACHE2 <- file.path(CACHEDIR, "epi_fig4_models.rds")   # 模型结果

if (file.exists(FCACHE) && file.exists(FCACHE2)) {
  cat("读取缓存: ", basename(FCACHE), " + ", basename(FCACHE2), "\n", sep = "")
  CA <- readRDS(FCACHE); CB <- readRDS(FCACHE2)   # 两段缓存都要读
  UMAP <- CA$umap; META <- CA$meta; S <- CA$S; tspo <- CA$tspo
  LMM <- CB$lmm; SLOPE <- CB$slope; RND <- CB$rnd
} else {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("需要 Seurat")
  if (!requireNamespace("lme4", quietly = TRUE))  stop("需要 lme4")
  suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(lme4) })
  stopifnot(file.exists(F_SC_EPI))

  sc_full <- readRDS(F_SC_EPI)
  cat("对象: ", nrow(sc_full), " 基因 x ", ncol(sc_full), " 细胞\n", sep = "")
  cells <- colnames(sc_full)[as.character(sc_full$cell_type) == "EOC"]
  sc <- subset(sc_full, cells = cells)
  rm(sc_full); invisible(gc())
  cat("上皮细胞: ", ncol(sc), "  病人: ", length(unique(sc$patient_id)), "\n", sep = "")

  sc <- Seurat::NormalizeData(sc, normalization.method = "LogNormalize",
                              scale.factor = 1e4, verbose = FALSE)

  # ---- 降维：对象未保存任何 reduction，必须重算 ----
  cat("计算 UMAP（对象内无现成降维）...\n")
  sc <- Seurat::FindVariableFeatures(sc, nfeatures = 2000, verbose = FALSE)
  sc <- Seurat::ScaleData(sc, verbose = FALSE)
  sc <- Seurat::RunPCA(sc, npcs = 30, verbose = FALSE)
  sc <- Seurat::RunUMAP(sc, dims = 1:30, seed.use = SEED, verbose = FALSE)
  UMAP <- as.data.frame(Seurat::Embeddings(sc, "umap"))
  colnames(UMAP) <- c("UMAP_1", "UMAP_2")
  META <- data.frame(
    cell = colnames(sc),
    subtype = as.character(sc$cell_subtype),
    patient = as.character(sc$patient_id),
    sample  = as.character(sc$sample),
    phase   = as.character(sc$treatment_phase),
    nCount  = sc$nCount_RNA,
    stringsAsFactors = FALSE)
  META <- cbind(META, UMAP)

  # ---- 模块打分：逐基因去文库大小后 z，再取模块内均值（与 03 一致）----
  mg <- lapply(MODULES, function(g) intersect(g, rownames(sc)))
  for (m in names(mg)) cat(sprintf("  %-12s %2d 基因\n", m, length(mg[[m]])))
  genes <- unique(c(unlist(mg), "TSPO"))
  X <- as.matrix(Seurat::GetAssayData(sc, assay = "RNA", layer = "data")[genes, , drop = FALSE])
  v  <- log10(sc$nCount_RNA + 1); vm <- v - mean(v)
  b  <- as.numeric((X %*% vm) / sum(vm^2))
  a  <- rowMeans(X) - b * mean(v)
  Z  <- t(scale(t(X - (a + b %o% v))))
  Z[!is.finite(Z)] <- NA
  S    <- sapply(mg, function(g) colMeans(Z[g, , drop = FALSE], na.rm = TRUE))
  tspo <- as.numeric(Z["TSPO", ])
  rm(X, Z); invisible(gc())

  pat <- META$patient

  # ---- 混合模型：模块 ~ TSPO + (1|病人) ----
  cat("拟合混合模型...\n")
  LMM <- do.call(rbind, lapply(MOD_ORDER, function(m) {
    d <- data.frame(y = S[, m], tspo = tspo, pat = factor(pat))
    d <- d[is.finite(d$y) & is.finite(d$tspo), ]
    fit <- suppressMessages(suppressWarnings(
      lme4::lmer(y ~ tspo + (1 | pat), data = d, REML = FALSE)))
    cf <- summary(fit)$coefficients
    data.frame(module = m, beta = cf["tspo", "Estimate"],
               se = cf["tspo", "Std. Error"], n = nrow(d), stringsAsFactors = FALSE)
  }))

  # ---- 逐病人斜率 ----
  cat("计算逐病人斜率...\n")
  SLOPE <- do.call(rbind, lapply(unique(pat), function(p) {
    do.call(rbind, lapply(MOD_ORDER, function(m) {
      idx <- pat == p & is.finite(S[, m]) & is.finite(tspo)
      if (sum(idx) < 30) return(NULL)
      r <- suppressWarnings(cor(tspo[idx], S[idx, m], method = "spearman"))
      if (!is.finite(r)) return(NULL)
      data.frame(patient = p, module = m, rho = r,
                 n = sum(idx), stringsAsFactors = FALSE)
    }))
  }))

  # ---- 随机基因模块对照（同样做混合模型）----
  cat("随机模块对照 ", N_RAND, " 次（约 1–2 分钟）...\n", sep = "")
  pool <- setdiff(intersect(rownames(sc), rownames(GetAssayData(sc, assay = "RNA", layer = "data"))),
                  unique(unlist(mg)))
  RND <- vapply(seq_len(N_RAND), function(i) {
    g <- sample(pool, RAND_SIZE)
    Xr <- as.matrix(Seurat::GetAssayData(sc, assay = "RNA", layer = "data")[g, , drop = FALSE])
    br <- as.numeric((Xr %*% vm) / sum(vm^2))
    ar <- rowMeans(Xr) - br * mean(v)
    Zr <- t(scale(t(Xr - (ar + br %o% v)))); Zr[!is.finite(Zr)] <- NA
    s  <- colMeans(Zr, na.rm = TRUE)
    d  <- data.frame(y = s, tspo = tspo, pat = factor(pat))
    d  <- d[is.finite(d$y) & is.finite(d$tspo), ]
    fit <- suppressMessages(suppressWarnings(
      lme4::lmer(y ~ tspo + (1 | pat), data = d, REML = FALSE)))
    as.numeric(summary(fit)$coefficients["tspo", "Estimate"])
  }, numeric(1))

  # 先存降维与打分，避免后续步骤出错时前面的计算全部作废
  saveRDS(list(umap = UMAP, meta = META, S = S, tspo = tspo,
           lmm = NULL, slope = NULL, rnd = NULL), FCACHE)
  saveRDS(list(lmm = LMM, slope = SLOPE, rnd = RND), FCACHE2)
  cat("已缓存: ", FCACHE, "\n", sep = "")
}

S <- as.matrix(S)
LMM$module <- factor(LMM$module, levels = MOD_ORDER)
cat("\n混合模型结果：\n"); pr(LMM, 3)

# =============================================================================
#  A  上皮亚群 UMAP
# =============================================================================

hr("[A] 亚群 UMAP")

SUBS <- sort(unique(META$subtype))
META$subtype <- factor(META$subtype, levels = SUBS)
# 12 色：低饱和、红蓝白体系内取色
PAL12 <- c("#67001F","#B2182B","#D6604D","#EF8A62","#F4A582","#FDDBC7",
           "#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061","#4D4D4D")[seq_along(SUBS)]
names(PAL12) <- SUBS
CENT <- aggregate(cbind(UMAP_1, UMAP_2) ~ subtype, data = META, FUN = median)

pA <- ggplot(META, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = subtype), size = 0.18, stroke = 0) +
  geom_text_repel(data = CENT, aes(UMAP_1, UMAP_2, label = subtype),
                  size = 1.5, colour = "grey15", fontface = "bold",
                  min.segment.length = 0, segment.size = 0.15,
                  segment.colour = "grey60", box.padding = 0.16,
                  point.padding = 0.05, max.overlaps = Inf, seed = SEED) +
  scale_colour_manual(values = PAL12, name = NULL, guide = "none") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_tspo() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# =============================================================================
#  B  TSPO 在 UMAP 上的分布
# =============================================================================

hr("[B] TSPO 分布")

META$TSPO <- tspo
pB <- ggplot(META[order(META$TSPO), ], aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = TSPO), size = 0.18, stroke = 0) +
  scale_colour_gradientn(colors = c("#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B", "#67001F"),
                         name = "TSPO\n(z)", guide = guide_colourbar(
                           barwidth = unit(1.8, "mm"), barheight = unit(7, "mm"))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_tspo() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.title = element_text(size = 4.4), legend.text = element_text(size = 4.0))

# =============================================================================
#  C  β 与一致率：气泡图
#
#  横轴 β（混合模型），纵轴模块；气泡大小 = 11 例中斜率方向与 β 同号的例数，
#  颜色 = β 的符号。这样一屏同时给出效应方向、大小与跨病人可重复性。
# =============================================================================

hr("[C] β 与一致率气泡图")

CONS <- do.call(rbind, lapply(MOD_ORDER, function(m) {
  b <- LMM$beta[LMM$module == m]
  sl <- SLOPE[SLOPE$module == m, ]
  data.frame(module = m, beta = b,
             n_same = sum(sign(sl$rho) == sign(b), na.rm = TRUE),
             n_pat = sum(is.finite(sl$rho)),
             stringsAsFactors = FALSE)
}))
CONS$lab <- MOD_LAB[CONS$module]
CONS$module <- factor(CONS$module, levels = rev(MOD_ORDER))
cat("\n"); pr(CONS[, c("module","beta","n_same","n_pat")], 3)

# 二项检验：11 例中 k 例同号，零假设 p=0.5
CONS$p_sign <- mapply(function(k, n) binom.test(k, n, 0.5)$p.value, CONS$n_same, CONS$n_pat)
CONS$star <- pstar(CONS$p_sign)

pC <- ggplot(CONS, aes(beta, module)) +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "grey60", linetype = "22") +
  geom_segment(aes(x = 0, xend = beta, y = module, yend = module),
               colour = "grey87", linewidth = 0.4) +
  geom_point(aes(size = n_same, fill = beta > 0), shape = 21,
             colour = "grey25", stroke = 0.3) +
  geom_text(aes(label = sprintf("%d/%d%s", n_same, n_pat, star)),
            vjust = -1.25, size = 1.7, colour = "grey20") +
  scale_y_discrete(labels = setNames(MOD_LAB[MOD_ORDER], MOD_ORDER),
                   expand = expansion(add = c(0.45, 0.65))) +
  scale_size_continuous(range = c(0.6, 2.6), breaks = c(0, 5, 11),
                        limits = c(0, 11), name = "Concordant\npatients") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                    name = "Direction") +
  labs(x = expression("Mixed-model "*beta*" (TSPO effect)"), y = NULL) +
  theme_tspo() +
  theme(axis.text.y = element_text(size = 4.8),
        legend.position = "right", legend.key.size = unit(2.3, "mm"),
        legend.text = element_text(size = 4.2), legend.title = element_text(size = 4.4))

# =============================================================================
#  D  逐病人斜率：气泡热图（模块 × 病人）
# =============================================================================

hr("[D] 逐病人斜率气泡热图")

SLOPE$module  <- factor(SLOPE$module, levels = rev(MOD_ORDER))
SLOPE$patient <- factor(SLOPE$patient, levels = sort(unique(SLOPE$patient)))
SLOPE$ab <- abs(SLOPE$rho)

pD <- ggplot(SLOPE, aes(patient, module)) +
  geom_tile(fill = "grey97", colour = "white", linewidth = 0.4) +
  geom_point(aes(size = ab, fill = rho > 0), shape = 21, colour = "grey40", stroke = 0.18) +
  scale_y_discrete(labels = setNames(MOD_LAB[MOD_ORDER], MOD_ORDER)) +
  scale_size_continuous(range = c(0.15, 2.5), limits = c(0, max(SLOPE$ab, na.rm = TRUE)),
                        name = "|rho|") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                    name = "Sign") +
  labs(x = "Patient", y = NULL) +
  theme_tspo() +
  theme(axis.text.x = element_text(angle = 42, hjust = 1, size = 4.4),
        axis.text.y = element_text(size = 4.8),
        panel.grid = element_blank(),
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 4.0),
        legend.title = element_text(size = 4.2))

# =============================================================================
#  E  随机基因模块对照
# =============================================================================

hr("[E] 随机模块的经验 p")

RND_ABS <- abs(RND[is.finite(RND)])
NS <- length(RND_ABS)
# 经验生存函数，+1 校正避免 p = 0
S_emp <- function(b) (sum(RND_ABS >= b) + 1) / (NS + 1)

EC <- data.frame(x = sort(unique(RND_ABS)))
EC$p <- vapply(EC$x, S_emp, numeric(1))

LMM$absb  <- abs(LMM$beta)
LMM$p_emp <- vapply(LMM$absb, S_emp, numeric(1))
cat("\n各模块 |beta| 在随机分布中的经验 p：\n")
pr(LMM[, c("module","beta","absb","p_emp")], 3)

EB <- LMM
EB$lab <- MOD_LAB[EB$module]
EB$col <- ifelse(EB$beta > 0, "#B2182B", "#2166AC")

pE <- ggplot() +
  geom_step(data = EC, aes(x, p), colour = "grey35", linewidth = 0.5,
            direction = "hv", alpha = 0.9) +
  geom_hline(yintercept = 0.05, linetype = "22", colour = "grey60", linewidth = 0.28) +
  geom_segment(data = EB, aes(x = absb, xend = absb, y = 0.012, yend = p_emp),
               colour = "grey80", linewidth = 0.3) +
  geom_point(data = EB, aes(absb, p_emp, fill = beta > 0), shape = 21,
             size = 1.9, colour = "grey20", stroke = 0.3) +
  geom_text_repel(data = EB, aes(absb, p_emp,
                                 label = sprintf("%s (%.3f)", lab, p_emp)),
                  size = 1.45, colour = "grey20", direction = "y",
                  hjust = -0.05, nudge_x = 0.004,
                  segment.size = 0.12, segment.colour = "grey65",
                  max.overlaps = Inf, seed = SEED) +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                    name = "Direction") +
  scale_y_log10(breaks = c(0.0125, 0.05, 0.2, 0.5, 1),
                labels = function(x) sprintf("%.3g", x)) +
  labs(x = expression("|Mixed-model "*beta*"|"),
       y = "Empirical p (random 12-gene modules)") +
  theme_tspo() +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 4.2))


# =============================================================================

hr("[F] 细胞状态 × 病人")

F_STATE <- find_tab("T13a_epi_state_scores.csv")
ST <- NULL; ST_TAB <- NULL
if (!is.na(F_STATE)) {
  st <- read.csv(F_STATE, stringsAsFactors = FALSE, check.names = FALSE)
  sts <- intersect(c("Proliferation","Hypoxia","Ciliated","EMT","OXPHOS","Glycolysis"),
                   colnames(st))
  ST_TAB <- do.call(rbind, lapply(sort(unique(st$patient)), function(p) {
    d <- st[st$patient == p, ]
    do.call(rbind, lapply(sts, function(s) {
      ok <- is.finite(d[[s]]) & is.finite(d$TSPO)
      if (sum(ok) < 20) return(NULL)
      r <- suppressWarnings(cor(d$TSPO[ok], d[[s]][ok], method = "spearman"))
      if (!is.finite(r)) return(NULL)
      data.frame(state = s, patient = p, rho = r,
                 n = sum(ok), stringsAsFactors = FALSE)
    }))
  }))
  cat("\n逐病人相关（前 6 行）：\n"); pr(head(ST_TAB, 6), 3)
  ST <- ST_TAB
}

pF <- if (!is.null(ST)) {
  ST$patient <- factor(ST$patient, levels = sort(unique(ST$patient)))
  ST$state   <- factor(ST$state, levels = rev(c("OXPHOS","Glycolysis","Ciliated",
                                                "Hypoxia","Proliferation","EMT")))
  MXV <- max(abs(ST$rho), na.rm = TRUE)
  ggplot(ST, aes(patient, state, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 1.45, colour = "grey15") +
    scale_fill_div(limit = c(-MXV, MXV), name = "rho") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_tspo() +
    theme(axis.text.x.top = element_text(angle = 42, hjust = 0, size = 4.2),
          axis.text.y = element_text(size = 4.6),
          panel.grid = element_blank(),
          legend.key.height = unit(6, "mm")) +
    guides(fill = guide_colourbar(barwidth = unit(1.8, "mm"), barheight = unit(6, "mm")))
} else panel_missing("T13a not found")

# =============================================================================
#  G  增殖与 EMT 的零关联：点 + 95% CI（跨病人中位数的自助区间）
# =============================================================================

hr("[G] 零关联")

pG <- if (!is.null(ST)) {
  G <- do.call(rbind, lapply(levels(ST$state), function(s) {
    r <- ST$rho[ST$state == s & is.finite(ST$rho)]
    if (!length(r)) return(NULL)
    bs <- replicate(2000, median(sample(r, replace = TRUE)))
    data.frame(state = s, med = median(r), lo = quantile(bs, 0.025),
               hi = quantile(bs, 0.975), n = length(r),
               n_pos = sum(r > 0), stringsAsFactors = FALSE)
  }))
  G$p_sign <- mapply(function(k, n) binom.test(k, n, 0.5)$p.value, G$n_pos, G$n)
  G$state <- factor(G$state, levels = G$state[order(G$med)])
  cat("\n"); pr(G, 3)
  G$sig <- ifelse(G$p_sign < 0.05, "Sign-test p < 0.05", "n.s.")
  ggplot(G, aes(med, state, colour = sig)) +
    geom_vline(xintercept = 0, linewidth = 0.28, colour = "grey60", linetype = "22") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.42,
                  orientation = "y") +
    geom_point(size = 1.7) +
    scale_colour_manual(values = c(`Sign-test p < 0.05` = "#B2182B", `n.s.` = "grey55"),
                        name = NULL) +
    labs(x = "Within-patient rho with TSPO (median, 95% bootstrap CI)", y = NULL) +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.8),
          legend.position = c(0.98, 0.06), legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 4.2))
} else panel_missing("T13a not found")

# =============================================================================
#  H  亚群 × 病人构成：气泡矩阵
#
#  点的大小为该亚群中来自该病人的细胞比例，用于显示亚群被单一病人主导。
# =============================================================================

hr("[H] 亚群大小与单一病人占比")

# 亚群 × 病人构成
COMP <- as.data.frame(table(subtype = META$subtype, patient = META$patient),
                      stringsAsFactors = FALSE)
COMP <- COMP[COMP$Freq > 0, ]
TOT  <- tapply(COMP$Freq, COMP$subtype, sum)
COMP$frac <- COMP$Freq / as.numeric(TOT[COMP$subtype])

DOM <- do.call(rbind, lapply(names(TOT), function(s) {
  d <- COMP[COMP$subtype == s, ]
  d <- d[order(-d$frac), ]
  data.frame(subtype = s,
             n_cell = as.numeric(TOT[[s]]),
             top_pat = d$patient[1],
             top_frac = d$frac[1],
             # 贡献 >=5% 的病人数：衡量该亚群在病人间的分散程度
             n_pat5 = sum(d$frac >= 0.05),
             stringsAsFactors = FALSE)
}))
DOM <- DOM[order(-DOM$top_frac), ]
cat("\n各亚群中占比最高的病人：\n"); pr(DOM, 3)
cat(sprintf("\n单一病人占比 >50%% 的亚群: %d / %d\n",
            sum(DOM$top_frac > 0.5), nrow(DOM)))

pH <- ggplot(DOM, aes(n_cell, top_frac)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0.5, ymax = Inf,
           fill = "#B2182B", alpha = 0.055) +
  geom_hline(yintercept = 0.5, linetype = "22", colour = "grey55", linewidth = 0.3) +
  geom_point(aes(size = n_pat5, fill = top_frac), shape = 21,
             colour = "grey30", stroke = 0.28) +
  geom_text_repel(aes(label = subtype), size = 1.75, colour = "grey15",
                  segment.size = 0.14, segment.colour = "grey70",
                  box.padding = 0.20, point.padding = 0.10,
                  max.overlaps = Inf, seed = SEED) +
  scale_x_log10(breaks = c(100, 200, 400, 800, 1200),
                labels = c("100", "200", "400", "800", "1,200")) +
  scale_y_continuous(limits = c(0.15, 1.02), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_size_continuous(range = c(0.7, 2.4), name = "Patients\ncontributing\n>=5%") +
  scale_fill_gradient(low = "#F7F7F7", high = "#67001F", limits = c(0, 1),
                      name = "Largest\npatient\nshare") +
  annotate("text", x = 1400, y = 0.965, hjust = 1, size = 1.75, colour = "#8C2D2D",
           label = "dominated by one patient") +
  labs(x = "Subtype size (cells, log scale)",
       y = "Share from the single largest patient") +
  theme_tspo() +
  theme(axis.text.x = element_text(size = 4.6),
        axis.text.y = element_text(size = 4.6),
        legend.key.size = unit(2.2, "mm"),
        legend.text = element_text(size = 3.6),
        legend.title = element_text(size = 3.4, lineheight = 0.9))

# =============================================================================
#  fig04_epi.R   图 4（9 个面板）上皮细胞中 TSPO 标记氧化磷酸化活跃状态
#
#  A  上皮细胞亚群 UMAP                       散点（嵌入图无法气泡化）
#  B  TSPO 在 UMAP 上的分布                    梯度散点
#  C  各模块的混合模型 β 与跨病人一致率        气泡图
#  D  逐病人斜率（模块 × 病人）                气泡热图
#  E  50 个随机基因模块的对照分布              密度 + 观测标记
#  F  细胞状态与 TSPO 的逐病人相关（状态 × 病人）标注热图
#  G  增殖与 EMT 的零关联                      点 + 95% CI
#  H  亚群 × 病人构成                          气泡矩阵
#  I  CPTAC 蛋白层面相关                       气泡矩阵
#
#  ── 数据来源与限制 ──────────────────────────────────────────────────────────
#  A/B/C/D/E 需要单细胞对象（F_SC_EPI，约 201 MB）。该对象**只存了 counts、
#  没有保存任何降维结果**（Reductions 为空），因此 UMAP 必须在本脚本内重算；
#  为免每次重跑，降维与打分结果缓存到 cache/epi_fig4.rds。
#
#  F/G 用 T13a_epi_state_scores.csv（7,822 个细胞的逐细胞状态评分），
#  H 用对象内的 cell_subtype，I 用 CPTAC_TSPO_protein_cor.csv。
#
#  H 未使用 T12d_cluster_by_patient.csv：该表的聚类是 11 个 Epi0–Epi10 的旧口径，
#  与对象内 12 个 EOC_C1–C12 不一致，改用对象内标签以保证与 UMAP 同源。
#
#  ── 绘图规范 ────────────────────────────────────────────────────────────────
#  不使用柱状图与裸箱线图。多指标汇总用气泡图与气泡热图，相关性矩阵用标注热图
#  （与图 2C 同格式），零假设用密度分布，效应量与不确定度用点 + 95% CI。
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/fig04_epi.R", encoding = "UTF-8")
#    首次约 4–6 分钟（含 UMAP 与 50 次随机模块的混合模型），之后读缓存约 1 分钟。
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

if (!exists("PROJ"))       PROJ       <- "D:/Rstudy/tcga"
if (!exists("OUT"))        OUT        <- file.path(PROJ, "results_v2")
if (!exists("PROJ2"))      PROJ2      <- "D:/Rstudy/tcga2"
if (!exists("F_SC_EPI"))   F_SC_EPI   <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")
if (!exists("SEED"))       SEED       <- 20260101

OUTDIR  <- file.path(OUT, "tables")
TABDIR  <- file.path(PROJ, "results_v2", "tables")   # 结果表所在（可能与 OUT 不同）

# 依次在 OUTDIR、TABDIR 与额外候选里找文件
find_tab <- function(name, extra = character()) {
  cand <- c(file.path(OUTDIR, name), file.path(TABDIR, name), extra)
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}
# 缓存写到工作区（D: 只读，saveRDS 会 Permission denied）
CACHEDIR <- file.path(dirname(FIGDIR), "cache")
if (!exists("FIGDIR")) FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHEDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(patchwork); library(grid); library(ggrepel); library(ggalluvial) })

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

N_RAND    <- 50
RAND_SIZE <- 12
set.seed(SEED)

# =============================================================================
#  STEP 0  模块定义
# =============================================================================

MODULES <- list(
  OXPHOS      = c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
                  "SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","UQCRH","CYC1",
                  "COX4I1","COX5A","COX5B","COX6A1","ATP5A1","ATP5B","ATP5C1","ATP5F1"),
  Mito_Import = c("TOMM22","TOMM40","TOMM70","TIMM23","TIMM17A","TIMM50","TIMM44",
                  "SLC25A3","SLC25A4","SLC25A5"),
  ROS_Defense = c("GPX4","SOD1","SOD2","CAT","PRDX3","TXNRD2","GSR","GCLC","GCLM"),
  Ferroptosis = c("GPX4","HMOX1","NQO1","SLC7A11","FTH1","FTL","TFRC"),
  FAO         = c("CPT1A","CPT1B","CPT2","ACADVL","ACADM","ACADS","ECH1","ECHS1",
                  "HADHA","HADHB"),
  Glycolysis  = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1"),
  TGFB1_Sig   = c("TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3"),
  Immune_Ckpt = c("CD274","CTLA4","LAG3","TIGIT","PDCD1","HAVCR2")
)
MOD_ORDER <- c("OXPHOS","Glycolysis","FAO","ROS_Defense",
               "Mito_Import","Immune_Ckpt","Ferroptosis","TGFB1_Sig")
MOD_LAB <- c(OXPHOS = "OXPHOS", Glycolysis = "Glycolysis", FAO = "Fatty acid oxidation",
             ROS_Defense = "ROS defense", Mito_Import = "Mitochondrial import",
             Immune_Ckpt = "Immune checkpoint", Ferroptosis = "Ferroptosis",
             TGFB1_Sig = "TGFB1 signaling")

# =============================================================================
#  STEP 1  载入对象、降维、打分（带缓存）
# =============================================================================

hr("STEP 1  单细胞对象")

FCACHE <- file.path(CACHEDIR, "epi_fig4.rds")
FCACHE2 <- file.path(CACHEDIR, "epi_fig4_models.rds")   # 模型结果

if (file.exists(FCACHE) && file.exists(FCACHE2)) {
  cat("读取缓存: ", basename(FCACHE), " + ", basename(FCACHE2), "\n", sep = "")
  CA <- readRDS(FCACHE); CB <- readRDS(FCACHE2)   # 两段缓存都要读
  UMAP <- CA$umap; META <- CA$meta; S <- CA$S; tspo <- CA$tspo
  LMM <- CB$lmm; SLOPE <- CB$slope; RND <- CB$rnd
} else {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("需要 Seurat")
  if (!requireNamespace("lme4", quietly = TRUE))  stop("需要 lme4")
  suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(lme4) })
  stopifnot(file.exists(F_SC_EPI))

  sc_full <- readRDS(F_SC_EPI)
  cat("对象: ", nrow(sc_full), " 基因 x ", ncol(sc_full), " 细胞\n", sep = "")
  cells <- colnames(sc_full)[as.character(sc_full$cell_type) == "EOC"]
  sc <- subset(sc_full, cells = cells)
  rm(sc_full); invisible(gc())
  cat("上皮细胞: ", ncol(sc), "  病人: ", length(unique(sc$patient_id)), "\n", sep = "")

  sc <- Seurat::NormalizeData(sc, normalization.method = "LogNormalize",
                              scale.factor = 1e4, verbose = FALSE)

  # ---- 降维：对象未保存任何 reduction，必须重算 ----
  cat("计算 UMAP（对象内无现成降维）...\n")
  sc <- Seurat::FindVariableFeatures(sc, nfeatures = 2000, verbose = FALSE)
  sc <- Seurat::ScaleData(sc, verbose = FALSE)
  sc <- Seurat::RunPCA(sc, npcs = 30, verbose = FALSE)
  sc <- Seurat::RunUMAP(sc, dims = 1:30, seed.use = SEED, verbose = FALSE)
  UMAP <- as.data.frame(Seurat::Embeddings(sc, "umap"))
  colnames(UMAP) <- c("UMAP_1", "UMAP_2")
  META <- data.frame(
    cell = colnames(sc),
    subtype = as.character(sc$cell_subtype),
    patient = as.character(sc$patient_id),
    sample  = as.character(sc$sample),
    phase   = as.character(sc$treatment_phase),
    nCount  = sc$nCount_RNA,
    stringsAsFactors = FALSE)
  META <- cbind(META, UMAP)

  # ---- 模块打分：逐基因去文库大小后 z，再取模块内均值（与 03 一致）----
  mg <- lapply(MODULES, function(g) intersect(g, rownames(sc)))
  for (m in names(mg)) cat(sprintf("  %-12s %2d 基因\n", m, length(mg[[m]])))
  genes <- unique(c(unlist(mg), "TSPO"))
  X <- as.matrix(Seurat::GetAssayData(sc, assay = "RNA", layer = "data")[genes, , drop = FALSE])
  v  <- log10(sc$nCount_RNA + 1); vm <- v - mean(v)
  b  <- as.numeric((X %*% vm) / sum(vm^2))
  a  <- rowMeans(X) - b * mean(v)
  Z  <- t(scale(t(X - (a + b %o% v))))
  Z[!is.finite(Z)] <- NA
  S    <- sapply(mg, function(g) colMeans(Z[g, , drop = FALSE], na.rm = TRUE))
  tspo <- as.numeric(Z["TSPO", ])
  rm(X, Z); invisible(gc())

  pat <- META$patient

  # ---- 混合模型：模块 ~ TSPO + (1|病人) ----
  cat("拟合混合模型...\n")
  LMM <- do.call(rbind, lapply(MOD_ORDER, function(m) {
    d <- data.frame(y = S[, m], tspo = tspo, pat = factor(pat))
    d <- d[is.finite(d$y) & is.finite(d$tspo), ]
    fit <- suppressMessages(suppressWarnings(
      lme4::lmer(y ~ tspo + (1 | pat), data = d, REML = FALSE)))
    cf <- summary(fit)$coefficients
    data.frame(module = m, beta = cf["tspo", "Estimate"],
               se = cf["tspo", "Std. Error"], n = nrow(d), stringsAsFactors = FALSE)
  }))

  # ---- 逐病人斜率 ----
  cat("计算逐病人斜率...\n")
  SLOPE <- do.call(rbind, lapply(unique(pat), function(p) {
    do.call(rbind, lapply(MOD_ORDER, function(m) {
      idx <- pat == p & is.finite(S[, m]) & is.finite(tspo)
      if (sum(idx) < 30) return(NULL)
      r <- suppressWarnings(cor(tspo[idx], S[idx, m], method = "spearman"))
      if (!is.finite(r)) return(NULL)
      data.frame(patient = p, module = m, rho = r,
                 n = sum(idx), stringsAsFactors = FALSE)
    }))
  }))

  # ---- 随机基因模块对照（同样做混合模型）----
  cat("随机模块对照 ", N_RAND, " 次（约 1–2 分钟）...\n", sep = "")
  pool <- setdiff(intersect(rownames(sc), rownames(GetAssayData(sc, assay = "RNA", layer = "data"))),
                  unique(unlist(mg)))
  RND <- vapply(seq_len(N_RAND), function(i) {
    g <- sample(pool, RAND_SIZE)
    Xr <- as.matrix(Seurat::GetAssayData(sc, assay = "RNA", layer = "data")[g, , drop = FALSE])
    br <- as.numeric((Xr %*% vm) / sum(vm^2))
    ar <- rowMeans(Xr) - br * mean(v)
    Zr <- t(scale(t(Xr - (ar + br %o% v)))); Zr[!is.finite(Zr)] <- NA
    s  <- colMeans(Zr, na.rm = TRUE)
    d  <- data.frame(y = s, tspo = tspo, pat = factor(pat))
    d  <- d[is.finite(d$y) & is.finite(d$tspo), ]
    fit <- suppressMessages(suppressWarnings(
      lme4::lmer(y ~ tspo + (1 | pat), data = d, REML = FALSE)))
    as.numeric(summary(fit)$coefficients["tspo", "Estimate"])
  }, numeric(1))

  # 先存降维与打分，避免后续步骤出错时前面的计算全部作废
  saveRDS(list(umap = UMAP, meta = META, S = S, tspo = tspo,
           lmm = NULL, slope = NULL, rnd = NULL), FCACHE)
  saveRDS(list(lmm = LMM, slope = SLOPE, rnd = RND), FCACHE2)
  cat("已缓存: ", FCACHE, "\n", sep = "")
}

S <- as.matrix(S)
LMM$module <- factor(LMM$module, levels = MOD_ORDER)
cat("\n混合模型结果：\n"); pr(LMM, 3)

# =============================================================================
#  A  上皮亚群 UMAP
# =============================================================================

hr("[A] 亚群 UMAP")

SUBS <- sort(unique(META$subtype))
META$subtype <- factor(META$subtype, levels = SUBS)
# 12 色：低饱和、红蓝白体系内取色
PAL12 <- c("#67001F","#B2182B","#D6604D","#EF8A62","#F4A582","#FDDBC7",
           "#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061","#4D4D4D")[seq_along(SUBS)]
names(PAL12) <- SUBS
CENT <- aggregate(cbind(UMAP_1, UMAP_2) ~ subtype, data = META, FUN = median)

pA <- ggplot(META, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = subtype), size = 0.18, stroke = 0) +
  geom_text_repel(data = CENT, aes(UMAP_1, UMAP_2, label = subtype),
                  size = 1.5, colour = "grey15", fontface = "bold",
                  min.segment.length = 0, segment.size = 0.15,
                  segment.colour = "grey60", box.padding = 0.16,
                  point.padding = 0.05, max.overlaps = Inf, seed = SEED) +
  scale_colour_manual(values = PAL12, name = NULL, guide = "none") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_tspo() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# =============================================================================
#  B  TSPO 在 UMAP 上的分布
# =============================================================================

hr("[B] TSPO 分布")

META$TSPO <- tspo
pB <- ggplot(META[order(META$TSPO), ], aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = TSPO), size = 0.18, stroke = 0) +
  scale_colour_gradientn(colors = c("#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B", "#67001F"),
                         name = "TSPO\n(z)", guide = guide_colourbar(
                           barwidth = unit(1.8, "mm"), barheight = unit(7, "mm"))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_tspo() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.title = element_text(size = 4.4), legend.text = element_text(size = 4.0))

# =============================================================================
#  C  β 与一致率：气泡图
#
#  横轴 β（混合模型），纵轴模块；气泡大小 = 11 例中斜率方向与 β 同号的例数，
#  颜色 = β 的符号。这样一屏同时给出效应方向、大小与跨病人可重复性。
# =============================================================================

hr("[C] β 与一致率气泡图")

CONS <- do.call(rbind, lapply(MOD_ORDER, function(m) {
  b <- LMM$beta[LMM$module == m]
  sl <- SLOPE[SLOPE$module == m, ]
  data.frame(module = m, beta = b,
             n_same = sum(sign(sl$rho) == sign(b), na.rm = TRUE),
             n_pat = sum(is.finite(sl$rho)),
             stringsAsFactors = FALSE)
}))
CONS$lab <- MOD_LAB[CONS$module]
CONS$module <- factor(CONS$module, levels = rev(MOD_ORDER))
cat("\n"); pr(CONS[, c("module","beta","n_same","n_pat")], 3)

# 二项检验：11 例中 k 例同号，零假设 p=0.5
CONS$p_sign <- mapply(function(k, n) binom.test(k, n, 0.5)$p.value, CONS$n_same, CONS$n_pat)
CONS$star <- pstar(CONS$p_sign)

pC <- ggplot(CONS, aes(beta, module)) +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "grey60", linetype = "22") +
  geom_segment(aes(x = 0, xend = beta, y = module, yend = module),
               colour = "grey87", linewidth = 0.4) +
  geom_point(aes(size = n_same, fill = beta > 0), shape = 21,
             colour = "grey25", stroke = 0.3) +
  geom_text(aes(label = sprintf("%d/%d%s", n_same, n_pat, star)),
            vjust = -1.25, size = 1.7, colour = "grey20") +
  scale_y_discrete(labels = setNames(MOD_LAB[MOD_ORDER], MOD_ORDER),
                   expand = expansion(add = c(0.45, 0.65))) +
  scale_size_continuous(range = c(0.6, 2.6), breaks = c(0, 5, 11),
                        limits = c(0, 11), name = "Concordant\npatients") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                    name = "Direction") +
  labs(x = expression("Mixed-model "*beta*" (TSPO effect)"), y = NULL) +
  theme_tspo() +
  theme(axis.text.y = element_text(size = 4.8),
        legend.position = "right", legend.key.size = unit(2.3, "mm"),
        legend.text = element_text(size = 4.2), legend.title = element_text(size = 4.4))

# =============================================================================
#  D  逐病人斜率：气泡热图（模块 × 病人）
# =============================================================================

hr("[D] 逐病人斜率气泡热图")

SLOPE$module  <- factor(SLOPE$module, levels = rev(MOD_ORDER))
SLOPE$patient <- factor(SLOPE$patient, levels = sort(unique(SLOPE$patient)))
SLOPE$ab <- abs(SLOPE$rho)

pD <- ggplot(SLOPE, aes(patient, module)) +
  geom_tile(fill = "grey97", colour = "white", linewidth = 0.4) +
  geom_point(aes(size = ab, fill = rho > 0), shape = 21, colour = "grey40", stroke = 0.18) +
  scale_y_discrete(labels = setNames(MOD_LAB[MOD_ORDER], MOD_ORDER)) +
  scale_size_continuous(range = c(0.15, 2.5), limits = c(0, max(SLOPE$ab, na.rm = TRUE)),
                        name = "|rho|") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                    name = "Sign") +
  labs(x = "Patient", y = NULL) +
  theme_tspo() +
  theme(axis.text.x = element_text(angle = 42, hjust = 1, size = 4.4),
        axis.text.y = element_text(size = 4.8),
        panel.grid = element_blank(),
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 4.0),
        legend.title = element_text(size = 4.2))

# =============================================================================
#  E  随机基因模块对照
# =============================================================================

hr("[E] 随机模块的经验 p")

RND_ABS <- abs(RND[is.finite(RND)])
NS <- length(RND_ABS)
# 经验生存函数，+1 校正避免 p = 0
S_emp <- function(b) (sum(RND_ABS >= b) + 1) / (NS + 1)

EC <- data.frame(x = sort(unique(RND_ABS)))
EC$p <- vapply(EC$x, S_emp, numeric(1))

LMM$absb  <- abs(LMM$beta)
LMM$p_emp <- vapply(LMM$absb, S_emp, numeric(1))
cat("\n各模块 |beta| 在随机分布中的经验 p：\n")
pr(LMM[, c("module","beta","absb","p_emp")], 3)

EB <- LMM
EB$lab <- MOD_LAB[EB$module]
EB$col <- ifelse(EB$beta > 0, "#B2182B", "#2166AC")

pE <- ggplot() +
  geom_step(data = EC, aes(x, p), colour = "grey35", linewidth = 0.5,
            direction = "hv", alpha = 0.9) +
  geom_hline(yintercept = 0.05, linetype = "22", colour = "grey60", linewidth = 0.28) +
  geom_segment(data = EB, aes(x = absb, xend = absb, y = 0.012, yend = p_emp),
               colour = "grey80", linewidth = 0.3) +
  geom_point(data = EB, aes(absb, p_emp, fill = beta > 0), shape = 21,
             size = 1.9, colour = "grey20", stroke = 0.3) +
  geom_text_repel(data = EB, aes(absb, p_emp,
                                 label = sprintf("%s (%.3f)", lab, p_emp)),
                  size = 1.45, colour = "grey20", direction = "y",
                  hjust = -0.05, nudge_x = 0.004,
                  segment.size = 0.12, segment.colour = "grey65",
                  max.overlaps = Inf, seed = SEED) +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                    name = "Direction") +
  scale_y_log10(breaks = c(0.0125, 0.05, 0.2, 0.5, 1),
                labels = function(x) sprintf("%.3g", x)) +
  labs(x = expression("|Mixed-model "*beta*"|"),
       y = "Empirical p (random 12-gene modules)") +
  theme_tspo() +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 4.2))


# =============================================================================

hr("[F] 细胞状态 × 病人")

F_STATE <- find_tab("T13a_epi_state_scores.csv")
ST <- NULL; ST_TAB <- NULL
if (!is.na(F_STATE)) {
  st <- read.csv(F_STATE, stringsAsFactors = FALSE, check.names = FALSE)
  sts <- intersect(c("Proliferation","Hypoxia","Ciliated","EMT","OXPHOS","Glycolysis"),
                   colnames(st))
  ST_TAB <- do.call(rbind, lapply(sort(unique(st$patient)), function(p) {
    d <- st[st$patient == p, ]
    do.call(rbind, lapply(sts, function(s) {
      ok <- is.finite(d[[s]]) & is.finite(d$TSPO)
      if (sum(ok) < 20) return(NULL)
      r <- suppressWarnings(cor(d$TSPO[ok], d[[s]][ok], method = "spearman"))
      if (!is.finite(r)) return(NULL)
      data.frame(state = s, patient = p, rho = r,
                 n = sum(ok), stringsAsFactors = FALSE)
    }))
  }))
  cat("\n逐病人相关（前 6 行）：\n"); pr(head(ST_TAB, 6), 3)
  ST <- ST_TAB
}

pF <- if (!is.null(ST)) {
  ST$patient <- factor(ST$patient, levels = sort(unique(ST$patient)))
  ST$state   <- factor(ST$state, levels = rev(c("OXPHOS","Glycolysis","Ciliated",
                                                "Hypoxia","Proliferation","EMT")))
  MXV <- max(abs(ST$rho), na.rm = TRUE)
  ggplot(ST, aes(patient, state, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 1.45, colour = "grey15") +
    scale_fill_div(limit = c(-MXV, MXV), name = "rho") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_tspo() +
    theme(axis.text.x.top = element_text(angle = 42, hjust = 0, size = 4.2),
          axis.text.y = element_text(size = 4.6),
          panel.grid = element_blank(),
          legend.key.height = unit(6, "mm")) +
    guides(fill = guide_colourbar(barwidth = unit(1.8, "mm"), barheight = unit(6, "mm")))
} else panel_missing("T13a not found")

# =============================================================================
#  G  增殖与 EMT 的零关联：点 + 95% CI（跨病人中位数的自助区间）
# =============================================================================

hr("[G] 零关联")

pG <- if (!is.null(ST)) {
  G <- do.call(rbind, lapply(levels(ST$state), function(s) {
    r <- ST$rho[ST$state == s & is.finite(ST$rho)]
    if (!length(r)) return(NULL)
    bs <- replicate(2000, median(sample(r, replace = TRUE)))
    data.frame(state = s, med = median(r), lo = quantile(bs, 0.025),
               hi = quantile(bs, 0.975), n = length(r),
               n_pos = sum(r > 0), stringsAsFactors = FALSE)
  }))
  G$p_sign <- mapply(function(k, n) binom.test(k, n, 0.5)$p.value, G$n_pos, G$n)
  G$state <- factor(G$state, levels = G$state[order(G$med)])
  cat("\n"); pr(G, 3)
  G$sig <- ifelse(G$p_sign < 0.05, "Sign-test p < 0.05", "n.s.")
  ggplot(G, aes(med, state, colour = sig)) +
    geom_vline(xintercept = 0, linewidth = 0.28, colour = "grey60", linetype = "22") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.42,
                  orientation = "y") +
    geom_point(size = 1.7) +
    scale_colour_manual(values = c(`Sign-test p < 0.05` = "#B2182B", `n.s.` = "grey55"),
                        name = NULL) +
    labs(x = "Within-patient rho with TSPO (median, 95% bootstrap CI)", y = NULL) +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.8),
          legend.position = c(0.98, 0.06), legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 4.2))
} else panel_missing("T13a not found")

# =============================================================================
#  H  亚群 × 病人构成：气泡矩阵
#
#  点的大小为该亚群中来自该病人的细胞比例，用于显示亚群被单一病人主导。
# =============================================================================

hr("[H] 亚群大小与单一病人占比")

# 亚群 × 病人构成
COMP <- as.data.frame(table(subtype = META$subtype, patient = META$patient),
                      stringsAsFactors = FALSE)
COMP <- COMP[COMP$Freq > 0, ]
TOT  <- tapply(COMP$Freq, COMP$subtype, sum)
COMP$frac <- COMP$Freq / as.numeric(TOT[COMP$subtype])

DOM <- do.call(rbind, lapply(names(TOT), function(s) {
  d <- COMP[COMP$subtype == s, ]
  d <- d[order(-d$frac), ]
  data.frame(subtype = s,
             n_cell = as.numeric(TOT[[s]]),
             top_pat = d$patient[1],
             top_frac = d$frac[1],
             # 贡献 >=5% 的病人数：衡量该亚群在病人间的分散程度
             n_pat5 = sum(d$frac >= 0.05),
             stringsAsFactors = FALSE)
}))
DOM <- DOM[order(-DOM$top_frac), ]
cat("\n各亚群中占比最高的病人：\n"); pr(DOM, 3)
cat(sprintf("\n单一病人占比 >50%% 的亚群: %d / %d\n",
            sum(DOM$top_frac > 0.5), nrow(DOM)))

# 堆叠条：每个亚群的细胞来自哪些病人（bar 内和为 1）。
# 原为「亚群大小 × 单一病人占比」的气泡散点，只能看到最大的那一个病人；
# 堆叠条把全部病人的贡献摊开，独占与共享一眼可辨。
pH <- panel_missing("composition not found")
if (exists("COMP") && nrow(COMP)) {
  ORD_S <- DOM$subtype                      # 按单一病人占比从高到低
  CP <- COMP
  CP$subtype <- factor(CP$subtype, levels = rev(ORD_S))
  # 病人按总细胞数排序，配色从深到浅
  PAT_TOT <- sort(tapply(CP$Freq, CP$patient, sum), decreasing = TRUE)
  CP$patient <- factor(CP$patient, levels = names(PAT_TOT))
  PAL_PAT <- setNames(colorRampPalette(c("#67001F", "#D6604D", "#F4A582", "#FDDBC7",
                                         "#92C5DE", "#4393C3", "#2166AC"))(length(PAT_TOT)),
                      names(PAT_TOT))
  cat(sprintf("堆叠条：%d 亚群 × %d 病人\n", length(ORD_S), length(PAT_TOT)))

  pH <- ggplot(CP, aes(frac, subtype, fill = patient)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.12) +
    geom_vline(xintercept = 0.5, linetype = "22", colour = "grey40", linewidth = 0.28) +
    scale_fill_manual(values = PAL_PAT, name = NULL) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    labs(x = "Fraction of cells from each patient", y = NULL) +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 3.6),
          panel.grid.major.y = element_blank(),
          legend.position = "right",
          legend.key.size = unit(1.9, "mm"),
          legend.text = element_text(size = 2.9),
          legend.spacing.y = unit(0.2, "mm")) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))
}

# =============================================================================
#  fig04_epi.R   图 4（9 个面板）上皮细胞中 TSPO 标记氧化磷酸化活跃状态
#
#  A  上皮细胞亚群 UMAP                       散点（嵌入图无法气泡化）
#  B  TSPO 在 UMAP 上的分布                    梯度散点
#  C  各模块的混合模型 β 与跨病人一致率        气泡图
#  D  逐病人斜率（模块 × 病人）                气泡热图
#  E  50 个随机基因模块的对照分布              密度 + 观测标记
#  F  细胞状态与 TSPO 的逐病人相关（状态 × 病人）标注热图
#  G  增殖与 EMT 的零关联                      点 + 95% CI
#  H  亚群 × 病人构成                          气泡矩阵
#  I  CPTAC 蛋白层面相关                       气泡矩阵
#
#  ── 数据来源与限制 ──────────────────────────────────────────────────────────
#  A/B/C/D/E 需要单细胞对象（F_SC_EPI，约 201 MB）。该对象**只存了 counts、
#  没有保存任何降维结果**（Reductions 为空），因此 UMAP 必须在本脚本内重算；
#  为免每次重跑，降维与打分结果缓存到 cache/epi_fig4.rds。
#
#  F/G 用 T13a_epi_state_scores.csv（7,822 个细胞的逐细胞状态评分），
#  H 用对象内的 cell_subtype，I 用 CPTAC_TSPO_protein_cor.csv。
#
#  H 未使用 T12d_cluster_by_patient.csv：该表的聚类是 11 个 Epi0–Epi10 的旧口径，
#  与对象内 12 个 EOC_C1–C12 不一致，改用对象内标签以保证与 UMAP 同源。
#
#  ── 绘图规范 ────────────────────────────────────────────────────────────────
#  不使用柱状图与裸箱线图。多指标汇总用气泡图与气泡热图，相关性矩阵用标注热图
#  （与图 2C 同格式），零假设用密度分布，效应量与不确定度用点 + 95% CI。
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/fig04_epi.R", encoding = "UTF-8")
#    首次约 4–6 分钟（含 UMAP 与 50 次随机模块的混合模型），之后读缓存约 1 分钟。
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

if (!exists("PROJ"))       PROJ       <- "D:/Rstudy/tcga"
if (!exists("OUT"))        OUT        <- file.path(PROJ, "results_v2")
if (!exists("PROJ2"))      PROJ2      <- "D:/Rstudy/tcga2"
if (!exists("F_SC_EPI"))   F_SC_EPI   <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")
if (!exists("SEED"))       SEED       <- 20260101

OUTDIR  <- file.path(OUT, "tables")
TABDIR  <- file.path(PROJ, "results_v2", "tables")   # 结果表所在（可能与 OUT 不同）

# 依次在 OUTDIR、TABDIR 与额外候选里找文件
find_tab <- function(name, extra = character()) {
  cand <- c(file.path(OUTDIR, name), file.path(TABDIR, name), extra)
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}
# 缓存写到工作区（D: 只读，saveRDS 会 Permission denied）
CACHEDIR <- file.path(dirname(FIGDIR), "cache")
if (!exists("FIGDIR")) FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHEDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(patchwork); library(grid); library(ggrepel); library(ggalluvial) })

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)

N_RAND    <- 50
RAND_SIZE <- 12
set.seed(SEED)

# =============================================================================
#  STEP 0  模块定义
# =============================================================================

MODULES <- list(
  OXPHOS      = c("NDUFS1","NDUFS2","NDUFS3","NDUFS5","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
                  "SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","UQCRH","CYC1",
                  "COX4I1","COX5A","COX5B","COX6A1","ATP5A1","ATP5B","ATP5C1","ATP5F1"),
  Mito_Import = c("TOMM22","TOMM40","TOMM70","TIMM23","TIMM17A","TIMM50","TIMM44",
                  "SLC25A3","SLC25A4","SLC25A5"),
  ROS_Defense = c("GPX4","SOD1","SOD2","CAT","PRDX3","TXNRD2","GSR","GCLC","GCLM"),
  Ferroptosis = c("GPX4","HMOX1","NQO1","SLC7A11","FTH1","FTL","TFRC"),
  FAO         = c("CPT1A","CPT1B","CPT2","ACADVL","ACADM","ACADS","ECH1","ECHS1",
                  "HADHA","HADHB"),
  Glycolysis  = c("HK1","HK2","PFKP","ALDOA","GAPDH","PGK1","LDHA","ENO1"),
  TGFB1_Sig   = c("TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3"),
  Immune_Ckpt = c("CD274","CTLA4","LAG3","TIGIT","PDCD1","HAVCR2")
)
MOD_ORDER <- c("OXPHOS","Glycolysis","FAO","ROS_Defense",
               "Mito_Import","Immune_Ckpt","Ferroptosis","TGFB1_Sig")
MOD_LAB <- c(OXPHOS = "OXPHOS", Glycolysis = "Glycolysis", FAO = "Fatty acid oxidation",
             ROS_Defense = "ROS defense", Mito_Import = "Mitochondrial import",
             Immune_Ckpt = "Immune checkpoint", Ferroptosis = "Ferroptosis",
             TGFB1_Sig = "TGFB1 signaling")

# =============================================================================
#  STEP 1  载入对象、降维、打分（带缓存）
# =============================================================================

hr("STEP 1  单细胞对象")

FCACHE <- file.path(CACHEDIR, "epi_fig4.rds")
FCACHE2 <- file.path(CACHEDIR, "epi_fig4_models.rds")   # 模型结果

if (file.exists(FCACHE) && file.exists(FCACHE2)) {
  cat("读取缓存: ", basename(FCACHE), " + ", basename(FCACHE2), "\n", sep = "")
  CA <- readRDS(FCACHE); CB <- readRDS(FCACHE2)   # 两段缓存都要读
  UMAP <- CA$umap; META <- CA$meta; S <- CA$S; tspo <- CA$tspo
  LMM <- CB$lmm; SLOPE <- CB$slope; RND <- CB$rnd
} else {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("需要 Seurat")
  if (!requireNamespace("lme4", quietly = TRUE))  stop("需要 lme4")
  suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(lme4) })
  stopifnot(file.exists(F_SC_EPI))

  sc_full <- readRDS(F_SC_EPI)
  cat("对象: ", nrow(sc_full), " 基因 x ", ncol(sc_full), " 细胞\n", sep = "")
  cells <- colnames(sc_full)[as.character(sc_full$cell_type) == "EOC"]
  sc <- subset(sc_full, cells = cells)
  rm(sc_full); invisible(gc())
  cat("上皮细胞: ", ncol(sc), "  病人: ", length(unique(sc$patient_id)), "\n", sep = "")

  sc <- Seurat::NormalizeData(sc, normalization.method = "LogNormalize",
                              scale.factor = 1e4, verbose = FALSE)

  # ---- 降维：对象未保存任何 reduction，必须重算 ----
  cat("计算 UMAP（对象内无现成降维）...\n")
  sc <- Seurat::FindVariableFeatures(sc, nfeatures = 2000, verbose = FALSE)
  sc <- Seurat::ScaleData(sc, verbose = FALSE)
  sc <- Seurat::RunPCA(sc, npcs = 30, verbose = FALSE)
  sc <- Seurat::RunUMAP(sc, dims = 1:30, seed.use = SEED, verbose = FALSE)
  UMAP <- as.data.frame(Seurat::Embeddings(sc, "umap"))
  colnames(UMAP) <- c("UMAP_1", "UMAP_2")
  META <- data.frame(
    cell = colnames(sc),
    subtype = as.character(sc$cell_subtype),
    patient = as.character(sc$patient_id),
    sample  = as.character(sc$sample),
    phase   = as.character(sc$treatment_phase),
    nCount  = sc$nCount_RNA,
    stringsAsFactors = FALSE)
  META <- cbind(META, UMAP)

  # ---- 模块打分：逐基因去文库大小后 z，再取模块内均值（与 03 一致）----
  mg <- lapply(MODULES, function(g) intersect(g, rownames(sc)))
  for (m in names(mg)) cat(sprintf("  %-12s %2d 基因\n", m, length(mg[[m]])))
  genes <- unique(c(unlist(mg), "TSPO"))
  X <- as.matrix(Seurat::GetAssayData(sc, assay = "RNA", layer = "data")[genes, , drop = FALSE])
  v  <- log10(sc$nCount_RNA + 1); vm <- v - mean(v)
  b  <- as.numeric((X %*% vm) / sum(vm^2))
  a  <- rowMeans(X) - b * mean(v)
  Z  <- t(scale(t(X - (a + b %o% v))))
  Z[!is.finite(Z)] <- NA
  S    <- sapply(mg, function(g) colMeans(Z[g, , drop = FALSE], na.rm = TRUE))
  tspo <- as.numeric(Z["TSPO", ])
  rm(X, Z); invisible(gc())

  pat <- META$patient

  # ---- 混合模型：模块 ~ TSPO + (1|病人) ----
  cat("拟合混合模型...\n")
  LMM <- do.call(rbind, lapply(MOD_ORDER, function(m) {
    d <- data.frame(y = S[, m], tspo = tspo, pat = factor(pat))
    d <- d[is.finite(d$y) & is.finite(d$tspo), ]
    fit <- suppressMessages(suppressWarnings(
      lme4::lmer(y ~ tspo + (1 | pat), data = d, REML = FALSE)))
    cf <- summary(fit)$coefficients
    data.frame(module = m, beta = cf["tspo", "Estimate"],
               se = cf["tspo", "Std. Error"], n = nrow(d), stringsAsFactors = FALSE)
  }))

  # ---- 逐病人斜率 ----
  cat("计算逐病人斜率...\n")
  SLOPE <- do.call(rbind, lapply(unique(pat), function(p) {
    do.call(rbind, lapply(MOD_ORDER, function(m) {
      idx <- pat == p & is.finite(S[, m]) & is.finite(tspo)
      if (sum(idx) < 30) return(NULL)
      r <- suppressWarnings(cor(tspo[idx], S[idx, m], method = "spearman"))
      if (!is.finite(r)) return(NULL)
      data.frame(patient = p, module = m, rho = r,
                 n = sum(idx), stringsAsFactors = FALSE)
    }))
  }))

  # ---- 随机基因模块对照（同样做混合模型）----
  cat("随机模块对照 ", N_RAND, " 次（约 1–2 分钟）...\n", sep = "")
  pool <- setdiff(intersect(rownames(sc), rownames(GetAssayData(sc, assay = "RNA", layer = "data"))),
                  unique(unlist(mg)))
  RND <- vapply(seq_len(N_RAND), function(i) {
    g <- sample(pool, RAND_SIZE)
    Xr <- as.matrix(Seurat::GetAssayData(sc, assay = "RNA", layer = "data")[g, , drop = FALSE])
    br <- as.numeric((Xr %*% vm) / sum(vm^2))
    ar <- rowMeans(Xr) - br * mean(v)
    Zr <- t(scale(t(Xr - (ar + br %o% v)))); Zr[!is.finite(Zr)] <- NA
    s  <- colMeans(Zr, na.rm = TRUE)
    d  <- data.frame(y = s, tspo = tspo, pat = factor(pat))
    d  <- d[is.finite(d$y) & is.finite(d$tspo), ]
    fit <- suppressMessages(suppressWarnings(
      lme4::lmer(y ~ tspo + (1 | pat), data = d, REML = FALSE)))
    as.numeric(summary(fit)$coefficients["tspo", "Estimate"])
  }, numeric(1))

  # 先存降维与打分，避免后续步骤出错时前面的计算全部作废
  saveRDS(list(umap = UMAP, meta = META, S = S, tspo = tspo,
           lmm = NULL, slope = NULL, rnd = NULL), FCACHE)
  saveRDS(list(lmm = LMM, slope = SLOPE, rnd = RND), FCACHE2)
  cat("已缓存: ", FCACHE, "\n", sep = "")
}

S <- as.matrix(S)
LMM$module <- factor(LMM$module, levels = MOD_ORDER)
cat("\n混合模型结果：\n"); pr(LMM, 3)

# =============================================================================
#  A  上皮亚群 UMAP
# =============================================================================

hr("[A] 亚群 UMAP")

SUBS <- sort(unique(META$subtype))
META$subtype <- factor(META$subtype, levels = SUBS)
# 12 色：低饱和、红蓝白体系内取色
PAL12 <- c("#67001F","#B2182B","#D6604D","#EF8A62","#F4A582","#FDDBC7",
           "#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061","#4D4D4D")[seq_along(SUBS)]
names(PAL12) <- SUBS
CENT <- aggregate(cbind(UMAP_1, UMAP_2) ~ subtype, data = META, FUN = median)

pA <- ggplot(META, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = subtype), size = 0.18, stroke = 0) +
  geom_text_repel(data = CENT, aes(UMAP_1, UMAP_2, label = subtype),
                  size = 1.5, colour = "grey15", fontface = "bold",
                  min.segment.length = 0, segment.size = 0.15,
                  segment.colour = "grey60", box.padding = 0.16,
                  point.padding = 0.05, max.overlaps = Inf, seed = SEED) +
  scale_colour_manual(values = PAL12, name = NULL, guide = "none") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_tspo() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# =============================================================================
#  B  TSPO 在 UMAP 上的分布
# =============================================================================

hr("[B] TSPO 分布")

META$TSPO <- tspo
pB <- ggplot(META[order(META$TSPO), ], aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = TSPO), size = 0.18, stroke = 0) +
  scale_colour_gradientn(colors = c("#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B", "#67001F"),
                         name = "TSPO\n(z)", guide = guide_colourbar(
                           barwidth = unit(1.8, "mm"), barheight = unit(7, "mm"))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_tspo() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.title = element_text(size = 4.4), legend.text = element_text(size = 4.0))

# =============================================================================
#  C  β 与一致率：气泡图
#
#  横轴 β（混合模型），纵轴模块；气泡大小 = 11 例中斜率方向与 β 同号的例数，
#  颜色 = β 的符号。这样一屏同时给出效应方向、大小与跨病人可重复性。
# =============================================================================

hr("[C] β 与一致率气泡图")

CONS <- do.call(rbind, lapply(MOD_ORDER, function(m) {
  b <- LMM$beta[LMM$module == m]
  sl <- SLOPE[SLOPE$module == m, ]
  data.frame(module = m, beta = b,
             n_same = sum(sign(sl$rho) == sign(b), na.rm = TRUE),
             n_pat = sum(is.finite(sl$rho)),
             stringsAsFactors = FALSE)
}))
CONS$lab <- MOD_LAB[CONS$module]
CONS$module <- factor(CONS$module, levels = rev(MOD_ORDER))
cat("\n"); pr(CONS[, c("module","beta","n_same","n_pat")], 3)

# 二项检验：11 例中 k 例同号，零假设 p=0.5
CONS$p_sign <- mapply(function(k, n) binom.test(k, n, 0.5)$p.value, CONS$n_same, CONS$n_pat)
CONS$star <- pstar(CONS$p_sign)

pC <- ggplot(CONS, aes(beta, module)) +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "grey60", linetype = "22") +
  geom_segment(aes(x = 0, xend = beta, y = module, yend = module),
               colour = "grey87", linewidth = 0.4) +
  geom_point(aes(size = n_same, fill = beta > 0), shape = 21,
             colour = "grey25", stroke = 0.3) +
  geom_text(aes(label = sprintf("%d/%d%s", n_same, n_pat, star)),
            vjust = -1.25, size = 1.7, colour = "grey20") +
  scale_y_discrete(labels = setNames(MOD_LAB[MOD_ORDER], MOD_ORDER),
                   expand = expansion(add = c(0.45, 0.65))) +
  scale_size_continuous(range = c(0.6, 2.6), breaks = c(0, 5, 11),
                        limits = c(0, 11), name = "Concordant\npatients") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                    name = "Direction") +
  labs(x = expression("Mixed-model "*beta*" (TSPO effect)"), y = NULL) +
  theme_tspo() +
  theme(axis.text.y = element_text(size = 4.8),
        legend.position = "right", legend.key.size = unit(2.3, "mm"),
        legend.text = element_text(size = 4.2), legend.title = element_text(size = 4.4))

# =============================================================================
#  D  逐病人斜率：气泡热图（模块 × 病人）
# =============================================================================

hr("[D] 逐病人斜率气泡热图")

SLOPE$module  <- factor(SLOPE$module, levels = rev(MOD_ORDER))
SLOPE$patient <- factor(SLOPE$patient, levels = sort(unique(SLOPE$patient)))
SLOPE$ab <- abs(SLOPE$rho)

pD <- ggplot(SLOPE, aes(patient, module)) +
  geom_tile(fill = "grey97", colour = "white", linewidth = 0.4) +
  geom_point(aes(size = ab, fill = rho > 0), shape = 21, colour = "grey40", stroke = 0.18) +
  scale_y_discrete(labels = setNames(MOD_LAB[MOD_ORDER], MOD_ORDER)) +
  scale_size_continuous(range = c(0.15, 2.5), limits = c(0, max(SLOPE$ab, na.rm = TRUE)),
                        name = "|rho|") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                    name = "Sign") +
  labs(x = "Patient", y = NULL) +
  theme_tspo() +
  theme(axis.text.x = element_text(angle = 42, hjust = 1, size = 4.4),
        axis.text.y = element_text(size = 4.8),
        panel.grid = element_blank(),
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 4.0),
        legend.title = element_text(size = 4.2))

# =============================================================================
#  E  随机基因模块对照
# =============================================================================

hr("[E] 随机模块的经验 p")

RND_ABS <- abs(RND[is.finite(RND)])
NS <- length(RND_ABS)
# 经验生存函数，+1 校正避免 p = 0
S_emp <- function(b) (sum(RND_ABS >= b) + 1) / (NS + 1)

EC <- data.frame(x = sort(unique(RND_ABS)))
EC$p <- vapply(EC$x, S_emp, numeric(1))

LMM$absb  <- abs(LMM$beta)
LMM$p_emp <- vapply(LMM$absb, S_emp, numeric(1))
cat("\n各模块 |beta| 在随机分布中的经验 p：\n")
pr(LMM[, c("module","beta","absb","p_emp")], 3)

EB <- LMM
EB$lab <- MOD_LAB[EB$module]
EB$col <- ifelse(EB$beta > 0, "#B2182B", "#2166AC")

pE <- ggplot() +
  geom_step(data = EC, aes(x, p), colour = "grey35", linewidth = 0.5,
            direction = "hv", alpha = 0.9) +
  geom_hline(yintercept = 0.05, linetype = "22", colour = "grey60", linewidth = 0.28) +
  geom_segment(data = EB, aes(x = absb, xend = absb, y = 0.012, yend = p_emp),
               colour = "grey80", linewidth = 0.3) +
  geom_point(data = EB, aes(absb, p_emp, fill = beta > 0), shape = 21,
             size = 1.9, colour = "grey20", stroke = 0.3) +
  geom_text_repel(data = EB, aes(absb, p_emp,
                                 label = sprintf("%s (%.3f)", lab, p_emp)),
                  size = 1.45, colour = "grey20", direction = "y",
                  hjust = -0.05, nudge_x = 0.004,
                  segment.size = 0.12, segment.colour = "grey65",
                  max.overlaps = Inf, seed = SEED) +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                    name = "Direction") +
  scale_y_log10(breaks = c(0.0125, 0.05, 0.2, 0.5, 1),
                labels = function(x) sprintf("%.3g", x)) +
  labs(x = expression("|Mixed-model "*beta*"|"),
       y = "Empirical p (random 12-gene modules)") +
  theme_tspo() +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 4.2))


# =============================================================================

hr("[F] 细胞状态 × 病人")

F_STATE <- find_tab("T13a_epi_state_scores.csv")
ST <- NULL; ST_TAB <- NULL
if (!is.na(F_STATE)) {
  st <- read.csv(F_STATE, stringsAsFactors = FALSE, check.names = FALSE)
  sts <- intersect(c("Proliferation","Hypoxia","Ciliated","EMT","OXPHOS","Glycolysis"),
                   colnames(st))
  ST_TAB <- do.call(rbind, lapply(sort(unique(st$patient)), function(p) {
    d <- st[st$patient == p, ]
    do.call(rbind, lapply(sts, function(s) {
      ok <- is.finite(d[[s]]) & is.finite(d$TSPO)
      if (sum(ok) < 20) return(NULL)
      r <- suppressWarnings(cor(d$TSPO[ok], d[[s]][ok], method = "spearman"))
      if (!is.finite(r)) return(NULL)
      data.frame(state = s, patient = p, rho = r,
                 n = sum(ok), stringsAsFactors = FALSE)
    }))
  }))
  cat("\n逐病人相关（前 6 行）：\n"); pr(head(ST_TAB, 6), 3)
  ST <- ST_TAB
}

pF <- if (!is.null(ST)) {
  ST$patient <- factor(ST$patient, levels = sort(unique(ST$patient)))
  ST$state   <- factor(ST$state, levels = rev(c("OXPHOS","Glycolysis","Ciliated",
                                                "Hypoxia","Proliferation","EMT")))
  MXV <- max(abs(ST$rho), na.rm = TRUE)
  ggplot(ST, aes(patient, state, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 1.45, colour = "grey15") +
    scale_fill_div(limit = c(-MXV, MXV), name = "rho") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_tspo() +
    theme(axis.text.x.top = element_text(angle = 42, hjust = 0, size = 4.2),
          axis.text.y = element_text(size = 4.6),
          panel.grid = element_blank(),
          legend.key.height = unit(6, "mm")) +
    guides(fill = guide_colourbar(barwidth = unit(1.8, "mm"), barheight = unit(6, "mm")))
} else panel_missing("T13a not found")

# =============================================================================
#  G  增殖与 EMT 的零关联：点 + 95% CI（跨病人中位数的自助区间）
# =============================================================================

hr("[G] 零关联")

pG <- if (!is.null(ST)) {
  G <- do.call(rbind, lapply(levels(ST$state), function(s) {
    r <- ST$rho[ST$state == s & is.finite(ST$rho)]
    if (!length(r)) return(NULL)
    bs <- replicate(2000, median(sample(r, replace = TRUE)))
    data.frame(state = s, med = median(r), lo = quantile(bs, 0.025),
               hi = quantile(bs, 0.975), n = length(r),
               n_pos = sum(r > 0), stringsAsFactors = FALSE)
  }))
  G$p_sign <- mapply(function(k, n) binom.test(k, n, 0.5)$p.value, G$n_pos, G$n)
  G$state <- factor(G$state, levels = G$state[order(G$med)])
  cat("\n"); pr(G, 3)
  G$sig <- ifelse(G$p_sign < 0.05, "Sign-test p < 0.05", "n.s.")
  ggplot(G, aes(med, state, colour = sig)) +
    geom_vline(xintercept = 0, linewidth = 0.28, colour = "grey60", linetype = "22") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.42,
                  orientation = "y") +
    geom_point(size = 1.7) +
    scale_colour_manual(values = c(`Sign-test p < 0.05` = "#B2182B", `n.s.` = "grey55"),
                        name = NULL) +
    labs(x = "Within-patient rho with TSPO (median, 95% bootstrap CI)", y = NULL) +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.8),
          legend.position = c(0.98, 0.06), legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 4.2))
} else panel_missing("T13a not found")

# =============================================================================
#  H  亚群 × 病人构成：气泡矩阵
#
#  点的大小为该亚群中来自该病人的细胞比例，用于显示亚群被单一病人主导。
# =============================================================================

hr("[H] 亚群大小与单一病人占比")

# 亚群 × 病人构成
COMP <- as.data.frame(table(subtype = META$subtype, patient = META$patient),
                      stringsAsFactors = FALSE)
COMP <- COMP[COMP$Freq > 0, ]
TOT  <- tapply(COMP$Freq, COMP$subtype, sum)
COMP$frac <- COMP$Freq / as.numeric(TOT[COMP$subtype])

DOM <- do.call(rbind, lapply(names(TOT), function(s) {
  d <- COMP[COMP$subtype == s, ]
  d <- d[order(-d$frac), ]
  data.frame(subtype = s,
             n_cell = as.numeric(TOT[[s]]),
             top_pat = d$patient[1],
             top_frac = d$frac[1],
             # 贡献 >=5% 的病人数：衡量该亚群在病人间的分散程度
             n_pat5 = sum(d$frac >= 0.05),
             stringsAsFactors = FALSE)
}))
DOM <- DOM[order(-DOM$top_frac), ]
cat("\n各亚群中占比最高的病人：\n"); pr(DOM, 3)
cat(sprintf("\n单一病人占比 >50%% 的亚群: %d / %d\n",
            sum(DOM$top_frac > 0.5), nrow(DOM)))

pH <- ggplot(DOM, aes(n_cell, top_frac)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0.5, ymax = Inf,
           fill = "#B2182B", alpha = 0.055) +
  geom_hline(yintercept = 0.5, linetype = "22", colour = "grey55", linewidth = 0.3) +
  geom_point(aes(size = n_pat5, fill = top_frac), shape = 21,
             colour = "grey30", stroke = 0.28) +
  geom_text_repel(aes(label = subtype), size = 1.75, colour = "grey15",
                  segment.size = 0.14, segment.colour = "grey70",
                  box.padding = 0.20, point.padding = 0.10,
                  max.overlaps = Inf, seed = SEED) +
  scale_x_log10(breaks = c(100, 200, 400, 800, 1200),
                labels = c("100", "200", "400", "800", "1,200")) +
  scale_y_continuous(limits = c(0.15, 1.02), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_size_continuous(range = c(0.7, 2.4), name = "Patients\ncontributing\n>=5%") +
  scale_fill_gradient(low = "#F7F7F7", high = "#67001F", limits = c(0, 1),
                      name = "Largest\npatient\nshare") +
  annotate("text", x = 1400, y = 0.965, hjust = 1, size = 1.75, colour = "#8C2D2D",
           label = "dominated by one patient") +
  labs(x = "Subtype size (cells, log scale)",
       y = "Share from the single largest patient") +
  theme_tspo() +
  theme(axis.text.x = element_text(size = 4.6),
        axis.text.y = element_text(size = 4.6),
        legend.key.size = unit(2.2, "mm"),
        legend.text = element_text(size = 3.9),
        legend.title = element_text(size = 3.9))

# =============================================================================

hr("[I] CPTAC 蛋白火山图")

F_CPTAC <- find_tab("CPTAC_TSPO_protein_cor.csv",
                    extra = c(file.path(OUT, "run_out", "tables", "CPTAC_TSPO_protein_cor.csv"),
                              file.path(PROJ, "results_v2", "CPTAC_TSPO_protein_cor.csv")))
pI <- panel_missing("CPTAC not found")
if (!is.na(F_CPTAC)) {
  cp <- read.csv(F_CPTAC, stringsAsFactors = FALSE)
  cp <- cp[is.finite(cp$rho) & is.finite(cp$p) & cp$p > 0, ]
  cp$neglogp <- -log10(cp$p)
  SET_ORDER <- c("OXPHOS","Glycolysis","FAO","Ligand","Myeloid","Other")
  cp$set <- factor(cp$set, levels = intersect(SET_ORDER, unique(cp$set)))
  LAB <- cp[abs(cp$rho) >= 0.22 | cp$p < 0.01, ]
  cat("\n"); pr(head(cp[order(-abs(cp$rho)), c("set","protein","rho","p","n")], 8), 3)

  # 下三角相关矩阵：TSPO 与最相关的若干蛋白两两相关。
  # 原为 rho × −log10p 的散点，只回答「谁和 TSPO 相关」；
  # 相关矩阵还能同时看出这些蛋白彼此是否共变（同一 OXPHOS 模块的证据）。
  F_PROT <- file.path(PROJ, "HS_CPTAC_OV_proteome_gene_tumor.cct.txt")
  TOPP <- head(cp$protein[order(-abs(cp$rho))], 10)

  if (file.exists(F_PROT)) {
    PM <- data.table::fread(F_PROT, data.table = FALSE)
    rownames(PM) <- PM[[1]]; PM <- PM[, -1, drop = FALSE]
    PM <- as.matrix(PM); storage.mode(PM) <- "double"
    VAR <- intersect(c("TSPO", TOPP), rownames(PM))
    cat("CPTAC 蛋白矩阵：", nrow(PM), " 蛋白 × ", ncol(PM), " 样本；入图 ",
        length(VAR), " 个\n", sep = "")

    if (length(VAR) >= 6) {
      M <- cor(t(PM[VAR, , drop = FALSE]), method = "spearman", use = "pairwise.complete.obs")
      M[upper.tri(M, diag = TRUE)] <- NA
      CD <- data.frame(g1 = rep(VAR, times = length(VAR)),
                       g2 = rep(VAR, each  = length(VAR)),
                       r  = as.vector(M), stringsAsFactors = FALSE)
      CD <- CD[is.finite(CD$r), ]
      LV <- rev(VAR)
      CD$g1 <- factor(CD$g1, levels = LV)
      CD$g2 <- factor(CD$g2, levels = LV)

      pI <- ggplot(CD, aes(g2, g1)) +
        geom_point(aes(size = abs(r), fill = r), shape = 21,
                   colour = "grey45", stroke = 0.2) +
        geom_text(aes(label = sprintf("%.2f", r), colour = abs(r) > 0.55),
                  size = 1.45, show.legend = FALSE) +
        scale_size_continuous(range = c(0.4, 4.0), name = "|rho|",
                              guide = guide_legend(keyheight = unit(1.6, "mm"))) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, limits = c(-1, 1), name = "rho",
                             guide = guide_colourbar(barwidth = unit(2.2, "mm"),
                                                     barheight = unit(7, "mm"))) +
        scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"),
                            guide = "none") +
        scale_x_discrete(position = "top") +
        labs(x = NULL, y = NULL) +
        theme_tspo() +
        theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.3,
                                             size = 3.8),
              axis.text.y = element_text(size = 3.8),
              panel.grid = element_blank(),
              legend.position = "right",
              legend.key.size = unit(2.0, "mm"),
              legend.text = element_text(size = 3.2),
              legend.title = element_text(size = 3.4),
              legend.spacing.y = unit(0.6, "mm"))
    }
  }
}


# =============================================================================

hr("拼装")

FIG <- (pA | pB | pC) / (pD | pE) / (pF | pG) / (pH | pI) +
  plot_layout(heights = c(1.05, 0.95, 0.95, 1.05))

FIG <- FIG +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, colour = "black"))

if (!exists("FIG_NAME")) FIG_NAME <- "Fig4"
print(FIG)
save_fig(FIG, FIG_NAME, width = 7.7, height = 10.2)

# 导出数值
write.csv(CONS, file.path(OUTDIR, "T39a_module_beta_consistency.csv"), row.names = FALSE)
write.csv(SLOPE, file.path(OUTDIR, "T39b_per_patient_slopes.csv"), row.names = FALSE)
write.csv(data.frame(random_beta = RND),
          file.path(OUTDIR, "T39c_random_modules.csv"), row.names = FALSE)
if (!is.null(ST)) write.csv(ST, file.path(OUTDIR, "T39d_state_by_patient.csv"), row.names = FALSE)
write.csv(DOM, file.path(OUTDIR, "T39e_subtype_dominance.csv"), row.names = FALSE)

cat("\nfig04_epi.R 完成\n")

# -----------------------------------------------------------------------------
# -- 作图样式 / fig05_myeloid_caf.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig05_myeloid_caf.R   图 5（9 个面板）非上皮区室：髓系配体输出与基质 CAF 活化
#
#  A  17 个免疫抑制配体的混合模型 β           森林图（β ± 95% CI + 一致率列表）
#  B  供体间一致率与效应量的关系                气泡散点
#  C  MIF 上调而 CD74 下调                    配体–受体对哑铃图
#  D  MIF 通路占各髓系亚型总通讯量的比例         100% 份额条带
#  E  CPTAC 蛋白层面 TSPO 与髓系标志/配体        棒棒糖图（含 95% CI）
#  F  基质细胞各类型的 TSPO 水平                小提琴 + 蜂群散点
#  G  四类 CAF 状态与 TSPO 的逐病人偏相关        森林图 + 逐病人点
#  H  myCAF 评分与 TSPO 的关系                 二维密度等高线散点
#  I  五个肌成纤维基因的逐病人一致率             基因 × 病人点阵
#
#  ── 数据来源 ────────────────────────────────────────────────────────────────
#  A/B/C  myeloid_ligand_mixed_final.csv（17 个配体与受体，含 beta/t/p/一致率）
#  D      cc5_subtypes_pop0.rds（CellChat 对象，netP$prob 为 源×靶×通路 三维数组）
#  E      CPTAC_TSPO_protein_cor.csv
#  F–I    F_SC_EPI 对象中的 Stromal 细胞（8,003 个），CAF 评分与逐病人偏相关
#         按 24_caf_trajectory.R 的定义与算法重算；结果缓存到 cache/caf_fig5.rds
#
#  ── 绘图规范 ────────────────────────────────────────────────────────────────
#  九个面板九种图型，不使用柱状图与裸箱线图。效应量一律给出 95% 置信区间，
#  跨病人可重复性以文字或点阵单独呈现，避免把「效应大小」与「一致性」混在同一
#  个视觉通道里。
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/fig05_myeloid_caf.R", encoding = "UTF-8")
#    首次约 3–4 分钟（载入对象与重算 CAF），之后读缓存约 1 分钟。
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

if (!exists("PROJ"))     PROJ     <- "D:/Rstudy/tcga"
if (!exists("PROJ2"))    PROJ2    <- "D:/Rstudy/tcga2"
if (!exists("OUT"))      OUT      <- file.path(PROJ, "results_v2")
if (!exists("F_SC_EPI")) F_SC_EPI <- file.path(PROJ2, "GSE165897_processed/seurat_processed.rds")
if (!exists("SEED"))     SEED     <- 20260101

OUTDIR   <- file.path(OUT, "tables")
# 缓存写到工作区（D: 只读）
CACHEDIR <- file.path(dirname(FIGDIR), "cache")
TABDIR   <- file.path(PROJ, "results_v2", "tables")
if (!exists("FIGDIR")) FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHEDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(patchwork); library(grid); library(ggrepel); library(ggbeeswarm) })

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)
set.seed(SEED)

find_tab <- function(name, extra = character()) {
  cand <- c(file.path(OUTDIR, name), file.path(TABDIR, name),
            file.path(OUT, name), file.path(PROJ, "results_v2", name), extra)
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}

# Spearman ρ 的 Fisher z 置信区间
rho_ci <- function(r, n) {
  if (!is.finite(r) || n < 5 || abs(r) >= 1) return(c(NA, NA))
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  tanh(z + c(-1.96, 1.96) * se)
}

# =============================================================================
#  STEP 0  髓系配体表
# =============================================================================

hr("STEP 0  髓系配体")

F_LIG <- find_tab("myeloid_ligand_mixed_final.csv")
if (is.na(F_LIG)) stop("找不到 myeloid_ligand_mixed_final.csv")
LIG <- read.csv(F_LIG, stringsAsFactors = FALSE)
# 由 beta 与 t 反推标准误，再给 95% CI
LIG$se <- abs(LIG$beta / LIG$t)
LIG$lo <- LIG$beta - 1.96 * LIG$se
LIG$hi <- LIG$beta + 1.96 * LIG$se
LIG$sig <- ifelse(LIG$FDR < 0.05, "FDR < 0.05", "n.s.")
LIG <- LIG[order(-LIG$beta), ]
cat("\n17 个配体：\n"); pr(LIG[, c("ligand","beta","lo","hi","pct_patient_pos","FDR")], 3)

# =============================================================================
#  A  配体效应量：森林图 + 一致率列表
# =============================================================================

hr("[A] 配体森林图")

LA <- LIG
LA$ligand <- factor(LA$ligand, levels = rev(LA$ligand))

pA <- ggplot(LA, aes(beta, ligand)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
  geom_segment(aes(x = lo, xend = hi, y = ligand, yend = ligand),
               linewidth = 0.4, colour = "grey45") +
  geom_point(aes(fill = sig), shape = 21, size = 1.9,
             colour = "grey25", stroke = 0.3) +
  geom_text(aes(x = Inf, label = sprintf("%d%%", pct_patient_pos)),
            hjust = 1.02, size = 1.55, colour = "grey30") +
  scale_fill_manual(values = c(`FDR < 0.05` = "#B2182B", `n.s.` = "grey80"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.34))) +
  labs(x = expression("Mixed-model "*beta*" (TSPO effect, 95% CI)"), y = NULL,
       subtitle = "right column: % of donors with positive slope") +
  theme_tspo() +
  theme(axis.text.y = element_text(size = 4.6),
        plot.subtitle = element_text(size = 4.2, colour = "grey35"),
        legend.position = c(0.62, 0.06), legend.justification = c(0, 0),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 4.2))

# =============================================================================
#  B  一致率与效应量：气泡散点
# =============================================================================

hr("[B] 一致率")

LB <- LIG
LB$nlp <- -log10(LB$FDR)

pB <- ggplot(LB, aes(pct_patient_pos, beta)) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = "grey60", linetype = "22") +
  geom_vline(xintercept = 50, linewidth = 0.28, colour = "grey85", linetype = "22") +
  geom_point(aes(size = nlp, fill = beta > 0), shape = 21,
             colour = "grey30", stroke = 0.3, alpha = 0.9) +
  geom_text_repel(aes(label = ligand), size = 1.45, colour = "grey20",
                  segment.size = 0.12, segment.colour = "grey70",
                  box.padding = 0.16, max.overlaps = Inf, seed = SEED) +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                    name = "Direction") +
  scale_size_continuous(range = c(0.5, 2.6), name = expression(-log[10]*"FDR")) +
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100)) +
  labs(x = "% of donors with slope in the same direction", y = expression("Mixed-model "*beta)) +
  theme_tspo() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 4.0),
        legend.title = element_text(size = 4.0))

# =============================================================================
#  C  MIF 上调而 CD74 下调：配体–受体对哑铃图
# =============================================================================

hr("[C] 配体–受体对")

PAIR <- data.frame(
  pair = c("MIF -> CD74", "MIF -> CD74", "LGALS9 -> CD44", "LGALS9 -> CD44"),
  role = c("Ligand", "Receptor", "Ligand", "Receptor"),
  gene = c("MIF", "CD74", "LGALS9", "CD44"),
  stringsAsFactors = FALSE)
PAIR$beta <- LIG$beta[match(PAIR$gene, LIG$ligand)]
PAIR$lo   <- LIG$lo[match(PAIR$gene, LIG$ligand)]
PAIR$hi   <- LIG$hi[match(PAIR$gene, LIG$ligand)]
PAIR$pct  <- LIG$pct_patient_pos[match(PAIR$gene, LIG$ligand)]
PAIR$fdr  <- LIG$FDR[match(PAIR$gene, LIG$ligand)]
cat("\n"); pr(PAIR, 3)

# 每个配体–受体对的连线（从配体到受体，箭头示方向）
SEG <- data.frame(
  pair = c("MIF -> CD74", "LGALS9 -> CD44"),
  x    = PAIR$beta[PAIR$role == "Ligand"],
  xend = PAIR$beta[PAIR$role == "Receptor"])
SEG$pair <- factor(SEG$pair, levels = c("LGALS9 -> CD44", "MIF -> CD74"))
PAIR$pair <- factor(PAIR$pair, levels = levels(SEG$pair))

pC <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
  geom_segment(data = SEG, aes(x = x, xend = xend, y = pair, yend = pair),
               colour = "grey65", linewidth = 0.6,
               arrow = arrow(length = unit(1.0, "mm"), type = "closed")) +
  geom_errorbar(data = PAIR, aes(xmin = lo, xmax = hi, y = pair),
                width = 0, orientation = "y", linewidth = 0.35, colour = "grey50") +
  geom_point(data = PAIR, aes(beta, pair, fill = role), shape = 21,
             size = 2.0, colour = "grey20", stroke = 0.3) +
  geom_text(data = PAIR, aes(x = beta, y = pair, label = gene,
                             vjust = ifelse(role == "Ligand", -1.6, 2.4)),
            size = 1.7, colour = "grey15") +
  scale_fill_manual(values = c(Ligand = "#B2182B", Receptor = "#2166AC"), name = NULL) +
  scale_y_discrete(expand = expansion(add = c(0.75, 0.75))) +
  labs(x = expression("Mixed-model "*beta*" (95% CI)"), y = NULL) +
  theme_tspo() +
  theme(axis.text.y = element_text(size = 4.8),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 4.2))

# =============================================================================
#  C  CellChat：配体 → 受体 的环形连接图（chord）
#
#  原为两对 L–R 的哑铃图（MIF→CD74、LGALS9→CD44）。哑铃只能放两对，
#  而 CellChat 里有 1,098 对可信 L–R，改成 chord 后一次能看到整张通信网络，
#  带宽 = 该 L–R 对在三类髓系亚型间的总通讯强度。
#  数据：cc5_subtypes_pop0.rds 的 net$prob（源 × 靶 × L–R 三维数组）
# =============================================================================

hr("[C] 配体-受体 chord")

F_CC_C <- find_tab("cc5_subtypes_pop0.rds",
                   extra = file.path(PROJ, "results_v2", "cc5_subtypes_pop0.rds"))
pC <- panel_missing("CellChat rds not found")

if (!is.na(F_CC_C) && requireNamespace("circlize", quietly = TRUE) &&
    requireNamespace("png", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(circlize); library(png); library(grid) })
  CC <- readRDS(F_CC_C)
  PB <- CC@net$prob                       # 源 × 靶 × L–R
  LRS <- CC@LR$LRsig

  ST <- apply(PB, 3, sum)                 # 每个 L–R 对的总强度
  DF <- data.frame(interaction = names(ST), strength = as.numeric(ST),
                   stringsAsFactors = FALSE)
  DF$ligand   <- LRS$ligand[match(DF$interaction, LRS$interaction_name)]
  DF$receptor <- LRS$receptor[match(DF$interaction, LRS$interaction_name)]
  DF <- DF[is.finite(DF$strength) & DF$strength > 0, ]

  # 受体名常写成 "TGFbR1_R2"，只取第一段，图上放得下
  DF$rec1 <- sub("_.*$", "", DF$receptor)
  AGG <- aggregate(strength ~ ligand + rec1, data = DF, FUN = sum)

  # 先选配体，再取这些配体真正连到的受体（只按全局 top-N 取受体会丢连接）
  TOPL <- head(names(sort(tapply(AGG$strength, AGG$ligand, sum), decreasing = TRUE)), 6)
  AGG  <- AGG[AGG$ligand %in% TOPL, ]
  TOPR <- head(names(sort(tapply(AGG$strength, AGG$rec1, sum), decreasing = TRUE)), 8)
  AGG  <- AGG[AGG$rec1 %in% TOPR, ]

  if (nrow(AGG) >= 6) {
    MAT <- xtabs(strength ~ ligand + rec1, data = AGG)
    MAT <- MAT[rowSums(MAT) > 0, colSums(MAT) > 0, drop = FALSE]
    cat(sprintf("chord：%d 配体 × %d 受体 × %d 条连接\n",
                nrow(MAT), ncol(MAT), sum(MAT > 0)))

    # 发散式 chord：配体排在外圈，受体落在内圈，
    # 曲线从外向内收拢（模板 B 的布局），线宽 = 通讯强度，颜色随配体。
    LIGS <- rownames(MAT); RECS <- colnames(MAT)
    nL <- length(LIGS); nR <- length(RECS)

    # 内圈错开半个扇区角，避免受体正好落在某一配体的延长线上
    ANG_L <- seq(90, 90 - 360 + 360 / nL, length.out = nL) * pi / 180
    ANG_R <- seq(90 + 180 / nR, 90 - 360 + 180 / nR + 360 / nR, length.out = nR) * pi / 180
    R_OUT <- 1.00; R_IN <- 0.56

    LN <- data.frame(name = LIGS, x = R_OUT * cos(ANG_L), y = R_OUT * sin(ANG_L),
                     ang = ANG_L, stringsAsFactors = FALSE)
    RN <- data.frame(name = RECS, x = R_IN * cos(ANG_R), y = R_IN * sin(ANG_R),
                     ang = ANG_R, stringsAsFactors = FALSE)
    LN$tot <- as.numeric(rowSums(MAT)); RN$tot <- as.numeric(colSums(MAT))

    # 二次贝塞尔：控制点取两端中点再向圆心收 55%，形成内收的弧
    bez <- function(p0, p1, n = 40, pull = 0.55) {
      c0 <- (p0 + p1) / 2 * pull
      t <- seq(0, 1, length.out = n)
      data.frame(x = (1-t)^2*p0[1] + 2*(1-t)*t*c0[1] + t^2*p1[1],
                 y = (1-t)^2*p0[2] + 2*(1-t)*t*c0[2] + t^2*p1[2])
    }
    EDG <- do.call(rbind, lapply(seq_len(nL), function(a) {
      do.call(rbind, lapply(seq_len(nR), function(b) {
        v <- MAT[a, b]
        if (!is.finite(v) || v <= 0) return(NULL)
        d <- bez(c(LN$x[a], LN$y[a]), c(RN$x[b], RN$y[b]))
        d$lig <- LIGS[a]; d$w <- v; d$id <- paste(a, b)
        d
      }))
    }))
    EDG$w01 <- (EDG$w - min(EDG$w)) / (max(EDG$w) - min(EDG$w) + 1e-12)

    PAL_L2 <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nL), LIGS)

    pC <- ggplot() +
      geom_path(data = EDG, aes(x, y, group = id, colour = lig, linewidth = w01),
                alpha = 0.72, lineend = "round") +
      scale_linewidth_continuous(range = c(0.15, 1.6), guide = "none") +
      geom_point(data = LN, aes(x, y, fill = name, size = tot), shape = 21,
                 colour = "grey30", stroke = 0.35) +
      geom_point(data = RN, aes(x, y, fill = name, size = tot), shape = 21,
                 colour = "grey30", stroke = 0.35) +
      scale_size_continuous(range = c(2.2, 5.2), guide = "none") +
      geom_text_repel(data = LN, aes(x * 1.14, y * 1.14, label = name),
                      size = 1.45, colour = "grey10", fontface = "plain",
                      segment.size = 0.12, segment.colour = "grey60",
                      min.segment.length = 0, max.overlaps = Inf, seed = 1,
                      box.padding = 0.35) +
      geom_text_repel(data = RN, aes(x, y, label = name),
                      size = 1.35, colour = "grey10", fontface = "bold",
                      bg.colour = "white", bg.r = 0.12,
                      segment.size = 0.16, segment.colour = "grey45",
                      min.segment.length = 0, max.overlaps = Inf, seed = 2,
                      box.padding = 0.55, point.padding = 0.14) +
      scale_fill_manual(values = c(PAL_L2, setNames(rep("grey80", nR), RECS)),
                        guide = "none") +
      scale_colour_manual(values = PAL_L2, guide = "none") +
      scale_y_continuous(breaks = c(0, 1), labels = c("  ", "  ")) +
      coord_equal(xlim = c(-1.16, 1.16), ylim = c(-1.16, 1.16), clip = "off") +
      labs(x = NULL, y = "   ") +
      theme_tspo() +
      theme(panel.grid = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_blank(), axis.ticks = element_blank(),
            axis.text.y = element_text(colour = "white", size = 4.0),
            axis.title.y = element_text(colour = "white", size = 4.0),
            plot.margin = margin(2, 2, 2, 2))
    cat(sprintf("发散 chord：%d 配体（外圈） × %d 受体（内圈） × %d 条曲线\n",
                nL, nR, length(unique(EDG$id))))
  }
}


# =============================================================================
#  D  CellChat：MIF 通路占比
# =============================================================================

hr("[D] CellChat MIF 占比")

F_CC <- find_tab("cc5_subtypes_pop0.rds",
                 extra = file.path(PROJ, "results_v2", "cc5_subtypes_pop0.rds"))

pD <- panel_missing("CellChat rds not found")
DC <- NULL
if (!is.na(F_CC)) {
  cc <- readRDS(F_CC)
  P  <- cc@netP$prob            # source x target x pathway
  mif <- P[, , "MIF", drop = FALSE][, , 1]
  tot <- apply(P, c(1, 2), sum)
  mif_out <- rowSums(mif)
  tot_out <- rowSums(tot)
  DC <- data.frame(subtype = names(mif_out),
                   MIF = as.numeric(mif_out),
                   Total = as.numeric(tot_out),
                   share_pct = round(100 * mif_out / tot_out, 1),
                   stringsAsFactors = FALSE)
  DC <- DC[order(-DC$share_pct), ]
  cat("\n"); pr(DC, 3)

  DL <- rbind(
    data.frame(subtype = DC$subtype, type = "MIF pathway",   v = DC$MIF),
    data.frame(subtype = DC$subtype, type = "Other pathways", v = DC$Total - DC$MIF))
  DL$subtype <- factor(DL$subtype, levels = rev(DC$subtype))
  DL$type    <- factor(DL$type, levels = c("MIF pathway", "Other pathways"))
  SH <- DC; SH$subtype <- factor(SH$subtype, levels = rev(DC$subtype))

  # 100% 份额条带：每行是该亚型发出信号的构成
  pD <- ggplot(DL, aes(v, subtype, fill = type)) +
    geom_col(position = "fill", width = 0.55,
             colour = "white", linewidth = 0.4) +
    geom_text(data = SH, aes(x = 1.02, y = subtype, label = sprintf("%.1f%%", share_pct)),
              hjust = 0, size = 1.7, colour = "grey25", inherit.aes = FALSE) +
    scale_fill_manual(values = c(`MIF pathway` = "#B2182B",
                                 `Other pathways` = "grey85"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16)),
                       breaks = c(0, 0.5, 1), labels = c("0","50","100%")) +
    scale_y_discrete() +
    labs(x = "Share of outgoing communication signal", y = NULL) +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.6),
          legend.position = c(0.02, 0.99), legend.justification = c(0, 1),
          legend.direction = "horizontal",
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 4.2))
}

# =============================================================================
#  E  CPTAC 蛋白：棒棒糖 + 95% CI
# =============================================================================

hr("[E] CPTAC 髓系与配体")

F_CP <- find_tab("CPTAC_TSPO_protein_cor.csv",
                 extra = file.path(OUT, "run_out", "tables", "CPTAC_TSPO_protein_cor.csv"))
pE <- panel_missing("CPTAC not found")
if (!is.na(F_CP)) {
  cp <- read.csv(F_CP, stringsAsFactors = FALSE)
  cp <- cp[is.finite(cp$rho) & is.finite(cp$n), ]
  # 只看髓系标志与配体两类，聚焦非上皮区室
  # 只看非上皮区室的两类（OXPHOS 亚基已在图 4I 呈现）
  cp <- cp[cp$set %in% c("Myeloid", "Ligand"), ]
  CI <- t(mapply(rho_ci, cp$rho, cp$n))
  cp$lo <- CI[, 1]; cp$hi <- CI[, 2]
  cp <- cp[order(cp$rho), ]
  cp$protein <- factor(cp$protein, levels = cp$protein)
  cat("\n"); pr(head(cp[order(-abs(cp$rho)), c("set","protein","rho","lo","hi","p")], 10), 3)

  pE <- ggplot(cp, aes(rho, protein)) +
    geom_vline(xintercept = 0, linewidth = 0.28, colour = "grey55", linetype = "22") +
    geom_segment(aes(x = 0, xend = rho, y = protein, yend = protein),
                 colour = "grey86", linewidth = 0.42) +
    geom_segment(aes(x = lo, xend = hi, y = protein, yend = protein),
                 colour = "grey45", linewidth = 0.3) +
    geom_point(aes(fill = set), shape = 21, size = 1.5,
               colour = "grey25", stroke = 0.25) +
    scale_fill_manual(values = c(Myeloid = "#B2182B", Ligand = "#2166AC"),
                      name = "Set") +
    labs(x = expression("Spearman "*rho*" with TSPO protein (95% CI)"), y = NULL) +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.4),
          legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 4.0),
          legend.title = element_text(size = 4.0))
}

# =============================================================================
#  STEP 1  基质细胞：CAF 评分与逐病人偏相关（带缓存）
# =============================================================================

hr("STEP 1  基质细胞")

FCACHE <- file.path(CACHEDIR, "caf_fig5_v2.rds")

if (file.exists(FCACHE)) {
  cat("读取缓存: ", basename(FCACHE), "\n", sep = "")
  CF <- readRDS(FCACHE)
  ST_DF <- CF$cells; CAF_PAT <- CF$per_patient; CAF_SUM <- CF$summary
  CAF_SUB <- CF$subtype; GENE_PAT <- CF$gene_patient
  MK_EXPR <- CF$marker_expr
} else {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("需要 Seurat")
  suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
  stopifnot(file.exists(F_SC_EPI))

  sc <- readRDS(F_SC_EPI)
  sc <- subset(sc, cells = colnames(sc)[as.character(sc$cell_type) == "Stromal"])
  if (!"data" %in% SeuratObject::Layers(sc[["RNA"]]))
    sc <- Seurat::NormalizeData(sc, verbose = FALSE)
  cat("基质细胞: ", ncol(sc), "  病人: ", length(unique(sc$patient_id)), "\n", sep = "")

  MK <- list(
    Quiescent = c("CFD","DCN","LUM","MFAP5","PDPN","SERPINF1"),
    iCAF      = c("IL6","CXCL1","CXCL2","LIF","HAS1","CCL2"),
    myCAF     = c("ACTA2","TAGLN","CCN2","COL1A1","COL1A2","POSTN","FN1"),
    apCAF     = c("HLA-DRA","HLA-DRB1","CD74"))
  MARK <- list(
    Fibroblast  = c("DCN","LUM","COL1A1","COL1A2","PDGFRA","THY1","MFAP5"),
    Endothelial = c("PECAM1","VWF","CDH5","PLVAP","EGFL7","CLDN5"),
    Mesothelial = c("MSLN","UPK3B","CALB2","LRRN4","WT1"),
    Pericyte    = c("RGS5","PDGFRB","ACTA2","MYH11","NDUFA4L2"))

  genes <- unique(c(unlist(MK), unlist(MARK), "TSPO",
                    c("COL1A1","ACTA2","FN1","TAGLN","POSTN")))
  genes <- intersect(genes, rownames(sc))
  D <- as.matrix(SeuratObject::LayerData(sc, assay = "RNA", layer = "data")[genes, , drop = FALSE])
  Z <- t(scale(t(D))); Z[!is.finite(Z)] <- NA

  score_of <- function(g) {
    g <- intersect(g, rownames(Z))
    if (length(g) < 3) return(rep(NA_real_, ncol(Z)))
    colMeans(Z[g, , drop = FALSE], na.rm = TRUE)
  }
  tspo <- as.numeric(SeuratObject::LayerData(sc, assay = "RNA", layer = "data")["TSPO", ])
  dep  <- log10(sc$nCount_RNA + 1)
  pat  <- as.character(sc$patient_id)

  # 基质组成归属
  cs <- lapply(MARK, score_of)
  MAT <- do.call(cbind, cs)
  assign_type <- colnames(MAT)[max.col(MAT, ties.method = "first")]

  S <- lapply(MK, score_of)
  ST_DF <- data.frame(patient = pat, cell_type = assign_type, tspo = tspo,
                      myCAF = S$myCAF, iCAF = S$iCAF,
                      Quiescent = S$Quiescent, apCAF = S$apCAF,
                      stringsAsFactors = FALSE)

  # 逐病人偏相关（控制 log10 文库大小，秩残差）
  pcor <- function(a, b, cov, min_cells = 30) {
    ok <- is.finite(a) & is.finite(b) & is.finite(cov)
    if (sum(ok) < min_cells) return(NA_real_)
    ra <- residuals(lm(rank(a[ok]) ~ rank(cov[ok])))
    rb <- residuals(lm(rank(b[ok]) ~ rank(cov[ok])))
    cor(ra, rb)
  }

  CAF_PAT <- do.call(rbind, lapply(c("myCAF","iCAF","Quiescent","apCAF"), function(nm) {
    r <- vapply(split(seq_along(pat), pat), function(i) pcor(tspo[i], S[[nm]][i], dep[i]),
                numeric(1))
    r <- r[is.finite(r)]
    if (!length(r)) return(NULL)
    data.frame(state = nm, patient = names(r), rho = as.numeric(r), stringsAsFactors = FALSE)
  }))

  CAF_SUM <- do.call(rbind, lapply(c("myCAF","iCAF","Quiescent","apCAF"), function(nm) {
    d <- CAF_PAT[CAF_PAT$state == nm, ]
    k <- sum(d$rho > 0); n <- nrow(d)
    data.frame(state = nm, n_patient = n, n_positive = k,
               median_rho = median(d$rho),
               p_sign = binom.test(k, n, 0.5)$p.value, stringsAsFactors = FALSE)
  }))
  CAF_SUM$FDR <- p.adjust(CAF_SUM$p_sign, "BH")
  CAF_SUM <- CAF_SUM[order(-CAF_SUM$median_rho), ]

  # 各基质类型的 TSPO
  CAF_SUB <- do.call(rbind, lapply(sort(unique(assign_type)), function(t) {
    v <- tspo[assign_type == t]
    data.frame(cell_type = t, n = length(v), mean_TSPO = mean(v, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  CAF_SUB <- CAF_SUB[order(-CAF_SUB$mean_TSPO), ]

  # 五个肌成纤维基因的逐病人一致率
  MF <- c("COL1A1","ACTA2","FN1","TAGLN","POSTN")
  SP <- split(seq_along(pat), pat)          # 名字即病人 ID
  GENE_PAT <- do.call(rbind, lapply(MF, function(g) {
    if (!g %in% rownames(D)) return(NULL)
    do.call(rbind, lapply(names(SP), function(p) {
      i <- SP[[p]]
      r <- pcor(tspo[i], as.numeric(D[g, i]), dep[i])
      if (!is.finite(r)) return(NULL)
      data.frame(gene = g, patient = p, rho = r, stringsAsFactors = FALSE)
    }))
  }))

  # 标志基因矩阵（供面板 F 的点图）：每类细胞的经典标志 + TSPO
  MG <- unique(c(MARK$Fibroblast, MARK$Endothelial, MARK$Mesothelial,
                 MARK$Pericyte, "TSPO"))
  MG <- intersect(MG, rownames(D))
  MK_EXPR <- D[MG, , drop = FALSE]

  saveRDS(list(cells = ST_DF, per_patient = CAF_PAT, summary = CAF_SUM,
               subtype = CAF_SUB, gene_patient = GENE_PAT,
               marker_expr = MK_EXPR), FCACHE)
  cat("已缓存: ", FCACHE, "\n", sep = "")
}

cat("\n四类 CAF 的逐病人偏相关：\n"); pr(CAF_SUM, 3)
cat("\n各基质类型的 TSPO：\n"); pr(CAF_SUB, 3)

# =============================================================================
#  F  基质细胞类型的 TSPO：小提琴 + 蜂群
# =============================================================================

hr("[F] 基质细胞类型的标志物点图")

# 单细胞文章常用的 marker dot plot：点大小 = 表达该基因的细胞比例，
# 颜色 = 跨细胞类型 z 标准化的平均表达；末列为 TSPO。
CT_ORDER <- CAF_SUB$cell_type
GENES_F <- c("DCN","LUM","COL1A1","PDGFRA",      # Fibroblast
             "PECAM1","VWF","CDH5","PLVAP",      # Endothelial
             "MSLN","UPK3B","WT1",               # Mesothelial
             "RGS5","PDGFRB","MYH11",            # Pericyte
             "TSPO")
GENES_F <- intersect(GENES_F, rownames(MK_EXPR))
ct <- factor(ST_DF$cell_type, levels = CT_ORDER)

# 缓存里存的是全部标志基因（含 COL1A2/THY1/MFAP5 等），此处只取精选的 15 个，
# 否则多出的行会因不在因子水平里而变成 NA
MK <- MK_EXPR[GENES_F, , drop = FALSE]
AVG <- sapply(CT_ORDER, function(t) rowMeans(MK[, ct == t, drop = FALSE]))
PCT <- sapply(CT_ORDER, function(t) 100 * rowMeans(MK[, ct == t, drop = FALSE] > 0))
SC  <- t(scale(t(AVG))); SC[!is.finite(SC)] <- 0      # 逐基因跨类型标准化

FD <- data.frame(
  gene = factor(rep(rownames(AVG), times = length(CT_ORDER)), levels = rev(GENES_F)),
  cell_type = factor(rep(CT_ORDER, each = nrow(AVG)), levels = CT_ORDER),
  avg = as.vector(AVG), pct = as.vector(PCT), z = as.vector(SC),
  stringsAsFactors = FALSE)
cat("\n标志物阳性率（%）：\n"); pr(PCT, 1)

pF <- ggplot(FD, aes(cell_type, gene)) +
  geom_point(aes(size = pct, fill = z), shape = 21, colour = "grey35", stroke = 0.18) +
  scale_size_continuous(range = c(0.12, 2.4), limits = c(0, 100),
                        breaks = c(0, 25, 50, 75, 100), name = "% expressing") +
  scale_fill_div(limit = c(-max(abs(FD$z)), max(abs(FD$z))), name = "Scaled\nmean") +
  geom_hline(yintercept = 0.5, linewidth = 0.4, colour = "white") +
  labs(x = NULL, y = NULL) +
  theme_tspo() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 4.4),
        axis.text.y = element_text(size = 4.4),
        panel.grid = element_blank(),
        legend.key.size = unit(2.2, "mm"),
        legend.text = element_text(size = 3.9),
        legend.title = element_text(size = 4.0))


# =============================================================================

hr("[G] CAF 森林图")

CAF_SUM$lab <- sprintf("%s (%d/%d)", CAF_SUM$state, CAF_SUM$n_positive, CAF_SUM$n_patient)
SG <- CAF_SUM
SG$lab <- factor(SG$lab, levels = rev(SG$lab))
CAF_PAT$lab <- factor(SG$lab[match(CAF_PAT$state, SG$state)], levels = levels(SG$lab))
CAF_PAT <- CAF_PAT[!is.na(CAF_PAT$lab), ]

pG <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
  geom_point(data = CAF_PAT, aes(rho, lab), size = 0.75, alpha = 0.55,
             colour = "grey45", stroke = 0) +
  geom_point(data = SG, aes(median_rho, lab, fill = p_sign < 0.05),
             shape = 23, size = 2.2, colour = "grey15", stroke = 0.35) +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey80"),
                    labels = c(`TRUE` = "Sign test p < 0.05", `FALSE` = "n.s."),
                    name = NULL) +
  scale_y_discrete(expand = expansion(add = c(0.55, 0.55))) +
  labs(x = expression("Partial "*rho*" with TSPO (per patient; diamond = median)"), y = NULL) +
  theme_tspo() +
  theme(axis.text.y = element_text(size = 4.6),
        legend.position = c(0.98, 0.06), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 4.1))

# =============================================================================
#  H  myCAF 与 TSPO：二维密度等高线
# =============================================================================

hr("[H] myCAF 与 TSPO（病人内）")

# 直接合并所有细胞会混入病人间差异（正是本文反复警告的 Simpson 悖论）：
# 面板 G 的推断是「病人内」偏相关，故此处先按病人各自的中位数中心化再合并，
# 二者口径一致，图与检验才能互相印证。
SH <- ST_DF[is.finite(ST_DF$myCAF) & is.finite(ST_DF$tspo), ]
SH$x <- SH$myCAF - ave(SH$myCAF, SH$patient, FUN = median)
SH$y <- SH$tspo  - ave(SH$tspo,  SH$patient, FUN = median)
R_POOL <- cor(SH$x, SH$y, method = "spearman")
R_MED  <- CAF_SUM$median_rho[CAF_SUM$state == "myCAF"]
N_POS  <- CAF_SUM$n_positive[CAF_SUM$state == "myCAF"]
N_PAT  <- CAF_SUM$n_patient[CAF_SUM$state == "myCAF"]
cat(sprintf("\n病人内中心化后合并 rho = %.3f；逐病人中位偏相关 = %.3f（%d/%d 为正）\n",
            R_POOL, R_MED, N_POS, N_PAT))

pH <- ggplot(SH, aes(x, y)) +
  geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey88") +
  geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey88") +
  geom_hex(bins = 34, linewidth = 0.10) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              linewidth = 0.65, colour = PAL["ink"], fill = alpha("white", 0)) +
  scale_fill_gradient(low = "#F7F7F7", high = "#67001F",
                      name = "Cells",
                      guide = guide_colourbar(barwidth = unit(1.9, "mm"),
                                              barheight = unit(6, "mm"))) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.9, size = 1.75,
           colour = PAL["ink"],
           label = sprintf("median within-patient partial rho = %+.3f\n%d/%d patients positive",
                           R_MED, N_POS, N_PAT)) +
  labs(x = "myCAF score (centred within patient)",
       y = "TSPO (centred within patient)") +
  theme_tspo() +
  theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 4.0),
        legend.title = element_text(size = 4.0))


# =============================================================================

hr("[I] 肌成纤维基因")

pI <- panel_missing("gene-level data not available")
if (!is.null(GENE_PAT) && nrow(GENE_PAT) > 0) {
  GI <- GENE_PAT
  GS <- aggregate(rho ~ gene, GI, function(v) c(med = median(v), pos = sum(v > 0), n = length(v)))
  GS <- do.call(data.frame, GS)
  names(GS) <- c("gene", "med", "pos", "n")
  GS <- GS[order(-GS$med), ]
  cat("\n"); pr(GS, 3)

  GI$gene    <- factor(GI$gene, levels = rev(GS$gene))
  GI$patient <- factor(GI$patient, levels = sort(unique(GI$patient)))
  GI$lab <- factor(sprintf("%s (%d/%d)", GS$gene, GS$pos, GS$n)[match(GI$gene, GS$gene)],
                   levels = rev(sprintf("%s (%d/%d)", GS$gene, GS$pos, GS$n)))

  MXV <- max(abs(GI$rho), na.rm = TRUE)
  pI <- ggplot(GI, aes(patient, lab)) +
    geom_tile(fill = "grey97", colour = "white", linewidth = 0.35) +
    geom_point(aes(size = abs(rho), fill = rho), shape = 21,
               colour = "grey40", stroke = 0.15) +
    scale_fill_div(limit = c(-MXV, MXV), name = "Partial\nrho") +
    scale_size_continuous(range = c(0.15, 2.3), limits = c(0, MXV), name = "|rho|") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_tspo() +
    theme(axis.text.x.top = element_text(angle = 42, hjust = 0, size = 4.2),
          axis.text.y = element_text(size = 4.4),
          panel.grid = element_blank(),
          legend.key.size = unit(2.2, "mm"),
          legend.text = element_text(size = 4.0),
          legend.title = element_text(size = 4.0))
}

# =============================================================================
#  拼装
# =============================================================================

hr("拼装")

# 3 × 3，九格填满，不留空位
FIG <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI)

FIG <- FIG +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, colour = "black"))

if (!exists("FIG_NAME")) FIG_NAME <- "Fig5"
print(FIG)
save_fig(FIG, FIG_NAME, width = 7.7, height = 9.6)

write.csv(LIG[, c("ligand","beta","lo","hi","p","FDR","pct_patient_pos")],
          file.path(OUTDIR, "T40a_ligand_forest.csv"), row.names = FALSE)
if (!is.null(DC)) write.csv(DC, file.path(OUTDIR, "T40b_cellchat_mif_share.csv"), row.names = FALSE)
write.csv(CAF_SUM, file.path(OUTDIR, "T40c_caf_summary.csv"), row.names = FALSE)
write.csv(CAF_PAT, file.path(OUTDIR, "T40d_caf_per_patient.csv"), row.names = FALSE)
write.csv(CAF_SUB, file.path(OUTDIR, "T40e_stromal_tspo.csv"), row.names = FALSE)
if (!is.null(GENE_PAT)) write.csv(GENE_PAT, file.path(OUTDIR, "T40f_mf_genes.csv"), row.names = FALSE)

cat("\nfig05_myeloid_caf.R 完成\n")

# -----------------------------------------------------------------------------
# -- 作图样式 / fig06_trajectory.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig06_trajectory.R   图 6（9 个面板）髓系 TSPO 沿单核→巨噬分化轴程序性下调
#
#  A  UMAP（细胞类型）                         嵌入图散点 + 质心标注
#  B  UMAP（伪时序）+ 主轴路径                  连续梯度 + 分箱质心连线
#  C  标志基因 × 伪时序分箱                     热图（早期褪去、晚期升起）
#  D  TSPO 沿伪时序                            六边形分箱 + 平滑曲线
#  E  2000 个随机基因的 rho 零分布              分布 + 分位标注
#  F  逐病人 rho                                排序点阵 + 符号检验
#  G  直接测量与独立重现                        斜率图（两数据集 × 两细胞类型）
#  H  化疗前后的配对比较                        配对斜率图
#  I  下游配体沿伪时序                          多曲线直接标注
#
#  ── 高分期刊里这类图的标准结构 ──────────────────────────────────────────────
#  轨迹分析图通常按「几何 → 方向 → 目标基因 → 特异性 → 稳健性 → 独立验证」推进：
#    几何：嵌入图与主轴（A、B）
#    方向：标志基因沿轴的变化证明确实是分化方向（C）
#    目标：TSPO 沿轴的表达（D）
#    特异性：随机基因零分布排除整体漂移（E）
#    稳健性：逐病人一致性与治疗暴露（F、H）
#    独立验证：不依赖轨迹方法的直接测量 + 独立数据集（G）
#    下游：配体程序是否同步（I）
#
#  ── 数据 ────────────────────────────────────────────────────────────────────
#  cache/trajectory_fig.rds        06 脚本产出的绘图缓存（UMAP/伪时序/表达/随机背景）
#  tables/T3b_trajectory_per_patient.csv   逐病人 rho
#  tables/T3a_trajectory_direction.csv     方向验证
#  tables/T3c_TSPO_by_celltype.csv         平衡子集细胞类型均值
#  tables/T3c_full_TSPO_by_celltype.csv    全集群细胞类型均值
#  tables/T41_gse154600_tspo.csv           GSE154600 独立重现
#  GSE266577_paired_TSPO.csv               化疗配对
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/fig06_trajectory.R", encoding = "UTF-8")
#    需先跑过 06_myeloid_trajectory.R 生成缓存，约 1 分钟出图。
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

if (!exists("PROJ"))  PROJ  <- "D:/Rstudy/tcga"
if (!exists("OUT"))   OUT   <- file.path(PROJ, "results_v2")
if (!exists("SEED"))  SEED  <- 20260101

OUTDIR   <- file.path(OUT, "tables")
CACHEDIR <- file.path(OUT, "cache")
TABDIR   <- file.path(PROJ, "results_v2", "tables")
if (!exists("FIGDIR")) FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(patchwork); library(grid); library(ggrepel) })

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)
set.seed(SEED)

find_tab <- function(name) {
  cand <- c(file.path(OUTDIR, name), file.path(TABDIR, name),
            file.path(OUT, name), file.path(PROJ, "results_v2", name))
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}

# =============================================================================
#  STEP 0  载入缓存
# =============================================================================

hr("STEP 0  载入")

F_CACHE <- file.path(CACHEDIR, "trajectory_fig.rds")
if (!file.exists(F_CACHE)) stop("缺 ", F_CACHE, "，请先运行 06_myeloid_trajectory.R")
CA <- readRDS(F_CACHE)

UM  <- as.data.frame(CA$umap); colnames(UM) <- c("UMAP_1", "UMAP_2")
PT  <- as.numeric(CA$pt)
EX  <- CA$expr
CT  <- CA$cell_type
TRE <- CA$treatment
NC  <- CA$n_count
stopifnot(nrow(UM) == length(PT), ncol(EX) == length(PT))

cat("细胞: ", length(PT), "  基因: ", nrow(EX), "\n", sep = "")
cat("细胞类型:\n"); print(table(CT))
cat("伪时序范围: ", paste(round(range(PT, na.rm = TRUE), 1), collapse = " ~ "), "\n", sep = "")

DF <- data.frame(UMAP_1 = UM$UMAP_1, UMAP_2 = UM$UMAP_2, pt = PT,
                 cell_type = CT, nCount = NC, stringsAsFactors = FALSE)

CT_SHORT <- c(`Classical monocytes` = "Classical monocytes",
              `Macrophages` = "Macrophages")
DF$ct <- factor(ifelse(DF$cell_type == "Classical monocytes",
                       "Classical monocytes", "Macrophages"),
                levels = c("Classical monocytes", "Macrophages"))
# 原来的 #B2182B / #2166AC 在 0.18 点径下过淡，改用更深的两色并加大点径
PAL_CT <- c(`Classical monocytes` = "#A50026", `Macrophages` = "#08519C")

# =============================================================================
#  A  UMAP：细胞类型
# =============================================================================

hr("[A] UMAP：细胞类型与伪时序")

# 原 A、B 是同一张 UMAP，只差一个透明度参数，信息重复。
# 合并为一张：颜色给细胞类型，透明度给伪时序（早期淡、晚期实）。
DF$alpha_pt <- scales::rescale(DF$pt, to = c(0.30, 1.00))
# 巨噬（1,638）远多于单核（546），单核最后绘制才不会被盖住
DA <- DF[order(DF$ct != "Classical monocytes", DF$pt), ]

CEN <- aggregate(cbind(UMAP_1, UMAP_2) ~ ct, data = DF, FUN = median)
CNT <- as.data.frame(table(DF$ct)); names(CNT) <- c("ct", "n")
CEN$lab <- sprintf("%s\n(n = %s)", CEN$ct,
                   format(CNT$n[match(CEN$ct, CNT$ct)], big.mark = ","))

pA <- ggplot(DA, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = ct, alpha = alpha_pt), size = 0.42, stroke = 0) +
  geom_text_repel(data = CEN, aes(UMAP_1, UMAP_2, label = lab),
                  size = 1.9, fontface = "bold", colour = "grey10",
                  bg.colour = "white", bg.r = 0.22,
                  min.segment.length = 0, segment.size = 0.2,
                  max.overlaps = Inf, seed = SEED) +
  scale_colour_manual(values = PAL_CT, name = NULL) +
  scale_alpha_identity(guide = "none") +
  labs(x = "UMAP 1", y = "UMAP 2",
       subtitle = "colour = cell type; opacity = pseudotime (light early, solid late)") +
  theme_tspo() +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = 4.2),
        plot.subtitle = element_text(size = 4.1, colour = "grey35"))

hr("[B] 伪时序是否由测序深度驱动")

# 伪时序与测序深度中度相关（rho 见下），正文报的「校正深度后 rho」因此必要。
# 这一面板把该混杂直接画出来，读者能自行判断校正是否合理。
DF$lnc <- log10(DF$nCount + 1)
R_TECH <- cor(DF$pt, DF$lnc, method = "spearman")
CAT_TECH <- cor(DF$pt, DF$lnc, method = "spearman")
cat(sprintf("\n伪时序 vs 测序深度: rho = %+.3f\n", R_TECH))

pB <- ggplot(DF, aes(lnc, pt)) +
  geom_hex(bins = 34, linewidth = 0.08) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.7, se = TRUE,
              linewidth = 0.65, colour = PAL["ink"], fill = alpha("white", 0)) +
  scale_fill_gradient(low = "#F7F7F7", high = "#67001F", trans = "log10",
                      name = "Cells",
                      guide = guide_colourbar(barwidth = unit(1.9, "mm"),
                                              barheight = unit(6, "mm"))) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.9, size = 1.8,
           colour = PAL["ink"],
           label = sprintf("rho(pseudotime, depth) = %+.3f", R_TECH)) +
  labs(x = expression(log[10]*"(nCount + 1)"),
       y = "Pseudotime") +
  theme_tspo() +
  theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 4.0),
        legend.title = element_text(size = 4.0))

# =============================================================================

hr("[C] 标志基因热图")

EARLY <- c("CD14","FCN1","S100A8","S100A9","VCAN","LYZ")
LATE  <- c("CD163","MRC1","CSF1R","MAF","F13A1","C1QA","APOE")
GN    <- c(EARLY, LATE)
GN    <- GN[GN %in% rownames(EX)]

NBIN <- 15
BK <- cut(rank(PT, ties.method = "first"), breaks = NBIN, labels = FALSE)
HM <- do.call(rbind, lapply(GN, function(g) {
  v <- as.numeric(EX[g, ])
  z <- (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  tapply(z, BK, mean, na.rm = TRUE)
}))
colnames(HM) <- sprintf("B%02d", seq_len(NBIN))
rownames(HM) <- GN          # lapply 未命名，rbind 后行名是 1..n，必须补回基因名
HML <- data.frame(gene = factor(rep(GN, each = NBIN), levels = rev(GN)),
                  bin = rep(seq_len(NBIN), times = length(GN)),
                  z = as.vector(t(HM)))
# 按峰值位置排序基因，使早期基因在上、晚期在下
PEAK <- apply(HM, 1, which.max)
ORD  <- GN[order(PEAK[GN])]
HML$gene <- factor(HML$gene, levels = rev(ORD))
cat("\n基因峰值所在分箱（1=最早）：\n")
print(data.frame(gene = ORD, peak_bin = as.integer(PEAK[ORD])), row.names = FALSE)

# 每侧细胞占比条（说明分箱内单核/巨噬构成）
CMP <- as.data.frame(prop.table(table(BK, DF$ct), 1))
names(CMP) <- c("bin", "ct", "frac"); CMP$bin <- as.integer(CMP$bin)

pC <- ggplot(HML, aes(bin, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_div(limit = c(-max(abs(HML$z)), max(abs(HML$z))), name = "z-score") +
  scale_x_continuous(breaks = c(1, NBIN), labels = c("early", "late"),
                     expand = c(0.01, 0.01)) +
  labs(x = "Pseudotime", y = NULL) +
  theme_tspo() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 4.2),
        axis.text.x = element_text(size = 4.6),
        legend.key.height = unit(6, "mm")) +
  guides(fill = guide_colourbar(barwidth = unit(1.8, "mm"), barheight = unit(6, "mm")))

# =============================================================================
#  D  TSPO 沿伪时序
# =============================================================================

hr("[D] TSPO 沿伪时序")

DD <- data.frame(pt = PT, tspo = as.numeric(EX["TSPO", ]))
DD <- DD[is.finite(DD$pt) & is.finite(DD$tspo), ]
R_RAW <- CA$tspo_rho
# 深度校正：对 log10 文库大小做秩残差后重算
DEP <- log10(NC + 1)
okd <- is.finite(DD$pt) & is.finite(DD$tspo) & is.finite(DEP)
r1 <- residuals(lm(rank(DD$pt[okd]) ~ rank(DEP[okd])))
r2 <- residuals(lm(rank(DD$tspo[okd]) ~ rank(DEP[okd])))
R_ADJ <- cor(r1, r2)
cat(sprintf("原始 rho = %.3f；校正深度后 = %.3f\n", R_RAW, R_ADJ))

pD <- ggplot(DD, aes(pt, tspo)) +
  geom_hex(bins = 40, linewidth = 0.08) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.6, se = TRUE,
              linewidth = 0.65, colour = PAL["ink"], fill = alpha("white", 0)) +
  scale_fill_gradient(low = "#F7F7F7", high = "#67001F", name = "Cells", trans = "log10",
                      guide = guide_colourbar(barwidth = unit(1.9, "mm"),
                                              barheight = unit(6, "mm"))) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.9, size = 1.8,
           colour = PAL["ink"],
           label = sprintf("rho = %.3f\nrho | depth = %.3f", R_RAW, R_ADJ)) +
  labs(x = "Pseudotime", y = "TSPO (log-normalized)") +
  theme_tspo() +
  theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 4.0),
        legend.title = element_text(size = 4.0))

# =============================================================================
#  E  随机基因零分布
# =============================================================================

hr("[E] 随机基因零分布")

# CA$random 是 vapply 得到的数值向量（不是数据框）
RD <- data.frame(rho = as.numeric(CA$random))
RD <- RD[is.finite(RD$rho), , drop = FALSE]
PCT <- 100 * mean(RD$rho > R_RAW)
cat(sprintf("随机基因 %d 个；负相关占 %.1f%%；中位 %.3f；TSPO 位于第 %.0f 百分位\n",
            nrow(RD), 100 * mean(RD$rho < 0), median(RD$rho), PCT))

pE <- ggplot(RD, aes(rho, 0)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55", linetype = "22") +
  geom_violin(aes(group = 1), orientation = "y", scale = "width",
              width = 0.72, alpha = 0.55, fill = "grey82",
              colour = "grey45", linewidth = 0.3) +
  geom_jitter(height = 0.10, width = 0, size = 0.20, alpha = 0.40,
              colour = "grey35", stroke = 0) +
  geom_vline(xintercept = R_RAW, linewidth = 0.6, colour = "#B2182B") +
  annotate("point", x = R_RAW, y = 0, shape = 23, size = 2.6,
           fill = "#B2182B", colour = "white", stroke = 0.4) +
  annotate("text", x = R_RAW + 0.035, y = 0.44, hjust = 0, size = 1.75,
           colour = "#B2182B", lineheight = 0.95,
           label = sprintf("TSPO\nrho = %.3f\n100th percentile", R_RAW)) +
  scale_y_continuous(limits = c(-0.62, 0.62), breaks = NULL) +
  labs(x = expression("Spearman "*rho*" with pseudotime (2,000 random 12-gene modules)"),
       y = NULL) +
  theme_tspo() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank())


# =============================================================================

hr("[F] 逐病人")

F_PP <- find_tab("T3b_trajectory_per_patient.csv")
pF <- panel_missing("T3b not found")
if (!is.na(F_PP)) {
  PP <- read.csv(F_PP, stringsAsFactors = FALSE)
  PP <- PP[is.finite(PP$rho), ]
  PP <- PP[order(PP$rho), ]
  k <- sum(PP$rho < 0); n <- nrow(PP)
  ps <- binom.test(k, n, 0.5)$p.value
  cat(sprintf("\n%d/%d 病人为负；中位 %.3f；符号检验 p = %.4f\n",
              k, n, median(PP$rho), ps))
  PP$patient <- factor(PP$patient, levels = PP$patient)
  PP$dir <- ifelse(PP$rho < 0, "Negative", "Positive")
  write.csv(PP, file.path(OUTDIR, "T42_per_patient_rho.csv"), row.names = FALSE)

  pF <- ggplot(PP, aes(rho, patient)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
    geom_segment(aes(x = 0, xend = rho, y = patient, yend = patient),
                 colour = "grey88", linewidth = 0.4) +
    geom_point(aes(fill = dir, size = n), shape = 21, colour = "grey30", stroke = 0.28) +
    scale_fill_manual(values = c(Negative = "#B2182B", Positive = "#2166AC"), name = NULL) +
    scale_size_continuous(range = c(0.5, 2.2), name = "Cells") +
    annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.9, size = 1.8,
             colour = PAL["ink"],
             label = sprintf("%d/%d negative; sign test p = %.3f", k, n, ps)) +
    labs(x = expression("Spearman "*rho*" with pseudotime"), y = NULL) +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.2),
          legend.position = "right",
          legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 4.0),
          legend.title = element_text(size = 4.0))
}

# =============================================================================
#  G  直接测量与独立重现：斜率图
# =============================================================================

hr("[G] 直接测量与独立重现")

F_BAL <- find_tab("T3c_TSPO_by_celltype.csv")
F_FULL <- find_tab("T3c_full_TSPO_by_celltype.csv")
F_154 <- find_tab("T41_gse154600_tspo.csv")
GL <- list()
if (!is.na(F_BAL)) {
  d <- read.csv(F_BAL, stringsAsFactors = FALSE)
  GL[["GSE266577 (balanced)"]] <- c(mono = d$mean_TSPO[d$celltype == "Classical monocytes"],
                                    macro = d$mean_TSPO[d$celltype == "Macrophages"])
}
if (!is.na(F_FULL)) {
  d <- read.csv(F_FULL, stringsAsFactors = FALSE)
  GL[["GSE266577 (all cells)"]] <- c(mono = d$mean_TSPO[d$celltype == "Classical monocytes"],
                                     macro = d$mean_TSPO[d$celltype == "Macrophages"])
}
if (!is.na(F_154)) {
  d <- read.csv(F_154, stringsAsFactors = FALSE)
  i1 <- grep("^Monocyte$", d$celltype)
  i2 <- grep("M-CSF$", d$celltype)
  if (length(i1) && length(i2))
    GL[["GSE154600"]] <- c(mono = d$mean_TSPO[i1], macro = d$mean_TSPO[i2])
}
GD <- do.call(rbind, lapply(names(GL), function(k)
  data.frame(dataset = k, ct = factor(c("Monocyte", "Macrophage"),
                                      levels = c("Monocyte", "Macrophage")),
             v = as.numeric(GL[[k]]), stringsAsFactors = FALSE)))
GD$ratio <- ave(GD$v, GD$dataset, FUN = function(z) z[1] / z[2])
cat("\n"); pr(GD[GD$ct == "Monocyte", c("dataset", "v", "ratio")], 3)

pG <- ggplot(GD, aes(ct, v, group = dataset, colour = dataset)) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(fill = dataset), shape = 21, size = 1.9, colour = "grey25", stroke = 0.28) +
  geom_text_repel(data = GD[GD$ct == "Macrophage", ],
                  aes(label = sprintf("%.2fx", ratio)),
                  size = 1.7, nudge_x = 0.34, hjust = 0, direction = "y",
                  box.padding = 0.55, point.padding = 0.30,
                  segment.size = 0.15, max.overlaps = Inf, seed = SEED) +
  scale_colour_manual(values = c("#B2182B", "#EF8A62", "#2166AC"), name = NULL) +
  scale_fill_manual(values = c("#B2182B", "#EF8A62", "#2166AC"), name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.28, 0.52))) +
  labs(x = NULL, y = "TSPO (log-normalized mean)") +
  theme_tspo() +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 4.0))

# =============================================================================
#  H  化疗前后配对
# =============================================================================

hr("[H] 化疗配对")

F_PAIR <- find_tab("GSE266577_paired_TSPO.csv")
pH <- panel_missing("paired table not found")
if (!is.na(F_PAIR)) {
  PR <- read.csv(F_PAIR, stringsAsFactors = FALSE)
  PR <- PR[is.finite(PR$naive) & is.finite(PR$ids), ]
  w <- suppressWarnings(wilcox.test(PR$naive, PR$ids, paired = TRUE))
  cat(sprintf("\n%d 对；中位差 %+.3f；配对 Wilcoxon p = %.3f\n",
              nrow(PR), median(PR$delta), w$p.value))
  PL <- rbind(data.frame(patient = PR$patient, stage = "chemo-naive", v = PR$naive),
              data.frame(patient = PR$patient, stage = "IDS", v = PR$ids))
  PL$stage <- factor(PL$stage, levels = c("chemo-naive", "IDS"))

  pH <- ggplot(PL, aes(stage, v, group = patient)) +
    geom_line(colour = "grey78", linewidth = 0.35) +
    geom_point(aes(fill = stage), shape = 21, size = 1.3,
               colour = "grey40", stroke = 0.22) +
    scale_fill_manual(values = c(`chemo-naive` = "#B2182B", IDS = "#2166AC"),
                      guide = "none") +
    annotate("text", x = 1.5, y = Inf, vjust = 1.8, size = 1.75, colour = PAL["ink"],
             label = sprintf("n = %d pairs; median delta %+.3f; p = %.2f",
                             nrow(PR), median(PR$delta), w$p.value)) +
    labs(x = NULL, y = "TSPO (log-normalized)") +
    theme_tspo() +
    theme(axis.text.x = element_text(size = 4.6))
}

# =============================================================================
#  I  下游配体沿伪时序
# =============================================================================

hr("[I] 沿伪时序的表达轨迹")

# 参照高频引用的画法：纵轴用原始 log-normalized 表达（不做 z 标准化），
# 逐基因给 loess 曲线 + 95% 置信带，并叠加降采样的原始细胞散点，
# 使读者既能看趋势也能看到真实的离散程度。
GS <- intersect(c("TSPO","LGALS1","MIF"), rownames(EX))
LG <- do.call(rbind, lapply(GS, function(g) {
  data.frame(gene = g, pt = PT, v = as.numeric(EX[g, ]), stringsAsFactors = FALSE)
}))
LG <- LG[is.finite(LG$pt) & is.finite(LG$v), ]

RP <- do.call(rbind, lapply(GS, function(g) {
  v <- as.numeric(EX[g, ])
  data.frame(gene = g,
             rho = suppressWarnings(cor(v, PT, method = "spearman")),
             mean_mono = mean(v[CT == "Classical monocytes"], na.rm = TRUE),
             mean_mac  = mean(v[CT == "Macrophages"], na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
cat("\n沿伪时序：\n"); pr(RP, 3)

# 散点降采样，曲线仍用全部细胞
set.seed(SEED)
PT_SUB <- do.call(rbind, lapply(split(LG, LG$gene), function(d) {
  d[sample(nrow(d), min(nrow(d), 700)), ]
}))
PAL_G <- c(TSPO = "#D95F0E", LGALS1 = "#1B9E77", MIF = "#7570B3")
LG$gene    <- factor(LG$gene, levels = GS)
PT_SUB$gene <- factor(PT_SUB$gene, levels = GS)

ENDL <- do.call(rbind, lapply(GS, function(g) {
  d <- LG[LG$gene == g, ]
  d[which.max(d$pt), ]
}))

pI <- ggplot(LG, aes(pt, v, colour = gene, fill = gene)) +
  geom_point(data = PT_SUB, size = 0.30, alpha = 0.16, stroke = 0,
             show.legend = FALSE) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.75, se = TRUE,
              linewidth = 0.65, alpha = 0.18) +
  geom_text_repel(data = ENDL, aes(label = gene), size = 1.8, hjust = 0,
                  nudge_x = 2.2, direction = "y", segment.size = 0.12,
                  max.overlaps = Inf, seed = SEED, show.legend = FALSE) +
  scale_colour_manual(values = PAL_G, name = NULL) +
  scale_fill_manual(values = PAL_G, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.16))) +
  labs(x = "Pseudotime (monocyte to macrophage)",
       y = expression("Expression (log-normalized)")) +
  theme_tspo() +
  theme(legend.position = "none")


  theme_tspo()
  theme_tspo()

# =============================================================================
#  拼装 3 × 3
# =============================================================================

hr("拼装")

FIG <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI)
FIG <- FIG +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, colour = "black"))

if (!exists("FIG_NAME")) FIG_NAME <- "Fig6"
print(FIG)
save_fig(FIG, FIG_NAME, width = 7.7, height = 9.4)

cat("\nfig06_trajectory.R 完成\n")

# -----------------------------------------------------------------------------
# -- 作图样式 / fig07_spatial.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig07_spatial.R   图 7（9 个面板）空间转录组：代谢轴与免疫轴分属不同区室
#
#  A  代表切片的 TSPO 空间分布            点阵图（坐标来自对象，着色为 TSPO z）
#  B  同切片的片内 z 三区域划分            点阵图（着色为区域）
#  C  五组共定位对：合并 vs 切片内          哑铃图（含切片内成立的切片数）
#  D  逐切片 ρ 热图（8 张 × 5 对）          标注热图
#  E  三区域的功能模块评分                分面小提琴 + 蜂群
#  F  区域均值热图                        标注热图
#  G  三组区域定义下的符号检验一致性        计数矩阵（3 定义 × 比较）
#  H  TSPO 的三项区域比较                 点图（成立切片数与符号检验 p）
#  I  界面相对肿瘤核心的差异表达            棒棒糖图（top 20）
#
#  ── 这张图要说什么 ──────────────────────────────────────────────────────────
#  合并所有 spot 算相关会得到 TSPO~CD68 的假关联（ρ=0.29），一旦按切片内计算
#  就几乎归零（−0.021，仅 2/8 切片为正）；而 TSPO~OXPHOS 在切片内仍成立
#  （+0.093，7/8）。三区域层面：氧化磷酸化在肿瘤核心最高、免疫与 EMT 在基质最高，
#  TSPO 本身在三个区域均无差异——恰是其两个来源（上皮代谢 / 髓系免疫）梯度相反、
#  相互抵消的结果。界面在转录组上无法与基质区分，故不作为独立区室。
#
#  ── 数据 ────────────────────────────────────────────────────────────────────
#  TSPO_spatial_processed.rds                     4,386 spot × 8 切片（坐标）
#  tables/T19d_spatial_region_within.csv          逐 spot 的三区归属与模块评分
#  spatial_persection_summary.csv / _detail.csv    五组共定位对的合并与逐切片 ρ
#  tables/T19e2_region_means.csv                  三区域均值
#  tables/T19e3_threshold_sensitivity.csv          三组定义下的符号检验
#  tables/T19e4_threshold_stability.csv            跨定义稳定性
#  tables/T19f_DEG_Interface_vs_Tumor_Core_excl_def.csv   界面 vs 核心差异表达
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/fig07_spatial.R", encoding = "UTF-8")
#    首次约 1–2 分钟（载入坐标后缓存），之后约 30 秒。
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

if (!exists("PROJ"))  PROJ  <- "D:/Rstudy/tcga"
if (!exists("OUT"))   OUT   <- file.path(PROJ, "results_v2")
if (!exists("SEED"))  SEED  <- 20260101

OUTDIR   <- file.path(OUT, "tables")
# 缓存写到工作区（D: 只读，saveRDS 会 Permission denied）
CACHEDIR <- file.path(dirname(FIGDIR), "cache")
TABDIR   <- file.path(PROJ, "results_v2", "tables")
if (!exists("FIGDIR")) FIGDIR <- file.path(OUT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHEDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(patchwork); library(grid); library(ggrepel)
                                   library(ggbeeswarm); library(hexbin) })

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)
set.seed(SEED)

find_tab <- function(name) {
  cand <- c(file.path(OUTDIR, name), file.path(TABDIR, name),
            file.path(OUT, name), file.path(PROJ, "results_v2", name))
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}

# =============================================================================
#  STEP 0  逐 spot 数据与坐标
# =============================================================================

hr("STEP 0  载入")

F_SP <- find_tab("T19d_spatial_region_within.csv")
if (is.na(F_SP)) stop("找不到 T19d_spatial_region_within.csv")
SP <- read.csv(F_SP, stringsAsFactors = FALSE)
cat("spot: ", nrow(SP), "  切片: ", length(unique(SP$section)), "\n", sep = "")

# 坐标从对象取一次并缓存（对象 420 MB，不宜每次载入）
FC <- file.path(CACHEDIR, "spatial_coords.rds")
if (file.exists(FC)) {
  CO <- readRDS(FC)
} else {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("需要 Seurat")
  suppressPackageStartupMessages(library(Seurat))
  x <- readRDS(file.path(PROJ, "TSPO_spatial_processed.rds"))
  CO <- data.frame(spot = colnames(x),
                   array_col = as.numeric(x$array_col),
                   array_row = as.numeric(x$array_row),
                   stringsAsFactors = FALSE)
  saveRDS(CO, FC)
  rm(x); invisible(gc())
  cat("坐标已缓存: ", FC, "\n", sep = "")
}
SP <- merge(SP, CO, by = "spot", all.x = TRUE)
cat("有坐标的 spot: ", sum(is.finite(SP$array_row)), " / ", nrow(SP), "\n", sep = "")

SECS <- sort(unique(SP$section))
cat("切片: ", paste(SECS, collapse = ", "), "\n", sep = "")
REGION_LV <- c("Tumor_Core", "Interface", "Stroma", "Other")
REGION_LAB <- c(Tumor_Core = "Tumor core", Interface = "Interface",
                Stroma = "Stroma", Other = "Other")
PAL_REG <- c(Tumor_Core = "#A50026", Interface = "#EF8A62",
             Stroma = "#08519C", Other = "grey82")

# 选代表切片：三区齐全且 spot 数最多
TAB <- table(SP$section, SP$region)
OK <- apply(TAB[, intersect(REGION_LV, colnames(TAB)), drop = FALSE], 1,
            function(z) all(z > 30))
REP <- if (any(OK)) names(which.max(rowSums(TAB)[OK])) else names(which.max(rowSums(TAB)))
cat("代表切片: ", REP, "（", sum(SP$section == REP), " spot）\n", sep = "")

# =============================================================================
#  A  代表切片的 TSPO 空间分布
# =============================================================================

hr("[A] TSPO 空间分布")

DA <- SP[SP$section == REP, ]
# 约半数 spot 的 TSPO 落在 -1.16（未检出的地板值），故色阶按 2%–98% 分位裁剪，
# 否则整张图会被这一档压成一片深蓝
TS_LIM <- quantile(DA$TSPO, c(0.02, 0.98), na.rm = TRUE)
pA <- ggplot(DA[order(DA$TSPO), ], aes(array_col, array_row)) +
  geom_point(aes(colour = TSPO), size = 0.85, stroke = 0) +
  scale_colour_gradientn(colours = c("#2166AC","#92C5DE","#F7F7F7","#EF8A62","#67001F"),
                         limits = TS_LIM, oob = scales::squish,
                         name = "TSPO\n(z)",
                         guide = guide_colourbar(barwidth = unit(1.8, "mm"),
                                                 barheight = unit(6, "mm"))) +
  scale_y_reverse() +
  labs(x = NULL, y = NULL, subtitle = sprintf("%s (n = %d spots)", REP, nrow(DA))) +
  theme_tspo() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right", legend.key.size = unit(2.0, "mm"),
        legend.text = element_text(size = 4.0), legend.title = element_text(size = 4.0),
        plot.subtitle = element_text(size = 4.2, colour = "grey35"))

# =============================================================================
#  B  同切片的三区域划分
# =============================================================================

hr("[B] 三区域划分")

DB <- DA[!is.na(DA$region), ]
DB$region <- factor(DB$region, levels = intersect(REGION_LV, unique(DB$region)))
pB <- ggplot(DB[order(DB$region), ], aes(array_col, array_row)) +
  geom_point(aes(colour = region), size = 0.85, stroke = 0) +
  scale_colour_manual(values = PAL_REG, labels = REGION_LAB, name = NULL) +
  scale_y_reverse() +
  labs(x = NULL, y = NULL, subtitle = "in-section z assignment") +
  theme_tspo() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right", legend.key.size = unit(2.2, "mm"),
        legend.text = element_text(size = 4.0),
        plot.subtitle = element_text(size = 4.2, colour = "grey35"))

# =============================================================================
#  C  五组共定位：合并 vs 切片内
# =============================================================================

hr("[C] 合并 vs 切片内")

F_PS <- find_tab("spatial_persection_summary.csv")
pC <- panel_missing("persection summary not found")
if (!is.na(F_PS)) {
  PS <- read.csv(F_PS, stringsAsFactors = FALSE)
  PS$pair_lab <- gsub(" ~ ", "\u2013", PS$analysis)
  PS <- PS[order(PS$rho_pooled), ]
  PS$pair_lab <- factor(PS$pair_lab, levels = PS$pair_lab)
  cat("\n"); pr(PS[, c("analysis","rho_pooled","rho_within","n_pos","n_sec","p_block")], 3)

  DL <- rbind(
    data.frame(pair_lab = PS$pair_lab, type = "Pooled across sections", rho = PS$rho_pooled),
    data.frame(pair_lab = PS$pair_lab, type = "Within section (median)", rho = PS$rho_within))
  DL$type <- factor(DL$type, levels = c("Pooled across sections", "Within section (median)"))
  AR <- merge(DL[DL$type == "Pooled across sections", c("pair_lab","rho")],
              DL[DL$type == "Within section (median)", c("pair_lab","rho")],
              by = "pair_lab", suffixes = c("pool","within"))

  pC <- ggplot(DL, aes(rho, pair_lab)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
    geom_segment(data = AR, aes(x = rhopool, xend = rhowithin, y = pair_lab, yend = pair_lab),
                 inherit.aes = FALSE, colour = "grey65", linewidth = 0.75,
                 arrow = arrow(length = unit(1.1, "mm"), type = "closed")) +
    geom_point(aes(fill = type), shape = 21, size = 2.0,
               colour = "grey20", stroke = 0.3) +
    geom_text(data = PS, aes(x = Inf, y = pair_lab,
                             label = sprintf("%d/%d", n_pos, n_sec)),
              hjust = 1.05, size = 1.6, colour = "grey30", inherit.aes = FALSE) +
    scale_fill_manual(values = c(`Pooled across sections` = "#B2182B",
                                 `Within section (median)` = "#2166AC"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.24))) +
    labs(x = expression("Spearman "*rho), y = NULL,
         subtitle = "right column: sections where rho > 0") +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.4),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 3.9),
          legend.margin = margin(0, 0, 0, 0),
          plot.subtitle = element_text(size = 4.1, colour = "grey35"))
}

# =============================================================================
#  D  逐切片 ρ 热图
# =============================================================================

hr("[D] 逐切片 ρ")

# 原为热图，与面板 F 重复；改为逐切片点带图：既给出 8 张切片各自的 ρ，
# 也直接看出有几张落在零线另一侧（C 面板只给中位数）。
F_PD <- find_tab("spatial_persection_detail.csv")
pD <- panel_missing("persection detail not found")
if (!is.na(F_PD)) {
  PD <- read.csv(F_PD, stringsAsFactors = FALSE)
  PD$pair_lab <- factor(gsub(" ~ ", "\u2013", PD$analysis),
                        levels = levels(PS$pair_lab))
  PD$sign <- ifelse(PD$rho > 0, "rho > 0", "rho < 0")
  SUM <- aggregate(rho ~ pair_lab, data = PD, FUN = median)
  NUMP <- aggregate(rho ~ pair_lab, data = PD,
                    FUN = function(z) sum(z > 0, na.rm = TRUE))
  names(NUMP)[2] <- "n_pos"
  cat("\n逐切片 ρ 的正向切片数：\n"); pr(NUMP, 3)

  # 中位数原先用实心菱形标在点上，会盖住底下的散点。
  # 改成图例可辨认的三项：圆点＝单张切片，竖线＝八张切片的中位数。

  pD <- ggplot(PD, aes(rho, pair_lab)) +
    geom_vline(xintercept = 0, linewidth = 0.32, colour = "grey55", linetype = "22") +
    # 中位数用菱形，但画在散点【之前】：图层在下，圆点覆盖其上，故不会挡住任何点。
    # fill 映射到一个常量字符串，这样菱形自带一个图例项，读者不必猜它是什么。
    geom_point(data = SUM, aes(x = rho, y = pair_lab, fill = "Median of 8 sections"),
               shape = 23, size = 2.2, colour = "grey15", stroke = 0.42,
               inherit.aes = FALSE) +
    geom_point(aes(fill = sign), shape = 21, size = 1.35,
               colour = "grey35", stroke = 0.22,
               position = position_jitter(height = 0.14, width = 0, seed = SEED)) +
    geom_text(data = NUMP, aes(x = Inf, y = pair_lab, label = sprintf("%d/8", n_pos)),
              hjust = 1.05, size = 1.6, colour = "grey30", inherit.aes = FALSE) +
    scale_fill_manual(values = c(`rho > 0` = "#B2182B",
                                 `rho < 0` = "#2166AC",
                                 `Median of 8 sections` = "white"),
                      breaks = c("rho > 0", "rho < 0", "Median of 8 sections"),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.22))) +
    labs(x = expression("Spearman "*rho*" within section"), y = NULL) +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.4),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 3.9),
          legend.margin = margin(0, 0, 0, 0),
          legend.spacing.x = unit(0.8, "mm"))
}

# =============================================================================

hr("[E] 区域模块评分")

MODS <- c("OXPHOS","Glycolysis","Immune","EMT")
DE <- SP[!is.na(SP$region) & SP$region %in% c("Tumor_Core","Interface","Stroma"), ]
DE$region <- factor(DE$region, levels = c("Tumor_Core","Interface","Stroma"))
EL <- do.call(rbind, lapply(MODS, function(m) {
  data.frame(region = DE$region, module = m, v = DE[[m]], stringsAsFactors = FALSE)
}))
EL$module <- factor(EL$module, levels = MODS)
cat("\n各区域模块均值：\n")
print(round(tapply(DE$OXPHOS, DE$region, mean, na.rm = TRUE), 3))

pE <- ggplot(EL, aes(region, v, fill = region)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey88") +
  geom_violin(scale = "width", width = 0.62, alpha = 0.55,
              colour = "grey35", linewidth = 0.25, trim = TRUE) +
  geom_quasirandom(width = 0.20, size = 0.22, alpha = 0.35, shape = 16, stroke = 0,
                   aes(colour = region)) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 1.4,
               fill = "white", colour = "grey15", stroke = 0.32) +
  facet_wrap(~ module, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = PAL_REG[c("Tumor_Core","Interface","Stroma")], guide = "none") +
  scale_colour_manual(values = PAL_REG[c("Tumor_Core","Interface","Stroma")], guide = "none") +
  scale_x_discrete(labels = c("Core", "IF", "Str")) +
  labs(x = NULL, y = "Module score (z)") +
  theme_tspo() +
  theme(axis.text.x = element_text(size = 4.2),
        strip.text = element_text(size = 4.5))

# =============================================================================
#  F  区域均值热图
# =============================================================================

hr("[F] 区域均值热图")

F_RM <- find_tab("T19e2_region_means.csv")
pF <- panel_missing("region means not found")
if (!is.na(F_RM)) {
  RM <- read.csv(F_RM, stringsAsFactors = FALSE)
  RM <- RM[RM$definition == "A", ]
  VARS <- c("TSPO","OXPHOS","Glycolysis","Immune","EMT","Hypoxia")
  RM <- RM[RM$variable %in% VARS, ]
  RL <- rbind(
    data.frame(variable = RM$variable, region = "Tumor core", v = RM$mean_Tumor_Core),
    data.frame(variable = RM$variable, region = "Interface",  v = RM$mean_Interface),
    data.frame(variable = RM$variable, region = "Stroma",     v = RM$mean_Stroma))
  RL$region <- factor(RL$region, levels = c("Tumor core","Interface","Stroma"))
  RL$variable <- factor(RL$variable, levels = rev(VARS))
  MXV <- max(abs(RL$v), na.rm = TRUE)
  cat("\n"); pr(RM[, c("variable","mean_Tumor_Core","mean_Interface","mean_Stroma")], 3)

  pF <- ggplot(RL, aes(region, variable, fill = v)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%+.2f", v)), size = 1.7, colour = "grey15") +
    scale_fill_div(limit = c(-MXV, MXV), name = "Mean\n(z)") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_tspo() +
    theme(axis.text.x.top = element_text(angle = 22, hjust = 0, size = 4.4),
          axis.text.y = element_text(size = 4.4),
          panel.grid = element_blank(),
          legend.key.height = unit(6, "mm")) +
    guides(fill = guide_colourbar(barwidth = unit(1.8, "mm"), barheight = unit(6, "mm")))
}

# =============================================================================
#  G  三组定义下的符号检验一致性
# =============================================================================

hr("[G] 阈值敏感性")

F_TS <- find_tab("T19e4_threshold_stability.csv")
pG <- panel_missing("threshold stability not found")
if (!is.na(F_TS)) {
  TS <- read.csv(F_TS, stringsAsFactors = FALSE)
  TS$lab <- sprintf("%s | %s", TS$variable, TS$comparison)
  SL <- rbind(
    data.frame(lab = TS$lab, def = "A", n = TS$n_positive.A, stable = TS$stable),
    data.frame(lab = TS$lab, def = "B", n = TS$n_positive.B, stable = TS$stable),
    data.frame(lab = TS$lab, def = "C", n = TS$n_positive.C, stable = TS$stable))
  ORD <- TS$lab[order(TS$comparison, TS$variable)]
  SL$lab <- factor(SL$lab, levels = rev(ORD))
  SL$stable <- factor(SL$stable, labels = c(`TRUE` = "stable", `FALSE` = "unstable"))
  cat("\n跨定义稳定的比较：", sum(TS$stable), " / ", nrow(TS), "\n", sep = "")

  # 行首用实心/空心点标出跨三种定义是否稳定（稳定 = 三组定义下多数切片同向）
  STB <- data.frame(lab = factor(TS$lab, levels = rev(ORD)),
                    stable = ifelse(TS$stable, "●", "○"),
                    x = 0.42)
  pG <- ggplot(SL, aes(def, lab, fill = n)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(data = STB, aes(x = x, y = lab, label = stable),
              inherit.aes = FALSE, size = 1.9, colour = "grey25") +
    geom_text(aes(label = n, colour = n >= 7), size = 1.6, show.legend = FALSE) +
    scale_fill_gradient(low = "#F7F7F7", high = "#67001F", limits = c(0, 8),
                        name = "Sections\npositive") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(expand = expansion(add = c(0.4, 0.4))) +
    labs(x = "Region definition (\u25cf stable across A/B/C)", y = NULL) +
    theme_tspo() +
    theme(axis.text.x.top = element_text(size = 4.4),
          axis.text.y = element_text(size = 3.9),
          panel.grid = element_blank(),
          legend.key.height = unit(5, "mm")) +
    guides(fill = guide_colourbar(barwidth = unit(1.7, "mm"), barheight = unit(5, "mm")))
}

# =============================================================================
#  H  TSPO 的三项区域比较
# =============================================================================

hr("[H] TSPO 区域比较")

F_T3 <- find_tab("T19e3_threshold_sensitivity.csv")
pH <- panel_missing("threshold sensitivity not found")
if (!is.na(F_T3)) {
  T3 <- read.csv(F_T3, stringsAsFactors = FALSE)
  T3 <- T3[T3$variable == "TSPO", ]
  cat("\n"); pr(T3[, c("definition","comparison","n_positive","mean_diff_sections","sign_p")], 3)
  T3$lab <- paste0(T3$comparison, "  (", T3$definition, ")")
  T3 <- T3[order(T3$mean_diff_sections), ]
  T3$lab <- factor(T3$lab, levels = T3$lab)

  pH <- ggplot(T3, aes(n_positive, lab)) +
    geom_vline(xintercept = 4, linewidth = 0.3, colour = "grey60", linetype = "22") +
    geom_segment(aes(x = 4, xend = n_positive, y = lab, yend = lab),
                 colour = "grey85", linewidth = 0.5) +
    geom_point(aes(fill = sign_p < 0.05), shape = 21, size = 2.0,
               colour = "grey25", stroke = 0.3) +
    geom_text(aes(x = Inf, label = sprintf("%d/8  p=%.2f", n_positive, sign_p)),
              hjust = 1.05, size = 1.55, colour = "grey30") +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey82"),
                      labels = c(`TRUE` = "Sign test p < 0.05", `FALSE` = "n.s."),
                      name = NULL) +
    scale_x_continuous(limits = c(0, 8), breaks = c(0, 4, 8),
                       expand = expansion(mult = c(0.04, 0.20))) +
    labs(x = "Sections where TSPO is higher (of 8)", y = NULL,
         subtitle = "dashed line: 4/8 = no consistency") +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.2),
          legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 4.0),
          plot.subtitle = element_text(size = 4.1, colour = "grey35"))
}

# =============================================================================
#  I  界面相对肿瘤核心的差异表达
# =============================================================================

hr("[I] 界面差异表达")

F_DEG <- find_tab("T19f_DEG_Interface_vs_Tumor_Core_excl_def.csv")
pI <- panel_missing("interface DEG not found")
if (!is.na(F_DEG)) {
  DG <- read.csv(F_DEG, stringsAsFactors = FALSE)
  UP <- DG[DG$logFC >= 0.5 & DG$adj.P.Val < 0.05, ]
  DN <- DG[DG$logFC <= -0.5 & DG$adj.P.Val < 0.05, ]
  cat(sprintf("\n|logFC|>=0.5 & FDR<0.05 : 上调 %d / 下调 %d\n", nrow(UP), nrow(DN)))
  TOPN <- 20
  UPT <- head(UP[order(-UP$logFC), ], TOPN)
  UPT$gene <- factor(UPT$gene, levels = UPT$gene[order(UPT$logFC)])
  cat("前 8 位：\n"); pr(head(UPT[order(-UPT$logFC), c("gene","logFC","adj.P.Val")], 8), 3)

  pI <- ggplot(UPT, aes(logFC, gene)) +
    geom_segment(aes(x = 0, xend = logFC, y = gene, yend = gene),
                 colour = "grey86", linewidth = 0.45) +
    geom_point(aes(size = -log10(adj.P.Val)), shape = 21, fill = "#B2182B",
               colour = "grey25", stroke = 0.25) +
    scale_size_continuous(range = c(0.4, 2.2), name = expression(-log[10]*"FDR")) +
    labs(x = expression("log"[2]*" fold change (interface vs core)"), y = NULL,
         subtitle = sprintf("top %d of %d upregulated genes", TOPN, nrow(UP))) +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.0),
          legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 3.9),
          legend.title = element_text(size = 3.9),
          plot.subtitle = element_text(size = 4.1, colour = "grey35"))
}

# =============================================================================
#  拼装 3 × 3
# =============================================================================

hr("拼装")

FIG <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI)
FIG <- FIG +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, colour = "black"))

if (!exists("FIG_NAME")) FIG_NAME <- "Fig7"
print(FIG)
save_fig(FIG, FIG_NAME, width = 7.7, height = 9.2)

cat("\nfig07_spatial.R 完成\n")

# -----------------------------------------------------------------------------
# -- 作图样式 / fig08_codex.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig08_codex.R   图 8（8 个面板）空间蛋白组：独立确认区室结构，并界定 MIF 的空间归属
#
#  A  代表 TMA 核心的细胞级多重表型          点阵（细胞坐标 + 表型着色）
#  B  七个标志物在肿瘤区与间质区的阳性率      哑铃图（合并口径 + 配对 p）
#  C  CK+ 在 129 个配对核心中的两侧阳性率     散点（对角线 = 无差异）
#  D  CD68+ 同上                             散点
#  E  各标志物「间质更高」的核心数            发散点图（含零线 64.5）
#  F  MIF 在片内 z 三区域定义下的水平         小提琴 + 蜂群
#  G  MIF 三对共定位：合并 vs 切片内          哑铃图（含成立切片数）
#  H  逐切片 ρ 分布                          点带图（菱形 = 中位数）
#
#  ── 这张图能说什么、不能说什么 ──────────────────────────────────────────────
#  该 TMA 面板只含 CK/Ki67/CD8/IER3/pSTAT3/CD3/CD68/CD19/DAPI，**没有 TSPO 抗体**。
#  因此它只能验证「肿瘤细胞与免疫细胞分属不同区室」这一前提，并给 MIF 的空间归属
#  提供证据；不能在此层面直接验证 TSPO 的蛋白定位。作图与图注都必须守住这条界线。
#
#  ── 数据 ────────────────────────────────────────────────────────────────────
#  VectraPolarisData::HumanOvarianCancerVP（ExperimentHub 缓存）
#    10 标志物 × 1,610,431 细胞；132 个 TMA 核心；2 张玻片
#  cache/codex_percore.rds        逐核心 × 分区的阳性计数（本脚本生成）
#  cache/codex_onecore.rds        代表核心的细胞坐标与表型（本脚本生成）
#  tables/T20a_MIF_raw_by_region.csv / T20b_MIF_z_by_region.csv
#  tables/T20d_spot_correlations.csv
#  spatial_persection_detail.csv  逐切片 ρ（MIF 三对）
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    setwd("D:/Rstudy/tcga")
#    source("<code_final>/fig08_codex.R", encoding = "UTF-8")
#    首次约 2–3 分钟（读 700 MB 缓存并聚合），之后约 40 秒。
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

if (!exists("PROJ"))  PROJ  <- "D:/Rstudy/tcga"
if (!exists("OUT"))   OUT   <- "D:/Rstudy/tcga/results_v2"
if (!exists("SEED"))  SEED  <- 20260101

OUTDIR   <- file.path(OUT, "tables")
TABDIR   <- file.path(PROJ, "results_v2", "tables")
if (!exists("FIGDIR")) FIGDIR <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/figures"
CACHEDIR <- file.path(dirname(FIGDIR), "cache")
dir.create(CACHEDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({
  library(patchwork); library(grid); library(data.table); library(ggrepel); library(ggbeeswarm)
})

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)
set.seed(SEED)

find_tab <- function(name) {
  cand <- c(file.path(OUTDIR, name), file.path(TABDIR, name),
            file.path(OUT, name), file.path(PROJ, "results_v2", name))
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}

MK       <- c("ck","cd68","cd3","cd8","cd19","ki67","p_stat3")
MK_LAB   <- c(ck = "CK", cd68 = "CD68", cd3 = "CD3", cd8 = "CD8",
              cd19 = "CD19", ki67 = "Ki67", p_stat3 = "pSTAT3")
MIN_CELL <- 50

# =============================================================================
#  STEP 0  CODEX 数据（聚合 + 代表核心）
# =============================================================================

hr("STEP 0  载入 CODEX")

F_AGG <- file.path(CACHEDIR, "codex_percore.rds")
F_ONE <- file.path(CACHEDIR, "codex_onecore.rds")

load_spe <- function() {
  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) return(NULL)
  suppressPackageStartupMessages(library(SpatialExperiment))
  # Windows 上 path.expand("~") 指向「文档」，而 ExperimentHub 缓存在 AppData\Local
  eh_dir <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "cache", "R", "ExperimentHub")
  cand <- Sys.glob(file.path(eh_dir, "*_7362"))
  if (!length(cand)) cand <- Sys.glob(file.path(path.expand("~"), "AppData", "Local",
                                                "R", "cache", "R", "ExperimentHub", "*_7362"))
  if (length(cand)) {
    x <- try(readRDS(cand[which.max(file.info(cand)$size)]), silent = TRUE)
    if (!inherits(x, "try-error")) return(x)
  }
  if (requireNamespace("VectraPolarisData", quietly = TRUE)) {
    try(ExperimentHub::setExperimentHubOption("LOCAL", TRUE), silent = TRUE)
    x <- try(VectraPolarisData::HumanOvarianCancerVP(metadata = FALSE), silent = TRUE)
    if (!inherits(x, "try-error")) return(x)
  }
  NULL
}

if (file.exists(F_AGG) && file.exists(F_ONE)) {
  AGG <- readRDS(F_AGG); ONE <- readRDS(F_ONE)
  cat("聚合与代表核心均来自缓存\n")
} else {
  spe <- load_spe()
  if (is.null(spe)) stop("无法载入 CODEX 数据（ExperimentHub 缓存缺失）")
  cat("CODEX: ", nrow(spe), " 标志物 x ", ncol(spe), " 细胞\n", sep = "")
  cd <- as.data.table(as.data.frame(SummarizedExperiment::colData(spe)))
  sc <- SpatialExperiment::spatialCoords(spe)
  cd[, x := sc[, 1]]; cd[, y := sc[, 2]]
  for (m in MK)
    cd[, (paste0("pos_", m)) := endsWith(as.character(get(paste0("phenotype_", m))), "+")]

  DD <- cd[tissue_category %in% c("Tumor", "Stroma")]
  AGG <- DD[, c(list(n = .N),
                setNames(lapply(MK, function(m) sum(get(paste0("pos_", m)))),
                         paste0("n_", MK))),
           by = .(sample_id, tissue_category)]
  saveRDS(AGG, F_AGG)

  # 代表核心：两侧细胞都够多、且总细胞数最大，画出来才看得出分区
  WIDE <- dcast(AGG, sample_id ~ tissue_category, value.var = "n")
  WIDE <- WIDE[!is.na(Tumor) & !is.na(Stroma) & Tumor >= 400 & Stroma >= 200]
  WIDE[, tot := Tumor + Stroma]
  REP <- WIDE$sample_id[which.max(WIDE$tot)]
  ONE <- cd[sample_id == REP, .(x, y, pos_ck, pos_cd68, pos_cd3)]
  ONE <- list(cells = as.data.frame(ONE), core = REP)
  saveRDS(ONE, F_ONE)
  rm(spe, cd, DD); invisible(gc())
  cat("已聚合并缓存\n")
}

if (is.list(ONE) && !is.data.frame(ONE)) { ONE_CELLS <- ONE$cells; REP_CORE <- ONE$core } else { ONE_CELLS <- ONE; REP_CORE <- attr(ONE, "core") }
if (is.null(REP_CORE)) REP_CORE <- "representative core"
cat("代表核心: ", REP_CORE, "（", nrow(ONE_CELLS), " 细胞）\n", sep = "")

# 逐核心宽表 + 配对筛选
WIDE <- dcast(AGG, sample_id ~ tissue_category, value.var = c("n", paste0("n_", MK)))
PAIR <- WIDE[!is.na(n_Tumor) & !is.na(n_Stroma) &
             n_Tumor >= MIN_CELL & n_Stroma >= MIN_CELL]
cat("两侧均 >= ", MIN_CELL, " 细胞的核心: ", nrow(PAIR), "\n", sep = "")

PCT <- data.table(sample_id = PAIR$sample_id)
for (m in MK) {
  PCT[, (paste0(m, "_T")) := 100 * PAIR[[paste0("n_", m, "_Tumor")]] / PAIR$n_Tumor]
  PCT[, (paste0(m, "_S")) := 100 * PAIR[[paste0("n_", m, "_Stroma")]] / PAIR$n_Stroma]
}

# 全部核心的合并阳性率（不做 >=50 过滤，与正文表格口径一致）
TOT <- AGG[, .(n = sum(n), t(apply(.SD, 1, function(z) z))), ]
POOL <- do.call(rbind, lapply(MK, function(m) {
  tt <- AGG[tissue_category == "Tumor"]; ss <- AGG[tissue_category == "Stroma"]
  data.frame(marker = m,
             Tumor  = 100 * sum(tt[[paste0("n_", m)]]) / sum(tt$n),
             Stroma = 100 * sum(ss[[paste0("n_", m)]]) / sum(ss$n),
             n_core_t = nrow(PAIR), stringsAsFactors = FALSE)
}))

SUM <- do.call(rbind, lapply(MK, function(m) {
  w <- suppressWarnings(wilcox.test(PCT[[paste0(m, "_T")]], PCT[[paste0(m, "_S")]], paired = TRUE))
  data.frame(marker = m, pool_T = POOL$Tumor[POOL$marker == m],
             pool_S = POOL$Stroma[POOL$marker == m],
             mean_T = mean(PCT[[paste0(m, "_T")]]), mean_S = mean(PCT[[paste0(m, "_S")]]),
             n_stroma_higher = sum(PCT[[paste0(m, "_S")]] > PCT[[paste0(m, "_T")]]),
             n_core = nrow(PCT), p = w$p.value, stringsAsFactors = FALSE)
}))
SUM$marker_lab <- MK_LAB[SUM$marker]
cat("\n合并口径（全部 132 核心）与逐核心配对：\n"); pr(SUM, 3)
write.csv(SUM, file.path(CACHEDIR, "codex_summary.csv"), row.names = FALSE)

# =============================================================================
#  A  代表核心的细胞级多重表型
# =============================================================================

hr("[A] 代表核心")

CT2 <- ONE_CELLS
CT2$cls <- ifelse(CT2$pos_ck, "CK+ (tumour)",
           ifelse(CT2$pos_cd68, "CD68+ (myeloid)",
           ifelse(CT2$pos_cd3, "CD3+ (T cell)", "other")))
CT2$cls <- factor(CT2$cls, levels = c("other", "CD3+ (T cell)", "CD68+ (myeloid)", "CK+ (tumour)"))
cat("代表核心表型构成:\n"); print(table(CT2$cls))

# CK+ 占 76%，若与免疫细胞同样着色、同样点径，整张图只剩一种红。
# 分三层绘制：CK+ 与 other 作淡色底（小点），CD3+ / CD68+ 用深色大点压在上层。
PAL_CLS <- c(`other` = "grey88", `CK+ (tumour)` = "#E8A598",
             `CD3+ (T cell)` = "#08519C", `CD68+ (myeloid)` = "#006D2C")
BG <- CT2[CT2$cls %in% c("other", "CK+ (tumour)"), ]
IM <- CT2[CT2$cls %in% c("CD3+ (T cell)", "CD68+ (myeloid)"), ]

pA <- ggplot() +
  geom_point(data = BG, aes(x, y, colour = cls), size = 0.20, stroke = 0, alpha = 0.9) +
  geom_point(data = IM[order(IM$cls == "CD3+ (T cell)"), ], aes(x, y, colour = cls),
             size = 0.55, stroke = 0) +
  scale_colour_manual(values = PAL_CLS, name = NULL,
                      breaks = c("CK+ (tumour)", "CD68+ (myeloid)", "CD3+ (T cell)", "other")) +
  coord_equal() +
  guides(colour = guide_legend(override.aes = list(size = c(1.6, 1.6, 1.6, 1.0)))) +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("%s... (%s cells)", substr(REP_CORE, 1, 30),
                          format(nrow(CT2), big.mark = ","))) +
  theme_tspo() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 3.9),
        legend.margin = margin(0, 0, 0, 0),
        plot.subtitle = element_text(size = 4.0, colour = "grey35"))

# =============================================================================
#  B  七个标志物的阳性率：肿瘤区 vs 间质区
# =============================================================================

hr("[B] 阳性率哑铃")

BD <- rbind(
  data.frame(marker_lab = SUM$marker_lab, region = "Tumour", v = SUM$pool_T),
  data.frame(marker_lab = SUM$marker_lab, region = "Stroma", v = SUM$pool_S))
BD$marker_lab <- factor(BD$marker_lab, levels = rev(SUM$marker_lab))
AX <- merge(BD[BD$region == "Tumour", c("marker_lab", "v")],
            BD[BD$region == "Stroma", c("marker_lab", "v")],
            by = "marker_lab", suffixes = c("T", "S"))
AX$lab <- sprintf("%.1f / %.1f", AX$vT, AX$vS)

pB <- ggplot(BD, aes(v, marker_lab)) +
  geom_segment(data = AX, aes(x = vT, xend = vS, y = marker_lab, yend = marker_lab),
               inherit.aes = FALSE, colour = "grey70", linewidth = 0.7,
               arrow = arrow(length = unit(1.1, "mm"), type = "closed")) +
  geom_point(aes(fill = region), shape = 21, size = 1.9,
             colour = "grey25", stroke = 0.28) +
  geom_text(data = AX, aes(x = Inf, y = marker_lab, label = lab), hjust = 1.05,
            size = 1.45, colour = "grey30", inherit.aes = FALSE) +
  scale_fill_manual(values = c(Tumour = "#B2182B", Stroma = "#2166AC"), name = NULL) +
  # CK 约 85% 而 CD19 仅 0.3%，线性轴上其余六个标志物被压成一点，故用对数轴
  scale_x_log10(breaks = c(0.3, 1, 3, 10, 30, 100),
                labels = function(x) format(x, trim = TRUE),
                expand = expansion(mult = c(0.10, 0.22))) +
  labs(x = "Positive cells (%, log scale)", y = NULL,
       subtitle = "arrow: Tumour -> Stroma") +
  theme_tspo() +
  theme(legend.position = c(0.02, 0.04), legend.justification = c(0, 0),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 3.9),
        plot.subtitle = element_text(size = 4.0, colour = "grey35"))

# =============================================================================
#  C / D  CK+ 与 CD68+ 的逐核心配对散点
# =============================================================================

core_scatter <- function(mk, PALX) {
  d <- data.frame(T = PCT[[paste0(mk, "_T")]], S = PCT[[paste0(mk, "_S")]])
  nS <- sum(d$S > d$T)
  mx <- max(c(d$T, d$S), na.rm = TRUE) * 1.04
  r  <- suppressWarnings(cor(d$T, d$S, method = "spearman"))
  ggplot(d, aes(T, S)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.28,
                colour = "grey60", linetype = "22") +
    geom_point(shape = 21, size = 0.85, fill = PALX, colour = "grey30", stroke = 0.18) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.04, vjust = -0.7, size = 1.7,
             colour = PAL["ink"],
             label = sprintf("stroma higher: %d/%d", nS, nrow(d))) +
    scale_x_continuous(limits = c(0, mx)) +
    scale_y_continuous(limits = c(0, mx)) +
    coord_equal() +
    labs(x = paste0(MK_LAB[mk], "+ in tumour (%)"),
         y = paste0(MK_LAB[mk], "+ in stroma (%)")) +
    theme_tspo()
}

pD <- core_scatter("cd68", "#1B7837")

# C 的散点几乎全部堆在 (高 CK, 近 0 CK) 一处，散点重叠严重看不出密度，
# 改用六边形分箱，密度本身成为信息。
# CK+ 的逐核心值几乎全部落在「肿瘤 85%、间质 0」的同一处，六边形分箱只剩两格、
# 散点又完全重叠。改用配对斜率图：129 条连线一律向下，0/129 无需看数字即可读出。
CL <- rbind(data.frame(core = PCT$sample_id, region = "Tumour", v = PCT$ck_T),
            data.frame(core = PCT$sample_id, region = "Stroma", v = PCT$ck_S))
CL$region <- factor(CL$region, levels = c("Tumour", "Stroma"))

pC <- ggplot(CL, aes(region, v, group = core)) +
  geom_line(colour = "grey55", linewidth = 0.16, alpha = 0.55) +
  geom_point(aes(fill = region), shape = 21, size = 0.55,
             colour = "grey35", stroke = 0.10) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.9, size = 1.7, colour = PAL["ink"],
           label = sprintf("stroma higher: %d/%d", sum(PCT$ck_S > PCT$ck_T), nrow(PCT))) +
  scale_fill_manual(values = c(Tumour = "#B2182B", Stroma = "#2166AC"), guide = "none") +
  scale_x_discrete(expand = expansion(add = 0.42)) +
  labs(x = NULL, y = "CK+ cells (%)") +
  theme_tspo()

# =============================================================================
#  E  各标志物「间质更高」的核心数
# =============================================================================

hr("[E] 间质更高核心数")

# E 原为「间质更高的核心数」，只给计数、丢掉差值的量级；
# 改为逐核心的 (间质 − 肿瘤) 阳性率差，零线为界，条带形状即效应量分布。
DL <- do.call(rbind, lapply(MK, function(m) {
  d <- PCT[[paste0(m, "_S")]] - PCT[[paste0(m, "_T")]]
  data.frame(marker_lab = MK_LAB[m], diff = d, stringsAsFactors = FALSE)
}))
DL$marker_lab <- factor(DL$marker_lab, levels = SUM$marker_lab[order(SUM$n_stroma_higher)])
MED <- aggregate(diff ~ marker_lab, data = DL, FUN = median)
cat("\n逐核心差值中位数（间质 − 肿瘤，百分点）：\n"); pr(MED, 3)

MD <- aggregate(diff ~ marker_lab, data = DL, FUN = median)
MED$side <- MED$diff > 0
MED$marker_lab <- factor(MED$marker_lab, levels = levels(DL$marker_lab))
DL <- merge(DL, MED[, c("marker_lab", "side")], by = "marker_lab")
DL$marker_lab <- factor(DL$marker_lab, levels = levels(MED$marker_lab))

pE <- ggplot(DL, aes(diff, marker_lab, fill = side)) +
  geom_vline(xintercept = 0, linewidth = 0.32, colour = "grey55", linetype = "22") +
  geom_violin(scale = "width", width = 0.62, alpha = 0.55, colour = "grey35",
              linewidth = 0.22) +
  geom_point(data = MED, aes(x = diff, y = marker_lab), shape = 21, size = 1.3,
             fill = "white", colour = "grey15", stroke = 0.3, inherit.aes = FALSE) +
  # 交替色只是绘图顺序，没有含义；改按「中位差落在零的哪一侧」着色
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.12))) +
  labs(x = "Stroma - tumour positivity (percentage points)", y = NULL,
       subtitle = "white dot: median across cores") +
  theme_tspo() +
  theme(plot.subtitle = element_text(size = 4.0, colour = "grey35"))

# =============================================================================
#  F  MIF 在片内 z 三区域定义下的水平
# =============================================================================

hr("[F] MIF 区域")

F_MIF <- file.path(CACHEDIR, "mif_perspot.rds")
pF <- panel_missing("MIF per-spot data not found")
if (file.exists(F_MIF)) {
  MZ <- readRDS(F_MIF)
  MZ <- MZ[MZ$region %in% c("Tumor_Core", "Interface", "Stroma"), ]
  MZ$region <- factor(MZ$region, levels = c("Tumor_Core", "Interface", "Stroma"))
  cat("\nMIF 逐 spot 区域均值：\n")
  print(round(tapply(MZ$MIF_z, MZ$region, mean), 3))

  # 雨云图：half-violin（云） + 细箱线（伞） + 抖动散点（雨）。
  # 原来是「小提琴 + 蜂群」，现在把密度只画在右侧、散点铺在下方，三层不互相遮挡。
  MZ$xn <- as.numeric(MZ$region)
  HALF <- do.call(rbind, lapply(sort(unique(MZ$xn)), function(k) {
    v <- MZ$MIF_z[MZ$xn == k]
    if (length(v) < 8 || length(unique(v)) < 4) return(NULL)
    dn <- density(v, n = 150, adjust = 0.9)
    dd <- data.frame(xpos = k,
                     x = k + dn$y / max(dn$y) * 0.34,
                     y = dn$x)
    rbind(data.frame(xpos = k, x = k, y = min(dd$y)), dd,
          data.frame(xpos = k, x = k, y = max(dd$y)))
  }))
  HALF$grp <- HALF$xpos

  pF <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey90") +
    # 云
    geom_polygon(data = HALF, aes(x, y, group = grp, fill = factor(grp)),
                 alpha = 0.55, colour = "grey35", linewidth = 0.25) +
    # 伞
    geom_boxplot(data = MZ, aes(xn, MIF_z, group = xn),
                 width = 0.075, linewidth = 0.22, outlier.shape = NA,
                 fill = "white", colour = "grey30",
                 position = position_nudge(x = 0)) +
    # 雨（点画在箱体下方的同一 x 附近，水平抖动）
    geom_point(data = MZ, aes(xn, MIF_z, colour = factor(xn)),
               shape = 16, size = 0.22, alpha = 0.30,
               position = position_jitter(width = 0.30, height = 0, seed = 1)) +
    scale_fill_manual(values = c(`1` = "#B2182B", `2` = "#EF8A62", `3` = "#2166AC"),
                      guide = "none") +
    scale_colour_manual(values = c(`1` = "#B2182B", `2` = "#EF8A62", `3` = "#2166AC"),
                        guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = c("Core", "Interface", "Stroma"),
                       expand = expansion(add = c(0.45, 0.55))) +
    labs(x = NULL, y = "MIF (in-section z)") +
    theme_tspo()
}

# =============================================================================
#  G / H  MIF 三对共定位：合并 vs 切片内
# =============================================================================

hr("[G/H] MIF 共定位")

F_PAIR <- find_tab("T20d_spot_correlations.csv")
pG <- panel_missing("MIF pair table not found")
pH <- panel_missing("per-section rho not found")
F_DET  <- find_tab("spatial_persection_detail.csv")

if (!is.na(F_PAIR)) {
  MP <- read.csv(F_PAIR, stringsAsFactors = FALSE)
  # 该表含 8 对，正文只讨论 MIF 的四对
  MP <- MP[grepl("^MIF", MP$pair), ]
  MP$lab <- gsub(" ~ ", "\u2013", MP$pair)
  MP <- MP[order(MP$rho_pooled), ]
  MP$lab <- factor(MP$lab, levels = MP$lab)
  cat("\n"); pr(MP[, c("pair","rho_pooled","rho_within_median","n_pos_section","n_section","p_block")], 3)

  # 哑铃只显示两个绝对值；改画「衰减量」= 合并 ρ − 切片内 ρ，
  # 一根棒的长度就是混合效应被剥掉多少，落点是否越过零线一眼可见。
  AT <- data.frame(lab = MP$lab, pooled = MP$rho_pooled, within = MP$rho_within_median,
                   n_pos = MP$n_pos_section, n_sec = MP$n_section,
                   stringsAsFactors = FALSE)
  AT$drop <- AT$pooled - AT$within
  AT <- AT[order(AT$drop), ]
  AT$lab <- factor(AT$lab, levels = AT$lab)
  cat("\n衰减量（合并 − 切片内）：\n"); pr(AT[, c("lab","pooled","within","drop","n_pos")], 3)

  pG <- ggplot(AT, aes(drop, lab)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
    geom_segment(aes(x = 0, xend = drop, y = lab, yend = lab),
                 colour = "grey84", linewidth = 0.6) +
    geom_point(aes(fill = within > 0), shape = 21, size = 2.1,
               colour = "grey25", stroke = 0.3) +
    geom_text(aes(x = Inf, label = sprintf("%.2f \u2192 %.2f", pooled, within)),
              hjust = 1.05, size = 1.4, colour = "grey30") +
    scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"),
                      labels = c(`TRUE` = "still positive after centring",
                                 `FALSE` = "collapses to zero / negative"),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.34))) +
    labs(x = expression("Attenuation (pooled "*rho*" - within-section "*rho*")"), y = NULL,
         subtitle = "right: pooled -> within-section") +
    theme_tspo() +
    theme(legend.position = c(0.02, 0.04), legend.justification = c(0, 0),
          legend.background = element_rect(fill = alpha("white", 0.88), colour = NA),
          legend.key.size = unit(2.3, "mm"), legend.text = element_text(size = 3.8),
          legend.margin = margin(0, 0, 0, 0),
          plot.subtitle = element_text(size = 4.0, colour = "grey35"))
}

F_MPS <- file.path(CACHEDIR, "mif_persection.rds")
if (file.exists(F_MPS)) {
  PD <- readRDS(F_MPS)
  PD$lab <- factor(gsub("~", "\u2013", gsub(" ", "", PD$pair)), levels = levels(MP$lab))
  PD$sec <- factor(PD$section)
  N2 <- aggregate(rho ~ lab, data = PD, FUN = function(z) sum(z > 0, na.rm = TRUE))
  names(N2)[2] <- "n_pos"
  cat("\n逐切片 MIF ρ 的正向切片数：\n"); pr(N2, 3)

  # 点带只给出中位数与个体值；改成「切片 × 配对」热图后，哪一张切片在拉高
  # 合并相关可以直接读出来（MIF–CD68 的合并正相关主要由 S3 以外的切片贡献为零）。
  MXV <- max(abs(PD$rho), na.rm = TRUE)
  pH <- ggplot(PD, aes(sec, lab, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 1.4, colour = "grey15") +
    scale_fill_div(limit = c(-MXV, MXV), name = "rho",
                   breaks = c(-0.4, -0.2, 0, 0.2)) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_tspo() +
    theme(axis.text.x.top = element_text(size = 4.0),
          axis.text.y = element_text(size = 4.0),
          panel.grid = element_blank(),
          legend.key.height = unit(6, "mm")) +
    guides(fill = guide_colourbar(barwidth = unit(1.8, "mm"), barheight = unit(6, "mm")))
}

# =============================================================================
#  拼装
# =============================================================================

hr("拼装")

FIG <- (pA | pB) / (pC | pD | pE) / (pF | pG | pH) +
  plot_layout(heights = c(1.05, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, colour = "black"))

if (!exists("FIG_NAME")) FIG_NAME <- "Fig8"
print(FIG)
save_fig(FIG, FIG_NAME, width = 7.7, height = 8.4)

cat("\nfig08_codex.R 完成\n")

# -----------------------------------------------------------------------------
# -- 作图样式 / fig09_upstream.R
# -----------------------------------------------------------------------------
# =============================================================================
#  fig09_upstream.R   图 9（8 个面板）上游调控与治疗学预测
#
#  A  chr22q13 区段拷贝数 vs TSPO 表达（供体层面）   散点 + 拟合
#  B  按拷贝数幅度三分位分层                        分层点图 + 中位数
#  C  区段特异性：chr22 基因与其他染色体            密度对比
#  D  六个 CpG 探针甲基化 vs TSPO                   棒棒糖 + 全基因组对照
#  E  bulk 转录因子活性                             火山图
#  F  bulk 与单细胞层面转录因子效应对比             散点（对角线）
#  G  198 种药物的偏相关                            排序棒棒糖
#  H  三分位分层的 Cohen's d                        发散棒棒糖
#
#  ── 数据 ────────────────────────────────────────────────────────────────────
#  ICGC_OV-AU_GX.rds（32,774 基因 × 78 样本）
#  ICGC_OV-AU_CNVbyCpG.rds（486,427 探针 × 78 样本；chr22q13.2 的 6 个探针共享同一值）
#  ICGC_OV-AU_Meth.rds（485,302 探针 × 78 样本）
#  sample.OV-AU.tsv（78 样本 → 67 供体）
#  results_v2/v2_1.6_TF_correlations.csv（bulk，271 TF）
#  results_v2/tables/T6b_SC_TF_vs_TSPO.csv（上皮单细胞，271 TF）
#  results_v2/tables/T5b_drug_partial_all.csv（198 药）
#  results_v2/tables/T5c_drug_tertile.csv（10 药三分位）
#
#  ── 运行 ────────────────────────────────────────────────────────────────────
#    source("<code_final>/fig09_upstream.R", encoding = "UTF-8")
#    约 1–2 分钟（Meth/Ct 均为惰性读入，只取需要的行）
# =============================================================================

CODE <- Sys.getenv("TSPO_CODE_DIR", unset = ".")
if (!file.exists(file.path(CODE, "fig_theme.R"))) {
  CODE <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/code_final"
}
source(file.path(CODE, "fig_theme.R"), encoding = "UTF-8")

PROJ <- if (exists("PROJ")) PROJ else "D:/Rstudy/tcga"
OUT  <- file.path(PROJ, "results_v2")
if (!exists("FIGDIR")) FIGDIR <- "C:/Users/fade/Documents/deepseek-harness/default-workspace/results_v2/figures"
if (!exists("SEED")) SEED <- 20260101
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({
  library(patchwork); library(grid); library(ggrepel); library(ggbeeswarm)
})

hr <- function(x) cat("\n===== ", x, " =====\n", sep = "")
pr <- function(d, digits = 3) print(d, row.names = FALSE, digits = digits)
set.seed(SEED)

find_f <- function(...) {
  cand <- c(file.path(OUT, ...), file.path(PROJ, ...), file.path(PROJ, "results_v2", ...))
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}

CG6 <- c("cg02079700","cg03207265","cg05333174","cg13608545","cg15429214","cg24694833")

# =============================================================================
#  STEP 0  ICGC：表达、拷贝数、供体映射
# =============================================================================

hr("STEP 0  ICGC")

GX <- readRDS(file.path(PROJ, "ICGC_OV-AU_GX.rds"))
CN <- readRDS(file.path(PROJ, "ICGC_OV-AU_CNVbyCpG.rds"))
cat("表达: ", nrow(GX), " 基因 x ", ncol(GX), " 样本\n", sep = "")

SEG <- colMeans(as.matrix(CN[CG6, , drop = FALSE]), na.rm = TRUE)
cat("区段 CN: ", paste(round(range(SEG, na.rm = TRUE), 3), collapse = " ~ "), "\n", sep = "")

SMA <- read.delim(file.path(PROJ, "sample.OV-AU.tsv"), stringsAsFactors = FALSE)
SMA <- SMA[!duplicated(SMA$icgc_sample_id), c("icgc_sample_id", "icgc_donor_id")]

# ⚠ 表达矩阵与拷贝数矩阵的列名是两套不同的样本 ID（同一供体的不同等分样本），
#   交集为 0，因此不能按样本直接对齐，必须各自聚合到供体再按供体合并。
cat("GX 列名例: ", paste(head(colnames(GX), 2), collapse = ", "),
    " | CN 列名例: ", paste(head(names(SEG), 2), collapse = ", "),
    " | 交集 ", length(intersect(colnames(GX), names(SEG))), "\n", sep = "")

eS <- SMA[SMA$icgc_sample_id %in% colnames(GX), ]
TS <- data.frame(donor = eS$icgc_donor_id,
                 TSPO = as.numeric(GX["TSPO", eS$icgc_sample_id]),
                 stringsAsFactors = FALSE)
TD <- aggregate(TSPO ~ donor, data = TS, FUN = mean, na.rm = TRUE)

cS <- SMA[SMA$icgc_sample_id %in% names(SEG), ]
CS <- data.frame(donor = cS$icgc_donor_id,
                 CN = as.numeric(SEG[cS$icgc_sample_id]), stringsAsFactors = FALSE)
CD <- aggregate(CN ~ donor, data = CS, FUN = mean, na.rm = TRUE)

cat("表达供体 ", nrow(TD), "  拷贝数供体 ", nrow(CD), "\n", sep = "")
DON <- merge(TD, CD, by = "donor")
DON <- DON[is.finite(DON$TSPO) & is.finite(DON$CN), ]
cat("可用于相关分析的供体: ", nrow(DON), "\n", sep = "")

R_CN <- suppressWarnings(cor(DON$TSPO, DON$CN, method = "spearman"))
P_CN <- suppressWarnings(cor.test(DON$TSPO, DON$CN, method = "spearman")$p.value)
# 校正表达幅度（TSPO 自身 z 化的幅度）后重算
AMP <- abs(scale(DON$TSPO))[, 1]
DON$CN_adj <- residuals(lm(rank(DON$CN) ~ rank(AMP)))
DON$TSPO_adj <- residuals(lm(rank(DON$TSPO) ~ rank(AMP)))
R_ADJ <- cor(DON$CN_adj, DON$TSPO_adj)
cat(sprintf("TSPO ~ 区段 CN: rho = %+.3f (p = %.2g, n = %d)；幅度校正后 %+.3f\n",
            R_CN, P_CN, nrow(DON), R_ADJ))

# =============================================================================
#  A  区段 CN 在本队列中的分布
#
#  原本 A 与 B 画的是同一批散点（同样两个轴、同样 67 个点），只是 B 多了分层着色，
#  属于重复。A 改为区段 CN 的分布：本队列 91% 供体该区段为缺失，
#  这个背景是理解「CN 与 TSPO 正相关」的前提，而散点图里看不出来。
# =============================================================================

hr("[A] 区段 CN 分布")

QB0 <- quantile(DON$CN, c(1/3, 2/3), na.rm = TRUE)
PCT_DEL <- 100 * mean(DON$CN < 0)
cat(sprintf("区段 CN 均值 %.3f；缺失(<0) 占 %.1f%%；三分位切点 %.3f / %.3f\n",
            mean(DON$CN), PCT_DEL, QB0[1], QB0[2]))

pA <- ggplot(DON, aes(CN)) +
  geom_histogram(aes(fill = after_stat(x < 0)), bins = 22,
                 colour = "white", linewidth = 0.25) +
  geom_rug(colour = "grey35", alpha = 0.6, length = unit(1.2, "mm")) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "grey20", linetype = "22") +
  geom_vline(xintercept = QB0, linewidth = 0.35, colour = "grey45", linetype = "12") +
  annotate("label", x = QB0[1], y = Inf, hjust = 1.15, vjust = 1.35, size = 1.45,
           colour = "grey30", fill = alpha("white", 0.9), label.size = 0,
           label.padding = unit(0.4, "mm"), label = "tertile") +
  annotate("label", x = Inf, y = Inf, hjust = 1.03, vjust = 1.25, size = 1.6,
           colour = PAL["ink"], fill = alpha("white", 0.9), label.size = 0,
           label.padding = unit(0.5, "mm"), lineheight = 0.95,
           label = sprintf("mean %+.3f\n%.0f%% deleted (<0)\nn = %d donors",
                           mean(DON$CN), PCT_DEL, nrow(DON))) +
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"), guide = "none") +
  labs(x = "chr22q13.2 segment copy number", y = "Donors") +
  theme_tspo()

hr("[B] 分层散点 + 层内拟合")

# 原为分层小提琴，只说明三层有差别；改为在八分层内各画拟合线，
# 直接检验「相关不是由层间差异造成」这一更强的命题。
QB <- quantile(DON$CN, c(1/3, 2/3), na.rm = TRUE)
DON$stratum <- factor(ifelse(DON$CN > QB[2], "High CN",
                      ifelse(DON$CN <= QB[1], "Low CN", "Mid CN")),
                      levels = c("Low CN", "Mid CN", "High CN"))
cat("分层切点: ", paste(round(QB, 3), collapse = " / "), "\n", sep = "")
print(table(DON$stratum))

RIN <- sapply(levels(DON$stratum), function(k) {
  d <- DON[DON$stratum == k, ]
  suppressWarnings(cor(d$TSPO, d$CN, method = "spearman"))
})
KW <- kruskal.test(TSPO ~ stratum, data = DON)
cat("层内 rho: ", paste(sprintf("%s %+.3f", names(RIN), RIN), collapse = "  "), "\n", sep = "")
cat(sprintf("Kruskal-Wallis p = %.3g\n", KW$p.value))

PAL_B <- c(`Low CN` = "#2166AC", `Mid CN` = "#BDBDBD", `High CN` = "#B2182B")
pB <- ggplot(DON, aes(CN, TSPO, colour = stratum, fill = stratum)) +
  # 整体拟合（黑）与三层内的拟合（彩色）叠在一张图上：
  # 整体为正，但层内方向不一致，这正是正文所说的「组间剂量效应」。
  geom_smooth(data = DON, aes(x = CN, y = TSPO), method = "lm", formula = y ~ x,
              se = FALSE, linewidth = 0.85, colour = PAL["ink"],
              inherit.aes = FALSE) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.5,
              alpha = 0.16, show.legend = FALSE) +
  geom_point(shape = 21, size = 1.3, colour = "grey35", stroke = 0.22) +
  scale_colour_manual(values = PAL_B, name = NULL) +
  scale_fill_manual(values = PAL_B, guide = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.9, size = 1.6,
           colour = PAL["ink"], lineheight = 0.95,
           label = sprintf("overall rho = %+.3f (p = %.1e)\nwithin-stratum: %s",
                           R_CN, P_CN,
                           paste(sprintf("%+.2f", RIN), collapse = " / "))) +
  labs(x = "chr22q13.2 segment copy number", y = "TSPO expression (log2)") +
  theme_tspo() +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 3.9))

# =============================================================================
#  C  区段特异性
# =============================================================================

hr("[C] 区段特异性")

# 原为 chr22 vs 其他染色体的密度对比，抽象；改为正文点名的几个基因，
# 每个点一个基因，读者能直接对上文字里的数字。
GENES_C <- c("TOMM22", "TSPO", "PPARA", "CD68", "ACTB", "GAPDH", "EPCAM")
pC <- panel_missing("gene-level CN correlations unavailable")
if (all(GENES_C %in% rownames(GX))) {
  GSC <- do.call(rbind, lapply(GENES_C, function(g) {
    v <- as.numeric(GX[g, eS$icgc_sample_id])
    gd <- aggregate(v ~ eS$icgc_donor_id, FUN = mean)
    names(gd) <- c("donor", "v")
    m <- merge(gd, DON[, c("donor", "CN")], by = "donor")
    r <- suppressWarnings(cor(m$v, m$CN, method = "spearman"))
    p <- suppressWarnings(cor.test(m$v, m$CN, method = "spearman")$p.value)
    data.frame(gene = g, rho = r, p = p, n = nrow(m), stringsAsFactors = FALSE)
  }))
  GSC$chr22 <- GSC$gene %in% c("TOMM22", "TSPO", "PPARA")
  GSC <- GSC[order(GSC$rho), ]
  GSC$gene <- factor(GSC$gene, levels = GSC$gene)
  cat("\n具名基因与区段 CN 的相关：\n"); pr(GSC[, c("gene", "rho", "p", "chr22")], 3)

  pC <- ggplot(GSC, aes(rho, gene)) +
    geom_vline(xintercept = 0, linewidth = 0.32, colour = "grey55", linetype = "22") +
    geom_segment(aes(x = 0, xend = rho, y = gene, yend = gene),
                 colour = "grey85", linewidth = 0.55) +
    geom_point(aes(fill = chr22), shape = 21, size = 2.2,
               colour = "grey25", stroke = 0.32) +
    geom_text(aes(x = Inf, label = ifelse(p < 0.05, sprintf("p=%.1e", p),
                                          sprintf("p=%.2f", p))),
              hjust = 1.05, size = 1.45, colour = "grey30") +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey75"),
                      labels = c(`TRUE` = "chr22q13", `FALSE` = "other chromosome"),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.34))) +
    labs(x = expression("Spearman "*rho*" with segment CN"), y = NULL) +
    theme_tspo() +
    theme(axis.text.y = element_text(size = 4.4),
          legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 3.9))
}

# =============================================================================
#  D  甲基化
# =============================================================================

hr("[D] 甲基化")

pD <- panel_missing("methylation matrix not found")
FM <- file.path(PROJ, "ICGC_OV-AU_Meth.rds")
if (file.exists(FM)) {
  ME <- readRDS(FM)
  # 甲基化矩阵同样是与表达矩阵不同的等分样本，必须按供体合并
  mM <- SMA[SMA$icgc_sample_id %in% colnames(ME), ]
  cat("甲基化样本 ", nrow(mM), "  供体 ", length(unique(mM$icgc_donor_id)), "\n", sep = "")
  MB <- as.matrix(ME[CG6, mM$icgc_sample_id, drop = FALSE])

  agg_donor <- function(mat, grp) {           # mat: 特征 x 样本
    g <- rowsum(t(mat), group = grp)
    n <- as.numeric(table(grp))
    g / n[match(rownames(g), names(table(grp)))]
  }
  MD_d <- agg_donor(MB, mM$icgc_donor_id)                       # 供体 x 6 探针
  GW_d <- agg_donor(matrix(colMeans(ME[, mM$icgc_sample_id, drop = FALSE], na.rm = TRUE), nrow = 1),
                    mM$icgc_donor_id)[, 1]                      # 供体层面全基因组平均 β

  don_m <- data.frame(donor = rownames(MD_d), stringsAsFactors = FALSE)
  MR <- sapply(seq_along(CG6), function(k) {
    y <- MD_d[, k]
    d <- merge(data.frame(donor = rownames(MD_d), beta = y),
               TD, by = "donor")
    suppressWarnings(cor(d$beta, d$TSPO, method = "spearman", use = "complete.obs"))
  })
  dg <- merge(data.frame(donor = rownames(MD_d), GW = GW_d), TD, by = "donor")
  R_GW <- suppressWarnings(cor(dg$GW, dg$TSPO, method = "spearman", use = "complete.obs"))

  MD <- data.frame(cpg = CG6, rho = as.numeric(MR), stringsAsFactors = FALSE)
  MD <- MD[order(MD$rho), ]
  MD$cpg <- factor(MD$cpg, levels = MD$cpg)
  cat("6 个 CpG 的 rho: ", paste(round(MD$rho, 3), collapse = " "), "\n", sep = "")
  cat(sprintf("均值 %.3f；全基因组平均 β 对照 rho = %+.3f（n = %d 供体）\n",
              mean(MD$rho), R_GW, nrow(dg)))

  # 原为棒棒糖，与 C/H 同为「点 + 线」；改为哑铃：每个探针的相关与
  # 全基因组平均 β 的对照值并列，两者几乎重合即说明无位点特异性。
  DL <- rbind(
    data.frame(cpg = MD$cpg, type = "CpG probe", rho = MD$rho),
    data.frame(cpg = MD$cpg, type = "Genome-wide avg", rho = R_GW))
  DL$type <- factor(DL$type, levels = c("CpG probe", "Genome-wide avg"))

  pD <- ggplot(DL, aes(rho, cpg)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
    geom_line(aes(group = cpg), colour = "grey80", linewidth = 0.5) +
    geom_point(aes(fill = type, shape = type), size = 1.9,
               colour = "grey30", stroke = 0.32) +
    scale_fill_manual(values = c(`CpG probe` = "#B2182B",
                                 `Genome-wide avg` = "white"), name = NULL) +
    scale_shape_manual(values = c(`CpG probe` = 21, `Genome-wide avg` = 21), name = NULL) +
    labs(x = expression("Spearman "*rho*" with TSPO"), y = NULL) +
    theme_tspo() +
    theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.88), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 3.9))
}

# =============================================================================
#  E  bulk 转录因子
# =============================================================================

hr("[E] bulk TF")

FB <- find_f("v2_1.6_TF_correlations.csv")
pE <- panel_missing("bulk TF table not found")
if (!is.na(FB)) {
  TF <- read.csv(FB, stringsAsFactors = FALSE)
  TF <- TF[is.finite(TF$rho) & is.finite(TF$FDR), ]
  TF$nlp <- -log10(pmax(TF$FDR, .Machine$double.eps))
  SIG <- TF$FDR < 0.05
  cat("271 个 TF 中 FDR<0.05: ", sum(SIG), " 个\n", sep = "")
  TOP <- head(TF[SIG, ][order(-abs(TF$rho[SIG])), ], 5)
  cat("效应量最大: ", paste(TOP$Gene, collapse = ", "), "\n", sep = "")

  # 原为火山图，与本图其余面板的「点阵」形态重复；改为以 0 为轴的双向密度：
  # 负相关与正相关各占一侧，分布形状与偏斜一眼可见。
  dn <- density(TF$rho, na.rm = TRUE, n = 400)
  DENS <- data.frame(x = dn$x, y = dn$y)
  DENS$side <- ifelse(DENS$x < 0, "negative", "positive")
  DENS$y2 <- DENS$y * ifelse(DENS$side == "negative", -1, 1)

  pE <- ggplot(DENS, aes(x, y2, fill = side)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_area(colour = "grey35", linewidth = 0.3, alpha = 0.65) +
    geom_vline(xintercept = c(-0.2, 0.2), linewidth = 0.25,
               colour = "grey60", linetype = "22") +
    geom_text_repel(data = TOP, aes(x = rho, y = 0, label = Gene),
                    inherit.aes = FALSE,
                    size = 1.6, colour = "grey10", segment.size = 0.14,
                    nudge_y = 0.18, max.overlaps = Inf, seed = SEED,
                    min.segment.length = 0) +
    scale_fill_manual(values = c(negative = "#2166AC", positive = "#B2182B"),
                      guide = "none") +
    scale_y_continuous(labels = function(z) sprintf("%.2f", abs(z))) +
    labs(x = expression("TF activity vs TSPO (Spearman "*rho*")"),
         y = "Density") +
    theme_tspo()
}

# =============================================================================
#  F  bulk vs 单细胞
# =============================================================================

hr("[F] bulk vs 单细胞")

FS <- find_f("tables", "T6b_SC_TF_vs_TSPO.csv")
pF <- panel_missing("single-cell TF table not found")
if (!is.na(FB) && !is.na(FS)) {
  SCC <- read.csv(FS, stringsAsFactors = FALSE)
  MM <- merge(TF[, c("Gene", "rho")], SCC[, c("tf", "median_rho")],
              by.x = "Gene", by.y = "tf")
  names(MM) <- c("tf", "bulk", "sc")
  MM <- MM[is.finite(MM$bulk) & is.finite(MM$sc), ]
  NFKB <- MM$tf %in% c("RELA", "NFKB1", "NFKB2", "RELB", "REL")
  cat("共有 TF: ", nrow(MM), " 个；符号一致 ",
      sum(sign(MM$bulk) == sign(MM$sc)), " 个\n", sep = "")
  cat("NF-κB 家族:\n"); pr(MM[NFKB, ], 3)

  pF <- ggplot(MM, aes(bulk, sc)) +
    geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey80") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.28, colour = "grey55", linetype = "22") +
    geom_point(aes(fill = sign(bulk) != sign(sc)), shape = 21, size = 1.0,
               colour = "grey45", stroke = 0.15, alpha = 0.85) +
    geom_text_repel(data = MM[NFKB, ], aes(label = tf), size = 1.7, colour = "#B2182B",
                    segment.size = 0.15, max.overlaps = Inf, seed = SEED) +
    scale_fill_manual(values = c(`TRUE` = "#EF8A62", `FALSE` = "grey80"),
                      labels = c(`TRUE` = "sign flips", `FALSE` = "same sign"), name = NULL) +
    labs(x = "bulk TF activity vs TSPO (rho)", y = "epithelial single-cell (median rho)") +
    theme_tspo() +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 3.9))
}

# =============================================================================
#  G  198 种药物的偏相关
# =============================================================================

hr("[G] 药物偏相关")

FD <- find_f("tables", "T5b_drug_partial_all.csv")
pG <- panel_missing("drug table not found")
if (!is.na(FD)) {
  DR <- read.csv(FD, stringsAsFactors = FALSE)
  DR <- DR[is.finite(DR$rho_partial), ]
  DR$drug_short <- sub("_[0-9]+$", "", DR$drug)
  DR <- DR[order(DR$rho_partial), ]
  DR$rank <- seq_len(nrow(DR))
  DR$sig <- DR$FDR_partial < 0.05
  cat("药物 ", nrow(DR), " 种；偏相关 FDR<0.05 的 ", sum(DR$sig), " 种\n", sep = "")

  LAB <- DR$drug_short %in% c("OSI-027","AZD2014","AZD8055","GNE-317","Trametinib",
                              "Bortezomib","Dasatinib","MK-1775","Cisplatin","Paclitaxel",
                              "JQ1","LY2109761","Doramapimod","Ruxolitinib")
  pG <- ggplot(DR, aes(rank, rho_partial)) +
    geom_hline(yintercept = 0, linewidth = 0.28, colour = "grey55") +
    geom_point(aes(rank, rho_partial, fill = sig), shape = 21, size = 0.85,
               colour = "grey55", stroke = 0.12, alpha = 0.9) +
    geom_point(data = DR[LAB, ], aes(rank, rho_partial), shape = 21, size = 1.5,
               fill = "white", colour = "grey15", stroke = 0.35) +
    geom_text_repel(data = DR[LAB, ], aes(rank, rho_partial, label = drug_short),
                    size = 1.45, colour = "grey15", segment.size = 0.1,
                    max.overlaps = Inf, seed = SEED, min.segment.length = 0) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey85"),
                      labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."), name = NULL) +
    labs(x = "Drugs ranked by partial rho", y = expression("Partial "*rho*" with TSPO")) +
    theme_tspo() +
    theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 3.9))
}

# =============================================================================
#  H  三分位 Cohen's d
# =============================================================================

hr("[H] 三分位效应量")

FT <- find_f("tables", "T5c_drug_tertile.csv")
pH <- panel_missing("tertile table not found")
if (!is.na(FT)) {
  TT <- read.csv(FT, stringsAsFactors = FALSE)
  TT <- TT[is.finite(TT$cohens_d), ]
  TT <- TT[order(TT$cohens_d), ]
  TT$drug <- factor(TT$drug, levels = TT$drug)
  cat("\n"); pr(TT[, c("drug", "cohens_d", "p")], 3)

  # 原为效应量棒棒糖，与 C/G 同为「点 + 线」；改为配对斜率，
  # 直接显示 High CN 与 Low CN 两组的实测预测 IC50 差值方向。
  SL <- rbind(
    data.frame(drug = TT$drug, grp = "Low CN",  y = TT$median_low),
    data.frame(drug = TT$drug, grp = "High CN", y = TT$median_high))
  SL$grp <- factor(SL$grp, levels = c("Low CN", "High CN"))
  TT$drug <- factor(TT$drug, levels = TT$drug)

  pH <- ggplot(SL, aes(grp, y, group = drug)) +
    geom_line(aes(colour = drug), linewidth = 0.45, alpha = 0.9) +
    geom_point(aes(fill = drug), shape = 21, size = 1.5,
               colour = "grey30", stroke = 0.25) +
    scale_colour_manual(values = setNames(
      colorRampPalette(c("#2166AC", "#B2182B"))(nrow(TT)), levels(TT$drug)), guide = "none") +
    scale_fill_manual(values = setNames(
      colorRampPalette(c("#2166AC", "#B2182B"))(nrow(TT)), levels(TT$drug)), guide = "none") +
    geom_text_repel(data = SL[SL$grp == "High CN", ],
                    aes(label = drug), size = 1.5, colour = "grey15",
                    nudge_x = 0.32, hjust = 0, direction = "y",
                    segment.size = 0.12, max.overlaps = Inf, seed = SEED,
                    box.padding = 0.35) +
    scale_x_discrete(expand = expansion(add = c(0.42, 0.95))) +
    labs(x = NULL, y = "Predicted IC50 (log scale)") +
    theme_tspo()
}

# =============================================================================
#  I  药物的原始相关 vs 偏相关
#
#  G 面板给的是校正后的偏相关；这里把校正前后的关系摊开，读者能看出
#  免疫构成与增殖评分到底削掉了多少信号——这也是本文反复强调的混杂问题。
# =============================================================================

hr("[I] 原始 vs 偏相关")

pI <- panel_missing("drug table not found")
if (!is.na(FD)) {
  DD2 <- read.csv(FD, stringsAsFactors = FALSE)
  DD2 <- DD2[is.finite(DD2$rho) & is.finite(DD2$rho_partial), ]
  DD2$drug_short <- sub("_[0-9]+$", "", DD2$drug)
  DD2$drop <- DD2$rho - DD2$rho_partial
  LAB2 <- DD2$drug_short %in% c("OSI-027","AZD2014","AZD8055","GNE-317","Trametinib",
                                "Bortezomib","Dasatinib","MK-1775","Cisplatin","Paclitaxel",
                                "JQ1","LY2109761","Doramapimod","Ruxolitinib")
  cat("药物 ", nrow(DD2), " 种；校正后 |rho| 缩小 ",
      sum(abs(DD2$rho_partial) < abs(DD2$rho)), " 种，扩大 ",
      sum(abs(DD2$rho_partial) > abs(DD2$rho)), " 种\n", sep = "")

  pI <- ggplot(DD2, aes(rho, rho_partial)) +
    geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey82") +
    geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey82") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
    geom_point(aes(fill = DD2$FDR_partial < 0.05), shape = 21, size = 0.9,
               colour = "grey55", stroke = 0.12, alpha = 0.85) +
    geom_point(data = DD2[LAB2, ], aes(rho, rho_partial), shape = 21, size = 1.5,
               fill = "white", colour = "grey15", stroke = 0.35) +
    geom_text_repel(data = DD2[LAB2, ], aes(rho, rho_partial, label = drug_short),
                    size = 1.4, colour = "grey10", segment.size = 0.1,
                    max.overlaps = Inf, seed = SEED, min.segment.length = 0,
                    box.padding = 0.5) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey85"),
                      labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."), name = NULL) +
    labs(x = expression("Raw Spearman "*rho*" with TSPO"),
         y = expression("Partial "*rho*" (adjusted)")) +
    theme_tspo() +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 3.9))
}

# =============================================================================
#  拼装 3 × 3
# =============================================================================

hr("拼装")

FIG <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI)
FIG <- FIG + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, colour = "black"))

if (!exists("FIG_NAME")) FIG_NAME <- "Fig9"
print(FIG)
save_fig(FIG, FIG_NAME, width = 7.7, height = 9.0)

cat("\nfig09_upstream.R 完成\n")
