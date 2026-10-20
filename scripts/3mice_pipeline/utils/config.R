# ============================================================================
# 3mice_pipeline/utils/config.R — 三品系（BALB/c vs C57BL/6）绘图流程配置
# 从仓库根目录运行；数据默认放在 data/raw/3mice/
# ============================================================================
ROOT <- normalizePath(Sys.getenv("MPXV_PROJECT_ROOT", "."), winslash = "/", mustWork = FALSE)
DATA_DIR_3MICE <- normalizePath(
  Sys.getenv("MPXV_DATA_DIR_3MICE", file.path(ROOT, "data", "raw", "3mice")),
  winslash = "/", mustWork = FALSE)
OUT_DIR <- normalizePath(
  Sys.getenv("MPXV_OUT_3MICE", file.path(ROOT, "results", "3mice")),
  winslash = "/", mustWork = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
