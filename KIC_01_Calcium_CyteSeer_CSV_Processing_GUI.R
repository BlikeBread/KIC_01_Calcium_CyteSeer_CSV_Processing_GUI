###############################################################################
# KIC CALCIUM ANALYSIS PIPELINE — STEP 1
# -----------------------------------------------------------------------------
#
# PURPOSE:
#   Interactive GUI tool to process CyteSeer-exported KIC calcium data (CSV).
#
#   The script:
#     1. Prompts the user to select input/output folders (tk_choose.dir)
#     2. Recursively scans the input directory for CyteSeerResults* exports
#     3. Reads PRIMARY traces (Primary_DataTable_CStats.csv) and parses KIC
#        waveform strings to reconstruct time series
#     4. Reads CELL peak tables (_Cell Peak Measurements_DataTable.csv) and
#        computes missing upstroke/downstroke timing metrics (10→25/50/75 and
#        90→75/50/25), using waveform-based estimation when available
#     5. Reads WHOLE-IMAGE traces (_Whole Image Measurements_DataTable_CStats.csv)
#        and parses KIC waveform strings per well
#     6. Reads WHOLE-IMAGE peak tables (_Whole-Image Peak Measurements_DataTable.csv)
#        and computes the same additional timing metrics
#     7. Exports per-well RAW trace plots (JPG) and last-10s trace windows (CSV)
#        for both CELL and WHOLE-IMAGE signals
#     8. Exports CELL peak tables (per-peak) plus per-cell mean summaries
#     9. Exports WHOLE-IMAGE peak tables (per-peak) plus per-well summary tables
#        (FULL stats + MEANS-only)
#
# OUTPUT:
#   An output folder containing:
#     /raw_traces                  (CELL trace plots, JPG)
#     /raw_csv                     (CELL last-10s traces, CSV)
#     /cell_peak_parameters        (CELL peak tables + per-cell means, CSV)
#     /whole_image_traces          (WHOLE-IMAGE trace plots, JPG)
#     /whole_image_csv             (WHOLE-IMAGE last-10s traces, CSV)
#     /whole_image_parameters      (WHOLE-IMAGE peak tables, CSV)
#     /summary                     (WHOLE-IMAGE summaries FULL + MEANS-only, CSV)
#
# NOTES:
#   - Uses tcltk (tk_choose.dir) for macOS folder pickers.
#   - Designed for interactive use (GUI folder selection + automated exports).
#   - Intended as Step 1 of a 3-step KIC calcium processing pipeline.
#
# AUTHORS:
#   Michele Buono
#   Talitha Spanjersberg
#   Nikki Scheen
#   Nina van der Wilt
#   Regenerative Medicine Center Utrecht (2026)
###############################################################################

# ===================================================================================================================
# 01) LOAD PACKAGES & BASIC SETUP
# -------------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)      # ensure dplyr::filter() masks stats::filter()
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tcltk)
})

prefix <- format(Sys.Date(), "%d.%m.%Y")

# ===================================================================================================================
# 02) CHOOSE INPUT/OUTPUT FOLDERS (tk_choose.dir)
# -------------------------------------------------------------------------------------------------------------------
raw_in <- tk_choose.dir(caption = "Select INPUT folder (must contain CyteSeerResults...)")
if (is.na(raw_in) || !nzchar(raw_in)) stop("❌ No input folder selected.")

output <- tk_choose.dir(caption = "Select OUTPUT folder (results destination)")
if (is.na(output) || !nzchar(output)) stop("❌ No output folder selected.")

input_dir  <- path.expand(raw_in)
output_dir <- path.expand(output)

message("▶ input_dir  = ", input_dir,  " (exists? ", dir.exists(input_dir),  ")")
message("▶ output_dir = ", output_dir, " (exists? ", dir.exists(output_dir), ")")

if (!dir.exists(input_dir)) stop("❌ input_dir not found – check the selected folder.")
if (!dir.exists(output_dir)) {
  message("⚠️ output_dir not found; creating it now.")
  dir.create(output_dir, recursive = TRUE)
}

# ===================================================================================================================
# 03) PRIMARY TRACES: FIND & READ A REFERENCE CSV (trace_data)
# -------------------------------------------------------------------------------------------------------------------
files <- list.files(
  path = input_dir,
  pattern = "Primary_DataTable_CStats.csv$",
  full.names = TRUE,
  recursive = TRUE
)
files <- files[str_detect(files, "CyteSeerResults[^/\\\\]*[/\\\\]")]

if (length(files) == 0) {
  stop("❌ No primary trace CSVs found in folders starting with 'CyteSeerResults'.")
}

# NOTE: The skip/nrows here are tuned to your export format. Adjust if your export changes.
trace_data <- read.csv(files[1], skip = 2747, nrows = 50, fill = TRUE, header = TRUE)

# ===================================================================================================================
# 04) CELL PEAKS: FIND, READ & AUGMENT
# -------------------------------------------------------------------------------------------------------------------
peak_files <- list.files(
  path = input_dir,
  pattern = "_Cell Peak Measurements_DataTable.csv$",
  full.names = TRUE,
  recursive = TRUE
)
peak_files <- peak_files[str_detect(peak_files, "CyteSeerResults[^/\\\\]*[/\\\\]")]

message("🔍 Found ", length(peak_files), " peak data file(s).")
print(peak_files)

if (length(peak_files) == 0) {
  stop("❌ No peak data files found in expected CyteSeerResults folders.")
}

safe_read_peak <- function(file, skip_lines = 82) {
  df <- tryCatch(
    read.csv(file,
             skip = skip_lines,
             fill = TRUE,
             header = TRUE,
             check.names = TRUE,     # converts spaces to dots for safe column access
             stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  if (nrow(df) == 0 || ncol(df) == 0) {
    message("⚠️ Skipping empty or unreadable file: ", basename(file))
    return(df)
  }
  df
}

peak_data <- purrr::map_dfr(peak_files, function(file) {
  well_name <- basename(dirname(file))
  df <- safe_read_peak(file, skip_lines = 82)
  if (nrow(df) == 0) return(tibble::tibble())
  df$Well <- rep(well_name, nrow(df))
  tibble::as_tibble(df, .name_repair = "minimal")
})

# -------- Helpers (parsers & math) ----
vnum <- function(x) suppressWarnings(as.numeric(x))

parse_kic_series <- function(s) {
  if (is.na(s) || !nzchar(s)) return(numeric(0))
  toks <- unlist(strsplit(gsub("[{}]", "", s), "[;|,\\s]+"))
  x <- suppressWarnings(as.numeric(toks))
  x[is.finite(x)]
}

safe_dt_ms <- function(x) {
  dt <- suppressWarnings(as.numeric(x))
  if (!is.finite(dt)) return(NA_real_)
  if (dt < 1) dt * 1000 else dt
}

# robust downstroke crossing using a monotone envelope after peak
cross_down <- function(yN, t, level, pk) {
  seg  <- yN[pk:length(yN)]
  segM <- cummin(pmin(seg, 1))
  ts   <- t[pk:length(yN)]
  idx  <- which(segM <= level)[1]
  if (is.na(idx)) return(NA_real_)
  if (idx == 1) return(ts[1])
  j <- idx - 1L
  ts[j] + (level - segM[j]) * (ts[j+1]-ts[j]) / (segM[j+1]-segM[j])
}

# interpolate among available (p->90) spans
predict_span_vec <- function(p, spans_list, perc_vec) {
  if (!length(perc_vec)) return(rep(NA_real_, length(spans_list[[1]])))
  key <- as.character(p)
  if (key %in% names(spans_list)) return(spans_list[[key]])
  percs <- sort(perc_vec); n <- length(percs)
  if (n < 2) return(rep(NA_real_, length(spans_list[[1]])))
  if (p <= percs[1])       { l <- percs[1];     u <- percs[2] }
  else if (p >= percs[n])  { l <- percs[n-1];   u <- percs[n] }
  else                     { l <- max(percs[percs < p]); u <- min(percs[percs > p]) }
  yL <- spans_list[[as.character(l)]]
  yU <- spans_list[[as.character(u)]]
  yL + (p - l) / (u - l) * (yU - yL)
}

# ---------- CELL peaks: compute 10->25/50/75 and 90->75/50/25 ----------
series_cell <- intersect(c("Cell.Peak.Smoothed.F.t.", "Cell.Peak.F.t."), names(peak_data))
series_cell <- if (length(series_cell)) series_cell[1] else NA_character_

ts_candidates_cell <- intersect(c("Timeslice", "Cell.Peak.Timeslice"), names(peak_data))
ts_cell <- if (length(ts_candidates_cell)) ts_candidates_cell[1] else NA_character_

if (is.na(series_cell)) message("⚠️ Cell peaks: no waveform column found (Cell.Peak.[Smoothed.]F.t.).")
if (is.na(ts_cell))     message("⚠️ Cell peaks: no timeslice column found (Timeslice / Cell.Peak.Timeslice).")

for (nm in c("Cell.Peak.Upstroke.10.25.Time","Cell.Peak.Upstroke.10.50.Time","Cell.Peak.Upstroke.10.75.Time",
             "Cell.Peak.Downstroke.90.75.Time","Cell.Peak.Downstroke.90.50.Time","Cell.Peak.Downstroke.90.25.Time")) {
  if (!nm %in% names(peak_data)) peak_data[[nm]] <- NA_real_
}

U1090 <- if ("Cell.Peak.Upstroke.10.90.Time" %in% names(peak_data)) vnum(peak_data[["Cell.Peak.Upstroke.10.90.Time"]]) else NA_real_

up_nms <- grep("^Cell\\.Peak\\.Upstroke\\.(\\d{2})\\.90\\.Time$", names(peak_data), value = TRUE)
if (length(up_nms)) {
  up_perc  <- as.numeric(sub("^.*Upstroke\\.(\\d{2})\\.90\\.Time$", "\\1", up_nms))
  up_spans <- setNames(lapply(up_nms, \(nm) vnum(peak_data[[nm]])), as.character(up_perc))
  up_perc  <- sort(unique(up_perc))
  
  U2590 <- predict_span_vec(25, up_spans, up_perc)
  U5090 <- predict_span_vec(50, up_spans, up_perc)
  U7590 <- predict_span_vec(75, up_spans, up_perc)
  
  peak_data[["Cell.Peak.Upstroke.10.25.Time"]] <- ifelse(is.finite(U1090), U1090 - U2590, NA_real_)
  peak_data[["Cell.Peak.Upstroke.10.50.Time"]] <- ifelse(is.finite(U1090), U1090 - U5090, NA_real_)
  peak_data[["Cell.Peak.Upstroke.10.75.Time"]] <- ifelse(is.finite(U1090), U1090 - U7590, NA_real_)
}

# Downstroke preferred: waveform monotone envelope
if (!is.na(series_cell) && !is.na(ts_cell)) {
  for (i in seq_len(nrow(peak_data))) {
    y <- parse_kic_series(peak_data[[series_cell]][i]); if (length(y) < 5) next
    dt <- safe_dt_ms(peak_data[[ts_cell]][i]);          if (!is.finite(dt) || dt <= 0) next
    base <- stats::quantile(y, 0.05, na.rm = TRUE); A <- max(y, na.rm = TRUE) - base
    if (!is.finite(A) || A <= 0) next
    yN <- (y - base) / A; t <- (seq_along(yN)-1) * dt; pk <- which.max(yN)
    
    t90 <- cross_down(yN, t, 0.90, pk)
    t75 <- cross_down(yN, t, 0.75, pk)
    t50 <- cross_down(yN, t, 0.50, pk)
    t25 <- cross_down(yN, t, 0.25, pk)
    
    if (is.finite(t90) && is.finite(t75)) peak_data[["Cell.Peak.Downstroke.90.75.Time"]][i] <- t75 - t90
    if (is.finite(t90) && is.finite(t50)) peak_data[["Cell.Peak.Downstroke.90.50.Time"]][i] <- t50 - t90
    if (is.finite(t90) && is.finite(t25)) peak_data[["Cell.Peak.Downstroke.90.25.Time"]][i] <- t25 - t90
  }
}

# Downstroke fallback: infer from existing spans
c_9020 <- if ("Cell.Peak.Downstroke.90.20.Time" %in% names(peak_data)) vnum(peak_data[["Cell.Peak.Downstroke.90.20.Time"]]) else NA_real_
c_9030 <- if ("Cell.Peak.Downstroke.90.30.Time" %in% names(peak_data)) vnum(peak_data[["Cell.Peak.Downstroke.90.30.Time"]]) else NA_real_
c_9040 <- if ("Cell.Peak.Downstroke.90.40.Time" %in% names(peak_data)) vnum(peak_data[["Cell.Peak.Downstroke.90.40.Time"]]) else NA_real_
c_7525 <- if ("Cell.Peak.Downstroke.75.25.Time" %in% names(peak_data)) vnum(peak_data[["Cell.Peak.Downstroke.75.25.Time"]]) else NA_real_

C_9075_est  <- ifelse(is.finite(c_9040), (15/50) * c_9040, NA_real_)
C_9050_est  <- ifelse(is.finite(c_9040), (40/50) * c_9040, NA_real_)
C_9025_est1 <- ifelse(is.finite(c_9030) & is.finite(c_9020), c_9020 + 0.5 * (c_9030 - c_9020), NA_real_)
C_9025_est2 <- ifelse(is.finite(c_7525) & is.finite(C_9075_est), C_9075_est + c_7525, NA_real_)

peak_data[["Cell.Peak.Downstroke.90.75.Time"]] <- dplyr::coalesce(vnum(peak_data[["Cell.Peak.Downstroke.90.75.Time"]]), C_9075_est)
peak_data[["Cell.Peak.Downstroke.90.50.Time"]] <- dplyr::coalesce(vnum(peak_data[["Cell.Peak.Downstroke.90.50.Time"]]), C_9050_est)
peak_data[["Cell.Peak.Downstroke.90.25.Time"]] <- dplyr::coalesce(vnum(peak_data[["Cell.Peak.Downstroke.90.25.Time"]]), C_9025_est1, C_9025_est2)

# ===================================================================================================================
# 05) WHOLE-IMAGE TRACES: FIND & READ (_Whole Image Measurements_DataTable_CStats.csv)
# -------------------------------------------------------------------------------------------------------------------
whole_image_files <- list.files(
  path = input_dir,
  pattern = "_Whole Image Measurements_DataTable_CStats.csv$",
  full.names = TRUE,
  recursive = TRUE
)
whole_image_files <- whole_image_files[str_detect(whole_image_files, "CyteSeerResults[^/\\\\]*[/\\\\]")]

message("🔍 Found ", length(whole_image_files), " whole image trace file(s).")
print(whole_image_files)

if (length(whole_image_files) == 0) {
  stop("❌ No whole image data files found in expected CyteSeerResults folders.")
}

whole_image_data <- purrr::map_dfr(whole_image_files, function(file) {
  run_name <- basename(dirname(file))
  n_lines <- length(readLines(file, warn = FALSE))
  if (n_lines <= 405) {
    warning(paste("File skipped (too few lines):", file))
    return(NULL)
  }
  
  df <- tryCatch(
    readr::read_csv(
      file,
      skip = 405,
      na = c("", "NA", "N/A", "—", "-"),
      show_col_types = FALSE,
      col_types = readr::cols(
        Whole.Image.Peak.Has.Peaks_Sum = readr::col_double(),
        .default = readr::col_guess()
      )
    ),
    error = function(e) {
      warning(paste("Could not read:", file, "->", e$message))
      return(NULL)
    }
  )
  if (is.null(df)) return(NULL)
  
  if (!"Well" %in% names(df)) {
    df$Well <- tools::file_path_sans_ext(basename(file))
  }
  
  df$Run <- run_name
  df$SourceFile <- basename(file)
  df
})

names(whole_image_data) <- make.names(names(whole_image_data), unique = TRUE)

# ===================================================================================================================
# 06) WHOLE-IMAGE PEAKS: FIND, READ & AUGMENT
# -------------------------------------------------------------------------------------------------------------------
whole_image_peak_files <- list.files(
  path = input_dir,
  pattern = "_Whole-Image Peak Measurements_DataTable.csv$",
  full.names = TRUE,
  recursive = TRUE
)

message("🔍 Found ", length(whole_image_peak_files), " whole image peak data file(s).")
if (length(whole_image_peak_files) == 0) stop("❌ No whole image peak data files found.")

whole_image_peak_data <- purrr::map_dfr(whole_image_peak_files, function(file) {
  well_name <- basename(dirname(file))
  
  total_lines <- length(readLines(file, warn = FALSE))
  if (total_lines <= 82) {
    message("⚠️ Skipping short file: ", file)
    return(NULL)
  }
  
  df <- tryCatch(read.csv(file, skip = 82, fill = TRUE, header = TRUE),
                 error = function(e) { message("⚠️ Error reading file: ", file, " — ", e$message); NULL })
  
  if (is.null(df) || nrow(df) == 0) {
    message("⚠️ Skipping empty or invalid file: ", file)
    return(NULL)
  }
  
  df$Well <- well_name
  df
})

# compute missing spans (same logic)
series_wh <- intersect(c("Whole.Image.Peak.Smoothed.F.t.", "Whole.Image.Peak.F.t."),
                       names(whole_image_peak_data))
series_wh <- if (length(series_wh)) series_wh[1] else NA_character_

ts_candidates <- intersect(c("Timeslice", "Whole.Image.Peak.Timeslice"),
                           names(whole_image_peak_data))
ts_wh <- if (length(ts_candidates)) ts_candidates[1] else NA_character_

for (nm in c("Whole.Image.Peak.Upstroke.10.25.Time","Whole.Image.Peak.Upstroke.10.50.Time","Whole.Image.Peak.Upstroke.10.75.Time",
             "Whole.Image.Peak.Downstroke.90.75.Time","Whole.Image.Peak.Downstroke.90.50.Time","Whole.Image.Peak.Downstroke.90.25.Time")) {
  if (!nm %in% names(whole_image_peak_data)) whole_image_peak_data[[nm]] <- NA_real_
}

U1090_WI <- if ("Whole.Image.Peak.Upstroke.10.90.Time" %in% names(whole_image_peak_data))
  vnum(whole_image_peak_data[["Whole.Image.Peak.Upstroke.10.90.Time"]]) else NA_real_

upWI_nms  <- grep("^Whole\\.Image\\.Peak\\.Upstroke\\.(\\d{2})\\.90\\.Time$", names(whole_image_peak_data), value = TRUE)
if (length(upWI_nms)) {
  upWI_perc  <- as.numeric(sub("^.*Upstroke\\.(\\d{2})\\.90\\.Time$", "\\1", upWI_nms))
  upWI_spans <- setNames(lapply(upWI_nms, \(nm) vnum(whole_image_peak_data[[nm]])), as.character(upWI_perc))
  upWI_perc  <- sort(unique(upWI_perc))
  
  U2590 <- predict_span_vec(25, upWI_spans, upWI_perc)
  U5090 <- predict_span_vec(50, upWI_spans, upWI_perc)
  U7590 <- predict_span_vec(75, upWI_spans, upWI_perc)
  
  whole_image_peak_data[["Whole.Image.Peak.Upstroke.10.25.Time"]] <- ifelse(is.finite(U1090_WI), U1090_WI - U2590, NA_real_)
  whole_image_peak_data[["Whole.Image.Peak.Upstroke.10.50.Time"]] <- ifelse(is.finite(U1090_WI), U1090_WI - U5090, NA_real_)
  whole_image_peak_data[["Whole.Image.Peak.Upstroke.10.75.Time"]] <- ifelse(is.finite(U1090_WI), U1090_WI - U7590, NA_real_)
}

if (!is.na(series_wh) && !is.na(ts_wh)) {
  for (i in seq_len(nrow(whole_image_peak_data))) {
    y  <- parse_kic_series(whole_image_peak_data[[series_wh]][i]); if (length(y) < 5) next
    dt <- safe_dt_ms(whole_image_peak_data[[ts_wh]][i]);           if (!is.finite(dt) || dt <= 0) next
    base <- stats::quantile(y, 0.05, na.rm=TRUE); A <- max(y, na.rm=TRUE) - base
    if (!is.finite(A) || A <= 0) next
    yN <- (y - base) / A; t <- (seq_along(yN)-1) * dt; pk <- which.max(yN)
    
    t90 <- cross_down(yN, t, 0.90, pk)
    t75 <- cross_down(yN, t, 0.75, pk)
    t50 <- cross_down(yN, t, 0.50, pk)
    t25 <- cross_down(yN, t, 0.25, pk)
    
    if (is.finite(t90) && is.finite(t75)) whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.75.Time"]][i] <- t75 - t90
    if (is.finite(t90) && is.finite(t50)) whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.50.Time"]][i] <- t50 - t90
    if (is.finite(t90) && is.finite(t25)) whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.25.Time"]][i] <- t25 - t90
  }
}

w_9020 <- if ("Whole.Image.Peak.Downstroke.90.20.Time" %in% names(whole_image_peak_data)) vnum(whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.20.Time"]]) else NA_real_
w_9030 <- if ("Whole.Image.Peak.Downstroke.90.30.Time" %in% names(whole_image_peak_data)) vnum(whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.30.Time"]]) else NA_real_
w_9040 <- if ("Whole.Image.Peak.Downstroke.90.40.Time" %in% names(whole_image_peak_data)) vnum(whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.40.Time"]]) else NA_real_
w_7525 <- if ("Whole.Image.Peak.Downstroke.75.25.Time" %in% names(whole_image_peak_data)) vnum(whole_image_peak_data[["Whole.Image.Peak.Downstroke.75.25.Time"]]) else NA_real_

W_9075_est  <- ifelse(is.finite(w_9040), (15/50) * w_9040, NA_real_)
W_9050_est  <- ifelse(is.finite(w_9040), (40/50) * w_9040, NA_real_)
W_9025_est1 <- ifelse(is.finite(w_9030) & is.finite(w_9020), w_9020 + 0.5 * (w_9030 - w_9020), NA_real_)
W_9025_est2 <- ifelse(is.finite(w_7525) & is.finite(W_9075_est), W_9075_est + w_7525, NA_real_)

whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.75.Time"]] <- dplyr::coalesce(vnum(whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.75.Time"]]), W_9075_est)
whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.50.Time"]] <- dplyr::coalesce(vnum(whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.50.Time"]]), W_9050_est)
whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.25.Time"]] <- dplyr::coalesce(vnum(whole_image_peak_data[["Whole.Image.Peak.Downstroke.90.25.Time"]]), W_9025_est1, W_9025_est2)

# ===================================================================================================================
# 07) PLOT RAW CELL TRACES (+ export last 10s CSV)
# -------------------------------------------------------------------------------------------------------------------
raw_dir <- file.path(output_dir, "raw_traces")
csv_dir <- file.path(output_dir, "raw_csv")
cell_dir <- file.path(output_dir, "cell_peak_parameters")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

columns <- c("Well", "KIC.API.on.Cell.F.t._Average")
missing_cols <- setdiff(columns, names(trace_data))
if (length(missing_cols)) stop("❌ Missing required columns in trace_data: ", paste(missing_cols, collapse = ", "))

calcium <- trace_data[, columns, drop = FALSE]

make_time <- function(n, freq = NA_real_) {
  if (!is.finite(freq) || freq <= 0) freq <- 1
  (seq_len(n) - 1) / freq
}

get_last10s_df <- function(text_data, well_label = "NA") {
  if (is.na(text_data) || !nzchar(text_data)) return(NULL)
  
  data_str   <- gsub("\\{|\\}", "", text_data)
  data_parts <- unlist(strsplit(data_str, ";|\\|"), use.names = FALSE)
  
  raw_freq_val <- suppressWarnings(as.numeric(data_parts[4]))
  freq <- 10 / (raw_freq_val / 100)
  if (!is.finite(freq) || freq <= 0) freq <- 1
  
  values <- suppressWarnings(as.numeric(data_parts[11:length(data_parts)]))
  values <- values[is.finite(values)]
  npts <- length(values)
  if (npts == 0L) return(NULL)
  
  time_seq <- make_time(npts, freq)
  data_df  <- data.frame(Time = time_seq, Value = values)
  
  n_keep <- min(npts, max(1L, as.integer(round(freq * 10))))
  out <- tail(data_df, n = n_keep)
  out$Time <- out$Time - min(out$Time, na.rm = TRUE)
  out
}

for (n in seq_len(nrow(calcium))) {
  well <- calcium$Well[n]
  df10 <- get_last10s_df(calcium$KIC.API.on.Cell.F.t._Average[n], well_label = well)
  
  if (is.null(df10)) {
    warning(sprintf("Well %s: empty/invalid KIC string, skipping.", well))
    next
  }
  
  plot_obj <- ggplot(df10, aes(x = Time, y = Value)) +
    geom_line() +
    labs(
      title = paste0("Raw transients - CELL - Well ", well),
      x = "Time (seconds)",
      y = "Intensity (AU)"
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path(raw_dir, paste0("Raw_Trace_", well, ".jpg")),
    plot = plot_obj, width = 8, height = 6, units = "in", dpi = 120
  )
  
  write.csv(
    df10,
    file = file.path(csv_dir, paste0("Well_", well, "_raw.csv")),
    row.names = FALSE
  )
}

# ===================================================================================================================
# 08) CELL-PEAK TABLES: SAVE PER-PEAK & MEANS PER CELL
# -------------------------------------------------------------------------------------------------------------------
columns_to_select <- c(
  "Well",
  "Cell.Peak.Measurements.ID",
  "Cell.Peak.Cell.ID",
  "Cell.Peak.Peak.Time",
  "Cell.Peak.Amplitude" ,
  "Cell.Peak.Upstroke.10.100.Time",
  "Cell.Peak.Upstroke.10.90.Time",
  "Cell.Peak.Downstroke.100.10.Time",
  "Cell.Peak.Downstroke.100.30.Time",
  "Cell.Peak.Upstroke.Time",
  "Cell.Peak.Downstroke.Time",
  "Cell.Peak.Max.Upstroke.Velocity",
  "Cell.Peak.Max.Downstroke.Velocity",
  "Cell.Peak.Decay.Time",
  "Cell.Peak.Upstroke.10.25.Time","Cell.Peak.Upstroke.10.50.Time","Cell.Peak.Upstroke.10.75.Time",
  "Cell.Peak.Downstroke.90.75.Time","Cell.Peak.Downstroke.90.50.Time","Cell.Peak.Downstroke.90.25.Time",
  "Cell.Peak.Downstroke.90.10.Time",
  "Cell.Peak.Downstroke.90.20.Time","Cell.Peak.Downstroke.90.30.Time","Cell.Peak.Downstroke.90.40.Time",
  "Cell.Peak.Downstroke.75.25.Time"
)

columns_to_select <- intersect(columns_to_select, names(peak_data))
all_peak_data <- peak_data[, columns_to_select, drop = FALSE]

write.csv(all_peak_data,
          file = file.path(cell_dir, paste0("Cell_Peak_Data_", prefix, ".csv")),
          row.names = FALSE)

metric_cols_cell <- setdiff(names(all_peak_data), c("Well","Cell.Peak.Measurements.ID","Cell.Peak.Cell.ID"))

mean_cellpeak_data <- all_peak_data %>%
  group_by(Cell.Peak.Cell.ID, Well) %>%
  summarise(
    Num.Peaks = n(),
    across(.cols = all_of(metric_cols_cell), .fns = mean, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(mean_cellpeak_data,
          file = file.path(cell_dir, paste0("Mean_Peak_Cell_Data_", prefix, ".csv")),
          row.names = FALSE)

# ===================================================================================================================
# 09) PLOT WHOLE-IMAGE TRACES (PER-WELL) + CSV EXPORT  [FIXED]
#     FIX vs your original: compute per-well window; DO NOT reuse 'final_10_seconds' from the CELL loop.
# -------------------------------------------------------------------------------------------------------------------
whole_image_dir <- file.path(output_dir, "whole_image_traces")
whole_image_csv_dir <- file.path(output_dir, "whole_image_csv")
whole_image_parameters <- file.path(output_dir, "whole_image_parameters")
dir.create(whole_image_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(whole_image_csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(whole_image_parameters, recursive = TRUE, showWarnings = FALSE)

# Whole-image peak table export (same as your intent)
columns_to_select3 <- c(
  "Well",
  "Whole.Image.Peak.Measurements.ID",
  "Whole.Image.Peak.Peak.Time",
  "Whole.Image.Peak.Amplitude",
  "Whole.Image.Peak.Upstroke.10.100.Time",
  "Whole.Image.Peak.Upstroke.10.90.Time",
  "Whole.Image.Peak.Downstroke.100.10.Time",
  "Whole.Image.Peak.Max.Upstroke.Velocity",
  "Whole.Image.Peak.Max.Downstroke.Velocity",
  "Whole.Image.Peak.Decay.Time",
  "Whole.Image.Peak.Max.Value",
  "Whole.Image.Peak.Upstroke.10.25.Time","Whole.Image.Peak.Upstroke.10.50.Time","Whole.Image.Peak.Upstroke.10.75.Time",
  "Whole.Image.Peak.Downstroke.90.75.Time","Whole.Image.Peak.Downstroke.90.50.Time","Whole.Image.Peak.Downstroke.90.25.Time",
  "Whole.Image.Peak.Downstroke.90.10.Time",
  "Whole.Image.Peak.Downstroke.90.20.Time","Whole.Image.Peak.Downstroke.90.30.Time","Whole.Image.Peak.Downstroke.90.40.Time",
  "Whole.Image.Peak.Downstroke.75.25.Time"
)
columns_to_select3 <- intersect(columns_to_select3, names(whole_image_peak_data))
whole_image_peaks <- whole_image_peak_data[, columns_to_select3, drop = FALSE]

write.csv(whole_image_peaks,
          file = file.path(whole_image_parameters, paste0("Whole_image_Peak_Data_", prefix, ".csv")),
          row.names = FALSE)

# Plot whole-image traces + export last 10s
kic_wh_col <- "KIC.API.on.Whole.Image.F.t._Average"
if (!kic_wh_col %in% names(whole_image_data)) {
  stop("❌ Missing required column in whole_image_data: ", kic_wh_col)
}

for (n in seq_len(nrow(whole_image_data))) {
  well <- whole_image_data$Well[n]
  df10 <- get_last10s_df(whole_image_data[[kic_wh_col]][n], well_label = well)
  
  if (is.null(df10)) {
    warning(sprintf("Well %s: empty/invalid WHOLE-IMAGE KIC string, skipping.", well))
    next
  }
  
  plot_obj <- ggplot(df10, aes(x = Time, y = Value)) +
    geom_line() +
    labs(
      title = paste0("Whole image transient - Well ", well),
      x = "Time (seconds)",
      y = "Intensity (AU)"
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path(whole_image_dir, paste0("Whole_Image_Trace_", well, ".jpg")),
    plot = plot_obj,
    width = 8, height = 6, units = "in", dpi = 120
  )
  
  write.csv(
    df10,
    file = file.path(whole_image_csv_dir, paste0("Well_", well, "_whole_image_raw_last10s.csv")),
    row.names = FALSE
  )
}

# ===================================================================================================================
# 10) WHOLE-IMAGE SUMMARIES (FULL & MEANS-ONLY)
# -------------------------------------------------------------------------------------------------------------------
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

group_keys <- intersect(c("Run", "Well"), names(whole_image_peaks))
metrics_whole <- setdiff(names(whole_image_peaks), c(group_keys, "Whole.Image.Peak.Measurements.ID"))

whole_image_summary_full <- whole_image_peaks %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) %>%
  dplyr::summarise(
    Num.Peaks = dplyr::n(),
    dplyr::across(
      .cols = dplyr::all_of(metrics_whole),
      .fns  = list(
        mean   = ~ mean(.x, na.rm = TRUE),
        sd     = ~ sd(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE),
        iqr    = ~ IQR(.x, na.rm = TRUE)
      ),
      .names = "{.col}.{.fn}"
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::across(dplyr::all_of(group_keys)))

readr::write_csv(
  whole_image_summary_full,
  file.path(summary_dir, paste0("Whole_Image_Summary_FULL_", prefix, ".csv"))
)

whole_image_summary_means <- whole_image_summary_full %>%
  dplyr::select(dplyr::all_of(group_keys), Num.Peaks, dplyr::ends_with(".mean"))

readr::write_csv(
  whole_image_summary_means,
  file.path(summary_dir, paste0("Whole_Image_Summary_MEANS_", prefix, ".csv"))
)

# Done message (optional)
if (requireNamespace("tcltk", quietly = TRUE)) {
  tryCatch(
    tcltk::tkmessageBox(
      title = "Done ✅",
      message = paste0("Analysis completed.\n\nOutputs written to:\n", output_dir),
      icon = "info",
      type = "ok"
    ),
    error = function(e) NULL
  )
}

message("✅ All done. Output folder: ", output_dir)