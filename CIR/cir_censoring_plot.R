################################################################################
################ CIR censoring-model comparison boxplots #######################
################################################################################

## This script compares three DR settings, all with the SOP model misspecified:
##   1. PS correct, censoring correct
##   2. PS correct, censoring wrong
##   3. PS wrong,   censoring correct
##


## Parse command-line options supplied as --name=value.
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) return(list())
  key_value <- strsplit(sub("^--", "", args), "=", fixed = TRUE)
  setNames(lapply(key_value, `[`, 2L), vapply(key_value, `[`, "", 1L))
}

get_arg <- function(args, name, default) {
  value <- args[[name]]
  if (is.null(value) || is.na(value) || !nzchar(value)) default else value
}

## Locate this script so the default Results path is reproducible.
infer_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0L) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
}

## Resolve the input and output paths.
args <- parse_args()
script_dir <- infer_script_dir()
result_dir <- normalizePath(
  get_arg(args, "result_dir", file.path(script_dir, "Results")),
  mustWork = TRUE
)
output_file <- get_arg(
  args,
  "output",
  file.path(result_dir, "cir_censoring_boxplots.pdf")
)
## Select the three DR robustness settings compared in the figure.
settings <- c(
  "sop_wrong_ps_correct_cen_correct",
  "sop_wrong_ps_correct_cen_wrong",
  "sop_wrong_ps_wrong_cen_correct"
)

setting_labels <- c(
  "PS correct\nCEN correct",
  "PS correct\nCEN wrong",
  "PS wrong\nCEN correct"
)

## Muted, colorblind-friendly palette for publication figures.
setting_colors <- c("#4E79A7", "#F28E2B", "#59A14F")
setting_border_colors <- c("#2F506F", "#A85B16", "#356A35")

## Find replicate-level simulation results across all run directories.
raw_files <- list.files(
  result_dir,
  pattern = "^raw_results\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(raw_files) == 0L) {
  stop("No raw_results.csv files found under: ", result_dir, call. = FALSE)
}

## Read one run, recover its design information, and retain figure inputs.
read_run <- function(path) {
  run_name <- basename(dirname(path))
  # Extract the example and sample size from the folder name.
  match_info <- regexec(
    "^cir_example([0-9]+)_n([0-9]+)$",
    run_name
  )
  fields <- regmatches(run_name, match_info)[[1L]]
  if (length(fields) != 3L) return(NULL)

  # Validate the input schema before filtering the simulation results.
  dat <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("method", "sim", "setting", "parameter", "truth", "bias")
  missing_columns <- setdiff(required, names(dat))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing column(s) in ", path, ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  # Keep DR results for rho0 and rho1 under the selected settings.
  dat <- dat[
    dat$method == "DR" &
      dat$setting %in% settings &
      dat$parameter %in% c("rho0", "rho1"),
    ,
    drop = FALSE
  ]
  dat$example <- as.integer(fields[[2L]])
  dat$n <- as.integer(fields[[3L]])
  dat
}

## Combine all recognized runs into one plotting data set.
plot_data <- do.call(rbind, Filter(Negate(is.null), lapply(raw_files, read_run)))

if (is.null(plot_data) || nrow(plot_data) == 0L) {
  stop("No matching DR results were found.", call. = FALSE)
}

## Apply readable labels and a fixed setting order across panels.
plot_data$setting <- factor(
  plot_data$setting,
  levels = settings,
  labels = setting_labels
)

## Create one panel for each example-by-sample-size combination.
panel_info <- unique(plot_data[c("example", "n")])
panel_info <- panel_info[order(panel_info$example, panel_info$n), , drop = FALSE]
parameters <- c("rho0", "rho1")

missing_parameters <- setdiff(parameters, unique(plot_data$parameter))
if (length(missing_parameters) > 0L) {
  stop(
    "Missing parameter(s): ", paste(missing_parameters, collapse = ", "),
    call. = FALSE
  )
}

## Open the PDF graphics device for the completed multi-panel figure.
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

grDevices::pdf(output_file, width = 13, height = 12, onefile = TRUE)
on.exit(grDevices::dev.off(), add = TRUE)

old_par <- graphics::par(no.readonly = TRUE)
on.exit(graphics::par(old_par), add = TRUE)

## Use a common y-axis range so panels can be compared directly.
y_values <- plot_data$bias
y_range <- range(y_values, finite = TRUE)
padding <- max(diff(y_range) * 0.06, 0.02)
y_limit <- y_range + c(-padding, padding)

## Define the panel layout and grouped positions for rho0 and rho1.
n_panels <- nrow(panel_info)
n_columns <- min(2L, n_panels)
n_rows <- ceiling(n_panels / n_columns)
box_positions <- c(1.15, 2.4, 3.65, 5.95, 7.2, 8.45)
group_centers <- c(mean(box_positions[1:3]), mean(box_positions[4:6]))
group_divider <- mean(box_positions[3:4])
x_limit <- c(0.25, 9.35)
short_labels <- rep(c("PS + CEN\ncorrect", "Only PS\ncorrect", "Only CEN\ncorrect"), 2L)

graphics::par(
  mfrow = c(n_rows, n_columns),
  mar = c(5.7, 5.0, 4.0, 1.0),
  oma = c(1.7, 0.5, 0.5, 0.5),
  las = 1
)

## Draw one panel for each simulation design.
for (panel_index in seq_len(n_panels)) {
  example_value <- panel_info$example[[panel_index]]
  n_value <- panel_info$n[[panel_index]]
  panel_data <- plot_data[
    plot_data$example == example_value & plot_data$n == n_value,
    ,
    drop = FALSE
  ]

  # Assemble three PS/censoring settings for each target parameter.
  box_values <- vector("list", length(box_positions))
  box_index <- 1L
  for (parameter in parameters) {
    for (setting_label in setting_labels) {
      box_values[[box_index]] <- panel_data[
        panel_data$parameter == parameter &
          as.character(panel_data$setting) == setting_label,
        "bias"
      ]
      box_index <- box_index + 1L
    }
  }

  # Draw three boxes for rho0 followed by three boxes for rho1.
  graphics::boxplot(
    box_values,
    at = box_positions,
    names = short_labels,
    col = rep(setting_colors, 2L),
    border = rep(setting_border_colors, 2L),
    boxwex = 0.78,
    outline = TRUE,
    xaxt = "n",
    xlim = x_limit,
    ylim = y_limit,
    ylab = "Bias",
    xlab = "",
    main = sprintf("Example %d, n = %d", example_value, n_value),
    cex.axis = 0.96,
    cex.lab = 1.18,
    cex.main = 1.25
  )

  # Add compact setting labels below the six boxes.
  original_mgp <- graphics::par("mgp")
  graphics::par(mgp = c(original_mgp[[1L]], 1.50, original_mgp[[3L]]))
  graphics::axis(
    side = 1,
    at = box_positions,
    labels = short_labels,
    cex.axis = 0.96,
    gap.axis = -1
  )
  graphics::par(mgp = original_mgp)

  # Mark zero and visually separate the two parameter groups.
  graphics::abline(h = 0, lty = 2, lwd = 1.1, col = "#6B7280")
  graphics::abline(v = group_divider, lty = 3, lwd = 1.0, col = "#CBD5E1")

  # Overlay Monte Carlo means and their numerical values.
  group_means <- vapply(box_values, mean, numeric(1L), na.rm = TRUE)
  mean_labels <- sprintf("%.2f", group_means)
  graphics::points(
    box_positions,
    group_means,
    pch = 18,
    cex = 1.20,
    col = "#111827"
  )

  # Keep each mean label inside the shared plotting range.
  plot_span <- diff(y_limit)
  label_x <- box_positions
  label_y <- group_means + 0.018 * plot_span
  upper_guard <- y_limit[[2L]] - 0.04 * plot_span
  lower_guard <- y_limit[[1L]] + 0.04 * plot_span
  label_y[label_y > upper_guard] <- group_means[label_y > upper_guard] - 0.080 * plot_span
  label_y <- pmax(lower_guard, pmin(upper_guard, label_y))

  label_cex <- 0.72
  graphics::text(
    label_x,
    label_y,
    labels = mean_labels,
    cex = label_cex,
    font = 2,
    col = "#1F2937"
  )
  # Label the rho0 and rho1 groups above the panel.
  graphics::mtext(
    c(expression(rho[0]), expression(rho[1])),
    side = 3,
    line = 0.25,
    at = group_centers,
    cex = 1.18,
    font = 2
  )
}

## Close the PDF and report its saved location.
grDevices::dev.off()
on.exit(NULL, add = FALSE)

message("Saved boxplots to: ", normalizePath(output_file, mustWork = TRUE))
