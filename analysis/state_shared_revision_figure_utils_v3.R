suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
})

figure_font <- "Arial"
figure_colours <- c(
  ink = "#24292E", muted = "#6B747C", line = "#C7CDD2",
  normal = "#3274A1", adenoma = "#D45D3F", gold = "#D6A23A",
  purple = "#7669A8", green = "#3F8A70", grey = "#7B858D",
  pale_blue = "#EDF4F8", pale_orange = "#FBF0EC",
  pale_gold = "#FBF6E9", pale_purple = "#F2F0F7", pale = "#F5F6F7"
)

read_tsv <- function(path) {
  connection <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection), add = TRUE)
  read.delim(connection, check.names = FALSE, quote = "", comment.char = "")
}

write_source <- function(data, directory, name) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  write.table(
    data, file.path(directory, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  )
}

theme_cb <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = figure_font) +
    theme(
      text = element_text(family = figure_font, colour = figure_colours[["ink"]]),
      axis.line = element_line(linewidth = 0.35, colour = figure_colours[["ink"]]),
      axis.ticks = element_line(linewidth = 0.3, colour = figure_colours[["ink"]]),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.6),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(2, 2, 2, 2, "mm"),
      plot.tag = element_text(size = 9, face = "bold"),
      plot.tag.location = "panel",
      plot.tag.position = c(0.012, 0.992)
    )
}

tagged <- function(plot) {
  plot & theme(
    plot.tag = element_text(family = figure_font, size = 9, face = "bold"),
    plot.tag.location = "panel",
    plot.tag.position = c(0.012, 0.992)
  )
}

export_cb_figure <- function(
    plot, output_directory, stem, width_mm = 178, height_mm = 195) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(output_directory, paste0(stem, ".pdf")), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = cairo_pdf, bg = "white"
  )
  ggsave(
    file.path(output_directory, paste0(stem, ".svg")), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = function(...) svglite::svglite(..., fix_text_size = FALSE),
    bg = "white"
  )
  ggsave(
    file.path(output_directory, paste0(stem, ".tiff")), plot,
    width = width_mm, height = height_mm, units = "mm", dpi = 600,
    device = ragg::agg_tiff, compression = "lzw", background = "white"
  )
  ggsave(
    file.path(output_directory, paste0(stem, ".png")), plot,
    width = width_mm, height = height_mm, units = "mm", dpi = 300,
    device = ragg::agg_png, background = "white"
  )
}
