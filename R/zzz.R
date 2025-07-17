# load Roboto font for use in plotting
.onLoad <- function(...) {
  fs::path_package(
    "roboto",
    package = "hptmUsage"
  ) |>
    paste0("/*.ttf") |>
    Sys.glob() |>
    systemfonts::add_fonts()
}
