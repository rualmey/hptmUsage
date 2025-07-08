# always attach devtools
if (interactive()) {
  suppressMessages(require(devtools))
}

# warn on partial matches
options(
  warnPartialMatchAttr = TRUE,
  warnPartialMatchDollar = TRUE,
  warnPartialMatchArgs = TRUE
)

# enable autocompletions for package names in
# `require()`, `library()`
utils::rc.settings(ipck = TRUE)

# fancy quotes are annoying and lead to
# 'copy + paste' bugs / frustrations
options(useFancyQuotes = FALSE)

# print libPaths on startup
if (length(.libPaths()) > 1) {
  msg <- "Using libraries at paths:\n"
} else {
  msg <- "Using library at path:\n"
}
libs <- paste("-", .libPaths(), collapse = "\n")
message(msg, libs, sep = "")

# stop flooding the terminal
options("max.print" = 100)
