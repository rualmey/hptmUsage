# helper functions -------------------------------------------------------------

## .paramcheck -----------------------------------------------------------------

test_that(".paramcheck passes with valid arguments", {
  # Create a dummy function to check against
  dummy_fun <- function(a, b, c) {}
  valid_params <- list(a = 1, b = 2)

  expect_silent(.paramcheck(valid_params, dummy_fun))
})

test_that(".paramcheck throws error for invalid argument names", {
  dummy_fun <- function(a, b) {}
  invalid_params <- list(a = 1, z = 3) # 'z' is not in dummy_fun

  expect_error(
    .paramcheck(invalid_params, dummy_fun),
    "Invalid argument in 'invalid_params': z"
  )
})

test_that(".paramcheck throws error for excluded arguments", {
  dummy_fun <- function(object, i, other) {}
  # 'object' and 'i' are typically excluded in wrapper functions
  params <- list(other = 1, object = 2)

  expect_error(
    .paramcheck(params, dummy_fun, exclusions = c("object", "i")),
    "Argument not allowed: object"
  )
})

test_that(".paramcheck handles multiple errors gracefully", {
  dummy_fun <- function(a, b) {}
  params <- list(x = 1, y = 2)

  expect_error(
    .paramcheck(params, dummy_fun),
    "Invalid arguments in 'params': x, y"
  )
})

# Main Function ----------------------------------------------------------------

## generateReport --------------------------------------------------------------

test_that("generateReport runs successfully (Happy Path)", {
  # Setup test environment
  output_dir <- withr::local_tempdir("output_dir")
  resource_dir <- withr::local_tempdir("resource_dir")

  # Create a dummy template file so fs::file_exists passes
  fs::file_create(fs::path(resource_dir, "report_template.qmd"))

  # Mock Data
  mock_dataset <- "Simulated QFeatures Object"

  # Mocks
  # 1. fs::path_package: Point to our temp resource dir instead of installed package
  mock_path_package <- function(...) return(resource_dir)

  # 2. quarto::quarto_available: Always TRUE
  mock_quarto_avail <- function(...) TRUE

  # 3. quarto::quarto_render:
  #    This mock must simulate the side effects of rendering:
  #    a) Creating the HTML output file
  #    b) Creating the processed dataset RDS
  mock_quarto_render <- function(input, output_file, execute_params, ...) {
    # Simulate HTML creation (input is the qmd, output_file is the target html name)
    # The function runs in a temp dir, so 'input' path contains that dir
    work_dir <- fs::path_dir(input)
    fs::file_create(fs::path(work_dir, output_file))

    # Simulate saving the processed dataset
    saveRDS(
      paste0(readRDS(execute_params$ds_path), " - Processed"),
      fs::path(work_dir, "ds_processed.rds")
    )
    return(TRUE)
  }

  # Apply mocks using mockery
  # We stub functions called *inside* generateReport
  mockery::stub(generateReport, "fs::path_package", mock_path_package)
  mockery::stub(generateReport, "quarto::quarto_available", mock_quarto_avail)
  mockery::stub(generateReport, "quarto::quarto_render", mock_quarto_render)

  # Execution
  res <- suppressMessages(generateReport(
    dataset = mock_dataset,
    output_dir = output_dir
  ))

  # Assertions
  # 1. Return value is the "processed" dataset
  expect_equal(res, "Simulated QFeatures Object - Processed")

  # 2. Output file exists in the user-specified output directory
  files <- fs::dir_ls(output_dir, glob = "*.html")
  expect_length(files, 1)
  expect_true(grepl("_hptmUsage.html$", files))
})

test_that("generateReport validates specific *_params arguments", {
  # We don't need to mock quarto here because it should fail before reaching it
  # provided .paramcheck works.

  # Invalid histone_params (hptmUsage::histonesFromUniprot doesn't have 'bad_arg')
  expect_error(
    generateReport("data", "out", histone_params = list(bad_arg = 1)),
    "Invalid argument in 'histone_params': bad_arg"
  )

  # Invalid msa_params (hptmUsage::alignHistones exclusion 'unaligned_histones')
  expect_error(
    generateReport("data", "out", msa_params = list(unaligned_histones = 1)),
    "Argument not allowed: unaligned_histones"
  )
})

test_that("generateReport fails if Quarto is not available", {
  # Mock quarto_available to throw error (as intended by min/max/error=TRUE args)
  mock_quarto_avail <- function(...) stop("Quarto not found")

  mockery::stub(generateReport, "quarto::quarto_available", mock_quarto_avail)

  expect_error(
    generateReport("data", "out"),
    "Quarto not found"
  )
})

test_that("generateReport fails if template is missing", {
  # Setup empty resource dir (no .qmd)
  resource_dir <- withr::local_tempdir("resource_dir")
  mock_path_package <- function(...) return(resource_dir)

  mockery::stub(generateReport, "fs::path_package", mock_path_package)
  mockery::stub(generateReport, "quarto::quarto_available", function(...) TRUE)

  expect_error(
    generateReport("data", "out"),
    "Template not found in package resources"
  )
})

test_that("generateReport handles Quarto rendering errors", {
  resource_dir <- withr::local_tempdir("resource_dir")
  fs::file_create(fs::path(resource_dir, "report_template.qmd"))

  mock_path_package <- function(...) return(resource_dir)
  mock_quarto_render <- function(...) stop("Pandoc error")

  mockery::stub(generateReport, "fs::path_package", mock_path_package)
  mockery::stub(generateReport, "quarto::quarto_available", function(...) TRUE)
  mockery::stub(generateReport, "quarto::quarto_render", mock_quarto_render)

  expect_error(
    suppressMessages(generateReport("data", "out")),
    "Quarto rendering failed: Pandoc error"
  )
})

test_that("generateReport prevents overwriting existing files unless specified", {
  output_dir <- withr::local_tempdir("output_dir")
  resource_dir <- withr::local_tempdir("resource_dir")
  fs::file_create(fs::path(resource_dir, "report_template.qmd"))

  # We need to predict the filename to create a collision
  # Since filename involves Sys.time(), we can mock Sys.time OR
  # we can just ensure the mock render produces a file that we *already* put in output

  fixed_time <- as.POSIXct("2025-01-01 12:00:00")
  expected_filename <- "250101_1200_hptmUsage.html"

  # Create the collision file
  fs::file_create(fs::path(output_dir, expected_filename))

  mock_path_package <- function(...) return(resource_dir)

  # Mock render: creates the file in temp, which will conflict when moved
  mock_quarto_render <- function(input, output_file, ...) {
    work_dir <- fs::path_dir(input)
    fs::file_create(fs::path(work_dir, output_file)) # Create file in temp
    # Also need processed data for function to finish before check
    saveRDS("data", fs::path(work_dir, "ds_processed.rds"))
    return(TRUE)
  }

  mockery::stub(generateReport, "fs::path_package", mock_path_package)
  mockery::stub(generateReport, "quarto::quarto_available", function(...) TRUE)
  mockery::stub(generateReport, "quarto::quarto_render", mock_quarto_render)
  # Mock Sys.time so the filename matches our collision file
  mockery::stub(generateReport, "Sys.time", fixed_time)

  # Case 1: overwrite = FALSE (Default)
  expect_error(
    suppressMessages(generateReport("data", output_dir)),
    "Output file already exists"
  )

  # Case 2: overwrite = TRUE
  expect_silent(
    suppressMessages(generateReport("data", output_dir, overwrite = TRUE))
  )
})

test_that("generateReport validates simple arguments", {
  expect_error(
    generateReport("data", "out", generate_usageplots = "invalid_option"),
    "'arg' should be one of"
  )
})
