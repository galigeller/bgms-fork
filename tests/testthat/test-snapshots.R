# test this file: testthat::test_file("tests/testthat/test-snapshots.R", reporter = "progress")

test_that("bgms snapshots match expected outputs", {
  skip_on_cran()
  set.seed(123)
  
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:5]
  
  tests <- list(
    list(
      label     = "single_bgm",
      fun_label = "bgm",
      fun       = bgms::bgm,
      args      = list(
        x                = dat,
        iter             = 1000,
        warmup           = 1000,
        chains           = 2,
        edge_selection   = TRUE,
        edge_prior       = "Bernoulli",
        na_action        = "listwise",
        update_method    = "adaptive-metropolis",
        display_progress = "none"
      )
    ),
    list(
      label     = "compare_bgm",
      fun_label = "bgmCompare",
      fun       = bgms::bgmCompare,
      args      = list(
        x                    = dat,
        group_indicator      = rep(1:2, each = 20),
        iter                 = 1000,
        warmup               = 1000,
        chains               = 2,
        difference_selection = FALSE,
        na_action            = "listwise",
        update_method        = "adaptive-metropolis",
        display_progress     = "none"
      )
    )
  )
  
  out <- list()
  
  for (spec in tests) {
    result <- do.call(spec$fun, spec$args)
    
    out[[spec$label]] <- list(
      fun_label        = spec$fun_label,
      result_class     = unclass(class(result)),
      result_names     = names(result),
      # full bgms object encoded as text so json2 can handle it
      full_result_dput = paste(capture.output(dput(result)), collapse = "\n")
    )
  }
  
  expect_snapshot_value(out, style = "json2")
})
