## Run with this command to view the extra info for failing tests: testthat::test_file("tests/testthat/test-cran.R", reporter = "progress") 

test_that("bgms functions satisfy user-visible return contract", {
  set.seed(123)
  
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:5]
  p   <- ncol(dat)
  
  tests <- list(
    list(
      label      = "single_bgm",
      fun_label  = "bgm",
      fun        = bgms::bgm,
      args       = list(
        x               = dat,
        iter            = 1000,
        warmup          = 1000,
        chains          = 2,
        edge_selection  = TRUE,
        edge_prior      = "Bernoulli",
        na_action       = "listwise",
        update_method   = "adaptive-metropolis",
        display_progress = "none"
      ),
      expected_fields = c(
        "posterior_summary_main",
        "posterior_summary_pairwise",
        "posterior_summary_indicator",
        "posterior_mean_main",
        "posterior_mean_pairwise",
        "posterior_mean_indicator",
        "raw_samples",
        "arguments"
      ),
      matrix_fields = c(
        "posterior_mean_pairwise",
        "posterior_mean_indicator"
      )
    ),
    list(
      label      = "compare_bgm",
      fun_label  = "bgmCompare",
      fun        = bgms::bgmCompare,
      args       = list(
        x                   = dat,
        group_indicator     = rep(1:2, each = 20),
        iter                = 1000,
        warmup              = 1000,
        chains              = 2,
        difference_selection = FALSE,
        na_action           = "listwise",
        update_method       = "adaptive-metropolis",
        display_progress    = "none"
      ),
      expected_fields = c(
        "posterior_summary_main_baseline",
        "posterior_summary_pairwise_baseline",
        "posterior_summary_main_differences",
        "posterior_summary_pairwise_differences",
        "posterior_mean_main_baseline",
        "posterior_mean_pairwise_baseline",
        "raw_samples",
        "arguments"
      ),
      matrix_fields = c(
        "posterior_mean_pairwise_baseline"
      )
    )
  )
  
  for (spec in tests) {
    
    result <- do.call(spec$fun, spec$args)
    
    missing <- setdiff(spec$expected_fields, names(result))
    expect_equal(
      length(missing), 0,
      info = sprintf(
        "[%s / %s] Missing fields: %s",
        spec$fun_label,
        spec$label,
        if (length(missing) == 0) "<none>" else paste(missing, collapse = ", ")
      )
    )
    
    for (fld in spec$matrix_fields) {
      if (!fld %in% names(result)) next
      
      actual_dim <- if (!is.null(result[[fld]])) {
        paste(dim(result[[fld]]), collapse = "x")
      } else {
        "NULL"
      }
      
      expect_equal(
        dim(result[[fld]]),
        c(p, p),
        info = sprintf(
          "[%s / %s / %s] Wrong dim: expected %ix%i, got %s",
          spec$fun_label,
          spec$label,
          fld,
          p, p,
          actual_dim
        )
      )
      
      expect_false(
        all(is.na(result[[fld]])),
        info = sprintf(
          "[%s / %s / %s] All entries are NA",
          spec$fun_label,
          spec$label,
          fld
        )
      )
    }
  }
})
