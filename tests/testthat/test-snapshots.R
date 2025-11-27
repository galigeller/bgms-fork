test_that("bgms user-visible structure matches snapshot", {
  skip_on_cran()
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
        x                = dat,
        iter             = 1000,
        warmup           = 1000,
        chains           = 2,
        edge_selection   = TRUE,
        edge_prior       = "Bernoulli",
        na_action        = "listwise",
        update_method    = "adaptive-metropolis",
        display_progress = "none"
      ),
      matrix_fields = c("posterior_mean_pairwise", "posterior_mean_indicator")
    ),
    list(
      label      = "compare_bgm",
      fun_label  = "bgmCompare",
      fun        = bgms::bgmCompare,
      args       = list(
        x                    = dat,
        group_indicator      = rep(1:2, each = 20),
        iter                 = 1000,
        warmup               = 1000,
        chains               = 2,
        difference_selection = FALSE,
        na_action            = "listwise",
        update_method        = "adaptive-metropolis",
        display_progress     = "none"
      ),
      matrix_fields = c("posterior_mean_pairwise_baseline")
    )
  )
  
  out <- list()
  
  for (spec in tests) {
    result <- do.call(spec$fun, spec$args)
    
    ## 🔥 THIS preserves the ORIGINAL snapshot structure you pasted
    out[[spec$label]] <- list(
      fun_label    = spec$fun_label,
      result_class = unclass(class(result)),
      result_names = names(result),
      matrix_meta  = lapply(spec$matrix_fields, function(fld) {
        if (!fld %in% names(result)) return(NULL)
        obj <- result[[fld]]
        list(
          dim      = dim(obj),
          rownames = rownames(obj),
          colnames = colnames(obj)
        )
      }),
      ## 🔥 NEW — adds the FULL raw bgms object without changing the structure
      full_result_dput = paste(capture.output(dput(result)), collapse = "\n")
      
    )
  }
  
  ## 🔥 EXACT SAME snapshot engine as before
  expect_snapshot_value(out, style = "json2")
})
