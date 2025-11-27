data("Wenchuan", package = "bgms")
dat <- na.omit(Wenchuan)[1:40, 1:5]

tests <- list(
  list(
    label = "single_bgm",
    fun   = bgms::bgm,
    args  = list(
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
    label = "compare_bgm",
    fun   = bgms::bgmCompare,
    args  = list(
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

for (spec in tests) {
  local_spec <- spec  # capture for the closure
  
  test_that(paste0("bgms snapshot: ", local_spec$label), {
    skip_on_cran()
    set.seed(123)
    
    result <- do.call(local_spec$fun, local_spec$args)
    
    obj <- list(
      fun_label        = local_spec$label,
      result_class     = unclass(class(result)),
      result_names     = names(result),
      full_result_dput = paste(capture.output(dput(result)), collapse = "\n")
    )
    
    expect_snapshot_value(obj, style = "json2")
  })
}
