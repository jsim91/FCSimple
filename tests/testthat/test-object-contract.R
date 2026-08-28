# Core object-contract tests. These exercise the documented FCSimple analysis
# pipeline on the bundled example FCS files. They are skipped when the heavy
# Bioconductor dependencies (flowCore, flowWorkspace, FlowSOM) are not
# installed, so they run in a standard check environment but degrade gracefully
# in a minimal one.

test_that("fcs_join returns the documented object contract", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")

  files <- fcs_example_files()
  expect_gt(length(files), 0)
  expect_true(all(file.exists(files)))

  obj <- fcs_join(files = files, apply_transform = TRUE)

  expect_type(obj, "list")
  expected <- c("data", "raw", "source", "run_date", "metadata",
                "collection_instrument", "object_history")
  for (nm in expected) {
    expect_true(nm %in% names(obj), info = nm)
  }

  expect_equal(nrow(obj$data), length(obj$source))
  expect_equal(ncol(obj$data), ncol(obj$raw))
  expect_equal(nrow(obj$metadata), length(unique(obj$source)))
})

test_that("pipeline steps augment the object rather than replacing it", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("igraph")

  obj <- fcs_join(files = fcs_example_files(), downsample_size = 500)

  audited <- fcs_audit(obj)
  expect_true("object_history" %in% names(audited))

  clustered <- fcs_cluster(audited, algorithm = "leiden")
  expect_true("leiden" %in% names(clustered))
  expect_equal(length(clustered$leiden$clusters), nrow(audited$data))
})

test_that("the internal temp-dir helper is writable and unique", {
  d1 <- FCSimple:::.fcs_temp_dir()
  d2 <- FCSimple:::.fcs_temp_dir()
  on.exit(unlink(c(d1, d2), recursive = TRUE, force = TRUE), add = TRUE)

  expect_true(dir.exists(d1))
  expect_true(dir.exists(d2))
  expect_false(identical(normalizePath(d1), normalizePath(d2)))

  probe <- file.path(d1, "probe.txt")
  writeLines("ok", probe)
  expect_true(file.exists(probe))
})
