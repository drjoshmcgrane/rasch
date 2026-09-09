factor_item_selector_rasch_data <- function() {
  set.seed(9062)
  X <- matrix(rbinom(240L, 1L, 0.5), 60L, 4L,
              dimnames = list(NULL, paste0("I0", 1:4)))
  list(X = X, data = as.data.frame(X))
}

test_that("ordinary external factors can share item names when items are explicit", {
  made <- factor_item_selector_rasch_data()
  X <- made$X
  d <- made$data
  different <- data.frame(I01 = factor(rep(c("A", "B"), length.out = nrow(X))))

  expect_error(
    rasch(d, factors = different),
    "share item-data names but contain different values"
  )

  by_name <- rasch(d, factors = different, items = colnames(X),
                   n_groups = 4, maxit = 60, tol = 1e-6)
  by_index <- rasch(d, factors = different,
                    items = match(colnames(X), names(d)),
                    n_groups = 4, maxit = 60, tol = 1e-6)
  from_matrix <- rasch(X, factors = different, items = colnames(X),
                       n_groups = 4, maxit = 60, tol = 1e-6)

  expect_equal(by_name$X, X)
  expect_equal(by_index$X, X)
  expect_equal(by_name$items$location, from_matrix$items$location,
               tolerance = 1e-8)
  expect_identical(names(by_name$factors), "I01")
  expect_identical(as.character(by_name$factors$I01),
                   as.character(different$I01))
  expect_identical(by_name$items$item, by_index$items$item)

  with_roles <- cbind(
    id = seq_len(nrow(X)), grp = rep(c("A", "B"), length.out = nrow(X)), d
  )
  expect_error(
    rasch(with_roles, id = "id", factors = "grp",
          items = c("id", colnames(X))),
    "items= includes id/factor column"
  )

  same <- data.frame(I01 = X[, "I01"])
  implicit <- rasch(d, factors = same, n_groups = 4,
                    maxit = 60, tol = 1e-6)
  explicit_same <- rasch(d, factors = same, items = colnames(X),
                         n_groups = 4, maxit = 60, tol = 1e-6)
  expect_false("I01" %in% implicit$items$item)
  expect_true("I01" %in% explicit_same$items$item)
})

test_that("EFRM external factors preserve item selection and frame roles", {
  d <- simulate_efrm(40, 4, 1, 2, seed = 882)
  truth <- attr(d, "truth")
  items <- unlist(truth$item_sets, use.names = FALSE)
  item_index <- match(items, names(d))
  different <- data.frame(S1I01 = factor(rep(c("A", "B"),
                                               length.out = nrow(d))))

  expect_error(
    rasch_efrm(d, item_sets = list(set1 = items), groups = "group",
               id = "id", factors = different, n_groups = 3,
               boot_reps = 0, maxit = 30, tol = 1e-5,
               se_method = "hybrid", workers = 1),
    "share item-data names but contain different values"
  )

  by_name <- rasch_efrm(
    d, item_sets = list(set1 = items), groups = "group", id = "id",
    factors = different, items = items, n_groups = 3, boot_reps = 0,
    maxit = 30, tol = 1e-5, se_method = "hybrid", workers = 1
  )
  by_index <- rasch_efrm(
    d, item_sets = list(set1 = items), groups = "group", id = "id",
    factors = different, items = item_index, n_groups = 3, boot_reps = 0,
    maxit = 30, tol = 1e-5, se_method = "hybrid", workers = 1
  )
  from_matrix <- rasch_efrm(
    as.matrix(d[, items, drop = FALSE]), item_sets = list(set1 = items),
    groups = d$group, id = d$id, factors = different, items = items,
    n_groups = 3, boot_reps = 0, maxit = 30, tol = 1e-5,
    se_method = "hybrid", workers = 1
  )

  expect_true(by_name$est$converged)
  expect_equal(by_name$X, by_index$X)
  expect_equal(by_name$X, from_matrix$X)
  expect_equal(by_name$items$location, from_matrix$items$location,
               tolerance = 1e-8)
  expect_identical(names(by_name$factors), c("group", "S1I01"))
  expect_identical(by_name$frame_group, "group")
  expect_identical(as.character(by_name$factors$S1I01),
                   as.character(different$S1I01))

  same <- data.frame(S1I01 = d[["S1I01"]])
  implicit <- rasch_efrm(
    d, item_sets = list(set1 = items[-1L]), groups = "group", id = "id",
    factors = same, n_groups = 3, boot_reps = 0, maxit = 30, tol = 1e-5,
    se_method = "hybrid", workers = 1
  )
  expect_false(any(grepl("S1I01", implicit$items$item, fixed = TRUE)))

  explicit_same <- rasch_efrm(
    d, item_sets = list(set1 = items), groups = "group", id = "id",
    factors = same, items = items, n_groups = 3, boot_reps = 0,
    maxit = 30, tol = 1e-5, se_method = "hybrid", workers = 1
  )
  expect_true(any(grepl("S1I01", explicit_same$items$item, fixed = TRUE)))

  role_named <- same
  names(role_named) <- "group"
  expect_error(
    rasch_efrm(
      d, item_sets = list(set1 = items), groups = "group", id = "id",
      factors = role_named, items = c("group", items), n_groups = 3,
      boot_reps = 0, maxit = 30, tol = 1e-5,
      se_method = "hybrid", workers = 1
    ),
    "includes id, group, or factor column"
  )
})
