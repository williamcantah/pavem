# pavem: simple Panel VAR and Panel VECM workflows.

chk_cols <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop("These columns are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

dec_p <- function(p, alpha = 0.05, null = "null hypothesis") {
  if (is.na(p)) return("No p-value available")
  if (p < alpha) {
    paste("Reject", null)
  } else {
    paste("Do not reject", null)
  }
}

first_num <- function(x) {
  if (is.null(x)) return(NA_real_)
  as.numeric(x)[1]
}

get_ht <- function(x) {
  if (inherits(x, "purtest")) x <- x$statistic
  list(
    statistic = first_num(x$statistic),
    p.value = first_num(x$p.value),
    parameter = first_num(x$parameter),
    method = if (!is.null(x$method)) x$method else NA_character_
  )
}

as_panel_list <- function(data, id, time, vars) {
  chk_cols(data, c(id, time, vars))
  data <- data[order(data[[id]], data[[time]]), , drop = FALSE]
  ids <- unique(data[[id]])
  out <- lapply(ids, function(one_id) {
    one <- data[data[[id]] == one_id, vars, drop = FALSE]
    stats::ts(one, start = min(data[data[[id]] == one_id, time]), frequency = 1)
  })
  names(out) <- ids
  out
}

companion_stab <- function(A, lags, vars) {
  k <- length(vars)
  lag_cols <- grep("\\.l[0-9]+$", colnames(A), value = TRUE)
  lag_mat <- as.matrix(A[, lag_cols, drop = FALSE])

  if (lags == 1) {
    comp <- lag_mat
  } else {
    comp <- rbind(
      lag_mat,
      cbind(
        diag(k * (lags - 1)),
        matrix(0, nrow = k * (lags - 1), ncol = k)
      )
    )
  }

  ev <- eigen(comp)$values
  data.frame(
    eigenvalue = as.character(round(ev, 5)),
    modulus = Mod(ev),
    stable = Mod(ev) < 1,
    row.names = NULL
  )
}

get_jstat <- function(model) {
  out <- tryCatch(suppressWarnings(panelvar::hansen_j_test(model)), error = function(e) NULL)
  nums <- suppressWarnings(as.numeric(unlist(out)))
  nums <- nums[is.finite(nums)]
  if (length(nums) == 0) return(NA_real_)
  nums[1]
}

fit_pvar_lag <- function(data, id, time, vars, lag,
                         transformation, steps, system,
                         max_instr, min_instr, collapse) {
  panelvar::pvargmm(
    dependent_vars = vars,
    lags = lag,
    transformation = transformation,
    data = data,
    panel_identifier = c(id, time),
    steps = steps,
    system_instruments = system,
    max_instr_dependent_vars = max_instr,
    min_instr_dependent_vars = min_instr,
    collapse = collapse
  )
}

choose_lag <- function(tab, prefer = "MBIC") {
  prefer <- toupper(prefer)
  if (!prefer %in% c("MBIC", "MAIC", "MQIC")) prefer <- "MBIC"
  ok <- is.finite(tab[[prefer]])
  if (!any(ok)) return(NA_integer_)
  tab$lag[which.min(ifelse(ok, tab[[prefer]], Inf))]
}

safe_varselect <- function(x, max_lag, type) {
  out <- tryCatch(vars::VARselect(x, lag.max = max_lag, type = type),
                  error = function(e) NULL)
  if (is.null(out)) return(NULL)
  out$selection
}

diag_dec <- function(p, stat, crit, alpha = 0.05, null = "null") {
  if (is.finite(p)) {
    if (p < alpha) paste("Reject", null) else paste("Do not reject", null)
  } else if (is.finite(stat) && is.finite(crit)) {
    if (stat > crit) paste("Reject", null) else paste("Do not reject", null)
  } else {
    "Decision not available"
  }
}

hansen_table <- function(model, alpha = 0.05) {
  h <- tryCatch(suppressWarnings(panelvar::hansen_j_test(model)), error = function(e) e)
  if (inherits(h, "error")) {
    return(data.frame(
      test = "Hansen J", statistic = NA_real_, df = NA_real_,
      p.value = NA_real_, critical.value = NA_real_,
      null = "Overidentifying restrictions are valid",
      decision = paste("Not available:", h$message),
      available = FALSE, row.names = NULL
    ))
  }
  stat <- first_num(h$statistic)
  df <- first_num(h$parameter)
  p <- first_num(h$p.value)
  crit <- if (is.finite(df)) stats::qchisq(1 - alpha, df = df) else NA_real_
  data.frame(
    test = "Hansen J",
    statistic = stat,
    df = df,
    p.value = p,
    critical.value = crit,
    null = "Overidentifying restrictions are valid",
    decision = diag_dec(p, stat, crit, alpha, "H0"),
    available = TRUE,
    row.names = NULL
  )
}

unavailable_iv_tests <- function(system = NA) {
  rbind(
    data.frame(
      test = "Sargan",
      statistic = NA_real_, df = NA_real_, p.value = NA_real_,
      critical.value = NA_real_,
      null = "Overidentifying restrictions are valid",
      decision = "Not available from panelvar::pvargmm object",
      available = FALSE, row.names = NULL
    ),
    data.frame(
      test = "Difference-in-Hansen",
      statistic = NA_real_, df = NA_real_, p.value = NA_real_,
      critical.value = NA_real_,
      null = "Instrument subset is valid",
      decision = if (isTRUE(system)) {
        "Not available from panelvar::pvargmm object"
      } else {
        "Not applicable because system instruments are not used"
      },
      available = FALSE, row.names = NULL
    )
  )
}

pvar_resid_list <- function(model) {
  lapply(model$residuals, function(x) as.matrix(x))
}

vecm_resid_list <- function(obj) {
  lapply(obj$model$L.varx, function(x) t(as.matrix(x$resid)))
}

serial_table <- function(resids, lags = c(1, 2), alpha = 0.05) {
  vars <- unique(unlist(lapply(resids, colnames)))
  rows <- list()
  idx <- 1
  for (v in vars) {
    pooled <- unlist(lapply(resids, function(m) {
      if (!v %in% colnames(m)) return(NULL)
      as.numeric(m[, v])
    }))
    pooled <- pooled[is.finite(pooled)]
    pooled <- pooled[abs(pooled) > .Machine$double.eps]
    for (lg in lags) {
      bt <- tryCatch(stats::Box.test(pooled, lag = lg, type = "Ljung-Box"),
                     error = function(e) e)
      if (inherits(bt, "error")) {
        rows[[idx]] <- data.frame(
          variable = v, test = paste0("Ljung-Box AR(", lg, ")"),
          statistic = NA_real_, df = lg, p.value = NA_real_,
          critical.value = stats::qchisq(1 - alpha, df = lg),
          null = paste("No serial correlation up to lag", lg),
          decision = paste("Not available:", bt$message),
          row.names = NULL
        )
      } else {
        stat <- first_num(bt$statistic)
        p <- first_num(bt$p.value)
        crit <- stats::qchisq(1 - alpha, df = lg)
        rows[[idx]] <- data.frame(
          variable = v, test = paste0("Ljung-Box AR(", lg, ")"),
          statistic = stat, df = lg, p.value = p,
          critical.value = crit,
          null = paste("No serial correlation up to lag", lg),
          decision = diag_dec(p, stat, crit, alpha, "H0"),
          row.names = NULL
        )
      }
      idx <- idx + 1
    }
  }
  do.call(rbind, rows)
}

plot_stability <- function(tab, title = "Stability condition", file = NULL) {
  if (!is.null(file)) grDevices::pdf(file, width = 7, height = 7)
  on.exit(if (!is.null(file)) grDevices::dev.off(), add = TRUE)
  theta <- seq(0, 2 * pi, length.out = 300)
  graphics::plot(cos(theta), sin(theta), type = "l", asp = 1,
                 xlab = "Real", ylab = "Imaginary", main = title,
                 col = "gray50")
  graphics::abline(h = 0, v = 0, col = "gray85")
  eig <- tab$eigenvalue
  if (!is.complex(eig)) eig <- suppressWarnings(as.complex(eig))
  graphics::points(Re(eig), Im(eig), pch = 19, col = ifelse(tab$stable, "steelblue", "firebrick"))
  invisible(tab)
}

pub_diag_table <- function(instrument = NULL, serial = NULL, stability = NULL) {
  out <- data.frame()

  if (!is.null(instrument) && nrow(instrument) > 0) {
    out <- rbind(out, data.frame(
      section = "Instrument validity",
      test = instrument$test,
      variable = NA_character_,
      statistic = instrument$statistic,
      df = instrument$df,
      p.value = instrument$p.value,
      critical.value = instrument$critical.value,
      decision = instrument$decision,
      row.names = NULL
    ))
  }

  if (!is.null(serial) && nrow(serial) > 0) {
    out <- rbind(out, data.frame(
      section = "Serial correlation",
      test = serial$test,
      variable = serial$variable,
      statistic = serial$statistic,
      df = serial$df,
      p.value = serial$p.value,
      critical.value = serial$critical.value,
      decision = serial$decision,
      row.names = NULL
    ))
  }

  if (!is.null(stability) && nrow(stability) > 0) {
    ev <- if ("eigenvalue" %in% names(stability)) stability$eigenvalue else stability$Eigenvalue
    mod <- if ("modulus" %in% names(stability)) stability$modulus else stability$Modulus
    out <- rbind(out, data.frame(
      section = "Stability condition",
      test = paste("Eigenvalue", seq_len(nrow(stability))),
      variable = as.character(ev),
      statistic = as.numeric(mod),
      df = NA_real_,
      p.value = NA_real_,
      critical.value = 1,
      decision = stability$decision,
      row.names = NULL
    ))
  }

  out
}

print.pav_diag <- function(x, ...) {
  cat("\nPublication-style post-estimation diagnostics\n")
  cat("------------------------------------------------\n")
  print(x$table, row.names = FALSE, ...)
  invisible(x)
}

#' Prepare panel data
#'
#' Sorts a panel data frame and optionally keeps selected variables.
#'
#' @param data Data frame.
#' @param id Cross-section identifier column name.
#' @param time Time identifier column name.
#' @param vars Optional vector of variable names to keep.
#' @return Sorted data frame.
#' @export
prep <- function(data, id, time, vars = NULL) {
  data <- as.data.frame(data)
  keep <- unique(c(id, time, vars))
  chk_cols(data, keep)
  out <- data[order(data[[id]], data[[time]]), keep, drop = FALSE]
  row.names(out) <- NULL
  out
}

#' Cross-sectional dependence tests
#'
#' Runs Pesaran CD or Breusch-Pagan LM tests for each variable using
#' `plm::pcdtest()`.
#'
#' @param data Data frame.
#' @param id Cross-section identifier column.
#' @param time Time identifier column.
#' @param vars Variables to test.
#' @param test Test type passed to `plm::pcdtest()`, usually `"cd"` or `"lm"`.
#' @param alpha Significance level.
#' @return Data frame with statistics, p-values, and decisions.
#' @export
csdt <- function(data, id, time, vars, test = "cd", alpha = 0.05) {
  chk_cols(data, c(id, time, vars))
  pdata <- plm::pdata.frame(data, index = c(id, time))
  ans <- lapply(vars, function(v) {
    fml <- stats::as.formula(paste(v, "~ 1"))
    z <- plm::pcdtest(fml, data = pdata, test = test)
    p <- first_num(z$p.value)
    data.frame(
      variable = v,
      test = test,
      statistic = first_num(z$statistic),
      p.value = p,
      reject = !is.na(p) && p < alpha,
      null = "No cross-sectional dependence",
      decision = dec_p(p, alpha, "H0"),
      row.names = NULL
    )
  })
  do.call(rbind, ans)
}

#' Panel unit-root tests
#'
#' Runs `plm::purtest()` for all requested variables and tests.
#'
#' @param data Data frame.
#' @param id Cross-section identifier column.
#' @param time Time identifier column.
#' @param vars Variables to test.
#' @param tests Unit-root tests, e.g. `"ips"`, `"madwu"`, `"levinlin"`,
#'   `"Pm"`, `"breitung"`, or `"hadri"`.
#' @param exo Deterministic component passed to `plm::purtest()`.
#' @param lags Lag choice passed to `plm::purtest()`.
#' @param alpha Significance level.
#' @return Data frame with statistics, p-values, and decisions.
#' @export
urts <- function(data, id, time, vars,
                 tests = c("ips", "madwu"),
                 exo = "intercept", lags = 1, alpha = 0.05) {
  chk_cols(data, c(id, time, vars))
  pdata <- plm::pdata.frame(data, index = c(id, time))

  rows <- list()
  k <- 1
  for (v in vars) {
    for (tt in tests) {
      z <- tryCatch(
        plm::purtest(pdata[[v]], test = tt, exo = exo, lags = lags),
        error = function(e) e
      )

      if (inherits(z, "error")) {
        rows[[k]] <- data.frame(
          variable = v, test = tt, statistic = NA_real_, p.value = NA_real_,
          reject = NA, null = "Unit root",
          decision = paste("Test failed:", z$message),
          row.names = NULL
        )
      } else {
        ht <- get_ht(z)
        p <- ht$p.value
        reject <- if (!is.na(p)) p < alpha else NA
        rows[[k]] <- data.frame(
          variable = v,
          test = tt,
          statistic = ht$statistic,
          p.value = p,
          reject = reject,
          null = if (tt == "hadri") "Stationarity" else "Unit root",
          decision = if (tt == "hadri") {
            dec_p(p, alpha, "H0 of stationarity")
          } else {
            dec_p(p, alpha, "H0 of unit root")
          },
          row.names = NULL
        )
      }
      k <- k + 1
    }
  }
  do.call(rbind, rows)
}

#' Panel cointegration tests
#'
#' Runs either Pedroni's test via `pco::pedroni99()` or the panel Johansen
#' rank test via `pvars::pcoint.JO()`.
#'
#' @param data Data frame.
#' @param id Cross-section identifier column.
#' @param time Time identifier column.
#' @param y Dependent variable for Pedroni test.
#' @param x Regressor variable for Pedroni test, or variables for Johansen.
#' @param vars Variables for Johansen test. If supplied, overrides `y` and `x`.
#' @param method `"ped"` for Pedroni or `"jo"` for panel Johansen.
#' @param lags Lag order.
#' @param type Deterministic case for `pvars::pcoint.JO()`.
#' @param alpha Significance level.
#' @return List with raw model output and a tidy decision table.
#' @export
coint <- function(data, id, time, y = NULL, x = NULL, vars = NULL,
                  method = c("ped", "jo"), lags = 1,
                  type = "Case3", alpha = 0.05) {
  method <- match.arg(method)
  if (method == "ped") {
    chk_cols(data, c(id, time, y, x))
    data <- data[order(data[[id]], data[[time]]), , drop = FALSE]
    ids <- unique(data[[id]])
    tlen <- length(unique(data[[time]]))
    n <- length(ids)
    ymat <- matrix(data[[y]], nrow = tlen, ncol = n)
    xmat <- matrix(data[[x]], nrow = tlen, ncol = n)
    raw <- pco::pedroni99(Y = ymat, X = xmat, kk = lags, type.stat = 1, ka = 2)
    stat <- as.data.frame(raw$STATISTIC)
    stat$test <- row.names(stat)
    row.names(stat) <- NULL
    stat$direction <- ifelse(grepl("ni", stat$test), "right-tail", "left-tail")
    stat$critical <- ifelse(stat$direction == "right-tail",
                            stats::qnorm(1 - alpha), stats::qnorm(alpha))
    stat$reject <- ifelse(stat$direction == "right-tail",
                          stat$standardized > stat$critical,
                          stat$standardized < stat$critical)
    stat$null <- "No cointegration"
    stat$decision <- ifelse(stat$reject, "Reject H0", "Do not reject H0")
    return(list(method = "pedroni99", raw = raw, table = stat))
  }

  if (is.null(vars)) vars <- c(y, x)
  chk_cols(data, c(id, time, vars))
  L.data <- as_panel_list(data, id, time, vars)
  raw <- pvars::pcoint.JO(L.data = L.data, lags = lags, type = type)
  pvals <- as.data.frame(raw$panel$pvals)
  stats <- as.data.frame(raw$panel$stats)
  long <- data.frame()
  for (i in seq_len(nrow(pvals))) {
    for (j in seq_len(ncol(pvals))) {
      p <- as.numeric(pvals[i, j])
      long <- rbind(long, data.frame(
        test = row.names(pvals)[i],
        rank_null = colnames(pvals)[j],
        statistic = as.numeric(stats[i, j]),
        p.value = p,
        reject = !is.na(p) && p < alpha,
        null = paste("Cointegration rank", colnames(pvals)[j]),
        decision = dec_p(p, alpha, "H0"),
        row.names = NULL
      ))
    }
  }
  list(method = "pcoint.JO", raw = raw, table = long)
}

#' Panel VAR lag selection
#'
#' Estimates each candidate lag separately and computes moment selection
#' criteria: MBIC, MAIC, and MQIC. The preferred default lag is the lag with
#' the most negative MBIC.
#'
#' @export
pvls <- function(data, id, time, vars, lags = 1:3,
                 transformation = "fod", steps = "onestep",
                 system = TRUE, max_instr = 4, min_instr = 2L,
                 collapse = TRUE, prefer = "MBIC") {
  data <- as.data.frame(data)
  chk_cols(data, c(id, time, vars))
  lags <- as.integer(lags)
  k <- length(vars)
  n <- length(unique(data[[id]]))
  tt <- length(unique(data[[time]]))
  nt <- n * tt

  rows <- list()
  models <- list()

  for (lag in lags) {
    fit <- tryCatch(
      fit_pvar_lag(data, id, time, vars, lag, transformation, steps,
                   system, max_instr, min_instr, collapse),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      rows[[as.character(lag)]] <- data.frame(
        lag = lag, J = NA_real_, K2p = k^2 * lag,
        MBIC = NA_real_, MAIC = NA_real_, MQIC = NA_real_,
        converged = FALSE, message = fit$message, row.names = NULL
      )
    } else {
      j <- get_jstat(fit)
      k2p <- k^2 * lag
      rows[[as.character(lag)]] <- data.frame(
        lag = lag, J = j, K2p = k2p,
        MBIC = j - k2p * log(nt),
        MAIC = j - 2 * k2p,
        MQIC = j - 2.01 * k2p * log(log(nt)),
        converged = TRUE, message = NA_character_, row.names = NULL
      )
      models[[as.character(lag)]] <- fit
    }
  }

  table <- do.call(rbind, rows)
  row.names(table) <- NULL
  selected <- choose_lag(table, prefer = prefer)

  structure(
    list(
      table = table,
      selected = selected,
      criterion = toupper(prefer),
      rule = "Choose the lag with the most negative criterion; MBIC is the conservative default.",
      models = models
    ),
    class = "pav_pvls"
  )
}

#' Estimate a Panel VAR
#'
#' @export
pvar <- function(data, id, time, vars, lags = 1,
                 transformation = "fod", steps = "onestep",
                 system = TRUE, max_instr = 4, min_instr = 2L,
                 collapse = TRUE, prefer = "MBIC") {
  data <- as.data.frame(data)
  lag_selection <- NULL
  if (is.character(lags) && tolower(lags[1]) == "auto") {
    lag_selection <- pvls(data, id, time, vars, lags = 1:3,
                          transformation = transformation, steps = steps,
                          system = system, max_instr = max_instr,
                          min_instr = min_instr, collapse = collapse,
                          prefer = prefer)
    lags <- lag_selection$selected
    if (is.na(lags)) stop("Automatic lag selection failed for all candidate lags.", call. = FALSE)
  }
  fit <- fit_pvar_lag(data, id, time, vars, as.integer(lags), transformation,
                      steps, system, max_instr, min_instr, collapse)
  structure(list(model = fit, vars = vars, id = id, time = time,
                 lags = as.integer(lags), lag_selection = lag_selection),
            class = "pav_pvar")
}

#' Panel VAR diagnostics
#'
#' Returns instrument validity tests, residual serial-correlation tests, and
#' stability results. Use `plot = TRUE` to draw the stability plot, or `file`
#' to save it as a PDF.
#'
#' @export
pvdi <- function(obj, alpha = 0.05, serial_lags = c(1, 2),
                 plot = TRUE, file = NULL) {
  model <- if (inherits(obj, "pav_pvar")) obj$model else obj
  stab_raw <- panelvar::stability(model)
  stab <- as.data.frame(stab_raw)
  names(stab) <- tolower(names(stab))
  stab$stable <- stab$modulus < 1
  stab$critical.value <- 1
  stab$decision <- ifelse(stab$stable,
                          "Stable: modulus below 1",
                          "Unstable: modulus at least 1")

  instrument <- rbind(
    hansen_table(model, alpha = alpha),
    unavailable_iv_tests(system = isTRUE(model$system_instruments))
  )
  serial <- serial_table(pvar_resid_list(model), lags = serial_lags, alpha = alpha)
  table <- pub_diag_table(instrument = instrument, serial = serial, stability = stab)

  if (isTRUE(plot) || !is.null(file)) {
    plot_stability(stab, title = "PVAR stability condition", file = file)
  }

  structure(list(
    table = table,
    instrument = instrument,
    hansen = instrument[instrument$test == "Hansen J", , drop = FALSE],
    serial = serial,
    stability = stab,
    stability_plot = if (isTRUE(plot) || !is.null(file)) "generated" else "not requested"
  ), class = "pav_diag")
}

#' Panel VAR impulse responses
#'
#' @export
pvir <- function(obj, n.ahead = 10, draws = 200, ci = 0.95, cores = 1) {
  model <- if (inherits(obj, "pav_pvar")) obj$model else obj
  point <- panelvar::oirf(model, n.ahead = n.ahead)
  bands <- panelvar::bootstrap_irf(
    model,
    typeof_irf = "OIRF",
    n.ahead = n.ahead,
    nof_Nstar_draws = draws,
    confidence.band = ci,
    mc.cores = cores
  )
  table <- data.frame()
  for (response in names(point)) {
    for (impulse in colnames(point[[response]])) {
      h <- seq_len(nrow(point[[response]])) - 1
      table <- rbind(table, data.frame(
        horizon = h,
        response = response,
        impulse = impulse,
        irf = point[[response]][, impulse],
        lower = bands$Lower[[response]][, impulse],
        upper = bands$Upper[[response]][, impulse],
        row.names = NULL
      ))
    }
  }
  list(point = point, bands = bands, table = table)
}

#' Panel VAR FEVD
#'
#' @export
pvfd <- function(obj, n.ahead = 10) {
  model <- if (inherits(obj, "pav_pvar")) obj$model else obj
  raw <- panelvar::fevd_orthogonal(model, n.ahead = n.ahead)
  table <- data.frame()
  for (response in names(raw)) {
    mat <- raw[[response]]
    for (impulse in colnames(mat)) {
      table <- rbind(table, data.frame(
        horizon = seq_len(nrow(mat)),
        response = response,
        impulse = impulse,
        share = mat[, impulse],
        row.names = NULL
      ))
    }
  }
  list(raw = raw, table = table)
}

#' Plot Panel VAR IRFs
#'
#' Plots PVAR impulse responses with confidence intervals.
#'
#' @param x Output from `pvir()` or a fitted object from `pvar()`.
#' @param n.ahead IRF horizon if `x` is a fitted PVAR object.
#' @param draws Bootstrap draws if `x` is a fitted PVAR object.
#' @param ci Confidence interval if `x` is a fitted PVAR object.
#' @param file Optional PDF file path.
#' @param col Line color.
#' @param shade Shaded confidence-band color.
#' @return Invisibly returns the IRF table used for plotting.
#' @export
pipl <- function(x, n.ahead = 10, draws = 200, ci = 0.95,
                 file = NULL, col = "steelblue",
                 shade = grDevices::rgb(0.2, 0.4, 0.8, 0.20)) {
  ir <- if (is.list(x) && !is.null(x$table) && all(c("irf", "lower", "upper") %in% names(x$table))) {
    x
  } else {
    pvir(x, n.ahead = n.ahead, draws = draws, ci = ci)
  }

  if (!is.null(file)) grDevices::pdf(file, width = 10, height = 7)
  on.exit(if (!is.null(file)) grDevices::dev.off(), add = TRUE)

  tab <- ir$table
  pairs <- unique(tab[, c("response", "impulse")])
  for (i in seq_len(nrow(pairs))) {
    one <- tab[tab$response == pairs$response[i] & tab$impulse == pairs$impulse[i], ]
    yr <- range(c(one$lower, one$upper, one$irf), na.rm = TRUE)
    graphics::plot(
      one$horizon, one$irf, type = "n", ylim = yr,
      xlab = "Years after shock", ylab = "Response",
      main = paste("PVAR:", pairs$response[i], "response to", pairs$impulse[i], "shock")
    )
    graphics::polygon(
      c(one$horizon, rev(one$horizon)),
      c(one$lower, rev(one$upper)),
      col = shade, border = NA
    )
    graphics::abline(h = 0, lty = 2, col = "gray50")
    graphics::lines(one$horizon, one$irf, lwd = 2, col = col)
    graphics::lines(one$horizon, one$lower, lty = 3, col = col)
    graphics::lines(one$horizon, one$upper, lty = 3, col = col)
  }
  invisible(tab)
}

#' Plot Panel VAR FEVD
#'
#' Plots PVAR forecast error variance decompositions.
#'
#' @param x Output from `pvfd()` or a fitted object from `pvar()`.
#' @param n.ahead FEVD horizon if `x` is a fitted PVAR object.
#' @param file Optional PDF file path.
#' @param cols Bar colors.
#' @return Invisibly returns the FEVD table used for plotting.
#' @export
pfpl <- function(x, n.ahead = 10, file = NULL,
                 cols = c("steelblue", "goldenrod", "darkseagreen",
                          "firebrick", "mediumpurple", "gray50")) {
  fd <- if (is.list(x) && !is.null(x$table) && "share" %in% names(x$table)) {
    x
  } else {
    pvfd(x, n.ahead = n.ahead)
  }

  if (!is.null(file)) grDevices::pdf(file, width = 10, height = 7)
  on.exit(if (!is.null(file)) grDevices::dev.off(), add = TRUE)

  tab <- fd$table
  for (response in unique(tab$response)) {
    one <- tab[tab$response == response, ]
    mat <- stats::xtabs(share ~ impulse + horizon, data = one)
    graphics::barplot(
      mat, beside = FALSE, col = cols[seq_len(nrow(mat))],
      main = paste("PVAR FEVD for", response),
      xlab = "Forecast horizon", ylab = "Proportion",
      legend.text = rownames(mat),
      args.legend = list(x = "topright", bty = "n")
    )
  }
  invisible(tab)
}

#' Panel VECM lag selection
#'
#' Selects the level-VAR lag used behind the Panel VECM. It applies
#' `vars::VARselect()` country by country and summarizes AIC, BIC/SC, HQ, and
#' FPE choices across the panel. The default selected lag is the rounded median
#' BIC/SC lag.
#'
#' @export
vels <- function(data, id, time, vars, lags = 1:4,
                 type = "const", prefer = "BIC") {
  data <- as.data.frame(data)
  chk_cols(data, c(id, time, vars))
  max_lag <- max(as.integer(lags))
  ids <- unique(data[[id]])
  rows <- list()

  for (one_id in ids) {
    one <- data[data[[id]] == one_id, vars, drop = FALSE]
    sel <- safe_varselect(one, max_lag = max_lag, type = type)
    if (is.null(sel)) {
      rows[[as.character(one_id)]] <- data.frame(
        id = as.character(one_id), AIC = NA_integer_, HQ = NA_integer_,
        SC = NA_integer_, FPE = NA_integer_, message = "VARselect failed",
        row.names = NULL
      )
    } else {
      rows[[as.character(one_id)]] <- data.frame(
        id = as.character(one_id),
        AIC = as.integer(sel[grep("AIC", names(sel))[1]]),
        HQ = as.integer(sel[grep("HQ", names(sel))[1]]),
        SC = as.integer(sel[grep("SC", names(sel))[1]]),
        FPE = as.integer(sel[grep("FPE", names(sel))[1]]),
        message = NA_character_,
        row.names = NULL
      )
    }
  }

  table <- do.call(rbind, rows)
  pref_col <- toupper(prefer)
  if (pref_col == "BIC") pref_col <- "SC"
  if (!pref_col %in% c("AIC", "HQ", "SC", "FPE")) pref_col <- "SC"
  vals <- table[[pref_col]]
  selected <- as.integer(round(stats::median(vals[is.finite(vals)], na.rm = TRUE)))
  if (!is.finite(selected)) selected <- NA_integer_

  structure(
    list(
      table = table,
      selected = selected,
      criterion = pref_col,
      rule = "Panel VECM uses the rounded median country-level selected lag; SC/BIC is the conservative default."
    ),
    class = "pav_vels"
  )
}

#' Estimate a Panel VECM
#'
#' @export
vecm <- function(data, id, time, vars, lags = 2, rank = 1,
                 type = "Case3", order = vars, prefer = "BIC") {
  data <- as.data.frame(data)
  lag_selection <- NULL
  if (is.character(lags) && tolower(lags[1]) == "auto") {
    lag_selection <- vels(data, id, time, vars, lags = 1:4,
                          type = "const", prefer = prefer)
    lags <- lag_selection$selected
    if (is.na(lags)) stop("Automatic Panel VECM lag selection failed.", call. = FALSE)
  }
  lags <- as.integer(lags)
  L.data <- as_panel_list(data, id, time, vars)
  rank_test <- pvars::pcoint.JO(L.data = L.data, lags = lags, type = type)
  fit <- pvars::pvarx.VEC(L.data = L.data, lags = lags, dim_r = rank, type = type)
  ident <- pvars::pid.chol(fit, order_k = order)
  structure(list(
    model = fit, identified = ident, rank_test = rank_test,
    L.data = L.data, vars = vars, id = id, time = time, lags = lags,
    rank = rank, type = type, order = order, lag_selection = lag_selection
  ), class = "pav_vecm")
}

#' Panel VECM diagnostics
#'
#' Returns rank-test information, residual serial-correlation tests, stability
#' results, cointegrating vectors, and the Cholesky impact matrix. Use
#' `plot = TRUE` to draw the stability plot, or `file` to save it as a PDF.
#'
#' @export
vedi <- function(obj, alpha = 0.05, serial_lags = c(1, 2),
                 plot = TRUE, file = NULL) {
  if (!inherits(obj, "pav_vecm")) stop("obj must come from vecm().", call. = FALSE)
  stab <- companion_stab(obj$model$A, obj$lags, obj$vars)
  stab$critical.value <- 1
  stab$decision <- ifelse(stab$stable,
                          "Stable: modulus below 1",
                          "Unstable: modulus at least 1")
  serial <- serial_table(vecm_resid_list(obj), lags = serial_lags, alpha = alpha)
  instrument <- data.frame(
    test = "Instrument validity",
    statistic = NA_real_, df = NA_real_, p.value = NA_real_,
    critical.value = NA_real_,
    null = "Instrument validity is not a VECM post-estimation test",
    decision = "Not applicable: Panel VECM is not estimated by GMM instruments",
    available = FALSE,
    row.names = NULL
  )
  table <- pub_diag_table(instrument = instrument, serial = serial, stability = stab)

  if (isTRUE(plot) || !is.null(file)) {
    plot_stability(stab, title = "Panel VECM stability condition", file = file)
  }

  structure(list(
    table = table,
    rank_test = obj$rank_test$panel,
    instrument = instrument,
    serial = serial,
    stability = stab,
    stable = all(stab$stable),
    beta = obj$model$beta,
    impact = obj$identified$B,
    stability_plot = if (isTRUE(plot) || !is.null(file)) "generated" else "not requested"
  ), class = "pav_diag")
}

boot_veir <- function(obj, n.ahead, draws, seed) {
  country_names <- names(obj$L.data)
  point <- vars::irf(obj$identified, n.ahead = n.ahead)
  irf_names <- names(point$irf)[-1]
  horizon <- point$irf[[1]]
  vals <- array(NA_real_, c(length(horizon), length(irf_names), draws),
                dimnames = list(horizon, irf_names, paste0("draw_", seq_len(draws))))
  set.seed(seed)
  for (d in seq_len(draws)) {
    samp <- sample(country_names, length(country_names), replace = TRUE)
    Ls <- obj$L.data[samp]
    names(Ls) <- paste0(samp, "_", seq_along(samp))
    draw <- tryCatch({
      fit <- pvars::pvarx.VEC(L.data = Ls, lags = obj$lags,
                              dim_r = obj$rank, type = obj$type)
      idfit <- pvars::pid.chol(fit, order_k = obj$order)
      vars::irf(idfit, n.ahead = n.ahead)
    }, error = function(e) NULL)
    if (!is.null(draw)) {
      for (nm in irf_names) vals[, nm, d] <- draw$irf[[nm]]
    }
  }
  list(
    lower = apply(vals, c(1, 2), stats::quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(vals, c(1, 2), stats::quantile, probs = 0.975, na.rm = TRUE),
    draws_used = sum(!is.na(vals[1, 1, ]))
  )
}

#' Panel VECM impulse responses
#'
#' @export
veir <- function(obj, n.ahead = 10, draws = 200, seed = 123) {
  if (!inherits(obj, "pav_vecm")) stop("obj must come from vecm().", call. = FALSE)
  point <- vars::irf(obj$identified, n.ahead = n.ahead)
  bands <- boot_veir(obj, n.ahead = n.ahead, draws = draws, seed = seed)
  table <- data.frame()
  for (response in obj$vars) {
    for (impulse in obj$vars) {
      nm <- paste0("epsilon[ ", impulse, " ] %->% ", response)
      table <- rbind(table, data.frame(
        horizon = point$irf[[1]],
        response = response,
        impulse = impulse,
        irf = point$irf[[nm]],
        lower = bands$lower[, nm],
        upper = bands$upper[, nm],
        row.names = NULL
      ))
    }
  }
  list(point = point, bands = bands, table = table)
}

#' Plot Panel VECM IRFs
#'
#' Plots Panel VECM impulse responses with confidence intervals.
#'
#' @param x Output from `veir()` or a fitted object from `vecm()`.
#' @param n.ahead IRF horizon if `x` is a fitted Panel VECM object.
#' @param draws Bootstrap draws if `x` is a fitted Panel VECM object.
#' @param seed Random seed for bootstrap if `x` is a fitted Panel VECM object.
#' @param file Optional PDF file path.
#' @param col Line color.
#' @param shade Shaded confidence-band color.
#' @return Invisibly returns the IRF table used for plotting.
#' @export
vipl <- function(x, n.ahead = 10, draws = 200, seed = 123,
                 file = NULL, col = "firebrick",
                 shade = grDevices::rgb(0.7, 0.3, 0.1, 0.20)) {
  ir <- if (is.list(x) && !is.null(x$table) && all(c("irf", "lower", "upper") %in% names(x$table))) {
    x
  } else {
    veir(x, n.ahead = n.ahead, draws = draws, seed = seed)
  }

  if (!is.null(file)) grDevices::pdf(file, width = 10, height = 7)
  on.exit(if (!is.null(file)) grDevices::dev.off(), add = TRUE)

  tab <- ir$table
  pairs <- unique(tab[, c("response", "impulse")])
  for (i in seq_len(nrow(pairs))) {
    one <- tab[tab$response == pairs$response[i] & tab$impulse == pairs$impulse[i], ]
    yr <- range(c(one$lower, one$upper, one$irf), na.rm = TRUE)
    graphics::plot(
      one$horizon, one$irf, type = "n", ylim = yr,
      xlab = "Years after shock", ylab = "Response",
      main = paste("Panel VECM:", pairs$response[i], "response to", pairs$impulse[i], "shock")
    )
    graphics::polygon(
      c(one$horizon, rev(one$horizon)),
      c(one$lower, rev(one$upper)),
      col = shade, border = NA
    )
    graphics::abline(h = 0, lty = 2, col = "gray50")
    graphics::lines(one$horizon, one$irf, lwd = 2, col = col)
    graphics::lines(one$horizon, one$lower, lty = 3, col = col)
    graphics::lines(one$horizon, one$upper, lty = 3, col = col)
  }
  invisible(tab)
}

#' Panel VECM FEVD
#'
#' @export
vefd <- function(obj, n.ahead = 10) {
  if (!inherits(obj, "pav_vecm")) stop("obj must come from vecm().", call. = FALSE)
  ir <- vars::irf(obj$identified, n.ahead = n.ahead)
  table <- data.frame()
  for (response in obj$vars) {
    sq <- matrix(NA_real_, nrow = length(ir$irf[[1]]), ncol = length(obj$vars))
    colnames(sq) <- obj$vars
    for (impulse in obj$vars) {
      nm <- paste0("epsilon[ ", impulse, " ] %->% ", response)
      sq[, impulse] <- ir$irf[[nm]]^2
    }
    cum <- apply(sq, 2, cumsum)
    den <- rowSums(cum)
    share <- sweep(cum, 1, den, "/")
    for (impulse in obj$vars) {
      table <- rbind(table, data.frame(
        horizon = ir$irf[[1]],
        response = response,
        impulse = impulse,
        share = share[, impulse],
        row.names = NULL
      ))
    }
  }
  table
}

#' Plot Panel VECM FEVD
#'
#' Plots Panel VECM forecast error variance decompositions.
#'
#' @param x Output from `vefd()` or a fitted object from `vecm()`.
#' @param n.ahead FEVD horizon if `x` is a fitted Panel VECM object.
#' @param file Optional PDF file path.
#' @param cols Bar colors.
#' @return Invisibly returns the FEVD table used for plotting.
#' @export
vfpl <- function(x, n.ahead = 10, file = NULL,
                 cols = c("steelblue", "goldenrod", "darkseagreen",
                          "firebrick", "mediumpurple", "gray50")) {
  tab <- if (is.data.frame(x) && "share" %in% names(x)) {
    x
  } else {
    vefd(x, n.ahead = n.ahead)
  }

  if (!is.null(file)) grDevices::pdf(file, width = 10, height = 7)
  on.exit(if (!is.null(file)) grDevices::dev.off(), add = TRUE)

  for (response in unique(tab$response)) {
    one <- tab[tab$response == response, ]
    mat <- stats::xtabs(share ~ impulse + horizon, data = one)
    graphics::barplot(
      mat, beside = FALSE, col = cols[seq_len(nrow(mat))],
      main = paste("Panel VECM FEVD for", response),
      xlab = "Forecast horizon", ylab = "Proportion",
      legend.text = rownames(mat),
      args.legend = list(x = "topright", bty = "n")
    )
  }
  invisible(tab)
}

coef_status <- function(p.value, alpha) {
  ifelse(is.na(p.value), "Not available",
         ifelse(p.value < alpha, "Significant*", "Insignificant"))
}

mg_table <- function(coef, var, dep, alpha, labels = NULL) {
  n <- rep(NA_integer_, length(coef))
  if (!is.null(dim(var)) && length(dim(var)) == 3) {
    n <- apply(var[dep, , , drop = FALSE], 2, function(z) sum(is.finite(z)))
  }
  if (length(n) != length(coef)) n <- rep(NA_integer_, length(coef))
  se <- sqrt(as.numeric(var) / n)
  se[!is.finite(se)] <- NA_real_
  t.stat <- as.numeric(coef) / se
  df <- pmax(n - 1, 1)
  p.value <- 2 * stats::pt(-abs(t.stat), df = df)
  p.value[!is.finite(p.value)] <- NA_real_
  data.frame(
    variable = if (is.null(labels)) names(coef) else labels,
    coefficient = as.numeric(coef),
    std.error = se,
    t.statistic = t.stat,
    p.value = p.value,
    status = coef_status(p.value, alpha),
    row.names = NULL
  )
}

norm_beta_tables <- function(beta_coef, dep, alpha) {
  vars <- dimnames(beta_coef)[[1]]
  dep_pos <- match(paste0(dep, ".l1"), vars)
  if (is.na(dep_pos)) dep_pos <- match(dep, vars)
  if (is.na(dep_pos)) stop("Could not find dependent variable in beta.", call. = FALSE)

  rows <- list()
  for (v in vars[-dep_pos]) {
    v_pos <- match(v, vars)
    vals <- -beta_coef[v_pos, 1, ] / beta_coef[dep_pos, 1, ]
    vals <- vals[is.finite(vals)]
    coef <- mean(vals, na.rm = TRUE)
    se <- stats::sd(vals, na.rm = TRUE) / sqrt(length(vals))
    if (!is.finite(se) || se == 0) se <- NA_real_
    t.stat <- coef / se
    p.value <- 2 * stats::pt(-abs(t.stat), df = max(length(vals) - 1, 1))
    rows[[v]] <- data.frame(
      variable = sub("\\.l[0-9]+$", "", v),
      coefficient = coef,
      std.error = se,
      t.statistic = t.stat,
      p.value = p.value,
      status = coef_status(p.value, alpha),
      row.names = NULL
    )
  }
  do.call(rbind, rows)
}

short_labels <- function(x) {
  out <- sub("^const$", "Constant", x)
  out <- sub("\\.l([0-9]+)$", " lag \\1", out)
  ifelse(out == "Constant", out, paste0("D(", out, ")"))
}

#' Extract Panel VECM Long-Run and Short-Run Results
#'
#' Creates publication-style coefficient tables from a fitted Panel VECM. A
#' separate long-run and short-run table is created for each dependent variable.
#'
#' @param obj Object returned by `vecm()`.
#' @param alpha Significance level used for the status column.
#' @return A list with `long_run`, `short_run`, and `combined` tables.
#' @export
vres <- function(obj, alpha = 0.05) {
  if (!inherits(obj, "pav_vecm")) stop("obj must come from vecm().", call. = FALSE)
  mg <- obj$model$MG_VECM
  vars <- obj$vars

  if (is.null(mg$beta$coef) || is.null(mg$alpha$coef) || is.null(mg$GAMMA$coef)) {
    stop("The fitted VECM object does not contain the required mean-group coefficients.",
         call. = FALSE)
  }

  long_run <- list()
  short_run <- list()
  combined <- data.frame()

  for (dep in vars) {
    long <- norm_beta_tables(mg$beta$coef, dep = dep, alpha = alpha)
    long$dependent <- dep
    long$scenario <- "Long-run"
    long_run[[dep]] <- long[, c("dependent", "scenario", "variable",
                                "coefficient", "std.error", "t.statistic",
                                "p.value", "status")]

    dep_pos <- match(dep, rownames(mg$alpha$mean))
    a_coef <- mg$alpha$mean[dep_pos, , drop = TRUE]
    a_var <- mg$alpha$var[dep_pos, , drop = TRUE]
    a_n <- sum(is.finite(mg$alpha$coef[dep_pos, 1, ]))
    a_se <- sqrt(as.numeric(a_var) / a_n)
    a_t <- as.numeric(a_coef) / a_se
    a_p <- 2 * stats::pt(-abs(a_t), df = max(a_n - 1, 1))
    alpha_tab <- data.frame(
      variable = "Cointegrating equation",
      coefficient = as.numeric(a_coef),
      std.error = a_se,
      t.statistic = a_t,
      p.value = a_p,
      status = coef_status(a_p, alpha),
      row.names = NULL
    )

    g_coef <- mg$GAMMA$mean[dep_pos, , drop = TRUE]
    g_var <- mg$GAMMA$var[dep_pos, , drop = TRUE]
    g_n <- apply(mg$GAMMA$coef[dep_pos, , , drop = FALSE], 2,
                 function(z) sum(is.finite(z)))
    g_se <- sqrt(as.numeric(g_var) / g_n)
    g_se[!is.finite(g_se)] <- NA_real_
    g_t <- as.numeric(g_coef) / g_se
    g_p <- 2 * stats::pt(-abs(g_t), df = pmax(g_n - 1, 1))
    gamma_tab <- data.frame(
      variable = short_labels(names(g_coef)),
      coefficient = as.numeric(g_coef),
      std.error = g_se,
      t.statistic = g_t,
      p.value = g_p,
      status = coef_status(g_p, alpha),
      row.names = NULL
    )

    short <- rbind(alpha_tab, gamma_tab)
    short$dependent <- dep
    short$scenario <- "Short-run"
    short_run[[dep]] <- short[, c("dependent", "scenario", "variable",
                                  "coefficient", "std.error", "t.statistic",
                                  "p.value", "status")]

    combined <- rbind(combined, long_run[[dep]], short_run[[dep]])
  }
  rownames(combined) <- NULL

  structure(
    list(long_run = long_run, short_run = short_run, combined = combined),
    class = "pav_vres"
  )
}

print.pav_vres <- function(x, ...) {
  for (dep in names(x$long_run)) {
    cat("\n", dep, " as the dependent variable\n", sep = "")
    cat("\nLong-run scenario\n")
    print(x$long_run[[dep]][, -c(1, 2)], row.names = FALSE)
    cat("\nShort-run scenario (error correction)\n")
    print(x$short_run[[dep]][, -c(1, 2)], row.names = FALSE)
  }
  invisible(x)
}

#' Publication-ready Panel VECM results
#'
#' @export
vpub <- function(obj, n.ahead = 10, draws = 200) {
  if (!inherits(obj, "pav_vecm")) stop("obj must come from vecm().", call. = FALSE)
  list(
    model_coefficients = obj$model$A,
    beta = obj$model$beta,
    vecm_results = vres(obj),
    diagnostics = vedi(obj),
    irf = veir(obj, n.ahead = n.ahead, draws = draws)$table,
    fevd = vefd(obj, n.ahead = n.ahead)
  )
}
