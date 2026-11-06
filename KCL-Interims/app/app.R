# Interim / futility design explorer for a time-to-event trial.
#
# Scale: d events, log HR estimate has variance 4/d. Z = -logHR_hat * sqrt(d) / 2, so
# benefit (HR < 1) gives positive Z.
#   * Final analysis at N events. Critical value from the 1-sided level alpha/2.
#   * Interim efficacy looks are 1-sided (benefit), boundaries via rpact
#     (Lan-DeMets O'Brien-Fleming, Lan-DeMets Pocock) or fixed nominal 1-sided p-values.
#   * Futility looks are 1-sided (stop if HR > threshold). Futility is NON-BINDING:
#     efficacy boundaries and the final critical value are computed as if futility
#     stopping did not exist.
#   * Success is HR <= boundary (equality counts as success); futility stop is HR > threshold.
#   * Design tab: either give the maximum events and solve for the design HR, or give a
#     target HR and solve for the number of events. Results tab: simulate and export.

library(shiny)

# --- shinylive / webR additions (no change to the app's logic) ---------------------
# Lets shinylive::export() detect rpact, and installs it from the webR repository at
# run time if the exported site does not already contain it. Harmless in native R.
if (FALSE) library(rpact)
if (identical(R.version$platform, "wasm32-unknown-emscripten") &&
    !requireNamespace("rpact", quietly = TRUE)) webr::install("rpact")
# -----------------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# ------------------------------------------------------------------ rpact availability

rpact_available <- function() {
  requireNamespace("rpact", quietly = TRUE) &&
    all(c("getDesignGroupSequential", "getGroupSequentialProbabilities") %in%
          getNamespaceExports("rpact"))
}
RPACT_MESSAGE <- paste0(
  "The rpact package is required but is not installed (or is too old to provide ",
  "getDesignGroupSequential and getGroupSequentialProbabilities). ",
  "Install it with install.packages(\"rpact\") and restart the app.")
require_rpact <- function() if (!rpact_available()) stop(RPACT_MESSAGE, call. = FALSE)

MAX_LOOKS <- 6
MAX_SIMS <- 10000
MAX_SCEN <- 30              # largest number of (effect, SD) scenarios
T_TEST_MAX_PER_ARM <- 50    # the t-test option is offered below this many subjects per arm
MIN_SIMS <- 1000
MAX_PATHS <- 30
BAR_LWD <- 5
PALETTE <- c("#0072B2", "#E69F00", "#CC79A7", "#56B4E9", "#7F7F7F",
             "#000000", "#8C564B", "#9467BD", "#17BECF", "#BCBD22")

# ---------------------------------------------------------------- input parsing

parse_hr_input <- function(txt, max_n = 10, positive = TRUE, what = "hazard ratio") {
  txt <- gsub("\\s+", "", txt %||% "")
  if (!nzchar(txt)) stop(sprintf("Enter at least one true %s.", what))
  vals <- numeric(0)
  for (p in strsplit(txt, ",", fixed = TRUE)[[1]]) {
    if (!nzchar(p)) next
    if (grepl(":", p, fixed = TRUE)) {
      r <- suppressWarnings(as.numeric(strsplit(p, ":", fixed = TRUE)[[1]]))
      if (length(r) != 3 || anyNA(r) || r[3] <= 0 || r[2] < r[1])
        stop("Write ranges as from:to:step, for example 0.6:1:0.1")
      if ((r[2] - r[1]) / r[3] > 1000) stop("Range has too many values.")
      vals <- c(vals, seq(r[1], r[2] + r[3] * 1e-8, by = r[3]))
    } else {
      v <- suppressWarnings(as.numeric(p))
      if (is.na(v)) stop(sprintf("Cannot read '%s' as a number.", p))
      vals <- c(vals, v)
    }
  }
  vals <- sort(unique(round(vals, 6)))
  if (!length(vals) || any(!is.finite(vals)) || (positive && any(vals <= 0)))
    stop(sprintf("True %ss must be %snumbers.", what, if (positive) "positive " else ""))
  if (length(vals) > max_n)
    stop(sprintf("%d %ss given; at most %d can be plotted.", length(vals), what, max_n))
  vals
}

check_events <- function(x, N, what, even = FALSE) {
  unit <- if (even) "total sample sizes" else "event counts"
  if (anyNA(x)) stop(sprintf("Enter a %s for every %s look.", if (even) "total sample size" else "event count", what))
  if (any(x != round(x)) || any(x < 1) || any(x >= N))
    stop(sprintf("%s %s must be whole numbers between 1 and N-1 (%d).", what, unit, N - 1))
  if (even && any(x %% 2 != 0))
    stop(sprintf("%s total sample sizes must be even numbers (whole subjects per arm).", what))
  if (anyDuplicated(x)) stop(sprintf("%s %s must be distinct.", what, unit))
  sort(x)
}

check_pct <- function(x, what) {
  if (anyNA(x)) stop(sprintf("Enter an information fraction for every %s look.", what))
  if (any(x <= 0) || any(x >= 100))
    stop(sprintf("%s information fractions must be greater than 0 and less than 100 (percent).", what))
  if (anyDuplicated(x)) stop(sprintf("%s information fractions must be distinct.", what))
  x
}

ev_from_pct <- function(pct, N) floor(pct / 100 * N + 0.5)   # whole events, halves round up
# Same for look sizes; with even = TRUE the result is an even total (whole subjects per arm).
look_from_pct <- function(pct, N, even = FALSE)
  if (even) 2 * floor(pct / 100 * N / 2 + 0.5) else ev_from_pct(pct, N)

whole_between <- function(x, lo, hi, what) {
  if (is.null(x) || length(x) != 1 || is.na(x) || x != round(x) || x < lo || x > hi)
    stop(sprintf("%s must be a whole number between %d and %d.", what, lo, hi))
  x
}

# ----------------------------------------------------------------- group sequential

# Exact probabilities for a Z-process with independent increments and information
# fractions `info`. Continuation region at look k is (lower[k], upper[k]); mu is the mean
# of Z_k. Returns P(stop below lower), P(stop above upper), P(reach look k).
gs_probs <- function(info, lower, upper, mu = 0) {
  K <- length(info)
  mu <- rep_len(mu, K)
  clip <- function(x) pmin(pmax(x, -40), 40)
  p <- rpact::getGroupSequentialProbabilities(rbind(clip(lower - mu), clip(upper - mu)), info)
  # rpact can return tiny negative values (about -1e-12) where the true probability is 0
  list(fut = pmax(p[1, ], 0), eff = pmax(p[3, ] - p[2, ], 0), reach = p[3, ])
}

# Efficacy boundaries from the interim information fractions alone (they do not depend on N).
efficacy_design_info <- function(info_int, method, nominal_p, alpha2) {
  require_rpact()
  a1 <- alpha2 / 2
  m <- length(info_int)
  if (m == 0) return(list(info = 1, z = qnorm(1 - a1)))
  info <- c(info_int, 1)
  if (method == "fixed") {
    zi <- qnorm(1 - nominal_p)
    spent <- sum(gs_probs(info[seq_len(m)], rep(-Inf, m), zi)$eff)
    if (spent >= a1)
      stop("The interim nominal p-values already use up the whole 1-sided alpha; lower them.")
    f <- function(cf) sum(gs_probs(info, rep(-Inf, m + 1), c(zi, cf))$eff) - a1
    z <- c(zi, uniroot(f, c(-5, 12), tol = 1e-10)$root)
  } else {
    des <- rpact::getDesignGroupSequential(kMax = m + 1, alpha = a1, sided = 1,
                                           informationRates = info,
                                           typeOfDesign = if (method == "asOF") "asOF" else "asP")
    z <- des$criticalValues
  }
  list(info = info, z = z)
}

efficacy_design <- function(N, int_ev, method, nominal_p, alpha2) {
  e <- efficacy_design_info(int_ev / N, method, nominal_p, alpha2)
  list(events = c(int_ev, N), z = e$z)
}

# Cumulative 1-sided alpha at each efficacy analysis: exact probability under HR = 1 of
# having crossed any efficacy boundary up to and including that analysis (futility ignored).
cum_alpha <- function(N, eff) {
  K <- length(eff$events)
  cumsum(gs_probs(eff$events / N, rep(-Inf, K), eff$z)$eff)
}

# Overall power (any efficacy crossing, futility ignored) when -log(HR) = theta.
overall_power <- function(events, z, N, theta) {
  K <- length(events)
  sum(gs_probs(events / N, rep(-Inf, K), z, mu = theta * sqrt(events) / 2)$eff)
}

# True HR whose overall power equals `power` (given the number of events N).
design_hr <- function(N, eff, power) {
  pw <- function(th) overall_power(eff$events, eff$z, N, th)
  hi <- 20 / sqrt(N)
  if (pw(hi) < power) stop("The target power cannot be reached with this design.")
  th <- uniroot(function(t) pw(t) - power, c(1e-8, hi), tol = 1e-12)$root
  exp(-th)
}

# Number of events needed for `power` at a target HR. Looks are given as information
# percentages (sorted ascending, interim p-values in the same order). The continuous solution
# uses the fact that the Z means are delta * sqrt(t_k) with delta = -log(HR) * sqrt(N) / 2.
# The integer solution is the smallest N at or above it for which the design with whole-event
# looks still reaches the power.
solve_events <- function(target_hr, alpha2, power, int_pct, int_p, method, fut_pct, even = FALSE) {
  require_rpact()
  if (is.null(target_hr) || length(target_hr) != 1 || !is.finite(target_hr) || target_hr <= 0 || target_hr >= 1)
    stop("The target HR must be above 0 and below 1.")
  ed <- efficacy_design_info(int_pct / 100, method, int_p, alpha2)
  K <- length(ed$info)
  pw <- function(delta) sum(gs_probs(ed$info, rep(-Inf, K), ed$z, mu = delta * sqrt(ed$info))$eff)
  if (pw(30) < power) stop("The target power cannot be reached with this design.")
  delta <- uniroot(function(d) pw(d) - power, c(1e-8, 30), tol = 1e-12)$root
  theta <- -log(target_hr)
  N_cont <- (2 * delta / theta)^2
  for (N in ceiling(N_cont) + 0:300) {
    if (even && N %% 2 != 0) next
    int_ev <- look_from_pct(int_pct, N, even); fut_ev <- look_from_pct(fut_pct, N, even)
    if (any(int_ev < 1 | int_ev >= N) || any(fut_ev < 1 | fut_ev >= N) ||
        anyDuplicated(int_ev) || anyDuplicated(fut_ev)) next
    eff <- efficacy_design(N, int_ev, method, int_p, alpha2)
    pwr <- overall_power(eff$events, eff$z, N, theta)
    if (pwr >= power - 1e-9)
      return(list(N = N, N_cont = N_cont, int_ev = int_ev, fut_ev = fut_ev, power_at_target = pwr))
  }
  stop("No whole number of events places the looks at distinct event counts with the requested power; ",
       "move the looks further apart.")
}

hr_bound <- function(z, d) exp(-z * 2 / sqrt(d))          # HR at which Z equals z
fut_z <- function(h, d) -log(h) * sqrt(d) / 2              # Z at which HR equals h

# ---- normally distributed endpoint
# With n subjects per arm, common SD s and 1:1 allocation, the difference in means has variance
# 4 s^2 / (2n). That is the form of the log HR variance 4/d with d = 2n (total sample size) and
# log HR replaced by -(benefit-positive difference)/s. The whole engine therefore runs on the
# "log HR equivalent" scale; these helpers convert to and from display units.
endpoint_of <- function(cfg) cfg$endpoint %||% "tte"
is_normal <- function(cfg) identical(endpoint_of(cfg), "normal")
dir_of <- function(cfg) cfg$dir %||% 1                      # +1 higher is better, -1 lower is better

# scenarios: eff (HR, or true treatment minus control difference), sd, mu (log HR equivalent mean)
cfg_scen <- function(cfg)
  if (is_normal(cfg)) cfg$scen else data.frame(eff = cfg$hrs, sd = NA_real_, mu = log(cfg$hrs))

# boundary in display units at an analysis of size d with critical value z (design SD for normal)
eff_bound <- function(cfg, z, d)
  if (is_normal(cfg)) dir_of(cfg) * z * 2 * cfg$sigma / sqrt(d) else hr_bound(z, d)
# Z at which the estimate equals a futility threshold (design SD for normal)
fut_z_cfg <- function(cfg, thr, d)
  if (is_normal(cfg)) dir_of(cfg) * thr / cfg$sigma * sqrt(d) / 2 else fut_z(thr, d)
# futility threshold on the log HR equivalent scale for a scenario with SD sd
fut_thr_log <- function(cfg, sd)
  if (is_normal(cfg)) -dir_of(cfg) * cfg$fut_hr / sd else log(cfg$fut_hr)

# ----------------------------------------------------------------------- simulation

# log HR estimates at the look event counts (exact joint law; independent increments)
sim_look_loghat <- function(loghr, looks, nsim) {
  dd <- diff(c(0, looks))
  K <- length(looks)
  inc <- matrix(rnorm(nsim * K, mean = rep(loghr * dd, each = nsim),
                      sd = rep(2 * sqrt(dd), each = nsim)), nrow = nsim)
  for (k in seq_len(K)[-1]) inc[, k] <- inc[, k] + inc[, k - 1]
  sweep(inc, 2, looks, "/")
}

# full running log HR estimate at every event 1..N
sim_paths <- function(loghr, N, npaths) {
  x <- matrix(rnorm(npaths * N, mean = loghr, sd = 2), nrow = npaths)
  for (j in seq_len(N)[-1]) x[, j] <- x[, j] + x[, j - 1]
  sweep(x, 2, seq_len(N), "/")
}

# Two-arm patient-level simulation on the benefit-positive scale (treatment mean delta_b, control 0,
# common SD sigma). Returns the difference in means and the t statistic at every per-arm size 1..M.
sim_t_arms <- function(delta_b, sigma, M, nrep) {
  cs <- function(a) { for (j in seq_len(ncol(a))[-1]) a[, j] <- a[, j] + a[, j - 1]; a }
  x <- matrix(rnorm(nrep * M, delta_b, sigma), nrep)
  y <- matrix(rnorm(nrep * M, 0, sigma), nrep)
  sx <- cs(x); sy <- cs(y); qx <- cs(x^2); qy <- cs(y^2)
  mm <- matrix(rep(seq_len(M), each = nrep), nrep)
  diff <- (sx - sy) / mm
  ss <- pmax((qx - sx^2 / mm) + (qy - sy^2 / mm), 0)
  tst <- diff / sqrt(ss / (2 * mm - 2) * 2 / mm)
  tst[, 1] <- NA
  list(diff = diff, t = tst)
}

# t statistics converted to the z with the same one-sided p-value
t_to_z <- function(tm, m) matrix(qnorm(pt(as.vector(tm), rep(2 * m - 2, each = nrow(tm)))), nrow = nrow(tm))

# Simulated look-level estimates for one scenario. loghat is on the log HR equivalent scale; z is
# NULL (the z statistic follows from loghat) or, for the t-test, the z with the same p-value as t.
sim_scen <- function(cfg, mu, sd, looks, nsim) {
  if (is_normal(cfg) && identical(cfg$analysis, "t")) {
    m <- looks / 2
    s <- sim_t_arms(-mu * sd, sd, max(m), nsim)
    return(list(loghat = -s$diff[, m, drop = FALSE] / sd, z = t_to_z(s$t[, m, drop = FALSE], m)))
  }
  list(loghat = sim_look_loghat(mu, looks, nsim), z = NULL)
}

# Full paths for one scenario: loghat (equivalent scale), y (display units), z at the looks (t-test).
sim_paths_scen <- function(cfg, mu, sd, N, npaths, looks) {
  if (!is_normal(cfg)) { lh <- sim_paths(mu, N, npaths); return(list(loghat = lh, y = exp(lh), z = NULL)) }
  if (identical(cfg$analysis, "t")) {
    m <- looks / 2
    s <- sim_t_arms(-mu * sd, sd, N / 2, npaths)
    idx <- pmax(1, floor(seq_len(N) / 2))
    d <- s$diff[, idx, drop = FALSE]
    return(list(loghat = -d / sd, y = dir_of(cfg) * d, z = t_to_z(s$t[, m, drop = FALSE], m)))
  }
  lh <- sim_paths(mu, N, npaths)
  list(loghat = lh, y = dir_of(cfg) * (-lh) * sd, z = NULL)
}

# Sequential decision rule applied at the look columns of `loghat`.
# Success: HR_hat <= boundary (Z >= critical value). Futility stop: HR_hat > threshold.
# Efficacy is checked before futility at a shared look. With follow_futility = FALSE the
# futility rule is ignored (needed to show the non-binding type I error and power).
apply_rules <- function(loghat, looks, N, fut_ev, fut_hr, eff_ev, eff_z, follow_futility = TRUE,
                        zmat = NULL, fut_thr = log(fut_hr)) {
  n <- nrow(loghat)
  outcome <- rep(NA_character_, n)
  stop_event <- rep(NA_real_, n)
  alive <- rep(TRUE, n)
  for (k in seq_along(looks)) {
    ev <- looks[k]
    hit <- rep(FALSE, n)
    ie <- match(ev, eff_ev)
    if (!is.na(ie)) {
      zk <- if (is.null(zmat)) -loghat[, k] * sqrt(ev) / 2 else zmat[, k]
      hit <- alive & zk >= eff_z[ie]
      outcome[hit] <- if (ev == N) "final_success" else "interim_success"
      stop_event[hit] <- ev
    }
    if (ev == N) {
      miss <- alive & !hit
      outcome[miss] <- "final_fail"
      stop_event[miss] <- ev
    } else if (follow_futility) {
      jf <- match(ev, fut_ev)
      if (!is.na(jf)) {
        fl <- alive & !hit & loghat[, k] > fut_thr[jf]
        outcome[fl] <- "futility"
        stop_event[fl] <- ev
        hit <- hit | fl
      }
    }
    alive <- alive & !hit
    if (!any(alive)) break
  }
  list(outcome = outcome, stop_event = stop_event)
}

run_design <- function(cfg) {
  require_rpact()
  set.seed(cfg$seed)
  N <- cfg$N
  eff <- efficacy_design(N, cfg$int_ev, cfg$method, cfg$int_p, cfg$alpha)
  dhr <- design_hr(N, eff, cfg$power)
  looks <- sort(unique(c(cfg$fut_ev, cfg$int_ev, N)))
  if (is_normal(cfg) && identical(cfg$analysis, "t") && min(looks) < 4)
    stop("With the t-test every analysis needs at least 2 subjects per arm (4 in total).")
  sc <- cfg_scen(cfg)
  rules <- function(lh, follow, sd, z = NULL)
    apply_rules(lh, looks, N, cfg$fut_ev, cfg$fut_hr, eff$events, eff$z, follow,
                zmat = z, fut_thr = fut_thr_log(cfg, sd))
  scen <- function(mu, sd) {
    s <- sim_scen(cfg, mu, sd, looks, cfg$nsim)
    list(follow = rules(s$loghat, TRUE, sd, s$z), ignore = rules(s$loghat, FALSE, sd, s$z))
  }
  oc <- lapply(seq_len(nrow(sc)), function(i) scen(sc$mu[i], sc$sd[i]))
  sd0 <- if (is_normal(cfg)) cfg$sigma else NA_real_       # null and design-effect scenarios use the design SD
  cmp <- list(null = scen(0, sd0), alt = scen(log(dhr), sd0))
  paths <- lapply(seq_len(nrow(sc)), function(i) {
    p <- sim_paths_scen(cfg, sc$mu[i], sc$sd[i], N, cfg$npaths, looks)
    c(p, rules(p$loghat[, looks, drop = FALSE], TRUE, sc$sd[i], p$z))
  })
  deff <- if (is_normal(cfg)) dir_of(cfg) * (-log(dhr)) * cfg$sigma else NA_real_
  list(cfg = cfg, eff = eff, dhr = dhr, deff = deff, mdd = eff_bound(cfg, eff$z[length(eff$z)], N),
       looks = looks, oc = oc, cmp = cmp, paths = paths, scen = sc)
}

# -------------------------------------------------------------------------- display

# How scenarios map to panels and colours. For the normal endpoint with several effects and several
# SDs, `facet` ("eff" or "sd") chooses the parameter that defines the panels and colour shows the other.
series_layout <- function(cfg, facet = "eff") {
  sc <- cfg_scen(cfg); n <- nrow(sc); one <- rep(1L, n)
  if (!is_normal(cfg))
    return(list(panel = one, col = seq_len(n), leg = sprintf("True HR %.3g", sc$eff), ptitle = ""))
  ue <- unique(sc$eff); us <- unique(sc$sd)
  if (length(ue) > 1 && length(us) > 1) {
    if (identical(facet, "sd"))
      list(panel = match(sc$sd, us), col = match(sc$eff, ue), leg = sprintf("True effect %.3g", ue),
           ptitle = sprintf("SD = %.3g", us))
    else
      list(panel = match(sc$eff, ue), col = match(sc$sd, us), leg = sprintf("SD %.3g", us),
           ptitle = sprintf("True effect = %.3g", ue))
  } else if (length(ue) > 1)
    list(panel = one, col = match(sc$eff, ue), leg = sprintf("True effect %.3g", ue), ptitle = sprintf("SD = %.3g", us))
  else if (length(us) > 1)
    list(panel = one, col = match(sc$sd, us), leg = sprintf("SD %.3g", us), ptitle = sprintf("True effect = %.3g", ue))
  else list(panel = one, col = one, leg = sprintf("True effect %.3g, SD %.3g", ue, us), ptitle = "")
}

# Rows of panels in the plot (two panels per row), used to size the plot.
panel_rows <- function(r, facet = "eff", ncols = 2) {
  P <- max(series_layout(r$cfg, facet)$panel)
  ceiling(P / (if (P > 1) ncols else 1))
}

# ncols: panels per row when faceted (the app stacks panels, exports use two per row)
plot_paths <- function(r, yscale = "log", axis_unit = "total", facet = "eff", ncols = 2) {
  cfg <- r$cfg
  N <- cfg$N
  nrm <- is_normal(cfg); dirn <- dir_of(cfg)
  islog <- !nrm && !identical(yscale, "linear")
  xdiv <- if (nrm && identical(axis_unit, "arm")) 2 else 1          # per-arm axis: totals divided by 2
  sc <- cfg_scen(cfg); sl <- series_layout(cfg, facet)
  P <- max(sl$panel); nc <- if (P > 1) ncols else 1; nr <- ceiling(P / nc)
  nser <- length(sl$leg)
  cols <- PALETTE[seq_len(nser)]
  start <- min(if (xdiv == 2) 2 * cfg$start else cfg$start, N)     # 'Paths start at' is in axis units
  eff_b <- eff_bound(cfg, r$eff$z, r$eff$events)
  path_y <- unlist(lapply(r$paths, function(p) p$y[, start:N]))
  null_line <- if (nrm) 0 else 1
  vals <- c(path_y, sc$eff, cfg$fut_hr, eff_b, null_line)
  vals <- vals[is.finite(vals)]
  ylim <- if (nrm) { rg <- range(vals); rg + c(-1, 1) * 0.04 * diff(rg) } else range(vals) * c(0.92, 1.08)
  worse_up <- !nrm || dirn < 0                 # does "worse" point up the y-axis?
  
  red <- adjustcolor("red", 0.55)
  green <- adjustcolor("forestgreen", 0.55)
  
  # Legend entries, ordered so they read left to right across each row. legend() fills down
  # columns, so the entries are permuted and padded to a full grid.
  leg_label <- c(sl$leg, "Futility stop", "Interim success", "Futility / fail region", "Success region")
  leg_col <- c(cols, "black", "black", red, green)
  leg_pch <- c(rep(16, nser), 4, 17, 15, 15)
  leg_lty <- c(rep(1, nser), NA, NA, NA, NA)
  n_leg <- length(leg_label)
  col_w <- max(strwidth(leg_label, units = "inches", cex = 0.85)) + 0.75   # text plus symbol and gaps
  ncl <- max(2, min(n_leg, floor(0.96 * dev.size("in")[1] / col_w)))
  nrows <- ceiling(n_leg / ncl)
  
  # Information fraction labels (percent of maximum sample size) for every futility or interim
  # analysis and the final one, drawn beneath the axis. Each label goes on the first line where it
  # does not touch the previous label on that line; the first line also holds the row title.
  pw <- dev.size("in")[1] / nc                                       # width of one panel in inches
  looks_all <- sort(unique(c(cfg$fut_ev, cfg$int_ev, N)))     # N is the final analysis (100%)
  lab <- paste0(vapply(round(100 * looks_all / N, 1), format, "", trim = TRUE), "%")
  csi <- par("csi")
  ev_per_in <- N * 1.02 / (pw - 6 * csi)                             # x units per inch of plot width
  title_end <- (0.005 * pw + strwidth("Information fraction", units = "inches", cex = 0.8) - 4.5 * csi) * ev_per_in
  half <- strwidth(lab, units = "inches", cex = 0.8) / 2 * ev_per_in
  gap <- 0.02 * N
  right <- title_end
  lvl <- integer(length(looks_all))
  for (k in seq_along(looks_all)) {
    l <- 1
    while (l <= length(right) && looks_all[k] - half[k] < right[l] + gap) l <- l + 1
    right[l] <- looks_all[k] + half[k]
    lvl[k] <- l - 1
  }
  last_line <- 2.4 + 0.9 * max(lvl)
  xlab_line <- last_line + 2.0
  ord <- as.vector(matrix(seq_len(nrows * ncl), nrow = nrows, byrow = TRUE))
  pad <- function(x, fill) c(x, rep(fill, nrows * ncl - n_leg))[ord]
  
  # Panels (two per row), then a legend panel underneath, so the legend can never sit on the data.
  lay <- matrix(c(seq_len(P), rep(0, nr * nc - P)), nr, nc, byrow = TRUE)
  graphics::layout(rbind(lay, P + 1), heights = c(rep(1, nr), 0.05 + 0.06 * nrows))
  on.exit({graphics::layout(1); par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))})
  par(oma = c(0, 0, 2.4, 0), cex = 1)
  
  xt <- pretty(c(0, N / xdiv)); xt <- xt[xt <= 1.02 * N / xdiv]
  xlab <- if (!nrm) "Number of events" else if (xdiv == 2) "Sample size per arm" else "Total sample size"
  ylab <- if (nrm) "Estimated treatment difference"
  else if (islog) "Hazard ratio (log scale)" else "Hazard ratio (linear scale)"
  ticks <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.5, 2, 3, 5)
  ticks <- ticks[ticks >= ylim[1] & ticks <= ylim[2]]
  
  for (pn in seq_len(P)) {
    rows <- which(sl$panel == pn)
    ptitle <- if (length(sl$ptitle) >= pn) sl$ptitle[pn] else ""
    par(mar = c(xlab_line + 1.3, 4.5, if (nzchar(ptitle)) 2 else 1, 1.5), xpd = FALSE)
    plot(NA, xlim = c(0, N * 1.02), ylim = ylim, log = if (islog) "y" else "", xaxt = "n", yaxt = "n",
         xlab = "", ylab = ylab)
    axis(1, at = xt * xdiv, labels = xt)
    title(xlab = xlab, line = xlab_line)
    if (nzchar(ptitle)) title(main = ptitle, cex.main = 0.9, font.main = 2, line = 0.4)
    if (!islog || length(ticks) < 2) axis(2, las = 1) else axis(2, at = ticks, las = 1)
    abline(h = null_line, lty = 2, col = "grey60")
    for (i in rows) lines(c(0, N), rep(sc$eff[i], 2), lty = 3, col = adjustcolor(cols[sl$col[i]], 0.5))
    
    par(xpd = NA)
    for (k in seq_along(looks_all))
      mtext(lab[k], side = 1, at = looks_all[k], line = 2.4 + 0.9 * lvl[k], cex = 0.8, col = "grey25")
    mtext("Information fraction", side = 1, at = grconvertX(0.005, "ndc", "user"), adj = 0,
          line = 2.4, cex = 0.8, col = "grey25")
    par(xpd = FALSE)
    
    for (i in rows) {
      p <- r$paths[[i]]; cl <- cols[sl$col[i]]
      for (j in seq_len(nrow(p$y))) {
        e <- p$stop_event[j]
        if (e >= start) lines(start:e, p$y[j, start:e], col = cl)
        if (p$outcome[j] == "futility") points(e, p$y[j, e], pch = 4, lwd = 2, col = cl, cex = 1.3)
        if (p$outcome[j] == "interim_success") points(e, p$y[j, e], pch = 17, col = cl, cex = 1.3)
      }
    }
    par(xpd = NA)
    points(rep(0, length(rows)), sc$eff[rows], pch = 16, col = cols[sl$col[rows]], cex = 2)
    par(xpd = FALSE)
    
    usr <- if (islog) 10^par("usr")[3:4] else par("usr")[3:4]
    to_red <- if (worse_up) usr[2] else usr[1]
    to_green <- if (worse_up) usr[1] else usr[2]
    for (j in seq_along(cfg$fut_ev))
      segments(cfg$fut_ev[j], cfg$fut_hr[j], cfg$fut_ev[j], to_red, lwd = BAR_LWD, col = red, lend = 1)
    K <- length(r$eff$events)
    for (k in seq_len(K - 1))
      segments(r$eff$events[k], eff_b[k], r$eff$events[k], to_green, lwd = BAR_LWD, col = green, lend = 1)
    segments(N, r$mdd, N, to_red, lwd = BAR_LWD, col = red, lend = 1)
    segments(N, r$mdd, N, to_green, lwd = BAR_LWD, col = green, lend = 1)
  }
  
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", ncol = ncl, bty = "n", cex = 0.85,
         legend = pad(leg_label, ""), col = pad(leg_col, NA),
         pch = pad(leg_pch, NA), lty = pad(leg_lty, NA))
  main <- if (!nrm)
    sprintf("Simulated running hazard ratio estimates: %d path%s per true HR, %d events at the final analysis",
            cfg$npaths, if (cfg$npaths == 1) "" else "s", N)
  else
    sprintf("Simulated running treatment difference estimates: %d path%s per scenario, %d subjects in total (%d per arm) at the final analysis",
            cfg$npaths, if (cfg$npaths == 1) "" else "s", N, N / 2)
  tcex <- min(0.95, 0.95 * 0.96 * dev.size("in")[1] / strwidth(main, units = "inches", cex = 0.95))   # shrink to fit
  mtext(main, side = 3, outer = TRUE, line = 0.6, cex = tcex)
}

pct <- function(x) sprintf("%.1f", 100 * x)

boundary_table <- function(r) {
  cfg <- r$cfg
  N <- cfg$N
  nrm <- is_normal(cfg); dirn <- dir_of(cfg)
  has_int <- length(cfg$int_ev) > 0
  ca <- cum_alpha(N, r$eff)
  rows <- list()
  for (j in seq_along(cfg$fut_ev)) {
    d <- cfg$fut_ev[j]
    rule <- if (!nrm) sprintf("Stop if HR > %.3f", cfg$fut_hr[j])
    else sprintf("Stop if effect %s %.3f", if (dirn > 0) "<" else ">", cfg$fut_hr[j])
    row <- data.frame(
      Analysis = "Futility (non-binding)", Events = as.integer(d), `Information %` = pct(d / N),
      `Z` = sprintf("%.3f", round(fut_z_cfg(cfg, cfg$fut_hr[j], d), 3) + 0), `1-sided p` = "",
      Rule = rule, check.names = FALSE)
    if (has_int) row[["Cumulative alpha, 1-sided (2-sided equivalent)"]] <- ""
    rows[[length(rows) + 1]] <- row
  }
  K <- length(r$eff$events)
  for (k in seq_len(K)) {
    d <- r$eff$events[k]
    z <- r$eff$z[k]
    b <- eff_bound(cfg, z, d)
    rule <- if (!nrm) sprintf("Success if HR <= %.3f%s", b, if (k == K) " (MDD)" else "")
    else sprintf("Success if effect %s %.3f%s", if (dirn > 0) ">=" else "<=", b, if (k == K) " (MDD)" else "")
    row <- data.frame(
      Analysis = if (k < K) "Interim efficacy" else "Final", Events = as.integer(d),
      `Information %` = pct(d / N), `Z` = sprintf("%.3f", z),
      `1-sided p` = sprintf("%.5f", 1 - pnorm(z)), Rule = rule, check.names = FALSE)
    if (has_int) row[["Cumulative alpha, 1-sided (2-sided equivalent)"]] <- sprintf("%.5f (%.5f)", ca[k], 2 * ca[k])
    rows[[length(rows) + 1]] <- row
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$Events, out$Analysis == "Final"), ]
  if (nrm) {
    out[["Per arm"]] <- as.integer(out$Events / 2)
    out <- out[, c(1, 2, ncol(out), 3:(ncol(out) - 1))]
    names(out)[2] <- "Total n"
  }
  out
}

oc_tables <- function(r) {
  cfg <- r$cfg
  nrm <- is_normal(cfg)
  sc <- cfg_scen(cfg)
  t1 <- list(); t2 <- list()
  for (i in seq_len(nrow(sc))) {
    h <- sc$eff[i]
    o <- r$oc[[i]]$follow$outcome
    se <- r$oc[[i]]$follow$stop_event
    fut <- mean(o == "futility"); ints <- mean(o == "interim_success")
    fs <- mean(o == "final_success"); ff <- mean(o == "final_fail")
    # mirrors the rules: success is at or beyond the boundary, a futility stop is strictly worse than the threshold
    benefit <- if (nrm) dir_of(cfg) * h >= 0 else h <= 1
    correct <- if (benefit) ints + fs else fut + ff
    ids <- if (nrm) data.frame(`True effect` = sprintf("%.3g", h), `True SD` = sprintf("%.3g", sc$sd[i]), check.names = FALSE)
    else data.frame(`True HR` = sprintf("%.3g", h), check.names = FALSE)
    truth <- if (nrm) { if (benefit) "benefit (better direction or 0)" else "no benefit (worse direction)" }
    else if (benefit) "benefit (HR<=1)" else "no benefit (HR>1)"
    t1[[i]] <- cbind(ids, data.frame(
      Truth = truth,
      `Simulated studies` = as.integer(cfg$nsim), `Futility stop` = pct(fut), `Interim success` = pct(ints), `Final success` = pct(fs),
      `Final not successful` = pct(ff), `Overall success` = pct(ints + fs),
      `Correct decision` = pct(correct), `Incorrect decision` = pct(1 - correct),
      check.names = FALSE))
    row <- cbind(ids, data.frame(`Simulated studies` = as.integer(cfg$nsim), check.names = FALSE))
    for (e in cfg$fut_ev) row[[sprintf("Futility @ %d", e)]] <- pct(mean(o == "futility" & se == e))
    for (e in cfg$int_ev) row[[sprintf("Interim success @ %d", e)]] <- pct(mean(o == "interim_success" & se == e))
    t2[[i]] <- row
  }
  list(by_hr = do.call(rbind, t1), by_look = do.call(rbind, t2))
}

compare_table <- function(r) {
  cfg <- r$cfg
  n <- cfg$nsim
  a1 <- cfg$alpha / 2
  succ <- function(x) mean(x$outcome %in% c("interim_success", "final_success"))
  se <- function(p) sqrt(p * (1 - p) / n)
  t1i <- succ(r$cmp$null$ignore); t1f <- succ(r$cmp$null$follow)
  pwi <- succ(r$cmp$alt$ignore);  pwf <- succ(r$cmp$alt$follow)
  verdict <- function(est, ref, s) {
    z <- (est - ref) / s
    if (abs(z) <= 3) sprintf("consistent with nominal (%.1f SE)", z) else sprintf("differs from nominal (%.1f SE)", z)
  }
  nrm <- is_normal(cfg)
  null_lab <- if (nrm) "effect = 0" else "HR = 1"
  des_lab <- if (nrm) sprintf("design effect %.3f (SD %.3g)", r$deff, cfg$sigma) else sprintf("design HR %.3f", r$dhr)
  data.frame(
    Quantity = c(sprintf("Type I error, 1-sided (%s), futility ignored", null_lab),
                 sprintf("Type I error, 1-sided (%s), futility followed", null_lab),
                 sprintf("Power at %s, futility ignored", des_lab),
                 sprintf("Power at %s, futility followed", des_lab)),
    Simulated = sprintf("%.4f", c(t1i, t1f, pwi, pwf)),
    `MC SE` = sprintf("%.4f", c(se(t1i), se(t1f), se(pwi), se(pwf))),
    Nominal = c(sprintf("%.4f (alpha/2)", a1), sprintf("<= %.4f", a1),
                sprintf("%.4f", cfg$power), sprintf("<= %.4f", cfg$power)),
    Comment = c(verdict(t1i, a1, se(t1i)),
                if (t1f <= a1 + 3 * se(t1f)) sprintf("within the nominal bound; %.4f below futility-ignored", t1i - t1f) else "EXCEEDS nominal bound",
                verdict(pwi, cfg$power, se(pwi)),
                sprintf("power loss from following futility: %.1f percentage points", 100 * (pwi - pwf))),
    check.names = FALSE)
}

# ------------------------------------------------------------------------------ export

# Every table shown in the app, as data frames, for CSV and PDF export.
result_tables <- function(r) {
  cfg <- r$cfg
  ft <- function(x) if (length(x)) paste(format(x, trim = TRUE), collapse = ", ") else "none"
  method <- c(asOF = "Lan-DeMets O'Brien-Fleming", asP = "Lan-DeMets Pocock", fixed = "Fixed nominal p-value per look")
  nrm <- is_normal(cfg); sc <- cfg_scen(cfg)
  dirtxt <- if (dir_of(cfg) > 0) "higher" else "lower"
  thr_lab <- if (!nrm) "Futility thresholds (stop if HR >)"
  else sprintf("Futility thresholds (stop if effect %s)", if (dir_of(cfg) > 0) "<" else ">")
  pairs <- list(
    if (nrm) c("Endpoint", "Normally distributed (difference in means)") else c("Endpoint", "Time-to-event (hazard ratio)"),
    if (nrm) c("Better is", dirtxt),
    if (nrm) c("Design SD", format(cfg$sigma)),
    if (nrm) c("Analysis of simulated studies", if (identical(cfg$analysis, "t")) "t-test, observed pooled SD" else "z-test, SD known"),
    c(if (nrm) "Maximum total sample size" else "Maximum events", cfg$N),
    c("2-sided alpha", cfg$alpha), c("Power", cfg$power),
    c(if (nrm) "Design effect (power, futility ignored)" else "Design HR (power, futility ignored)",
      sprintf("%.4f", if (nrm) r$deff else r$dhr)),
    c(if (nrm) "Final critical value (MDD, at the design SD)" else "Final critical value (MDD, success if HR <=)", sprintf("%.4f", r$mdd)),
    c(if (nrm) "Futility analyses, total sample size (non-binding)" else "Futility analyses, events (non-binding)", ft(cfg$fut_ev)),
    c(thr_lab, ft(cfg$fut_hr)),
    c(if (nrm) "Interim efficacy analyses, total sample size" else "Interim efficacy analyses, events", ft(cfg$int_ev)),
    c("Interim boundary", if (length(cfg$int_ev)) method[[cfg$method]] else "none"),
    if (nrm) c("True treatment effects simulated", ft(unique(sc$eff))) else c("True HRs simulated", ft(cfg$hrs)),
    if (nrm) c("True SDs simulated", ft(unique(sc$sd))),
    c(if (nrm) "Simulations per scenario" else "Simulations per HR", cfg$nsim), c("Seed", cfg$seed))
  pairs <- Filter(Negate(is.null), pairs)
  settings <- data.frame(Setting = vapply(pairs, `[`, "", 1), Value = vapply(pairs, `[`, "", 2),
                         check.names = FALSE, stringsAsFactors = FALSE)
  settings$Value <- as.character(settings$Value)
  oc <- oc_tables(r)
  out <- list(`Settings` = settings, `Analyses and boundaries` = boundary_table(r))
  out[[sprintf("Operating characteristics by %s (%% of simulated studies, futility followed)",
               if (nrm) "true effect and SD" else "true HR")]] <- oc$by_hr
  if (ncol(oc$by_look) > (if (nrm) 3 else 2)) out[["Stopping by analysis (% of simulated studies)"]] <- oc$by_look
  out[["Comparison with nominal alpha and power"]] <- compare_table(r)
  out
}

write_results_csv <- function(r, file) {
  con <- file(file, "w", encoding = "UTF-8"); on.exit(close(con))
  for (nm in names(tl <- result_tables(r))) {
    writeLines(sprintf('"%s"', gsub('"', '""', nm)), con)
    suppressWarnings(write.table(tl[[nm]], con, sep = ",", row.names = FALSE, qmethod = "double", na = ""))
    writeLines("", con)
  }
}

# Simple paged table drawing in inches (landscape letter), base graphics only.
write_results_pdf <- function(r, file, yscale = "log", axis_unit = "total", facet = "eff") {
  W <- 11; H <- 8.5; M <- 0.5
  pdf(file, width = W, height = H, onefile = TRUE)
  on.exit(dev.off())
  new_page <- function(title = NULL) {
    par(mar = c(0, 0, 0, 0), xpd = NA)
    plot.new(); plot.window(xlim = c(0, W), ylim = c(H, 0))
    text(M, 0.35, "Interim and futility analyses: design and results", adj = 0, font = 2, cex = 1.1)
    text(W - M, 0.35, format(Sys.Date()), adj = 1, cex = 0.8, col = "grey40")
    export_footer()
    M + 0.5
  }
  y <- new_page()
  tl <- result_tables(r)
  for (nm in names(tl)) {
    df <- tl[[nm]]
    head <- vapply(names(df), function(h) paste(strwrap(h, 16), collapse = "\n"), "")
    cells <- as.matrix(df); cells[is.na(cells)] <- ""
    cex <- 0.75
    cw <- function(cex) pmax(
      vapply(seq_along(head), function(j) max(strwidth(strsplit(head[j], "\n")[[1]], cex = cex, font = 2)), 0),
      vapply(seq_len(ncol(cells)), function(j) max(strwidth(cells[, j], cex = cex), 0), 0)) + 0.18
    while (sum(cw(cex)) > W - 2 * M && cex > 0.35) cex <- cex - 0.02
    wd <- cw(cex)
    lh <- abs(strheight("Mg", cex = cex)) * 1.9
    hh <- lh * max(lengths(strsplit(head, "\n")))
    need <- 0.4 + hh + lh * min(nrow(cells), 3)
    if (y + need > H - M) y <- new_page()
    text(M, y, nm, adj = c(0, 1), font = 2, cex = 0.85)
    y <- y + 0.28
    x <- M + c(0, cumsum(wd))[seq_along(wd)]
    for (j in seq_along(head)) text(x[j], y, head[j], adj = c(0, 1), cex = cex, font = 2)
    y <- y + hh + 0.04
    segments(M, y, M + sum(wd), y, lwd = 0.6)
    y <- y + 0.06
    for (i in seq_len(nrow(cells))) {
      if (y + lh > H - M) y <- new_page()
      for (j in seq_along(head)) text(x[j], y, cells[i, j], adj = c(0, 1), cex = cex)
      y <- y + lh
    }
    y <- y + 0.35
  }
  # plot on its own page
  par(xpd = FALSE, mar = c(5.1, 4.1, 4.1, 2.1), oma = c(1.8, 0, 0, 0))
  plot_paths(r, yscale, axis_unit, facet)
  export_footer()
  invisible(file)
}

write_plot_png <- function(r, file, width = 1400, height = NULL, res = 130, yscale = "log",
                           axis_unit = "total", facet = "eff") {
  height <- height %||% round(1000 * (150 + 450 * panel_rows(r, facet)) / 600)
  png(file, width = width, height = height, res = res)
  on.exit(dev.off())
  par(oma = c(1.8, 0, 0, 0))
  plot_paths(r, yscale, axis_unit, facet)
  export_footer()
  invisible(file)
}

# ------------------------------------------------------------- reading and building

# Reads and checks the Design tab inputs. Look timing is in events or in information
# percent; information percent is always used when the number of events is calculated.
read_design <- function(input) {
  endpoint <- input$endpoint %||% "tte"
  if (!endpoint %in% c("tte", "normal")) stop("Unknown endpoint.")
  nrm <- endpoint == "normal"
  mode <- (if (nrm) input$mode_n else input$mode) %||% "hr"
  if (!mode %in% c("events", "hr")) stop("Unknown calculation mode.")
  alpha <- input$alpha; power <- input$power
  if (is.null(alpha) || is.na(alpha) || alpha <= 0 || alpha >= 1) stop("Two-sided alpha must be between 0 and 1.")
  if (is.null(power) || is.na(power) || power <= alpha / 2 || power >= 1) stop("Power must be between alpha/2 and 1.")
  nf <- input$n_fut %||% 0; ni <- input$n_int %||% 0
  if (nf < 0 || nf > MAX_LOOKS || ni < 0 || ni > MAX_LOOKS || nf != round(nf) || ni != round(ni))
    stop(sprintf("Use between 0 and %d looks of each type.", MAX_LOOKS))
  timing <- if (mode == "events") ((if (nrm) input$timing_n else input$timing) %||% "events") else "info"
  key <- if (timing == "events") (if (nrm) "nn" else "ev") else "if"
  getv <- function(prefix, n) vapply(seq_len(n), function(i) as.numeric(input[[paste0(prefix, i)]] %||% NA), numeric(1))
  fut_t <- getv(paste0("fut_", key), nf); fut_hr <- getv(if (nrm) "fut_d" else "fut_hr", nf)
  int_t <- getv(paste0("int_", key), ni); int_p <- getv("int_p", ni)
  method <- input$method %||% "asOF"
  if (nf > 0 && !nrm && (anyNA(fut_hr) || any(fut_hr <= 0))) stop("Futility HR thresholds must be positive.")
  if (nf > 0 && nrm && anyNA(fut_hr)) stop("Enter a futility threshold (treatment difference) for every futility analysis.")
  if (ni > 0 && method == "fixed" && (anyNA(int_p) || any(int_p <= 0) || any(int_p >= 0.5)))
    stop("Nominal 1-sided p-values must be between 0 and 0.5.")
  dirn <- as.numeric(input$dir %||% 1); sigma <- input$sd_d
  if (nrm) {
    if (!dirn %in% c(1, -1)) stop("Choose whether higher or lower values are better.")
    if (is.null(sigma) || is.na(sigma) || sigma <= 0) stop("The design SD must be positive.")
    if (mode == "hr" && (is.null(input$target_eff) || is.na(input$target_eff))) stop("Enter the target treatment effect.")
  }
  list(endpoint = endpoint, dir = dirn, sigma = sigma, target_eff = input$target_eff,
       mode = mode, timing = timing, N = if (nrm) input$N_n else input$N, target_hr = input$target_hr,
       alpha = alpha, power = power, fut_t = fut_t, fut_hr = fut_hr, int_t = int_t, int_p = int_p, method = method)
}

# Turns a design specification into event counts (or total sample sizes), boundaries and the design HR
# (or design effect).
build_design <- function(spec) {
  require_rpact()
  nrm <- identical(spec$endpoint, "normal")
  o <- order(spec$fut_t); fut_t <- spec$fut_t[o]; fut_hr <- spec$fut_hr[o]
  o <- order(spec$int_t); int_t <- spec$int_t[o]; int_p <- spec$int_p[o]
  nf <- length(fut_t); ni <- length(int_t)
  sol <- NULL
  if (spec$mode == "events") {
    N <- spec$N
    if (nrm) {
      if (is.null(N) || length(N) != 1 || is.na(N) || N != round(N) || N < 10 || N %% 2 != 0)
        stop("Maximum total sample size must be an even whole number of at least 10.")
    } else if (is.null(N) || length(N) != 1 || is.na(N) || N != round(N) || N < 10)
      stop("Maximum events must be a whole number of at least 10.")
    if (spec$timing == "info") {
      if (nf > 0) fut_t <- look_from_pct(check_pct(fut_t, "Futility"), N, nrm)
      if (ni > 0) int_t <- look_from_pct(check_pct(int_t, "Interim"), N, nrm)
    }
    fut_ev <- if (nf > 0) check_events(fut_t, N, "Futility", nrm) else numeric(0)
    int_ev <- if (ni > 0) check_events(int_t, N, "Interim", nrm) else numeric(0)
  } else {
    if (nf > 0) fut_t <- check_pct(fut_t, "Futility")
    if (ni > 0) int_t <- check_pct(int_t, "Interim")
    target_eq <- spec$target_hr
    if (nrm) {
      th <- spec$dir * spec$target_eff / spec$sigma
      if (!is.finite(th) || th <= 0)
        stop("The target effect must be in the better direction (positive if higher is better, negative if lower is better).")
      target_eq <- exp(-th)
    }
    sol <- solve_events(target_eq, spec$alpha, spec$power, int_t, int_p, spec$method, fut_t, even = nrm)
    N <- sol$N; fut_ev <- sol$fut_ev; int_ev <- sol$int_ev
    sol$target_hr <- if (nrm) spec$target_eff else spec$target_hr
  }
  cfg <- list(N = N, alpha = spec$alpha, power = spec$power, fut_ev = fut_ev, fut_hr = fut_hr,
              int_ev = int_ev, int_p = int_p, method = spec$method)
  if (nrm) cfg <- c(cfg, list(endpoint = "normal", dir = spec$dir, sigma = spec$sigma))
  eff <- efficacy_design(N, int_ev, spec$method, int_p, spec$alpha)
  dhr <- design_hr(N, eff, spec$power)
  list(cfg = cfg, eff = eff, dhr = dhr, mdd = eff_bound(cfg, tail(eff$z, 1), N), mode = spec$mode, sol = sol,
       deff = if (nrm) spec$dir * (-log(dhr)) * spec$sigma else NA_real_)
}

# Reads the Results tab inputs. `dcfg` is the design configuration (it decides the endpoint).
read_oc <- function(input, N, dcfg = list()) {
  common <- list(nsim = whole_between(input$nsim, MIN_SIMS, MAX_SIMS, "Simulations per scenario"),
                 npaths = whole_between(input$npaths, 1, MAX_PATHS, "Plotted paths per scenario"),
                 start = whole_between(input$start, 1, N, "The path start point"),
                 seed = if (is.null(input$seed) || is.na(input$seed)) stop("Enter a seed.") else input$seed)
  if (!is_normal(dcfg)) return(c(list(hrs = parse_hr_input(input$hrs)), common))
  effs <- parse_hr_input(input$effs, positive = FALSE, what = "treatment effect")
  sds <- parse_hr_input(input$sds, what = "SD")
  if (length(effs) * length(sds) > MAX_SCEN)
    stop(sprintf("%d treatment effects times %d SDs is %d scenarios; at most %d are allowed.",
                 length(effs), length(sds), length(effs) * length(sds), MAX_SCEN))
  sc <- data.frame(eff = rep(effs, times = length(sds)), sd = rep(sds, each = length(effs)))
  sc$mu <- -dcfg$dir * sc$eff / sc$sd
  analysis <- if (N / 2 < T_TEST_MAX_PER_ARM) (input$analysis %||% "t") else "z"
  if (analysis == "t" && min(c(dcfg$fut_ev, dcfg$int_ev, N)) < 4)
    stop("With the t-test every analysis needs at least 2 subjects per arm (4 in total).")
  c(list(scen = sc, analysis = analysis), common)
}

read_cfg <- function(input) {
  d <- build_design(read_design(input))
  c(d$cfg, read_oc(input, d$cfg$N, d$cfg))
}

# ------------------------------------------------------------------------------ help

HELP <- list(
  mode = paste("Choose what the app calculates. 'Design HR' means you enter the maximum number of events and the app finds",
               "the true HR that gives the stated power. 'Number of events' means you enter a target HR and the app finds the",
               "smallest number of events that gives the stated power. Both allow for efficacy interims and ignore futility",
               "(futility is non-binding), so a less efficient boundary such as Pocock needs more events than O'Brien-Fleming for the",
               "same HR and power, or has a more extreme design HR at the same number of events."),
  N = "Total events at the final analysis. You enter it in this mode and the app never changes it.",
  target_hr = paste("The true HR the study is designed to detect with the stated power. It must be above 0 and below 1",
                    "(HR below 1 favours treatment). The look timing is then entered as information fractions, because the",
                    "event counts are not known until the calculation is done."),
  endpoint = paste("Time-to-event uses the log hazard ratio with events as the information. Normally distributed uses a difference in",
                   "means with a common SD in both arms and 1:1 allocation, with total sample size as the information. Both use the same",
                   "group sequential machinery (rpact boundaries, non-binding futility), because the difference in means with total",
                   "sample size n has variance 4 SD^2 / n, the same form as the log HR with n events."),
  dir = paste("Which direction of the treatment minus control difference is better. All effects and thresholds are entered and shown as",
              "treatment minus control. If higher is better the target effect is positive and a futility stop is an estimate below the",
              "threshold; if lower is better the target effect is negative and a futility stop is an estimate above the threshold."),
  mode_n = paste("'Sample size' means you enter a target effect and SD and the app finds the smallest even total sample size that gives the",
                 "stated power. 'Design effect' means you enter the maximum total sample size and the app finds the effect that gives the",
                 "stated power. Both allow for efficacy interims and ignore futility, which is non-binding."),
  N_n = "Total subjects at the final analysis, both arms together. It must be even (whole subjects per arm). The app never changes it.",
  target_eff = paste("True treatment minus control difference the study is designed to detect with the stated power, in the units of the",
                     "endpoint. Type one number with a decimal point (for example 5 or 2.5); a list or range is not accepted here, because the",
                     "Design tab describes a single design (the Results tab takes lists and ranges). It must be in the better direction:",
                     "positive if higher is better, negative if lower is better. Use the same units as the SD and the thresholds."),
  sd_d = paste("Standard deviation assumed in the design, common to both arms. Type one positive number with a decimal point (for example 10",
               "or 7.5), in the same units as the treatment effect; a list or range is not accepted here. It sets the sample size and the",
               "minimum detectable difference (MDD) shown in the table and on the plot. The Results tab can simulate other SDs, entered as a",
               "list (8, 10, 12) or a range written from:to:step (8:12:2)."),
  timing_n = paste("Enter the timing of each analysis as a total sample size (even) or as an information fraction (sample size at the analysis",
                   "divided by the maximum, as a percentage). Information fractions are used when the sample size is calculated; the",
                   "resulting totals are rounded to even whole numbers."),
  look_nn = "Total subjects at this analysis, both arms together: an even whole number from 2 to the maximum total sample size minus 2.",
  fut_d = paste("Stop for futility if the estimated treatment minus control difference at this analysis is worse than this value (below it if",
                "higher is better, above it if lower is better). A value of 0 stops when the estimate shows no benefit. Equality does not stop."),
  effs = paste("True treatment minus control differences to simulate: values separated by commas (0, 3), a range from:to:step (0:5:1), or both.",
               "Negative values are allowed. A value of 0 is the null. At most 10 values, and at most", MAX_SCEN, "scenarios once combined with the SDs."),
  sds = paste("True SDs to simulate, entered like the effects. The test statistic assumes the scenario SD, so each effect and SD pair has its own",
              "standardised effect (effect / SD). Using SDs different from the design SD shows how sensitive power is to a wrong SD."),
  facet_by = paste("Appears when you give several effects and several SDs. One panel is drawn for each value of the chosen parameter, and colours",
                   "show the other one. Changing it redraws the plot without a new Run."),
  analysis = paste("Offered when the maximum sample size is below", 2 * T_TEST_MAX_PER_ARM, "in total (fewer than", T_TEST_MAX_PER_ARM, "per arm).",
                   "The t-test option simulates patient-level data, estimates the pooled SD at every analysis and converts the t statistic to the",
                   "same one-sided p-value scale as the z boundary. Group sequential theory is exact only for known variance, so with the t-test",
                   "type I error and power are approximate, and the comparison table shows the simulated values with their Monte Carlo standard",
                   "errors. The design table, cumulative alpha and sample size calculation stay z-based. The MDD and the success bars are drawn",
                   "at the design SD, so they are approximate when the t-test is used or the true SD differs. The z-test option treats the scenario",
                   "SD as known. Futility uses the estimated difference only, so the analysis choice does not affect it."),
  axis_unit = paste("Show the plot x-axis as total sample size (default) or sample size per arm. Only the plot and its exports change.",
                    "'Paths start at' is read in the same units."),
  alpha = "Two-sided significance level. Interim and final boundaries use alpha/2 one-sided, so alpha/2 is the nominal 1-sided type I error.",
  power = "Overall probability of success at any efficacy analysis (interim or final) at the design HR, with futility ignored.",
  timing = paste("Enter the timing of each analysis as an event count or as an information fraction (events at the analysis",
                 "divided by maximum events, shown as a percentage). Information fractions are used automatically when the number",
                 "of events is calculated; the resulting event counts are rounded to whole events."),
  n_fut = paste("Number of futility analyses (0 to 6). Futility is 1-sided and non-binding: a study stops if the estimated HR is",
                "strictly above the threshold. Efficacy boundaries and the final critical value ignore futility stopping, so type I",
                "error is controlled whether or not a futility boundary is followed. A futility analysis may share an event count",
                "with an interim analysis."),
  look_events = "Number of events at which this analysis is done: a whole number from 1 to maximum events minus 1.",
  look_info = paste("Events at this analysis as a percentage of maximum events (above 0 and below 100). The event count is",
                    "rounded to a whole number."),
  fut_hr = "Stop for futility if the estimated HR at this analysis is strictly greater than this value. A value of 1 stops when the estimate shows no benefit.",
  n_int = paste("Number of interim efficacy analyses (0 to 6). They are 1-sided for benefit: success if the estimated HR is at or",
                "below the boundary HR."),
  method = paste("Lan-DeMets O'Brien-Fleming spends little alpha early and keeps the final critical value close to the unadjusted",
                 "one. Lan-DeMets Pocock spends alpha more evenly, so interim boundaries are easier to cross and the final one is",
                 "harder. Fixed nominal uses the 1-sided p-value you enter for each interim. In every case the final critical value",
                 "is set so that the overall 1-sided type I error equals alpha/2."),
  int_p = "1-sided p-value that must be reached at this interim to declare success, for example 0.001. Used only with the fixed nominal method.",
  hrs = paste("Values separated by commas (0.6, 0.75), a range written from:to:step (0.6:1:0.1), or both. At most 10 distinct values.",
              "HR below 1 favours treatment. The HR = 1 and design HR scenarios used for type I error and power are always simulated",
              "in addition."),
  npaths = sprintf("Number of simulated paths drawn for each true HR or scenario (1 to %d). Affects the plot only.", MAX_PATHS),
  start = paste("Where the plotted paths begin, in the units of the x-axis. For a time-to-event endpoint this is an event count. For a normal",
                "endpoint it is a total sample size if 'Plot x-axis' is set to total sample size, and a sample size per arm if it is set to",
                "per arm (so 30 per arm starts the paths at 60 subjects in total). Early estimates are very variable, so paths start",
                "part-way through. Affects the plot only."),
  yscale = paste("Scale of the hazard ratio axis. Log (default) spaces equal ratios equally, so 0.5 to 1 takes the same space as 1 to 2.",
                 "Over a narrow range of HRs a log axis looks close to linear. Changing the scale redraws the plot and the PNG and PDF",
                 "exports without needing a new Run."),
  nsim = sprintf("Number of simulated studies for each true HR or scenario (%d to %d) behind every percentage in the tables. More simulations reduce Monte Carlo error.", MIN_SIMS, MAX_SIMS),
  seed = "Random number seed. The same seed and settings reproduce the same results.",
  design_tbl = paste("Z is the standardised statistic, -log(HR) times sqrt(events) / 2, so benefit gives a positive Z; with d events the",
                     "log HR estimate has variance 4/d (1:1 allocation, Schoenfeld). Each boundary is also given as the HR at which it is",
                     "met. Cumulative alpha (shown when there are interim efficacy analyses) is the exact probability, when the true HR is",
                     "1, of having crossed an efficacy boundary at or before that analysis with futility ignored. The final value equals",
                     "alpha/2 and the 2-sided equivalent is twice the 1-sided value. Success is HR at or below the boundary; a futility",
                     "stop is HR strictly above the threshold. The final boundary is the minimum detectable difference (MDD)."),
  plot = paste("Dotted lines and left-hand dots are the true HRs (or true treatment effects for a normal endpoint, with the estimated difference on the y-axis). Red bars: futility stop region, or not successful at the final.",
               "Green bars: success region (estimated HR at or below the boundary). The final bars sit at the minimum detectable",
               "difference (MDD). Percentages beneath the axis are the information fraction at each futility or interim analysis and at the final analysis (100%).",
               "Paths stop at the first analysis where a rule fires (success is checked before futility at a shared analysis)."),
  oc_tbl = paste("Percentages of simulated studies with the futility rule followed. Truth mirrors the rules: a true HR at or below 1 is",
                 "classed as benefit and a true HR above 1 as no benefit. A success is correct for benefit; a futility stop or a",
                 "not-successful final is correct for no benefit. A true HR of exactly 1 is therefore classed as benefit, so a success",
                 "at HR = 1 is counted as correct here; the type I error at HR = 1 is compared with alpha in the last table. The table",
                 "after this one splits futility and interim stops by analysis."),
  cmp_tbl = paste("Type I error is the probability of any efficacy success when the true HR is 1, compared with alpha/2. Power is the",
                  "same probability at the design HR, compared with the power you entered. The futility-ignored rows check alpha and",
                  "power; following a non-binding futility rule can only lower both. Monte Carlo standard errors show how far the",
                  "simulated values can be expected to differ from the exact ones.")
)

with_help <- function(label, id) {
  tagList(label, tags$a(
    href = "#", class = "help-q", role = "button", title = "Click for help",
    onclick = sprintf(paste0("var e=document.getElementById('help_%s');",
                             "e.style.display=(getComputedStyle(e).display==='none')?'block':'none';return false;"), id),
    "?"))
}
help_box <- function(id, text) tags$div(id = paste0("help_", id), class = "help-box", text)

help_css <- "
.help-q{display:inline-block;width:1.4em;height:1.4em;line-height:1.4em;text-align:center;border-radius:50%;
  background:#e9ecef;color:#495057;font-size:.75em;font-weight:bold;text-decoration:none;margin-left:.5em;cursor:pointer}
.help-q:hover{background:#ced4da;text-decoration:none}
.help-box{display:none;background:#f8f9fa;border-left:3px solid #adb5bd;padding:6px 10px;margin:0 0 12px 0;
  font-size:.9em;color:#495057}
.look-row{margin-bottom:4px}
"

# One futility or interim analysis block: timing (events, total n, or information %) and a second field.
look_block <- function(kind, i, label, ev_default, nn_default, pct_default, thr, thr2 = NULL) {
  n_id <- if (kind == "fut") "n_fut" else "n_int"
  ev_id <- sprintf("%s_ev%d", kind, i)
  nn_id <- sprintf("%s_nn%d", kind, i)
  pc_id <- sprintf("%s_if%d", kind, i)
  th_id <- sprintf("%s_%s%d", kind, thr$key, i)
  ev_cond <- "input.endpoint == 'tte' && input.mode == 'events' && input.timing == 'events'"
  nn_cond <- "input.endpoint == 'normal' && input.mode_n == 'events' && input.timing_n == 'events'"
  pc_cond <- paste0("!(", ev_cond, ") && !(", nn_cond, ")")
  thr_ui <- tagList(numericInput(th_id, with_help(thr$label, th_id), thr$default, step = thr$step))
  thr_help <- help_box(th_id, thr$help)
  if (!is.null(thr$show_if)) {
    thr_ui <- conditionalPanel(thr$show_if, thr_ui)
    thr_help <- conditionalPanel(thr$show_if, thr_help)
  }
  if (!is.null(thr2)) {      # a separate threshold field for the normal endpoint
    th2_id <- sprintf("%s_%s%d", kind, thr2$key, i)
    thr_ui <- tagList(conditionalPanel("input.endpoint == 'tte'", thr_ui),
                      conditionalPanel("input.endpoint == 'normal'",
                                       numericInput(th2_id, with_help(thr2$label, th2_id), thr2$default, step = thr2$step)))
    thr_help <- tagList(conditionalPanel("input.endpoint == 'tte'", thr_help),
                        conditionalPanel("input.endpoint == 'normal'", help_box(th2_id, thr2$help)))
  }
  conditionalPanel(sprintf("input.%s >= %d", n_id, i),
                   tags$div(class = "look-row",
                            tags$strong(sprintf("%s analysis %d", label, i)),
                            fluidRow(
                              column(6,
                                     conditionalPanel(ev_cond, numericInput(ev_id, with_help("Events", ev_id), ev_default, min = 1, step = 1)),
                                     conditionalPanel(nn_cond, numericInput(nn_id, with_help("Total n", nn_id), nn_default, min = 2, step = 2)),
                                     conditionalPanel(pc_cond,
                                                      numericInput(pc_id, with_help("Info. fraction (%)", pc_id), pct_default, min = 0, max = 100, step = 0.1))),
                              column(6, thr_ui)),
                            conditionalPanel(ev_cond, help_box(ev_id, HELP$look_events)),
                            conditionalPanel(nn_cond, help_box(nn_id, HELP$look_nn)),
                            conditionalPanel(pc_cond, help_box(pc_id, HELP$look_info)),
                            thr_help))
}

design_panel <- sidebarLayout(
  sidebarPanel(width = 4,
               radioButtons("endpoint", with_help("Endpoint", "endpoint"),
                            c("Time-to-event (hazard ratio)" = "tte", "Normally distributed (difference in means)" = "normal")),
               help_box("endpoint", HELP$endpoint),
               conditionalPanel("input.endpoint == 'tte'",
                                radioButtons("mode", with_help("Calculate", "mode"),
                                             c("Number of events, given a target HR" = "hr",
                                               "Design HR, given the number of events" = "events"), selected = "hr"),
                                help_box("mode", HELP$mode),
                                conditionalPanel("input.mode == 'events'",
                                                 numericInput("N", with_help("Max events", "N"), 624, min = 10, step = 1), help_box("N", HELP$N)),
                                conditionalPanel("input.mode == 'hr'",
                                                 numericInput("target_hr", with_help("Target HR", "target_hr"), 0.75, min = 0.01, max = 0.99, step = 0.01),
                                                 help_box("target_hr", HELP$target_hr))),
               conditionalPanel("input.endpoint == 'normal'",
                                radioButtons("dir", with_help("Better is", "dir"), c("Higher values" = "1", "Lower values" = "-1"), inline = TRUE),
                                help_box("dir", HELP$dir),
                                radioButtons("mode_n", with_help("Calculate", "mode_n"),
                                             c("Sample size, given a target effect" = "hr",
                                               "Design effect, given the sample size" = "events"), selected = "hr"),
                                help_box("mode_n", HELP$mode_n),
                                conditionalPanel("input.mode_n == 'events'",
                                                 numericInput("N_n", with_help("Max total sample size", "N_n"), 170, min = 10, step = 2), help_box("N_n", HELP$N_n)),
                                conditionalPanel("input.mode_n == 'hr'",
                                                 numericInput("target_eff", with_help("Target treatment effect", "target_eff"), 5, step = 0.5),
                                                 help_box("target_eff", HELP$target_eff)),
                                numericInput("sd_d", with_help("Design SD", "sd_d"), 10, min = 0.0001, step = 0.5),
                                help_box("sd_d", HELP$sd_d)),
               fluidRow(
                 column(6, numericInput("alpha", with_help("2-sided alpha", "alpha"), 0.05, min = 0.0001, max = 0.5, step = 0.005)),
                 column(6, numericInput("power", with_help("Power", "power"), 0.9, min = 0.5, max = 0.999, step = 0.01))),
               help_box("alpha", HELP$alpha), help_box("power", HELP$power),
               hr(),
               conditionalPanel("input.endpoint == 'tte' && input.mode == 'events'",
                                radioButtons("timing", with_help("Enter analysis timing as", "timing"),
                                             c("Events" = "events", "Information fraction (%)" = "info"), inline = TRUE),
                                help_box("timing", HELP$timing)),
               conditionalPanel("input.endpoint == 'normal' && input.mode_n == 'events'",
                                radioButtons("timing_n", with_help("Enter analysis timing as", "timing_n"),
                                             c("Total sample size" = "events", "Information fraction (%)" = "info"), inline = TRUE),
                                help_box("timing_n", HELP$timing_n)),
               numericInput("n_fut", with_help(sprintf("Futility analyses (0 to %d)", MAX_LOOKS), "n_fut"), 1, min = 0, max = MAX_LOOKS, step = 1),
               help_box("n_fut", HELP$n_fut),
               lapply(seq_len(MAX_LOOKS), function(i)
                 look_block("fut", i, "Futility", c(312, 406, 468, 530, 562, 593)[i], c(86, 110, 128, 144, 154, 162)[i],
                            c(50, 65, 75, 85, 90, 95)[i],
                            list(key = "hr", label = "Stop if HR >", default = 1, step = 0.01, help = HELP$fut_hr),
                            list(key = "d", label = "Stop if worse than", default = 0, step = 0.1, help = HELP$fut_d))),
               hr(),
               numericInput("n_int", with_help(sprintf("Interim efficacy analyses (0 to %d)", MAX_LOOKS), "n_int"), 0, min = 0, max = MAX_LOOKS, step = 1),
               help_box("n_int", HELP$n_int),
               conditionalPanel("input.n_int >= 1",
                                selectInput("method", with_help("Interim boundary (1-sided, benefit)", "method"),
                                            c("Lan-DeMets O'Brien-Fleming" = "asOF", "Lan-DeMets Pocock" = "asP",
                                              "Fixed nominal p-value per look" = "fixed")),
                                help_box("method", HELP$method)),
               lapply(seq_len(MAX_LOOKS), function(i)
                 look_block("int", i, "Interim", c(312, 406, 468, 530, 562, 593)[i], c(86, 110, 128, 144, 154, 162)[i],
                            c(50, 65, 75, 85, 90, 95)[i],
                            list(key = "p", label = "Nominal 1-sided p", default = 0.001, step = 0.0005, help = HELP$int_p,
                                 show_if = "input.method == 'fixed'")))
  ),
  mainPanel(width = 8,
            h4("Design"),
            htmlOutput("design"),
            h4(with_help("Analyses and boundaries", "design_tbl")),
            help_box("design_tbl", HELP$design_tbl),
            tableOutput("bounds"),
            tags$p(class = "text-muted", style = "font-size:.9em",
                   "Boundaries, event counts and design HRs are calculated with the R package ",
                   tags$a(href = "https://www.rpact.org", target = "_blank", rel = "noopener", "rpact"),
                   " (Wassmer and Pahlke, Confirmatory Adaptive Clinical Trial Design and Analysis; ",
                   tags$a(href = "https://cran.r-project.org/package=rpact", target = "_blank", rel = "noopener",
                          "CRAN"), ")."),
            hr(),
            actionButton("show_results", "Show Results", class = "btn-primary")
  )
)

oc_panel <- sidebarLayout(
  sidebarPanel(width = 4,
               conditionalPanel("input.endpoint == 'tte'",
                                textInput("hrs", with_help("True hazard ratios", "hrs"), "0.6, 0.75, 1"),
                                help_box("hrs", HELP$hrs)),
               conditionalPanel("input.endpoint == 'normal'",
                                textInput("effs", with_help("True treatment effects", "effs"), "0, 3, 5"), help_box("effs", HELP$effs),
                                textInput("sds", with_help("True SDs", "sds"), "10"), help_box("sds", HELP$sds),
                                radioButtons("axis_unit", with_help("Plot x-axis", "axis_unit"),
                                             c("Total sample size" = "total", "Sample size per arm" = "arm"), inline = TRUE),
                                help_box("axis_unit", HELP$axis_unit),
                                uiOutput("facet_ui"), uiOutput("analysis_ui")),
               fluidRow(
                 column(6, numericInput("npaths", with_help("Plotted paths per scenario", "npaths"), 8, min = 1, max = MAX_PATHS, step = 1)),
                 column(6, numericInput("start", with_help("Paths start at", "start"), 50, min = 1, step = 1))),
               help_box("npaths", HELP$npaths), help_box("start", HELP$start),
               conditionalPanel("input.endpoint == 'tte'",
                                radioButtons("yscale", with_help("Plot y-axis scale", "yscale"), c("Log" = "log", "Linear" = "linear"), inline = TRUE),
                                help_box("yscale", HELP$yscale)),
               fluidRow(
                 column(6, numericInput("nsim", with_help("Simulations per scenario", "nsim"), MAX_SIMS, min = MIN_SIMS, max = MAX_SIMS, step = 1000)),
                 column(6, numericInput("seed", with_help("Seed", "seed"), 123, step = 1))),
               help_box("nsim", HELP$nsim), help_box("seed", HELP$seed),
               actionButton("run", "Run", class = "btn-primary")
  ),
  mainPanel(width = 8,
            uiOutput("stale"),
            div(style = "margin-bottom:10px",
                actionButton("exp_png", "Plot (PNG)", icon = icon("download")),
                actionButton("exp_csv", "All results (CSV)", icon = icon("download")),
                actionButton("exp_pdf", "All results (PDF)", icon = icon("download"))),
            h4(with_help("Simulated paths", "plot")),
            help_box("plot", HELP$plot),
            uiOutput("plot_ui"),
            h4(with_help("Operating characteristics by scenario (% of simulated studies, futility followed)", "oc_tbl")),
            help_box("oc_tbl", HELP$oc_tbl),
            tableOutput("oc1"),
            h4("Stopping by analysis (% of simulated studies)"), tableOutput("oc2"),
            h4(with_help("Comparison with nominal alpha and power", "cmp_tbl")),
            help_box("cmp_tbl", HELP$cmp_tbl),
            tableOutput("cmp")
  )
)

# ---- footer for exported PNG and PDF files: notice, app address, UTC date and time ----
APP_URL <- "https://kinnersleyn.github.io/kcl-apps/Interims/"
export_footer <- function() {
  txt <- sprintf("For educational purposes only, not for decision-making. %s Generated %s UTC",
                 APP_URL, format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"))
  fs <- 7
  w <- grid::convertWidth(grid::grobWidth(grid::textGrob(txt, gp = grid::gpar(fontsize = fs))), "npc", valueOnly = TRUE)
  if (w > 0.96) fs <- fs * 0.96 / w
  grid::grid.text(txt, y = grid::unit(2.5, "mm"), gp = grid::gpar(fontsize = fs, col = "grey35"))
}
# viewport that leaves room at the bottom of the page for the footer
footer_vp <- function() grid::viewport(y = grid::unit(0.25, "in"), height = grid::unit(1, "npc") - grid::unit(0.25, "in"), just = "bottom")

# ---- educational-use banner shown at the top of the app ----
edu_banner <- tagList(
  tags$style(HTML(".edu-banner{background:#fff4d6;border:1px solid #e0b64a;border-radius:6px;padding:7px 12px;margin:6px 0 10px;font-size:13px;line-height:1.35;color:#5c4400;}")),
  tags$div(class = "edu-banner", role = "note",
           tags$strong("For educational purposes only. "),
           "This app uses simulated or illustrative data and has not been validated. It must not be used for decision-making, including the design, conduct or interpretation of any real clinical trial or any decision about patient care.")
)

ui <- fluidPage(
  edu_banner,
  tags$head(tags$style(HTML(help_css)),
            # Saves a file sent from R (base64) through the browser, instead of Shiny's download link.
            # Shiny's own download links do not work in shinylive/webR.
            tags$script(HTML(paste0("$(document).on('shiny:connected', function(){",
                                    "Shiny.addCustomMessageHandler('save_file', function(m){",
                                    "var bin=atob(m.b64), n=bin.length, bytes=new Uint8Array(n);",
                                    "for (var i=0;i<n;i++) bytes[i]=bin.charCodeAt(i);",
                                    "var url=URL.createObjectURL(new Blob([bytes],{type:m.mime}));",
                                    "var a=document.createElement('a'); a.href=url; a.download=m.name;",
                                    "document.body.appendChild(a); a.click(); document.body.removeChild(a);",
                                    "setTimeout(function(){URL.revokeObjectURL(url);},10000);});});")))),
  titlePanel("Interim and futility analyses: event-driven trial paths and operating characteristics"),
  tabsetPanel(id = "main_tabs",
              tabPanel("Design", design_panel),
              tabPanel("Results", oc_panel))
)

# ----------------------------------------------------------------------------- server

server <- function(input, output, session) {
  as_validation <- function(e) {
    if (inherits(e, "shiny.silent.error")) stop(e)
    validate(need(FALSE, conditionMessage(e)))
  }
  # Design updates as soon as an input changes (no simulation involved).
  des <- reactive(tryCatch(build_design(read_design(input)), error = as_validation))
  oc_cfg <- reactive({
    d <- des()
    tryCatch(c(d$cfg, read_oc(input, d$cfg$N, d$cfg)), error = as_validation)
  })
  # Simulations run only when Run is pressed.
  res <- eventReactive(input$run, tryCatch(run_design(oc_cfg()), error = as_validation), ignoreNULL = FALSE)
  
  output$design <- renderUI({
    d <- des()
    N <- d$cfg$N
    pw <- sprintf("%.0f%%", 100 * d$cfg$power)
    if (is_normal(d$cfg)) {
      ge <- if (dir_of(d$cfg) > 0) "&ge;" else "&le;"
      mdd <- sprintf("Final critical value (MDD, at the design SD): success if the estimated effect is %s <b>%.3f</b> (z = %.3f).",
                     ge, d$mdd, tail(d$eff$z, 1))
      if (d$mode == "events") {
        HTML(sprintf("<p>Maximum total sample size: <b>%d</b> (%d per arm, entered). Design effect for %s power with SD %s, futility ignored: <b>%.3f</b>. %s</p>",
                     N, N / 2, pw, format(d$cfg$sigma), d$deff, mdd))
      } else {
        s <- d$sol
        HTML(sprintf(paste0("<p>Target effect <b>%s</b> (SD %s) needs <b>%d</b> subjects in total (<b>%d</b> per arm) for %s power, futility ignored ",
                            "(%.1f before rounding up to an even total and whole-subject analyses). ",
                            "Power at the target effect with %d subjects: <b>%.4f</b>. Design effect for %s power with %d subjects: <b>%.3f</b>. %s</p>"),
                     format(s$target_hr), format(d$cfg$sigma), N, N / 2, pw, s$N_cont, N, s$power_at_target, pw, N, d$deff, mdd))
      }
    } else {
      mdd <- sprintf("Final critical value (MDD): success if HR &le; <b>%.3f</b> (z = %.3f).", d$mdd, tail(d$eff$z, 1))
      if (d$mode == "events") {
        HTML(sprintf("<p>Maximum events: <b>%d</b> (entered). Design HR for %s power, futility ignored: <b>%.3f</b>. %s</p>",
                     N, pw, d$dhr, mdd))
      } else {
        s <- d$sol
        HTML(sprintf(paste0("<p>Target HR <b>%.3f</b> needs <b>%d</b> events for %s power, futility ignored ",
                            "(%.1f before rounding up to whole events and whole-event analyses). ",
                            "Power at the target HR with %d events: <b>%.4f</b>. Design HR for %s power with %d events: <b>%.3f</b>. %s</p>"),
                     s$target_hr, N, pw, s$N_cont, N, s$power_at_target, pw, N, d$dhr, mdd))
      }
    }
  })
  output$bounds <- renderTable(boundary_table(des()), striped = TRUE, spacing = "s")
  
  # Normal endpoint: several effects and several SDs need a choice of how to split the plot.
  both_multi <- reactive({
    if (!identical(input$endpoint, "normal")) return(FALSE)
    tryCatch(length(parse_hr_input(input$effs, positive = FALSE, what = "treatment effect")) > 1 &&
               length(parse_hr_input(input$sds, what = "SD")) > 1, error = function(e) FALSE)
  })
  facet_rv <- reactiveVal("eff")
  observeEvent(input$facet_by, facet_rv(input$facet_by))
  output$facet_ui <- renderUI({
    if (!both_multi()) return(NULL)
    tagList(radioButtons("facet_by", with_help("Facet plot by", "facet_by"),
                         c("Treatment effect" = "eff", "SD" = "sd"), selected = isolate(facet_rv()), inline = TRUE),
            help_box("facet_by", HELP$facet_by))
  })
  observeEvent(input$run, {
    if (both_multi())
      showModal(modalDialog(
        title = "How should the plot be split?",
        "You gave several treatment effects and several SDs, so the plot needs panels. Make one panel for each:",
        radioButtons("facet_pick", NULL, c("Treatment effect (colours show the SD)" = "eff",
                                           "SD (colours show the treatment effect)" = "sd"), selected = facet_rv()),
        footer = actionButton("facet_ok", "Show plot", class = "btn-primary"), easyClose = FALSE))
  }, ignoreInit = TRUE)
  observeEvent(input$facet_ok, {
    facet_rv(input$facet_pick %||% "eff")
    updateRadioButtons(session, "facet_by", selected = facet_rv())
    removeModal()
  })
  # The t-test option is offered when there are fewer than T_TEST_MAX_PER_ARM subjects per arm.
  output$analysis_ui <- renderUI({
    d <- tryCatch(des(), error = function(e) NULL)
    if (is.null(d) || !is_normal(d$cfg) || d$cfg$N / 2 >= T_TEST_MAX_PER_ARM) return(NULL)
    tagList(radioButtons("analysis", with_help("Analysis of simulated studies", "analysis"),
                         c("t-test, observed pooled SD" = "t", "z-test, SD known" = "z"),
                         selected = isolate(input$analysis) %||% "t"),
            help_box("analysis", HELP$analysis))
  })
  
  output$stale <- renderUI({
    r <- res()
    cur <- tryCatch(oc_cfg(), error = function(e) NULL)
    if (!is.null(cur) && !identical(cur, r$cfg))
      div(class = "alert alert-warning",
          "The design or the options on this tab have changed since the last run. Press Run to update these results.")
  })
  output$plot_ui <- renderUI(plotOutput("plot", height = sprintf("%dpx", 150 + 450 * panel_rows(res(), facet_rv(), 1))))
  output$plot <- renderPlot(plot_paths(res(), input$yscale %||% "log", input$axis_unit %||% "total", facet_rv(), 1), res = 100)
  output$oc1 <- renderTable(oc_tables(res())$by_hr, striped = TRUE, spacing = "s")
  output$oc2 <- renderTable({
    r <- res()
    t <- oc_tables(r)$by_look
    validate(need(ncol(t) > (if (is_normal(r$cfg)) 3 else 2), "No futility or interim analyses specified."))
    t
  }, striped = TRUE, spacing = "s")
  output$cmp <- renderTable(compare_table(res()), striped = TRUE, spacing = "s")
  
  observeEvent(input$show_results, updateTabsetPanel(session, "main_tabs", selected = "Results"))
  
  # Results are current only if the design and Results tab inputs still match the last Run.
  results_current <- function() {
    r <- tryCatch(res(), error = function(e) NULL)
    cur <- tryCatch(oc_cfg(), error = function(e) NULL)
    !is.null(r) && !is.null(cur) && identical(cur, r$cfg)
  }
  STALE_MSG <- "The design or options have changed since the last run, so nothing was exported. Press Run to update the results, then export."
  # Writes the file on the server side and sends its bytes to the browser, which saves it.
  send_file <- function(name, mime, writer) {
    f <- tempfile(); on.exit(unlink(f))
    writer(f)
    session$sendCustomMessage("save_file", list(name = name, mime = mime,
      b64 = jsonlite::base64_enc(readBin(f, "raw", file.size(f)))))
  }
  guarded <- function(btn, name, mime, writer) observeEvent(input[[btn]], {
    if (!results_current()) {
      showNotification(STALE_MSG, type = "warning", duration = 8)
      return()
    }
    tryCatch(send_file(name, mime, writer),
             error = function(e) showNotification(paste("Export failed:", conditionMessage(e)),
                                                  type = "error", duration = NULL))
  })
  guarded("exp_png", "interim_paths.png", "image/png", function(f)
    write_plot_png(res(), f, yscale = input$yscale %||% "log", axis_unit = input$axis_unit %||% "total",
                   facet = facet_rv()))
  guarded("exp_csv", "interim_results.csv", "text/csv", function(f) write_results_csv(res(), f))
  guarded("exp_pdf", "interim_results.pdf", "application/pdf", function(f)
    write_results_pdf(res(), f, input$yscale %||% "log", input$axis_unit %||% "total", facet_rv()))
}

missing_rpact_ui <- fluidPage(
  titlePanel("Interim and futility analyses"),
  div(class = "alert alert-danger", role = "alert", RPACT_MESSAGE))

app_object <- function() {
  if (rpact_available()) shinyApp(ui, server)
  else shinyApp(missing_rpact_ui, function(input, output, session) NULL)
}

app_object()
