# Proportional and non-proportional hazards demonstrator
# Two-arm trial, time-to-event endpoint

library(shiny)
library(ggplot2)
library(dplyr)
library(survival)
if (FALSE) library(gtable)   # lets shinylive::export() detect gtable, used via gtable::
if (FALSE) library(jsonlite) # used via jsonlite:: for downloads

arm_cols <- c(Control = "#1b6ca8", Experimental = "#d9561f")

# ---- helpers ---------------------------------------------------------------

step_val <- function(t, cuts, vals) {
  if (length(cuts) == 0) rep(vals[1], length(t)) else vals[findInterval(t, cuts) + 1]
}

# Piecewise-constant representation of both hazards on a common grid.
# The grid contains every user cut time, so step scenarios are exact; smooth
# presets (waning) are approximated by a fine grid.
make_cells <- function(tmax, c_cuts, c_vals, hr_fun, hr_cuts, n_fine = 400) {
  bp <- c(seq(0, tmax, length.out = n_fine + 1), c_cuts, hr_cuts)
  bp <- sort(unique(round(bp[bp >= 0 & bp <= tmax], 10)))
  a <- bp[-length(bp)]
  b <- bp[-1]
  mid <- (a + b) / 2
  hc <- step_val(mid, c_cuts, c_vals)
  hr <- hr_fun(mid)
  data.frame(a = a, b = b, hc = hc, hr = hr, he = hc * hr)
}

# Exponential dropout hazard from a dropout percentage per interval (rpact convention:
# phi = -log(1 - dropoutRate) / dropoutTime)
dropout_rate <- function(pct, interval) -log(1 - pct / 100) / interval + 0

# Event times by inversion of the piecewise cumulative hazard
simulate_arm <- function(n, a, b, h, lam, tmax, arm) {
  Hc <- c(0, cumsum(h * (b - a)))
  e <- rexp(n)
  idx <- findInterval(e, Hc)
  t_event <- rep(Inf, n)
  ok <- idx <= length(a)
  t_event[ok] <- a[idx[ok]] + (e[ok] - Hc[idx[ok]]) / h[idx[ok]]
  cens <- if (lam > 0) rexp(n, lam) else rep(Inf, n)
  data.frame(arm = arm,
             time = pmin(t_event, cens, tmax),
             event = as.integer(t_event <= cens & is.finite(t_event)))
}

stack_plots <- function(plots, heights) {
  grobs <- lapply(plots, ggplot2::ggplotGrob)
  g <- Reduce(function(x, y) rbind(x, y, size = "max"), grobs)
  rows <- sort(g$layout$t[grepl("^panel", g$layout$name)])
  g$heights[rows] <- grid::unit(heights, "null")
  g
}

hr_breaks <- function(lim) {
  cand <- c(0.05, 0.1, 0.2, 0.25, 0.33, 0.5, 0.67, 0.75, 1, 1.25, 1.5, 2, 3, 4, 5, 10, 20)
  b <- cand[cand >= lim[1] & cand <= lim[2]]
  if (length(b) < 3) b <- unique(signif(exp(seq(log(lim[1]), log(lim[2]), length.out = 5)), 2))
  b
}

haz_grob_fn <- function(cells, tmax) {
  xs <- c(cells$a, tmax)
  ext <- function(v) c(v, tail(v, 1))
  d <- rbind(data.frame(time = xs, haz = ext(cells$hc), arm = "Control"),
             data.frame(time = xs, haz = ext(cells$he), arm = "Experimental"))
  d$arm <- factor(d$arm, levels = names(arm_cols))
  r <- data.frame(time = xs, hr = ext(cells$hr))
  xscale <- scale_x_continuous(limits = c(0, tmax), expand = expansion(mult = 0.02))
  
  p1 <- ggplot(d, aes(time, haz, colour = arm)) +
    geom_step(linewidth = 1) +
    scale_colour_manual(values = arm_cols, name = NULL) +
    xscale + expand_limits(y = 0) +
    labs(x = NULL, y = "Hazard rate") +
    theme_bw(base_size = 13) +
    theme(legend.position = "top")
  
  rng <- range(c(r$hr, 1))
  p2 <- ggplot(r, aes(time, hr)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_step(linewidth = 1, colour = "black") +
    scale_y_log10(limits = c(rng[1] / 1.1, rng[2] * 1.1), breaks = hr_breaks,
                  labels = function(x) format(x, trim = TRUE, drop0trailing = TRUE)) +
    xscale +
    labs(x = "Time", y = "Hazard ratio (log scale)") +
    theme_bw(base_size = 13)
  
  stack_plots(list(p1, p2), c(2, 1))
}

km_grob_fn <- function(sim, cells, tmax, step, ci = FALSE) {
  fit <- survfit(Surv(time, event) ~ arm, data = sim, conf.type = "log-log", conf.int = 0.95)
  s <- summary(fit, censored = TRUE)
  km <- data.frame(time = s$time, surv = s$surv, lo = s$lower, hi = s$upper, ncens = s$n.censor,
                   arm = factor(sub("^arm=", "", as.character(s$strata)), levels = names(arm_cols)))
  km$lo <- ifelse(is.na(km$lo), km$surv, km$lo)
  km$hi <- ifelse(is.na(km$hi), km$surv, km$hi)
  km <- bind_rows(data.frame(time = 0, surv = 1, lo = 1, hi = 1, ncens = 0,
                             arm = factor(names(arm_cols), levels = names(arm_cols))), km) |>
    arrange(arm, time)
  
  # step-shaped band: each estimate holds until the next time point
  ribbon <- NULL
  if (ci) {
    ci_df <- km |>
      group_by(arm) |>
      group_modify(~ {
        n <- nrow(.x)
        data.frame(time = as.vector(rbind(.x$time, c(.x$time[-1], .x$time[n]))),
                   lo = rep(.x$lo, each = 2), hi = rep(.x$hi, each = 2))
      }) |>
      ungroup()
    ribbon <- geom_ribbon(data = ci_df, aes(time, ymin = lo, ymax = hi, fill = arm),
                          alpha = 0.15, colour = NA, show.legend = FALSE)
  }
  
  theo <- function(h, arm) {
    H <- c(0, cumsum(h * (cells$b - cells$a)))
    data.frame(time = c(cells$a, tmax), surv = exp(-H), arm = arm)
  }
  th <- bind_rows(theo(cells$hc, "Control"), theo(cells$he, "Experimental"))
  th$arm <- factor(th$arm, levels = names(arm_cols))
  
  times <- seq(0, tmax, by = step)
  xscale <- scale_x_continuous(limits = c(0, tmax), breaks = times, expand = expansion(mult = 0.02))
  
  p1 <- ggplot() +
    ribbon +
    geom_line(data = th, aes(time, surv, colour = arm), linewidth = 2.4, alpha = 0.3,
              show.legend = FALSE) +
    geom_step(data = km, aes(time, surv, colour = arm), linewidth = 0.8) +
    geom_point(data = filter(km, ncens > 0), aes(time, surv, colour = arm),
               shape = 124, size = 2.5, show.legend = FALSE) +
    scale_colour_manual(values = arm_cols, name = NULL) +
    scale_fill_manual(values = arm_cols, guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0.02)) +
    xscale +
    labs(x = NULL, y = "Survival probability",
         subtitle = paste0("Solid: simulated KM, ticks mark censored times. Faint line: true survival.",
                           if (ci) "\nShaded: pointwise 95% CI (Greenwood variance, log-log transformation)." else "")) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.subtitle = element_text(size = 10))
  
  rk <- expand.grid(time = times, arm = names(arm_cols), stringsAsFactors = FALSE)
  rk$n <- mapply(function(t, a) sum(sim$time[sim$arm == a] >= t), rk$time, rk$arm)
  rk$arm <- factor(rk$arm, levels = rev(names(arm_cols)))
  
  p2 <- ggplot(rk, aes(time, arm, label = n, colour = arm)) +
    geom_text(size = 4) +
    scale_colour_manual(values = arm_cols, guide = "none") +
    xscale +
    labs(x = "Time", y = NULL, title = "Number at risk") +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.ticks.y = element_blank(), axis.line.x = element_line(),
          plot.title = element_text(size = 11, face = "bold"))
  
  stack_plots(list(p1, p2), c(5, 1))
}

cox_text <- function(sim, true_hr = NULL) {
  fit <- tryCatch(suppressWarnings(coxph(Surv(time, event) ~ arm, data = sim)),
                  error = function(e) NULL)
  if (is.null(fit) || !is.finite(coef(fit))) {
    return("Cox model could not be fitted (too few events in an arm).")
  }
  s <- summary(fit)
  paste0(
    sprintf("Cox HR (experimental vs control): %.2f, 95%% CI %.2f to %.2f, Wald p = %s.",
            s$conf.int[1, "exp(coef)"], s$conf.int[1, "lower .95"], s$conf.int[1, "upper .95"],
            format.pval(s$coefficients[1, "Pr(>|z|)"], digits = 3, eps = 1e-4)),
    if (!is.null(true_hr)) sprintf(" True HR: %.2f.", true_hr) else ""
  )
}

save_grob <- function(g, file, w, h) {
  png(file, width = w, height = h + 0.25, units = "in", res = 300)
  on.exit(dev.off())
  grid::grid.newpage()
  grid::pushViewport(footer_vp())
  grid::grid.draw(g)
  grid::popViewport()
  export_footer()
}

fmt <- function(x) vapply(x, function(v) format(signif(v, 3), trim = TRUE, scientific = FALSE), "")

# Add wrapped text under a gtable; note_pt is the height reserved for it in points
add_note <- function(g, txt, note_pt) {
  tg <- grid::textGrob(paste(txt, collapse = "\n"),
                       x = grid::unit(6, "pt"), y = grid::unit(1, "npc") - grid::unit(4, "pt"),
                       hjust = 0, vjust = 1,
                       gp = grid::gpar(fontsize = 9, col = "grey25", lineheight = 1.15))
  g <- gtable::gtable_add_rows(g, grid::unit(note_pt, "pt"))
  gtable::gtable_add_grob(g, tg, t = nrow(g), l = 1, r = ncol(g), clip = "off", name = "note")
}

# PNG of a plot with the parameter description written underneath it
save_km_png <- function(g, lines, file, w = 8, h = 8, width = 130) {
  txt <- unlist(lapply(lines, strwrap, width = width, exdent = 2))
  note_pt <- length(txt) * 12.5 + 12
  save_grob(add_note(g, txt, note_pt), file, w, h + note_pt / 72)
}

# ---- help ------------------------------------------------------------------

help_btn <- function(id) {
  actionButton(id, "?", class = "btn-info btn-xs", title = "Help",
               style = "margin-left:6px; padding:0 7px; border-radius:50%; font-weight:bold;")
}
lab <- function(text, id) tagList(text, help_btn(id))

help_content <- list(
  help_ctrl = list(
    title = "Control arm: median survival time or hazard rate",
    body = tagList(
      p("Choose 'Median survival time' to enter the control arm on the survival scale. The app converts ",
        "it to a hazard rate assuming an exponential distribution: hazard = ln(2) / median, ",
        "equivalently median = ln(2) / hazard. For example, a median of 12 gives a hazard of 0.0578 ",
        "per unit time. The conversion is shown in a table under the inputs."),
      p("With a piecewise control arm the conversion is applied within each interval, so each median ",
        "describes the constant hazard in that interval and not the median of the overall survival curve."),
      p("The experimental arm hazard is always the control hazard multiplied by the hazard ratio. ",
        "A zero hazard cannot be entered as a median, so use the hazard rate option for that case.")
    )
  ),
  help_ctrl_type = list(
    title = "Control hazard over time: constant or piecewise",
    body = tagList(
      p("Constant means one hazard applies for the whole follow-up, from time 0 to the maximum time."),
      p("Piecewise constant means follow-up is split into 2 to 5 intervals at cut times you choose. The ",
        "hazard is constant within each interval and jumps to a new value at each cut time, giving a step ",
        "function. For example, with cut times 12 and 24 and a maximum time of 36 you enter three values: ",
        "one for 0 to 12, one for 12 to 24 and one for 24 to 36. Cut times must increase strictly and lie ",
        "between 0 and the maximum time."),
      p("Whenever the control hazard changes, the experimental hazard changes with it, because it is the ",
        "control hazard multiplied by the hazard ratio at that time.")
    )
  ),
  help_hr = list(
    title = "Hazard ratio",
    body = tagList(
      p("The hazard ratio (HR) is the experimental hazard divided by the control hazard at the same time, ",
        "so values below 1 favour the experimental arm. The lower panel of the plot shows the HR on a log ",
        "axis with a dashed line at 1."),
      p("Proportional: the HR is constant over time and equals the target HR."),
      p("Non-proportional, piecewise HR: choose 2 to 5 intervals, each with its own HR, so the HR changes in ",
        "steps at the cut times you set. The target HR is not used."),
      p("Non-proportional, preset: delayed effect (HR = 1 until the delay time, then the target HR); ",
        "waning effect (HR starts at the target and rises linearly to 1 at the end of the waning period, ",
        "then stays at 1); crossing hazards (an early HR before the crossing time, usually above 1, then ",
        "the target HR).")
    )
  ),
  help_cens = list(
    title = "Censoring (dropout)",
    body = tagList(
      p("As in rpact, the percentage applies to every interval of the stated length. Dropout times are ",
        "exponential with rate -ln(1 - p) / interval, where p is the proportion. So 10% per 12 means that, ",
        "in each 12 time units, 10% of patients still under observation drop out (ignoring events). ",
        "The same rate is used in both arms."),
      p("Censoring at the maximum time is administrative and comes on top of dropout, so the censored ",
        "percentage reported in the results table is higher than the dropout percentage.")
    )
  ),
  help_sim = list(
    title = "Generate new simulation",
    body = p("The button replaces the random seed with a new random value and redraws the simulated data and ",
             "Kaplan-Meier curve. Enter that seed again to reproduce a draw.")
  ),
  help_ci = list(
    title = "Confidence interval method",
    body = p("Pointwise 95% confidence intervals from Greenwood's variance estimate on the complementary ",
             "log-log scale (survfit conf.type = 'log-log'), which keeps the limits within 0 and 1. ",
             "They are pointwise, not simultaneous bands.")
  )
)

# ---- UI --------------------------------------------------------------------

tab_css <- "
.navbar-default { background-color:#1f3a5f; border-color:#1f3a5f; border-radius:0; }
.navbar-default .navbar-brand { color:#ffffff; font-weight:600; }
.navbar-default .navbar-nav > li > a { color:#dce6f2; font-weight:700; font-size:18px; padding:15px 32px; }
.navbar-default .navbar-nav > li > a:hover { background-color:#2d5086; color:#ffffff; }
.navbar-default .navbar-nav > .active > a,
.navbar-default .navbar-nav > .active > a:hover,
.navbar-default .navbar-nav > .active > a:focus { background-color:#d9561f; color:#ffffff; }
"

# ---- footer for exported PNG and PDF files: notice, app address, UTC date and time ----
APP_URL <- "https://kinnersleyn.github.io/kcl-apps/Hazards/"
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

ui <- navbarPage(
  "Proportional and non-proportional hazards",
  id = "tabs",
  windowTitle = "Proportional and non-proportional hazards",
  header = tagList(edu_banner, tags$head(tags$style(HTML(tab_css)),
    tags$script(HTML(paste0("$(document).on('shiny:connected', function(){",
                                "Shiny.addCustomMessageHandler('save_file', function(m){",
                                "var bin=atob(m.b64), n=bin.length, bytes=new Uint8Array(n);",
                                "for (var i=0;i<n;i++) bytes[i]=bin.charCodeAt(i);",
                                "var url=URL.createObjectURL(new Blob([bytes],{type:m.mime}));",
                                "var a=document.createElement('a'); a.href=url; a.download=m.name;",
                                "document.body.appendChild(a); a.click(); document.body.removeChild(a);",
                                "setTimeout(function(){URL.revokeObjectURL(url);},10000);});});"))))),
  
  tabPanel(
    "Design",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Time and control arm"),
        numericInput("tmax", "Maximum time", 36, min = 1, step = 1),
        radioButtons("ctrl_type", lab("Control hazard over time", "help_ctrl_type"),
                     c("Constant" = "constant", "Piecewise constant" = "piecewise")),
        radioButtons("ctrl_scale", lab("Specify control arm as", "help_ctrl"),
                     c("Median survival time" = "median", "Hazard rate" = "hazard"),
                     selected = "median"),
        conditionalPanel(
          "input.ctrl_type == 'constant' && input.ctrl_scale == 'hazard'",
          numericInput("h0", "Control hazard rate (per unit time)", 0.0866, min = 0, step = 0.005)
        ),
        conditionalPanel(
          "input.ctrl_type == 'constant' && input.ctrl_scale == 'median'",
          numericInput("med0", "Control median survival time", 8, min = 0.001, step = 0.5)
        ),
        conditionalPanel(
          "input.ctrl_type == 'piecewise'",
          selectInput("ctrl_k", "Number of intervals", 2:5, selected = 3),
          uiOutput("ctrl_ui")
        ),
        conditionalPanel(
          "input.ctrl_scale == 'median'",
          tableOutput("ctrl_conv")
        ),
        h4("Hazard ratio"),
        radioButtons("hr_type", lab("Pattern over time", "help_hr"),
                     c("Proportional" = "ph",
                       "Non-proportional: piecewise HR" = "pw",
                       "Non-proportional: preset" = "preset")),
        conditionalPanel(
          "input.hr_type != 'pw'",
          numericInput("hr_t", "Target hazard ratio (experimental vs control)", 0.7, min = 0.01, step = 0.05)
        ),
        conditionalPanel(
          "input.hr_type == 'pw'",
          selectInput("hr_k", "Number of HR intervals", 2:5, selected = 3),
          uiOutput("hr_ui")
        ),
        conditionalPanel(
          "input.hr_type == 'preset'",
          selectInput("preset", "Preset",
                      c("Delayed effect" = "delay", "Waning effect" = "wane", "Crossing hazards" = "cross")),
          conditionalPanel(
            "input.preset == 'delay'",
            numericInput("p_delay", "Delay (HR = 1 before this time)", 6, min = 0.1)
          ),
          conditionalPanel(
            "input.preset == 'wane'",
            numericInput("p_wane", "Waning period (HR moves linearly from target to 1 by this time)", 24, min = 0.1)
          ),
          conditionalPanel(
            "input.preset == 'cross'",
            numericInput("p_cross", "Crossing time (target HR applies after this time)", 12, min = 0.1),
            numericInput("p_early", "Early HR (before crossing)", 1.3, min = 0.01, step = 0.05)
          )
        )
      ),
      mainPanel(
        width = 8,
        div(style = "margin-bottom:10px; display:flex; flex-wrap:wrap; gap:6px; align-items:center;",
            tags$b("Save:"),
            actionButton("dl_haz", "Hazard plot (PNG)", icon = icon("download"))),
        plotOutput("haz_plot", height = "560px")
      )
    ),
    div(style = "margin:10px 0 30px 0; text-align:center;",
        actionButton("show_results", "Show Results", class = "btn-primary btn-lg"))
  ),
  
  tabPanel(
    "Results",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Simulation"),
        numericInput("n_arm", "Patients per arm (maximum 5000)", 300, min = 10, max = 5000, step = 50),
        numericInput("cens", lab("Censoring (dropout) %, identical in each arm", "help_cens"),
                     10, min = 0, max = 99, step = 1),
        numericInput("cens_int", "Censoring interval (time units)", 12, min = 0.001),
        numericInput("risk_step", "Numbers-at-risk interval", 6, min = 0.1),
        numericInput("seed", "Random seed", 123, step = 1),
        div(style = "margin-bottom:15px;",
            actionButton("resim", "Generate new simulation"), help_btn("help_sim")),
        div(style = "display:flex; align-items:center;",
            checkboxInput("km_ci", "Show 95% confidence intervals", FALSE), help_btn("help_ci"))
      ),
      mainPanel(
        width = 8,
        div(style = "margin-bottom:10px; display:flex; flex-wrap:wrap; gap:6px; align-items:center;",
            tags$b("Save:"),
            actionButton("dl_km", "KM plot (PNG)", icon = icon("download")),
            actionButton("dl_dat", "Simulated data (CSV)", icon = icon("download"))),
        plotOutput("km_plot", height = "620px"),
        uiOutput("sim_desc"),
        tableOutput("sim_tbl"),
        textOutput("cox_txt")
      )
    )
  )
)

# ---- server ----------------------------------------------------------------

server <- function(input, output, session) {
  
  num_inputs <- function(prefix, n) {
    vapply(seq_len(n), function(i) {
      x <- input[[paste0(prefix, i)]]
      if (is.null(x) || is.na(x)) NA_real_ else as.numeric(x)
    }, numeric(1))
  }
  
  # Last valid piecewise control hazards (hazard scale), kept so that switching
  # between hazard and median input carries the values across.
  ctrl_store <- reactiveVal(NULL)
  observe({
    req(input$ctrl_type == "piecewise")
    k <- as.integer(input$ctrl_k)
    v <- if (identical(input$ctrl_scale, "median")) log(2) / num_inputs("cm_", k) else num_inputs("ch_", k)
    if (!anyNA(v) && all(is.finite(v)) && all(v >= 0)) ctrl_store(v)
  })
  
  observeEvent(input$ctrl_scale, {
    if (identical(input$ctrl_scale, "median")) {
      if (isTRUE(input$h0 > 0)) updateNumericInput(session, "med0", value = signif(log(2) / input$h0, 4))
    } else if (isTRUE(input$med0 > 0)) {
      updateNumericInput(session, "h0", value = signif(log(2) / input$med0, 4))
    }
  }, ignoreInit = TRUE)
  
  lapply(names(help_content), function(id) {
    observeEvent(input[[id]], {
      showModal(modalDialog(title = help_content[[id]]$title, help_content[[id]]$body,
                            easyClose = TRUE, footer = modalButton("Close")))
    }, ignoreInit = TRUE)
  })
  
  observeEvent(input$show_results, {
    updateNavbarPage(session, "tabs", selected = "Results")
  })
  
  observeEvent(input$resim, {
    set.seed(NULL)
    updateNumericInput(session, "seed", value = sample.int(999999L, 1L))
  })
  
  output$ctrl_ui <- renderUI({
    k <- as.integer(input$ctrl_k)
    scale <- input$ctrl_scale
    tm <- isolate(input$tmax)
    base <- isolate(if (identical(scale, "median")) log(2) / input$med0 else input$h0)
    req(k, !is.na(tm), !is.na(base), is.finite(base))
    cuts <- round(tm * seq_len(k - 1) / k, 2)
    mult <- c(1.5, 1, 0.7, 0.5, 0.4)[seq_len(k)]
    st <- isolate(ctrl_store())
    hz <- if (!is.null(st) && length(st) == k && all(st > 0)) st else base * mult
    cut_inputs <- lapply(seq_len(k - 1), function(i)
      numericInput(paste0("cc_", i), paste("Control cut time", i), cuts[i], min = 0))
    val_inputs <- if (identical(scale, "median")) {
      lapply(seq_len(k), function(i)
        numericInput(paste0("cm_", i), paste("Control median, interval", i),
                     signif(log(2) / hz[i], 4), min = 0.001, step = 0.5))
    } else {
      lapply(seq_len(k), function(i)
        numericInput(paste0("ch_", i), paste("Control hazard, interval", i),
                     signif(hz[i], 3), min = 0, step = 0.005))
    }
    tagList(cut_inputs, val_inputs)
  })
  
  output$ctrl_conv <- renderTable({
    req(input$ctrl_scale == "median")
    if (input$ctrl_type == "constant") {
      m <- input$med0
      req(isTRUE(m > 0))
      data.frame(Interval = "All times", Median = format(m, digits = 4),
                 `Hazard rate (ln 2 / median)` = sprintf("%.5f", log(2) / m), check.names = FALSE)
    } else {
      k <- as.integer(input$ctrl_k)
      m <- num_inputs("cm_", k)
      cuts <- num_inputs("cc_", k - 1)
      req(!anyNA(m), !anyNA(cuts), all(m > 0), isTRUE(input$tmax > 0))
      edges <- format(c(0, cuts, input$tmax), trim = TRUE, digits = 4)
      data.frame(Interval = paste(edges[-length(edges)], "to", edges[-1]),
                 Median = format(m, digits = 4, trim = TRUE),
                 `Hazard rate (ln 2 / median)` = sprintf("%.5f", log(2) / m), check.names = FALSE)
    }
  })
  
  output$hr_ui <- renderUI({
    k <- as.integer(input$hr_k)
    tm <- isolate(input$tmax)
    tg <- isolate(input$hr_t)
    req(k, !is.na(tm), !is.na(tg), tg > 0)
    cuts <- round(tm * seq_len(k - 1) / k, 2)
    hrs <- signif(exp(log(tg) * seq(0, 1, length.out = k)), 3)
    tagList(
      lapply(seq_len(k - 1), function(i)
        numericInput(paste0("rc_", i), paste("HR cut time", i), cuts[i], min = 0)),
      lapply(seq_len(k), function(i)
        numericInput(paste0("rv_", i), paste("HR, interval", i), hrs[i], min = 0.01, step = 0.05))
    )
  })
  
  spec <- reactive({
    tmax <- input$tmax
    validate(need(isTRUE(tmax > 0), "Maximum time must be positive."))
    
    use_median <- identical(input$ctrl_scale, "median")
    if (input$ctrl_type == "constant") {
      c_cuts <- numeric(0)
      if (use_median) {
        validate(need(isTRUE(input$med0 > 0), "Control median must be positive."))
        c_vals <- log(2) / input$med0
      } else {
        validate(need(isTRUE(input$h0 >= 0), "Control hazard must be non-negative."))
        c_vals <- input$h0
      }
    } else {
      k <- as.integer(input$ctrl_k)
      c_cuts <- num_inputs("cc_", k - 1)
      raw <- if (use_median) num_inputs("cm_", k) else num_inputs("ch_", k)
      validate(
        need(!anyNA(c_cuts) && !anyNA(raw), "Complete all control cut times and control values."),
        need(all(diff(c(0, c_cuts, tmax)) > 0),
             "Control cut times must increase strictly and lie between 0 and the maximum time."),
        need(if (use_median) all(raw > 0) else all(raw >= 0),
             if (use_median) "Control medians must be positive." else "Control hazards must be non-negative.")
      )
      c_vals <- if (use_median) log(2) / raw else raw
    }
    edges <- c(0, c_cuts, tmax)
    ctrl_desc <- if (length(c_vals) == 1) {
      sprintf("constant hazard %s (median %s under an exponential distribution)",
              fmt(c_vals), fmt(log(2) / c_vals))
    } else {
      paste0("piecewise constant hazard, ",
             paste(sprintf("[%s, %s): %s (median %s)", fmt(edges[-length(edges)]), fmt(edges[-1]),
                           fmt(c_vals), fmt(log(2) / c_vals)), collapse = "; "))
    }
    
    tg <- NULL
    if (input$hr_type == "pw") {
      k <- as.integer(input$hr_k)
      r_cuts <- num_inputs("rc_", k - 1)
      r_vals <- num_inputs("rv_", k)
      validate(
        need(!anyNA(r_cuts) && !anyNA(r_vals), "Complete all HR cut times and hazard ratios."),
        need(all(diff(c(0, r_cuts, tmax)) > 0),
             "HR cut times must increase strictly and lie between 0 and the maximum time."),
        need(all(r_vals > 0), "Hazard ratios must be positive.")
      )
      hr_fun <- function(t) step_val(t, r_cuts, r_vals)
      hr_cuts <- r_cuts
      r_edges <- c(0, r_cuts, tmax)
      hr_desc <- paste0("non-proportional, piecewise HR, ",
                        paste(sprintf("[%s, %s): %s", fmt(r_edges[-length(r_edges)]), fmt(r_edges[-1]),
                                      fmt(r_vals)), collapse = "; "))
    } else {
      tg <- input$hr_t
      validate(need(isTRUE(tg > 0), "Target hazard ratio must be positive."))
      if (input$hr_type == "ph") {
        hr_fun <- function(t) rep(tg, length(t))
        hr_cuts <- numeric(0)
        hr_desc <- sprintf("proportional, HR = %s", fmt(tg))
      } else if (input$preset == "delay") {
        d <- input$p_delay
        validate(need(isTRUE(d > 0), "Delay must be positive."))
        hr_fun <- function(t) ifelse(t < d, 1, tg)
        hr_cuts <- d
        hr_desc <- sprintf("non-proportional, delayed effect (HR = 1 before time %s, then %s)", fmt(d), fmt(tg))
      } else if (input$preset == "wane") {
        w <- input$p_wane
        validate(need(isTRUE(w > 0), "Waning period must be positive."))
        hr_fun <- function(t) tg + (1 - tg) * pmin(t / w, 1)
        hr_cuts <- w
        hr_desc <- sprintf("non-proportional, waning effect (HR rises linearly from %s at time 0 to 1 at time %s)",
                           fmt(tg), fmt(w))
      } else {
        cx <- input$p_cross
        early <- input$p_early
        validate(need(isTRUE(cx > 0), "Crossing time must be positive."),
                 need(isTRUE(early > 0), "Early HR must be positive."))
        hr_fun <- function(t) ifelse(t < cx, early, tg)
        hr_cuts <- cx
        hr_desc <- sprintf("non-proportional, crossing hazards (HR = %s before time %s, then %s)",
                           fmt(early), fmt(cx), fmt(tg))
      }
    }
    list(tmax = tmax, c_cuts = c_cuts, c_vals = c_vals, hr_fun = hr_fun,
         hr_cuts = hr_cuts, hr_type = input$hr_type, tg = tg,
         ctrl_desc = ctrl_desc, hr_desc = hr_desc)
  })
  
  cells <- reactive({
    s <- spec()
    make_cells(s$tmax, s$c_cuts, s$c_vals, s$hr_fun, s$hr_cuts)
  })
  
  haz_grob <- reactive(haz_grob_fn(cells(), spec()$tmax))
  
  sim <- reactive({
    cl <- cells()
    tmax <- spec()$tmax
    cens_ok <- isTRUE(input$cens >= 0) && isTRUE(input$cens <= 99)
    validate(
      need(isTRUE(input$n_arm >= 2) && isTRUE(input$n_arm <= 5000),
           "Patients per arm must be between 2 and 5000."),
      need(cens_ok, "Censoring percentage must be between 0 and 99."),
      need(isTRUE(input$cens_int > 0), "Censoring interval must be positive."),
      need(isTRUE(!is.na(input$seed)), "Enter a random seed.")
    )
    n <- as.integer(input$n_arm)
    lam <- dropout_rate(input$cens, input$cens_int)
    set.seed(input$seed)
    d <- bind_rows(
      simulate_arm(n, cl$a, cl$b, cl$hc, lam, tmax, "Control"),
      simulate_arm(n, cl$a, cl$b, cl$he, lam, tmax, "Experimental")
    )
    d$arm <- factor(d$arm, levels = names(arm_cols))
    list(data = d, lambda = lam)
  })
  
  desc_lines <- reactive({
    s <- spec()
    sm <- sim()
    ci_txt <- if (isTRUE(input$km_ci)) {
      "Pointwise 95% confidence intervals are shown (Greenwood variance, complementary log-log transformation)."
    } else {
      "Confidence intervals are not shown."
    }
    c(
      sprintf("Design: maximum time %s. Control arm: %s. Hazard ratio: %s.",
              fmt(s$tmax), s$ctrl_desc, s$hr_desc),
      sprintf(paste0("Simulation: %d patients per arm, random seed %s. Event times are drawn from the ",
                     "piecewise exponential model defined by these hazards."),
              as.integer(input$n_arm), format(input$seed, scientific = FALSE, trim = TRUE)),
      sprintf(paste0("Censoring: dropout of %s%% per %s time units in both arms (exponential, rate %.5f), ",
                     "plus administrative censoring at time %s."),
              fmt(input$cens), fmt(input$cens_int), sm$lambda, fmt(s$tmax)),
      sprintf(paste0("Display: Kaplan-Meier estimates with censoring ticks; faint lines are the true survival ",
                     "functions. %s Numbers at risk are shown every %s time units."),
              ci_txt, fmt(input$risk_step))
    )
  })
  
  output$sim_desc <- renderUI({
    div(style = "font-size:90%; color:#333; margin:6px 0 14px 0;",
        strong("Simulation parameters"),
        lapply(desc_lines(), tags$p, style = "margin:2px 0;"))
  })
  
  km_grob <- reactive({
    tmax <- spec()$tmax
    validate(need(isTRUE(input$risk_step > 0) && isTRUE(tmax / input$risk_step <= 40),
                  "Numbers-at-risk interval must be positive and give at most 40 time points."))
    km_grob_fn(sim()$data, cells(), tmax, input$risk_step, isTRUE(input$km_ci))
  })
  
  output$haz_plot <- renderPlot({
    g <- haz_grob()
    grid::grid.newpage()
    grid::grid.draw(g)
  }, res = 110)
  
  output$km_plot <- renderPlot({
    g <- km_grob()
    grid::grid.newpage()
    grid::grid.draw(g)
  }, res = 110)
  
  output$sim_tbl <- renderTable({
    s <- sim()
    s$data |>
      group_by(arm) |>
      summarise(N = n(), Events = sum(event), Censored = N - Events,
                `Censored (%)` = round(100 * Censored / N, 1), .groups = "drop") |>
      mutate(Arm = as.character(arm),
             `Dropout rate` = sprintf("%.5f", s$lambda)) |>
      select(Arm, N, Events, Censored, `Censored (%)`, `Dropout rate`)
  }, digits = 3)
  
  output$cox_txt <- renderText({
    s <- spec()
    cox_text(sim()$data, if (s$hr_type == "ph") s$tg else NULL)
  })
  
  # Saves files through the browser (base64 + Blob). Shiny's own download links do not
  # work in shinylive/webR.
  send_file <- function(name, mime, writer) {
    f <- tempfile(); on.exit(unlink(f))
    writer(f)
    session$sendCustomMessage("save_file", list(name = name, mime = mime,
      b64 = jsonlite::base64_enc(readBin(f, "raw", file.size(f)))))
  }
  on_save <- function(btn, name, mime, writer) observeEvent(input[[btn]], {
    tryCatch(send_file(name(), mime, writer),
             error = function(e) showNotification(paste("Download failed:", conditionMessage(e)),
                                                  type = "error", duration = NULL))
  })
  on_save("dl_haz", function() "hazard_plot.png", "image/png",
          function(f) save_grob(haz_grob(), f, 8, 6.5))
  on_save("dl_km", function() "km_plot.png", "image/png",
          function(f) save_km_png(km_grob(), desc_lines(), f))
  on_save("dl_dat", function() "simulated_data.csv", "text/csv", function(f) {
    d <- sim()$data
    write.csv(cbind(id = seq_len(nrow(d)), d), f, row.names = FALSE)
  })
}

shinyApp(ui, server)
