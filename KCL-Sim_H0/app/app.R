# Long-run behaviour of the 1-sided test under H0: simulated clinical trials
# Shiny app: single file (app.R)
# Display conventions follow the Type I and II Error (Power) Explorer.

library(shiny)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(plotly)
if (FALSE) library(jsonlite)  # used via jsonlite:: for downloads

`%||%` <- function(a, b) if (is.null(a)) b else a

max_trials   <- 10000  # maximum number of simulated trials
marker_limit <- 2000   # per-trial blocks up to this many trials, bin-level bars above

# ── Core calculations ─────────────────────────────────────────────────────────
# Two-arm trial, equal allocation, common SD, true effect 0 (H0), 1-sided test.
# Patient-level data are simulated for every trial. t-test: pooled-variance t
# statistic (2n - 2 df). Normal approximation: z statistic using the true SD.

sim_trials <- function(nsim, n, sd, normal = FALSE) {
  se_fac <- sqrt(2 / n)
  vapply(seq_len(nsim), function(i) {
    x <- rnorm(n, 0, sd)   # control arm
    y <- rnorm(n, 0, sd)   # treatment arm
    d <- mean(y) - mean(x)
    if (normal) d / (sd * se_fac)
    else d / (sqrt((var(x) + var(y)) / 2) * se_fac)
  }, numeric(1))
}

pval_fun <- function(stat, df, normal = FALSE) {
  if (normal) pnorm(stat, lower.tail = FALSE) else pt(stat, df, lower.tail = FALSE)
}

crit_fun <- function(alpha, df, normal = FALSE) {
  if (normal) qnorm(1 - alpha) else qt(1 - alpha, df)
}

dens_fun <- function(x, df, normal = FALSE) {
  if (normal) dnorm(x) else dt(x, df)
}

# Freedman-Diaconis width, 2 x IQR x n^(-1/3), with the IQR of the null
# distribution of the statistic, rounded to 2 significant figures.
fd_width <- function(nsim, df, normal = FALSE) {
  q <- if (normal) qnorm(0.75) else qt(0.75, df)
  signif(2 * (2 * q) * nsim^(-1 / 3), 2)
}

# Monte Carlo error line for the summary box. Under H0 with an exact test the
# number of rejections is Binomial(k, alpha), so the SE of the observed
# proportion is sqrt(alpha (1 - alpha) / k).
mc_line <- function(n_ge, k, alpha, dp) {
  if (k < 1) return(NULL)
  se <- sqrt(alpha * (1 - alpha) / k)
  lo <- max(0, alpha - 1.96 * se)
  hi <- min(1, alpha + 1.96 * se)
  z  <- (n_ge / k - alpha) / se
  sprintf(paste0(
    "Monte Carlo error: SE of the proportion = %s. If the test is exact, 95%% ",
    "of simulations of this size give a proportion between %s and %s. ",
    "Observed proportion is %+.1f SE from alpha."),
    formatC(se, format = "f", digits = dp + 1),
    formatC(lo, format = "f", digits = dp),
    formatC(hi, format = "f", digits = dp), z)
}

# Decimal places for proportions: set by the size of the simulation, so that
# 10 trials gives 1 place, 100 gives 2 and 10000 gives 4
dp_fun <- function(nsim) max(1, ceiling(log10(nsim)))

stat_label <- function(s) if (s$normal) "z-statistic" else "t-statistic"

# "Nice" major tick spacing (1, 2, 5 x 10^k) giving roughly `target` ticks
nice_dtick <- function(width, target = 10) {
  raw <- width / target
  mag <- 10^floor(log10(raw))
  m   <- raw / mag
  mag * if (m < 1.5) 1 else if (m < 3.5) 2 else if (m < 7.5) 5 else 10
}

# Colours. red_col matches alpha_col in the power explorer, without transparency.
red_col  <- "rgb(177, 69, 53)"
grey_col <- "rgb(110, 110, 110)"

# ── Help tooltip helper (same approach as the power explorer) ─────────────────
help_icon <- function(text) {
  tags$span(
    "\u24d8",
    class = "help-icon", tabindex = "0", role = "button",
    `data-toggle` = "popover", `data-bs-toggle` = "popover",
    `data-trigger` = "focus", `data-bs-trigger` = "focus",
    `data-placement` = "right", `data-bs-placement` = "right",
    `data-content` = text, `data-bs-content` = text,
    style = "cursor: pointer; color: #337ab7; margin-left: 6px; font-weight: bold;"
  )
}

labelled <- function(label, help) {
  tags$div(
    style = "display:flex; align-items:center;",
    tags$label(label, style = "margin-bottom: 0;"),
    help_icon(help)
  )
}

n_help <- paste0(
  "Patients per arm (one value). Each simulated trial draws N patients per arm ",
  "from a normal distribution with mean 0 (H0 true) and the SD below."
)
test_help <- paste0(
  "t-test: pooled-variance t statistic with 2n - 2 degrees of freedom, compared ",
  "with the central t distribution. Normal approximation: z statistic using the ",
  "true SD, compared with the standard normal distribution."
)
alpha_help <- paste0(
  "1-sided significance level. The dashed vertical line is the critical value ",
  "for this alpha. It can be changed without re-simulating."
)
ntr_help <- paste0("Number of simulated trials, 1 to ", max_trials, ".")
seed_help <- paste0(
  "Random seed. The same seed and settings reproduce the same trials. Change ",
  "the seed, then click Simulate, for a fresh set."
)
units_help <- paste0(
  "Bar height as a count of trials or as a percentage of all simulated trials. ",
  "The default is counts below 20 trials and percentages from 20 trials upwards; ",
  "it is switched when a simulation crosses that threshold."
)
bw_help <- paste0(
  "Bar width in units of the test statistic. Leave empty for the ",
  "Freedman-Diaconis default, 2 x IQR x trials^(-1/3), with the IQR taken from ",
  "the null distribution and rounded to 2 significant figures."
)

# ── UI ────────────────────────────────────────────────────────────────────────

# ---- footer for exported PNG and PDF files: notice, app address, UTC date and time ----
APP_URL <- "https://kinnersleyn.github.io/kcl-apps/Sim_H0/"
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
  tags$head(
    # Tell the server when plotly.js and crosstalk are loaded, so plots are not drawn too early
    tags$script(HTML(paste0("$(document).on('shiny:connected', function(){",
      "var t=setInterval(function(){",
      "if (typeof Plotly!=='undefined' && typeof crosstalk!=='undefined') {",
      "clearInterval(t); Shiny.setInputValue('plotly_ready', true);}}, 100);});"))),
    tags$script(HTML(paste0("$(document).on('shiny:connected', function(){",
                            "Shiny.addCustomMessageHandler('save_file', function(m){",
                            "var bin=atob(m.b64), n=bin.length, bytes=new Uint8Array(n);",
                            "for (var i=0;i<n;i++) bytes[i]=bin.charCodeAt(i);",
                            "var url=URL.createObjectURL(new Blob([bytes],{type:m.mime}));",
                            "var a=document.createElement('a'); a.href=url; a.download=m.name;",
                            "document.body.appendChild(a); a.click(); document.body.removeChild(a);",
                            "setTimeout(function(){URL.revokeObjectURL(url);},10000);});});"))),
    tags$script(HTML(paste0(
      "$(function(){ $('.help-icon').popover({container: 'body'}); });\n",
      "function downloadPlotPng(id, fname, w, h) {\n",
      "  var el = document.getElementById(id);\n",
      "  var gd = el && (el.classList.contains('js-plotly-plot') ? el :\n",
      "                  el.querySelector('.js-plotly-plot'));\n",
      "  if (!gd || typeof Plotly === 'undefined') {\n",
      "    alert('The plot is not ready yet, so it cannot be downloaded.');\n",
      "    return;\n",
      "  }\n",
      "  var url = 'https://kinnersleyn.github.io/kcl-apps/Sim_H0/';\n",
      "  var d = new Date(), p2 = function(n){ return (n < 10 ? '0' : '') + n; };\n",
      "  var stamp = d.getUTCFullYear() + '-' + p2(d.getUTCMonth() + 1) + '-' + p2(d.getUTCDate()) + ' ' + p2(d.getUTCHours()) + ':' + p2(d.getUTCMinutes()) + ':' + p2(d.getUTCSeconds()) + ' UTC';\n",
      "  var txt = 'For educational purposes only, not for decision-making. ' + url + ' Generated ' + stamp;\n",
      "  var m = gd._fullLayout.margin, old = gd.layout.annotations || [], extra = 22;\n",
      "  var lg = gd.layout.legend, ly = lg && lg.y;\n",
      "  var size = Math.min(11, 0.95 * w / (txt.length * 0.56));\n",
      "  var ann = {xref: 'paper', yref: 'paper', x: 0.5, y: 0, xanchor: 'center', yanchor: 'bottom', showarrow: false,\n",
      "             xshift: (m.r - m.l) / 2, yshift: -(m.b + extra - 5), text: txt, font: {size: size, color: '#555555'}};\n",
      "  var on = {'margin.b': m.b + extra, annotations: old.concat([ann])};\n",
      "  var off = {'margin.b': m.b, annotations: old};\n",
      "  if (lg && lg.yref === 'container' && typeof ly === 'number') { on['legend.y'] = ly + extra / h; off['legend.y'] = ly; }\n",
      "  Plotly.relayout(gd, on).then(function() {\n",
      "    return Plotly.downloadImage(gd, {format: 'png', filename: fname, width: w, height: h, scale: 2});\n",
      "  }).then(function() { return Plotly.relayout(gd, off); });\n",
      "}"
    )))
  ),
  # Load plotly.js and its dependencies (crosstalk, typedarray, css) with the page. Under shinylive/webR it otherwise loads on demand and the
  # first plot can render before it is ready ("Plotly is not defined").
  htmltools::attachDependencies(
    tags$div(style = "display:none"),
    htmltools::findDependencies(plotly::as_widget(plotly::plotly_build(plotly::plot_ly(x = 1, y = 1))))),
  uiOutput("dyn_css"),
  uiOutput("dyn_plot_css"),
  titlePanel("Long-run behaviour of the 1-sided test under H0: simulated trials"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      tags$b("Display"),
      fluidRow(
        column(6, selectInput(
          "text_size", "Option text size",
          choices = c("100%" = 100, "90%" = 90, "80%" = 80, "70%" = 70),
          selected = 100
        )),
        column(6, sliderInput("plot_scale", "Plot size (%)",
                              min = 40, max = 100, value = 100, step = 5))
      ),
      checkboxInput("gridlines", "Show gridlines", value = TRUE),
      
      tags$hr(),
      tags$b("Trial design"),
      labelled("Patients per arm:", n_help),
      numericInput("n", label = NULL, value = 175, min = 2, step = 1),
      
      labelled("Test:", test_help),
      radioButtons(
        "test", label = NULL,
        choices = c("t-test" = "t", "Normal approximation" = "normal"),
        selected = "t", inline = TRUE
      ),
      
      fluidRow(
        column(6, numericInput("sd", "SD", value = 230, min = 0.0001)),
        column(6, tagList(
          labelled("alpha (1-sided)", alpha_help),
          numericInput("alpha", label = NULL, value = 0.025, min = 0.0001,
                       max = 0.5, step = 0.005)
        ))
      ),
      
      tags$hr(),
      tags$b("Simulation"),
      fluidRow(
        column(6, tagList(
          labelled("Number of trials", ntr_help),
          numericInput("n_sims", label = NULL, value = 10, min = 1,
                       max = max_trials, step = 1)
        )),
        column(6, tagList(
          labelled("Seed", seed_help),
          numericInput("seed", label = NULL, value = 123, step = 1)
        ))
      ),
      actionButton("simulate", "Simulate", class = "btn-primary"),
      checkboxInput("animate", "Animate (add trials progressively)", value = TRUE),
      tags$small(tags$em(paste0(
        "Patients per arm, test, SD, number of trials and seed are applied when ",
        "Simulate is clicked. alpha and the display options update immediately."
      ))),
      
      tags$hr(),
      tags$b("Bars"),
      labelled("Bar height:", units_help),
      radioButtons(
        "units", label = NULL,
        choices = c("Counts" = "count", "% of trials" = "pct"),
        selected = "count", inline = TRUE
      ),
      labelled("Bar width override:", bw_help),
      numericInput("bw", label = NULL, value = NA, min = 0.001, step = 0.05),
      uiOutput("bw_used"),
      checkboxInput("show_curve", "Show expected counts under H0 (curve)",
                    value = TRUE),
      
    ),
    mainPanel(
      width = 9,
      tags$div(
        style = "margin-bottom:10px; display:flex; flex-wrap:wrap; gap:6px; align-items:center;",
        tags$b("Save:"),
        tags$button(
          type = "button", class = "btn btn-default",
          onclick = "downloadPlotPng('plot', 'simulated_trials_under_H0', 1400, 800);",
          icon("download"), " Plot (PNG)"
        ),
        actionButton("dl_csv", "Values (CSV)", icon = icon("download")),
        actionButton("dl_pdf", "All results (PDF)", icon = icon("download"))
      ),
      uiOutput("summary_box"),
      plotlyOutput("plot", height = "620px")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  rv <- reactiveValues(sim = NULL, k = 0, running = FALSE, side = NULL)
  
  # ── Simulate ────────────────────────────────────────────────────────────────
  observeEvent(input$simulate, ignoreNULL = FALSE, {
    n    <- input$n
    sd   <- input$sd
    nsim <- input$n_sims
    seed <- input$seed
    
    msg <- NULL
    if (anyNA(c(n, sd, nsim, seed)) || length(c(n, sd, nsim, seed)) < 4) {
      msg <- "Patients per arm, SD, number of trials and seed must all be numbers."
    } else if (n < 2) {
      msg <- "Patients per arm must be at least 2."
    } else if (sd <= 0) {
      msg <- "SD must be positive."
    } else if (nsim < 1 || nsim > max_trials) {
      msg <- paste0("Number of trials must be between 1 and ", max_trials, ".")
    }
    if (!is.null(msg)) {
      showNotification(msg, type = "error")
      return()
    }
    
    n      <- as.integer(round(n))
    nsim   <- as.integer(round(nsim))
    normal <- input$test == "normal"
    df     <- 2L * n - 2L
    
    set.seed(as.integer(seed))
    ts <- sim_trials(nsim, n, sd, normal)
    
    rv$sim <- list(t = ts, p = pval_fun(ts, df, normal), n = n, sd = sd,
                   nsim = nsim, seed = as.integer(seed), normal = normal, df = df)
    
    # Bar height default: counts below 20 trials, percentages from 20 upwards.
    # Switched only when a simulation crosses the threshold.
    side <- nsim >= 20
    if (!identical(side, rv$side)) {
      updateRadioButtons(session, "units", selected = if (side) "pct" else "count")
      rv$side <- side
    }
    
    if (isTRUE(input$animate)) {
      rv$k <- 0
      rv$running <- TRUE
    } else {
      rv$k <- nsim
      rv$running <- FALSE
    }
  })
  
  # ── Animation: add trials in steps until all are shown ──────────────────────
  observe({
    req(rv$running)
    s <- isolate(rv$sim)
    frames <- min(s$nsim, 60)
    step   <- ceiling(s$nsim / frames)
    invalidateLater(if (frames <= 20) 300 else 100, session)
    isolate({
      k <- min(rv$k + step, s$nsim)
      rv$k <- k
      if (k >= s$nsim) rv$running <- FALSE
    })
  })
  
  # ── Bar width ───────────────────────────────────────────────────────────────
  bw_override <- reactive({
    bw <- input$bw
    if (!is.null(bw) && length(bw) == 1 && !is.na(bw) && bw > 0) bw else NULL
  })
  
  bar_w <- reactive({
    s <- rv$sim
    req(s)
    bw_override() %||% fd_width(s$nsim, s$df, s$normal)
  })
  
  output$bw_used <- renderUI({
    s <- rv$sim
    req(s)
    tags$small(sprintf(
      "Bar width in use: %s (%s).", format(bar_w()),
      if (is.null(bw_override())) "Freedman-Diaconis default" else "user override"
    ))
  })
  
  # ── Display size ────────────────────────────────────────────────────────────
  output$dyn_css <- renderUI({
    s <- suppressWarnings(as.numeric(input$text_size)) / 100
    if (length(s) != 1 || is.na(s)) s <- 1
    px <- function(x) paste0(round(x, 1), "px")
    tags$style(HTML(paste0(
      ".well, .well label, .well .control-label, .well .checkbox, .well .radio, ",
      ".well .btn, .well .form-control, .well .selectize-input, ",
      ".well .selectize-dropdown, .well .irs { font-size: ", px(14 * s), "; }\n",
      ".well .form-control { height: ", px(34 * s), "; padding: ",
      px(6 * s), " ", px(12 * s), "; }\n",
      ".well .form-group { margin-bottom: ", px(6 * s), "; }\n",
      ".well .checkbox, .well .radio { margin-top: ", px(2 * s),
      "; margin-bottom: ", px(2 * s), "; }\n",
      ".well hr { margin-top: ", px(6 * s), "; margin-bottom: ", px(6 * s), "; }\n",
      ".well { padding: ", px(10 * s), "; }"
    )))
  })
  
  plot_scale_r <- debounce(
    reactive(if (is.null(input$plot_scale)) 100 else input$plot_scale), 400
  )
  plot_f <- reactive(max(plot_scale_r() / 100, 0.7))
  
  output$dyn_plot_css <- renderUI({
    s <- plot_scale_r() / 100
    tagList(
      tags$style(HTML(paste0("#plot { height: ", round(620 * s), "px !important; }"))),
      tags$script(HTML(paste0(
        "setTimeout(function(){ window.dispatchEvent(new Event('resize')); }, 100);",
        "setTimeout(function(){ window.dispatchEvent(new Event('resize')); }, 500);"
      )))
    )
  })
  
  alpha_ok <- reactive({
    !is.null(input$alpha) && !is.na(input$alpha) &&
      input$alpha > 0 && input$alpha < 0.5
  })
  
  # ── Summary box ─────────────────────────────────────────────────────────────
  output$summary_box <- renderUI({
    s <- rv$sim
    req(s, alpha_ok())
    k     <- min(rv$k, s$nsim)
    alpha <- input$alpha
    crit  <- crit_fun(alpha, s$df, s$normal)
    slab  <- stat_label(s)
    n_ge  <- sum(s$t[seq_len(k)] >= crit)
    
    line1 <- sprintf("Trials simulated: %d of %d%s.", k, s$nsim,
                     if (isTRUE(rv$running)) " (animation running)" else "")
    line2 <- if (k == 0) {
      sprintf("Trials with %s \u2265 critical value (%s): none simulated yet.",
              slab, formatC(crit, format = "f", digits = 3))
    } else {
      sprintf(
        "Trials with %s \u2265 critical value (%s): frequency %d, proportion %s (%s%%).",
        slab, formatC(crit, format = "f", digits = 3), n_ge,
        formatC(n_ge / k, format = "f", digits = dp_fun(s$nsim)),
        formatC(100 * n_ge / k, format = "f", digits = max(0, dp_fun(s$nsim) - 2))
      )
    }
    line3 <- sprintf(
      "Expected under H0 (1-sided alpha = %s): proportion %s, about %s trials.",
      format(alpha, scientific = FALSE), format(alpha, scientific = FALSE),
      formatC(alpha * k, format = "f", digits = 1)
    )
    
    tags$div(
      style = paste0("margin: 10px 0; padding: 8px 12px; border: 1px solid #ccc; ",
                     "border-radius: 4px; background: #fafafa;"),
      line1, tags$br(), tags$b(line2), tags$br(), line3,
      if (k > 0) tagList(tags$br(), mc_line(n_ge, k, alpha, dp_fun(s$nsim)))
    )
  })
  
  # ── Plot ────────────────────────────────────────────────────────────────────
  output$plot <- renderPlotly({
    req(isTRUE(input$plotly_ready))
    s <- rv$sim
    req(s)
    validate(need(alpha_ok(), "alpha (1-sided) must be between 0 and 0.5."))
    k     <- min(rv$k, s$nsim)   # guard against a stale k after re-simulating
    alpha <- input$alpha
    f     <- plot_f()
    w     <- bar_w()
    crit  <- crit_fun(alpha, s$df, s$normal)
    pct   <- identical(input$units, "pct")
    u     <- if (pct) 100 / s$nsim else 1       # height of one trial
    slab  <- stat_label(s)
    f2    <- function(x) formatC(x, format = "f", digits = 2)
    f3    <- function(x) formatC(x, format = "f", digits = 3)
    fp    <- function(x) formatC(x, format = "g", digits = 3)
    
    # Axis ranges come from the full set of trials, so nothing moves while animating
    q    <- if (s$normal) qnorm(0.9995) else qt(0.9995, s$df)
    m    <- max(abs(c(s$t, crit, q))) + w
    xr   <- c(-m, m)
    xs   <- seq(-m, m, length.out = 501)
    ycv  <- (if (pct) 100 else s$nsim) * w * dens_fun(xs, s$df, s$normal)
    ymax <- 1.12 * max(max(table(round(s$t / w))) * u, ycv)
    
    # Trials shown so far, assigned to bins centred on multiples of w
    d <- tibble(t = s$t[seq_len(k)], p = s$p[seq_len(k)]) |>
      mutate(bin = round(t / w), ge = t >= crit)
    
    if (s$nsim <= marker_limit) {
      # One block per trial, stacked within its bin (trials below the critical
      # value first, so those at or above it sit on top)
      d <- d |>
        arrange(bin, ge, t) |>
        group_by(bin) |>
        mutate(pos = row_number()) |>
        ungroup() |>
        mutate(x = bin * w, base = (pos - 1) * u,
               txt = sprintf("%s: %s<br>1-sided p-value: %s",
                             slab, f3(t), fp(p)))
      mk <- function(df) list(x = df$x, y = rep(u, nrow(df)), base = df$base,
                              text = df$txt)
      lo <- mk(filter(d, !ge))
      hi <- mk(filter(d, ge))
      # White outlines only when blocks are tall enough to show them; with
      # many trials each block is a thin sliver and the outline hides the fill
      lw <- if (s$nsim <= 100) 1 else if (s$nsim <= 400) 0.5 else 0
    } else {
      # Bin-level bars
      b <- d |>
        group_by(bin) |>
        summarise(n_lo = sum(!ge), n_hi = sum(ge), .groups = "drop") |>
        mutate(n = n_lo + n_hi, x = bin * w,
               txt = sprintf(
                 "Bin: %s to %s<br>Count: %d<br>1-sided p-value at bin centre: %s",
                 f3(x - w / 2), f3(x + w / 2), n,
                 fp(pval_fun(x, s$df, s$normal))))
      b_lo <- filter(b, n_lo > 0)
      b_hi <- filter(b, n_hi > 0)
      lo <- list(x = b_lo$x, y = b_lo$n_lo * u, base = rep(0, nrow(b_lo)),
                 text = b_lo$txt)
      hi <- list(x = b_hi$x, y = b_hi$n_hi * u, base = b_hi$n_lo * u,
                 text = b_hi$txt)
      lw <- 0.5
    }
    
    add_blocks <- function(p, z, name, color) {
      empty <- length(z$x) == 0
      if (empty) {
        # Zero-height placeholder keeps the legend entry without NA values,
        # which make plotly warn "Ignoring 1 observations"
        z <- list(x = 0, y = 0, base = 0, text = "")
      }
      add_trace(
        p, type = "bar", x = z$x, y = z$y, base = z$base, width = w,
        name = name, text = z$text, textposition = "none",
        hoverinfo = if (empty) "skip" else "text",
        marker = list(color = color, line = list(color = "white", width = lw)),
        hoverlabel = list(bgcolor = color, font = list(color = "white"))
      )
    }
    
    p <- plot_ly() |>
      add_blocks(lo, paste0(slab, " < critical value"), grey_col) |>
      add_blocks(hi, paste0(slab, " \u2265 critical value"), red_col)
    
    if (isTRUE(input$show_curve)) {
      p <- p |>
        add_trace(
          type = "scatter", mode = "lines", x = xs, y = ycv,
          name = if (s$normal) "Expected under H0 (standard normal)"
          else paste0("Expected under H0 (central t, df = ", s$df, ")"),
          line = list(color = "black", width = 2), hoverinfo = "skip"
        )
    }
    
    ttl <- paste0(
      "Simulated ", slab, "s under H0: ", s$nsim, " trials, N = ", s$n,
      " per arm<br><sup>SD = ", s$sd, ", 1-sided alpha = ", alpha,
      if (s$normal) ", normal approximation" else ", t-test",
      ", seed = ", s$seed, ", bar width = ", format(w), "</sup>"
    )
    
    p |>
      layout(
        font = list(size = 12 * f),
        title = list(text = ttl, font = list(size = 14 * f)),
        hovermode = "closest",
        barmode = "overlay", bargap = 0,
        xaxis = list(
          title = list(text = slab, standoff = 10),
          range = xr, zeroline = FALSE, automargin = TRUE,
          showline = TRUE, linecolor = "black", linewidth = 1, ticks = "outside",
          showgrid = isTRUE(input$gridlines), gridcolor = "#D0D0D0",
          tick0 = 0, dtick = nice_dtick(diff(xr)),
          minor = list(dtick = nice_dtick(diff(xr)) / 2,
                       showgrid = isTRUE(input$gridlines), gridcolor = "#EDEDED")
        ),
        yaxis = list(
          title = list(text = if (pct) "% of trials" else "Number of trials",
                       standoff = 10),
          range = c(0, ymax), zeroline = FALSE, automargin = TRUE,
          showline = TRUE, linecolor = "black", linewidth = 1, ticks = "outside",
          showgrid = isTRUE(input$gridlines), gridcolor = "#D0D0D0",
          tick0 = 0,
          dtick = if (pct) nice_dtick(ymax, 8) else max(1, nice_dtick(ymax, 8)),
          ticksuffix = if (pct) "%" else ""
        ),
        shapes = list(list(
          type = "line", xref = "x", x0 = crit, x1 = crit,
          yref = "paper", y0 = 0, y1 = 1,
          line = list(dash = "dash", color = "black", width = 1)
        )),
        annotations = list(list(
          xref = "x", yref = "paper", x = crit, xanchor = "center",
          y = 1, yanchor = "bottom", align = "center", showarrow = FALSE,
          text = paste0("Critical value = ", f2(crit),
                        "<br>1-sided alpha = ", alpha),
          font = list(size = 11 * f)
        )),
        legend = list(orientation = "h", x = 0.5, xanchor = "center",
                      yref = "container", y = 0.01, yanchor = "bottom"),
        margin = list(t = 115 * f, b = 135 * f, r = 30)
      ) |>
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = list("toImage"),
        toImageButtonOptions = list(format = "png",
                                    filename = "simulated_trials_under_H0",
                                    width = 1400, height = 800, scale = 2)
      )
  })
  
  # ── CSV download ────────────────────────────────────────────────────────────
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
  on_save("dl_csv",
          function() paste0("simulated_trials_under_H0_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          "text/csv", function(f) {
    s <- rv$sim
    if (is.null(s)) stop("Run a simulation first.")
    if (!alpha_ok()) stop("alpha (1-sided) must be between 0 and 0.5.")
    alpha <- input$alpha
    crit  <- crit_fun(alpha, s$df, s$normal)
    tibble(
      Trial            = seq_along(s$t),
      Statistic        = s$t,
      P_value_1sided   = s$p,
      Ge_critical      = s$t >= crit,
      Critical_value   = crit,
      Alpha_1sided     = alpha,
      N_per_arm        = s$n,
      SD               = s$sd,
      Test             = if (s$normal) "normal approximation" else "t-test",
      Seed             = s$seed
    ) |>
      write.csv(f, row.names = FALSE)
  })

  text_page <- function(lines) {
    lines <- unlist(lapply(lines, function(l) if (nzchar(l)) strwrap(l, width = 105) else ""))
    grid::grid.newpage()
    grid::grid.text(paste(lines, collapse = "\n"), x = grid::unit(0.5, "in"),
                    y = grid::unit(1, "npc") - grid::unit(0.5, "in"), just = c("left", "top"),
                    gp = grid::gpar(fontsize = 10))
    export_footer()
  }

  # Static version of the plot for the PDF (plotly cannot be written to PDF from R)
  sim_gg <- function() {
    s <- rv$sim
    if (is.null(s)) stop("Run a simulation first.")
    if (!alpha_ok()) stop("alpha (1-sided) must be between 0 and 0.5.")
    alpha <- input$alpha
    w     <- bar_w()
    crit  <- crit_fun(alpha, s$df, s$normal)
    pct   <- identical(input$units, "pct")
    u     <- if (pct) 100 / s$nsim else 1
    slab  <- stat_label(s)
    q     <- if (s$normal) qnorm(0.9995) else qt(0.9995, s$df)
    m     <- max(abs(c(s$t, crit, q))) + w
    xs    <- seq(-m, m, length.out = 501)
    ycv   <- (if (pct) 100 else s$nsim) * w * dens_fun(xs, s$df, s$normal)
    d <- tibble(t = s$t) |>
      mutate(bin = round(t / w), ge = t >= crit) |>
      count(bin, ge) |>
      mutate(x = bin * w, y = n * u,
             grp = factor(ifelse(ge, paste0(slab, " >= critical value"),
                                 paste0(slab, " < critical value")),
                          levels = c(paste0(slab, " < critical value"),
                                     paste0(slab, " >= critical value"))))
    p <- ggplot(d, aes(x, y, fill = grp)) +
      geom_col(width = w, colour = NA) +
      scale_fill_manual(values = c("#6E6E6E", "#B14535"), name = NULL, drop = FALSE) +
      geom_vline(xintercept = crit, linetype = "dashed") +
      annotate("text", x = crit, y = Inf, vjust = 1.4, hjust = -0.05, size = 3.3,
               label = paste0("Critical value = ", formatC(crit, format = "f", digits = 2))) +
      coord_cartesian(xlim = c(-m, m)) +
      labs(x = slab, y = if (pct) "% of trials" else "Number of trials",
           title = paste0("Simulated ", slab, "s under H0: ", s$nsim, " trials, N = ", s$n, " per arm"),
           subtitle = paste0("SD = ", s$sd, ", 1-sided alpha = ", alpha,
                             if (s$normal) ", normal approximation" else ", t-test",
                             ", seed = ", s$seed, ", bar width = ", format(w))) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
    if (isTRUE(input$show_curve)) {
      p <- p + geom_line(data = data.frame(x = xs, y = ycv), aes(x, y), inherit.aes = FALSE, linewidth = 0.9)
    }
    p
  }

  sim_lines <- function() {
    s <- rv$sim
    alpha <- input$alpha
    k     <- s$nsim
    crit  <- crit_fun(alpha, s$df, s$normal)
    slab  <- stat_label(s)
    n_ge  <- sum(s$t >= crit)
    c(sprintf("Long-run behaviour of the 1-sided test under H0: %d simulated trials.", k),
      sprintf("N = %d per arm, SD = %s, %s, 1-sided alpha = %s, seed = %d.", s$n, format(s$sd),
              if (s$normal) "normal approximation" else "t-test", format(alpha, scientific = FALSE), s$seed),
      "",
      sprintf("Trials with %s >= critical value (%s): frequency %d, proportion %s (%s%%).",
              slab, formatC(crit, format = "f", digits = 3), n_ge,
              formatC(n_ge / k, format = "f", digits = dp_fun(k)),
              formatC(100 * n_ge / k, format = "f", digits = max(0, dp_fun(k) - 2))),
      sprintf("Expected under H0 (1-sided alpha = %s): proportion %s, about %s trials.",
              format(alpha, scientific = FALSE), format(alpha, scientific = FALSE),
              formatC(alpha * k, format = "f", digits = 1)),
      mc_line(n_ge, k, alpha, dp_fun(k)))
  }
  on_save("dl_pdf", function() "simulated_trials_under_H0.pdf", "application/pdf", function(f) {
    g <- sim_gg(); l <- sim_lines()
    pdf(f, width = 10, height = 7); on.exit(dev.off())
    print(g, vp = footer_vp()); export_footer(); text_page(l)
  })
}

# ── Run ───────────────────────────────────────────────────────────────────────

shinyApp(ui = ui, server = server)