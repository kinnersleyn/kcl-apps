# Type I and II Error (Power) Explorer: 2-arm trial, continuous endpoint, 1-sided test
# Shiny app: single file (app.R)
# Adapted from the Type I and II Error Demo QMD (R-CONIS MDD vs target concept)

library(shiny)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(plotly)
if (FALSE) library(jsonlite)  # used via jsonlite:: for downloads

`%||%` <- function(a, b) if (is.null(a)) b else a

# ── Core calculations ─────────────────────────────────────────────────────────
# Two-sample comparison of means, equal allocation, common SD, 1-sided test.
# The t-test power uses the non-central t distribution, which is what rpact's
# getPowerMeans(normalApproximation = FALSE) uses for a fixed (kMax = 1) design.

std_err <- function(n, sd) sd * sqrt(2 / n)

# Power at a given true effect size. ncp is clipped to +/-37 because pt() with
# a non-central parameter is not reliable beyond about 37.6 (power there is
# 0 or 1 to machine precision anyway).
power_fun <- function(eff, n, sd, alpha, normal = FALSE) {
  ncp <- eff / std_err(n, sd)
  if (normal) return(pnorm(ncp - qnorm(1 - alpha)))
  df    <- 2 * n - 2
  tcrit <- qt(1 - alpha, df)
  suppressWarnings(
    pt(tcrit, df, ncp = pmax(pmin(ncp, 37), -37), lower.tail = FALSE)
  )
}

# 1-sided p-value for an observed effect size
pval_fun <- function(eff, n, sd, normal = FALSE) {
  z <- eff / std_err(n, sd)
  if (normal) pnorm(z, lower.tail = FALSE) else pt(z, 2 * n - 2, lower.tail = FALSE)
}

# Minimum detectable difference: smallest observed effect with 1-sided p <= alpha
mdd_fun <- function(n, sd, alpha, normal = FALSE) {
  q <- if (normal) qnorm(1 - alpha) else qt(1 - alpha, 2 * n - 2)
  q * std_err(n, sd)
}

# True effect size at which power equals p (vectorised over p)
effect_at_power <- function(p, n, sd, alpha, normal = FALSE) {
  se <- std_err(n, sd)
  if (normal) return((qnorm(1 - alpha) + qnorm(p)) * se)
  df    <- 2 * n - 2
  tcrit <- qt(1 - alpha, df)
  vapply(p, function(pi) {
    f <- function(ncp) {
      suppressWarnings(pt(tcrit, df, ncp = ncp, lower.tail = FALSE)) - pi
    }
    r <- tryCatch(uniroot(f, c(-30, 37), tol = 1e-10)$root,
                  error = function(e) NA_real_)
    r * se
  }, numeric(1))
}

# "Nice" major tick spacing (1, 2, 5 x 10^k) giving roughly `target` ticks
nice_dtick <- function(width, target = 10) {
  raw <- width / target
  mag <- 10^floor(log10(raw))
  m   <- raw / mag
  mag * if (m < 1.5) 1 else if (m < 3.5) 2 else if (m < 7.5) 5 else 10
}

# Read plotly events from the power plot (source "tm"). The server observers
# start before the plot has rendered and registered its events, which makes
# plotly warn that the event is not registered; the events work once the plot
# exists, so the warning is suppressed here.
ev_data <- function(event, priority = "input") {
  suppressWarnings(event_data(event, source = "tm", priority = priority))
}

# Parse a comma/space separated list of numbers
parse_nums <- function(x) {
  if (is.null(x)) return(numeric(0))
  v <- suppressWarnings(as.numeric(unlist(strsplit(trimws(x), "[,;[:space:]]+"))))
  v[!is.na(v)]
}

# ── Distribution plot helpers ─────────────────────────────────────────────────
# Colours for plot 2. alpha_col is the dark red sampled from the reference site
# (about rgb(177, 69, 53)); change it here if you want a different shade.
alpha_col <- "rgba(177, 69, 53, 0.55)"   # dark red: alpha (Type I error) area
beta_col  <- "rgba(27, 94, 32, 0.55)"    # dark green: beta (Type II error) area
# The last number in alpha_col and beta_col is the opacity (1 = solid).
grey_col  <- "#6E6E6E"                   # null (solid) and alternative (dotted)

# Densities on the effect size scale. Null: central t (or normal) centred at 0.
# Alternative: non-central t (or normal) centred at eff. The t versions match
# the power calculation used in plot 1. ncp is clipped as in power_fun().
dens_null <- function(x, se, df, normal = FALSE) {
  if (normal) dnorm(x / se) / se else dt(x / se, df) / se
}
dens_alt <- function(x, eff, se, df, normal = FALSE) {
  if (normal) return(dnorm((x - eff) / se) / se)
  suppressWarnings(dt(x / se, df, ncp = max(min(eff / se, 37), -37))) / se
}

# Scalar summaries for a selection (effect size and N)
dist_stats <- function(s, d) {
  se   <- std_err(s$n, d$sd)
  df   <- 2 * s$n - 2
  crit <- if (d$normal) qnorm(1 - d$alpha) * se else qt(1 - d$alpha, df) * se
  pow  <- power_fun(s$eff, s$n, d$sd, d$alpha, d$normal)
  list(se = se, df = df, crit = crit, power = pow, beta = 1 - pow)
}

src_label <- function(src) {
  switch(src,
         locked  = "locked",
         hover   = "last hovered",
         default = "default: effect size giving the target power, first N")
}

# Everything plot 2 needs for one selection: trace data, names, shapes,
# annotations and title. Used for both the initial render and proxy updates.
dist_calc <- function(s, d, rng, f = 1) {
  st <- dist_stats(s, d)
  xs <- seq(rng$x[1], rng$x[2], length.out = 1001)
  yn <- dens_null(xs, st$se, st$df, d$normal)
  ya <- dens_alt(xs, s$eff, st$se, st$df, d$normal)
  
  xa  <- c(st$crit, xs[xs > st$crit])
  yar <- dens_null(xa, st$se, st$df, d$normal)
  xb  <- c(xs[xs < st$crit], st$crit)
  ybr <- dens_alt(xb, s$eff, st$se, st$df, d$normal)
  
  f2 <- function(x) formatC(x, format = "f", digits = 2)
  fmt1 <- function(x) formatC(x, format = "f", digits = 1)
  eff_txt <- f2(s$eff)
  # Beta rounded first, so the displayed beta and power add to exactly 100%
  b_r <- round(100 * st$beta, 1)
  p_r <- 100 - b_r
  
  shapes <- list(list(
    type = "line", xref = "x", x0 = st$crit, x1 = st$crit,
    yref = "paper", y0 = 0, y1 = 1,
    line = list(dash = "dash", color = "black", width = 1)
  ))
  annots <- list(
    list(xref = "x", yref = "y", x = 0,
         y = dens_null(0, st$se, st$df, d$normal), yanchor = "bottom",
         text = "H<sub>0</sub>", showarrow = FALSE, font = list(size = 14 * f)),
    list(xref = "x", yref = "y", x = s$eff, y = max(ya), yanchor = "bottom",
         text = "H<sub>1</sub>", showarrow = FALSE, font = list(size = 14 * f)),
    # Placed just above the plot area, on top of the critical value line, so
    # it cannot overlap the H0 / H1 labels inside the plot.
    list(xref = "x", yref = "paper", x = st$crit, xanchor = "center",
         y = 1, yanchor = "bottom", align = "center", showarrow = FALSE,
         text = paste0("Critical value = ", f2(st$crit),
                       "<br>Power = ", fmt1(p_r), "%"),
         font = list(size = 11 * f))
  )
  
  list(
    xs = xs, yn = yn, ya = ya, xa = xa, yar = yar, xb = xb, ybr = ybr,
    names = c(
      "Null (H0): effect size 0",
      paste0("Alternative (H1): effect size ", eff_txt),
      paste0("alpha (Type I error) = ", format(100 * d$alpha, digits = 3), "%"),
      paste0("beta (Type II error) = ", fmt1(b_r), "%, power = ", fmt1(p_r), "%")
    ),
    shapes = shapes, annots = annots,
    title = paste0(
      "H1 effect size = ", eff_txt, ", N = ", s$n, " per arm (", src_label(s$src), ")",
      "<br><sup>1-sided alpha = ", d$alpha, ", SD = ", d$sd,
      if (d$normal) ", normal approximation" else ", t-test", "</sup>"
    )
  )
}

# ── Help tooltip helper (same approach as the Bayesian Binomial app) ──────────
help_icon <- function(text) {
  # Click-to-show popover (Bootstrap), dismissed by clicking elsewhere.
  # Both data-* and data-bs-* attributes are set so it works with Bootstrap 3 or 5.
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
  "Patients per arm. Enter one number (e.g. 175) or two separated by a comma ",
  "(e.g. 175, 113). The second value is overlaid with dotted lines. ",
  "Default vertical reference lines and markers always refer to the first value."
)
test_help <- paste0(
  "t-test uses the non-central t distribution with 2n - 2 degrees of freedom ",
  "(as rpact with normalApproximation = FALSE). Normal approximation uses the ",
  "standard normal distribution."
)
alpha_help <- "1-sided significance level. A 1-sided 0.025 equals a 2-sided 0.05."
target_help <- paste0(
  "Target power as a percentage (e.g. enter 90 for 90%). With markers ",
  "switched on, a marker is placed at the effect size that gives this power ",
  "for the first N."
)
href_help <- paste0(
  "Horizontal dashed reference lines, as percentages on the primary y-axis, ",
  "comma-separated (e.g. 2.5, 50, 90). The same heights read as 1-sided ",
  "p-values of 0.025, 0.5 and 0.9 on the secondary axis."
)
vauto_help <- paste0(
  "Vertical dashed lines placed at the true effect size where the power curve ",
  "for the first N reaches each of these power levels (%)."
)
vextra_help <- "Extra vertical dashed lines at these effect sizes, comma-separated."
plots_help <- paste0(
  "Choose which plots to display. At least one must stay selected. The ",
  "distributions plot follows the effect size and N of the curve you hover over ",
  "in the power plot. Click the power plot to lock the distributions at that ",
  "point, and click again to release."
)
marker_help <- paste0(
  "For the first N: power at effect size 0 (equals alpha), the minimum ",
  "detectable difference (MDD; observed effect with 1-sided p = alpha) shown ",
  "on both curves, and the effect size giving the target power."
)

# ── UI ────────────────────────────────────────────────────────────────────────

# ---- footer for exported PNG and PDF files: notice, app address, UTC date and time ----
APP_URL <- "https://kinnersleyn.github.io/kcl-apps/T12E/"
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
    tags$link(rel = "icon", type = "image/jpeg", href = "icon.jpg"),
    tags$link(rel = "apple-touch-icon", href = "icon.jpg"),
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
      "  var url = 'https://kinnersleyn.github.io/kcl-apps/T12E/';\n",
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
  titlePanel("Type I and II Error (Power) Explorer: 2-arm trial, 1-sided test"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
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
      
      labelled("Plots shown:", plots_help),
      checkboxGroupInput(
        "show_plots", label = NULL,
        choices = c("Power and p-value plot" = "p1",
                    "Null and alternative distributions plot" = "p2"),
        selected = c("p1", "p2")
      ),
      
      labelled("Patients per arm (1 or 2 values):", n_help),
      textInput("n_text", label = NULL, value = "175, 113"),
      
      labelled("Test:", test_help),
      radioButtons(
        "test", label = NULL,
        choices = c("t-test" = "t", "Normal approximation" = "normal"),
        selected = "t", inline = TRUE
      ),
      
      fluidRow(
        column(4, numericInput("sd", "SD", value = 230, min = 0.0001)),
        column(4, tagList(
          labelled("alpha (1-sided)", alpha_help),
          numericInput("alpha", label = NULL, value = 0.025, min = 0.0001,
                       max = 0.5, step = 0.005)
        )),
        column(4, tagList(
          labelled("Target power (%)", target_help),
          numericInput("target_power", label = NULL, value = 90, min = 1,
                       max = 99.9, step = 5)
        ))
      ),
      
      tags$hr(),
      tags$b("Effect size range (x-axis data)"),
      fluidRow(
        column(4, numericInput("eff_min", "Min", value = 0)),
        column(4, numericInput("eff_max", "Max", value = 100)),
        column(4, numericInput("eff_step", "Step", value = 0.5, min = 0.0001))
      ),
      
      tags$hr(),
      tags$b("Reference lines (dashed)"),
      labelled("Horizontal, primary y-axis (%):", href_help),
      textInput("h_ref", label = NULL, value = "2.5, 50, 90"),
      checkboxInput("auto_v", "Vertical lines at power levels (%) for first N:",
                    value = TRUE),
      textInput("v_power", label = NULL, value = "50, 90"),
      labelled("Additional vertical lines (effect size):", vextra_help),
      textInput("v_extra", label = NULL, value = ""),
      tags$div(
        style = "display:flex; align-items:center;",
        checkboxInput("markers", "Show markers (first N)", value = TRUE),
        help_icon(marker_help)
      ),
      checkboxInput("gridlines", "Show gridlines", value = TRUE),
      
      tags$hr(),
      tags$b("Axis ranges (also updated by zooming on the plot)"),
      fluidRow(
        column(6, numericInput("x_min", "x min", value = 0)),
        column(6, numericInput("x_max", "x max", value = 100))
      ),
      fluidRow(
        column(6, numericInput("y_min", "y min (%)", value = 0)),
        column(6, numericInput("y_max", "y max (%)", value = 100))
      ),
      actionButton("reset_axes", "Reset axes to data range"),
      
    ),
    mainPanel(
      width = 8,
      tags$div(
        style = "margin-bottom:10px; display:flex; flex-wrap:wrap; gap:6px; align-items:center;",
        tags$b("Save:"),
        tags$button(
          type = "button", class = "btn btn-default",
          onclick = "downloadPlotPng('plot', 'type1_type2_error', 1400, 800);",
          icon("download"), " Power plot (PNG)"
        ),
        tags$button(
          type = "button", class = "btn btn-default",
          onclick = paste0("downloadPlotPng('dist_plot', ",
                           "'null_alternative_distributions', 1400, 700);"),
          icon("download"), " Distributions plot (PNG)"
        ),
        actionButton("dl_csv", "Values (CSV)", icon = icon("download")),
        actionButton("dl_pdf", "All results (PDF)", icon = icon("download"))
      ),
      conditionalPanel(
        "input.show_plots && input.show_plots.indexOf('p1') > -1",
        plotlyOutput("plot", height = "620px")
      ),
      conditionalPanel(
        "input.show_plots && input.show_plots.indexOf('p2') > -1",
        uiOutput("dist_status"),
        plotlyOutput("dist_plot", height = "460px")
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  # plotly checks, after the first flush, that every event read with
  # event_data() was registered by a rendered plot. The power plot cannot render
  # until the debounced inputs have settled, so that check ran too early and
  # produced "event ... is not registered" warnings. The plot does register
  # these events (event_register() in the plot code) once it renders, so list
  # them as registered up front.
  session$userData$plotlyShinyEventIDs <- unique(c(
    session$userData$plotlyShinyEventIDs,
    paste(c("plotly_relayout", "plotly_hover", "plotly_click"), "tm", sep = "-")
  ))
  
  # Debounce the free-text inputs so partially typed values are not evaluated
  ns_r     <- debounce(reactive(parse_nums(input$n_text)), 500)
  href_r   <- debounce(reactive(parse_nums(input$h_ref)), 500)
  vpower_r <- debounce(reactive(parse_nums(input$v_power)), 500)
  vextra_r <- debounce(reactive(parse_nums(input$v_extra)), 500)
  
  normal_r <- reactive(input$test == "normal")
  
  # ── Validated design parameters ─────────────────────────────────────────────
  design <- reactive({
    ns <- ns_r()
    validate(
      need(length(ns) %in% 1:2,
           "Enter one or two values for the number of patients per arm."),
      need(all(ns >= 2), "Patients per arm must be at least 2."),
      need(!is.na(input$sd) && input$sd > 0, "SD must be positive."),
      need(!is.na(input$alpha) && input$alpha > 0 && input$alpha < 0.5,
           "alpha (1-sided) must be between 0 and 0.5."),
      need(!is.na(input$target_power) && input$target_power > 0 &&
             input$target_power < 100,
           "Target power (%) must be a percentage between 0 and 100 (e.g. 90)."),
      need(!anyNA(c(input$eff_min, input$eff_max, input$eff_step)) &&
             input$eff_max > input$eff_min && input$eff_step > 0,
           "Effect size range: max must exceed min and step must be positive."),
      need((input$eff_max - input$eff_min) / input$eff_step <= 5000,
           "Too many grid points (limit 5000). Increase the step.")
    )
    list(ns = ns, sd = input$sd, alpha = input$alpha, tp = input$target_power / 100,
         normal = normal_r())
  })
  
  # ── Curves: one row per (N, effect size) ────────────────────────────────────
  curves <- reactive({
    d   <- design()
    eff <- seq(input$eff_min, input$eff_max, by = input$eff_step)
    map_dfr(d$ns, \(n) tibble(
      N      = n,
      Effect = eff,
      Power  = power_fun(eff, n, d$sd, d$alpha, d$normal),
      PValue = pval_fun(eff, n, d$sd, d$normal)
    ))
  })
  
  # ── Vertical reference lines ────────────────────────────────────────────────
  vlines <- reactive({
    d <- design()
    auto <- numeric(0)
    if (isTRUE(input$auto_v)) {
      p <- vpower_r() / 100
      p <- p[p > 0 & p < 1]
      if (length(p) > 0) {
        auto <- effect_at_power(p, d$ns[1], d$sd, d$alpha, d$normal)
      }
    }
    v <- c(auto, vextra_r())
    sort(unique(round(v[!is.na(v)], 6)))
  })
  
  # Axis ranges used when the plot is (re)built. Isolated so that zooming
  # (which updates the numeric inputs) does not rebuild the plot.
  axis_ranges <- function() {
    xr <- c(input$x_min, input$x_max)
    yr <- c(input$y_min, input$y_max) / 100
    if (anyNA(xr) || xr[1] >= xr[2]) xr <- c(input$eff_min, input$eff_max)
    if (anyNA(yr) || yr[1] >= yr[2]) yr <- c(0, 1)
    list(x = xr, y = yr)
  }
  
  # ── Plot ────────────────────────────────────────────────────────────────────
  output$plot <- renderPlotly({
    req(isTRUE(input$plotly_ready))
    d    <- design()
    cv   <- curves()
    vl   <- vlines()
    hl   <- href_r() / 100
    ax   <- isolate(axis_ranges())
    f    <- plot_f()   # font / margin scale from the plot-size slider
    blue <- "#0000FF"
    red  <- "#FF0000"
    ns   <- d$ns
    
    p <- plot_ly(source = "tm")
    
    for (i in seq_along(ns)) {
      n    <- ns[i]
      di   <- cv[cv$N == n, ]
      dash <- c("solid", "dot")[i]
      lab  <- paste0("N = ", n, " per arm")
      
      p <- p |>
        add_trace(
          type = "scatter", mode = "lines", x = di$Effect, y = di$Power,
          yaxis = "y", name = paste0("Power, ", lab),
          line = list(color = blue, dash = dash, width = 2),
          hovertemplate = paste0(
            lab, "<br>Effect size: %{x:.2f}<br>Power: %{y:.2%}<extra></extra>"
          ),
          hoverlabel = list(bgcolor = blue, font = list(color = "white"))
        ) |>
        add_trace(
          type = "scatter", mode = "lines", x = di$Effect, y = di$PValue,
          yaxis = "y2", name = paste0("1-sided p-value, ", lab),
          line = list(color = red, dash = dash, width = 2),
          hovertemplate = paste0(
            lab, "<br>Effect size: %{x:.2f}<br>1-sided p-value: %{y:.4f}<extra></extra>"
          ),
          hoverlabel = list(bgcolor = red, font = list(color = "white"))
        )
    }
    
    # Markers for the first N
    if (isTRUE(input$markers)) {
      n1     <- ns[1]
      mdd    <- mdd_fun(n1, d$sd, d$alpha, d$normal)
      eff_tp <- effect_at_power(d$tp, n1, d$sd, d$alpha, d$normal)
      pw_mdd <- power_fun(mdd, n1, d$sd, d$alpha, d$normal)
      
      bx <- c(0, mdd)
      by <- c(power_fun(0, n1, d$sd, d$alpha, d$normal), pw_mdd)
      bt <- c(
        sprintf("N = %s: effect size 0, power = %.2f%% (equals alpha)", n1, 100 * by[1]),
        sprintf("N = %s: MDD = %.2f, power = %.2f%%", n1, mdd, 100 * pw_mdd)
      )
      if (!is.na(eff_tp)) {
        bx <- c(bx, eff_tp)
        by <- c(by, d$tp)
        bt <- c(bt, sprintf("N = %s: effect size %.2f gives target power %.1f%%",
                            n1, eff_tp, 100 * d$tp))
      }
      
      p <- p |>
        add_trace(
          type = "scatter", mode = "markers", x = bx, y = by, yaxis = "y",
          marker = list(color = blue, size = 9), showlegend = FALSE,
          hoverinfo = "text", text = bt,
          hoverlabel = list(bgcolor = blue, font = list(color = "white"))
        ) |>
        add_trace(
          type = "scatter", mode = "markers", x = mdd, y = d$alpha, yaxis = "y2",
          marker = list(color = red, size = 9), showlegend = FALSE,
          hoverinfo = "text",
          text = sprintf("N = %s: MDD = %.2f, 1-sided p-value = %.4f", n1, mdd, d$alpha),
          hoverlabel = list(bgcolor = red, font = list(color = "white"))
        )
    }
    
    # Reference lines as dashed shapes with value labels
    shapes <- list()
    annots <- list()
    for (h in hl) {
      shapes[[length(shapes) + 1]] <- list(
        type = "line", xref = "paper", x0 = 0, x1 = 1, yref = "y", y0 = h, y1 = h,
        line = list(dash = "dash", color = "black", width = 1)
      )
      # Primary axis units (blue), left side
      annots[[length(annots) + 1]] <- list(
        xref = "paper", x = 0.005, xanchor = "left", yref = "y", y = h,
        yanchor = "bottom", text = paste0(signif(100 * h, 4), "%"),
        showarrow = FALSE, font = list(size = 10 * f, color = blue)
      )
      # Secondary axis units (red), right side
      annots[[length(annots) + 1]] <- list(
        xref = "paper", x = 0.995, xanchor = "right", yref = "y2", y = h,
        yanchor = "bottom", text = format(signif(h, 4), scientific = FALSE),
        showarrow = FALSE, font = list(size = 10 * f, color = red)
      )
    }
    for (v in vl) {
      shapes[[length(shapes) + 1]] <- list(
        type = "line", xref = "x", x0 = v, x1 = v, yref = "paper", y0 = 0, y1 = 1,
        line = list(dash = "dash", color = "black", width = 1)
      )
      annots[[length(annots) + 1]] <- list(
        xref = "x", x = v, xanchor = "left", yref = "paper", y = 0.99,
        yanchor = "top", text = format(round(v, 2), nsmall = 2),
        showarrow = FALSE, font = list(size = 10 * f)
      )
    }
    
    ttl <- paste0(
      "1-sided alpha = ", d$alpha, ", SD = ", d$sd,
      if (d$normal) " (normal approximation)" else " (t-test)",
      "<br>N = ", ns[1], " per arm (solid)",
      if (length(ns) == 2) paste0("; N = ", ns[2], " per arm (dotted)") else ""
    )
    
    p |>
      layout(
        font = list(size = 12 * f),
        title = list(text = ttl, font = list(size = 14 * f)),
        hovermode = "closest",
        xaxis = list(
          title = list(text = "Effect size", standoff = 10),
          range = ax$x, zeroline = FALSE, automargin = TRUE,
          showline = TRUE, linecolor = "black", linewidth = 1, ticks = "outside",
          showgrid = isTRUE(input$gridlines), gridcolor = "#D0D0D0",
          tick0 = 0, dtick = nice_dtick(diff(ax$x)),
          minor = list(dtick = nice_dtick(diff(ax$x)) / 2,
                       showgrid = isTRUE(input$gridlines), gridcolor = "#EDEDED")
        ),
        yaxis = list(
          title = list(text = "Rejection Probability (Power)",
                       font = list(color = blue), standoff = 10),
          tickfont = list(color = blue), tickformat = ".1~%", tick0 = 0, dtick = 0.1,
          range = ax$y, zeroline = FALSE, automargin = TRUE,
          showline = TRUE, linecolor = blue, linewidth = 1, ticks = "outside",
          showgrid = isTRUE(input$gridlines), gridcolor = "#D0D0D0",
          minor = list(dtick = 0.05, showgrid = isTRUE(input$gridlines),
                       gridcolor = "#EDEDED")
        ),
        yaxis2 = list(
          title = list(text = "1-sided p-value", font = list(color = red),
                       standoff = 10),
          tickfont = list(color = red), tickformat = ".3f", tick0 = 0, dtick = 0.1,
          overlaying = "y", side = "right", range = ax$y,
          showgrid = FALSE, zeroline = FALSE, automargin = TRUE,
          showline = TRUE, linecolor = red, linewidth = 1, ticks = "outside"
        ),
        shapes = shapes,
        annotations = annots,
        legend = list(orientation = "h", x = 0.5, xanchor = "center",
                      yref = "container", y = 0.01, yanchor = "bottom"),
        margin = list(t = 80 * f, b = 135 * f, r = 30)
      ) |>
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = list("toImage"),
        toImageButtonOptions = list(format = "png", filename = "type1_type2_error",
                                    width = 1400, height = 800, scale = 2)
      ) |>
      event_register("plotly_relayout") |>
      event_register("plotly_hover") |>
      event_register("plotly_click")
  })
  
  # ── Axis synchronisation: numeric inputs <-> plot zoom ──────────────────────
  
  # Reset numeric axis inputs when the data range changes or on request
  observeEvent(list(input$eff_min, input$eff_max, input$reset_axes), {
    req(!is.na(input$eff_min), !is.na(input$eff_max))
    updateNumericInput(session, "x_min", value = input$eff_min)
    updateNumericInput(session, "x_max", value = input$eff_max)
    updateNumericInput(session, "y_min", value = 0)
    updateNumericInput(session, "y_max", value = 100)
  }, ignoreInit = TRUE)
  
  # Numeric inputs -> plot
  observeEvent(list(input$x_min, input$x_max, input$y_min, input$y_max), {
    req(!anyNA(c(input$x_min, input$x_max, input$y_min, input$y_max)),
        input$x_min < input$x_max, input$y_min < input$y_max)
    yr <- c(input$y_min, input$y_max) / 100
    plotlyProxy("plot", session) |>
      plotlyProxyInvoke("relayout", list(
        `xaxis.range`  = c(input$x_min, input$x_max),
        `xaxis.dtick`  = nice_dtick(input$x_max - input$x_min),
        `xaxis.minor.dtick` = nice_dtick(input$x_max - input$x_min) / 2,
        `yaxis.range`  = yr,
        `yaxis2.range` = yr
      ))
  }, ignoreInit = TRUE)
  
  # Plot zoom -> numeric inputs
  observeEvent(ev_data("plotly_relayout"), {
    ev <- ev_data("plotly_relayout")
    if (!is.null(ev[["xaxis.autorange"]]) || !is.null(ev[["yaxis.autorange"]])) {
      updateNumericInput(session, "x_min", value = input$eff_min)
      updateNumericInput(session, "x_max", value = input$eff_max)
      updateNumericInput(session, "y_min", value = 0)
      updateNumericInput(session, "y_max", value = 100)
      return()
    }
    if (!is.null(ev[["xaxis.range[0]"]])) {
      updateNumericInput(session, "x_min", value = round(ev[["xaxis.range[0]"]], 3))
      updateNumericInput(session, "x_max", value = round(ev[["xaxis.range[1]"]], 3))
    }
    if (!is.null(ev[["yaxis.range[0]"]])) {
      updateNumericInput(session, "y_min", value = round(100 * ev[["yaxis.range[0]"]], 2))
      updateNumericInput(session, "y_max", value = round(100 * ev[["yaxis.range[1]"]], 2))
    }
  })
  
  # ── Display size: option text size (CSS) and plot heights ───────────────────
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
      ".well .form-group { margin-bottom: ", px(15 * s), "; }\n",
      ".well hr { margin-top: ", px(12 * s), "; margin-bottom: ", px(12 * s), "; }\n",
      ".well { padding: ", px(19 * s), "; }"
    )))
  })
  
  plot_scale_r <- debounce(
    reactive(if (is.null(input$plot_scale)) 100 else input$plot_scale), 400
  )
  # Font and margin scale for the plots: follows the plot-size slider, but never
  # below 0.7 so text stays legible. Equals 1 at 100%, so the plots look
  # exactly as before at full size.
  plot_f <- reactive(max(plot_scale_r() / 100, 0.7))
  
  # Plot heights are set with CSS on the existing plot containers (not by
  # recreating them), then a window resize event makes plotly redraw each plot
  # to fit its container.
  output$dyn_plot_css <- renderUI({
    s <- plot_scale_r() / 100
    tagList(
      tags$style(HTML(paste0(
        "#plot { height: ", round(620 * s), "px !important; }\n",
        "#dist_plot { height: ", round(460 * s), "px !important; }"
      ))),
      tags$script(HTML(paste0(
        "setTimeout(function(){ window.dispatchEvent(new Event('resize')); }, 100);",
        "setTimeout(function(){ window.dispatchEvent(new Event('resize')); }, 500);"
      )))
    )
  })
  
  # ── Plot visibility: at least one plot must remain selected ─────────────────
  rv <- reactiveValues(sel = NULL, locked = FALSE, last_show = c("p1", "p2"))
  
  observeEvent(input$show_plots, ignoreNULL = FALSE, {
    if (length(input$show_plots) == 0) {
      updateCheckboxGroupInput(session, "show_plots", selected = rv$last_show)
    } else {
      rv$last_show <- input$show_plots
    }
  })
  
  # ── Selection driving plot 2: hover on plot 1, click to lock / release ──────
  # sel = list(eff, n, src), src in "hover" / "locked". The N is that of the
  # curve hovered (power and p-value traces come in pairs per N; the marker
  # traces after them refer to the first N).
  sel_from_event <- function(ev, d, src) {
    if (is.null(ev) || is.null(ev$x) || length(ev$x) == 0) return(NULL)
    cn  <- ev$curveNumber[1]
    idx <- if (!is.null(cn) && cn < 2 * length(d$ns)) floor(cn / 2) + 1 else 1
    list(eff = as.numeric(ev$x[1]), n = d$ns[idx], src = src)
  }
  
  safe_design <- function() tryCatch(design(), error = function(e) NULL)
  
  observe({
    ev <- ev_data("plotly_hover", priority = "event")
    isolate({
      d <- safe_design()
      if (!isTRUE(rv$locked) && !is.null(d)) {
        s <- sel_from_event(ev, d, "hover")
        if (!is.null(s)) rv$sel <- s
      }
    })
  })
  
  observe({
    ev <- ev_data("plotly_click", priority = "event")
    isolate({
      if (isTRUE(rv$locked)) {
        # Second click releases the lock
        rv$locked <- FALSE
        if (!is.null(rv$sel)) rv$sel$src <- "hover"
      } else {
        d <- safe_design()
        s <- if (!is.null(d)) sel_from_event(ev, d, "locked") else NULL
        if (!is.null(s)) {
          rv$sel    <- s
          rv$locked <- TRUE
        }
      }
    })
  })
  
  # Release a lock whose N is no longer one of the entered values
  observe({
    d <- safe_design()
    if (!is.null(d) && isTRUE(rv$locked) && !is.null(rv$sel) &&
        !(rv$sel$n %in% d$ns)) {
      rv$locked <- FALSE
      rv$sel    <- NULL
    }
  })
  
  sel_thr <- throttle(reactive(rv$sel), 100)
  
  sel_r <- reactive({
    d <- design()
    s <- sel_thr()
    if (is.null(s) || !(s$n %in% d$ns)) {
      n1 <- d$ns[1]
      e  <- effect_at_power(d$tp, n1, d$sd, d$alpha, d$normal)
      if (is.na(e)) e <- mdd_fun(n1, d$sd, d$alpha, d$normal)
      s <- list(eff = e, n = n1, src = "default")
    }
    s
  })
  
  # Fixed axis ranges so the plot does not jump while hovering
  dist_rng <- reactive({
    d     <- design()
    se_n  <- std_err(d$ns, d$sd)
    peak  <- max(mapply(function(n) dens_null(0, std_err(n, d$sd), 2 * n - 2, d$normal),
                        d$ns))
    list(
      x    = c(min(0, input$eff_min) - 4 * max(se_n),
               max(0, input$eff_max) + 4 * max(se_n)),
      ymax = 1.18 * peak
    )
  })
  
  output$dist_status <- renderUI({
    d  <- design()
    s  <- sel_r()
    st <- dist_stats(s, d)
    f2 <- function(x) formatC(x, format = "f", digits = 2)
    badge_style <- function(bg) paste0(
      "background:", bg, "; color:white; padding:2px 8px; border-radius:3px; ",
      "font-weight:bold; margin-right:8px;")
    badge <- switch(
      s$src,
      locked  = tags$span("LOCKED", style = badge_style("#9B1C1C")),
      hover   = tags$span("LAST HOVERED", style = badge_style("#555555")),
      default = tags$span("DEFAULT", style = badge_style("#999999"))
    )
    hint <- if (!("p1" %in% input$show_plots)) {
      "Show the power and p-value plot to hover over it or click it."
    } else if (s$src == "locked") {
      "Locked. Click the power and p-value plot again to release."
    } else if (s$src == "hover") {
      "Click the power and p-value plot to lock this position."
    } else {
      "Hover over a curve in the power and p-value plot to change this; click to lock."
    }
    tags$div(
      style = "margin: 10px 0 0 0;",
      badge,
      sprintf(
        paste0("Effect size under H1 = %s, N = %s per arm, SD = %s, SE = %s, ",
               "critical value = %s, alpha = %s%%, beta = %s%%, power = %s%%."),
        f2(s$eff), s$n, format(d$sd), f2(st$se), f2(st$crit),
        format(100 * d$alpha, digits = 3),
        formatC(round(100 * st$beta, 1), format = "f", digits = 1),
        formatC(100 - round(100 * st$beta, 1), format = "f", digits = 1)
      ),
      tags$br(),
      tags$em(hint)
    )
  })
  
  output$dist_plot <- renderPlotly({
    req(isTRUE(input$plotly_ready))
    d   <- design()
    rng <- dist_rng()
    f   <- plot_f()   # font / margin scale from the plot-size slider
    cc  <- dist_calc(isolate(sel_r()), d, rng, f)
    
    plot_ly() |>
      add_trace(
        type = "scatter", mode = "lines", x = cc$xs, y = cc$yn, name = cc$names[1],
        line = list(color = grey_col, dash = "solid", width = 2), hoverinfo = "skip"
      ) |>
      add_trace(
        type = "scatter", mode = "lines", x = cc$xs, y = cc$ya, name = cc$names[2],
        line = list(color = grey_col, dash = "dot", width = 2), hoverinfo = "skip"
      ) |>
      add_trace(
        type = "scatter", mode = "lines", x = cc$xa, y = cc$yar, name = cc$names[3],
        fill = "tozeroy", fillcolor = alpha_col,
        line = list(color = "rgba(0,0,0,0)", width = 0), hoverinfo = "skip"
      ) |>
      add_trace(
        type = "scatter", mode = "lines", x = cc$xb, y = cc$ybr, name = cc$names[4],
        fill = "tozeroy", fillcolor = beta_col,
        line = list(color = "rgba(0,0,0,0)", width = 0), hoverinfo = "skip"
      ) |>
      layout(
        font = list(size = 12 * f),
        title = list(text = cc$title, font = list(size = 14 * f)),
        hovermode = FALSE,
        xaxis = list(
          title = list(text = "Effect size", standoff = 10),
          range = rng$x, fixedrange = TRUE, zeroline = FALSE,
          showline = TRUE, linecolor = "black", linewidth = 1, ticks = "outside",
          showgrid = isTRUE(input$gridlines), gridcolor = "#D0D0D0",
          tick0 = 0, dtick = nice_dtick(diff(rng$x)),
          minor = list(dtick = nice_dtick(diff(rng$x)) / 2,
                       showgrid = isTRUE(input$gridlines), gridcolor = "#EDEDED")
        ),
        yaxis = list(
          title = list(text = "Density", standoff = 10),
          range = c(0, rng$ymax), fixedrange = TRUE, zeroline = FALSE,
          showgrid = FALSE, showticklabels = FALSE,
          showline = TRUE, linecolor = "black", linewidth = 1
        ),
        shapes = cc$shapes,
        annotations = cc$annots,
        legend = list(orientation = "h", x = 0.5, xanchor = "center",
                      yref = "container", y = 0.01, yanchor = "bottom"),
        margin = list(t = 115 * f, b = 135 * f, r = 30)
      ) |>
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = list("toImage"),
        toImageButtonOptions = list(format = "png", filename = "null_alternative_distributions",
                                    width = 1400, height = 700, scale = 2)
      )
  })
  
  # Update plot 2 in place when the hovered / locked selection changes
  observeEvent(sel_thr(), ignoreInit = TRUE, {
    req(sel_thr())
    d <- safe_design()
    req(d)
    cc <- dist_calc(sel_r(), d, dist_rng(), plot_f())
    plotlyProxy("dist_plot", session) |>
      plotlyProxyInvoke(
        "restyle",
        list(
          x    = list(I(cc$xs), I(cc$xs), I(cc$xa), I(cc$xb)),
          y    = list(I(cc$yn), I(cc$ya), I(cc$yar), I(cc$ybr)),
          name = as.list(cc$names)
        ),
        list(0L, 1L, 2L, 3L)
      ) |>
      plotlyProxyInvoke(
        "relayout",
        list(`title.text` = cc$title, shapes = cc$shapes, annotations = cc$annots)
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
          function() paste0("type1_type2_error_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          "text/csv", function(f) {
    d <- design()
    out <- curves() |>
      transmute(
        N_per_arm        = N,
        Effect_size      = Effect,
        Power            = Power,
        P_value_1sided   = PValue,
        Test             = if (d$normal) "normal approximation" else "t-test",
        SD               = d$sd,
        Alpha_1sided     = d$alpha
      )
    write.csv(out, f, row.names = FALSE)
  })

  text_page <- function(lines) {
    lines <- unlist(lapply(lines, function(l) if (nzchar(l)) strwrap(l, width = 105) else ""))
    grid::grid.newpage()
    grid::grid.text(paste(lines, collapse = "\n"), x = grid::unit(0.5, "in"),
                    y = grid::unit(1, "npc") - grid::unit(0.5, "in"), just = c("left", "top"),
                    gp = grid::gpar(fontsize = 10))
    export_footer()
  }

  # Static versions of the plots for the PDF (plotly cannot be written to PDF from R)
  pdf_power_gg <- function() {
    d  <- design(); cv <- curves(); vl <- vlines(); hl <- href_r() / 100
    ns <- d$ns
    xr <- c(input$x_min, input$x_max); yr <- c(input$y_min, input$y_max) / 100
    if (anyNA(xr) || xr[1] >= xr[2]) xr <- c(input$eff_min, input$eff_max)
    if (anyNA(yr) || yr[1] >= yr[2]) yr <- c(0, 1)
    nl   <- paste0("N = ", ns, " per arm")
    long <- bind_rows(mutate(cv, Series = "Power", Value = Power),
                      mutate(cv, Series = "1-sided p-value", Value = PValue)) |>
      mutate(Nlab = factor(paste0("N = ", N, " per arm"), levels = nl))
    p <- ggplot(long, aes(Effect, Value, colour = Series, linetype = Nlab)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = c(Power = "#0000FF", `1-sided p-value` = "#FF0000"), name = NULL) +
      scale_linetype_manual(values = c("solid", "dotted")[seq_along(ns)], name = NULL) +
      scale_y_continuous(labels = scales::percent, name = "Rejection probability (power)",
                         sec.axis = sec_axis(~ ., name = "1-sided p-value")) +
      coord_cartesian(xlim = xr, ylim = yr) +
      labs(x = "Effect size",
           title = paste0("1-sided alpha = ", d$alpha, ", SD = ", d$sd,
                          if (d$normal) " (normal approximation)" else " (t-test)"),
           subtitle = paste0("N = ", ns[1], " per arm (solid)",
                             if (length(ns) == 2) paste0("; N = ", ns[2], " per arm (dotted)") else "")) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      guides(colour   = guide_legend(order = 1, override.aes = list(linetype = "solid")),
             linetype = guide_legend(order = 2, override.aes = list(colour = "black")))
    if (length(hl)) p <- p + geom_hline(yintercept = hl, linetype = "dashed", colour = "grey30")
    if (length(vl)) p <- p + geom_vline(xintercept = vl, linetype = "dashed", colour = "grey30")
    if (isTRUE(input$markers)) {
      n1  <- ns[1]
      mdd <- mdd_fun(n1, d$sd, d$alpha, d$normal)
      et  <- effect_at_power(d$tp, n1, d$sd, d$alpha, d$normal)
      mk  <- data.frame(Effect = c(0, mdd, et),
                        Value  = c(power_fun(0, n1, d$sd, d$alpha, d$normal),
                                   power_fun(mdd, n1, d$sd, d$alpha, d$normal), d$tp))
      mk <- mk[!is.na(mk$Effect), ]
      p <- p + geom_point(data = mk, aes(Effect, Value), inherit.aes = FALSE, colour = "#0000FF", size = 2.5) +
        geom_point(data = data.frame(Effect = mdd, Value = d$alpha), aes(Effect, Value),
                   inherit.aes = FALSE, colour = "#FF0000", size = 2.5)
    }
    p
  }

  pdf_dist_gg <- function() {
    d <- design(); s <- sel_r(); rng <- dist_rng()
    cc <- dist_calc(s, d, rng, 1); st <- dist_stats(s, d)
    ttl <- gsub("</?sup>", "", gsub("<br>", "\n", cc$title))
    ln  <- rbind(data.frame(x = cc$xs, y = cc$yn, g = cc$names[1]),
                 data.frame(x = cc$xs, y = cc$ya, g = cc$names[2]))
    ln$g <- factor(ln$g, levels = cc$names[1:2])
    ggplot() +
      geom_area(data = data.frame(x = cc$xa, y = cc$yar, f = cc$names[3]), aes(x, y, fill = f), alpha = 0.55) +
      geom_area(data = data.frame(x = cc$xb, y = cc$ybr, f = cc$names[4]), aes(x, y, fill = f), alpha = 0.55) +
      geom_line(data = ln, aes(x, y, linetype = g), colour = "#6E6E6E", linewidth = 0.9) +
      geom_vline(xintercept = st$crit, linetype = "dashed") +
      scale_fill_manual(values = setNames(c("#B14535", "#1B5E20"), cc$names[3:4]), name = NULL) +
      scale_linetype_manual(values = setNames(c("solid", "dotted"), cc$names[1:2]), name = NULL) +
      coord_cartesian(xlim = rng$x, ylim = c(0, rng$ymax)) +
      labs(x = "Effect size", y = "Density", title = ttl,
           subtitle = sprintf("Critical value = %.2f, power = %.1f%%", st$crit, 100 * st$power)) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom", legend.box = "vertical") +
      guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 2))
  }

  pdf_lines <- function() {
    d <- design(); s <- sel_r(); st <- dist_stats(s, d)
    f2 <- function(x) formatC(x, format = "f", digits = 2)
    per_n <- vapply(d$ns, function(n) {
      sprintf("N = %s per arm: minimum detectable difference %s; effect size for %s%% power %s.",
              n, f2(mdd_fun(n, d$sd, d$alpha, d$normal)), format(100 * d$tp),
              f2(effect_at_power(d$tp, n, d$sd, d$alpha, d$normal)))
    }, "")
    c("Type I and II error (power) explorer: 2-arm trial, 1-sided test.",
      sprintf("SD = %s, 1-sided alpha = %s, %s, target power %s%%.", format(d$sd), format(d$alpha),
              if (d$normal) "normal approximation" else "t-test", format(100 * d$tp)),
      "", per_n, "",
      "Distributions plot selection:",
      sprintf("Effect size under H1 = %s, N = %s per arm, SE = %s, critical value = %s, alpha = %s%%, beta = %s%%, power = %s%%.",
              f2(s$eff), s$n, f2(st$se), f2(st$crit), format(100 * d$alpha, digits = 3),
              formatC(round(100 * st$beta, 1), format = "f", digits = 1),
              formatC(100 - round(100 * st$beta, 1), format = "f", digits = 1)))
  }
  on_save("dl_pdf", function() "type1_type2_error.pdf", "application/pdf", function(f) {
    show <- input$show_plots
    g1 <- if ("p1" %in% show) pdf_power_gg()
    g2 <- if ("p2" %in% show) pdf_dist_gg()
    l  <- pdf_lines()
    pdf(f, width = 10, height = 7); on.exit(dev.off())
    if (!is.null(g1)) { print(g1, vp = footer_vp()); export_footer() }
    if (!is.null(g2)) { print(g2, vp = footer_vp()); export_footer() }
    text_page(l)
  })
}

# ── Run ───────────────────────────────────────────────────────────────────────

shinyApp(ui = ui, server = server)