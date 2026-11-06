# Bayesian Distribution Explorer for Response Endpoint (Beta-Binomial)
# Shiny app — single file (app.R)

library(shiny)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggridges)
if (FALSE) library(jsonlite)  # used via jsonlite:: for downloads
library(scales)

# ── ESS helper ────────────────────────────────────────────────────────────────
# For a Beta(alpha, beta) distribution the ESS is alpha + beta.
# This is the standard conjugate-prior result: the prior carries the same
# information as (alpha + beta) Binomial observations.
ess_beta <- function(alpha, beta) alpha + beta

# ── Credible interval helpers ──────────────────────────────────────────────────

# Equal-tailed interval: cuts (1-conf)/2 probability from each tail.
eti_beta <- function(alpha, beta, conf) {
  tail_p <- (1 - conf) / 2
  c(lower = qbeta(tail_p, alpha, beta), upper = qbeta(1 - tail_p, alpha, beta))
}

# Highest posterior density interval: shortest interval containing `conf`
# probability mass. For a unimodal Beta (alpha > 1 and beta > 1) this is found
# by minimising the interval width over all starting quantiles p_low such that
# CDF(p_high) - CDF(p_low) = conf. Falls back to the equal-tailed interval for
# non-unimodal cases (alpha <= 1 or beta <= 1), where a single HPD interval
# need not exist.
hpd_beta <- function(alpha, beta, conf) {
  if (alpha <= 1 || beta <= 1) {
    return(eti_beta(alpha, beta, conf))
  }
  width <- function(p_low) {
    p_high <- p_low + conf
    if (p_high > 1) return(Inf)
    qbeta(p_high, alpha, beta) - qbeta(p_low, alpha, beta)
  }
  opt <- optimize(width, interval = c(0, 1 - conf))
  p_low <- opt$minimum
  c(lower = qbeta(p_low, alpha, beta), upper = qbeta(p_low + conf, alpha, beta))
}

# ── Successes input parser ─────────────────────────────────────────────────────
# Accepts either a comma-separated list (e.g. "4, 6, 8") or a sequence
# expression of the form "a to b by c" (e.g. "4 to 8 by 2").
parse_successes <- function(x) {
  x <- trimws(x)
  seq_match <- regmatches(
    x,
    regexec(
      "^([0-9]*\\.?[0-9]+)\\s+to\\s+([0-9]*\\.?[0-9]+)\\s+by\\s+([0-9]*\\.?[0-9]+)$",
      x, ignore.case = TRUE
    )
  )[[1]]

  if (length(seq_match) == 4) {
    from <- as.numeric(seq_match[2])
    to   <- as.numeric(seq_match[3])
    by   <- as.numeric(seq_match[4])
    successes <- seq(from, to, by = by)
  } else {
    successes <- suppressWarnings(as.numeric(unlist(strsplit(x, ","))))
  }

  successes <- successes[!is.na(successes)]
  if (length(successes) == 0) successes <- 0
  successes
}

# ── Help tooltip helper ────────────────────────────────────────────────────────
# Lightweight hover help using the native HTML title attribute — no extra
# package needed. Shows a small "?" icon; the tooltip text appears on hover.
help_icon <- function(text) {
  tags$span(
    "\u24d8",
    title = text,
    style = "cursor: help; color: #337ab7; margin-left: 6px; font-weight: bold;"
  )
}

alpha_beta_help_text <- paste0(
  "The response rate is modelled as a Beta(alpha, beta) distribution. ",
  "Alpha is the prior number of successes and beta the prior number of ",
  "failures (pseudo-observations). For example, alpha = 1 and beta = 1 gives ",
  "a Beta(1, 1) distribution, which is a Uniform distribution between 0 and 1 ",
  "\u2014 i.e. a flat, non-informative prior."
)

posterior_dist_help_text <- paste0(
  "The posterior distribution updates the prior Beta(alpha, beta) using the ",
  "observed data. With n trials and s observed successes, the posterior is ",
  "Beta(alpha + s, beta + n - s)."
)

posterior_mean_ci_help_text <- paste0(
  "The posterior mean is alpha / (alpha + beta). The credible interval is the ",
  "range containing the chosen percentage of posterior probability, computed ",
  "using the interval type selected below (HPD or ETI)."
)

hpd_help_text <- paste0(
  "Highest Posterior Density (HPD) interval: the narrowest interval that ",
  "contains the chosen probability mass, so every point inside the interval ",
  "has higher density than any point outside it."
)

eti_help_text <- paste0(
  "Equal-Tailed Interval (ETI): cuts an equal probability from each tail of ",
  "the posterior distribution, e.g. the 2.5th and 97.5th percentiles for a ",
  "95% interval."
)

response_range_help_text <- paste0(
  "Set the response rate range you want the posterior probability calculated ",
  "for. For example, if a response rate greater than or equal to 40% is of ",
  "interest (e.g. a \u201cGo Criteria\u201d), specify the region of interest as 40 to 100."
)

# ── UI ────────────────────────────────────────────────────────────────────────

# ---- footer for exported PNG and PDF files: notice, app address, UTC date and time ----
APP_URL <- "https://kinnersleyn.github.io/kcl-apps/Bayesian/"
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
    tags$link(rel = "icon", type = "image/jpeg", href = "icon.jpg"),
    tags$link(rel = "apple-touch-icon", href = "icon.jpg"),
    # Saves a file sent from R (base64) through the browser. Shiny's own download
    # links do not work in shinylive/webR.
    tags$script(HTML(paste0("$(document).on('shiny:connected', function(){",
      "Shiny.addCustomMessageHandler('save_file', function(m){",
      "var bin=atob(m.b64), n=bin.length, bytes=new Uint8Array(n);",
      "for (var i=0;i<n;i++) bytes[i]=bin.charCodeAt(i);",
      "var url=URL.createObjectURL(new Blob([bytes],{type:m.mime}));",
      "var a=document.createElement('a'); a.href=url; a.download=m.name;",
      "document.body.appendChild(a); a.click(); document.body.removeChild(a);",
      "setTimeout(function(){URL.revokeObjectURL(url);},10000);});});")))
  ),
  titlePanel("Bayesian Distribution Explorer for Response Endpoint (Beta-Binomial)"),
  sidebarLayout(
    sidebarPanel(
      tags$div(
        style = "display:flex; align-items:center;",
        tags$label("Alpha (prior):", style = "margin-bottom: 0;"),
        help_icon(alpha_beta_help_text)
      ),
      numericInput("alpha", label = NULL, value = 1, min = 0.1, step = 0.1),

      tags$div(
        style = "display:flex; align-items:center;",
        tags$label("Beta (prior):", style = "margin-bottom: 0;"),
        help_icon(alpha_beta_help_text)
      ),
      numericInput("beta_par", label = NULL, value = 1, min = 0.1, step = 0.1),

      checkboxInput("plot_prior", "Plot Prior (Beta Distribution)?", value = TRUE),
      numericInput("n", "Trial size (n):", value = 10, min = 1, step = 1),
      textInput(
        "successes",
        "Observed successes (comma-separated, or 'a to b by c'):",
        value = "5"
      ),
      tags$div(
        style = "display:flex; align-items:center;",
        tags$label("Response Region of Interest (%):", style = "margin-bottom: 0;"),
        help_icon(response_range_help_text)
      ),
      sliderInput(
        "response_range",
        label = NULL,
        min = 0, max = 100, value = c(40, 100), step = 1
      ),
      checkboxInput("shade", "Shade Response Region of Interest?", value = TRUE),
      radioButtons(
        "plot_type", "Choose Plot Type:",
        choices  = c("Panel Plot" = "panel", "Ridgeline Plot" = "ridgeline"),
        selected = "panel"
      ),
      numericInput(
        "ci_level",
        "Credible Interval Level (%):",
        value = 95, min = 50, max = 99.9, step = 0.1
      ),
      radioButtons(
        "ci_type",
        "Credible Interval Type:",
        choiceNames = list(
          tagList("Equal-Tailed Interval (ETI)", help_icon(eti_help_text)),
          tagList("Highest Posterior Density (HPD)", help_icon(hpd_help_text))
        ),
        choiceValues = c("eti", "hpd"),
        selected = "eti"
      ),
      actionButton("update", "Update Plot and Probability")
    ),
    mainPanel(
      tags$div(
        style = "margin-bottom:10px; display:flex; flex-wrap:wrap; gap:6px; align-items:center;",
        tags$b("Save (after Update):"),
        actionButton("dl_png", "Plots (PNG)", icon = icon("download")),
        actionButton("dl_csv", "Table (CSV)", icon = icon("download")),
        actionButton("dl_pdf", "All results (PDF)", icon = icon("download"))
      ),
      # Prior plot + ESS — hidden entirely when checkbox is off
      conditionalPanel(
        condition = "input.plot_prior == true",
        plotOutput("priorPlot", height = "200px"),
        verbatimTextOutput("priorESS")
      ),
      plotOutput("distPlot"),
      uiOutput("posteriorTable"),
      verbatimTextOutput("posteriorHeader")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Shared x grid — defined once, reused everywhere
  x_grid <- seq(0, 1, length.out = 1000)

  # Debounce the response-range slider: while dragging, values fire on every
  # tick, and re-running the shaded plot / posterior probability on each tick
  # feels sluggish (especially on mobile). Debouncing waits for 300ms of no
  # movement before propagating the new range downstream.
  response_range_db <- debounce(reactive(input$response_range), 300)

  # ── Reactive: validated prior parameters ──────────────────────────────────
  prior_params <- reactive({
    list(
      alpha = if (is.na(input$alpha)    || input$alpha    <= 0) 1 else input$alpha,
      beta  = if (is.na(input$beta_par) || input$beta_par <= 0) 1 else input$beta_par
    )
  })

  # ── Reactive: prior density data ──────────────────────────────────────────
  prior_data <- reactive({
    p <- prior_params()
    tibble(
      Probability = x_grid,
      Density     = dbeta(x_grid, p$alpha, p$beta)
    )
  })

  # ── Reactive: posterior parameters (updated on button click only) ─────────
  posterior_data <- eventReactive(input$update, {
    p <- prior_params()
    n <- input$n

    # Parse and validate successes (comma-separated or "a to b by c")
    successes <- parse_successes(input$successes)

    validate(
      need(
        all(successes >= 0 & successes <= n),
        paste0(
          "Observed successes must be between 0 and the trial size (n = ", n,
          "). Please check the values entered."
        )
      )
    )

    data.frame(
      Successes = successes,
      Alpha     = p$alpha + successes,
      Beta      = p$beta  + n - successes
    )
  })

  # ── Prior plot ────────────────────────────────────────────────────────────
  prior_plot <- reactive({
    p <- prior_params()
    ggplot(prior_data(), aes(x = Probability, y = Density)) +
      geom_line(linewidth = 1.2, colour = "black", linetype = "dashed") +
      scale_x_continuous(labels = label_percent()) +
      labs(
        title = paste0("Prior Distribution, Beta(", p$alpha, ", ", p$beta, ")"),
        x     = "Response Rate",
        y     = "Density"
      ) +
      theme_minimal()
  })

  output$priorPlot <- renderPlot(prior_plot())

  # ── Prior ESS ─────────────────────────────────────────────────────────────
  output$priorESS <- renderText({
    p   <- prior_params()
    ess <- ess_beta(p$alpha, p$beta)
    paste0("Prior Effective Sample Size (ESS) = ", round(ess, 2))
  })

  # ── Posterior plot ────────────────────────────────────────────────────────
  dist_plot <- reactive({
    pd <- posterior_data()

    # Build full density curves — one row per (Successes x x_grid point)
    plot_data <- pd |>
      rowwise() |>
      reframe(
        Probability = x_grid,
        Density     = dbeta(x_grid, Alpha, Beta),
        Successes   = Successes
      )

    # Shading: filter existing plot_data rather than recomputing dbeta()
    rng <- response_range_db()
    range_low_p  <- rng[1] / 100
    range_high_p <- rng[2] / 100

    shading_data <- if (input$shade) {
      plot_data |>
        filter(Probability >= range_low_p, Probability <= range_high_p)
    } else {
      NULL
    }

    if (input$plot_type == "panel") {

      base_plot <- ggplot(
        plot_data,
        aes(x = Probability, y = Density, colour = as.factor(Successes))
      ) +
        geom_line(linewidth = 1.2, show.legend = TRUE) +
        scale_x_continuous(labels = label_percent()) +
        facet_wrap(~ Successes, scales = "free_y", labeller = label_both) +
        labs(
          title  = paste0("Posterior Distributions by Observed Successes (n = ", input$n, ")"),
          x      = "Response Rate",
          y      = "Density",
          colour = "Successes"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom") +
        guides(col = guide_legend(ncol = 8, theme = theme(legend.byrow = TRUE)))

      if (!is.null(shading_data)) {
        base_plot <- base_plot +
          geom_ribbon(
            data        = shading_data,
            aes(x = Probability, ymin = 0, ymax = Density,
                fill = as.factor(Successes)),
            alpha       = 0.3,
            inherit.aes = FALSE
          ) +
          scale_fill_discrete(name = "Response Region of Interest") +
          guides(fill = "none")
      }

    } else {

      base_plot <- ggplot(
        plot_data,
        aes(x = Probability, y = Successes, height = Density,
            group = Successes, colour = as.factor(Successes))
      ) +
        geom_ridgeline(stat = "identity", alpha = 0) +
        scale_x_continuous(labels = label_percent()) +
        scale_y_continuous(labels = label_number(accuracy = 1)) +
        labs(
          title  = paste0("Posterior Distributions by Observed Successes (n = ", input$n, ")"),
          x      = "Response Rate",
          y      = "Successes",
          colour = "Successes"
        ) +
        theme_ridges(center_axis_labels = TRUE, grid = TRUE) +
        theme(legend.position = "bottom") +
        guides(col = guide_legend(ncol = 8, theme = theme(legend.byrow = TRUE)))

      if (!is.null(shading_data)) {
        base_plot <- base_plot +
          geom_ridgeline(
            data        = shading_data,
            aes(x = Probability, y = Successes, height = Density,
                group = Successes, colour = as.factor(Successes),
                fill  = as.factor(Successes)),
            stat        = "identity",
            alpha       = 0.3,
            inherit.aes = FALSE,
            show.legend = FALSE
          ) +
          theme_ridges(center_axis_labels = TRUE, grid = TRUE) +
          theme(legend.position = "bottom")
      }
    }

    base_plot
  })

  output$distPlot <- renderPlot(dist_plot())

  # ── Posterior mean & credible interval ─────────────────────────────────────
  # Depends only on posterior_data(), ci_level and ci_type — NOT on the
  # response-range slider — so dragging the slider never triggers the HPD
  # optimisation or ETI calculation again.
  posterior_ci <- reactive({
    pd <- posterior_data()
    conf <- input$ci_level / 100
    ci_fun <- if (input$ci_type == "hpd") hpd_beta else eti_beta

    ci <- Map(ci_fun, pd$Alpha, pd$Beta, MoreArgs = list(conf = conf))

    pd$Mean    <- pd$Alpha / (pd$Alpha + pd$Beta)
    pd$CI_low  <- sapply(ci, `[`, "lower")
    pd$CI_high <- sapply(ci, `[`, "upper")
    pd
  })

  # ── Posterior probability in region of interest ───────────────────────────
  # Depends on the debounced response range, so it only recalculates once the
  # user stops moving the slider, not on every intermediate tick.
  posterior_prob <- reactive({
    pd  <- posterior_data()
    rng <- response_range_db()
    range_low_p  <- rng[1] / 100
    range_high_p <- rng[2] / 100

    pd$Prob <- pbeta(range_high_p, pd$Alpha, pd$Beta) -
      pbeta(range_low_p, pd$Alpha, pd$Beta)
    pd
  })

  # ── Header text: ESS + region of interest ─────────────────────────────────
  output$posteriorHeader <- renderText({
    pd  <- posterior_data()
    rng <- response_range_db()
    ess <- ess_beta(pd$Alpha[1], pd$Beta[1])

    paste0(
      "Response Region of Interest: ",
      rng[1], "% - ", rng[2], "%\n",
      "Posterior Effective Sample Size (ESS) = ", round(ess, 2)
    )
  })

  # ── Output table: one column per successes scenario ───────────────────────
  output$posteriorTable <- renderUI({
    ci   <- posterior_ci()
    prob <- posterior_prob()
    conf <- input$ci_level / 100
    ci_label <- toupper(input$ci_type)

    header_cells <- lapply(ci$Successes, function(s) {
      tags$th(paste0("Successes = ", s))
    })

    dist_cells <- lapply(seq_len(nrow(ci)), function(i) {
      tags$td(paste0(
        "Beta(", round(ci$Alpha[i], 2), ", ", round(ci$Beta[i], 2), ")"
      ))
    })

    mean_ci_cells <- lapply(seq_len(nrow(ci)), function(i) {
      tags$td(paste0(
        scales::percent(ci$Mean[i], accuracy = 0.1), " (",
        scales::percent(ci$CI_low[i], accuracy = 0.1), " - ",
        scales::percent(ci$CI_high[i], accuracy = 0.1), ")"
      ))
    })

    prob_cells <- lapply(seq_len(nrow(prob)), function(i) {
      tags$td(sprintf("%.3f", prob$Prob[i]))
    })

    n_cols <- nrow(ci) + 1

    tags$table(
      class = "table table-bordered",
      tags$tr(
        tags$th(
          colspan = n_cols,
          style = "text-align: center;",
          paste0("Trial Size (n) = ", input$n)
        )
      ),
      tags$tr(tags$th("Metric"), header_cells),
      tags$tr(
        tags$td(tagList("Posterior Distribution", help_icon(posterior_dist_help_text))),
        dist_cells
      ),
      tags$tr(
        tags$td(tagList(
          paste0("Posterior Mean & ", scales::percent(conf, accuracy = 0.1), " CrInt (", ci_label, ")"),
          help_icon(posterior_mean_ci_help_text)
        )),
        mean_ci_cells
      ),
      tags$tr(
        tags$td("Posterior Probability in Region of Interest"),
        prob_cells
      )
    )
  })

  # ── Saving results ────────────────────────────────────────────────────────
  # Table as a data frame, one row per observed-successes scenario
  results_df <- function() {
    ci   <- posterior_ci()
    prob <- posterior_prob()
    rng  <- response_range_db()
    p    <- prior_params()
    data.frame(
      Trial_size_n       = input$n,
      Prior_alpha        = p$alpha,
      Prior_beta         = p$beta,
      Successes          = ci$Successes,
      Posterior_alpha    = ci$Alpha,
      Posterior_beta     = ci$Beta,
      Posterior_mean     = ci$Mean,
      CrI_level_pct      = input$ci_level,
      CrI_type           = toupper(input$ci_type),
      CrI_lower          = ci$CI_low,
      CrI_upper          = ci$CI_high,
      Region_lower_pct   = rng[1],
      Region_upper_pct   = rng[2],
      Prob_in_region     = prob$Prob
    )
  }

  # Prior plot (if shown) above the posterior plot
  draw_plots <- function() {
    gp <- list(if (isTRUE(input$plot_prior)) prior_plot(), dist_plot())
    gp <- Filter(Negate(is.null), gp)
    h  <- if (length(gp) == 2) c(1, 3) else 1
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(y = grid::unit(0.25, "in"), height = grid::unit(1, "npc") - grid::unit(0.25, "in"), just = "bottom",
                                      layout = grid::grid.layout(length(gp), 1, heights = grid::unit(h, "null"))))
    for (i in seq_along(gp)) {
      grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
      grid::grid.draw(ggplot2::ggplotGrob(gp[[i]]))
      grid::popViewport()
    }
    grid::popViewport()
    export_footer()
  }

  result_lines <- function() {
    df  <- results_df()
    tab <- data.frame(
      Successes = df$Successes,
      Posterior = sprintf("Beta(%s, %s)", round(df$Posterior_alpha, 2), round(df$Posterior_beta, 2)),
      `Mean (CrI)` = sprintf("%s (%s - %s)", scales::percent(df$Posterior_mean, accuracy = 0.1),
                             scales::percent(df$CrI_lower, accuracy = 0.1),
                             scales::percent(df$CrI_upper, accuracy = 0.1)),
      `P(region)` = sprintf("%.3f", df$Prob_in_region),
      check.names = FALSE)
    c(sprintf("Bayesian beta-binomial: prior Beta(%s, %s), trial size n = %s.",
              df$Prior_alpha[1], df$Prior_beta[1], df$Trial_size_n[1]),
      sprintf("Posterior mean with %s%% %s credible interval (CrI).", df$CrI_level_pct[1], df$CrI_type[1]),
      sprintf("Response region of interest: %s%% to %s%%.", df$Region_lower_pct[1], df$Region_upper_pct[1]),
      sprintf("Posterior effective sample size (ESS) = %s.", round(df$Posterior_alpha[1] + df$Posterior_beta[1], 2)),
      "", capture.output(print(tab, row.names = FALSE)))
  }

  send_file <- function(name, mime, writer) {
    f <- tempfile(); on.exit(unlink(f))
    writer(f)
    session$sendCustomMessage("save_file", list(name = name, mime = mime,
      b64 = jsonlite::base64_enc(readBin(f, "raw", file.size(f)))))
  }
  on_save <- function(btn, name, mime, writer) observeEvent(input[[btn]], {
    if (input$update < 1) {
      showNotification("Press Update Plot and Probability first, then save.", type = "warning")
      return()
    }
    tryCatch(send_file(name, mime, writer),
             error = function(e) showNotification(paste("Save failed:", conditionMessage(e)),
                                                  type = "error", duration = NULL))
  })
  on_save("dl_png", "posterior_plots.png", "image/png", function(f) {
    png(f, width = 10, height = if (isTRUE(input$plot_prior)) 10 else 7, units = "in", res = 150)
    on.exit(dev.off()); draw_plots()
  })
  on_save("dl_csv", "posterior_table.csv", "text/csv", function(f)
    write.csv(results_df(), f, row.names = FALSE))
  on_save("dl_pdf", "bayesian_results.pdf", "application/pdf", function(f) {
    pdf(f, width = 10, height = if (isTRUE(input$plot_prior)) 10 else 7)
    on.exit(dev.off()); draw_plots()
    grid::grid.newpage()
    grid::grid.text(paste(result_lines(), collapse = "\n"), x = grid::unit(0.5, "in"),
                    y = grid::unit(1, "npc") - grid::unit(0.5, "in"), just = c("left", "top"),
                    gp = grid::gpar(fontfamily = "mono", fontsize = 9))
    export_footer()
  })
}

# ── Run ───────────────────────────────────────────────────────────────────────

shinyApp(ui = ui, server = server)
