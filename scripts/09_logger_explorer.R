# ══════════════════════════════════════════════════════════════════════════════
#  09_logger_explorer.R — interactive version of 08_logger_check.R
#
#  Click a logger on the reef map, then move (or play) the hour slider to watch
#  the shadow move across the zoomed-in view. For each hour it shows the shading
#  re-extracted now next to the value stage 02 stored, so a MISMATCH or a
#  logger sitting on a shadow edge is easy to spot.
#
#  Start it from a terminal in the project folder:
#     Rscript scripts/09_logger_explorer.R
#  or in RStudio:  shiny::runApp("scripts/09_logger_explorer.R")
#  It opens in your browser. Stop it with Ctrl+C (terminal) or the red Stop
#  button (RStudio). Needs the shiny package once:  install.packages("shiny")
#
#  Speed: each new hour is one ray-trace of the whole reef (a few seconds, the
#  same call as stage 02). Hours already seen are kept in memory, so once a
#  day is computed ("Ray-trace the whole day") the slider and Play are instant.
#
#  PROVENANCE: [NEW] Claude, 22 Oct 2026. Same shading code as 02
#    (make_shadow, refracted sun, BUFFER_M extraction).
# ══════════════════════════════════════════════════════════════════════════════

if (!requireNamespace("shiny", quietly = TRUE))
  stop("Install shiny once with: install.packages(\"shiny\")", call. = FALSE)
source(here::here("config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "functions_raytrace.R"))
suppressPackageStartupMessages({ library(dplyr); library(shiny) })

if (!exists(".SHADE_FIX", envir = globalenv())) init_geometry_checks()

EXPLORER_MODEL <- FIG_MAIN_MODEL
OVERVIEW_FACT  <- 4        # overview map is drawn at 1/4 resolution (speed)
TOL            <- 0.01     # |now − stage 02| above this is flagged MISMATCH

geom  <- readRDS(in_path("geom.rds"))
cache <- readRDS(in_path("cache.rds"))
preds <- if (file.exists(in_path("predictions.rds")))
  readRDS(in_path("predictions.rds")) %>% filter(model == EXPLORER_MODEL) else NULL
EVENTS <- sort(intersect(unique(cache$sitekey),
                         sub("\\.rds$", "", list.files(in_path("dems"), "\\.rds$"))))
LT <- function(t) format(t, "%H:%M", tz = LOCAL_TZ)

# ── per-event data, loaded once ───────────────────────────────────────────────
.events <- new.env()
get_event <- function(key) {
  if (!is.null(.events[[key]])) return(.events[[key]])
  dem  <- load_dem(in_path("dems", paste0(key, ".rds")))
  hobo <- load_hobo_locations(geom, key)
  e    <- terra::ext(dem$rast)
  mask <- matrix_to_rast(dem$reef_mask * 1, dem$rast)
  elev <- terra::mask(dem$rast, mask, maskvalues = 0)
  amb  <- dem$ambient; amb[!dem$reef_mask] <- NA
  ev <- list(
    key = key, dem = dem, hobo = hobo, xy = sf::st_coordinates(hobo),
    ext = as.vector(e),                                      # xmin xmax ymin ymax
    res = terra::res(dem$rast)[1],
    cen = get_site_centroid(dem$rast),
    elev_m = terra::as.matrix(elev, wide = TRUE),              # map layout, row 1 = north
    amb_r  = matrix_to_rast(amb, dem$rast),
    amb_now = extract_at_loggers(amb, dem$rast, hobo, BUFFER_M),
    shadows = new.env())
  .events[[key]] <- ev
  ev
}

# ── one hour of shading, computed once and kept ──────────────────────────────
#  Stored compactly: the full-resolution map as 1-byte values (0–254, 255 = NA),
#  a 1/4-resolution overview, and the exact values extracted at every logger.
get_shadow <- function(ev, t) {
  id <- format(t, "%Y%m%d%H", tz = LOCAL_TZ)
  if (!is.null(ev$shadows[[id]])) return(ev$shadows[[id]])
  sun <- sun_at_midhour(t, ev$cen$lat, ev$cen$lon)
  if (sun$altitude <= 0) {
    ev$shadows[[id]] <- list(sun = sun, night = TRUE); return(ev$shadows[[id]])
  }
  r <- matrix_to_rast(make_shadow(ev$dem, sun, refract = TRUE), ev$dem$rast)
  m <- terra::as.matrix(r, wide = TRUE)
  code <- ifelse(is.finite(m), round(pmin(1, pmax(0, m)) * 254), 255)
  ev$shadows[[id]] <- list(
    sun = sun, night = FALSE, angles = raytrace_angles(sun, refract = TRUE),
    full = matrix(as.raw(code), nrow(m)),
    overview = terra::as.matrix(terra::aggregate(r, OVERVIEW_FACT, mean, na.rm = TRUE),
                                wide = TRUE),
    vals = terra::extract(r, terra::vect(sf::st_buffer(ev$hobo, BUFFER_M)),
                          fun = mean, na.rm = TRUE)[, 2])
  ev$shadows[[id]]
}
decode <- function(raw_m) { v <- as.integer(raw_m); v[v == 255L] <- NA; matrix(v / 254, nrow(raw_m)) }

# rows/cols of the map-layout matrix covering a window around (x, y)
window_idx <- function(ev, x, y, half) {
  nr <- nrow(ev$elev_m); nc <- ncol(ev$elev_m)
  c1 <- max(1, floor((x - half - ev$ext["xmin"]) / ev$res) + 1)
  c2 <- min(nc, floor((x + half - ev$ext["xmin"]) / ev$res) + 1)
  r1 <- max(1, floor((ev$ext["ymax"] - (y + half)) / ev$res) + 1)
  r2 <- min(nr, floor((ev$ext["ymax"] - (y - half)) / ev$res) + 1)
  list(rows = r1:r2, cols = c1:c2,
       xlim = ev$ext["xmin"] + c(c1 - 1, c2) * ev$res,
       ylim = ev$ext["ymax"] - c(r2, r1 - 1) * ev$res)
}
grey_img <- function(m) {
  col <- gray(pmin(1, pmax(0, ifelse(is.finite(m), 0.05 + 0.93 * m, 1))))
  col[!is.finite(m)] <- "white"
  as.raster(matrix(col, nrow(m)))
}
cached_val <- function(key, logger, t, col) {
  v <- cache[[col]][cache$sitekey == key & cache$hoboID_env == logger & cache$datetime_pos == t]
  if (length(v)) v[1] else NA_real_
}
verdict <- function(now, then) {
  if (!is.finite(then)) "no stage-02 value" else if (abs(now - then) <= TOL) "MATCH" else "MISMATCH"
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style("body{font-size:13px} .well{padding:10px} h4{margin-top:4px}")),
  titlePanel("Logger shadow explorer"),
  sidebarLayout(
    sidebarPanel(width = 3,
      selectInput("event", "Event", EVENTS, selected = if ("HS19" %in% EVENTS) "HS19" else EVENTS[1]),
      selectInput("day", "Day", character()),
      sliderInput("hour", "Hour (local, sun at mid-hour)", min = 6, max = 18, value = 12,
                  step = 1, animate = animationOptions(interval = 1200, loop = TRUE)),
      sliderInput("zoom", "Zoom half-width (cm)", min = 5, max = 100, value = 25, step = 5),
      checkboxInput("contours", "Elevation contours on zoom", TRUE),
      actionButton("whole_day", "Ray-trace the whole day (for Play)"),
      hr(),
      selectInput("logger", "Logger (or click one on the map)", character()),
      helpText("Yellow = selected logger, red = others, blue square = zoom window.",
               "Values: 'now' is re-extracted from this map; 'stage 02' is what the",
               "models used (cache.rds). Hours with the sun below", SUN_ALT_MIN_DEG,
               "deg are drawn but were not modelled.")
    ),
    mainPanel(width = 9,
      fluidRow(
        column(6, h4(textOutput("map_title")),
               plotOutput("overview", height = "460px", click = "map_click")),
        column(6, h4(textOutput("zoom_title")), plotOutput("zoom", height = "460px"))),
      fluidRow(
        column(6, plotOutput("series", height = "300px")),
        column(6, plotOutput("bars", height = "300px"))),
      tableOutput("table")
    )
  )
)

# ── server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  ev <- reactive(get_event(input$event))

  observeEvent(input$event, {
    k <- input$event
    days <- sort(unique(format(cache$datetime_pos[cache$sitekey == k], "%Y-%m-%d", tz = LOCAL_TZ)))
    n    <- table(format(cache$datetime_pos[cache$sitekey == k], "%Y-%m-%d", tz = LOCAL_TZ))
    updateSelectInput(session, "day", choices = days, selected = names(which.max(n)))
    ids  <- sort(ev()$hobo$hoboID)
    with_data <- intersect(ids, unique(cache$hoboID_env[cache$sitekey == k]))
    lab  <- ifelse(ids %in% with_data, ids, paste(ids, "(no light data)"))
    updateSelectInput(session, "logger", choices = setNames(ids, lab),
                      selected = if (length(with_data)) with_data[1] else ids[1])
  })
  observeEvent(input$map_click, {
    xy <- ev()$xy; d <- sqrt((xy[, 1] - input$map_click$x)^2 + (xy[, 2] - input$map_click$y)^2)
    if (min(d) < 0.6) updateSelectInput(session, "logger", selected = ev()$hobo$hoboID[which.min(d)])
  })

  t_now <- reactive({
    req(input$day)
    as.POSIXct(sprintf("%s %02d:00:00", input$day, input$hour), tz = LOCAL_TZ)
  })
  sh <- reactive({
    e <- ev(); t <- t_now()
    id <- format(t, "%Y%m%d%H", tz = LOCAL_TZ)
    if (is.null(e$shadows[[id]]))
      withProgress(message = paste("Ray-tracing", e$key, LT(t)), value = 0.5, get_shadow(e, t))
    else get_shadow(e, t)
  })
  observeEvent(input$whole_day, {
    e <- ev(); hrs <- 6:18
    withProgress(message = paste("Ray-tracing", e$key, input$day), value = 0, {
      for (h in hrs) {
        get_shadow(e, as.POSIXct(sprintf("%s %02d:00:00", input$day, h), tz = LOCAL_TZ))
        incProgress(1 / length(hrs), detail = sprintf("%02d:00", h))
      }
    })
  })
  sel <- reactive({ req(input$logger); which(ev()$hobo$hoboID == input$logger) })

  output$map_title <- renderText({
    s <- sh(); t <- t_now()
    if (s$night) return(sprintf("%s %s: sun below horizon", input$event, LT(t)))
    sprintf("%s %s  sun %.0f deg (in water %.0f deg), bearing %.0f deg%s", input$event, LT(t),
            s$sun$altitude * 180 / pi, s$angles$alt_deg, s$angles$bearing,
            if (s$sun$altitude * 180 / pi <= SUN_ALT_MIN_DEG) "  [not modelled]" else "")
  })

  output$overview <- renderPlot({
    e <- ev(); s <- sh(); i <- sel()
    par(mar = c(0.5, 0.5, 0.5, 0.5))
    plot(NA, xlim = e$ext[1:2], ylim = e$ext[3:4], asp = 1, axes = FALSE, xlab = "", ylab = "")
    if (!s$night) rasterImage(grey_img(s$overview), e$ext["xmin"], e$ext["ymin"],
                              e$ext["xmax"], e$ext["ymax"], interpolate = FALSE)
    points(e$xy[, 1], e$xy[, 2], pch = 21, col = "red", cex = 1.4, lwd = 1.6)
    text(e$xy[, 1], e$xy[, 2], e$hobo$hoboID, pos = 3, cex = 0.6, col = "red3")
    half <- input$zoom / 100
    rect(e$xy[i, 1] - half, e$xy[i, 2] - half, e$xy[i, 1] + half, e$xy[i, 2] + half,
         border = "dodgerblue", lwd = 2)
    points(e$xy[i, 1], e$xy[i, 2], pch = 23, bg = "yellow", cex = 2)
  })

  output$zoom_title <- renderText({
    s <- sh(); i <- sel(); if (s$night) return("")
    then <- cached_val(input$event, input$logger, t_now(), "shadow_refr")
    sprintf("%s  direct light: now %.2f | stage 02 %s  [%s]   sky view %.2f",
            input$logger, s$vals[i], if (is.finite(then)) sprintf("%.2f", then) else "-",
            verdict(s$vals[i], then), ev()$amb_now[i])
  })

  output$zoom <- renderPlot({
    e <- ev(); s <- sh(); i <- sel(); x <- e$xy[i, 1]; y <- e$xy[i, 2]
    w <- window_idx(e, x, y, input$zoom / 100)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(1.8, 0.5, 0), cex.axis = 0.7)
    plot(NA, xlim = w$xlim, ylim = w$ylim, asp = 1, axes = FALSE,
         xlab = "cm east of logger", ylab = "cm north of logger")
    tx <- pretty((w$xlim - x) * 100); ty <- pretty((w$ylim - y) * 100)
    axis(1, at = x + tx / 100, labels = tx); axis(2, at = y + ty / 100, labels = ty); box()
    if (!s$night)
      rasterImage(grey_img(decode(s$full[w$rows, w$cols, drop = FALSE])),
                  w$xlim[1], w$ylim[1], w$xlim[2], w$ylim[2], interpolate = FALSE)
    if (input$contours) {
      z <- e$elev_m[w$rows, w$cols, drop = FALSE]
      if (any(is.finite(z)))
        contour(x = seq(w$xlim[1], w$xlim[2], length.out = ncol(z)),
                y = seq(w$ylim[1], w$ylim[2], length.out = nrow(z)),
                z = t(z[nrow(z):1, , drop = FALSE]), add = TRUE, col = "#E07B39",
                nlevels = 8, labcex = 0.6, drawlabels = TRUE)
    }
    a <- seq(0, 2 * pi, length.out = 100)
    polygon(x + BUFFER_M * cos(a), y + BUFFER_M * sin(a), border = "dodgerblue",
            col = adjustcolor("dodgerblue", 0.2), lwd = 2)
    others <- setdiff(seq_len(nrow(e$xy)), i)
    points(e$xy[others, 1], e$xy[others, 2], pch = 21, col = "red", cex = 1.4, lwd = 1.6)
    points(x, y, pch = 23, bg = "yellow", cex = 1.8)
    if (!s$night) {   # arrow pointing towards the sun
      b <- s$angles$bearing * pi / 180; L <- input$zoom / 100 * 0.8
      arrows(x, y, x + L * sin(b), y + L * cos(b), col = "gold3", lwd = 2, length = 0.1)
    }
  })

  output$series <- renderPlot({
    k <- input$event; lg <- input$logger; req(input$day)
    d <- cache %>% filter(sitekey == k, hoboID_env == lg,
                          format(datetime_pos, "%Y-%m-%d", tz = LOCAL_TZ) == input$day) %>%
      arrange(datetime_pos)
    p <- if (!is.null(preds)) preds %>% filter(sitekey == k, hoboID_env == lg) else NULL
    par(mar = c(3.5, 4, 2, 4), mgp = c(2, 0.6, 0))
    if (nrow(d) == 0) { plot.new(); title(paste(lg, "has no light data on", input$day)); return() }
    h  <- as.numeric(format(d$datetime_pos, "%H", tz = LOCAL_TZ)) + 0.5
    yl <- c(0, max(c(d$light_mol, d$surface_par * 0.6, p$pred), na.rm = TRUE))
    plot(h, d$light_mol, type = "b", pch = 16, col = "steelblue", ylim = yl, xlim = c(6, 19),
         xlab = "Hour", ylab = expression("PAR (" * mu * "mol m"^-2 * " s"^-1 * ")"),
         main = paste(lg, "on", input$day))
    if (!is.null(p) && nrow(p)) {
      j <- match(d$datetime_pos, p$datetime_pos)
      lines(h, p$pred[j], type = "b", pch = 17, col = "firebrick", lty = 2)
    }
    abline(v = input$hour + 0.5, col = "gold3", lwd = 2)
    par(new = TRUE)
    plot(h, d$shadow_refr, type = "b", pch = 1, col = "grey30", ylim = c(0, 1), xlim = c(6, 19),
         axes = FALSE, xlab = "", ylab = "")
    e <- ev(); i <- sel()
    now <- vapply(6:18, function(hh) {
      id <- sprintf("%s%02d", gsub("-", "", input$day), hh)
      s <- e$shadows[[id]]; if (is.null(s) || s$night) NA_real_ else s$vals[i]
    }, numeric(1))
    points(6:18 + 0.5, now, pch = 4, col = "black", cex = 1.2, lwd = 2)
    axis(4); mtext("direct-light fraction", side = 4, line = 2)
    legend("topleft", bty = "n", cex = 0.8,
           c("measured", paste("predicted", EXPLORER_MODEL), "direct light (stage 02)",
             "direct light (now)"),
           col = c("steelblue", "firebrick", "grey30", "black"), pch = c(16, 17, 1, 4),
           lty = c(1, 2, 1, NA))
  })

  output$bars <- renderPlot({
    req(input$day); t <- t_now()
    p <- if (!is.null(preds)) preds %>% filter(sitekey == input$event, datetime_pos == t) else NULL
    par(mar = c(4, 4, 2, 0.5), mgp = c(2, 0.6, 0))
    if (is.null(p) || nrow(p) == 0) {
      plot.new(); title(paste("No", EXPLORER_MODEL, "predictions at", LT(t))); return()
    }
    p  <- p[order(p$hoboID_env), ]
    yl <- c(0, max(c(p$light_mol, p$pred)) * 1.08)
    bp <- barplot(rbind(p$pred, p$light_mol), beside = TRUE, col = c("firebrick", "steelblue"),
                  border = NA, names.arg = p$hoboID_env, las = 2, cex.names = 0.6, ylim = yl,
                  ylab = "PAR", main = paste("All loggers at", LT(t)))
    j <- match(input$logger, p$hoboID_env)
    if (!is.na(j)) rect(bp[1, j] - 0.6, 0, bp[2, j] + 0.6, yl[2], border = "dodgerblue", lwd = 2)
    legend("topright", c(paste("predicted", EXPLORER_MODEL), "measured"), bty = "n",
           fill = c("firebrick", "steelblue"), border = NA, cex = 0.8)
  })

  output$table <- renderTable({
    e <- ev(); i <- sel(); req(input$day)
    rows <- lapply(6:18, function(hh) {
      id <- sprintf("%s%02d", gsub("-", "", input$day), hh)
      s <- e$shadows[[id]]; if (is.null(s) || s$night) return(NULL)
      t <- as.POSIXct(sprintf("%s %02d:00:00", input$day, hh), tz = LOCAL_TZ)
      then <- cached_val(input$event, input$logger, t, "shadow_refr")
      data.frame(hour = sprintf("%02d:00", hh), now = s$vals[i], stage02 = then,
                 check = verdict(s$vals[i], then))
    })
    do.call(rbind, rows)
  }, digits = 2, caption = "Hours ray-traced so far for this logger")
}

app <- shinyApp(ui, server)
if (sys.nframe() == 0L) runApp(app, launch.browser = TRUE) else app
