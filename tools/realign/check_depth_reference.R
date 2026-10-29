# ══════════════════════════════════════════════════════════════════════════════
#  check_depth_reference.R
#  Are the new (2016-aligned) DEMs on a tide datum, on each site's 2016
#  survey-day water surface, or not aligned in depth at all?
#
#  For every logger:
#     z_new   = new DEM height at the realigned position (median in a 6 cm circle)
#     dpHB    = old logger height (old DEM, which had that year's capture tide baked in)
#     offset  = z_new - dpHB - capture_tide_m(old, that site-year)
#
#  How to read the offsets:
#     ~0 everywhere                      -> new DEMs are on a tide datum: use capture tide 0
#     constant within a site, differs
#       between sites                    -> each site is on its own 2016 survey-day surface:
#                                           capture tide for that site = -offset
#     differs between years at a site    -> years NOT aligned in depth: depth did not carry over
#
#  Expect ~10-20 cm of noise: the old capture tides come from a harmonic fit
#  (RMSE 8-15 cm) and the old dpHB was a single-pixel value.
#
#  Run from the LightNiches folder. Reads only; writes three files to output/:
#     depth_reference_loggers.csv   one row per logger (includes z_new = new dpHB)
#     depth_reference_summary.csv   one row per site-year
#     depth_reference.png           offsets by site and year
# ══════════════════════════════════════════════════════════════════════════════

library(terra)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/LightNiches")

GEOM     <- "data/geom_new_realigned.csv"
DEM_DIR  <- "data/dems_new"
CAPTURE  <- "output/capture_tide.csv"     # old capture tides (from the last full run)
BUFFER_M <- 0.06                          # same as BUFFER_M in config.R
TMERC    <- paste0("+proj=tmerc +lat_0=-14.6989736557007 +lon_0=145.448257446289",
                   " +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# ── one row per logger ────────────────────────────────────────────────────────
loggers <- read.csv(GEOM) %>%
  distinct(rec, year, new_unit, .keep_all = TRUE) %>%
  transmute(site = rec, year = as.integer(year), sitekey = sprintf("%s%02d", rec, year),
            hoboID = new_unit, lon, lat, dpHB = suppressWarnings(as.numeric(dpHB)))

capture <- read.csv(CAPTURE) %>% select(sitekey, capture_tide_m)

# ── sample each new DEM at its loggers ────────────────────────────────────────
z_new <- lapply(split(loggers, loggers$sitekey), function(d) {
  f <- file.path(DEM_DIR, paste0(d$sitekey[1], ".tif"))
  if (!file.exists(f)) { message("no DEM for ", d$sitekey[1]); return(NULL) }
  dem  <- rast(f)
  pts  <- project(vect(d, geom = c("lon", "lat"), crs = TMERC), crs(dem))
  circ <- buffer(pts, BUFFER_M)
  d$z_new <- extract(dem, circ, fun = median, na.rm = TRUE)[, 2]
  message(sprintf("%s: %d loggers, %d with a DEM value", d$sitekey[1], nrow(d), sum(is.finite(d$z_new))))
  d
}) %>% bind_rows()

res <- z_new %>%
  left_join(capture, by = "sitekey") %>%
  mutate(diff   = z_new - dpHB,                 # how much higher the new DEM sits
         offset = diff - capture_tide_m)        # what's left after the old capture tide

# ── per site-year ─────────────────────────────────────────────────────────────
summ <- res %>%
  group_by(site, year, sitekey) %>%
  summarise(n = sum(is.finite(diff)),
            capture_tide_old = first(capture_tide_m),
            median_diff   = median(diff, na.rm = TRUE),
            spread_diff   = IQR(diff, na.rm = TRUE),        # small = consistent within the plot
            median_offset = median(offset, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 3)))

cat("\n── per site-year ─────────────────────────────────────────────\n")
cat("median_diff = new DEM height - old dpHB (m); spread = IQR between loggers\n")
cat("median_offset = median_diff - old capture tide (needs a capture time)\n\n")
print(as.data.frame(summ), row.names = FALSE)

# ── per site: do the years agree? ─────────────────────────────────────────────
by_site <- summ %>%
  filter(is.finite(median_offset)) %>%
  group_by(site) %>%
  summarise(years = paste(year, collapse = ","),
            offset_mean  = round(mean(median_offset), 3),
            offset_range = round(diff(range(median_offset)), 3),   # big = years disagree
            .groups = "drop")
cat("\n── per site (years with a capture time) ──────────────────────\n")
cat("offset_range small (< ~0.2 m) = the years agree, depth carried over\n\n")
print(as.data.frame(by_site), row.names = FALSE)

overall <- median(summ$median_offset, na.rm = TRUE)
cat(sprintf("\nOverall median offset: %.2f m  (near 0 suggests a tide datum)\n", overall))

# ── save + plot ───────────────────────────────────────────────────────────────
write.csv(res,  "output/depth_reference_loggers.csv", row.names = FALSE)
write.csv(summ, "output/depth_reference_summary.csv", row.names = FALSE)

p <- ggplot(filter(res, is.finite(offset)),
            aes(x = factor(year), y = offset, colour = factor(year))) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~ site, nrow = 1) +
  labs(x = NULL, y = "new DEM height - old dpHB - old capture tide (m)",
       title = "Depth reference of the new DEMs",
       subtitle = "On 0 = tide datum | flat but off 0 = 2016 survey surface | steps between years = not aligned in depth") +
  theme_minimal() + theme(legend.position = "none")
ggsave("output/depth_reference.png", p, width = 11, height = 4.5, dpi = 150)
print(p)
message("\nWrote output/depth_reference_loggers.csv, _summary.csv and depth_reference.png")
