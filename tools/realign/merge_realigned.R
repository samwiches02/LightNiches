# ══════════════════════════════════════════════════════════════════════════════
#  merge_realigned.R
#  Copies the realigned logger positions from QGIS into a NEW copy of geom_new.
#  geom_new.csv itself is only read, never written.
#
#  Input : geom_new.csv            original table (lon/lat = local tmerc metres)
#          hobos_realign.csv       exported from QGIS with Geometry = AS_XY,
#                                  CRS = Lizard local tmerc (columns X, Y, rec,
#                                  year, new_unit, status, ...)
#  Output: geom_new_realigned.csv  same rows and columns as geom_new, with lon/lat
#                                  replaced by the new positions, plus:
#                                    lon_orig, lat_orig  the original values
#                                    realign_status      manual / auto / unmoved
#                                    shift_m             distance moved
#          realign_summary.csv     one row per site/year
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)

GEOM_IN  <- "geom_new.csv"
NEW_IN   <- "hobos_realign.csv"
GEOM_OUT <- "geom_new_realigned.csv"
stopifnot(normalizePath(GEOM_OUT, mustWork = FALSE) != normalizePath(GEOM_IN))

geom <- read.csv(GEOM_IN, stringsAsFactors = FALSE, check.names = FALSE)
new  <- read.csv(NEW_IN,  stringsAsFactors = FALSE, check.names = FALSE)

# QGIS writes the coordinates as X / Y when exported with GEOMETRY = AS_XY
if (!all(c("X", "Y") %in% names(new))) {
  stop("hobos_realign.csv has no X/Y columns - export again with Layer options ",
       "GEOMETRY = AS_XY (see step 4).")
}
if (!"status" %in% names(new)) new$status <- NA_character_

new <- new %>%
  transmute(rec = as.character(rec), year = as.integer(year),
            new_unit = as.character(new_unit),
            lon_new = X, lat_new = Y, realign_status = status)

# [28 Oct 2026] guard: an export made while the QGIS layer was still filtered
#   only contains one site, and every other site silently keeps its ORIGINAL
#   position. Stop if the export covers less than 80% of the site-years.
n_geom <- nrow(dplyr::distinct(geom, rec, year))
n_new  <- nrow(dplyr::distinct(new, rec, year))
if (n_new < 0.8 * n_geom) {
  stop(sprintf(paste0("hobos_realign.csv only has %d of %d site-years (%s). ",
                      "Was a filter still on in QGIS? Clear the filter on hobos_realign, ",
                      "export again, then re-run this script."),
               n_new, n_geom,
               paste(unique(paste0(new$rec, sprintf("%02d", new$year %% 100))), collapse = ", ")),
       call. = FALSE)
}

# one position per logger, or the join would duplicate rows
dups <- new %>% count(rec, year, new_unit) %>% filter(n > 1)
if (nrow(dups) > 0) {
  print(dups)
  stop("Some loggers appear more than once in hobos_realign.csv (listed above).")
}

out <- geom %>%
  mutate(rec = as.character(rec), year = as.integer(year),
         new_unit = as.character(new_unit),
         .row = row_number()) %>%
  left_join(new, by = c("rec", "year", "new_unit")) %>%
  arrange(.row) %>%
  mutate(lon_orig = lon, lat_orig = lat,
         lon = coalesce(lon_new, lon),
         lat = coalesce(lat_new, lat),
         shift_m = sqrt((lon - lon_orig)^2 + (lat - lat_orig)^2)) %>%
  select(all_of(names(geom)), lon_orig, lat_orig, realign_status, shift_m)

stopifnot(nrow(out) == nrow(geom))          # same rows as the original

# ── checks ────────────────────────────────────────────────────────────────────
summary_tab <- out %>%
  distinct(rec, year, new_unit, .keep_all = TRUE) %>%
  group_by(rec, year) %>%
  summarise(loggers   = n(),
            manual    = sum(realign_status %in% "manual"),
            auto      = sum(realign_status %in% "auto"),
            unmoved   = sum(realign_status %in% "unmoved"),
            not_in_qgis = sum(is.na(realign_status)),
            median_shift_m = round(median(shift_m, na.rm = TRUE), 3),
            max_shift_m    = round(max(shift_m, na.rm = TRUE), 3),
            .groups = "drop")
print(as.data.frame(summary_tab))

if (any(summary_tab$unmoved > 0 | summary_tab$not_in_qgis > 0)) {
  message("\nNote: site/years with 'unmoved' or 'not_in_qgis' > 0 keep their ",
          "ORIGINAL coordinates for those loggers.")
}

write.csv(out, GEOM_OUT, row.names = FALSE)
write.csv(summary_tab, "realign_summary.csv", row.names = FALSE)
message("\nWrote ", GEOM_OUT, " (", nrow(out), " rows). ", GEOM_IN, " was not changed.")
