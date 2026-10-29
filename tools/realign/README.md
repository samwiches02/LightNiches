# Realigning HOBO logger positions onto new DEMs (QGIS)

Moves every logger from its position on the **old** DEM to the same spot on the
**new** DEM. You place about 10 loggers per site by hand ("anchors"), and a script
places the rest.

| File | Runs in | What it does |
|---|---|---|
| `hobo_finder.py` | QGIS Python console | Window that shows the old and new maps side by side. You click once to place a logger. |
| `autofill_hobos_qgis.py` | QGIS Python console | Fits a transform (shift, rotation, scale, skew) per site from your anchors and moves every other logger. |
| `merge_realigned.R` | R | Writes the new positions into a **copy** of `geom_new.csv`. |
| `check_depth_reference.R` | R | Samples the new DEM height at each logger (gives the new `dpHB`). |

All four live in `tools/realign/` in the LightNiches repo.

**To run a `.py` script in QGIS:** Plugins → Python Console → **Show Editor** (notepad
icon) → **Open Script** (folder icon) → choose the file → **Run ▶**.

---

## 1. One-time setup

1. **Custom CRS.** Settings → Custom Projections → **+**. Name: `Lizard local tmerc`.
   Format: Proj String. Paste the line below, then **Validate** → OK.
   ```
   +proj=tmerc +lat_0=-14.6989736557007 +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
   ```
2. **Load the loggers.** Layer → Add Layer → Add Delimited Text Layer → `data/geom_new.csv`.
   Set X = `lon`, Y = `lat`, CRS = **Lizard local tmerc**.
3. **One row per logger.**
   - Right-click the layer → Filter → `"L" = 0.25`.
   - Processing Toolbox → **Delete duplicates by attribute** (under *Vector general*,
     not WhiteboxTools), fields `rec`, `year`, `new_unit` → Run.
4. **Save the working layer.** Right-click the result → Export → Save Features As →
   GeoPackage, file `hobos_realign.gpkg`, layer `hobos_realign`, CRS **Lizard local tmerc**.
   Remove every other copy from the Layers panel; **keep exactly one `hobos_realign`**.
5. **Add the working fields.** Run `autofill_hobos_qgis.py` once. It stores each
   logger's original position (`lon_old`, `lat_old`) and adds `status`, `loo_err_m`,
   `auto_x`, `auto_y`.
6. **Style** (right-click `hobos_realign` → Properties):
   - **Symbology:**
     - Click **Simple Marker** → Fill color → data-defined button → Edit, and paste:
       ```
       CASE WHEN "status"='manual' THEN '#2ca02c' WHEN "status"='auto' THEN '#ff7f0e' ELSE '#d62728' END
       ```
       Green = anchor, orange = placed by the script, red = not moved.
     - Click **Marker** → **+** → Geometry Generator (LineString), and paste:
       ```
       make_line(make_point("lon_old", "lat_old"), $geometry)
       ```
       Then change the Simple Line under it to **Arrow** (yellow, thin).
   - **Labels:** Single Labels, value `new_unit`, tick Draw text buffer.
   - Style (bottom-left) → **Save as Default**, so the style survives reloading.
7. **Tidy up.** Drag `hobos_realign` to the top of the Layers panel. Project →
   Snapping Options → snapping **off**. Save the project.

**Layer groups.** Keep one group per site/year, each containing the old DEM and the
new DEM. Name the old DEM like the site/year (`CB17`, `TB19`, ...) so the finder can
pick it automatically.

---

## 2. For each site/year

1. **Show only this site.** Tick only that site/year's group.
2. **Filter the loggers.** Right-click `hobos_realign` → Filter, e.g.
   ```
   "rec" = 'TM' AND "year" = 17
   ```
   Site codes: `CB GT HS L1 N3 RS SE TB TM`. Years: `17 18 19`.
3. **Open the finder.** Run `hobo_finder.py`. Check that **Old map** and **New map**
   are the right way round.
4. **Place the anchors.** Set **Step through** to **Anchors to add**, up to **10**
   per site, then click **Next ▶**. For each logger:
   - Find the feature under the **red cross** (old map, close-up and overview).
   - Move your mouse over the main map until the **cursor view** shows the same
     feature under the **magenta cross** (or tick **Overlay**: Red/cyan turns grey
     when they line up).
   - **Left-click once.** The logger moves exactly there and the finder moves on to
     the next one. Cmd+Z undoes.
   - No clear feature? Click Next to skip it, and place a nearby logger instead.
5. **Save Layer Edits** (disk icon).
6. **Fill in the rest.** Run `autofill_hobos_qgis.py`. The console prints one line per
   site, e.g.
   ```
   TM17: 10 anchors, 40 auto — rmse 0.040 m, worst anchor H40 0.090 m
   ```
7. **Fix bad anchors.** Set **Step through** to **Bad anchors**, error > **15 cm**.
   Next ▶ goes worst first; re-place each one, then save and run the autofill again.
   Repeat until there's nothing above 15 cm (or what's left looks right on the map).
8. **Check the auto-placed loggers.** Set **Step through** to **All loggers**. The
   orange points should sit on the right features, and the arrows should change
   smoothly across the plot. Click to fix any that are off, then save and run the
   autofill again.

**Just one logger?** Type its ID in **Go to logger** (e.g. `H40`) → Enter → click to
re-place it.

---

## 3. When all sites are done

1. Run `autofill_hobos_qgis.py` one last time. Clear the filter, then Select by
   Expression `"status" = 'unmoved'`. Anything selected still has its old position.
2. **Export.** Right-click `hobos_realign` → Export → Save Features As:
   - Format **CSV**, file `data/hobos_realign.csv`
   - CRS **Lizard local tmerc**
   - Layer Options → **GEOMETRY = AS_XY**
3. **Merge.** In R:
   ```r
   setwd("~/Desktop/LightNiches/data"); source("../tools/realign/merge_realigned.R")
   ```
   This writes `data/geom_new_realigned.csv`. `geom_new.csv` itself is never changed.
4. **New logger heights.** In R:
   ```r
   source("~/Desktop/LightNiches/tools/realign/check_depth_reference.R")
   ```
   Then rebuild `data/dphb_newdem.csv` from `output/depth_reference_loggers.csv`
   (columns `site, year, unit, dpHB` = `site, year, hoboID, z_new`).
5. **Re-run the pipeline** from stage 01. Any site whose positions changed must be
   ray-traced again: delete its file in `output_newdem/raytrace/` first.

---

## Finder window reference

| Part | What it does |
|---|---|
| Old map, overview | Old DEM zoomed out. Red = this logger, yellow = the others, green circles = anchors to add, orange circles = bad anchors, red box = close-up area |
| Old map, close-up | Old DEM on the logger's original spot |
| New map, your cursor | New DEM under your mouse on the main map (magenta cross = mouse, cyan cross = original spot) |
| Overlay (tick to show) | Old map around the logger with the new map around your cursor on top. Red/cyan: grey = lined up |
| `1:___` boxes | Zoom of each view. Scrolling inside a view also changes it |
| Pan main map to logger | Next/Prev moves the main map, keeping its zoom |
| Click main map to place logger / then Next | One click places the logger; then moves on to the next |
| Step through | Anchors to add / Bad anchors / All loggers |
| Go to logger | Jump to any logger by ID |

## Reading the autofill output

- `rmse`: how well the transform fits your anchors. A few cm is normal hand error.
- `worst anchor Hxx 0.31 m  ⚠ CHECK THAT ANCHOR`: that anchor disagrees with the
  others by more than 10 cm. It's usually clicked on the wrong feature. Re-place it.
- Several anchors at one site all at 20–30 cm: matching is hard there. Add more
  anchors, or place every logger by hand.
- Each point's `status`: `manual` = your anchor, `auto` = placed by the script,
  `unmoved` = the site has fewer than 3 anchors.

## Troubleshooting

| Message | Fix |
|---|---|
| `No loggers in the current filter` | The filter matches nothing, or there are two `hobos_realign` layers. Keep one, check the filter with **Test**, and re-run the finder. |
| `All sites in the filter have their anchors` | That site already has enough. Run the autofill, or raise the "up to" number. |
| `Null geometry` / `skipping N feature(s) with no location` | A point with no position (usually an accidental Add Feature). It's skipped; delete it via its `fid` if it's empty. |
| `wrapped C/C++ object ... has been deleted` | The layer was removed while the finder was open. Close the finder and run it again. |
| `NetCDF: HDF error` (satellite Kd) | A `.nc` file is only partly downloaded. Delete it and download it again. |
