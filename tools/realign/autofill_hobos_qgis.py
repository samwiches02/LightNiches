# ══════════════════════════════════════════════════════════════════════════════
#  autofill_hobos_qgis.py  —  run inside QGIS (Plugins → Python Console → Show
#  Editor → open this file → Run ▶)
#
#  For every site/year in the layer "hobos_realign":
#    • loggers you moved by hand are ANCHORS (they stay exactly where you put them)
#    • with ≥ 3 anchors, an affine transform (shift + rotation + scale + skew)
#      is fitted from the anchors' old → new positions and applied to every
#      other logger of that site/year
#
#  Fields it uses / creates on the layer:
#    lon_old, lat_old  original position (created from the current position on
#                      the very first run if missing — so run it once BEFORE
#                      moving anything)
#    status            manual | auto | unmoved
#    loo_err_m         anchors only: miss distance when that anchor is left out
#                      of the fit (large → that anchor was probably misplaced)
#    auto_x, auto_y    where the script last put an auto point. If you nudge an
#                      auto point by hand it becomes an anchor on the next run.
#    snap_x, snap_y    where every logger was at the last run (the finder uses it
#                      to drop re-placed anchors from its "Bad anchors" list)
#
#  Safe to run as often as you like. Works on all sites at once, even if the
#  layer is filtered (the filter is restored afterwards).
# ══════════════════════════════════════════════════════════════════════════════

import math
import numpy as np
from qgis.core import QgsProject, QgsField, QgsGeometry, QgsPointXY
from qgis.PyQt.QtCore import QVariant

LAYER_NAME = "hobos_realign"
MIN_MOVE   = 0.001   # m — moved more than this = moved by hand
MAX_LOO    = 0.10    # m — flag anchors that disagree with the others by more than this
MIN_ANCHORS = 3


def fit_affine(src, dst):
    cen = src.mean(axis=0)
    A = np.column_stack([np.ones(len(src)), src - cen])
    B, *_ = np.linalg.lstsq(A, dst, rcond=None)
    return cen, B


def apply_affine(tr, xy):
    cen, B = tr
    return np.column_stack([np.ones(len(xy)), xy - cen]) @ B


def run(layer_name=LAYER_NAME):
    layers = QgsProject.instance().mapLayersByName(layer_name)
    if not layers:
        raise Exception(f"No layer called '{layer_name}' in the project.")
    lyr = layers[0]

    # save any pending hand edits first
    if lyr.isEditable():
        if not lyr.commitChanges():
            raise Exception("Could not save pending edits: " + "; ".join(lyr.commitErrors()))

    subset = lyr.subsetString()
    lyr.setSubsetString("")

    try:
        # ── make sure the fields exist ────────────────────────────────────────
        names = lyr.fields().names()
        first_run = "lon_old" not in names or "lat_old" not in names
        new = []
        for n, t in [("lon_old", QVariant.Double), ("lat_old", QVariant.Double),
                     ("status", QVariant.String), ("loo_err_m", QVariant.Double),
                     ("auto_x", QVariant.Double), ("auto_y", QVariant.Double),
                     ("snap_x", QVariant.Double), ("snap_y", QVariant.Double)]:
            if n not in names:
                new.append(QgsField(n, t, len=20 if t == QVariant.String else 0,
                                    prec=0 if t == QVariant.String else 6))
        if new:
            lyr.dataProvider().addAttributes(new)
            lyr.updateFields()

        idx = {n: lyr.fields().indexOf(n) for n in
               ["rec", "year", "new_unit", "lon_old", "lat_old",
                "status", "loo_err_m", "auto_x", "auto_y", "snap_x", "snap_y"]}
        if idx["rec"] < 0 or idx["year"] < 0:
            raise Exception("Layer needs 'rec' and 'year' fields.")

        feats = list(lyr.getFeatures())
        empty = [f for f in feats if not f.hasGeometry() or f.geometry().isNull() or f.geometry().isEmpty()]
        if empty:
            print(f"Skipping {len(empty)} feature(s) with no location, fid(s): "
                  + ", ".join(str(f.id()) for f in empty[:20]))
            feats = [f for f in feats if f not in empty]
        lyr.startEditing()

        if first_run:
            for f in feats:
                p = f.geometry().asPoint()
                lyr.changeAttributeValue(f.id(), idx["lon_old"], p.x())
                lyr.changeAttributeValue(f.id(), idx["lat_old"], p.y())
                lyr.changeAttributeValue(f.id(), idx["status"], "unmoved")
            lyr.commitChanges()
            print("First run: stored the current positions as lon_old/lat_old for",
                  len(feats), "loggers. Now move ~6 anchors per site and run again.")
            return

        def num(v):
            try:
                return None if v is None or (isinstance(v, float) and math.isnan(v)) else float(v)
            except (TypeError, ValueError):
                return None

        groups = {}
        for f in feats:
            groups.setdefault((str(f[idx["rec"]]), str(f[idx["year"]])), []).append(f)

        report = []
        for (rec, year), fs in sorted(groups.items()):
            now = np.array([[f.geometry().asPoint().x(), f.geometry().asPoint().y()] for f in fs])
            old = np.array([[num(f[idx["lon_old"]]), num(f[idx["lat_old"]])] for f in fs], dtype=float)
            if np.isnan(old).any():
                print(f"{rec}{year}: some lon_old/lat_old empty — skipped")
                continue
            ax = np.array([num(f[idx["auto_x"]]) if num(f[idx["auto_x"]]) is not None else np.nan for f in fs])
            ay = np.array([num(f[idx["auto_y"]]) if num(f[idx["auto_y"]]) is not None else np.nan for f in fs])

            moved = np.hypot(*(now - old).T) > MIN_MOVE
            untouched_auto = ~np.isnan(ax) & (np.hypot(now[:, 0] - ax, now[:, 1] - ay) <= MIN_MOVE)
            anchor = moved & ~untouched_auto
            a = np.where(anchor)[0]
            rest = np.where(~anchor)[0]

            status = ["manual" if anchor[i] else "unmoved" for i in range(len(fs))]
            loo = [None] * len(fs)
            new_xy = now.copy()
            new_xy[rest] = old[rest]          # anything not anchored starts from old
            auto_xy = [None] * len(fs)
            note = ""

            if len(a) >= MIN_ANCHORS:
                sv = np.linalg.svd(old[a] - old[a].mean(axis=0), compute_uv=False)
                if sv[1] / sv[0] < 0.1:
                    note = "anchors almost in a line — spread them out; not filled"
                else:
                    tr = fit_affine(old[a], now[a])
                    if len(rest):
                        new_xy[rest] = apply_affine(tr, old[rest])
                        for i in rest:
                            status[i] = "auto"
                            auto_xy[i] = new_xy[i]
                    res = np.hypot(*(apply_affine(tr, old[a]) - now[a]).T)
                    if len(a) >= 4:
                        for k, i in enumerate(a):
                            keep = np.delete(a, k)
                            t2 = fit_affine(old[keep], now[keep])
                            loo[i] = float(np.hypot(*(apply_affine(t2, old[[i]]) - now[[i]]).ravel()))
                    worst = max(a, key=lambda i: loo[i] if loo[i] is not None else -1)
                    wl = loo[worst]
                    note = (f"rmse {np.sqrt(np.mean(res**2)):.3f} m"
                            + (f", worst anchor {fs[worst][idx['new_unit']]} {wl:.3f} m" if wl is not None else "")
                            + ("  ⚠ CHECK THAT ANCHOR" if wl is not None and wl > MAX_LOO else ""))
            elif len(a):
                note = f"only {len(a)} anchor(s), need {MIN_ANCHORS}"

            for i, f in enumerate(fs):
                # absolute tolerance (np.allclose is relative: ~2.6 cm at these coordinates)
                if np.hypot(*(new_xy[i] - now[i])) > 1e-6:
                    lyr.changeGeometry(f.id(), QgsGeometry.fromPointXY(QgsPointXY(*new_xy[i])))
                lyr.changeAttributeValue(f.id(), idx["status"], status[i])
                lyr.changeAttributeValue(f.id(), idx["loo_err_m"], None if loo[i] is None else round(loo[i], 4))
                lyr.changeAttributeValue(f.id(), idx["auto_x"], None if auto_xy[i] is None else float(auto_xy[i][0]))
                lyr.changeAttributeValue(f.id(), idx["auto_y"], None if auto_xy[i] is None else float(auto_xy[i][1]))
                # snapshot of where every logger was when the errors were computed
                lyr.changeAttributeValue(f.id(), idx["snap_x"], float(new_xy[i][0]))
                lyr.changeAttributeValue(f.id(), idx["snap_y"], float(new_xy[i][1]))

            if len(a):
                report.append(f"{rec}{year}: {len(a)} anchors, {sum(s == 'auto' for s in status)} auto — {note}")

        if not lyr.commitChanges():
            raise Exception("Could not save: " + "; ".join(lyr.commitErrors()))

        print("\n".join(report) if report else "No anchors found yet — move some loggers first.")
    finally:
        lyr.setSubsetString(subset)
        lyr.triggerRepaint()


run()
