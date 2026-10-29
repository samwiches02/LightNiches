# ══════════════════════════════════════════════════════════════════════════════
#  hobo_finder.py  —  run inside QGIS (Python Console → Show Editor → open → Run ▶)
#
#  One tall floating window with three views stacked top to bottom, each with
#  its own zoom box:
#
#    1 "Old map — overview":  the old map zoomed out around the selected logger
#      (red cross), other loggers of the site as yellow crosses, and a red box
#      showing the area of the close-up view.
#    2 "Old map — close-up":  the old map zoomed in on the logger's ORIGINAL
#      position (red cross).
#    3 "New map — your cursor": a magnifier on the new map that follows your
#      mouse on the main QGIS map (magenta crosshair = your mouse; cyan cross =
#      the logger's original spot).
#
#    4 "Overlay" (optional, tick to show): the old map around the logger with
#      the new map around your cursor laid on top, both at the overlay's zoom.
#      When the two images line up, your cursor is on the logger's new spot.
#      Modes: Blend (slider = new-map opacity) or Red/cyan (old = red, new =
#      cyan; aligned → grey, misaligned → coloured fringes).
#
#  Leave the main map zoomed out, move the mouse until view 3 shows the same
#  feature as view 2 (or the overlay lines up), then click to place.
#
#  • Select a logger (click it / attribute table) or use ◀ Prev / Next ▶.
#  • "Step through:" picks what Prev/Next walks through:
#      Anchors to add  – loggers still needed to reach N anchors per site (below)
#      Bad anchors     – anchors whose error (loo_err_m, from the last autofill
#                        run) is above the cm threshold, worst first; click to
#                        re-place, then re-run the autofill
#      All loggers     – everything in the filter, in ID order
#  • "Go to logger": type an ID (e.g. H40) + Enter to jump to it and re-place it.
#  • Anchors to add, "up to N per site" (N = 10): Prev/Next
#    steps only through loggers you still need to place. Loggers you've already
#    placed by hand are skipped; each new suggestion is the logger farthest from
#    every placed anchor (fills the biggest gap). A site with none yet starts
#    with 5 around the edge + 1 near the centre, then fills gaps. Suggestions
#    are green circles on the overview. Untick it to step through every logger
#    (e.g. to check the auto points).
#  • Zoom each view with its own 1:___ box, or scroll inside the view (the box
#    follows). Drag inside a view to pan it.
#  • Old map / New map dropdowns: guessed from your layer groups (old = raster
#    named like the site/year, e.g. "CB17"; new = the other raster in the same
#    group). Change them if the guess is wrong.
#  • "Pan main map": moves the main map to each logger without changing its zoom.
#  • "Click to place" (on by default): one left-click on the main map puts the
#    selected logger EXACTLY under your cursor (editing is switched on for you;
#    Cmd/Ctrl+Z undoes). With "then Next" ticked it jumps to the next logger.
#    Untick it to get your normal QGIS tools back.
# ══════════════════════════════════════════════════════════════════════════════

from qgis.core import (QgsProject, QgsRasterLayer, QgsCoordinateTransform,
                       QgsPointXY, QgsFeatureRequest, QgsRectangle, QgsGeometry,
                       QgsMapSettings, QgsMapRendererSequentialJob)
import numpy as np
import math
from qgis.PyQt import sip
from qgis.gui import (QgsMapCanvas, QgsVertexMarker, QgsMapToolPan, QgsRubberBand,
                      QgsMapToolEmitPoint)
from qgis.PyQt.QtCore import Qt, QTimer, QSize
from qgis.PyQt.QtGui import QColor, QImage, QPixmap, QPainter, QPen
from qgis.PyQt.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton,
                                 QComboBox, QLabel, QSpinBox, QCheckBox, QLineEdit, QSlider,
                                 QWidget, QSizePolicy)
try:
    from qgis.core import Qgis
    POLYGON = Qgis.GeometryType.Polygon
except Exception:                                   # older QGIS
    from qgis.core import QgsWkbTypes
    POLYGON = QgsWkbTypes.PolygonGeometry

LAYER_NAME = "hobos_realign"
SCALE_OVERVIEW = 250        # starting zooms — change them in the window
SCALE_CLOSEUP = 40
SCALE_CURSOR = 40
N_ANCHORS = 10            # target anchors per site (change in the window)
BAD_CM = 15               # "Bad anchors" mode: anchors with error above this (cm)


class HoboFinder(QDialog):
    # The layer is looked up again if it was removed/reloaded while the window
    # is open (e.g. after removing a duplicate layer), instead of crashing.
    @property
    def lyr(self):
        L = getattr(self, "_lyr", None)
        try:
            if L is not None and not sip.isdeleted(L):
                return L
        except Exception:
            pass
        found = QgsProject.instance().mapLayersByName(self._layer_name)
        if not found:
            raise Exception(f"No layer called '{self._layer_name}' in the project.")
        self._lyr = found[0]
        try:
            self._lyr.selectionChanged.connect(self.on_selection)
        except Exception:
            pass
        self.fid = None
        print(f"HOBO finder: re-attached to layer '{self._layer_name}'")
        return self._lyr

    @lyr.setter
    def lyr(self, value):
        self._lyr = value

    def __init__(self, iface, layer_name=LAYER_NAME):
        super().__init__(iface.mainWindow())
        self.iface = iface
        self.main = iface.mapCanvas()
        lyrs = QgsProject.instance().mapLayersByName(layer_name)
        if not lyrs:
            raise Exception(f"No layer called '{layer_name}' in the project.")
        self._layer_name = layer_name
        self.lyr = lyrs[0]
        self.markers = []
        self.fid = None
        self.old_pt = None          # logger's original position (canvas CRS)
        self.mouse_pt = None
        self._prog = False          # True while the script itself is zooming

        self.setWindowTitle("HOBO finder")
        self.setWindowFlags(self.windowFlags() | Qt.Tool | Qt.WindowStaysOnTopHint)
        self.resize(420, 1000)

        crs = self.main.mapSettings().destinationCrs()
        self.over_cv, self.s_over, over_box = self.make_view("Old map — overview", SCALE_OVERVIEW, crs)
        self.old_cv, self.s_old, old_box = self.make_view("Old map — close-up", SCALE_CLOSEUP, crs)
        self.new_cv, self.s_new, new_box = self.make_view("New map — your cursor", SCALE_CURSOR, crs)

        self.s_over.valueChanged.connect(lambda _: self.recentre(self.over_cv, self.s_over))
        self.s_old.valueChanged.connect(lambda _: self.recentre(self.old_cv, self.s_old))
        self.s_new.valueChanged.connect(lambda _: self.update_new() if self.mouse_pt
                                        else self.recentre(self.new_cv, self.s_new))
        self.old_cv.extentsChanged.connect(self.update_box)

        # red box on the overview = what the close-up shows
        self.box = QgsRubberBand(self.over_cv, POLYGON)
        self.box.setStrokeColor(QColor("#ff0000"))
        self.box.setFillColor(QColor(255, 0, 0, 30))
        self.box.setWidth(2)

        # crosshair at the centre of the cursor view = your mouse position
        self.cross = QgsVertexMarker(self.new_cv)
        self.cross.setColor(QColor("#ff00ff")); self.cross.setIconSize(22)
        self.cross.setIconType(QgsVertexMarker.ICON_CROSS); self.cross.setPenWidth(2)
        self.cross.hide()

        self.rasters = sorted([l for l in QgsProject.instance().mapLayers().values()
                               if isinstance(l, QgsRasterLayer)],
                              key=lambda l: l.name().lower())
        self.old_combo, self.new_combo = QComboBox(), QComboBox()
        for c in (self.old_combo, self.new_combo):
            for l in self.rasters:
                c.addItem(l.name(), l.id())
        self.old_combo.currentIndexChanged.connect(self.set_layers)
        self.new_combo.currentIndexChanged.connect(self.set_layers)

        self.pan_main = QCheckBox("Pan main map to logger")
        self.pan_main.setChecked(True)

        # click-to-place: one click on the main map = logger goes exactly there
        self.place_tool = QgsMapToolEmitPoint(self.main)
        self.place_tool.canvasClicked.connect(self.place)
        self.prev_tool = None
        self.place_cb = QCheckBox("Click main map to place logger")
        self.next_cb = QCheckBox("then Next")
        self.next_cb.setChecked(True)
        self.place_cb.toggled.connect(self.toggle_place)
        self.main.mapToolSet.connect(self.on_tool_changed)

        self.info = QLabel("Select a logger or press Next ▶")
        prev_b, next_b = QPushButton("◀ Prev"), QPushButton("Next ▶")
        prev_b.clicked.connect(lambda: self.step(-1))
        next_b.clicked.connect(lambda: self.step(1))

        r1 = QHBoxLayout(); r1.addWidget(QLabel("Old map:")); r1.addWidget(self.old_combo, 1)
        r2 = QHBoxLayout(); r2.addWidget(QLabel("New map:")); r2.addWidget(self.new_combo, 1)
        r3 = QHBoxLayout(); r3.addWidget(self.place_cb); r3.addWidget(self.next_cb); r3.addStretch()
        # what Prev/Next steps through
        self.mode = QComboBox()
        self.mode.addItems(["Anchors to add", "Bad anchors (worst first)", "All loggers"])
        self.mode.currentIndexChanged.connect(lambda _: self.on_mode())
        self.n_anchor = QSpinBox(); self.n_anchor.setRange(3, 60); self.n_anchor.setValue(N_ANCHORS)
        self.n_anchor.setPrefix("up to "); self.n_anchor.setSuffix(" per site")
        self.n_anchor.valueChanged.connect(lambda _: self.refresh())
        self.bad_cm = QSpinBox(); self.bad_cm.setRange(1, 200); self.bad_cm.setValue(BAD_CM)
        self.bad_cm.setPrefix("error > "); self.bad_cm.setSuffix(" cm")
        self.bad_cm.valueChanged.connect(lambda _: self.refresh())
        r3b = QHBoxLayout(); r3b.addWidget(QLabel("Step through:")); r3b.addWidget(self.mode)
        r3b.addWidget(self.n_anchor); r3b.addWidget(self.bad_cm); r3b.addStretch()
        # jump straight to any logger by its ID
        self.goto = QLineEdit(); self.goto.setPlaceholderText("logger ID, e.g. H40")
        self.goto.returnPressed.connect(self.go_to)
        go_b = QPushButton("Go"); go_b.clicked.connect(self.go_to)
        r3c = QHBoxLayout(); r3c.addWidget(QLabel("Go to logger:")); r3c.addWidget(self.goto, 1); r3c.addWidget(go_b)
        top = QVBoxLayout(); top.addLayout(r1); top.addLayout(r2); top.addWidget(self.pan_main)
        top.addLayout(r3); top.addLayout(r3b); top.addLayout(r3c)
        views = QVBoxLayout()
        for b in (over_box, old_box, new_box):
            views.addLayout(b, 1)
        views.addLayout(self.make_overlay(), 0)
        bot = QHBoxLayout(); bot.addWidget(prev_b); bot.addWidget(self.info, 1); bot.addWidget(next_b)
        lay = QVBoxLayout(self); lay.addLayout(top); lay.addLayout(views, 1); lay.addLayout(bot)

        # mouse moves fire very often: coalesce them so the view doesn't lag
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.setInterval(60)
        self.timer.timeout.connect(self.update_new)
        # the overlay is heavier: only redraw it once the mouse pauses
        self.ov_timer = QTimer(self)
        self.ov_timer.setSingleShot(True)
        self.ov_timer.setInterval(250)
        self.ov_timer.timeout.connect(self.update_overlay)

        self.lyr.selectionChanged.connect(self.on_selection)
        self.main.xyCoordinates.connect(self.on_mouse)
        self.set_layers()
        if self.lyr.selectedFeatureIds():
            self.on_selection()
        self.place_cb.setChecked(True)
        self.on_mode()

    # ── building blocks ───────────────────────────────────────────────────────
    def make_view(self, title, scale, crs):
        cv = QgsMapCanvas(self)
        cv.setCachingEnabled(True)            # keep rendered layers between redraws
        cv.setParallelRenderingEnabled(True)
        cv.setCanvasColor(QColor("white"))
        cv.setDestinationCrs(crs)
        cv.setMapTool(QgsMapToolPan(cv))
        spin = QSpinBox()
        spin.setRange(2, 20000)
        spin.setPrefix("1:")
        spin.setValue(scale)
        # scrolling inside the view updates its box (without re-centring)
        def sync(s, spin=spin):
            if self._prog or self.old_pt is None:   # ignore start-up noise
                return
            spin.blockSignals(True); spin.setValue(int(round(s))); spin.blockSignals(False)
        cv.scaleChanged.connect(sync)
        box = QVBoxLayout()
        head = QHBoxLayout(); head.addWidget(QLabel(f"<b>{title}</b>")); head.addStretch(); head.addWidget(spin)
        box.addLayout(head); box.addWidget(cv, 1)
        return cv, spin, box

    def to_crs(self, x, y, canvas):
        tr = QgsCoordinateTransform(self.lyr.crs(), canvas.mapSettings().destinationCrs(),
                                    QgsProject.instance())
        return tr.transform(QgsPointXY(x, y))

    def main_to(self, pt, canvas):
        tr = QgsCoordinateTransform(self.main.mapSettings().destinationCrs(),
                                    canvas.mapSettings().destinationCrs(),
                                    QgsProject.instance())
        return tr.transform(pt)

    def old_xy(self, f):
        try:
            return float(f["lon_old"]), float(f["lat_old"])
        except (KeyError, TypeError, ValueError):
            p = f.geometry().asPoint()
            return p.x(), p.y()

    def ordered(self):
        fs = list(self.lyr.getFeatures())     # respects the layer filter
        # skip points with no location (e.g. an accidental "Add Feature")
        empty = [f for f in fs if not f.hasGeometry() or f.geometry().isNull() or f.geometry().isEmpty()]
        if empty and not getattr(self, "_warned_empty", False):
            self._warned_empty = True
            print(f"HOBO finder: skipping {len(empty)} feature(s) with no location, fid(s): "
                  + ", ".join(str(f.id()) for f in empty[:20]))
        fs = [f for f in fs if f not in empty]
        return sorted(fs, key=lambda f: (str(f["rec"]), str(f["year"]), str(f["new_unit"])))

    def clear_markers(self):
        for m, c in self.markers:
            c.scene().removeItem(m)
        self.markers = []

    def marker(self, canvas, pt, colour, size, width, icon=QgsVertexMarker.ICON_CROSS):
        m = QgsVertexMarker(canvas)
        m.setCenter(pt); m.setColor(QColor(colour)); m.setIconSize(size)
        m.setIconType(icon); m.setPenWidth(width)
        self.markers.append((m, canvas))

    def centre_on(self, canvas, pt, scale):
        # explicit extent first: an empty/new canvas can't be re-centred
        self._prog = True
        try:
            canvas.setExtent(QgsRectangle(pt.x() - 1, pt.y() - 1, pt.x() + 1, pt.y() + 1))
            canvas.zoomScale(scale)
            canvas.setCenter(pt)
        finally:
            self._prog = False
        canvas.refresh()

    def recentre(self, canvas, spin):
        if self.old_pt is not None:
            self.centre_on(canvas, self.old_pt, spin.value())

    def update_box(self, *args):
        if self.old_pt is not None:
            self.box.setToGeometry(QgsGeometry.fromRect(self.old_cv.extent()), None)

    def guess_rasters(self, f):
        key = f"{f['rec']}{f['year']}".lower()
        old = next((l for l in self.rasters if l.name().lower() == key), None)
        if old is None:
            return
        self.old_combo.setCurrentIndex(self.rasters.index(old))
        # new map = another raster in the same layer-tree group as the old one
        node = QgsProject.instance().layerTreeRoot().findLayer(old.id())
        if node is not None and node.parent() is not None:
            for child in node.parent().findLayers():
                l = child.layer()
                if isinstance(l, QgsRasterLayer) and l.id() != old.id() and l in self.rasters:
                    self.new_combo.setCurrentIndex(self.rasters.index(l))
                    return

    # ── actions ───────────────────────────────────────────────────────────────
    def set_layers(self, *args):
        prj = QgsProject.instance()
        old = prj.mapLayer(self.old_combo.currentData()) if self.old_combo.count() else None
        new = prj.mapLayer(self.new_combo.currentData()) if self.new_combo.count() else None
        for cv in (self.over_cv, self.old_cv):
            cv.setLayers([old] if old else [])
            cv.refresh()
        self.new_cv.setLayers([self.lyr] + ([new] if new else []))   # points drawn on top
        self.new_cv.refresh()

    def on_selection(self, *args):
        ids = self.lyr.selectedFeatureIds()
        if ids:
            self.fid = ids[-1]
            self.refresh()

    def is_placed(self, f):
        """Already placed by hand: status 'manual', or moved since the last autofill run."""
        names = f.fields().names()
        st = f["status"] if "status" in names else None
        if st == "manual":
            return True
        p = f.geometry().asPoint()
        ox, oy = self.old_xy(f)
        moved = math.hypot(p.x() - ox, p.y() - oy) > 0.001
        if st == "auto" and "auto_x" in names:
            try:   # an auto point counts only if you nudged it since
                return math.hypot(p.x() - float(f["auto_x"]), p.y() - float(f["auto_y"])) > 0.001
            except (TypeError, ValueError):
                return False
        return moved

    def suggested_anchors(self, fs):
        """Loggers still to place so the site reaches the target number of anchors.
        With no anchors yet: 5 around the edge + 1 near the centre, then more.
        With some already placed: each new one is the logger FARTHEST from all
        placed anchors, so the new ones fill the biggest gaps."""
        target = self.n_anchor.value()
        placed = [i for i, f in enumerate(fs) if self.is_placed(f)]
        need = target - len(placed)
        if need <= 0:
            return []
        P = np.array([self.old_xy(f) for f in fs], dtype=float)
        if placed:
            chosen = []
            free = [i for i in range(len(fs)) if i not in placed]
            for _ in range(min(need, len(free))):
                ref = P[placed + chosen]
                i = max(free, key=lambda j: np.min(np.hypot(*(ref - P[j]).T)))
                chosen.append(i); free.remove(i)
            return [fs[i] for i in chosen]
        if len(fs) <= 6:
            return fs
        P = np.array([self.old_xy(f) for f in fs], dtype=float)
        c = P.mean(axis=0)
        D = P - c
        chosen = []
        for k in range(5):
            th = np.radians(90 - k * 72)
            u = np.array([np.cos(th), np.sin(th)])
            along = D @ u
            perp = np.abs(D @ np.array([-u[1], u[0]]))
            score = along - 0.5 * perp            # far out in that direction, near the ray
            for i in np.argsort(-score):
                if i not in chosen:
                    chosen.append(int(i)); break
        rest = [i for i in range(len(fs)) if i not in chosen]
        chosen.append(min(rest, key=lambda i: np.hypot(*D[i])))
        free = [i for i in range(len(fs)) if i not in chosen]
        while len(chosen) < min(target, len(fs)):      # beyond 6: fill the gaps
            i = max(free, key=lambda j: np.min(np.hypot(*(P[chosen] - P[j]).T)))
            chosen.append(i); free.remove(i)
        return [fs[i] for i in chosen[:target]]

    # ── modes ─────────────────────────────────────────────────────────────────
    def adding(self):
        return self.mode.currentIndex() == 0

    def fixing(self):
        return self.mode.currentIndex() == 1

    def on_mode(self):
        self.n_anchor.setVisible(self.adding())
        self.bad_cm.setVisible(self.fixing())
        self.refresh()

    def err(self, f):
        """Leave-one-out error of an anchor (m) from the last autofill run, or None."""
        if "loo_err_m" not in f.fields().names():
            return None
        try:
            v = float(f["loo_err_m"])
            return v if math.isfinite(v) else None
        except (TypeError, ValueError):
            return None

    def replaced_since_autofill(self, f):
        """True if this logger was moved after the last autofill run."""
        names = f.fields().names()
        if "snap_x" not in names:
            return False
        try:
            p = f.geometry().asPoint()
            return math.hypot(p.x() - float(f["snap_x"]), p.y() - float(f["snap_y"])) > 0.001
        except (TypeError, ValueError):
            return False

    def bad_anchors(self, fs):
        """Anchors above the error threshold that you haven't re-placed yet."""
        thr = self.bad_cm.value() / 100
        bad = [f for f in fs if self.err(f) is not None and self.err(f) > thr
               and not self.replaced_since_autofill(f)]
        return sorted(bad, key=lambda f: -self.err(f))

    def go_to(self):
        want = self.goto.text().strip().lower()
        if not want:
            return
        hit = [f for f in self.ordered() if str(f["new_unit"]).lower() == want]
        if not hit:
            self.info.setText(f"{self.goto.text().strip()} is not in the current filter")
            return
        self.lyr.selectByIds([hit[0].id()])   # → on_selection → refresh

    def sequence(self):
        fs = self.ordered()
        if self.fixing():
            return self.bad_anchors(fs)
        if not self.adding():
            return fs
        seq, groups = [], {}
        for f in fs:
            groups.setdefault((str(f["rec"]), str(f["year"])), []).append(f)
        for key in sorted(groups):
            seq += self.suggested_anchors(groups[key])
        return seq

    def step(self, d):
        fs = self.sequence()
        if not fs:
            if not self.ordered():
                msg = "No loggers in the current filter"
            elif self.adding():
                msg = "All sites in the filter have their anchors — run the autofill"
            else:
                fixed = sum(self.replaced_since_autofill(f) for f in self.ordered())
                msg = (f"Done here: {fixed} re-placed. Save + run the autofill to update the errors, "
                       "then change the filter to the next site"
                       if fixed else
                       f"No anchors with error > {self.bad_cm.value()} cm in the filter. "
                       "Change the filter to the next site")
            self.info.setText(msg)
            return
        ids = [f.id() for f in fs]
        i = ids.index(self.fid) + d if self.fid in ids else (0 if d > 0 else len(ids) - 1)
        self.lyr.selectByIds([ids[max(0, min(i, len(ids) - 1))]])   # → on_selection

    def on_mouse(self, pt):
        self.mouse_pt = pt
        if not self.timer.isActive():
            self.timer.start()

    def update_new(self):
        if self.mouse_pt is None:
            return
        c = self.main_to(self.mouse_pt, self.new_cv)
        cv, want = self.new_cv, self.s_new.value()
        if cv.extent().width() > 0 and abs(cv.scale() - want) <= 0.01 * want:
            cv.setCenter(c)                   # fast path: one redraw, zoom unchanged
        else:
            self.centre_on(cv, c, want)
        self.cross.setCenter(c)
        self.cross.show()
        if self.ov_on.isChecked():
            self.ov_timer.start()             # restarts on every move -> draws on pause

    # ── overlay ───────────────────────────────────────────────────────────────
    def make_overlay(self):
        self.ov_on = QCheckBox("Overlay")
        self.ov_mode = QComboBox(); self.ov_mode.addItems(["Red/cyan", "Blend"])
        self.ov_alpha = QSlider(Qt.Horizontal); self.ov_alpha.setRange(0, 100); self.ov_alpha.setValue(50)
        self.ov_alpha.setToolTip("New-map opacity (Blend mode)")
        self.s_ov = QSpinBox(); self.s_ov.setRange(2, 20000); self.s_ov.setPrefix("1:")
        self.s_ov.setValue(SCALE_CLOSEUP)
        self.ov_img = QLabel(); self.ov_img.setAlignment(Qt.AlignCenter)
        self.ov_img.setMinimumHeight(200)
        self.ov_img.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored)
        self.ov_img.setStyleSheet("background: white;")
        self.ov_old = None            # cached render of the old map at the logger
        self.ov_old_key = None

        head = QHBoxLayout()
        head.addWidget(self.ov_on); head.addWidget(self.ov_mode); head.addWidget(self.ov_alpha, 1)
        head.addWidget(self.s_ov)
        self.ov_box = QWidget(); vb = QVBoxLayout(self.ov_box); vb.setContentsMargins(0, 0, 0, 0)
        vb.addWidget(self.ov_img, 1)
        self.ov_box.setVisible(False)
        outer = QVBoxLayout(); outer.addLayout(head); outer.addWidget(self.ov_box, 1)

        self.ov_on.toggled.connect(self.toggle_overlay)
        for sig in (self.ov_mode.currentIndexChanged, self.ov_alpha.valueChanged,
                    self.s_ov.valueChanged):
            sig.connect(lambda *_: self.update_overlay())
        return outer

    def toggle_overlay(self, on):
        self.ov_box.setVisible(on)
        if on:
            self.resize(self.width(), max(self.height(), 1150))
        self.update_overlay()

    def render(self, layer, centre, scale, w, h):
        ms = QgsMapSettings()
        ms.setLayers([layer])
        ms.setDestinationCrs(self.old_cv.mapSettings().destinationCrs())
        ms.setOutputSize(QSize(w, h))
        ms.setOutputDpi(96)
        ms.setBackgroundColor(QColor("white"))
        half_w = scale * w / (96 / 0.0254) / 2          # metres per half-width
        half_h = half_w * h / w
        ms.setExtent(QgsRectangle(centre.x() - half_w, centre.y() - half_h,
                                  centre.x() + half_w, centre.y() + half_h))
        job = QgsMapRendererSequentialJob(ms)
        job.start(); job.waitForFinished()
        return job.renderedImage()

    @staticmethod
    def grey(img):
        img = img.convertToFormat(QImage.Format_RGB32)
        w, h = img.width(), img.height()
        ptr = img.constBits(); ptr.setsize(img.sizeInBytes() if hasattr(img, "sizeInBytes") else img.byteCount())
        a = np.frombuffer(ptr, np.uint8).reshape(h, img.bytesPerLine() // 4, 4)[:, :w, :3]
        return a.astype(np.float32).mean(axis=2)

    def update_overlay(self):
        if not self.ov_on.isChecked() or self.old_pt is None:
            return
        prj = QgsProject.instance()
        old = prj.mapLayer(self.old_combo.currentData()) if self.old_combo.count() else None
        new = prj.mapLayer(self.new_combo.currentData()) if self.new_combo.count() else None
        if old is None or new is None:
            return
        w, h = max(self.ov_img.width(), 50), max(self.ov_img.height(), 50)
        scale = self.s_ov.value()
        key = (old.id(), self.old_pt.x(), self.old_pt.y(), scale, w, h)
        if key != self.ov_old_key:                    # old view only changes per logger
            self.ov_old = self.grey(self.render(old, self.old_pt, scale, w, h))
            self.ov_old_key = key
        c = self.main_to(self.mouse_pt, self.old_cv) if self.mouse_pt else self.old_pt
        g_new = self.grey(self.render(new, c, scale, w, h))
        g_old = self.ov_old

        if self.ov_mode.currentText() == "Blend":
            a = self.ov_alpha.value() / 100
            mix = (1 - a) * g_old + a * g_new
            rgb = np.dstack([mix, mix, mix])
        else:                                         # red = old, cyan = new
            rgb = np.dstack([g_old, g_new, g_new])
        rgb = np.ascontiguousarray(np.clip(rgb, 0, 255).astype(np.uint8))
        img = QImage(rgb.data, w, h, 3 * w, QImage.Format_RGB888).copy()

        p = QPainter(img)
        p.setPen(QPen(QColor("#ff00ff"), 2))
        cx, cy = w // 2, h // 2
        p.drawLine(cx - 12, cy, cx + 12, cy); p.drawLine(cx, cy - 12, cx, cy + 12)
        p.drawEllipse(cx - 15, cy - 15, 30, 30)
        p.end()
        self.ov_img.setPixmap(QPixmap.fromImage(img))

    def resizeEvent(self, e):
        super().resizeEvent(e)
        self.update_overlay()

    def refresh(self):
        if self.fid is None:
            return
        f = next(self.lyr.getFeatures(QgsFeatureRequest(self.fid)), None)
        if f is None:
            return
        self.guess_rasters(f)
        self.clear_markers()
        ox, oy = self.old_xy(f)
        self.old_pt = self.to_crs(ox, oy, self.old_cv)

        # yellow: the site's other loggers (overview + close-up), red: this one
        for g in self.ordered():
            if g.id() != f.id() and g["rec"] == f["rec"] and g["year"] == f["year"]:
                gx, gy = self.old_xy(g)
                for cv in (self.over_cv, self.old_cv):
                    self.marker(cv, self.to_crs(gx, gy, cv), "#ffd400", 10, 2)
        site = [g for g in self.ordered() if g["rec"] == f["rec"] and g["year"] == f["year"]]
        if self.fixing():                        # orange circles: bad anchors at this site
            for g in self.bad_anchors(site):
                gx, gy = self.old_xy(g)
                self.marker(self.over_cv, self.to_crs(gx, gy, self.over_cv), "#ff6d00", 16, 2,
                            QgsVertexMarker.ICON_CIRCLE)
        if self.adding():                        # green circles: suggested anchors
            for g in self.suggested_anchors(site):
                gx, gy = self.old_xy(g)
                self.marker(self.over_cv, self.to_crs(gx, gy, self.over_cv), "#00c853", 16, 2,
                            QgsVertexMarker.ICON_CIRCLE)
        self.marker(self.over_cv, self.old_pt, "#ff0000", 14, 3)
        self.marker(self.old_cv, self.old_pt, "#ff0000", 18, 3)
        self.marker(self.old_cv, self.old_pt, "#ff0000", 30, 2, QgsVertexMarker.ICON_CIRCLE)

        self.centre_on(self.old_cv, self.old_pt, self.s_old.value())
        self.centre_on(self.over_cv, self.old_pt, self.s_over.value())
        self.update_box()

        # cyan cross: the logger's original spot, on the cursor view and main map
        self.marker(self.new_cv, self.to_crs(ox, oy, self.new_cv), "#00e5ff", 14, 2)
        m = self.to_crs(ox, oy, self.main)
        self.marker(self.main, m, "#00e5ff", 18, 3)
        if self.pan_main.isChecked():
            self.main.setCenter(m)            # pan only, keep your zoom
        self.main.refresh()
        if self.mouse_pt is None:             # no mouse yet: show the same spot
            self.centre_on(self.new_cv, self.to_crs(ox, oy, self.new_cv), self.s_new.value())

        fs = self.sequence()
        ids = [g.id() for g in fs]
        tag = "to place " if self.adding() else ("bad " if self.fixing() else "")
        pos = f"{tag}{ids.index(f.id()) + 1}/{len(ids)}" if f.id() in ids else ""
        st = f["status"] if "status" in f.fields().names() and f["status"] else ""
        e = self.err(f)
        es = f"  error {e * 100:.0f} cm" if e is not None else ""
        if e is not None and self.replaced_since_autofill(f):
            es += " (re-placed — run the autofill)"
        self.info.setText(f"<b>{f['new_unit']}</b>  {f['rec']}{f['year']}  {pos}  {st}{es}")
        self.update_overlay()

    # ── click to place ────────────────────────────────────────────────────────
    def toggle_place(self, on):
        if on:
            cur = self.main.mapTool()
            if cur is not self.place_tool:
                self.prev_tool = cur
            self.main.setMapTool(self.place_tool)
        elif self.main.mapTool() is self.place_tool:
            if self.prev_tool is not None:
                self.main.setMapTool(self.prev_tool)
            else:
                self.main.unsetMapTool(self.place_tool)

    def on_tool_changed(self, new, old=None):
        # you picked another QGIS tool → untick, so the box shows the truth
        if new is not self.place_tool and self.place_cb.isChecked():
            self.place_cb.blockSignals(True)
            self.place_cb.setChecked(False)
            self.place_cb.blockSignals(False)

    def place(self, pt, button):
        if button != Qt.LeftButton:
            return
        if self.fid is None:
            self.info.setText("Select a logger first (Next ▶)")
            return
        tr = QgsCoordinateTransform(self.main.mapSettings().destinationCrs(),
                                    self.lyr.crs(), QgsProject.instance())
        p = tr.transform(pt)
        if not self.lyr.isEditable():
            self.lyr.startEditing()
        self.lyr.beginEditCommand("Place logger")
        self.lyr.changeGeometry(self.fid, QgsGeometry.fromPointXY(QgsPointXY(p)))
        self.lyr.endEditCommand()
        self.lyr.triggerRepaint()
        self.new_cv.refresh()
        if self.next_cb.isChecked():
            self.step(1)

    def closeEvent(self, e):
        self.timer.stop()
        self.ov_timer.stop()
        self.place_cb.setChecked(False)       # give the main map its tool back
        for sig, slot in [(self.lyr.selectionChanged, self.on_selection),
                          (self.main.xyCoordinates, self.on_mouse),
                          (self.main.mapToolSet, self.on_tool_changed)]:
            try:
                sig.disconnect(slot)
            except Exception:
                pass
        self.clear_markers()
        self.main.refresh()
        super().closeEvent(e)


try:
    _hobo_finder.close()
except Exception:
    pass
_hobo_finder = HoboFinder(iface)
_hobo_finder.show()
