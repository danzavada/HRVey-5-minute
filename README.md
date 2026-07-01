# HRVey 5-min — R-peak detection for 5-minute ECG PDFs

A browser port of the **ECG Peak Detector** desktop app (`MonkECG-5`). It opens a
5-minute vector ECG PDF, **auto-extracts the waveform straight from the PDF's vector
paths**, **auto-detects R-peaks** (Pan–Tompkins), lets you edit and annotate beats, and
**exports R-R intervals to Excel** — all **entirely in your browser**, offline. No upload,
no server, no install.

This is a faithful web reimplementation of the Python/PyQt desktop tool. The signal
pipeline is a line-for-line port of the desktop app's numpy/scipy code and has been
validated to reproduce it to within the PDF format's floating-point floor (see
[Fidelity](#fidelity)).

---

## What it does

1. **Open** a 5-minute ECG PDF. The ECG lives on **page 2** as thin vector strokes
   (14 stacked strips ≈ 22 s each).
2. It renders page 2 and, in parallel, reads the raw vector paths, reconstructs the
   continuous voltage signal, resamples it to 2000 Hz, and runs R-peak detection.
3. Peaks are drawn on the trace, colour-coded, with R-R interval labels when you zoom in.
4. You **review and correct**: add missed peaks, delete false ones, annotate ectopics.
5. **Export** the R-R intervals (ms) + annotations to an `.xlsx` file.

### Beat colours

| Colour | Meaning |
|--------|---------|
| 🟢 Green | Normal beat |
| 🔵 Blue | Long R-R (> 20 % above the median) |
| 🟠 Orange | Short R-R (> 20 % below the median) |
| 🔴 Red | PVC / PAC (manually tagged) |
| ⚪ Gray | "Can't detect accurately" (manually tagged) |

---

## Usage

Open `index.html` — either double-click it, or serve the folder and browse to it:

```bash
# any static server works; e.g.
python -m http.server 8765
# then open http://localhost:8765/
```

Then **Open PDF** (or drag-and-drop, or `Ctrl+O`) and pick your ECG PDF. There's a
ready-made synthetic example in [`EXAMPLE/`](EXAMPLE/) with its expected `.xlsx` output.

### Modes & shortcuts

| Key | Action |
|-----|--------|
| `S` | **Select** — click a beat to select it (then annotate) |
| `A` | **Add peak** — click on a beat; snaps to the local max. `Shift+click` = exact |
| `D` | **Delete peak** — click a peak to remove it |
| `P` / `Q` / `U` | Tag the selected beat **PVC** / **PAC** / **can't-detect** (toggles) |
| `Delete` | Remove the selected peak |
| `Ctrl+Z` | Undo |
| `Z` | Fit page to window |
| Drag | Pan · `Ctrl`+wheel | Zoom · wheel | Pan |

Mode keys and detection **sensitivity** (1–10) are configurable in **⚙ Settings**.
Sensitivity re-runs detection; a higher value finds more (smaller) peaks.

---

## Fidelity

The hard part of the port is reproducing the desktop app's exact signal processing in the
browser. Both the PDF vector extraction and the DSP were validated numerically against the
original Python (`scipy`/`numpy`) implementation:

- **PDF path extraction** — `pdf.js` operator lists are transformed through the CTM +
  viewport to reproduce PyMuPDF's `get_drawings()` coordinates and stroke widths **exactly**
  (pixel-identical point sets).
- **Butterworth band-pass + `filtfilt`, `uniform_filter1d`, `find_peaks`
  (height / distance / prominence)** — reimplemented from scratch and checked **bit-for-bit**
  against `scipy` fixtures (Δ = 0).
- **Full pipeline** on a synthetic 5-minute ECG: **identical R-peak count**, peak times
  within one sample (0.5 ms), voltage within 0.2 µV. The residual is the PDF format's
  inherent float32 coordinate storage — a hard floor, not a divergence.

See [`ecg-core.js`](ecg-core.js) for the ported pipeline (pure, dependency-free, runs in the
browser and in Node).

---

## The expected PDF format

- **Two pages**; the ECG signal is on **page 2** (page index 1).
- The trace is drawn as **vector strokes 0.4–0.8 pt wide** (thicker/thinner lines such as
  grid lines are ignored).
- **14 strips** stacked vertically with blank gaps between them; paper speed **12.5 mm/s**,
  gain **2.5 mm/mV**. Strips are stitched back into one continuous signal.

Scanned/raster ECGs (no vector paths) are **not** supported by this tool — use the
image-based sibling app for those.

---

## Tech stack

- **Vanilla JS + Canvas** — no framework.
- [**pdf.js**](https://mozilla.github.io/pdf.js/) 3.32.2 — PDF rendering + vector path extraction (vendored).
- [**SheetJS**](https://sheetjs.com/) — `.xlsx` export (vendored).
- PWA (installable, offline via service worker).

Everything runs client-side; your ECG never leaves the machine.

---

## Project layout

```
index.html                  app shell + UI + all interaction logic
ecg-core.js                 ported signal pipeline (extraction + DSP + detection)
manifest.webmanifest, sw.js PWA manifest + offline service worker
vendor/                     pdf.js + SheetJS (pinned)
icons/                      app icons
EXAMPLE/                    synthetic 5-minute ECG PDF + its exported .xlsx
```

Ported from the `MonkECG-5` desktop app (Python / PyQt6 / pyqtgraph / PyMuPDF / scipy).
