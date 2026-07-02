# HRVey 5-min

R-peak detection for **5-minute vector ECG PDFs**, entirely in the browser — offline, no upload, no install.

## Use

Open `index.html` (or serve the folder and browse to it), then **Open PDF** (or drag-and-drop). The app extracts the ECG from the PDF's vector paths, detects R-peaks automatically, and lets you review, edit, and export.

- **Sensitivity** slider (toolbar) re-detects live.
- **Select / Add / Delete** peaks (`S` / `A` / `D`); `Shift+click` adds exactly; annotate a selected beat **PVC** / **PAC** (`P` / `Q`); `Ctrl+Z` undo.
- **Export .pdf / .xlsx** (right panel), named `<Subject ID>_YYYY_MM_DD` (date read from the PDF). The `.xlsx` includes an **Analysis** tab; **Load analysis** restores those beats onto the matching PDF.
- **`N / 14 rows`** badge + **Rows** overlay flag any strip the extractor misses.

A synthetic example ECG is in [`EXAMPLE/`](EXAMPLE/).

## Expected PDF

Two pages, ECG on **page 2** as thin (0.4–0.8 pt) vector strokes in 14 stacked strips (12.5 mm/s, 2.5 mm/mV). Scanned/raster ECGs are not supported.

## Stack

Vanilla JS + Canvas · [pdf.js](https://mozilla.github.io/pdf.js/) (vector extraction) · [SheetJS](https://sheetjs.com/) (xlsx) · PWA. Everything runs client-side. Pipeline: [`ecg-core.js`](ecg-core.js).
