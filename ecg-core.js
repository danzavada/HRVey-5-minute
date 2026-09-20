/* ecg-core.js — faithful JS port of ecg_viewer.py's signal pipeline.
 *
 * Pure, dependency-free (no DOM, no pdf.js). Works in the browser and in Node.
 * The only PDF-specific piece is extractStrokePaths(), which consumes an
 * already-decoded pdf.js operator list (fnArray/argsArray) + viewport transform.
 *
 * Every numeric routine mirrors the exact scipy/numpy behaviour the desktop app
 * relies on (butter/filtfilt, uniform_filter1d reflect alignment, find_peaks
 * height/distance/prominence, interp1d linear). Validated against scipy fixtures
 * and the full Python reference — see test.js.
 */
(function (root) {
  "use strict";

  // ── Physical constants (verbatim from ecg_viewer.py) ──────────────────── //
  // PAPER_SPEED / VOLT_SCALE are the defaults for the 5-minute report; a
  // document may override them (e.g. the 1-minute report prints at
  // 25 mm/s 10 mm/mV) — see parseScale() and extractEcg(paths, opts).
  const PAPER_SPEED = 12.5;
  const VOLT_SCALE  = 2.5;
  const PT_PER_MM   = 2.834645669291339;
  const N_ROWS      = 14;          // strips in the 5-minute (2-page) report
  const FS          = 2000;
  const S_PER_PT    = 1.0 / (PAPER_SPEED * PT_PER_MM);
  const MV_PER_PT   = 1.0 / (VOLT_SCALE  * PT_PER_MM);
  const ABNORMAL_THRESHOLD = 0.2;

  // Read "12.5mm/s 2.5mm/mV" / "25mm/s 10mm/mV" from the page text. Returns
  // nulls for whatever is not printed so the caller can fall back to defaults.
  function parseScale(text) {
    const t = String(text || "");
    const sp = /(\d+(?:[.,]\d+)?)\s*mm\s*\/\s*s(?![a-z])/i.exec(t);
    const vs = /(\d+(?:[.,]\d+)?)\s*mm\s*\/\s*mv/i.exec(t);
    const num = (m) => (m ? parseFloat(m[1].replace(",", ".")) : NaN);
    const a = num(sp), b = num(vs);
    return { paperSpeed: a > 0 ? a : null, voltScale: b > 0 ? b : null };
  }

  // Fixed Butterworth band-pass (order 2, 5–15 Hz @ 2000 Hz) + filtfilt zi.
  // Precomputed with scipy so results are bit-comparable to the desktop app.
  const B  = [0.00024135904904198078, 0.0, -0.00048271809808396156, 0.0, 0.00024135904904198078];
  const A  = [1.0, -3.954114210492194, 5.86480480132174, -3.867233731430246, 0.9565436765112041];
  const ZI = [-0.00024135904902100982, -0.0002413590491039314, 0.0002413590491030208, 0.00024135904902192114];
  const PADLEN = 15; // scipy filtfilt default: 3 * max(len(a), len(b)) = 3*5

  // ── small numeric helpers ─────────────────────────────────────────────── //
  function median(arr) {
    const n = arr.length;
    if (n === 0) return 0;
    const s = Float64Array.from(arr).sort();
    const m = n >> 1;
    return (n & 1) ? s[m] : (s[m - 1] + s[m]) / 2;
  }
  function mean(arr) { let s = 0; for (let i = 0; i < arr.length; i++) s += arr[i]; return s / arr.length; }
  function std(arr) { const m = mean(arr); let s = 0; for (let i = 0; i < arr.length; i++) { const d = arr[i] - m; s += d * d; } return Math.sqrt(s / arr.length); }
  function argmaxRange(arr, lo, hi) { // inclusive hi; returns absolute index, first max (numpy argmax)
    let bi = lo, bv = arr[lo];
    for (let i = lo + 1; i <= hi; i++) { if (arr[i] > bv) { bv = arr[i]; bi = i; } }
    return bi;
  }
  // numpy searchsorted(a, v) side='left': first i with a[i] >= v
  function searchsortedLeft(a, v) {
    let lo = 0, hi = a.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (a[mid] < v) lo = mid + 1; else hi = mid; }
    return lo;
  }

  // scipy.ndimage.uniform_filter1d, mode='reflect' (half-sample: d c b a | a b c d | d c b a),
  // default origin=0 → window [i - size//2 , i - size//2 + size - 1].
  function uniformFilter1d(x, size) {
    const n = x.length, out = new Float64Array(n), lo = size >> 1;
    const refl = (j) => { // reflect index into [0, n-1]
      if (n === 1) return 0;
      while (j < 0 || j >= n) { if (j < 0) j = -j - 1; else if (j >= n) j = 2 * n - 1 - j; }
      return j;
    };
    // prefix sums over a reflected virtual array via direct windowed sum with running update
    let acc = 0;
    for (let k = 0; k < size; k++) acc += x[refl(0 - lo + k)];
    out[0] = acc / size;
    for (let i = 1; i < n; i++) {
      acc -= x[refl(i - 1 - lo)];
      acc += x[refl(i - lo + size - 1)];
      out[i] = acc / size;
    }
    return out;
  }

  // scipy.signal.lfilter, Direct Form II Transposed; a[0] assumed 1, len(a)==len(b).
  function lfilter(b, a, x, zi) {
    const n = x.length, na = a.length, y = new Float64Array(n);
    const z = zi ? zi.slice() : new Array(na - 1).fill(0);
    for (let i = 0; i < n; i++) {
      const xi = x[i];
      const yi = b[0] * xi + z[0];
      y[i] = yi;
      for (let k = 1; k < na - 1; k++) z[k - 1] = b[k] * xi + z[k] - a[k] * yi;
      z[na - 2] = b[na - 1] * xi - a[na - 1] * yi;
    }
    return y;
  }
  // scipy.signal.filtfilt, method='pad', padtype='odd', default padlen.
  function filtfilt(b, a, x) {
    const n = x.length;
    if (n <= PADLEN) return filtfiltShort(b, a, x); // guard tiny inputs
    const ext = new Float64Array(n + 2 * PADLEN);
    for (let i = 0; i < PADLEN; i++) ext[i] = 2 * x[0] - x[PADLEN - i];              // odd left
    ext.set(x, PADLEN);
    for (let i = 0; i < PADLEN; i++) ext[n + PADLEN + i] = 2 * x[n - 1] - x[n - 2 - i]; // odd right
    const scale = (v) => ZI.map((z) => z * v);
    let y = lfilter(b, a, ext, scale(ext[0]));
    y.reverse();
    y = lfilter(b, a, y, scale(y[0]));
    y.reverse();
    return y.subarray(PADLEN, PADLEN + n);
  }
  function filtfiltShort(b, a, x) { // rare fallback: no padding
    const scale = (v) => ZI.map((z) => z * v);
    let y = lfilter(b, a, x, scale(x[0])); y.reverse();
    y = lfilter(b, a, y, scale(y[0])); y.reverse();
    return y;
  }

  // scipy.signal._local_maxima_1d: strict local maxima, plateau → midpoint.
  function localMaxima1d(x) {
    const n = x.length, out = [];
    let i = 1;
    const iMax = n - 1;
    while (i < iMax) {
      if (x[i - 1] < x[i]) {
        let iAhead = i + 1;
        while (iAhead < iMax && x[iAhead] === x[i]) iAhead++;
        if (x[iAhead] < x[i]) { out.push((i + iAhead - 1) >> 1); i = iAhead; }
        else i = iAhead;
      } else i++;
    }
    return out;
  }
  // scipy.signal._peak_prominences (wlen=None → full search).
  function peakProminences(x, peaks) {
    const n = x.length, prom = new Float64Array(peaks.length);
    for (let pi = 0; pi < peaks.length; pi++) {
      const p = peaks[pi], h = x[p];
      let i = p, leftMin = h;
      while (i >= 0 && x[i] <= h) { if (x[i] < leftMin) leftMin = x[i]; i--; }
      i = p; let rightMin = h;
      while (i < n && x[i] <= h) { if (x[i] < rightMin) rightMin = x[i]; i++; }
      prom[pi] = h - Math.max(leftMin, rightMin);
    }
    return prom;
  }
  // scipy.signal._select_by_peak_distance
  function selectByDistance(peaks, priority, distance) {
    const np_ = peaks.length, keep = new Uint8Array(np_).fill(1);
    const dist = Math.ceil(distance);
    const order = Array.from(peaks.keys()).sort((a, b) => priority[a] - priority[b]); // ascending priority
    for (let idx = np_ - 1; idx >= 0; idx--) {
      const j = order[idx];
      if (!keep[j]) continue;
      let k = j - 1;
      while (k >= 0 && peaks[j] - peaks[k] < dist) { keep[k] = 0; k--; }
      k = j + 1;
      while (k < np_ && peaks[k] - peaks[j] < dist) { keep[k] = 0; k++; }
    }
    const res = [];
    for (let m = 0; m < np_; m++) if (keep[m]) res.push(peaks[m]);
    return res;
  }
  // scipy.signal.find_peaks — supports {height, distance, prominence} (order: height→distance→prominence).
  function findPeaks(x, opts) {
    opts = opts || {};
    let peaks = localMaxima1d(x);
    if (opts.height != null) peaks = peaks.filter((p) => x[p] >= opts.height);
    if (opts.distance != null && peaks.length) {
      const prio = peaks.map((p) => x[p]);
      peaks = selectByDistance(peaks, prio, opts.distance);
    }
    if (opts.prominence != null && peaks.length) {
      const prom = peakProminences(x, peaks);
      peaks = peaks.filter((_, i) => prom[i] >= opts.prominence);
    }
    return peaks;
  }

  // numpy.histogram(a, bins=n) over [min,max]; last bin closed on the right.
  function histogram(a, nbins) {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < a.length; i++) { const v = a[i]; if (v < lo) lo = v; if (v > hi) hi = v; }
    const counts = new Float64Array(nbins), width = (hi - lo) / nbins;
    for (let i = 0; i < a.length; i++) {
      let b = width > 0 ? Math.floor((a[i] - lo) / width) : 0;
      if (b >= nbins) b = nbins - 1; if (b < 0) b = 0;
      counts[b]++;
    }
    const edges = new Float64Array(nbins + 1);
    for (let i = 0; i <= nbins; i++) edges[i] = lo + i * width;
    return { counts, edges };
  }

  // ── strip-boundary detection ──────────────────────────────────────────── //
  // Strips are separated by horizontal bands that no trace enters. Build a
  // 1-pt histogram of every stroke vertex's y, cut it at runs of ≥ MIN_GAP
  // (near-)empty bins, then fold fragments back into their neighbours: a
  // fragment that is too light to be a strip (an artifact plateau carved off
  // its own trace), or one sitting closer than half the strip pitch to its
  // neighbour (a plateau cannot lie between two strips that ARE a pitch
  // apart). The number of strips is whatever survives — 14 for the 5-minute
  // report, 4 or 7 for the 1-minute reports (3/6 full strips + a 1.2 s stub).
  // The previous density-valley approach needed a fixed strip count and lost
  // the stub: with ~2 % of the vertices it never produced a deep enough valley.
  // Returns the boundary y's (ascending, in page points), i.e. n_rows - 1 values.
  function findStripBoundaries(pathYs) {
    let total = 0, yLo = Infinity, yHi = -Infinity;
    for (const ys of pathYs) for (let i = 0; i < ys.length; i++) {
      const v = ys[i]; if (v < yLo) yLo = v; if (v > yHi) yHi = v; total++;
    }
    if (!(yHi > yLo)) return [];
    const nBins = Math.ceil(yHi - yLo) + 1;
    const counts = new Float64Array(nBins);
    for (const ys of pathYs) for (let i = 0; i < ys.length; i++) counts[Math.floor(ys[i] - yLo)]++;

    const EPS = Math.max(1, total * 1e-4); // a bin this sparse is an artifact trail, not a strip
    const MIN_GAP = 3;                     // bins; a steep upstroke may skip 1–2 bins of its own strip
    const MIN_MASS = 0.01;                 // fraction of all vertices a strip must hold (a stub has ~2 %)

    // maximal runs of non-empty bins, split wherever ≥ MIN_GAP empty bins intervene
    const segs = [];                       // {lo, hi}: inclusive bin range
    let lo = -1, last = -1;
    for (let b = 0; b < nBins; b++) {
      if (counts[b] <= EPS) continue;
      if (lo < 0 || b - last - 1 >= MIN_GAP) { if (lo >= 0) segs.push({ lo, hi: last }); lo = b; }
      last = b;
    }
    if (lo >= 0) segs.push({ lo, hi: last });
    const weigh = (s) => { let m = 0, c = 0; for (let b = s.lo; b <= s.hi; b++) { m += counts[b]; c += counts[b] * (b + 0.5); } s.mass = m; s.c = m ? c / m : (s.lo + s.hi) / 2; return s; };
    segs.forEach(weigh);
    const fold = (i, into) => { const a = segs[i], b = segs[into]; b.lo = Math.min(a.lo, b.lo); b.hi = Math.max(a.hi, b.hi); weigh(b); segs.splice(i, 1); };
    const nearer = (i) => (i === 0 ? 1 : i === segs.length - 1 ? i - 1
      : (segs[i].c - segs[i - 1].c <= segs[i + 1].c - segs[i].c ? i - 1 : i + 1));

    while (segs.length > 1) {
      let li = 0; for (let i = 1; i < segs.length; i++) if (segs[i].mass < segs[li].mass) li = i;
      if (segs[li].mass < MIN_MASS * total) { fold(li, nearer(li)); continue; }
      if (segs.length < 3) break;
      const gaps = []; for (let i = 1; i < segs.length; i++) gaps.push(segs[i].c - segs[i - 1].c);
      const pitch = median(gaps);
      let gi = 0; for (let i = 1; i < gaps.length; i++) if (gaps[i] < gaps[gi]) gi = i;
      if (gaps[gi] < 0.5 * pitch) { const a = gi, b = gi + 1; if (segs[a].mass < segs[b].mass) fold(a, b); else fold(b, a); continue; }
      break;
    }
    const out = [];
    for (let i = 1; i < segs.length; i++) out.push(yLo + (segs[i - 1].hi + 1 + segs[i].lo) / 2);
    return out;
  }

  // ── main extraction (mirrors extract_ecg) ─────────────────────────────── //
  // ECG-candidate paths: thin strokes (the trace is 0.567 pt; grid lines are
  // fills or 0/1.134 pt) with more than one line SEGMENT (≥3 vertices), matching
  // PyMuPDF's `len(d['items']) > 1`. Tiny 2-vertex stub paths carry outlier
  // y-values that mis-assigned whole ECG rows on real device PDFs. See [[row-detection-bug]].
  function isEcgPath(p) { return p.width > 0.4 && p.width < 0.8 && p.xs.length > 2; }
  // Vertex count of the ECG-candidate paths — used to pick the page that carries
  // the trace (a 2-page report's summary page has a few hundred from the R-R graph).
  function countEcgVertices(paths) {
    let n = 0; for (const p of paths) if (isEcgPath(p)) n += p.xs.length; return n;
  }

  // paths: [{width:Number, xs:Float64Array, ys:Float64Array}] — all stroked paths.
  // opts:  {paperSpeed (mm/s), voltScale (mm/mV)} — defaults 12.5 / 2.5.
  // The strip count is detected, not assumed: 14 for the 5-minute report, 4 or
  // 7 for the 1-minute ones. Row indices in raw_r / cum_t / baseline_y run
  // 0..n_rows-1 top to bottom.
  function extractEcg(paths, opts) {
    opts = opts || {};
    const paperSpeed = opts.paperSpeed > 0 ? opts.paperSpeed : PAPER_SPEED;
    const voltScale  = opts.voltScale  > 0 ? opts.voltScale  : VOLT_SCALE;
    const sPerPt = 1.0 / (paperSpeed * PT_PER_MM), mvPerPt = 1.0 / (voltScale * PT_PER_MM);
    const empty = { time: new Float64Array(0), voltage: new Float64Array(0),
      cum_t: new Float64Array(0), baseline_y: new Float64Array(0), n_rows: 0,
      x_left: 0, y_min: 0, band: 1, raw_r: new Int8Array(0),
      paper_speed: paperSpeed, volt_scale: voltScale, s_per_pt: sPerPt, mv_per_pt: mvPerPt };

    const pd = [];
    for (const p of paths) {
      if (!isEcgPath(p)) continue;
      pd.push({ xs: p.xs, ys: p.ys, med: median(p.ys) });
    }
    if (!pd.length) return empty;

    pd.sort((a, b) => a.med - b.med);
    let xLeft = Infinity;
    for (const p of pd) for (let i = 0; i < p.xs.length; i++) if (p.xs[i] < xLeft) xLeft = p.xs[i];
    const meds = pd.map((p) => p.med);
    const yMin = meds[0], yMax = meds[meds.length - 1];

    const boundaries = findStripBoundaries(pd.map((p) => p.ys));
    const nRows = boundaries.length + 1;
    const band = nRows > 1 ? (yMax - yMin) / (nRows - 1) : 1.0;
    const rows = meds.map((m) => searchsortedLeft(boundaries, m));

    const rowXs = Array.from({ length: nRows }, () => []);
    const rowYs = Array.from({ length: nRows }, () => []);
    for (let i = 0; i < pd.length; i++) { rowXs[rows[i]].push(pd[i].xs); rowYs[rows[i]].push(pd[i].ys); }

    const xEnds = new Float64Array(nRows), baseYs = new Float64Array(nRows), cumT = new Float64Array(nRows);
    const tParts = [], vParts = [], rParts = [];
    let lastValid = -1;
    for (let r = 0; r < nRows; r++) {
      if (!rowXs[r].length) continue;
      let rx = concatF64(rowXs[r]), ry = concatF64(rowYs[r]);
      const order = Array.from(rx.keys()).sort((a, b) => rx[a] - rx[b]);
      const sx = new Float64Array(rx.length), sy = new Float64Array(rx.length);
      for (let i = 0; i < order.length; i++) { sx[i] = rx[order[i]]; sy[i] = ry[order[i]]; }
      rx = sx; ry = sy;
      let mx = -Infinity; for (let i = 0; i < rx.length; i++) if (rx[i] > mx) mx = rx[i];
      xEnds[r] = mx; baseYs[r] = median(ry);
      if (lastValid >= 0) cumT[r] = cumT[lastValid] + (xEnds[lastValid] - xLeft) * sPerPt;
      lastValid = r;
      const t = new Float64Array(rx.length), v = new Float64Array(rx.length), rr = new Int8Array(rx.length);
      for (let i = 0; i < rx.length; i++) {
        t[i] = cumT[r] + (rx[i] - xLeft) * sPerPt;
        v[i] = (baseYs[r] - ry[i]) * mvPerPt;
        rr[i] = r;
      }
      tParts.push(t); vParts.push(v); rParts.push(rr);
    }

    let t = concatF64(tParts), v = concatF64(vParts);
    // sort by time
    const ord = Array.from(t.keys()).sort((a, b) => t[a] - t[b]);
    const ts = new Float64Array(t.length), vs = new Float64Array(t.length);
    for (let i = 0; i < ord.length; i++) { ts[i] = t[ord[i]]; vs[i] = v[ord[i]]; }
    // unique by round(t,5), keep first occurrence (numpy.unique returns first index)
    const ut = [], uv = [];
    let prevKey = NaN;
    for (let i = 0; i < ts.length; i++) {
      const key = Math.round(ts[i] * 1e5) / 1e5;
      if (key !== prevKey) { ut.push(ts[i]); uv.push(vs[i]); prevKey = key; }
    }
    const utA = Float64Array.from(ut), uvA = Float64Array.from(uv);
    // resample to uniform grid: np.arange(t0, t_last, 1/FS)
    const t0 = utA[0], tLast = utA[utA.length - 1], step = 1 / FS;
    const nUni = Math.max(0, Math.ceil((tLast - t0) / step));
    const tUni = new Float64Array(nUni), vUni = new Float64Array(nUni);
    let j = 0;
    for (let i = 0; i < nUni; i++) {
      const tt = t0 + i * step;
      tUni[i] = tt;
      while (j < utA.length - 1 && utA[j + 1] < tt) j++;
      if (tt <= utA[0]) vUni[i] = uvA[0];
      else if (tt >= utA[utA.length - 1]) vUni[i] = uvA[utA.length - 1];
      else {
        const x0 = utA[j], x1 = utA[j + 1];
        const f = (tt - x0) / (x1 - x0);
        vUni[i] = uvA[j] + f * (uvA[j + 1] - uvA[j]);
      }
    }
    return { time: tUni, voltage: vUni, cum_t: cumT, baseline_y: baseYs, n_rows: nRows,
      x_left: xLeft, y_min: yMin, band, raw_r: concatI8(rParts),
      paper_speed: paperSpeed, volt_scale: voltScale, s_per_pt: sPerPt, mv_per_pt: mvPerPt };
  }

  // How many strips the page LAYOUT has room for, from the baselines that were
  // found: strips sit on a regular pitch, so span/pitch + 1 counts a strip the
  // extractor merged into a neighbour (the pitch is the median baseline gap,
  // which survives one missing strip). Only strips with a baseline are counted.
  function expectedRows(meta) {
    const ys = [];
    const seen = new Set(meta.raw_r);
    for (let r = 0; r < meta.n_rows; r++) if (seen.has(r)) ys.push(meta.baseline_y[r]);
    if (ys.length < 2) return ys.length;
    ys.sort((a, b) => a - b);
    const gaps = []; for (let i = 1; i < ys.length; i++) gaps.push(ys[i] - ys[i - 1]);
    const pitch = median(gaps);
    return pitch > 0 ? Math.round((ys[ys.length - 1] - ys[0]) / pitch) + 1 : ys.length;
  }

  function concatF64(list) {
    let n = 0; for (const a of list) n += a.length;
    const out = new Float64Array(n); let o = 0;
    for (const a of list) { out.set(a, o); o += a.length; }
    return out;
  }
  function concatI8(list) {
    let n = 0; for (const a of list) n += a.length;
    const out = new Int8Array(n); let o = 0;
    for (const a of list) { out.set(a, o); o += a.length; }
    return out;
  }

  // ── peak detection (mirrors detect_peaks) ─────────────────────────────── //
  function detectPeaks(voltage, sensitivity) {
    if (sensitivity == null) sensitivity = 8;
    const n = voltage.length;
    if (n < 20) return [];
    const filt = filtfilt(B, A, voltage);
    const deriv = new Float64Array(n); // np.diff(filt, prepend=filt[0]) → deriv[0]=0
    for (let i = 1; i < n; i++) deriv[i] = filt[i] - filt[i - 1];
    const sq = new Float64Array(n); for (let i = 0; i < n; i++) sq[i] = deriv[i] * deriv[i];
    const win = Math.floor(0.15 * FS);
    const mwi = uniformFilter1d(sq, win);
    const coeff = (10 - sensitivity) * 0.2 + 0.1;
    const thr = mean(mwi) + coeff * std(mwi);
    const peaks = findPeaks(mwi, { height: thr, distance: Math.floor(0.3 * FS) });
    const ref = Math.floor(0.05 * FS);
    const refined = new Set();
    for (const p of peaks) {
      const lo = Math.max(0, p - ref), hi = Math.min(n - 1, p + ref);
      refined.add(argmaxRange(voltage, lo, hi));
    }
    return Array.from(refined).sort((a, b) => a - b);
  }

  // ── beat classification (mirrors _classify_beats) ─────────────────────── //
  function classifyBeats(time, peakIdx) {
    const longSet = new Set(), shortSet = new Set();
    if (peakIdx.length < 3) return { long: longSet, short: shortSet };
    const rr = [];
    for (let i = 1; i < peakIdx.length; i++) rr.push((time[peakIdx[i]] - time[peakIdx[i - 1]]) * 1000);
    const med = median(rr);
    for (let i = 0; i < rr.length; i++) {
      if (rr[i] > med * (1 + ABNORMAL_THRESHOLD)) longSet.add(peakIdx[i + 1]);
      else if (rr[i] < med * (1 - ABNORMAL_THRESHOLD)) shortSet.add(peakIdx[i + 1]);
    }
    return { long: longSet, short: shortSet };
  }

  // ── pdf.js operator-list → stroked paths (page space, top-left origin) ─── //
  // Reproduces PyMuPDF get_drawings() coordinates + stroke widths. Key details
  // learned empirically from pdf.js 3.x output:
  //   • constructPath is emitted BEFORE setLineWidth within a drawing, so the
  //     width that matters is the one in effect at the *paint* (stroke) op.
  //   • several constructPath ops may accumulate into one path that a single
  //     paint op then strokes — they must be merged.
  // Coordinates are transformed on construction by (viewport ∘ CTM); the width
  // is scaled by the CTM at paint time — matching the on-page point/width units.
  function extractStrokePaths(fnArray, argsArray, viewportTransform, OPS) {
    const compose = (A, B) => [
      A[0] * B[0] + A[2] * B[1], A[1] * B[0] + A[3] * B[1],
      A[0] * B[2] + A[2] * B[3], A[1] * B[2] + A[3] * B[3],
      A[0] * B[4] + A[2] * B[5] + A[4], A[1] * B[4] + A[3] * B[5] + A[5],
    ];
    const scaleOf = (M) => Math.sqrt(Math.abs(M[0] * M[3] - M[1] * M[2]));
    const STROKE_OPS = new Set([OPS.stroke, OPS.closeStroke, OPS.fillStroke,
      OPS.eoFillStroke, OPS.closeFillStroke, OPS.closeEOFillStroke].filter((v) => v != null));
    let ctm = [1, 0, 0, 1, 0, 0];
    const stack = [];
    let lineWidth = 1;
    let bufXs = [], bufYs = []; // accumulated path since last paint op
    const out = [];
    const flush = (stroked) => {
      if (stroked && bufXs.length > 1) {
        out.push({ width: lineWidth * scaleOf(ctm), xs: Float64Array.from(bufXs), ys: Float64Array.from(bufYs) });
      }
      bufXs = []; bufYs = [];
    };
    for (let i = 0; i < fnArray.length; i++) {
      const fn = fnArray[i], args = argsArray[i];
      if (fn === OPS.save) stack.push(ctm.slice());
      else if (fn === OPS.restore) ctm = stack.pop() || [1, 0, 0, 1, 0, 0];
      else if (fn === OPS.transform) ctm = compose(ctm, args);
      else if (fn === OPS.setLineWidth) lineWidth = args[0];
      else if (fn === OPS.constructPath) {
        const M = compose(viewportTransform, ctm);
        const ops = args[0], co = args[1];
        let k = 0;
        const push = (x, y) => { bufXs.push(M[0] * x + M[2] * y + M[4]); bufYs.push(M[1] * x + M[3] * y + M[5]); };
        for (const op of ops) {
          if (op === OPS.moveTo || op === OPS.lineTo) { push(co[k], co[k + 1]); k += 2; }
          else if (op === OPS.curveTo) { push(co[k + 4], co[k + 5]); k += 6; }
          else if (op === OPS.curveTo2) { push(co[k + 2], co[k + 3]); k += 4; }
          else if (op === OPS.curveTo3) { push(co[k + 2], co[k + 3]); k += 4; }
          else if (op === OPS.rectangle) { k += 4; }
          // closePath consumes no coordinates
        }
      } else if (STROKE_OPS.has(fn)) flush(true);
      else if (fn === OPS.fill || fn === OPS.eoFill || fn === OPS.endPath) flush(false);
    }
    return out;
  }

  root.ECGCore = {
    PAPER_SPEED, VOLT_SCALE, PT_PER_MM, N_ROWS, FS, S_PER_PT, MV_PER_PT, ABNORMAL_THRESHOLD,
    // numeric (exported for tests)
    median, mean, std, uniformFilter1d, lfilter, filtfilt, findPeaks, peakProminences,
    localMaxima1d, histogram, searchsortedLeft,
    // pipeline
    parseScale, isEcgPath, countEcgVertices, findStripBoundaries, extractEcg, expectedRows,
    detectPeaks, classifyBeats, extractStrokePaths,
  };
})(typeof window !== "undefined" ? window : (typeof globalThis !== "undefined" ? globalThis : this));
if (typeof module !== "undefined") module.exports = (typeof window !== "undefined" ? window : globalThis).ECGCore;
