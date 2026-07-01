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
  const PAPER_SPEED = 12.5;
  const VOLT_SCALE  = 2.5;
  const PT_PER_MM   = 2.834645669291339;
  const N_ROWS      = 14;
  const FS          = 2000;
  const S_PER_PT    = 1.0 / (PAPER_SPEED * PT_PER_MM);
  const MV_PER_PT   = 1.0 / (VOLT_SCALE  * PT_PER_MM);
  const ABNORMAL_THRESHOLD = 0.2;

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

  // ── strip-boundary detection (mirrors _find_strip_boundaries) ─────────── //
  function findStripBoundaries(pathYs) {
    let all = [];
    for (const ys of pathYs) all = all.concat(Array.from(ys));
    const arr = Float64Array.from(all);
    let yLo = Infinity, yHi = -Infinity;
    for (let i = 0; i < arr.length; i++) { if (arr[i] < yLo) yLo = arr[i]; if (arr[i] > yHi) yHi = arr[i]; }
    if (yHi <= yLo) return [];
    const nBins = Math.max(100, Math.floor(yHi - yLo));
    const { counts, edges } = histogram(arr, nBins);
    const win = Math.max(3, Math.floor(nBins / (N_ROWS * 4)));
    const density = uniformFilter1d(counts, win);
    let dmax = -Infinity; for (let i = 0; i < density.length; i++) if (density[i] > dmax) dmax = density[i];
    const neg = new Float64Array(density.length); for (let i = 0; i < density.length; i++) neg[i] = -density[i];
    const valleys = findPeaks(neg, { prominence: dmax * 0.25 });
    if (!valleys.length) return [];
    const binW = edges[1] - edges[0];
    const out = valleys.map((v) => edges[v] + binW / 2);
    out.sort((a, b) => a - b);
    return out;
  }

  // ── main extraction (mirrors extract_ecg) ─────────────────────────────── //
  // paths: [{width:Number, xs:Float64Array, ys:Float64Array}] — all stroked paths.
  function extractEcg(paths) {
    const empty = { time: new Float64Array(0), voltage: new Float64Array(0),
      cum_t: new Float64Array(N_ROWS), baseline_y: new Float64Array(N_ROWS),
      x_left: 0, y_min: 0, band: 1, raw_r: new Int8Array(0) };

    const pd = [];
    for (const p of paths) {
      if (!(p.width > 0.4 && p.width < 0.8)) continue;
      if (p.xs.length <= 1) continue;
      pd.push({ xs: p.xs, ys: p.ys, med: median(p.ys) });
    }
    if (!pd.length) return empty;

    pd.sort((a, b) => a.med - b.med);
    let xLeft = Infinity;
    for (const p of pd) for (let i = 0; i < p.xs.length; i++) if (p.xs[i] < xLeft) xLeft = p.xs[i];
    const meds = pd.map((p) => p.med);
    const yMin = meds[0], yMax = meds[meds.length - 1];
    const band = yMax > yMin ? (yMax - yMin) / Math.max(N_ROWS - 1, 1) : 1.0;

    const boundaries = findStripBoundaries(pd.map((p) => p.ys));
    const rows = meds.map((m) => {
      let r;
      if (boundaries.length) r = searchsortedLeft(boundaries, m);
      else r = Math.round((m - yMin) / band);
      return Math.min(Math.max(r, 0), N_ROWS - 1);
    });

    const rowXs = Array.from({ length: N_ROWS }, () => []);
    const rowYs = Array.from({ length: N_ROWS }, () => []);
    for (let i = 0; i < pd.length; i++) { rowXs[rows[i]].push(pd[i].xs); rowYs[rows[i]].push(pd[i].ys); }

    const xEnds = new Float64Array(N_ROWS), baseYs = new Float64Array(N_ROWS), cumT = new Float64Array(N_ROWS);
    const tParts = [], vParts = [], rParts = [];
    let lastValid = -1;
    for (let r = 0; r < N_ROWS; r++) {
      if (!rowXs[r].length) continue;
      let rx = concatF64(rowXs[r]), ry = concatF64(rowYs[r]);
      const order = Array.from(rx.keys()).sort((a, b) => rx[a] - rx[b]);
      const sx = new Float64Array(rx.length), sy = new Float64Array(rx.length);
      for (let i = 0; i < order.length; i++) { sx[i] = rx[order[i]]; sy[i] = ry[order[i]]; }
      rx = sx; ry = sy;
      let mx = -Infinity; for (let i = 0; i < rx.length; i++) if (rx[i] > mx) mx = rx[i];
      xEnds[r] = mx; baseYs[r] = median(ry);
      if (lastValid >= 0) cumT[r] = cumT[lastValid] + (xEnds[lastValid] - xLeft) * S_PER_PT;
      lastValid = r;
      const t = new Float64Array(rx.length), v = new Float64Array(rx.length), rr = new Int8Array(rx.length);
      for (let i = 0; i < rx.length; i++) {
        t[i] = cumT[r] + (rx[i] - xLeft) * S_PER_PT;
        v[i] = (baseYs[r] - ry[i]) * MV_PER_PT;
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
    return { time: tUni, voltage: vUni, cum_t: cumT, baseline_y: baseYs,
      x_left: xLeft, y_min: yMin, band, raw_r: concatI8(rParts) };
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
    findStripBoundaries, extractEcg, detectPeaks, classifyBeats, extractStrokePaths,
  };
})(typeof window !== "undefined" ? window : (typeof globalThis !== "undefined" ? globalThis : this));
if (typeof module !== "undefined") module.exports = (typeof window !== "undefined" ? window : globalThis).ECGCore;
