/* ERA5 proxy tracks: decode and filter only. No canvas, no DOM.

   window.HF_ERA5 (data/era5-tracks.js, built by tools/build_era5_tracks.py)
   holds pipeline A's hurricane-force-EQUIVALENT tracks: cyclones whose ERA5
   fields look like the ones OPC warned for as hurricane force (800 km gust
   index >= 71.7 kt). It is a proxy, not the archive. Counts before 2001-02
   are for comparison within that era only (STATUS.md, decision 1), so the
   page keeps them off unless asked.

   decode(raw) -> [{hf0, hf1 (first and last HF-equivalent fix index), key, basin, season, start, month, pHf, peakGust, minP,
                    fixes: [{lat, lon, pres, gust, hf}]}]
   filter(events, f) -> the subset matching the page's filters:
     f.basin   'all' | 'atl' | 'pac'
     f.season0, f.season1   inclusive season start years
     f.months  {1: true, ...} or empty/undefined for every month (the month
               of the track's first fix, as for the archive step view)
   decodePeaks(raw) / filterPeaks(P, f): the same for the all-cyclone peak
   positions (parallel typed arrays; filterPeaks returns indices).

   The class filter, the pressure cut and the text search describe archive
   fields the proxy does not have, so they are not applied. */

window.HF = window.HF || {};

(function (HF) {
  'use strict';

  function decode(raw) {
    var thr = raw.meta.threshold_kt;
    return raw.ev.map(function (r) {
      var flat = r[7], fixes = [];
      for (var i = 0; i < flat.length; i += 4) {
        var lat = flat[i] / 4, lon = flat[i + 1] / 4;
        // sp, cp, lam: sin/cos latitude and longitude in radians, so the
        // globe can project ~80,000 fixes a frame without redoing the trig.
        fixes.push({ lat: lat, lon: lon, pres: flat[i + 2], gust: flat[i + 3], hf: flat[i + 3] >= thr,
                     sp: Math.sin(lat * Math.PI / 180), cp: Math.cos(lat * Math.PI / 180), lam: lon * Math.PI / 180 });
      }
      var hf0 = -1, hf1 = -1;
      for (var k = 0; k < fixes.length; k++) if (fixes[k].hf) { if (hf0 < 0) hf0 = k; hf1 = k; }
      if (hf0 < 0) { hf0 = 0; hf1 = fixes.length - 1; }        // no HF fix (not in this catalog): treat the whole track as the stretch
      return {
        hf0: hf0, hf1: hf1,
        key: 'era5:' + r[0], basin: r[1] === 0 ? 'atl' : 'pac', season: r[2], start: r[3],
        month: Math.floor(r[3] / 10000) % 100, pHf: r[4] / 1000, peakGust: r[5] / 10, minP: r[6] / 10,
        fixes: fixes
      };
    });
  }

  function filter(events, f) {
    var months = f.months || {}, anyMonth = false, k;
    for (k in months) if (months[k]) { anyMonth = true; break; }
    return events.filter(function (e) {
      if (f.basin && f.basin !== 'all' && e.basin !== f.basin) return false;
      if (e.season < f.season0 || e.season > f.season1) return false;
      if (anyMonth && !months[e.month]) return false;
      return true;
    });
  }

  /** The all-cyclone population: one peak position per cyclone (the full
      tracks were never stored), as parallel arrays so 75,087 points stay
      cheap. thr is the pipeline A threshold from raw.meta. */
  function decodePeaks(raw) {
    var ev = raw.ev, n = ev.length, thr = raw.meta.threshold_kt;
    var P = { n: n, thr: thr, basin: new Uint8Array(n), season: new Uint16Array(n), month: new Uint8Array(n),
              hf: new Uint8Array(n), sp: new Float32Array(n), cp: new Float32Array(n), lam: new Float32Array(n),
              gust: new Float32Array(n), minP: new Float32Array(n) };
    for (var i = 0; i < n; i++) {
      var r = ev[i], lat = r[4] / 4, lon = r[5] / 4;
      P.basin[i] = r[1]; P.season[i] = r[2]; P.month[i] = Math.floor(r[3] / 10000) % 100;
      P.gust[i] = r[6] / 10; P.hf[i] = r[6] / 10 >= thr ? 1 : 0; P.minP[i] = r[7] / 10;
      P.sp[i] = Math.sin(lat * Math.PI / 180); P.cp[i] = Math.cos(lat * Math.PI / 180); P.lam[i] = lon * Math.PI / 180;
    }
    return P;
  }

  /** Indices of the peaks that pass the same filter as filter(). */
  function filterPeaks(P, f) {
    var months = f.months || {}, anyMonth = false, k, out = [];
    for (k in months) if (months[k]) { anyMonth = true; break; }
    var b = f.basin === 'atl' ? 0 : f.basin === 'pac' ? 1 : -1;
    for (var i = 0; i < P.n; i++) {
      if (b >= 0 && P.basin[i] !== b) continue;
      if (P.season[i] < f.season0 || P.season[i] > f.season1) continue;
      if (anyMonth && !months[P.month[i]]) continue;
      out.push(i);
    }
    return out;
  }

  /** P(HF within 24 h) classes for colouring: 0 = not scored, 1..5 by the
      upper edges below (probabilities, not percentages). */
  var P_EDGES = [0.05, 0.2, 0.5, 0.8];
  function pClass(p) {
    if (p < 0) return 0;
    var c = 1;
    while (c <= P_EDGES.length && p >= P_EDGES[c - 1]) c++;
    return c;
  }

  /** Attach data/era5-prob.js (raw = window.HF_ERA5_PROB) to decoded tracks:
      track.pfix (P per fix, -1 not scored) and track.ec, the class of each
      edge j-1 -> j at index j (0 for j = 0). Scoring is at 00/12 UTC, so
      every edge has a scored end; the edge takes the larger of its two ends'
      P, since P looks 24 h ahead. Tracks without scores keep ec = null.
      Returns how many tracks got scores. */
  function attachProb(tracks, raw) {
    var n = 0;
    tracks.forEach(function (t) {
      var row = raw.p[t.key.slice(5)];
      if (!row) { t.pfix = null; t.ec = null; return; }
      t.pfix = row.map(function (v) { return v < 0 ? -1 : v / 1000; });
      t.ec = new Uint8Array(row.length);
      for (var j = 1; j < row.length; j++) t.ec[j] = pClass(Math.max(t.pfix[j - 1], t.pfix[j]));
      n++;
    });
    return n;
  }

  /** Every track the P(HF) tables score (data/era5-allprob.js), as parallel
      arrays: one entry per track (id, basin, season, kind, month, first fix
      offset, count) and one per fix (trig for projection, P in 0..1000 and its
      class, time as 12-hour steps since 1979-01-01 00 UTC). See
      tools/build_era5_allprob.py for the wire format. */
  var EPOCH_MS = Date.UTC(1979, 0, 1);
  function decodeAllP(raw) {
    var T = raw.t, n = T.length, nf = 0, i, j;
    for (i = 0; i < n; i++) nf += T[i][5].length / 4;
    var D = { n: n, nfix: nf, id: new Int32Array(n), basin: new Uint8Array(n), season: new Uint16Array(n),
              kind: new Uint8Array(n), month: new Uint8Array(n), off: new Int32Array(n + 1),
              sp: new Float32Array(nf), cp: new Float32Array(nf), lam: new Float32Array(nf),
              lat: new Float32Array(nf), lon: new Float32Array(nf),
              p: new Uint16Array(nf), cls: new Uint8Array(nf), step: new Int32Array(nf) };
    var k = 0;
    for (i = 0; i < n; i++) {
      var r = T[i], f = r[5], lat4 = 0, lon4 = 0, st = r[4];
      D.id[i] = r[0]; D.basin[i] = r[1]; D.season[i] = r[2]; D.kind[i] = r[3];
      D.month[i] = new Date(EPOCH_MS + r[4] * 43200000).getUTCMonth() + 1;
      D.off[i] = k;
      for (j = 0; j < f.length; j += 4, k++) {
        st += f[j]; lat4 += f[j + 1]; lon4 += f[j + 2];
        var lat = lat4 / 4, lon = lon4 / 4;
        D.lat[k] = lat; D.lon[k] = lon;
        D.sp[k] = Math.sin(lat * Math.PI / 180); D.cp[k] = Math.cos(lat * Math.PI / 180); D.lam[k] = lon * Math.PI / 180;
        D.p[k] = f[j + 3]; D.cls[k] = pClass(f[j + 3] / 1000); D.step[k] = st;
      }
    }
    D.off[n] = k;
    return D;
  }

  /** Track indices that pass basin ('all' | 'atl' | 'pac'), season range and
      the month chips (month of the track's first fix), as an Int32Array. */
  function filterAllP(D, f) {
    var months = f.months || {}, anyMonth = false, k, out = [];
    for (k in months) if (months[k]) { anyMonth = true; break; }
    var b = f.basin === 'atl' ? 0 : f.basin === 'pac' ? 1 : -1;
    for (var i = 0; i < D.n; i++) {
      if (b >= 0 && D.basin[i] !== b) continue;
      if (D.season[i] < f.season0 || D.season[i] > f.season1) continue;
      if (anyMonth && !months[D.month[i]]) continue;
      out.push(i);
    }
    return Int32Array.from(out);
  }

  /** By month shading: 12-hourly fixes of every kind-0 low (2004-05 on) per
      DENS_LAT x DENS_LON cell for each calendar month, in fixes per season.
      b = 0 atl, 1 pac, -1 both; seasons s0..s1 inclusive. Returns
      {byMonth: [12 lists of {lat0, lon0, dlat, dlon, h}], max, nSeasons}. */
  var DENS_LAT = 5, DENS_LON = 10, DENS_NLAT = 18, DENS_NLON = 36;
  function monthDensity(D, b, s0, s1) {
    var fmon = D.fmon;
    if (!fmon) {
      fmon = D.fmon = new Uint8Array(D.nfix);
      for (var q = 0; q < D.nfix; q++) fmon[q] = new Date(EPOCH_MS + D.step[q] * 43200000).getUTCMonth();
    }
    var n = DENS_NLAT * DENS_NLON, grids = [], m, i, k;
    for (m = 0; m < 12; m++) grids.push(new Float32Array(n));
    for (i = 0; i < D.n; i++) {
      if (D.kind[i] !== 0 || (b >= 0 && D.basin[i] !== b) || D.season[i] < s0 || D.season[i] > s1) continue;
      for (k = D.off[i]; k < D.off[i + 1]; k++) {
        var li = Math.floor(D.lat[k] / DENS_LAT), lj = Math.floor((((D.lon[k] + 180) % 360 + 360) % 360) / DENS_LON);   // tracks that cross 180 carry on past it
        if (li < 0 || li >= DENS_NLAT || lj < 0 || lj >= DENS_NLON) continue;
        grids[fmon[k]][li * DENS_NLON + lj] += 1;
      }
    }
    var ns = s1 - s0 + 1, max = 0, byMonth = [];
    for (m = 0; m < 12; m++) {
      var cells = [];
      for (i = 0; i < n; i++) {
        if (!grids[m][i]) continue;
        var v = grids[m][i] / ns;
        if (v > max) max = v;
        cells.push({ lat0: Math.floor(i / DENS_NLON) * DENS_LAT, lon0: (i % DENS_NLON) * DENS_LON - 180,
                     dlat: DENS_LAT, dlon: DENS_LON, h: v });
      }
      byMonth.push(cells);
    }
    return { byMonth: byMonth, max: max, nSeasons: ns };
  }

  /** UTC date and hour of fix k, "24 Dec 2013 12 UTC". */
  var MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
  function fixTimeLabel(D, k) {
    var d = new Date(EPOCH_MS + D.step[k] * 43200000);
    return d.getUTCDate() + ' ' + MON[d.getUTCMonth()] + ' ' + d.getUTCFullYear() + ' ' + (d.getUTCHours() < 10 ? '0' : '') + d.getUTCHours() + ' UTC';
  }

  HF.era5 = { decodeAllP: decodeAllP, filterAllP: filterAllP, monthDensity: monthDensity, fixTimeLabel: fixTimeLabel, decode: decode, filter: filter, decodePeaks: decodePeaks, filterPeaks: filterPeaks,
              pClass: pClass, attachProb: attachProb, P_EDGES: P_EDGES };
})(window.HF);
