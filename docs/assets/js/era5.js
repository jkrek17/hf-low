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

  HF.era5 = { decode: decode, filter: filter };
})(window.HF);
