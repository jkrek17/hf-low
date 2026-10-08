#!/usr/bin/env node
/* Tests for docs/assets/js/era5.js and the committed docs/data/era5-tracks.js.
     node tests/era5/test_era5.js                                              */
var fs = require('fs'), path = require('path'), vm = require('vm'), assert = require('assert');
var ROOT = path.join(__dirname, '..', '..');
var ctx = { window: {} };
vm.createContext(ctx);
['assets/js/era5.js', 'data/era5-tracks.js', 'data/era5-null.js', 'data/era5-peaks.js', 'data/era5-prob.js', 'data/era5-lists.js', 'data/atlas.js'].forEach(function (f) {
  vm.runInContext(fs.readFileSync(path.join(ROOT, 'docs', f), 'utf8'), ctx, { filename: f });
});
var HF = ctx.window.HF, RAW = ctx.window.HF_ERA5;
var passed = 0, failed = 0;
function test(name, fn) {
  try { fn(); passed++; console.log('  ok   ' + name); }
  catch (e) { failed++; console.log('  FAIL ' + name + '\n       ' + (e && e.message)); }
}
var small = {
  meta: { threshold_kt: 71.7 },
  ev: [
    [1, 0, 2004, 2004120100, 800, 750, 9800, [200, -200, 1000, 60, 204, -196, 990, 75]],
    [2, 1, 2010, 2011010100, 500, 720, 9900, [180, 640, 1000, 50, 184, 644, 995, 72]],
    [3, 0, 1985, 1985020100, 700, 730, 9700, [210, -100, 1000, 60]]
  ]
};
var ev = HF.era5.decode(small);
test('decode: scaled fields, basin, month, and the HF flag from the gust', function () {
  assert.strictEqual(ev.length, 3);
  assert.strictEqual(ev[0].basin, 'atl'); assert.strictEqual(ev[1].basin, 'pac');
  assert.strictEqual(ev[0].fixes[0].lat, 50); assert.strictEqual(ev[0].fixes[0].lon, -50);
  assert.strictEqual(ev[1].fixes[0].lon, 160);
  assert.strictEqual(ev[0].pHf, 0.8); assert.strictEqual(ev[0].peakGust, 75); assert.strictEqual(ev[0].minP, 980);
  assert.strictEqual(ev[0].month, 12); assert.strictEqual(ev[1].month, 1);
  assert.strictEqual(ev[0].fixes[0].hf, false); assert.strictEqual(ev[0].fixes[1].hf, true);
  assert.strictEqual(ev[1].fixes[1].hf, true);                   // 72 >= 71.7
});
test('filter: basin, season range, months', function () {
  var f = function (o) { return HF.era5.filter(ev, o).map(function (e) { return e.key; }).join(); };
  assert.strictEqual(f({ basin: 'all', season0: 2004, season1: 2025 }), 'era5:1,era5:2');
  assert.strictEqual(f({ basin: 'all', season0: 1979, season1: 2025 }), 'era5:1,era5:2,era5:3');
  assert.strictEqual(f({ basin: 'pac', season0: 1979, season1: 2025 }), 'era5:2');
  assert.strictEqual(f({ basin: 'all', season0: 1979, season1: 2025, months: { 1: true } }), 'era5:2');
  assert.strictEqual(f({ basin: 'all', season0: 1979, season1: 2025, months: { 1: false } }), 'era5:1,era5:2,era5:3');
});
test('committed data: 4,157 events, 79,476 fixes, fixes in range, none lost in encoding', function () {
  assert.strictEqual(RAW.meta.events, 4157); assert.strictEqual(RAW.ev.length, 4157);
  var all = HF.era5.decode(RAW), n = 0;
  all.forEach(function (e) {
    n += e.fixes.length;
    e.fixes.forEach(function (p) { assert.ok(p.lat > 0 && p.lat < 90 && p.lon >= -180 && p.lon <= 180); });
    assert.ok(e.peakGust >= 71.7, e.key + ' peak gust below threshold');
    assert.ok(e.fixes.some(function (p) { return p.hf; }) , e.key + ' has no HF-equivalent fix');
  });
  assert.strictEqual(n, 79476); assert.strictEqual(RAW.meta.fixes, 79476);
  var pre = all.filter(function (e) { return e.season < 2001; }).length;
  assert.ok(pre > 0 && pre < all.length);
  assert.strictEqual(HF.era5.filter(all, { basin: 'all', season0: 1979, season1: 2025 }).length, 4157);
});
test('peaks: decode and filter; 75,087 cyclones, the HF-equivalent ones are exactly the 4,157 events', function () {
  var raw = ctx.window.HF_ERA5_PEAKS, P = HF.era5.decodePeaks(raw);
  assert.strictEqual(P.n, 75087); assert.strictEqual(raw.meta.events, 75087);
  var hf = 0; for (var i = 0; i < P.n; i++) hf += P.hf[i];
  assert.strictEqual(hf, 4157);
  var all = HF.era5.filterPeaks(P, { basin: 'all', season0: 0, season1: 9999 });
  assert.strictEqual(all.length, 75087);
  var atl = HF.era5.filterPeaks(P, { basin: 'atl', season0: 0, season1: 9999 }).length;
  var pac = HF.era5.filterPeaks(P, { basin: 'pac', season0: 0, season1: 9999 }).length;
  assert.strictEqual(atl + pac, 75087);
  var post = HF.era5.filterPeaks(P, { basin: 'all', season0: 2004, season1: 2025 }).length;
  assert.ok(post > 0 && post < 75087);
  var jan = HF.era5.filterPeaks(P, { basin: 'all', season0: 0, season1: 9999, months: { 1: true } }).length;
  assert.ok(jan > 0 && jan < 75087);
});
test('null cases: 4,154 tracks, none reaches the threshold, 6-hourly fixes', function () {
  var nul = HF.era5.decode(ctx.window.HF_ERA5_NULL);
  assert.strictEqual(nul.length, 4154);
  nul.forEach(function (e) { assert.ok(e.peakGust < 71.7, e.key + ' null case at or above threshold'); });
});
test('pClass: edges 5%, 20%, 50%, 80%; -1 is unscored', function () {
  var c = HF.era5.pClass;
  assert.deepStrictEqual([c(-1), c(0), c(0.049), c(0.05), c(0.199), c(0.2), c(0.5), c(0.79), c(0.8), c(1)].join(), '0,1,1,2,2,3,4,4,5,5');
});
test('attachProb: an edge takes the larger P of its two ends; unmatched tracks get none', function () {
  var t = HF.era5.decode(small);
  var n = HF.era5.attachProb(t, { p: { '1': [100, -1, 600], '2': [-1, 30] } });
  assert.strictEqual(n, 2);
  assert.strictEqual(t[0].pfix.join(), '0.1,-1,0.6');
  assert.strictEqual(Array.prototype.join.call(t[0].ec), '0,2,4');       // max(0.1, none) -> 5-20%; max(none, 0.6) -> 50-80%
  assert.strictEqual(Array.prototype.join.call(t[1].ec), '0,1');                 // 3% -> under 5%
  assert.strictEqual(t[2].ec, null);
});
test('P(HF) sidecar: every scored fix is a 00/12 UTC fix of a catalog track, P in [0, 1]', function () {
  var raw = ctx.window.HF_ERA5_PROB, ev = HF.era5.decode(ctx.window.HF_ERA5), nu = HF.era5.decode(ctx.window.HF_ERA5_NULL);
  var by = {}; ev.concat(nu).forEach(function (e) { by[e.key.slice(5)] = e; });
  var tracks = Object.keys(raw.p), scored = 0, pre = 0;
  assert.strictEqual(tracks.length, raw.meta.tracks);
  tracks.forEach(function (id) {
    var e = by[id], row = raw.p[id];
    assert.ok(e, 'track ' + id + ' is not a mapped track');
    assert.strictEqual(row.length, e.fixes.length, 'track ' + id + ' length');
    var hh = e.start % 100;
    row.forEach(function (v, i) {
      if (v < 0) return;
      assert.ok(v <= 1000, 'P out of range');
      assert.ok((hh + 6 * i) % 24 === 0 || (hh + 6 * i) % 24 === 12, 'track ' + id + ' fix ' + i + ' is not 00/12 UTC');
      scored++; if (e.season < 2004) pre++;
    });
  });
  assert.strictEqual(scored, raw.meta.scored_fixes);
  assert.ok(pre > 0 && pre < scored);
});
test('strongest-storm lists: 150 rows, every track on the map, gust ranks within an era only', function () {
  var L = ctx.window.HF_ERA5_LISTS, ev = HF.era5.decode(ctx.window.HF_ERA5), nu = HF.era5.decode(ctx.window.HF_ERA5_NULL);
  var ids = {}; ev.concat(nu).forEach(function (e) { ids[e.key.slice(5)] = 1; });
  assert.strictEqual(L.rows.length, 150);
  L.rows.forEach(function (r) { assert.ok(ids[r[3]], 'track ' + r[3] + ' not on the map'); });
  var g = L.rows.filter(function (r) { return r[0] === 'gust'; });
  assert.strictEqual(g.length, 50);
  ['1979-2000', '2004-2025'].forEach(function (era) {
    var e = g.filter(function (r) { return r[10] === era; }).sort(function (a, b) { return a[2] - b[2]; });
    assert.strictEqual(e.length, 25);
    e.forEach(function (r, i) {
      assert.strictEqual(r[2], i + 1);
      if (i) assert.ok(r[8] <= e[i - 1][8], 'gust not descending in ' + era);
      assert.ok(era === '1979-2000' ? r[4] < 2001 : r[4] >= 2004, 'season outside its era');
    });
  });
  ['minp', 'depth'].forEach(function (k) {
    assert.strictEqual(L.rows.filter(function (r) { return r[0] === k; }).length, 50);
  });
});
test('atlas boxes reproduce the atlas report: busiest 5 x 10 degree box per basin and source', function () {
  var A = ctx.window.HF_ATLAS;
  function busiest(src, b) {
    var best = null;
    A.box[src][b].cells.forEach(function (c) {
      var n = 0; for (var m = 2; m < 14; m++) n += c[m];
      if (!best || n > best.n) best = { n: n, i: c[0], j: c[1] };
    });
    var g = A.box[src][b];
    return { h: Math.round(best.n * 6 / A.meta.seasons * 10) / 10, lat0: g.lat0 + 5 * best.i, lon0: g.lon0 + 10 * best.j };
  }
  assert.deepStrictEqual(busiest('archive', 'atl'), { h: 96.8, lat0: 60, lon0: -40 });
  assert.deepStrictEqual(busiest('proxy', 'atl'), { h: 88.9, lat0: 60, lon0: -40 });
  assert.deepStrictEqual(busiest('archive', 'pac'), { h: 43.1, lat0: 40, lon0: 160 });
  assert.deepStrictEqual(busiest('proxy', 'pac'), { h: 47.7, lat0: 40, lon0: 160 });
  assert.strictEqual(A.meta.seasons, 22);
});
test('atlas motion boxes have at least 25 steps; historic candidates are pre-2004-05 and on mapped tracks', function () {
  var A = ctx.window.HF_ATLAS, ev = HF.era5.decode(ctx.window.HF_ERA5), ids = {};
  ev.forEach(function (e) { ids[e.key.slice(5)] = e; });
  ['archive', 'proxy'].forEach(function (src) { ['atl', 'pac'].forEach(function (b) {
    assert.ok(A.motion[src][b].length > 10);
    A.motion[src][b].forEach(function (r) { assert.ok(r[4] >= 25); });
  }); });
  assert.strictEqual(A.historic.length, 60);
  A.historic.forEach(function (r) {
    assert.ok(r[3] < 2004, 'season ' + r[3]);
    assert.ok(ids[r[2]], 'track ' + r[2] + ' is not an event track');
  });
});
console.log('\n' + passed + ' passed, ' + failed + ' failed');
process.exit(failed ? 1 : 0);
