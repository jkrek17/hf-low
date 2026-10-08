#!/usr/bin/env node
/* Tests for docs/assets/js/era5.js and the committed docs/data/era5-tracks.js.
     node tests/era5/test_era5.js                                              */
var fs = require('fs'), path = require('path'), vm = require('vm'), assert = require('assert');
var ROOT = path.join(__dirname, '..', '..');
var ctx = { window: {} };
vm.createContext(ctx);
['assets/js/era5.js', 'data/era5-tracks.js'].forEach(function (f) {
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
console.log('\n' + passed + ' passed, ' + failed + ' failed');
process.exit(failed ? 1 : 0);
