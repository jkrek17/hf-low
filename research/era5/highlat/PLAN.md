# High-latitude Atlantic HF lows: hypotheses and planned tests

Written 2026-10-08 before the main analysis was run, so the commit time shows
the plan came first. Results go in `results.txt` next to this file; null
results are reported in full.

## Question

In the archive, 37.5% of Atlantic HF lows (2004-05 to 2025-26) have an HF fix
at or north of 60N, against 1.4% in the Pacific. Pipeline A's ERA5 proxy gives
37.5% of Atlantic events peaking at or north of 60N (1979-2025; 34.0% from
2004-05). Why, and can the proxy be trusted there?

## What had been looked at before this note

- The catalog shares above, and the count of HF-strength times to size the pull.
- `gustloc.py` was run on two times (one track, 2026-04-30, centre 63.5N 38-39W)
  to check that it reproduces pipeline A's index exactly. It did. In both, the
  gust maximum sat 450-530 km south-southwest of the centre near Cape Farewell,
  in westerly flow, and masking ocean within 300 km of Greenland took the index
  below threshold. That is one storm and is not evidence about the population.
- Nothing else: no regional statistics, no archive comparison.

## Hypotheses

- **H1 (terrain jets).** At Atlantic HF-strength fixes north of 60N, the pipeline A
  gust maximum lies far from the centre (> 400 km) and close to Greenland
  (<= 300 km) much more often than at fixes south of 60N.
- **H0 for H1.** The maximum sits at similar distances from the centre in both
  regions, and the high-latitude share comes from deep storms whose own wind
  field carries the index.
- **H2 (orographic signature in the archive).** Archive HF fixes near Greenland
  carry weaker (higher) central pressure than HF fixes south of 60N. A storm
  needs less depth to produce HF wind where terrain helps.
- **H3 (proxy artefact).** Pipeline A over-predicts relative to the archive north
  of 60N: higher false-alarm ratio there than south of 60N.

## Planned tests (2004-05 onward; pipeline A; ERA5 is a proxy)

1. Location of the gust maximum at every in-domain HF-strength fix (index >= 71.7 kt)
   of every Atlantic pipeline A event: distance and bearing from the centre,
   distance to Greenland, Iceland and any land, sea-ice concentration, and 10 m
   wind direction there. Compare north and south of 60N. Unit for inference:
   season (block bootstrap by season, 22 seasons).
2. Terrain-mask sensitivity of the track index: the share of events that fall
   below 71.7 kt when ocean within 100 km of Greenland, within 300 km of
   Greenland, within 50 km of any land, or under sea ice > 0.15 is removed, and
   when only points within 400 km of the centre count. By peak latitude.
3. Archive pressure at HF fixes by region (south of 60N; Cape Farewell and
   Irminger Sea; Denmark Strait and Iceland; Norwegian Sea), one value per event
   (median over its HF fixes in the region), season-bootstrap CI on differences.
4. Pipeline A against the archive by region: POD and FAR, matching as the
   calibration does (a fix within 400 km at the same time; 800 km for pressure-less
   events), split by whether the event has a fix at or north of 60N. Note that
   pipeline A's Atlantic domain stops at 67N, so archive HF fixes north of 67N
   cannot be matched there.
5. In situ: Prins Christian Sund (ISD 043900, the Cape Farewell station Moore
   and Renfrew used), wind and gust at times when the ERA5 maximum is in the
   Cape Farewell jet region against other HF-strength times.

## How the answer will be read

- Terrain jets are real wind, and OPC warns for them (the archive has a
  tip-jet class). A large terrain share is therefore an explanation, not by
  itself an artefact.
- An artefact is (a) an index carried by grid points within 100 km of
  Greenland, where 0.25-degree cells mix land, ice and steep terrain, or by
  ice-covered points, or (b) a higher false-alarm ratio north of 60N.
- Recommend a terrain or ice mask for pipeline A if more than 10% of
  Atlantic events (2004-05 on) depend on points within 100 km of Greenland or
  under sea ice, or if FAR north of 60N exceeds FAR south of 60N by more than
  0.10 with a season-bootstrap 95% interval that excludes zero.
