# Heating sources: `tabl`, `rpro`, and coupling an external network

## `in_srctype` values that actually exist

The case statement (`inputparmod.f` ~line 663) accepts `tabl` and `rpro`.
The inline comment on the declaration (line 89) lists only
`none|heav|strt|manu|surf` and is **stale**. Ni56 machinery is bypassed under
`tabl`/`rpro`, so there is no double counting with the decay defaults.

## `tabl`: the 8-column table format

Read by `tbsrcmod.f` from `Source/hrate_dyn.dat` and/or `Source/hrate_wnd.dat`.
Reference examples ship with the code (`Data/Source/hrate_*.dat`, Korobkin
WinNet tables) — copy their layout exactly.

```
<line 1>  comment
<line 2>  comment
<line 3>  (blank)
<line 4>  comment
<line 5>  comment
<line 6>  thermalization OPTION for cols 4-8   (5 ints: 0=const eff, 1=Barnes-Kasen)
<line 7>  thermalization PARAMETER for cols 4-8 (5 floats)
<line 8>  (blank)
<line 9>  # <n_rows>
<line 10> comment
<line 11> comment
<line 12+> the table, 8 columns:
   1:time[s] 2:Qdot_total 3:Qdot_radiation 4:alpha 5:beta 6:gamma 7:electrons 8:fission
```

**The reader is positional.** It performs exactly five bare `read(4,*)` before
line 6. Add or remove one header line and it parses garbage — typically ending
with a zero source, no warning. Emit this file programmatically and assert on
the structure afterwards (re-read it and check line 6 has 5 fields, line 9
parses as `# n`, and the body rows have 8 columns).

### How the columns are used

- Columns **4, 5, 7, 8** (alpha, beta, electrons, fission) are thermalized into
  the material source, **each with its own** Barnes–Kasen parameter from the
  header: `f = ln(1 + 2p/(t rho)) / (2p/(t rho))`.
- Column **6** (gamma) becomes a grey gamma-transport source (`gas_decaygamma`,
  `kappa_gamma = in_opcapgam * rho`).
- Column 2 is not used as the source; it is bookkeeping.
- **Hand it RAW rates.** SuperNu applies `f(t)` itself. Pre-applying it
  double-counts. Consequence: lumping everything into one channel is not
  merely cosmetic — channels you leave empty never receive their own (often
  much more efficient) thermalization.
- Time lookup is **nearest-neighbour with edge clamping** (`binsrch`), not
  interpolation. The table must densely span the whole simulation window; a run
  outside it silently clamps to the edge rate.
- Rates are per gram, multiplied by `rho`, and weighted by `dynfr` /
  `(1-dynfr)` per cell. See the `dyn_fr` trap in `gotchas.md`.

## `rpro`: the built-in fit

`SOURCE/rprocess_fitfrm_source.f90` + `rprocmod.f90` implement a Rosswog &
Korobkin Ye-dependent analytic heating fit. Useful as an independent
cross-check on your own heating handoff.

When comparing against it, note its **hardcoded partition**:
`X_alpha=0.05, X_beta=0.20, X_gamma=0.35, X_ff=0.00` — summing to **0.60**.
So `rpro` delivers only 60% of its own `heating_rate`. If your table delivers a
different fraction, part of any discrepancy is bookkeeping, not physics.
Its `vexp` is derived per cell from `gas_vol/gas_mass` against a hardcoded
`M_ref = 0.05 Msun`, and it is only calibrated over its r-process `Ye` range —
extrapolating to high `Ye` (where heating is negligible anyway) diverges.

## Coupling an external network — the trap that cost us most

If you drive SuperNu from a reaction network (e.g. SkyNet), **check what the
network's exported heating rate actually contains.**

SkyNet's `NetworkOutput::HeatingRateVsTime` is

```
-dotEpsilonNu - sumDYMassExcess/dt
```

`sumDYMassExcess` is the **full** nuclear mass-excess release — for beta decay
that includes the energy the antineutrino carries away. `dotEpsilonNu` is
nonzero **only if a `NeutrinoReactionLibrary` is loaded**; a plain
REACLIB strong+weak+fission setup does not load one, so nothing is subtracted.

Feeding that straight to SuperNu as deliverable heating **over-heats by
~1/0.65 = 1.54x**. It is not obviously wrong in any single output.

Three ways to catch it:
1. Read the network source for what the rate includes.
2. Reconstruct it independently: sum `lambda_i Y_i(t) Q_i` over spontaneous
   decays from the same rate library. If it reproduces the exported rate to
   ~1%, nothing was subtracted (we measured 0.994–1.001 across 121 grid
   points); if neutrinos *had* been removed you would find ~1.54.
3. End-to-end: correct the rate and confirm the bolometric peak moves by the
   ratio you expect. (Ours moved 0.649 against 0.650 predicted — in a *grey*
   run. With real wavelength-dependent opacity the peak ratio is not the
   heating ratio, because the heating level also shifts escape timing and
   photospheric temperature; we saw 0.720 there.)

The fix is a per-decay channel decomposition of `Y(A,Z;t)` — which also gives
you the alpha/fission channels that matter at late times in very
neutron-rich ejecta (they grew from ~2% at 1 d to ~16% at 10 d at `Ye=0.05`)
and which a single lumped channel would silently mis-thermalize.
