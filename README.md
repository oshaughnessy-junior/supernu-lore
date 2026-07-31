# SuperNu Lore

Working knowledge of **SuperNu** (LANL Monte-Carlo radiation transport,
<https://github.com/lanl/SuperNu>, GPL-3) as actually used for kilonova
radiative transfer: the input interfaces that matter, the opacity-table
machinery, the heating-source coupling, how to run it at scale on the OSPool,
and — mostly — the **silent failure modes**.

**Why this exists (and why it is NOT in a project repo):** this is
cross-cutting, slowly-changing lore about SuperNu's *interfaces and traps*.
It should survive project churn and transfer between projects that all drive
the same code (KN grids, adaptive inference, conditions inference, ...).
Project-specific pipelines stay in their own repos and link here.

Maintained by Richard O'Shaughnessy's group + Claude agents. Correct anything
you verify; mark unverified claims `(unverified)`. Cite `file:line` against a
named commit when you can.

Source citations below are against upstream commit
**`92816eb6a2b075fcb476a5199496284455b12635`** unless noted.

## GOLDEN RULES (read before running)

1. **SuperNu fails SILENTLY more often than it crashes.** Four separate
   input-handling mistakes we have hit all produce a run that completes,
   prints no warning, and is wrong. Always run an independent energy-closure
   check on any new handoff — see `energy-accounting.md`. This is the single
   most important thing on this page.
2. **`in_srctype` accepts `tabl` and `rpro`.** The inline comment on the
   declaration (`inputparmod.f:89`) lists only `none|heav|strt|manu|surf` and
   is **stale**; the case statement (~line 663) is authoritative.
3. **A wind-only `tabl` source needs an explicit `dyn_fr=0` column in
   `input.str`.** Without the column `gas_dynfr` defaults to **1**
   (`GAS/gas_setup.f:101`) and your `hrate_wnd.dat` source is multiplied by
   `(1-dynfr) = 0`. Zero heating, no warning.
4. **`Source/hrate_*.dat` has a RIGID header**: exactly 5 lines before the
   thermalization-option line. One extra comment line shifts the parse and
   again injects zero energy, silently. See `heating-sources.md`.
5. **Hand the tabular source RAW rates.** SuperNu applies Barnes–Kasen
   thermalization itself, per channel, from the file header. Pre-applying
   `f(t)` double-counts.
6. **The opacity table can be PRE-COARSENED, exactly.** `coarsen_tbxs` has no
   Planck/Rosseland weighting and no dependence on the local gas temperature,
   so the group-collapsed table is a reusable product. This turns a ~30 GB,
   multi-minute startup into ~0.3 GB and <1 s. Patch + recipe in
   `opacity-tables.md`. This is usually the difference between "runs on a few
   fat nodes" and "runs anywhere".
7. **Check element coverage per composition, never globally.** A table set
   that covers 99.9% of one composition can cover 8% of another (iron-peak vs
   lanthanide). Partial coverage biases brightness far more than colour.

## Contents

- `gotchas.md` — the traps, symptom -> cause. Start here.
- `opacity-tables.md` — Fontes table structure, what `coarsen_tbxs` actually
  computes, and the pre-coarsening patch/recipe.
- `heating-sources.md` — `tabl` / `rpro`, the 8-column format, thermalization,
  and coupling an external network (SkyNet) without double-counting.
- `energy-accounting.md` — how to prove SuperNu injected what you handed it.
- `inputs-reference.md` — `input.str` / `input.par` fields that actually bite.
- `building.md` — build recipes and runtime data files.
- `osg-deployment.md` — containers, HTCondor, measured resource footprints.

## Reusable artifacts

- `patches/supernu-precoarsen.patch` — adds `in_tbxs_dump` / `in_tbxs_coarse`
  (GPL-3, derivative of SuperNu; see `LICENSE-NOTES.md`).
- `tools/make_coarse_opacity.py` — builds the coarsened table using SuperNu
  itself, so generator and consumer can never disagree.
- `tools/verify_opacity_table.py` — acceptance test for a *raw* `opacities.h5`
  (h5py + numpy, analysis-machine): structure, grid conventions, C-order
  layout, per-composition coverage, and physics invariants that catch a
  silently transposed/mis-scaled table. Run it on any new table handoff.

**Not here, deliberately:** any LANL opacity data. The tables are not ours to
redistribute. This repo only tells you what to do with a set you already have.
