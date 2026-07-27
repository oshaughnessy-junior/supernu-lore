# Patches

## `supernu-precoarsen.patch`

**Base commit:** upstream `lanl/SuperNu`
`92816eb6a2b075fcb476a5199496284455b12635`.
**License:** GPL-3.0 (derivative of SuperNu) — see `../LICENSE-NOTES.md`.

Adds two namelist logicals, **both default `.false.`**, so with them off the
read/coarsen path is byte-identical to upstream:

| flag | effect |
|---|---|
| `in_tbxs_dump` | after `coarsen_tbxs`, write `opacities_coarse.h5` (datasets `cap`, `sig`, `temp`, `rho`, `elem_codes`; attributes `ng`, `nelem`, `wlmin`, `wlmax`) |
| `in_tbxs_coarse` | read that file directly and subset to the elements in `input.str`, skipping the ~30 GB raw read |

Touches `inputparmod.f` (declare + namelist + `insertl`), `supernu.f90` (branch
around the read/coarsen block), `tbxsmod.f` (the two new subroutines).

Apply:

```bash
git -C <supernu-src> checkout 92816eb6a2b075fcb476a5199496284455b12635
git -C <supernu-src> apply /path/to/supernu-precoarsen.patch
```

Composes cleanly with an independent checkpointing overlay we also carry (they
touch disjoint regions of `supernu.f90`), but **verify any overlay stack
applies in sequence on a clean checkout** before trusting a build.

Rationale, sizing, and the validation procedure: `../opacity-tables.md`.

### Safety property worth preserving if you re-derive this

`read_tbxs_coarse` hard-stops when the stored `ng`/`wlmin`/`wlmax` disagree
with the run's group structure. Keep that. Silently mixing a mismatched grid
would be exactly the kind of quiet wrongness this code is otherwise prone to
(`../gotchas.md`).
