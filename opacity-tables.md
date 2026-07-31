# Opacity tables: structure, what coarsening does, and pre-coarsening

## The Fontes per-element tables

SuperNu's tabulated-opacity mode (`in_notbopac=f`) uses **per-element** tables
(C. Fontes et al.), mixed per cell weighted by the actual mass fractions
(partial-density interpolation). Composition sensitivity is therefore
**continuous in `X_el`** — there is no `Ye`-binning anywhere in the opacity
step. (Language elsewhere about "N Ye-binned tables" refers to N *composition*
choices fed to this same machinery, not to the table structure.)

Two on-disk forms, same content:

- ASCII `Table/op_<el>_1Em<rho>gcc.table`, one file per (element, density).
- Consolidated HDF5 `opacities.h5` (`in_opfmthdf5=t`), one dataset per element
  plus `temperature_grid`, `density_grid`, `wavelength_grid`.

The three grids, and the two ways to misread them if you touch the h5 in Python
(SuperNu itself reads them correctly — these bite *your* tooling):

- `temperature_grid` — **eV**, ascending (a June-2025 set spans 0.01–5.0 eV, 27
  pts). `T[K] = T[eV] × 11604.5`.
- `density_grid` — **actual mass density in g cm⁻³**, ascending (e.g.
  `1e-20 … 1e-4`, 17 pts). It is **not** an exponent. A handoff note may
  describe it as "`10^-(value)`" — that is how it was *constructed*
  (`ρ_i = 10^(i-20)`), not an operation to apply. Running `10**-x` on the stored
  array gives `≈1.0` everywhere and **silently flattens the density axis**.
- `wavelength_grid` — dimensionless `u = E/kT` (see the column note below);
  physical energy is `wavelength_grid * temperature_grid[i]`.

Per element the array is **(17 rho, 27 T, 14900 wavelength, 6 col)** as h5py
sees it; Fortran reads it transposed. Column mapping (`tbxsmod.f`,
`tb_raw(2:7,:,:,:,l) = data4`), all columns in **cm² g⁻¹**:

| h5 col (0-based) | `tb_raw` row | meaning | used by SuperNu? |
|---|---|---|---|
| 0 | 2 | **total absorption** = col2+col3+col4 | no |
| 1 | 3 | **total extinction** = col0+col5 | no |
| 2 | 4 | bound-bound | yes (bb) |
| 3 | 5 | bound-free | yes (bf) |
| 4 | 6 | free-free | yes (ff) |
| 5 | 7 | scattering (**grey** — flat in wavelength) | yes (`tb_sig`) |

Cols 0,1 are **pre-summed totals, not "auxiliary noise"**: SuperNu ignores them
and rebuilds absorption from the bb/bf/ff components, but the two identities
hold pointwise to `< 1e-8` relative (we checked fe/nd/sm/u; residual is float32
round-off inherited from the ASCII source). That makes them a **free
internal-consistency check** on a table you did not build yourself — the
verifier below asserts on them. Do **not** add col0 and col1 together: col1
already contains col0.

`tb_raw(1,...)` is not read from file: it is built as
`tb_temp(itemp) * wavelength_grid`, i.e. the `wavelength_grid` is a
**non-dimensional energy** `u = E/T`, so the raw grid in physical energy
**shifts with the table temperature index**. (The dataset name is a legacy
misnomer — it is neither a wavelength nor a physical energy. Do not rename it;
`tbxsmod.f:199` looks it up by that exact string.) Confirmed independently from
the data: the free-free knee in `κ_ff·u³ ∝ (1−e^{−u})` sits at `u = 0.6925` at
every non-neutral temperature (spread `1.6e-16`) — fixed in grid units only if
the grid is already in units of `kT`.

## Verifying a raw table before you trust it

`tools/verify_opacity_table.py` (h5py + numpy; analysis-machine, not a worker)
runs the full acceptance set on a raw `opacities.h5`: structure, both grid
conventions, the C-order layout the reader demands, per-element population, the
component-sum identities above, and two physics invariants (the free-free knee;
scattering falling with density) that catch a **silently transposed or
mis-scaled** file which every structural check would still pass. Use `--require`
to check coverage for a specific composition, not just globally (golden rule 7):

```bash
python tools/verify_opacity_table.py opacities.h5 --require fe,nd,sm,u
```

Run it on any new table handoff, and always after rebuilding one from the ASCII
`op_*.table` files — the transpose trap (`gotchas.md`) aborts SuperNu, but
nothing warns you until it does.

**Constraint:** every element named in `input.str` must exist in the table set
or the reader hard-stops. Truncate/renormalise `X_el` to the tabulated list.

**Check coverage per composition.** A lanthanide/heavy-only set can cover
99.9% of a mid-`Ye` composition and **8%** of an iron-peak one. Partial
coverage biases **brightness** far more than colour: in one measured case
restoring the missing actinides (11–13% of mass, and the most opaque material
present) moved peak *g* by **2–3 mag** while the colour moved only 0.2–0.3 mag.
A colour-only sanity check will not catch it.

## What `coarsen_tbxs` actually computes

This is the key fact for scaling. Reading `tbxsmod.f`:

- For each (element, rho, temp) and each raw wavelength interval, it does a
  **trapezoid integral in inverse-energy space**, distributing the interval's
  bb+bf+ff opacity over whichever output groups it overlaps, then divides by
  the group's inverse-energy width.
- Scattering: `tb_sig = sum(tb_raw(7,:))/ngr`, a plain mean.
- **There is no Planck or Rosseland weighting, and no dependence on the local
  gas temperature.** It depends only on the table's own (T, rho) nodes and on
  the output group boundaries, which are fixed by
  `in_grp_ng` / `in_grp_wlmin` / `in_grp_wlmax`.

Therefore **the coarsened table is a well-defined, reusable product** for a
given group structure. Pre-computing it is *exact*, not an approximation.

`ngr = 14900` is a hardcoded private parameter, so you cannot simply supply a
smaller file — the pre-coarsen path below adds a separate reader instead.

## Why you want to pre-coarsen

`tb_raw` is allocated for **all elements at once**:

```
ncol(7) x ngr(14900) x ntemp(27) x nrho(17) x nelem x 8 bytes
```

= **~32 GB for 85 elements** — held only until coarsening, then freed. What
survives is `tb_cap(ng, ntemp, nrho, nelem)` in `real*4`: **156 MB** at
ng=1024, nelem=83.

Measured, same physics either way:

| | raw table | pre-coarsened |
|---|---|---|
| file staged | 27.9 GB | **0.32 GB** (87x) |
| peak RSS | ~30 GB | **269 MB** |
| startup | ~2.5–3 min | negligible |

On a shared pool that is the difference between matching ~0.2% of slots and
matching essentially all of them.

## The pre-coarsen patch

`patches/supernu-precoarsen.patch` adds two namelist logicals, both default
`.false.` (so with them off the read/coarsen path is byte-identical to
upstream):

- `in_tbxs_dump` — after coarsening, write `opacities_coarse.h5`
  (`cap`, `sig`, `temp`, `rho`, `elem_codes`, plus `ng`/`nelem`/`wlmin`/`wlmax`
  attributes).
- `in_tbxs_coarse` — read that file directly and subset it to the elements in
  `input.str`, skipping the raw read entirely.

**Design choice worth keeping:** the coarsened table is produced *by SuperNu
itself*, not by a reimplementation of the trapezoid logic in Python. Generator
and consumer are then the same code and cannot drift apart.

The reader hard-stops if the stored `ng`/`wlmin`/`wlmax` disagree with the run,
so a table built for another group structure cannot silently corrupt results.

### Recipe

1. Build SuperNu with `-DUSE_HDF5=ON` and the patch applied.
2. Run once with `in_tbxs_dump=t`, using an `input.str` listing **every**
   element in your table set, so the product is composition-independent.
   `tools/make_coarse_opacity.py` does this. Needs the ~32 GB transient once.
3. Ship the resulting `opacities_coarse.h5` and run with `in_tbxs_coarse=t`.

### Validating a coarsened table

Do **not** expect whole-run bit-identity — threaded MC diverges between runs
regardless. The sharp test is:

- **`output.grd_capgrey` at the FIRST timestep**, before any MC feedback,
  must be **bit-identical** to a full-table run (we measured `max|diff| = 0.0`
  across all cells). That isolates the table from the transport.
- Then light curves within the MC noise floor. Measure that floor rather than
  assuming it: at ng=1024 and 2^17 particles we saw ~0.01–0.03 mag in the
  optical, ~0.06 mag in K, 0.15% in bolometric — larger than a grey-mode
  figure would suggest, because photons are spread over 1024 groups.
