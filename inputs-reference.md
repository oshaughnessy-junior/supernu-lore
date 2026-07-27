# `input.str` / `input.par`: the fields that actually bite

## `input.str`

```
# <free-form comment>
# <nx> <ny> <nz> <ncol> <nabund>
# <col labels...>
<body: ncol numbers per cell, nx*ny*nz rows>
```

Labels are matched case-insensitively (`inputstrmod.f`). Recognised optional
per-cell columns: `temp`, `ye`, `cap`, `dyn_fr`; `mass` is required.
Abundance columns are labelled by **element symbol** and converted by
`elnam2elcode`; duplicates are rejected; `ni56`/`co56`/`fe52`/`mn52`/`cr48`/`v48`
are special-cased.

Traps:

- **`dyn_fr` is not optional in practice.** Absent -> `gas_dynfr = 1` -> a
  wind-only `tabl` source is silently zeroed. Write it explicitly.
- **Element symbols only go to Z=111.** Cap or lump anything higher.
- **1D spherical is `in_grd_igeom=11`**; for it the first column is the cell's
  right edge and the left edge of cell 1 is taken as 0.
- With `in_gas_gastempinit > 0` the `temp` column is **overridden** (you get a
  warning). Set `in_gas_gastempinit = 0d0` to honour `input.str`.
- `in_gas_cvcoef = -1d0` selects the physical heat capacity instead of the
  power-law form.

## `input.par` fragments worth memorising

Grey analytic opacity driven from the `cap` column (useful as a cheap
plumbing/contrast mode — it is achromatic, so it cannot produce line-blanketing
colour):

```
 in_opacanaltype = 'grey'
 in_gas_captpwr = 0d0
 in_gas_caprpwr = 1d0
 in_notbopac = t
 in_nobbopac = t / in_nobfopac = t / in_noffopac = t / in_nothmson = t
```

Per-element tabulated opacity (the science mode):

```
 in_notbopac = f
 in_notbbbopac = f / in_notbbfopac = f / in_notbffopac = f / in_notbthmson = f
 in_opfmthdf5 = t          ! consolidated opacities.h5 rather than ASCII Table/
 in_opacmixrossel = 1d0    ! Rosseland-weighted group mixing
 in_grp_ng = 1024          ! you need groups to have colour at all
```

External heating:

```
 in_novolsrc = f
 in_srctype = 'tabl'       ! or 'rpro'; 'none' + in_novolsrc=t for a null run
 in_opcapgam = 0.1d0
```

## Sanity runs worth having in your back pocket

- **null**: `in_srctype='none'`, `in_novolsrc=t`. Gives the floor from the
  thermal initial condition alone; anything heated must sit far above it.
- **grey vs tabulated at fixed everything else**: bolometric peak should agree
  in the optically-thin/blue limit (opacity redistributes, it does not create
  energy) but *not* at the red end, where wavelength-dependent opacity changes
  the escape timing. If they agree everywhere, suspect your opacity mode never
  engaged.
- **particle-count repeat**: same config at 2^N and 2^(N+2) — the difference is
  your MC noise floor, and you need it before believing any small difference.
