# Energy accounting: proving SuperNu injected what you handed it

Because the input-handling failures are silent (`gotchas.md`), a closure test
is not optional hygiene — it is the only thing standing between you and a
plausible-looking wrong grid. Run it after **any** change to the handoff.

## Reading `output.tot_energy`

16 columns, header written by `OUTPUT/output_grid.f`:

```
eerror erad emat eext eout evelo sfluxgamma sflux
sthermal smanufac sanalvol sanalsurf samp sdecaygamma sdecaybeta sdeposgamma
```

Two things to know before using it:

- **`eext` is CUMULATIVE**; most `s*` columns are per-step tallies.
- **`sthermal` is NOT your injected heating.** It is the fictitious thermal
  re-emission term (`sum(gas_emit)` before the material source is added,
  `SOURCE/sourceenergy.f90:49`). Comparing against it will mislead you badly.

Per step, `tot_eext` gains
`dt*sum(vol*matsrc)` (`sourceenergy.f90:58`) + gamma deposition
(`sourceenergy_misc.f:47`) + particle amplification, plus the initial energy on
the first step. So the injected material source is

```
E_matsrc = sum_{it>=2} [ eext(it) - eext(it-1) - sdeposgamma(it) - samp(it) ]
```

and the emitted gamma energy is `sum_{it>=2} sdecaygamma(it)`. Skip step 1: its
`eext` bundles the initial energy.

Parsing note: use an exponent-repair regex — Fortran drops the `E` for 3-digit
exponents (`gotchas.md`).

## The input side

Re-implement what `tabular_source.f90` does, from your own table:

```
E_therm = sum_it dt * sum_cells [ m_cell * f_BK(param_ch; t, rho_cell) ] * q_ch(t)
E_gamma = sum_it dt * sum_cells   m_cell                                 * q_gamma(t)
```

with **nearest-neighbour** lookup into your table (matching `binsrch`), the
per-channel Barnes–Kasen parameter from your own header, and homologous
densities `rho_i(t) = m_i / (vol0_i t^3)`.

## What "passing" looks like

We hold the handoff to **10%**; a correct one lands near **2%**, with the
residual explained by nearest-neighbour lookup versus midpoint integration on
an exponential time grid. SuperNu's own `eerror` should be ~1e-5 or smaller.

Also worth recording per run: the gamma **deposition fraction**
(`sdeposgamma/sdecaygamma`, ~0.8 in our optically-thick KN configs) and the
per-channel breakdown of `E_therm`, which is what tells you whether
alpha/fission are being thermalized separately or silently lumped.

## Failure signatures

| symptom | cause |
|---|---|
| `rel_diff = -1.0`, `eerror = 1.0`, `eext` flat | zero source: missing `dyn_fr` column, or a shifted `hrate` header |
| input side ~1.54x the SuperNu side | network rate still contains neutrino energy (`heating-sources.md`) |
| ~2x, only vs `rpro` | partition bookkeeping: `rpro` delivers 0.60 of its own rate |
| a few % | expected: nearest-neighbour vs midpoint integration |
