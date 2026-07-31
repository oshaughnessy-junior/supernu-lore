# Gotchas (symptom -> cause)

Ordered roughly by how much time each one costs before you find it. The first
four all produce a **completed run with no warning and a wrong answer**.

## Silent: the run finishes and is wrong

- **No `dyn_fr` column in `input.str` -> a wind-only `tabl` source is zeroed.**
  `gas_dynfr` defaults to **1** (`GAS/gas_setup.f:101`); `SOURCE/tabular_source.f90`
  weights `hrate_wnd.dat` by `(1-dynfr)` and `hrate_dyn.dat` by `dynfr`. Symptom:
  run completes, light curve is a bare cooling tail, `output.tot_energy` shows
  `eext` flat and `eerror` = 1. Fix: write an explicit `dyn_fr` column (0 for a
  pure wind, 1 for pure dynamical).

- **An extra comment line in `Source/hrate_*.dat` -> zero injected energy.**
  `read_tbsrc` (`tbsrcmod.f`) is positional, not token-driven: it does exactly
  five bare `read(4,*)` before the thermalization-option line, one more after
  the parameter line, then `read dmy, nt`, then two more, then the table.
  Symptom: identical to the `dyn_fr` failure. Fix: emit the header
  programmatically and assert on it — see `heating-sources.md`.

- **`tot_sthermal` is NOT the injected heating.** It is the *fictitious thermal
  re-emission* term, `sum(gas_emit)` before the material source is added
  (`SOURCE/sourceenergy.f90:49`). Comparing your input integral against it will
  mislead you by orders of magnitude. Use the per-step differences of
  cumulative `tot_eext` instead (`energy-accounting.md`).

- **Pre-applying thermalization `f(t)` to the source table double-counts.**
  SuperNu applies Barnes–Kasen itself, per channel, using the file's own
  header. Hand it raw rates.

## Loud, but the message is misleading

- **`unknown chemical element name:z104`** — SuperNu's element table
  (`elemdatamod.f`) stops at **Z=111**. Networks that reach further (SkyNet goes
  to Z=112) will emit labels it cannot parse. Only bites compositions with
  superheavy content, i.e. the low-`Ye` end; a blue-composition test will not
  reveal it. Fix: cap/lump at Z=111 before writing `input.str`.

- **`Array bound mismatch for dimension 1 of array 'tb_raw' (6/17)`** at
  `tbxsmod.f:209` — the raw `opacities.h5` is **transposed**. Each element
  dataset must be `(17 rho, 27 T, 14900 wl, 6 col)` in C/h5py order (component
  axis fastest); Fortran then sees `(6,14900,27,17)`. A writer that emits the
  axes in the other order produces this abort at load. Only bites tables you
  rebuilt yourself from the ASCII `op_*.table` files — run
  `tools/verify_opacity_table.py` after any such rebuild; it catches the
  transpose (and a *silent* mis-scaling that this abort would not) before you
  submit. See `opacity-tables.md`.

- **`STOP read_tbxs_coarse: ng /= in_grp_ng`** (with the pre-coarsen patch) —
  the group structure is baked into a coarsened table. `in_grp_ng`,
  `in_grp_wlmin`, `in_grp_wlmax` are then **not free parameters**; changing any
  of them requires regenerating the table. Working as designed: it refuses
  rather than silently mixing a mismatched grid.

- **HTCondor: "reading from file .../sim_out: No such file or directory"** —
  usually an *input* problem, not an output one. Your job script exited before
  creating the output dir. Read `job_N.err`, not the condor hold reason.

## Build / parse

- **gfortran >= 10 refuses to build without `-fallow-argument-mismatch`.**
  `mpimod.f` type-puns MPI buffer arguments. Not optional.

- **`Testsuite/first`'s shipped reference output is stale upstream.** It will
  "fail" against a correct build (the reference heats; the current code with
  that `input.par` has no active source). Do not use it as your build check —
  use a `tabl` smoke run with the shipped `Data/Source/hrate_*.dat` and confirm
  `output.flx_luminos` is nonzero.

- **Fortran `1p,e12.4` drops the `E` for 3-digit exponents.**
  `output.tot_energy` can contain literal `-2.5351-278`, which `numpy.loadtxt`
  rejects. Any parser needs `re.sub(r"(\d)([+-]\d\d\d)", r"\1E\2", line)`.
  Shows up once your source has genuinely tiny channels (alpha/fission).

- **Runtime data files must be in the run directory**: `data.ion`,
  `data.bf_verner`, `data.ff_sutherland`. Missing `data.ion` gives
  `STOP ion_read_data: cannot read data.ion` immediately.

## Scaling

- **The raw opacity read is the memory wall, and it is transient.**
  `tb_raw(ncol=7, ngr=14900, ntemp=27, nrho=17, nelem)` in float64: ~32 GB for
  85 elements, held only until `coarsen_tbxs` collapses it to
  `tb_cap(ng,27,17,nelem)` in float32 (~156 MB at ng=1024). Sizing a cluster
  request off the transient wastes an enormous amount of the pool. Pre-coarsen
  instead (`opacity-tables.md`).

- **Only ONE table load fits on a typical node.** ~30 GB each. Two concurrent
  loads on a 62 GB box will OOM. If you parallelise a slice locally, serialise
  the *load* phase and let transport (<1 GB) overlap freely.

- **Threaded runs are not reproducible.** With OpenMP the particle processing
  order — and hence RNG consumption — varies, so two identical configs give
  independent MC realisations. Do not expect bit-identity between runs; expect
  agreement within the MC noise floor, and measure that floor.
