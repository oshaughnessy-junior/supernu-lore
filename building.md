# Building SuperNu

Upstream: <https://github.com/lanl/SuperNu> (GPL-3, Fortran + CMake, hybrid
MPI/OpenMP). Pin a commit; the repo moves.

## Minimal build

```bash
cmake -DUSE_MPI=ON -DUSE_HDF5=OFF \
      -DCMAKE_Fortran_FLAGS="-O2 -fallow-argument-mismatch" <srcdir>
make -j
```

- **`-fallow-argument-mismatch` is REQUIRED with gfortran >= 10.** `mpimod.f`
  calls the same MPI routine with different buffer types/ranks, which newer
  gfortran rejects as an error. Without it the build simply fails.
- Add `-ffree-line-length-none` if you patch free-form sources.

## HDF5 build (needed for `opacities.h5`)

```bash
cmake -DUSE_MPI=ON -DUSE_HDF5=ON \
      -DCMAKE_Fortran_FLAGS="-O2 -fallow-argument-mismatch" <srcdir>
```

Requires HDF5 **Fortran** bindings (`libhdf5_fortran`), not just the C library.
Both a system `libhdf5-*-fortran` and a conda-forge `hdf5` work.

## Runtime data

`Data/data.ion`, `Data/data.bf_verner`, `Data/data.ff_sutherland` ship in-tree
and must be present **in the run directory** (link them). `Data/Atoms/` ships
bound-bound line data for **H and C only** — everything else would have to be
fetched from Kurucz, whose lanthanide/actinide lists are severely incomplete.
That is why the inline-atomic-data path is a sensitivity bracket, not a science
mode, for r-process material.

## Smoke test that means something

Do **not** use `Testsuite/first` against its shipped reference: that reference
is stale upstream and "fails" on a correct build.

Use instead: the shipped `Data/Source/hrate_*.dat` with `in_srctype='tabl'`,
and check that energy deposits, gas temperatures rise, and
`output.flx_luminos` is nonzero.

## Editing fixed-form sources

`tbxsmod.f`, `inputparmod.f`, `tbsrcmod.f` etc. are **fixed-form Fortran**:
code must live in columns 7–72. Long `use ... , only: ...` lines silently
truncate at column 72 and produce baffling errors like
`Symbol 'hid' ... not found in module 'h5aux'` (from a truncated `HID_T`).
Keep patched lines short, or use a bare `use h5aux`.
