#!/usr/bin/env python3
"""Build a PRE-COARSENED SuperNu opacity table.

Why: SuperNu holds tb_raw(7, 14900, ntemp, nrho, nelem) in float64 while
coarsening -- ~32 GB for 85 elements -- then immediately collapses it to
tb_cap(ng, ntemp, nrho, nelem) in float32 (~156 MB at ng=1024). Shipping the
COARSENED table instead makes memory and disk trivial on a batch system.

Why it is exact, not an approximation: coarsen_tbxs (tbxsmod.f) is a pure
trapezoid average in inverse-energy space evaluated at the table's OWN
(T, rho) nodes. No Planck/Rosseland weighting, no dependence on the local gas
temperature -- only on the group boundaries, fixed by in_grp_ng / wlmin /
wlmax. So the coarsened table is a well-defined, reusable product.

How: rather than reimplementing that trapezoid logic (and risking a silent
divergence), we let SuperNu generate what SuperNu consumes. The patch in
patches/supernu-precoarsen.patch adds
    in_tbxs_dump    -> write opacities_coarse.h5 after coarsening
    in_tbxs_coarse  -> read it back, skipping the raw read
The dump uses an input.str listing EVERY element in the table, so the product
is composition-independent; the reader subsets to what a given run needs.

The file records ng/wlmin/wlmax and the reader hard-stops on a mismatch, so a
table built for a different group structure cannot silently corrupt a run.

Requires: the patched, HDF5-enabled SuperNu; ~32 GB RAM once; h5py.

Usage:
    make_coarse_opacity.py --supernu ./supernu --opacities opacities.h5 \
                           --data-dir /path/to/SuperNu/Data [--ng 1024]
"""

import argparse
import os
import shutil
import subprocess
import sys

# Z -> symbol, 1..111 (SuperNu's elemdatamod.f stops at 111)
ELEMENT_SYMBOLS = (
    "h he li be b c n o f ne na mg al si p s cl ar k ca sc ti v cr mn fe co "
    "ni cu zn ga ge as se br kr rb sr y zr nb mo tc ru rh pd ag cd in sn sb "
    "te i xe cs ba la ce pr nd pm sm eu gd tb dy ho er tm yb lu hf ta w re "
    "os ir pt au hg tl pb bi po at rn fr ra ac th pa u np pu am cm bk cf es "
    "fm md no lr rf db sg bh hs mt ds rg").split()

PAR = """&inputpars
 in_comment = "coarse-opacity-dump"
 in_nomp = {nomp}
 in_noreadstruct = f
 in_noeos = t
 in_io_nogridgroupdump = t
 in_tsp_tfirst = 1d4
 in_tsp_tlast = 1.2d4
 in_tsp_gridtype = 'expo'
 in_tsp_nt = 2
 in_grd_igeom = 11
 in_ndim = {nx}, 1, 1
 in_isvelocity = t
 in_gas_gastempinit = 0d0
 in_gas_cvcoef = -1d0
 in_flx_ndim = {ng}, 1, 1
 in_grp_ng = {ng}
 in_novolsrc = t
 in_srctype = 'none'
 in_notbopac = f
 in_notbbbopac = f
 in_notbbfopac = f
 in_notbffopac = f
 in_notbthmson = f
 in_opfmthdf5 = t
 in_tbxs_dump = t
 in_nobbopac = t
 in_nobfopac = t
 in_noffopac = t
 in_nothmson = t
 in_noplanckweighting = t
 in_opacmixrossel = 1d0
 in_src_n2s = 10
 in_prt_n2max = 18
/
"""

MSUN, CLIGHT = 1.9891e33, 2.99792458e10


def write_str(path, els, nx=8, mej=0.035, vej=0.15):
    """Minimal 1D spherical structure listing every element in the table.

    Composition is irrelevant to the coarsening (it is per-element); all that
    matters is that every element is PRESENT so tb_nelem covers the table.
    """
    # gas temperature is irrelevant here: coarsen_tbxs evaluates at the
    # TABLE's own (T, rho) nodes, not at the cell state. A flat 1e4 K is fine.
    x = 1.0 / len(els)
    vr = [(i + 1) / nx * vej * CLIGHT for i in range(nx)]
    m = mej * MSUN / nx
    ncol = 6 + len(els)
    with open(path, "w") as f:
        f.write("# all-element structure for a coarse-opacity dump\n")
        f.write(f"# {nx} 1 1 {ncol} {len(els)}\n")
        f.write("# " + "  ".join(f"{c:>10s}" for c in
                ["rightvel", "mass", "temp", "ye", "cap", "dyn_fr"] + els) + "\n")
        for i in range(nx):
            row = [vr[i], m, 1.0e4, 0.25, 1.0, 0.0] + [x] * len(els)
            f.write("  ".join(f"{v: .6e}" for v in row) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--supernu", required=True, help="patched, HDF5-enabled binary")
    p.add_argument("--opacities", required=True, help="raw opacities.h5")
    p.add_argument("--data-dir", required=True, help="SuperNu Data/ directory")
    p.add_argument("--ng", type=int, default=1024)
    p.add_argument("--nomp", type=int, default=8)
    p.add_argument("--workdir", default="./coarse_gen")
    args = p.parse_args()

    import h5py
    with h5py.File(args.opacities, "r") as f:
        els = sorted(k for k in f.keys() if k.isalpha() and k not in
                     ("density_grid", "temperature_grid", "wavelength_grid"))
    unknown = [e for e in els if e not in ELEMENT_SYMBOLS]
    if unknown:
        sys.exit(f"table contains elements SuperNu cannot label: {unknown}")
    print(f"{len(els)} elements in {args.opacities}")

    d = args.workdir
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    write_str(os.path.join(d, "input.str"), els)
    open(os.path.join(d, "input.par"), "w").write(
        PAR.format(ng=args.ng, nomp=args.nomp, nx=8))
    os.symlink(os.path.abspath(args.opacities), os.path.join(d, "opacities.h5"))
    for f in ("data.ion", "data.bf_verner", "data.ff_sutherland"):
        os.symlink(os.path.join(os.path.abspath(args.data_dir), f),
                   os.path.join(d, f))

    print(f"running SuperNu in {d} (needs ~32 GB RAM once) ...")
    env = dict(os.environ, OMP_NUM_THREADS=str(args.nomp))
    with open(os.path.join(d, "gen.log"), "w") as log:
        rc = subprocess.call([os.path.abspath(args.supernu)], cwd=d,
                             stdout=log, stderr=log, env=env)
    out = os.path.join(d, "opacities_coarse.h5")
    if not os.path.exists(out):
        sys.exit(f"FAILED (rc={rc}); see {d}/gen.log")
    print(f"wrote {out} ({os.path.getsize(out)/1e9:.3f} GB, "
          f"raw was {os.path.getsize(args.opacities)/1e9:.1f} GB)")
    print("validate: a run using it must reproduce output.grd_capgrey at the "
          "FIRST timestep bit-identically vs the raw table.")


if __name__ == "__main__":
    main()
