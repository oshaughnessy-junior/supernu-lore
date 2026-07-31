#!/usr/bin/env python3
"""Acceptance test for a SuperNu Fontes-format RAW HDF5 opacity table.

SuperNu fails silently more often than it crashes (README golden rule 1). A
mis-ordered, mis-scaled, or partially-populated opacity table is one of the
ways: several of the failures below either abort with an opaque Fortran message
or, worse, complete and are wrong. This runs the checks that catch them BEFORE
a campaign, against the raw ``opacities.h5`` (``in_opfmthdf5=t``) — structure,
the two grid conventions that are easy to misread, the C-order layout the
Fortran reader requires, per-element population, the component-sum identities,
and two physics invariants that catch a silently transposed/mis-scaled file
that every structural check would pass.

This is an ANALYSIS-MACHINE tool: it needs ``h5py``, which worker containers
deliberately do not ship (osg-deployment.md). Run it where you stage the table,
not inside a job.

Scope: the RAW per-element table. The pre-coarsened product
(``opacities_coarse.h5``, opacity-tables.md) has a different schema and its own
built-in ``ng``/``wlmin``/``wlmax`` guard; this tool does not check that one.

Usage
-----
    python verify_opacity_table.py /path/to/opacities.h5
    # --require checks coverage for a SPECIFIC composition (golden rule 7):
    # global coverage is not per-composition coverage.
    python verify_opacity_table.py /path/to/opacities.h5 --require fe,nd,sm,u

Exit status is 0 if every check passes, 1 otherwise. numpy 1.x / 2.x safe.
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np

# np.trapz was renamed np.trapezoid in NumPy 2.0; support both.
_trapz = getattr(np, "trapezoid", None) or np.trapz

# Compile-time constants in tbxsmod.f — the file must match these exactly.
N_RHO, N_TEMP, N_GROUP, N_COMP = 17, 27, 14900, 6

GRIDS = ("density_grid", "temperature_grid", "wavelength_grid")

# Trailing-axis component order (see opacity-tables.md, "The Fontes per-element
# tables"): 0=absorption total, 1=extinction total, 2=bb, 3=bf, 4=ff, 5=scatter.
ABS, EXT, BB, BF, FF, SCAT = range(6)

LANTHANIDES = "la ce pr nd pm sm eu gd tb dy ho er tm yb lu".split()
ACTINIDES = "ac th pa u np pu am cm bk cf es fm md no".split()
IRON_GROUP = "sc ti v cr mn fe co ni cu zn".split()

SIGMA_T = 6.6524587321e-25  # cm^2
N_A = 6.02214076e23


class Checker:
    def __init__(self) -> None:
        self.failed = 0

    def check(self, ok: bool, label: str, detail: str = "") -> bool:
        mark = "PASS" if ok else "FAIL"
        if not ok:
            self.failed += 1
        print(f"  [{mark}] {label}" + (f" — {detail}" if detail else ""))
        return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument(
        "--require",
        default="",
        help="comma-separated element symbols that must be present and populated",
    )
    args = ap.parse_args()

    c = Checker()
    f = h5py.File(args.path, "r")
    keys = list(f.keys())
    elements = sorted(k for k in keys if f[k].ndim == 4)

    print("\n== 1. Structure ==")
    for g in GRIDS:
        c.check(g in keys, f"grid dataset {g!r} present")
    c.check(
        len(elements) > 0, "element datasets present", f"{len(elements)} found"
    )
    c.check(
        not set(keys) - set(GRIDS) - set(elements),
        "no unexpected datasets",
        str(sorted(set(keys) - set(GRIDS) - set(elements))),
    )

    print("\n== 2. Grids ==")
    rho = f["density_grid"][:]
    T = f["temperature_grid"][:]
    u = f["wavelength_grid"][:]

    c.check(rho.shape == (N_RHO,), "density_grid length", f"{rho.shape} vs ({N_RHO},)")
    c.check(T.shape == (N_TEMP,), "temperature_grid length", f"{T.shape} vs ({N_TEMP},)")
    c.check(u.shape == (N_GROUP,), "wavelength_grid length", f"{u.shape} vs ({N_GROUP},)")

    c.check(np.all(np.diff(rho) > 0), "density_grid strictly ascending")
    c.check(np.all(np.diff(T) > 0), "temperature_grid strictly ascending")
    c.check(np.all(np.diff(u) > 0), "wavelength_grid strictly ascending")

    # density_grid holds ACTUAL densities in g/cm^3, not exponents.
    c.check(
        bool(rho.min() > 0 and rho.max() < 1.0),
        "density_grid values are densities in g/cm^3 (not exponents)",
        f"[{rho.min():.1e}, {rho.max():.1e}]",
    )
    c.check(
        bool(np.allclose(np.log10(rho), np.arange(N_RHO) - 20)),
        "density_grid follows rho_i = 10^(i-20) g/cm^3",
    )
    c.check(
        bool(T.min() >= 0.005 and T.max() <= 50.0),
        "temperature_grid plausible in eV",
        f"[{T.min()}, {T.max()}] eV",
    )

    print("\n== 3. Element dataset shape (C/h5py order) ==")
    want = (N_RHO, N_TEMP, N_GROUP, N_COMP)
    bad_shape = [e for e in elements if f[e].shape != want]
    c.check(
        not bad_shape,
        f"every element dataset is {want}",
        f"offenders: {bad_shape[:5]}" if bad_shape else "Fortran will see (6,14900,27,17)",
    )

    required = [s.strip().lower() for s in args.require.split(",") if s.strip()]
    if required:
        missing = [e for e in required if e not in elements]
        c.check(not missing, "required elements present", f"missing: {missing}")

    print("\n== 4. Group coverage ==")
    for group, name in (
        (LANTHANIDES, "lanthanides"),
        (ACTINIDES, "actinides"),
        (IRON_GROUP, "iron group"),
    ):
        have = [e for e in group if e in elements]
        c.check(
            len(have) == len(group),
            f"all {len(group)} {name} present",
            f"{len(have)}/{len(group)}; missing {sorted(set(group)-set(have))}"
            if len(have) != len(group)
            else "",
        )

    print("\n== 5. Per-element data sanity (sampled slab) ==")
    irho, itemp = N_RHO // 2, N_TEMP // 2
    empty, nonfinite, negative = [], [], []
    for e in elements:
        s = f[e][irho, itemp, :, :]
        if not np.all(np.isfinite(s)):
            nonfinite.append(e)
        if np.any(s < 0):
            negative.append(e)
        # A real table has a bound-bound line forest and a non-flat total.
        if s[:, BB].max() <= 0 or np.allclose(s[:, EXT], s[0, EXT]):
            empty.append(e)
    c.check(not nonfinite, "all sampled values finite", f"offenders: {nonfinite[:5]}")
    c.check(not negative, "no negative opacities", f"offenders: {negative[:5]}")
    c.check(
        not empty,
        "every element has real (non-flat, line-bearing) data",
        f"offenders: {empty[:10]}",
    )

    print("\n== 6. Component decomposition identities ==")
    tol = 1e-6
    for e in [x for x in ("fe", "nd", "sm", "u") if x in elements] or elements[:2]:
        s = f[e][irho, itemp, :, :]
        a, x = s[:, ABS], s[:, EXT]
        bb, bf, ff, sc = s[:, BB], s[:, BF], s[:, FF], s[:, SCAT]

        def relerr(p, q):
            d = np.abs(p - q)
            m = np.maximum(np.abs(p), np.abs(q))
            nz = m > 0
            return float(np.max(d[nz] / m[nz])) if np.any(nz) else 0.0

        c.check(relerr(a, bb + bf + ff) < tol, f"{e}: col0 == bb+bf+ff",
                f"rel {relerr(a, bb+bf+ff):.2e}")
        c.check(relerr(x, a + sc) < tol, f"{e}: col1 == col0+scatter",
                f"rel {relerr(x, a+sc):.2e}")
        c.check(bool(np.allclose(sc, sc[0])), f"{e}: scattering column is gray")

    print("\n== 7. Physics: energy grid is dimensionless u = h.nu/kT ==")
    # kappa_ff ~ nu^-3 (1 - exp(-u)).  y = kappa_ff * u^3 ~ (1 - e^-u), whose
    # half-plateau knee sits at u ~ 0.69 -- FIXED in grid units for every T if
    # (and only if) the grid is already in units of kT.
    probe = "fe" if "fe" in elements else elements[0]
    knees, skipped = [], []
    for it in range(0, N_TEMP, 2):
        y = f[probe][irho, it, :, FF] * u**3
        plateau = np.median(y[(u > 5) & (u < 11)])
        if plateau <= 0:
            # Cold enough that the gas is neutral: no free electrons, so the
            # free-free column is identically zero and the knee is undefined.
            skipped.append(float(T[it]))
            continue
        reg = u < 12
        knees.append(u[reg][np.argmin(np.abs(y[reg] - 0.5 * plateau))])
    knees = np.array(knees)
    spread = knees.std() / knees.mean() if len(knees) else np.inf
    c.check(
        len(knees) >= 5 and spread < 1e-3,
        "free-free knee is T-independent in grid units",
        f"knee={knees.mean():.4f} over {len(knees)} temperatures, spread={spread:.2e} "
        f"(=> E[eV] = wavelength_grid * T[eV])"
        + (f"; skipped neutral T={skipped} eV" if skipped else ""),
    )

    print("\n== 8. Physics: density axis orientation ==")
    # Mean ionisation inferred from the gray scattering column must DECREASE
    # with increasing density at fixed T (Saha recombination). If the rho axis
    # were reversed this trend inverts.
    it = int(np.argmin(np.abs(T - 1.0)))
    zbar = np.array([f[probe][i, it, 0, SCAT] for i in range(N_RHO)])
    c.check(
        bool(zbar[0] > zbar[-1] and np.mean(np.diff(zbar) <= 0) > 0.8),
        "scattering (∝ ionisation) falls with increasing density index",
        f"kappa_es {zbar[0]:.3e} -> {zbar[-1]:.3e} over rho {rho[0]:.0e}->{rho[-1]:.0e}",
    )

    print("\n== 9. Physics: lanthanide/actinide opacity contrast ==")
    # Rosseland mean of the total extinction. Because the grid is u = h.nu/kT,
    # the weighting function is temperature-independent.
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        uc = np.clip(u, 1e-30, 500.0)
        g = np.where(u < 500, (uc**4 * np.exp(-uc)) / np.expm1(-uc) ** 2, 0.0)
    g = np.abs(np.nan_to_num(g))

    def rosseland(el, i, it):
        k = f[el][i, it, :, EXT]
        return float(_trapz(g, u) / _trapz(g / np.maximum(k, 1e-99), u))

    i13 = int(np.argmin(np.abs(np.log10(rho) + 13)))
    ihalf = int(np.argmin(np.abs(T - 0.5)))
    heavy = [e for e in ("nd", "sm", "dy") if e in elements]
    light = [e for e in ("fe", "ni") if e in elements]
    if heavy and light:
        kh = np.mean([rosseland(e, i13, ihalf) for e in heavy])
        kl = np.mean([rosseland(e, i13, ihalf) for e in light])
        print(f"        rho=1e-13 g/cm3, T=0.5 eV: "
              f"kappa_R(lanthanide)={kh:.3g}, kappa_R(iron-group)={kl:.3g} cm^2/g")
        c.check(kh > 50 * kl, "lanthanide Rosseland opacity >> iron-group",
                f"ratio {kh/kl:.0f}x")

    print("\n== 10. Runtime footprint ==")
    per_elem = 7 * N_GROUP * N_TEMP * N_RHO * 8
    print(f"        tb_raw is {per_elem/2**30:.2f} GiB per element in the composition")
    print(f"        (SuperNu reads only elements named in input.str / input.compo)")

    f.close()
    print(f"\n{'ALL CHECKS PASSED' if not c.failed else f'{c.failed} CHECK(S) FAILED'}\n")
    return 1 if c.failed else 0


if __name__ == "__main__":
    sys.exit(main())
