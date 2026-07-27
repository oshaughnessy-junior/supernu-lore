# Licensing / redistribution

**Prose (`*.md`)** — group working notes. Reuse freely with attribution.

**`patches/supernu-precoarsen.patch`** — a derivative work of
[SuperNu](https://github.com/lanl/SuperNu), which is **GPL-3.0**. The patch is
therefore distributed under **GPL-3.0**. It is a diff, not a redistribution of
SuperNu itself; apply it to your own checkout of the upstream commit named in
the patch header. SuperNu is
"© Triad National Security, LLC ... released under the terms of the GNU GPLv3";
see the upstream `COPYING`.

**`tools/make_coarse_opacity.py`** — ours; same terms as the prose. It drives
SuperNu but contains none of its code.

**Opacity data is NOT here and must not be added.** The Fontes/LANL tables
(ASCII `Table/` or consolidated `opacities.h5`) carry unresolved redistribution
terms. Treat them as internal-use-only: this repo describes what to *do* with a
set you already have access to, and deliberately contains none of it. The same
applies to any derived coarsened table — it is a lossy but direct transform of
the licensed data, so it inherits the restriction.

**Cite the codes you use.** SuperNu: Wollaeger & van Rossum. Coupled networks
have their own citation conditions (e.g. SkyNet asks that you cite
arXiv:1706.06198 as a condition of use in publications).
