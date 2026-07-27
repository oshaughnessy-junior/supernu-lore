# Running SuperNu at scale (OSPool / HTCondor)

Numbers below are measured for 1D spherical KN runs (nx=32, ng=1024, nt=96);
treat the shapes as transferable and re-measure the magnitudes.

## Container

Build the binary into an image and pin the upstream commit inside the
definition, so the container records exactly what it contains. Two things
worth doing:

- Write provenance files into the image (`/opt/SUPERNU_COMMIT.txt`, one per
  applied overlay). Then a job's MANIFEST can record which code produced it.
- If you apply multiple patch overlays, verify they apply **in sequence** on a
  clean checkout before trusting the build. Ours (checkpointing + pre-coarsen)
  touch disjoint regions of `supernu.f90` and compose cleanly, but that is
  luck, not design.

Run with `universe = container` and `container_image = <osdf url>`; the job
then executes *inside* the image, so point your wrapper at the in-container
binary path rather than nesting apptainer.

## Data staging

OSDF caches **immutably by file name**. Never overwrite a staged object —
bump a version suffix (`-v1` -> `-v2`) and update the submit file. Verify by
checksum on the far end after transfer; it is cheap and it has caught a real
truncation for us before.

Stage the **pre-coarsened** opacity table (`opacity-tables.md`), not the raw
one. That single change moved us from a ~40 GB memory request (matching ~0.2%
of the pool) to **2 GB** (matching essentially all of it).

## Measured footprint

| quantity | value |
|---|---|
| MemoryUsage | **269 MB** (pre-coarsened; dominated by `tb_cap`, 156 MB) |
| DiskUsage | 750 MB |
| wall (2^12 particles) | 190–248 s |
| output per sim | **3.7 MB** |

**Output size scales with `nt x ng x ncells`, not with particle count** —
raising particles for MC convergence does not inflate storage. Budgeting a
1,536-sim campaign at a per-sim figure borrowed from a higher-resolution
project overestimated our storage by ~60x (6 GB actual vs 384 GB assumed).
Measure one run before sizing scratch.

## The submit-file trap that does not fail

```
transfer_output_files  = sim_out
```

Outputs return to the **submit directory**, so every process in a cluster
writes to the same `./sim_out` and they **clobber one another**. This does not
error: an N-job campaign "succeeds" and returns one directory, silently
discarding N-1 results. Always remap:

```
transfer_output_remaps = "sim_out = results/$(BATCH)/sims/sim_$(Process)"
```

## Job-wrapper conventions that pay off

- **Write the MANIFEST last.** Anything that sweeps completed sandboxes (a
  streaming cron, a collector) can then treat its presence as "this directory
  is complete and safe to move".
- **Reserve a distinct exit code for "checkpointed, migrate me"** (we use 85,
  matching `checkpoint_exit_code`) so a genuine failure is never mistaken for
  an eviction.
- **Hard-fail on missing inputs rather than letting a writer fall back to a
  default.** Our wrapper exits non-zero if the per-channel heating file is
  absent, because the fallback silently produced over-heated runs.
- **Reduce in-job.** Post-process to photometry inside the job; ship kB, not
  the full spectral arrays, unless you need them.
- Do **not** assume a worker container has your analysis stack. Ours ships
  `python3-numpy` only — an `import h5py` in the input writer failed every job.
  Prefer plain-text sidecars or hardcoded constants for anything a worker
  needs, and keep the h5-reading paths on the analysis machine.

## Smoke before scale

A 6-job smoke batch found four distinct bugs that no amount of local testing
would have: a tarball packed without a top-level directory, the missing h5py,
the output clobbering above, and a smoke setting that changed `in_grp_ng`
against a table built for a different one. Budget a smoke round; it is
cheaper than a wrong campaign.
