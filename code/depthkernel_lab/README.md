# Depth Kernel Lab

This folder isolates the current DSurfTomo depth-kernel calculation from the full inversion.

The baseline driver reproduces the existing finite-difference logic for one 1D column:

```text
1 baseline surfdisp96 run + 6 * nz perturbed surfdisp96 runs
```

It computes Rayleigh phase velocity by default and writes the local kernels:

- `dC/dVs`
- `dC/dVp`
- `dC/drho`

## Build And Run

```bash
make
make run
```

Output is written to `baseline_kernel.out`.

To compare the same node-perturbation finite-difference kernel against disba:

```bash
make compare
```

This uses `/opt/miniconda3/envs/dispa/bin/python`, where `disba` is available on this machine. Prefer this `dispa` environment, or `/opt/miniconda3/envs/seisloc/bin/python` when a newer package stack is needed, over the system `python3`. For Numba/disba/matplotlib scripts, set writable caches when needed:

```bash
NUMBA_CACHE_DIR=/private/tmp/numba_cache
MPLCONFIGDIR=/private/tmp/mplconfig
```

It writes:

- `disba_compare.out`: direct comparison between Fortran `surfdisp96` and disba finite-difference kernels using the same DSurfTomo node perturbation and sublayer mapping.
- `disba_layer_kernels.out`: disba's native layer-parameter sensitivity kernels. These are useful diagnostics but are not the same parameterization as DSurfTomo's node kernels.

To test the pDSurfTomo-style CPU batch finite-difference route:

```bash
make batch-compare
```

This keeps DSurfTomo's finite-difference kernel definition, but builds all `6*nz` perturbation
models first and runs the independent disba forward jobs with a thread pool. It is a prototype
for moving the same task decomposition back into the production Fortran code.

To test the actual pDSurfTomo Numba `surf96_batch` backend:

```bash
make pdsurf-numba-compare
```

Note that pDSurfTomo calls this path `disba`, but the local implementation is
`pDSurfTomo-main/bin/surf96_disba.py`: a Numba-vectorized surf96 translation,
not the pip package named `disba`.

To benchmark memory-bounded tile processing over real `MOD_Ref` columns:

```bash
make pdsurf-tile-benchmark
```

The script `pdsurf_tile_senk_solver.py` uses the same tile strategy but exposes
Fortran stream binary I/O for DAzi integration.

The preferred accelerated DAzi backend is:

```bash
DAZI_DEPTHKERNEL_BACKEND=pdsurf_numba_tight_root_fortran_pv
DAZI_PDSURF_ENGINE=tight_root
```

It uses pDSurf Numba for depth sensitivities, tightens the pDSurf root
convergence, and recomputes `pvRc` with the Fortran phase-velocity forward
solver, so FMM ray tracing stays aligned with the default backend.

The validated TI process-batch backend is:

```bash
DAZI_DEPTHKERNEL_TI_BACKEND=process_batch
```

It launches independent Fortran worker processes around the original
`surfdisp96/tregn96` path, avoiding unsafe OpenMP over `COMMON` blocks.

Fast combined run:

```bash
DAZI_DEPTHKERNEL_BACKEND=pdsurf_numba_tight_root_fortran_pv \
DAZI_PDSURF_ENGINE=tight_root \
DAZI_DEPTHKERNEL_TI_BACKEND=process_batch \
DAZI_DEPTHKERNEL_TI_WORKERS=4 \
DAZI_DEPTHKERNEL_TI_TILE_COLUMNS=64 \
DAZI_PDSURF_TILE_COLUMNS=64 \
DAZI_PYTHON=/opt/miniconda3/envs/dispa/bin/python \
/Users/liuxin/Desktop/DSurf_test/DAzi_vorotomo_migration/bin/DAzimSurfTomo para.in
```

Prefer `run_dazi_configured.py` for normal runs. It also exports the repo-local
`pdsurf_tile_senk_solver.py`, `depthkernel_ti_process_batch.py`, and
`depthkernel_ti_worker` paths so a run does not silently pick up stale helper
scripts from another checkout.

For `ani_voro` runs where strict equivalence to the original Fortran baseline
matters, use the stable anisotropic profile instead:

```bash
python3 run_dazi_configured.py /path/to/case \
  --profile strict_ani_voro \
  --threads 4 \
  --log /tmp/dazi_strict_ani.log \
  --metadata /tmp/dazi_strict_ani.json
```

This keeps the original isotropic finite-difference sensitivity kernel and
uses only the exact-validated acceleration pieces:

```text
depthkernel_backend = default
depthkernel_ti_backend = process_batch
voro_gp_backend = direct
```

On the `ani_voro, nreal=10, p10_s50` gate this matched the original baseline
exactly and reduced runtime only from `161.1 s` to `149.7 s`; it is a stable
baseline profile, not the high-speed profile.

For this Mac, `DAZI_DEPTHKERNEL_TI_WORKERS=8` or `10` is faster in standalone
p10_s50 tests and still matches the serial TI dump exactly.

Convenience runner:

```bash
./run_dazi_hybrid_backend.sh /path/to/case 4
```

For the server handoff package, keep normal user-facing runtime controls in
the tail of each case's `para.in`.  The unified runner also supports explicit
legacy key/value configs when `--config` is passed.  In that case, copy one of
these files to `/path/to/case/dazi_run.conf`:

```text
dazi_run.strict_ani_voro.conf.example
dazi_run.fast_exploratory.conf.example
dazi_run.fast_tight_root.conf.example
dazi_run.conf.example
```

Set either one default thread count:

```text
threads = 8
```

or separate counts when CPU pressure needs tighter control:

```text
omp_threads = 8
ti_workers = 4
```

Then run:

```bash
python3 run_dazi_configured.py /path/to/case --config /path/to/case/dazi_run.conf --metadata /tmp/dazi_run_meta.json
```

Use `--dry-run` to print the exact binary, input file, and environment mapping
without starting the inversion.

Available backend profiles:

| profile | backend meaning | intended use |
|---|---|---|
| `legacy_baseline` | default depthkernel, default TI, legacy GP | original-reference comparisons |
| `strict_ani_voro` | default depthkernel, process-batch TI, direct GP | stable `ani_voro` production checks |
| `fast_exploratory` | previous pDSurf sensitivity backend, process-batch TI, direct GP | old-engine regression and iso checks |
| `fast_tight_root` | tight-root pDSurf sensitivity backend, process-batch TI, direct GP | default fast line for development |

To run the four acceptance combinations with the same CPU/backend settings,
either copy `acceptance_matrix.conf.example` and fill in the four case paths,
or generate four matched cases from one template:

```bash
python3 make_four_combo_cases.py \
  --source-case /private/tmp/dazi_hybrid_p10s50_ab_20260506/template \
  --out-root /private/tmp/dazi_four_combo_cases_p10s50 \
  --maxiter 1 \
  --nrealizations 1 \
  --normal-cells 10 \
  --adaptive-cells 5 \
  --nzrand 4 \
  --force
```

Then run:

```bash
python3 run_four_combo_acceptance.py \
  --matrix-config /private/tmp/dazi_four_combo_cases_p10s50/acceptance_matrix.conf \
  --out-root /private/tmp/dazi_four_combo_acceptance \
  --max-rss-gb 32 \
  --dry-run
```

Remove `--dry-run` to start the actual jobs:

```bash
python3 run_four_combo_acceptance.py \
  --matrix-config /private/tmp/dazi_four_combo_cases_p10s50/acceptance_matrix.conf \
  --out-root /private/tmp/dazi_four_combo_acceptance \
  --max-rss-gb 32 \
  --force
```

Add `--with-baseline` when you want each combo to run both the legacy/default
backend and the optimized backend and report output differences.

For anisotropic output, use the dedicated angle-aware comparator instead of raw
per-column differences:

```bash
python3 compare_anisotropy_model.py \
  /path/to/baseline/run \
  /path/to/optimized/run \
  --out-md /tmp/anisotropy_diff.md \
  --out-json /tmp/anisotropy_diff.json
```

It compares the fast-axis angle modulo 180 degrees and reports `Gc/Gs` vector
differences separately from formatted angle/amplitude columns. The four-combo
acceptance summary automatically includes these anisotropy-aware metrics when
`--with-baseline` is used.

To build a bounded 5-period medium A/B case from the Anning real-data example:

```bash
/opt/miniconda3/envs/dispa/bin/python make_medium_ab_case.py \
  --source-case /Users/liuxin/Desktop/DSurf_test/DAzi_large_data_260128/example/anning_real_SRL_ani_cell30_15_nreal050_iaratio0p5 \
  --out-case /private/tmp/dazi_hybrid_medium_ab_20260506/template \
  --max-period 5 \
  --sources-per-period 30
```

The same generator can build a larger 10-period stress case:

```bash
/opt/miniconda3/envs/dispa/bin/python make_medium_ab_case.py \
  --source-case /Users/liuxin/Desktop/DSurf_test/DAzi_large_data_260128/example/anning_real_SRL_ani_cell30_15_nreal050_iaratio0p5 \
  --out-case /private/tmp/dazi_hybrid_p10s50_ab_20260506/template \
  --max-period 10 \
  --sources-per-period 50
```

In this larger case, compare both final output files and dumped G triplets.
The final anisotropy model can be sensitive to small Vs-kernel changes even
when Gc/Gs triplet columns are identical.

To compare two completed DAzi runs produced by different depth-kernel backends:

```bash
/opt/miniconda3/envs/dispa/bin/python compare_dazi_backend_outputs.py \
  --base /private/tmp/dazi_pdsurf_backend_ab_20260506/default \
  --test /private/tmp/dazi_pdsurf_backend_ab_20260506/pdsurf_numba \
  --outdir /private/tmp/dazi_pdsurf_backend_ab_20260506
```

This writes `log_summary.csv` and `file_diff_summary.csv`.

To compare dumped depth-kernel arrays before G assembly:

```bash
/opt/miniconda3/envs/dispa/bin/python compare_depthkernel_dumps.py \
  --base /private/tmp/dazi_pdsurf_backend_ab_20260506/default/depthkernel_dump \
  --test /private/tmp/dazi_pdsurf_backend_ab_20260506/pdsurf_numba/depthkernel_dump \
  --nx 28 --ny 38 --nz 20 --kmax 1 \
  --vel-bin /private/tmp/dazi_pdsurf_backend_ab_20260506/pdsurf_numba/.dazi_senk_tile/vel_senK.bin \
  --out /private/tmp/dazi_pdsurf_backend_ab_20260506/depthkernel_dump_diff_summary.csv
```

To compare dumped sparse G triplets from the joint anisotropic assembly:

```bash
/opt/miniconda3/envs/dispa/bin/python compare_g_triplet_dumps.py \
  --base /private/tmp/dazi_pdsurf_backend_ab_20260506/default/g_triplet_dump \
  --test /private/tmp/dazi_pdsurf_backend_ab_20260506/pdsurf_numba/g_triplet_dump \
  --nparpi 17784 \
  --out /private/tmp/dazi_pdsurf_backend_ab_20260506/g_triplet_diff_summary.csv
```

To run a complete iso depth-kernel A/B gate:

```bash
/opt/miniconda3/envs/dispa/bin/python run_iso_depthkernel_ab.py \
  --case /private/tmp/dazi_hybrid_p10s50_ab_20260506/template \
  --out-root /private/tmp/dazi_iso_closure_20260506/p10s50 \
  --case-label p10_s50 \
  --omp-threads 10 \
  --tile-columns 64 \
  --force
```

To run a complete TI depth-kernel A/B gate:

```bash
/opt/miniconda3/envs/dispa/bin/python run_ti_depthkernel_ab.py \
  --case /private/tmp/dazi_hybrid_p10s50_ab_20260506/template \
  --out-root /private/tmp/dazi_ti_closure_20260506/p10s50 \
  --case-label p10_s50 \
  --omp-threads 10 \
  --ti-workers 4 \
  --ti-tile-columns 64 \
  --force
```

Closure notes:

- `/Users/liuxin/Desktop/DSurf_test/iso_depthkernel_closure.md`
- `/Users/liuxin/Desktop/DSurf_test/ti_depthkernel_closure.md`
- `/Users/liuxin/Desktop/DSurf_test/lsmr_vorotomo_backend_notes.md`
- `/Users/liuxin/Desktop/DSurf_test/voro_gp_direct_closure.md`

To run a Fortran/SciPy LSMR backend A/B check:

```bash
/opt/miniconda3/envs/dispa/bin/python run_lsmr_backend_ab.py \
  --case /private/tmp/dazi_ti_closure_20260506/mini/process_batch \
  --out-root /private/tmp/dazi_lsmr_ab_20260506/mini_voro \
  --case-label mini_voro \
  --omp-threads 4 \
  --ti-workers 4 \
  --backend-scope voro \
  --voro-timer \
  --force
```

Use this as a solver sandbox, not as the default production backend yet. The
mini A/B showed that one Python process per Vorotomo realization adds more
overhead than native Fortran LSMR on small cases.

The optional `--voro-timer` flag sets `DAZI_VORO_TIMER=1` and writes per-realization
breakdown lines for cell generation, grid-to-cell mapping, projected `GP`
assembly, LSMR, and back projection.

To run the legacy/direct Vorotomo `GP` assembly A/B:

```bash
/opt/miniconda3/envs/dispa/bin/python run_voro_gp_backend_ab.py \
  --case /private/tmp/dazi_ti_closure_20260506/p10s50/process_batch \
  --out-root /private/tmp/dazi_voro_gp_ab_20260506/p10s50 \
  --case-label p10s50 \
  --omp-threads 4 \
  --ti-workers 4 \
  --voro-timer \
  --force
```

The current validated default backend is direct GP. To force it explicitly:

```bash
DAZI_VORO_GP_BACKEND=direct
```

To compare against the old dense `grow + aprod` projection path:

```bash
DAZI_VORO_GP_BACKEND=legacy
```

The real `nrealizations=50` one-iteration A/B is recorded at:

```text
/private/tmp/dazi_voro_gp_ab_20260506/real_ani_nreal50_iter1/summary_voro_gp_backend_ab.md
```

The full 5-iteration direct-GP production-style run is recorded at:

```text
/Users/liuxin/Desktop/DSurf_test/production_direct_gp_run_20260506.md
```

To sweep finite-difference step size and sublayer refinement:

```bash
make sweep
```

The first sweep scans:

- `dln = 0.02, 0.01, 0.005, 0.0025`
- `minthk = 1, 2, 4, 8`

It writes `sweep_summary.csv` plus per-case baseline/compare files.

To extract representative columns from the Anning real-data MOD file and run the same sweep:

```bash
make sweep-real
```

This writes results under `real_column_sweep/`, including:

- `columns_manifest.csv`
- `real_columns_sweep_summary.csv`

To diagnose where the pDSurf sensitivity backend differs from the original
Fortran finite-difference backend:

```bash
/opt/miniconda3/envs/dispa/bin/python diagnose_pdsurf_sensitivity_diff.py \
  --base-depth /path/to/default/depthkernel_dump \
  --test-depth /path/to/hybrid/depthkernel_dump \
  --base-triplets /path/to/default/g_triplet_dump \
  --test-triplets /path/to/hybrid/g_triplet_dump \
  --vel-bin /path/to/hybrid/.dazi_senk_tile/vel_senK.bin \
  --depz-bin /path/to/hybrid/.dazi_senk_tile/depz.bin \
  --periods-bin /path/to/hybrid/.dazi_senk_tile/tRc.bin \
  --nx 28 --ny 38 --nz 20 --kmax 10 \
  --out-dir /tmp/pdsurf_sensitivity_diagnostics
```

This writes period/depth bucket summaries for depth kernels and block/depth
bucket summaries for sparse G triplets, including near-`ftol` counts.

To isolate whether the mismatch already exists in the forward engine, build the
production `surfdisp96` batch driver:

```bash
make surfdisp96_forward_batch
```

Then compare pDSurf `surf96_batch` with the production Fortran `surfdisp96`
over the same refined baseline and perturbation models:

```bash
/opt/miniconda3/envs/dispa/bin/python compare_pdsurf_forward_engine.py \
  --vel-bin /path/to/hybrid/.dazi_senk_tile/vel_senK.bin \
  --depz-bin /path/to/hybrid/.dazi_senk_tile/depz.bin \
  --periods-bin /path/to/hybrid/.dazi_senk_tile/tRc.bin \
  --nx 28 --ny 38 --nz 20 \
  --columns 1,532,1064 \
  --minthk 3.0 --iflsph 1 --iwave 2 --igr 0 --dln 0.01 \
  --dc 0.005 --dt 0.005 \
  --fortran-exe ./surfdisp96_forward_batch \
  --out-dir /tmp/pdsurf_forward_parity
```

This writes per-job and grouped forward-curve difference tables plus
`FORWARD_ENGINE_PARITY.md`. On the p10_s50 diagnostic columns, pDSurf and
Fortran differed by about `5e-6 km/s` in phase velocity. That is small for a
single forward curve but large enough to create `~1e-4` sensitivity differences
after finite differencing.

By default the comparator also tests precision variants:

```text
default
model_f32
output_f32
model_output_f32
strict_fastmath
tight_root
strict_fastmath_tight_root
```

These do not change the production backend. They only test whether pDSurf's
double-precision path differs from DAzi's Fortran `real*4/sngl` precision path.
The standalone tile solver exposes the same diagnostic controls as
`--model-precision`, `--output-precision`, and `--engine`; the default engine
is unchanged:

```text
--engine default
```

The comparator also exposes `--iflsph 0/1` to separate flat-earth forward
differences from spherical-transform differences. The optional
`strict_fastmath` and `tight_root` variants write temporary pDSurf module
copies in the output directory with source-level diagnostic transforms; they
are diagnostic only.

For a meaningful horizontally varying model, run the same tool against `MOD_Ref`:

```bash
NUMBA_CACHE_DIR=/private/tmp/numba_cache /opt/miniconda3/envs/dispa/bin/python sweep_real_columns.py \
  --mod /Users/liuxin/Desktop/DSurf_test/DAzi_large_data_260128/example/anning_real_SRL_ani_cell30_15_nreal050_iaratio0p5/MOD_Ref \
  --outdir real_column_sweep_modref
```

To summarize and plot a sweep table:

```bash
make summarize
```

This writes `best_by_column.csv`, `best_overall.txt`, and heatmaps under `real_column_sweep_modref/summary_plots/`.

## Input Format

`sample_column.in` is intentionally simple:

```text
nz kmax
depz(1:nz)
vs(1:nz)
periods(1:kmax)
minthk iwave igr
```

`Vp` and `rho` are computed from the same empirical relations used in DSurfTomo.

## Next Checks

1. Compare this output with the full `CalSurfG.depthkernel` for the same column.
2. Add a Python/disba comparator for phase velocity and sensitivity shape.
3. Add step-size tests for `dlnVs`, `dlnVp`, and `dlnrho`.
4. Add profiling around `refineGrid2LayerMdl` and `surfdisp96` call counts.
