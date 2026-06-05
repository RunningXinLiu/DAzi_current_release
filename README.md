# DAzi Current Release

Stable DAzi surface-wave azimuthal anisotropy tomography package with Voronoi parameterization, forward modeling utilities, and runnable Yunnan examples.

This repository is intended as a clean public-facing release snapshot: it keeps the core Fortran inversion/forward code, compact examples, and practical run notes in one place.

## What This Code Does

DAzi estimates 3-D isotropic shear-wave velocity perturbations and azimuthal anisotropy from Rayleigh-wave phase/travel-time observations.

The current release supports:

- Joint inversion for isotropic `dVs` and azimuthal anisotropy terms `Gc/Gs`
- NoVoro fixed-grid inversion examples
- Voronoi-cell based inversion examples
- Forward modeling with `MODVs.true`, `MODGc.true`, and `MODGs.true`
- Expanded `para.in` runtime controls for Voro, FMM, LSMR, and thread settings
- Direct legacy-style execution through `DAzimSurfTomo para.in`

## Repository Layout

```text
code/
  src_inv_iso_joint/       main DAzi inversion source
  src_forward/             forward modeling source
  depthkernel_lab/         depth-kernel and backend utilities
  bin/                     prebuilt executables included for convenience

examples/
  test4_Yunnan/
    smooth/                current-release NoVoro example
    voro/                  current-release Voro example
  test4_Yunnan_legacy_old_para/
    ...                    old compact para.in examples, kept only as reference

docs/
  README_run.md            build and run guide
  VORO_PARAMETER_GUIDE.md  Voro parameter notes
  PARA_IN_DIRECT_MODE.md   direct para.in execution notes

scripts/
  convert_inv_to_modtrue.py
```

## Quick Start

Clone and build:

```bash
git clone https://github.com/RunningXinLiu/DAzi_current_release.git
cd DAzi_current_release
chmod +x build.sh run_one_case.sh
./build.sh
```

Run the current-release Voro example:

```bash
./run_one_case.sh examples/test4_Yunnan/voro 4
```

Or run directly in the original DAzi style:

```bash
cd examples/test4_Yunnan/voro
../../../code/bin/DAzimSurfTomo para.in
```

The output log is written to `para.in_inv.log` when using `run_one_case.sh`.

## Example Cases

| Case | Purpose | Notes |
|---|---|---|
| `examples/test4_Yunnan/smooth` | NoVoro baseline | Fixed-grid reference case with expanded current-release `para.in` |
| `examples/test4_Yunnan/voro` | Voro reference | Voronoi parameterized current-release case |
| `examples/test4_Yunnan_legacy_old_para` | Historical reference | Old compact `para.in` format, not recommended as the default entry point |

The current example `para.in` files include explicit blocks for:

- output controls
- Voronoi parameters
- `gp_spfra`
- runtime profile
- ray backend
- LSMR backend
- unified thread count
- depth-kernel tile controls

## Voro Notes

The current Voro workflow does **not** apply the traditional fixed-grid weight in the same way as the NoVoro workflow. During Voro inversion, the program prints:

```text
Please note that the voronoi cell scheme does not use weight.
```

This matters for uneven ray coverage. If long-period paths are sparse, first inspect path coverage and consider removing weakly constrained periods before increasing the number of Voro cells.

Practical first-pass advice:

- Start with conservative Voro cell counts
- Use fewer iterations first, for example 3-5
- Remove periods with very poor path coverage
- Increase `ncell`, `acell`, or `nreal` only after the residual and model are stable
- Treat poorly covered areas as low-resolution or null-space regions

## Important `para.in` Parameters

Typical current-release Voro block:

```text
1                                    c: voronoi mode (1: use; 0: not use)
100                                  c: Vs normal/random voro cell number per layer
0                                    c: Vs adaptive voro cell number per layer
100                                  c: number of realizations
9                                    c: voronoi cell layer number (nzrand)
0.5                                  c: iaratio = ANI cell count / Vs cell count
```

Optional runtime controls:

```text
fast_tight_root                      c: profile [legacy_baseline|strict_ani_voro|fast_exploratory|fast_tight_root]
fmm_parallel                         c: ray_backend [serial|fmm_parallel]
scipy                                c: lsmr_backend [fortran|scipy|unset]
4                                    c: threads, unified OpenMP/process CPU budget; DAZI_THREADS overrides
64                                   c: pdsurf_tile_columns
64                                   c: ti_tile_columns
1                                    c: voro_timer
```

These controls are designed to make the run configuration explicit in `para.in` while still allowing direct Fortran execution.

## Forward Modeling Utility

To convert inversion outputs into forward-model truth files:

```bash
python3 scripts/convert_inv_to_modtrue.py examples/test4_Yunnan/voro --force
```

This generates:

```text
MODVs.true
MODGc.true
MODGs.true
```

These files can be used with the forward executable `code/bin/SurfAAForward`.

## Recommended Workflow For A New Dataset

1. Prepare dispersion/travel-time data and an initial `MOD`.
2. Start with a NoVoro or conservative Voro case.
3. Check period-by-period path coverage.
4. Remove extremely sparse long-period data before interpreting anisotropy.
5. Run a short iteration test.
6. Sweep a small number of Voro settings only after the baseline is stable.
7. Plot residuals and mask low-coverage regions before geological interpretation.

## Version

Current public snapshot:

```text
v2026.06.05-current-para
```

This tag corresponds to the release package with expanded current-release `para.in` examples.

## Notes

This repository is a research code release. It is organized for reproducible collaboration and example-driven use, not as a polished software library.

