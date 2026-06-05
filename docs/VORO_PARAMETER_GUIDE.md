# Voro Parameter Guide

The Voro parameters are in `para.in`.

```text
voronoi mode
Vs normal/random voro cell number per layer
Vs adaptive voro cell number per layer
number of realizations
voronoi cell layer number (nzrand)
iaratio = ANI cell count / Vs cell count
```

## Recommended Conservative Usage

For a small or unevenly covered dataset, begin conservatively:

```text
voronoi mode = 1
ncell        = 20-60
acell        = 0 or about 10-20% of ncell
nreal        = 20-50
nzrand       = number of effective depth intervals
iaratio      = 0.5-0.7 for anisotropic inversion
```

If long-period path counts are very small, do not compensate by increasing cell count. Remove or down-select those periods first.

## How To Read Residual Behavior

If the first iteration reduces residual but later iterations rise slightly or plateau, it usually means the useful linear update has already happened and later iterations are fitting weakly constrained structure or noise.

If residual increases strongly from the first iteration, check:

```text
period coverage
data format and units
initial MOD
too many Voro cells
too many poorly constrained long-period paths
anisotropic smoothing too weak
```

## About Weight

The current Voro path does not apply the traditional fixed-grid weight as a simple spatial smoothing mask. Treat coverage QC and period selection as part of the inversion setup.

