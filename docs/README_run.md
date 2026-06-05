# DAzi Current Release Run Notes

This package contains a stable DAzi current-release code snapshot, a small Yunnan example, and simple shell scripts for compiling and running cases.

## 1. Build

From the package root:

```bash
chmod +x build.sh run_one_case.sh
./build.sh
```

The main executables are:

```text
code/bin/DAzimSurfTomo
code/bin/SurfAAForward
```

## 2. Run A Case Directly

The Fortran program still supports the original direct `para.in` workflow:

```bash
cd examples/test4_Yunnan/voro
../../code/bin/DAzimSurfTomo para.in
```

From the package root, the helper script does the same thing and only adds thread setup plus logging:

```bash
./run_one_case.sh examples/test4_Yunnan/voro 4
```

The second argument is the thread count. If it is omitted, the script uses `DAZI_THREADS`, then falls back to `4`.

## 3. Included Examples

```text
examples/test4_Yunnan/smooth   fixed-grid / no-Voro reference
examples/test4_Yunnan/voro     Voro example
examples/test4_Yunnan/voro50   lower-cell Voro example
examples/test4_Yunnan/voro200  higher-cell Voro example
```

Each example contains the needed `para.in`, `MOD`, and dispersion data file.

## 4. Important Voro Note

The current Voro inversion does not use the traditional fixed-grid weight in the same way as the NoVoro/smooth workflow. The program prints:

```text
Please note that the voronoi cell scheme does not use weight.
```

For uneven ray coverage, especially when long-period paths are sparse, first check the path distribution and consider removing poorly constrained periods before interpreting anisotropy.

## 5. Practical First Test

For a new dataset:

1. Start with NoVoro or conservative Voro settings.
2. Use only periods with enough path coverage.
3. Use fewer iterations first, for example 3-5.
4. Increase Voro cell count only after the residual and model are stable.

