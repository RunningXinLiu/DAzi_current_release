# Direct `para.in` Mode

The automatic wrappers are optional. They are useful for:

```text
preflight memory checks
thread setup
automatic Voro cell recommendation
case generation
consistent logs
batch execution
```

But the inversion executable itself still accepts the original command:

```bash
DAzimSurfTomo para.in
```

This means a user can manually edit `para.in` and run a case exactly like the early DAzi versions.

The helper script in this package:

```bash
./run_one_case.sh examples/test4_Yunnan/voro 4
```

is only a convenience wrapper around:

```bash
cd examples/test4_Yunnan/voro
OMP_NUM_THREADS=4 ../../code/bin/DAzimSurfTomo para.in
```

