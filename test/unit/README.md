# Unit Tests

Build and run:

```
cd test/unit
make runtest
```

Tests:

- `test_read_bgchg_from_gjf`: validates background charge parsing.
- `test_check_kywd_compatible_pass1`: valid default keyword set.
- `test_check_kywd_compatible_pass2`: X2C with CASSCF on ORCA, HF_prog auto-switch.
- `test_check_kywd_compatible_pass3`: valid CASCI with ORCA.
- `test_check_kywd_compatible_xfail`: invalid NewMult with ist!=5.
- `test_check_kywd_compatible_quote_xfail`: single quotes in mcpdft_prog.

Expected-fail tests:

- Any test binary named with the suffix `_xfail` is treated as expected to fail.

Choose a toolchain:

```
make runtest -f Makefile.gnu_openblas
make runtest -f Makefile.gnu_openblas_ci
make runtest -f Makefile.gnu_mkl
```
