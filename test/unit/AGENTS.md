# Unit Test Coding Guide for AI Agents

Scope

- This guide is only for unit tests under `test/unit/`.
- Follow existing Fortran test patterns and naming conventions.

Where to add tests

- Put tests for `mr_keyword` under `test/unit/mr_keyword_check_kywd/` or `test/unit/mr_keyword_read_bgchg/`.
- Add new test directories if needed, and keep file names clear and specific.

Naming conventions

- Fortran tests: `test_<feature>_<case>.f90`
- Expected-fail Fortran tests: `test_<feature>_<case>_xfail.f90`
- Python tests: `test_<feature>.py`

Runner behavior

- `run_tests.sh` treats `*_xfail` binaries as expected failures.
- Non-`*_xfail` Fortran tests must exit 0 to pass.
- Python tests are run with pytest.

Test utilities

- Use `test_utils.f90` for common helpers.
- For keyword tests, call `reset_mr_keyword_defaults()` to avoid duplicated setup.

Makefile integration

- Add new test binaries to `TESTS` in `test/unit/Makefile.main`.
- Add object lists if needed (see `OBJ_check_kywd_common`).
- Ensure any new helper modules are linked in the relevant `OBJ_*` list.

Running tests

Fortran tests:
```
cd test/unit
make runtest
```

Python tests:
```
cd test/unit
pytest test_*.py
```

Toolchains

```
make runtest -f Makefile.gnu_openblas
make runtest -f Makefile.gnu_mkl
```
