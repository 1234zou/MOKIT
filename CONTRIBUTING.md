## A few suggestions for contributors:

Please make sure that the diff stat is clean and all irrelevant files are ignored by `.gitignore`. Tips:
* It's fine to compile in `src` and send merge request, since `*.o`, `*.so`, etc are ignored in `.gitignore`.
* If the current Makefile ( a link to `Makefile.intel`) does not satisfy you, you can try `Makefile.gnu_openblas` or `Makefile.gnu_openblas_arch`. But do not modify them directly.
* If these three still don't work for you, you can modify one of them
```
cp Makefile.intel xxx.make
# after modifying xxx.make
make all -f xxx.make
```
since `*.make` is also ignored.
* do not run any examples in `examples` since the `*.gbw`, `*.gms`, etc are not ignored yet. Please run the examples in any other place.

## CI/CD
* any push or pr will trigger Gitlab CI build. The produced artifacts (binary) can be downloaded at pipelines.
* To skip CI (when updating docs or something not important), add [ci skip] in commit message. For example `git commit -m "xxx [ci skip]"`. 
* CI minutes is limited to 400min per user per repo per month. One pass of MOKIT CI will consume around 30 min.
