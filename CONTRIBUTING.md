## Guides for submitting an issue 
* Bug report: please provide following information, if applicable and available
  + version of MOKIT (1.x.y rc n), type of installation (conda, source, etc.)
  + error message or unexpected output
  + input files, command line input, steps for reproduce
* Feature requests, suggestions, proposals, documentation requests are also welcome.

## A few suggestions for contributors

Please make sure that the diff stat is clean and all irrelevant files are ignored by `.gitignore`. Tips:
* It's fine to compile in `src` and send merge request, since `*.o`, `*.so`, etc are ignored in `.gitignore`.
* If the current Makefile ( a link to `Makefile.intel`) does not satisfy you, you can try `Makefile.gnu_openblas` or `Makefile.gnu_mkl`. But do not modify them directly.
* If all of them still don't work for you, you can modify one of them
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
* CI minutes is limited to 400min per user per repo per month. One pass of MOKIT CI without conda build will consume around 20 min.
* Please don't modify the conda recipes, e.g. the rc number or build number, unless you do want to improve this part. The conda packaging will be taken care of by maintainers.
