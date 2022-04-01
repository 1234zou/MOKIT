A few suggestions for contributors:
* it's fine to compile in `src` and send merge request, since `*.o`, `*.so`, etc are ignored in `.gitignore`.
* if you want to modify `src/Makefile` and contribute, you can do
```
cp Makefile xxx.make
# after modifying xxx.make
make all -f xxx.make
```
since `*.make` is also ignored.
* do not run any examples in `examples` since the `*.gbw`, `*.gms`, etc are not ignored yet. Please run the examples in any other place.

