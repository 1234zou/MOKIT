
# install doxygen and graphviz first

doxygen -g

# Edit Doxyfile to enable Fortran and call graphs
sed -i 's/OPTIMIZE_FOR_FORTRAN.*= NO/OPTIMIZE_FOR_FORTRAN = YES/' Doxyfile
sed -i 's/EXTRACT_ALL.*= NO/EXTRACT_ALL = YES/' Doxyfile
sed -i 's/CALL_GRAPH.*= NO/CALL_GRAPH = YES/' Doxyfile
sed -i 's/CALLER_GRAPH.*= NO/CALLER_GRAPH = YES/' Doxyfile
sed -i 's/RECURSIVE.*= NO/RECURSIVE = YES/' Doxyfile
sed -i 's/HAVE_DOT.*= NO/HAVE_DOT = YES/' Doxyfile
sed -i 's/GENERATE_LATEX.*= YES/GENERATE_LATEX = NO/' Doxyfile

# Run doxygen on your source files
doxygen Doxyfile

# firefox html/index.html 
