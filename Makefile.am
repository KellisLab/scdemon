# Makefile.am

src/RcppExports.cpp:
	Rscript -e "Rcpp::compileAttributes()"
R/RcppExports.R:
	Rscript -e "Rcpp::compileAttributes()"

document: src/RcppExports.cpp R/RcppExports.R
	Rscript -e "devtools::document(roclets = c('rd', 'collate', 'namespace'))"

buildbin: src/RcppExports.cpp R/RcppExports.R
	Rscript -e "devtools::build(binary = TRUE, args = c('--preclean'))"

buildsrc: src/RcppExports.cpp R/RcppExports.R
	Rscript -e "devtools::build(args=c('--no-build-vignettes'))"

install: src/RcppExports.cpp R/RcppExports.R
	Rscript -e "devtools::install_local(force=TRUE, upgrade=FALSE)"

test:
	Rscript -e "devtools::test()"

site:
	Rscript -e "pkgdown::build_site(preview=FALSE)"

AM_CXXFLAGS = @EIGEN_CXXFLAGS@ @OPENMP_CXXFLAGS@ @HDF5_CXXFLAGS@

if ENABLE_CXX_TEST
	bin_PROGRAMS += bin/scdemon
	bin_scdemon_SOURCES = bin/scdemon.cpp
	bin_scdemon_LDADD = @EIGEN_LIBS@ @OPENMP_LIBS@ @HDF5_LIBS@
endif

all-local: document buildsrc buildbin install test site

clean-local:
	rm -f src/*.o src/*.so
	rm -f src/RcppExports.cpp R/RcppExports.R
	rm -f bin/scdemon bin/scdemon.o
