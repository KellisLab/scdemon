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


all-local: loadall document buildsrc buildbin install test site

clean-local:
	rm -f src/*.o src/*.so
	rm -f src/RcppExports.cpp R/RcppExports.R
