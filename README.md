# mixstock

This is the source code repository for the `mixstock` package for mixed stock analysis in R, based on work by Bolker and Okuyama and others.

The package is not on CRAN, and hasn't been for some time, because of fairly trivial problems with some of the C code. (In particular, there are some indexing problems somewhere in the C code for unconditional and conditional maximum likelihood (CML/UML) that trigger `valgrind` warnings. If this means nothing to you, don't worry about it.)

If you want to install the package, at present you will need to have the `devtools` package, *and the R compilation tools*, installed. Then you can either use `devtools::install_version("mixstock","0.9.5.1")` (to install the last archived version from CRAN) or `devtools::install_github("bbolker/mixstock")` (to install the version from this repository, which may be more up-to-date).

If you really need the package and are unable to install compilation tools for some reason, you'll need to contact me and convince me to build a binary version of the package for you.

If anyone is interested in helping to bring this package up-to-date and getting it back on CRAN, please let me know! This could be done either by (1) fixing the problems with the C code or (2) removing the C code (and the corresponding R functions) from the package.


