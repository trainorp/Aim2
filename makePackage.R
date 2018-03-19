library(devtools)
library(roxygen2)
library(statmod)

setwd("~/gdrive/Dissertation/Aim2/")

# devtools::create("BayesianGLasso")
# devtools::create_description("BayesianGLasso")
devtools::use_gpl3_license(pkg="BayesianGLasso")
devtools::use_package("statmod",type="Imports",pkg="BayesianGLasso")
devtools::use_package("MASS",type="Imports",pkg="BayesianGLasso")
devtools::use_package("glasso",type="Imports",pkg="BayesianGLasso")
devtools::use_rcpp("BayesianGLasso")
# devtools::use_package("Rcpp",type="Imports",pkg="BayesianGLasso")
# devtools::use_package("Rcpp",type="LinkingTo",pkg="BayesianGLasso")
devtools::use_package("RcppArmadillo",type="LinkingTo",pkg="BayesianGLasso")
devtools::document("BayesianGLasso")
devtools::use_cran_comments("BayesianGLasso")

currentCode<-as.package("BayesianGLasso")

devtools::load_all(currentCode)
devtools::build(currentCode)
devtools::build_win(currentCode)

devtools::release("BayesianGLasso")
