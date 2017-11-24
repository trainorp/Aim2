library(devtools)
library(roxygen2)
library(statmod)

setwd("~/gdrive/Dissertation/Aim2/")

devtools::create("BayesianGLasso")
devtools::use_gpl3_license("BayesianGLasso")
devtools::use_package("statmod",type="Imports",pkg="BayesianGLasso")
devtools::use_package("MASS",type="Imports",pkg="BayesianGLasso")
devtools::use_package("glasso",type="Imports",pkg="BayesianGLasso")
devtools::use_cran_comments("BayesianGLasso")
roxygen2::roxygenise(package.dir="~/gdrive/Dissertation/Aim2/BayesianGLasso")

currentCode<-as.package("BayesianGLasso")
devtools::load_all(currentCode)
devtools::document(currentCode)
devtools::build(currentCode)
devtools::build_win(currentCode)
devtools::release("BayesianGLasso")