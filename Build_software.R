### build the website
library("workflowr")
wflow_status()
wflow_view()
wflow_build()
wflow_publish()
# wflow_git_commit(all=TRUE)   You can use gitkraken anyway


### build the R package
require(devtools)
#devtools::missing_s3()  #list all S3 methods that you’ve forgotten to export
devtools::document()
# devtools::build_vignettes()
devtools::build()
devtools::install(build_vignettes = TRUE)

library("SlopeHunter")
data("data_example")
head(data_example)
?data_example
citation("SlopeHunter")

Sh.Model <- slopehunter(dat = data_example, xbeta_col="xbeta", xse_col="xse",
                        ybeta_col="ybeta", yse_col="yse", yp_col="yp",
                        comp.size = seq(0.03, 0.10, 0.01), xp.thresh = 0.1, coef.diff = 1,
                        correct.reg.dill = TRUE, show_adjustments = TRUE, seed=2019)

Sh.Model$Sh.b
Adj <- Sh.Model$Estimates
head(Adj)
head(Adj$ybeta)
head(Adj$ybeta.Adj)

Adj_sh <- SHadj(Sh.Model, dat = data_example, xbeta_col = "xbeta", xse_col = "xse",
             ybeta_col = "ybeta", yse_col = "yse")

head(Adj_sh)
?slopehunter
plot(Sh.Model)

# Install it from Github
devtools::install_github("Osmahmoud/SlopeHunter")
require(SlopeHunter)
