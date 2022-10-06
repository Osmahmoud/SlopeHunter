### build the website
library("workflowr")
wflow_status()
wflow_view()
wflow_build()
# wflow_publish()  # No need if you have commited your changes e.g. using gitkraken
# wflow_git_commit(all=TRUE)   You can use gitkraken anyway


### build the R package
require(devtools)
devtools::missing_s3()  #list all S3 methods that you’ve forgotten to export
devtools::document()
devtools::build_vignettes()
devtools::build()
devtools::build(binary = TRUE)
devtools::install(build_vignettes = TRUE)
devtools::check()


library("SlopeHunter")
utils::packageVersion("SlopeHunter")
data("data_example")
head(data_example)
?data_example
help("data_example")
citation("SlopeHunter")

Sh.Model <- hunt(dat = data_example, xbeta_col="xbeta", xse_col="xse",
                 ybeta_col="ybeta", yse_col="yse", yp_col="yp",
                 xp_thresh = 0.001, Bootstrapping = TRUE, show_adjustments = TRUE, seed=2021)


Sh.Model$b
Adj <- Sh.Model$Fit
head(Adj)
head(Adj$ybeta)

head(Sh.Model$est)

# Adj_sh <- SHadj(Sh.Model, dat = data_example, xbeta_col = "xbeta", xse_col = "xse",
#             ybeta_col = "ybeta", yse_col = "yse")
# head(Adj_sh)

?hunt
require(ggplot2)
require(plotly)
ggplotly(Sh.Model$plot)

# Install it from Github
devtools::install_github("Osmahmoud/SlopeHunter")
require(SlopeHunter)
