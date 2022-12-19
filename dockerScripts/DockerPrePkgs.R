focalPk = 'https://packagemanager.rstudio.com/all/__linux__/focal/latest'
locaLib = '/usr/local/lib/R/site-library/'
basePkg <- c(
    'shiny',
    'rmarkdown',
    'Hmisc',
    'rjson',
    'caret',
    'DBI',
    'RPostgres',
    'curl',
    'httr',
    'xml2'
)
bioCpre <- c('devtools', 'BiocVersion', 'BiocManager')
baseTab <- c('xtable', 'prettyunits', 'munsell', 'ggplot2', 'remotes')
baseAde <- c('nleqslv', 'ade4', 'pixmap', "stringi")

update.packages(contriburl = contrib.url(focalPk), ask = F)
if (!require('rJava')) {install.packages(lib = locaLib, 'rJava', dependencies=T, ask = F)}
install.packages(lib = locaLib, 'pak',
                 repos = 'https://r-lib.github.io/p/pak/dev/', dependencies = T)

pak::pkg_install(basePkg, dependencies = T, ask = F, upgrade = F)
pak::pkg_install(bioCpre, dependencies = T, ask = F, upgrade = T)
pak::pkg_install(baseTab, dependencies = T, ask = F, upgrade = F)
remotes::install_github("SymbolixAU/googlePolylines", dependencies = T)
install.packages(lib = locaLib, baseAde, dependencies = T, ask = F)
devtools::install_github("mdsumner/ncdf4", dependencies = T, upgrade ="always")
if(!require("arrow")){install.packages("arrow", dependencies = T, verbose = T, ask = F)}

