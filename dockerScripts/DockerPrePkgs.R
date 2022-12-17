update.packages(
    ask = FALSE,
    contriburl = contrib.url(
        'https://packagemanager.rstudio.com/all/__linux__/focal/latest'
    )
)

if (!require('rJava')) {
    install.packages(lib = '/usr/local/lib/R/site-library/', 'rJava')
}

install.packages(lib = '/usr/local/lib/R/site-library/', 'pak',
                 repos = 'https://r-lib.github.io/p/pak/dev/', dependencies = T)

options(repos = 'https://packagemanager.rstudio.com/all/__linux__/focal/latest')
pak::pkg_install(
    c(
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
    ), dependencies = T, ask = F 
)
pak::pkg_install(
    c('devtools', 'BiocVersion', 'BiocManager'),
    dependencies = T,
    ask = F,
    upgrade = T
)
pak::pkg_install(
    c('xtable', 'prettyunits', 'munsell', 'ggplot2', 'remotes'), dependencies = T, ask = F
)
remotes::install_github("SymbolixAU/googlePolylines", dependencies = T)

install.packages(lib = '/usr/local/lib/R/site-library/',
                 c('nleqslv', 'ade4', 'pixmap', 'stringi'),
                 dependencies = T)

