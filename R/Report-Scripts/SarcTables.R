#!/usr/bin/env Rscript
library("base")
gb <- globalenv(); assign("gb", gb)
options(width = 320, scipen = 5, knitr.kable.NA = '')
library(verbose=F, warn.conflicts = F, quietly = T, package= "dplyr")


PlotSuppInfo <- function(suppInfo){
    names(suppinfo)<- c("ID","Sentrix ID","Array type","Material type","Gender")
    suppinfo <- as.data.frame(suppinfo)
    colnames(suppinfo) <- " "
    suppinfo <- as.data.frame(t(suppinfo))
    rownames(suppinfo) <- NULL
    xtraCss1="border-radius:0px;border-width:1px;border-style:solid;border-color:rgb(75,0,130);"
    kgb <- c("striped",font_size = 9, bootstrap_options = c("bordered"), position = "float_left")
    suppTab <- suppinfo %>% 
    	knitr::kable("html",c(booktabs = T, escape = F, linesep = ""),align='clc') %>%
    	kableExtra::kable_styling(kgb, full_width = F, position="left") %>% 
        column_spec(column = c(1:5), width = "200px") %>%
        column_spec(1, border_left = "3px solid white") %>%
        column_spec(5, border_right = "3px solid white") %>%
        kableExtra::row_spec(row = 0, font_size = 12,  color = "black") %>%
        kableExtra::row_spec(row = 1, font_size = 11)
    return(suppTab)
}
