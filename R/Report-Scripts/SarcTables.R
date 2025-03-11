#!/usr/bin/env Rscript
## Script name: SarcTables.R
## Purpose: Source Global scripts for Sarcoma analysis report tables
## Date Created: January 19, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

gb <- globalenv()
assign("gb", gb)
options(width = 320, scipen = 5, knitr.kable.NA = '')

PredTable <- function(predLi){
    kgb <- c("striped", font_size = 9, bootstrap_options = c("bordered"))
    pt <- predLi$pred %>% 
        knitr::kable("html",c(booktabs = T, escape = F, linesep = ""), align = 'clc') %>% 
        kableExtra::kable_styling(kgb, full_width = F, position = "left") %>% 
        kableExtra::row_spec(row = 0, font_size = 12, background = "rgb(220, 208, 255)", color = "black") %>%
        kableExtra::row_spec(row = 1, font_size = 11)
    return(pt)
}


PredTopTable <- function(predLi){
    kgb <- c("striped", font_size = 9, bootstrap_options = c("bordered"))
    pt <- predLi$predTop %>% 
        knitr::kable("html",c(booktabs = T, escape = F, linesep = ""), align = 'clc') %>% 
        kableExtra::kable_styling(kgb, full_width = F, position = "left") %>% 
        kableExtra::row_spec(row = 0, font_size = 12, background = "rgb(211, 211, 211)", color = "black") %>%
        kableExtra::row_spec(row = c(1:3), font_size = 11)
    return(pt)
}
