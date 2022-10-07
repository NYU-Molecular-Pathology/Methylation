#!/usr/bin/env Rscript
library("base")
gb <- globalenv(); assign("gb", gb)
options(width = 320, scipen = 5, knitr.kable.NA = '')
library(verbose=F, warn.conflicts = F, quietly = T, package= "dplyr")
library(verbose=F, warn.conflicts = F, quietly = T, package= "sarc.v12b6")
library(verbose=F, warn.conflicts = F, quietly = T, package= "htmltools")

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
        kableExtra::row_spec(row = 0, font_size = 12, background="rgb(255, 235, 205)", color = "black") %>%
      kableExtra::row_spec(row = 0, background="rgb(255, 235, 205)") %>%
        kableExtra::row_spec(row = 1, font_size = 11)
    return(suppTab)
}

GetSarcPred <- function(predRaw){
    pred <-predRaw[,1:2]
    pred$abbrevation <- sarc.v12b6::reflist[match(pred$predicted, sarc.v12b6::reflist$internal_identifier),"abbreviation"]
    pred$maxscore <- round(pmax(pmin(pred$maxscore,1-1e-4),1e-4),4)
    pred$name <- sarc.v12b6::reflist[match(pred$predicted, sarc.v12b6::reflist$internal_identifier),"name"]
    pred <- as.data.frame(pred)
    pred$maxscore <- as.character(pred$maxscore)
    rownames(pred) <- NULL
    colnames(pred) <- c("Predicted", "Max Score", "Abbrevation", "Class Name")
    
    predNex <- predRaw[1,3:ncol(predRaw)]
    predTop <- as.data.frame(t(predNex))
    predTop$Class <- rownames(predTop)
    oo <- order(predTop[,1], decreasing = T)
    predTop <- predTop[oo,1:2]
    predTop <- predTop[2:6,1:2]
    predTop[1:3,1] <- as.character(round(pmax(pmin(predTop[1:3,1],1-1e-4),1e-4),4))
    colnames(predTop) <- c("Other Top Scores", "Class")
    predTop$Name <- sarc.v12b6::reflist[match(predTop$Class, sarc.v12b6::reflist$internal_identifier),"name"]
    rownames(predTop) <- NULL
    des <- sarc.v12b6::reflist[match(pred$Predicted, sarc.v12b6::reflist$internal_identifier),"description"]
    return(list("des"=des, "pred"=pred, "predTop"=predTop))
}

PredTable <- function(predLi){
    kgb <- c("striped", font_size = 9, bootstrap_options = c("bordered"))
    pt <- predLi$pred %>% 
        knitr::kable("html",c(booktabs = T, escape = F, linesep = ""),align='clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position="float_left") %>%
        kableExtra::row_spec(row = 0, font_size = 12, background = "rgb(220, 208, 255)", color = "black") %>%
        kableExtra::row_spec(row = 1, font_size = 11)
    return(pt)
}

PredTopTable <- function(predLi){
    kgb <- c("striped", font_size = 9, bootstrap_options = c("bordered"))
    pt <- predLi$predTop %>% 
        knitr::kable("html",c(booktabs = T, escape = F, linesep = ""),align='clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position="left") %>%
        kableExtra::row_spec(row = 0, font_size = 12, background = "rgb(211, 211, 211)", color = "black") %>%
        kableExtra::row_spec(row = c(1:3), font_size = 11)
    return(pt)
}

