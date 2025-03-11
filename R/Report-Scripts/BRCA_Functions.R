#!/usr/bin/env Rscript
## Script name: BRCA_Functions.R
## Purpose: Functions for methylation report.Rmd BRCA1 and BRCA2 loci status
## Date Created: January 19, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

gb <- globalenv(); assign("gb", gb)
# DEBUG & TEST ---------------------------------------------------------
# if(!require("minfiDataEPIC")){BiocManager::install("minfiDataEPIC")}
# RGset <- minfiDataEPIC::RGsetEPIC
# RGset <- RGset[,1]

GetBRCA_data <- function(ratioSet, RGset, brcaProbes){
    #Output is DF with sample information, betas for relevant BRCA loci
    brcaBetas = minfi::getBeta(ratioSet)[brcaProbes,]
    brcaBetas <- as.data.frame(brcaBetas)
    for (bal in 1:nrow(brcaBetas)) {
        if (is.nan(brcaBetas[bal, 1])) {brcaBetas[bal, 1] = 0.00}
        brcaBetas[bal, 1] <- round(brcaBetas[bal, 1], digits = 3)
        if (brcaBetas[bal, 1] < 0.09) {
            brcaBetas[bal, 1] = 0.00
        }
        if (brcaBetas[bal, 1] > 1.00) {
            brcaBetas[bal, 1] = 1.00
        } else {
            (brcaBetas[bal, 1] = brcaBetas[bal, 1])
        }
        }
    #Calculate detection of Pvalues using minfi
    detection.pvals = minfi::detectionP(RGset)
    l.pvals = detection.pvals[c(brcaProbes),]
    l.pvals <- as.data.frame(l.pvals)
    for (lpv in 1:nrow(l.pvals)) {
        if (is.nan(l.pvals[lpv,1])) {l.pvals[lpv,1] = 0.00}
        l.pvals[lpv,1] <- round(l.pvals[lpv,1], digits = 3)
        if (l.pvals[lpv, 1] < 0.09) {
            l.pvals[lpv, 1] = 0.00
        }
        if (l.pvals[lpv, 1] > 1.00) {
            l.pvals[lpv, 1] = 1.00
        } else {
            l.pvals[lpv, 1] = l.pvals[lpv, 1]
        }
    }
    #calculate number of methylated loci with beta values above a cutoff
    pos.loci = brcaBetas > 0.50
    m.reslt = rep("", nrow(brcaBetas))
    m.reslt[pos.loci == 5] = "POSITIVE"
    m.reslt[pos.loci >= 0 & pos.loci < 4] = "NEGATIVE"
    m.reslt[pos.loci == 4] = "INDETERMINATE"

    clin.res = data.frame(
        "Sample_Name" = 1,
        "betaValues" = brcaBetas[,1],
        "BRCA.pos.loci" = pos.loci[,1],
        "m.reslt" = m.reslt,
        "Pvalue" = l.pvals,
        stringsAsFactors = F
    )
    return(clin.res)
}


sanitizeDense <- function(beta.matrix){
    beta.matrix <- as.data.frame(beta.matrix)
    beta.matrix[is.na(beta.matrix)] <- 0
    beta.matrix[beta.matrix < 0] <- 0
    beta.matrix[beta.matrix > 1] <- 1
    beta.matrix <- round(beta.matrix, digits = 3)
    density <- density(beta.matrix[, 1])
    density_df <- data.frame(x = density$x, y = density$y)
    return(density_df)
}


renderBRCA_plot <- function(beta.matrix, clin.res, brca = "BRCA1") {
    clinVals <- as.data.frame(t(clin.res[2:5]))
    col_vect <- c(
        '#e6194b',
        '#3cb44b',
        '#ffe119',
        '#4363d8',
        '#f58231',
        '#911eb4',
        '#46f0f0',
        '#f032e6',
        '#bcf60c',
        '#fabebe'
    )
    col_vect <- col_vect[1:nrow(clin.res)]
    Probe.Names = paste(colnames(clinVals))
    names(Probe.Names) <- col_vect
    df <-
        data.frame(
            x = as.numeric(paste(clinVals[1, ])),
            y = seq(0.1, 2, length.out = nrow(clin.res)),
            Probe.Names = Probe.Names
        )
    density <- sanitizeDense(beta.matrix)
    fig <- ggplot(data = df, aes(x, y)) +
        geom_bar(
            aes(fill = Probe.Names),
            stat = 'identity',
            width = 0.05,
            color = "black"
        ) +
        scale_fill_manual("legend", values = names(Probe.Names)) +
        geom_line(density,
                  mapping = aes(x = x, y = y),
                  color = "black") +
        geom_vline(xintercept = 0.5,
                   color = "red",
                   linetype = 2) +
        xlab('Beta Values') + ylab('Density') + theme_bw() +
        guides(fill = guide_legend(title = paste(brca, "Probes"))) +
        labs(
            caption = "The black line shows overall beta value density for the sample between 0 to 1.
Beta values close to zero < 0.5 indicate low methylation for that probe loci. Betas closer to 1.0 > 0.5
indicate hypermethylation. Color bars show probe beta values on the x-axis.  Bar height is irrelevant."
        ) +
        theme(
            plot.caption = element_text(
                size = 9,
                face = "italic",
                color = "blue",
                hjust = 0
            ),
            plot.caption.position = 'plot'
        )
    return(suppressWarnings(fig))
}


drawPlotTab <- function(brcaTab){
    btso = c("bordered")
    kgb <- c("striped",font_size = 12, bootstrap_options = btso, position = "left")
    txtc = "text-align:center;"
    be = c(booktabs = T, escape = F, linesep = "")
    theBRTab <-
        brcaTab %>% knitr::kable("html", be, align = 'clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position = "left") %>%
        kableExtra::column_spec(column = 2, background = "rgb(204, 255, 204)", extra_css = txtc) %>%
        kableExtra::column_spec(column = 2, bold = T, border_right = T)
    return(theBRTab)
}


getBRCA_table <- function(clin.res, brca) {
    title.string = c(
        "BRCA Total Probe Loci" = paste(nrow(clin.res)),
        "Total Methylated Loci" = length(which(clin.res$BRCA.pos.loci))
    )
    title.string <- as.data.frame(title.string)
    colnames(title.string) = NULL
    return(title.string)
}


BRCAstatus <- function(RGset, brca = "BRCA1") {
    brca1Probes <- c(
        "cg22435081",
        "cg20130344",
        "cg26081003",
        "cg24900425",
        "cg15109179",
        "cg10527302",
        "cg14780845",
        "cg05368731",
        "cg26370022",
        "cg11529738"
    )

    brca2Probes <- c(
        "cg09698981",
        "cg19634711",
        "cg12836863",
        "cg05153735",
        "cg26458617",
        "cg24674189",
        "cg15856230"
    )
    if (brca == "BRCA1") {
        brcaProbes <- brca1Probes
    } else{
        brcaProbes <- brca2Probes
    }
    Mset = minfi::preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls")
    ratioSet = ratioConvert(Mset, what = "both", keepCN = FALSE)
    clin.res <- GetBRCA_data(ratioSet, RGset, brcaProbes)

    #old_par <- par()
    par(mar = c(5, 6, 4, 1) + .1)

    beta.matrix = minfi::getBeta(ratioSet)
    brcaPlot1 <- renderBRCA_plot(beta.matrix, clin.res, brca)
    brcaTab <- getBRCA_table(clin.res, brca)

    return(list(
        "brcaPlot1" = brcaPlot1,
        "brcaTab" = drawPlotTab(brcaTab),
        "theValue" = clin.res
    ))
}
