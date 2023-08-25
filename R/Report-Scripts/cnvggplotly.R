gb <- globalenv(); assign("gb", gb)

MNPcnv450kNew <- function(Mset, sex = NULL, ...) {
    library(
        verbose = F,
        warn.conflicts = F,
        quietly = T,
        package = "mnp.v11b4"
    )
    path <- file.path(path.package('mnp.v11b4'), 'ext')
    if (is.null(sex)) {
        Rset <- minfi::ratioConvert(Mset, what = "both", keepCN = TRUE)
        sex <-
            ifelse(MNPgetSex(Rset)$predictedSex == "M", "Male", "Female")
    }
    load(file.path(path, "CNanalysis4_conumee_ANNO.vh20150715.RData"))
    cndata <- conumee::CNV.load(Mset)
    if (sex == "Male") {
        load(file.path(path, "CNanalysis4_conumee_REF-M.vh20150715.RData"))
        x <- conumee::CNV.fit(cndata, refM.data, annoXY)
    }
    if (sex == "Female") {
        load(file.path(path, "CNanalysis4_conumee_REF-F.vh20150715.RData"))
        x <- conumee::CNV.fit(cndata, refF.data, annoXY)
    }
    x <- conumee::CNV.bin(x)
    x <- conumee::CNV.detail(x)
    x <- conumee::CNV.segment(x)
    return(x)
}

GetCNxx <- function(Mset, sex, sampleID) {
    is450k <-  Mset@annotation[["array"]] != "IlluminaHumanMethylationEPIC"
    if (is450k == TRUE) {
        load(paste(path.package('mnp.v11b4'),'/ext/ovgenes.RData',sep=''))
        xx <- MNPcnv450kNew(Mset, sex = sex, main = sampleID)
    } else {
        xx <- supM(mnp.v11b6::MNPcnv(Mset, sex = sex, main = sampleID))
    }
    return(xx)
}

GetOvAnnot <- function() {
    cnvUrl <-
        "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/cnvggplotly.R"
    githubURL <-
        "https://github.com/NYU-Molecular-Pathology/Methylation/raw/main/Rdata/newOvGenes.rds"
    gb <- globalenv()
    assign("gb", gb)
    if (!file.exists("newOvGenes.rds")) {
        tryCatch(
            expr = {
                newOvGenes <- readRDS(url(githubURL, method = "libcurl"))
            },
            error = function(e) {
                if (!file.exists("newOvGenes.rds")) {
                    utils::download.file(githubURL, file.path(getwd(), "newOvGenes.rds"), method = "libcurl")
                }
                newOvGenes <- readRDS(file.path(getwd(), "newOvGenes.rds"))
            }
        )
    } else{
        newOvGenes <- readRDS("newOvGenes.rds")
    }
    assign("newOvGenes", newOvGenes, gb)
    return(newOvGenes)
}

GetCNVTables <- function(dra){
    gainDf <- t(data.frame(Gains = c(rownames(dra[(dra$Gain == T),]))))
    lossDf <- t(data.frame(Loss = c(rownames(dra[(dra$Loss == T),]))))
    btso = c("bordered")
    kgb <- c("striped",font_size = 14, bootstrap_options = btso, position = "left")
    kgh <- c(booktabs = T, escape = F, linesep = "")
    extra_css1 = "border-radius:0px;border-width:1px;border-style:solid;border-color:rgb(192,192,192);"
    txtc = "text-align:center;"
    gainTab <- 
    	gainDf %>% knitr::kable("html", kgh, align = 'clc') %>% 
    	kableExtra::kable_styling(kgb, full_width = F, position="left") %>%
    	kableExtra::row_spec(row = 1, extra_css = extra_css1) %>% 
        kableExtra::column_spec(column = 1, background = "palegreen", extra_css = txtc)
    
    lossTab <- 
    	lossDf %>% knitr::kable("html", kgh, align = 'clc') %>% 
    	kableExtra::kable_styling(kgb, full_width = F, position="left") %>%
    	kableExtra::row_spec(row = 1, extra_css = extra_css1) %>% 
        kableExtra::column_spec(column = 1, background = "#CD5C5C", extra_css = txtc)
    return(list("gainDf"=gainDf, "lossDf"=lossDf, "gainTab"=gainTab, "lossTab"=lossTab))
}

gb$newOvGenes <- newOvGenes <- gb$GetOvAnnot()

NewGgplotly <- function (msetDat, dat) {
    Mset<-msetDat$Mset
    sex <- msetDat$sex
    sampleID <- dat$sampleID  
    options(warn = -1)
    newOvGenes <- gb$newOvGenes
    xx <- gb$GetCNxx(Mset, sex, sampleID)
    ylim = c(-2, 2)
    bin.ratio <- xx@bin$ratio - xx@bin$shift
    bin.ratio[bin.ratio < ylim[1]] <- ylim[1]
    bin.ratio[bin.ratio > ylim[2]] <- ylim[2]
    .cumsum0 <- function(x, left = TRUE, right = FALSE, n = NULL) {
        xx <- c(0, cumsum(as.numeric(x)))
        if (!left) {xx <- xx[-1]}
        if (!right) {xx <- head(xx, -1)}
        names(xx) <- n
        xx
    }
    chrX = TRUE
    chrY <- ifelse(sex == "male", T, F)
    chr <- xx@anno@genome$chr
    chr.cumsum0 <- .cumsum0(xx@anno@genome[chr, "size"], n = chr)
    if (!chrX & is.element("chrX", names(chr.cumsum0))) {chr.cumsum0["chrX"] <- NA}
    if (!chrY & is.element("chrY", names(chr.cumsum0))) {chr.cumsum0["chrY"] <- NA}
    x <- chr.cumsum0[as.vector(seqnames(xx@anno@bins))] + (xx@anno@bins)$midpoint
    
    chrs <- .cumsum0(xx@anno@genome[chr, "size"], right = TRUE)
    chrspq <- .cumsum0(xx@anno@genome[chr, "size"]) + xx@anno@genome[chr, "pq"]
    tickl <- .cumsum0(xx@anno@genome[chr, "size"]) + xx@anno@genome[chr, "size"] / 2
    cols2 = c("firebrick4", "firebrick2", "red", "lightpink", "darkgrey", "darkgrey",
              "lightgreen", "green", "green3", "darkgreen")
    colRamOut <- colorRamp(cols2)((bin.ratio + max(abs(ylim))) / (2 * max(abs(ylim))))
    bin.ratio.cols <- apply(colRamOut, 1, function(x){rgb(x[1], x[2], x[3], maxColorValue = 255)})
    df <- data.frame(x, bin.ratio, bin.ratio.cols)
    xs <- xx@seg$summary$loc.start + chr.cumsum0[xx@seg$summary$chrom]
    xe <- xx@seg$summary$loc.end + chr.cumsum0[xx@seg$summary$chrom]
    ys <- xx@seg$summary$seg.median - xx@bin$shift
    ye <- xx@seg$summary$seg.median - xx@bin$shift
    df2 <- data.frame(xs, xe, ys, ye)
    detail.ratio <- xx@detail$ratio - xx@bin$shift
    detail.ratio[detail.ratio < ylim[1]] <- ylim[1]
    detail.ratio[detail.ratio > ylim[2]] <- ylim[2]
    detail.ratio.above <- (detail.ratio > 0.15 & detail.ratio < 1) | detail.ratio < -0.15
    xAnno <- xx@anno@detail
    detail.x <- start( xAnno) + (end( xAnno) - start( xAnno)) / 2 + chr.cumsum0[as.vector(seqnames(xAnno))]
    df3 <- data.frame(detail.ratio, detail.x, names = (xAnno)$name)
    dra <- detail.ratio.above
    et14 <- element_text(size = 14)
    
    p <- 
        ggplot(df, aes(x = x, y = bin.ratio)) + geom_point(colour = bin.ratio.cols, size = 1) +
        geom_vline(xintercept = chrs, color = "grey44", size = 1, alpha = 0.5) +
        geom_vline(xintercept = chrspq, color = "black", size = 0.5, linetype = "dashed", alpha = 0.8) +
        geom_segment(aes(x = xs, y = ys, xend = xe, yend = ye), size = 0.5, data = df2, color = "blue") +
        xlab("Chromosomes") + ylab("Bin Ratios") +
        geom_point(aes(x = detail.x, y = detail.ratio), size = ifelse(test = dra, yes = 2, no = 1.5), 
                   alpha = 0.9, data = df3, color = ifelse(test = dra, yes = "blue", no = "black")) +
        scale_x_continuous(breaks = tickl, labels = c(chr)) +
        scale_y_continuous(breaks = seq(ylim[1], ylim[2], by = 0.25)) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 12, angle = 90), axis.text.y = et14,
              axis.title.y = et14, axis.title.x = et14)
    
    ggpb <- plotly::plotly_build(plotly::ggplotly(p)) 
    
    # Annotate Hover Probes
    if (length(xx@bin$ratio) >= 25666) {
        ggpb$x$data[[1]]$text <- paste0(
            seqnames(xx@anno@bins),
            "<br>", "start: ", start(xx@anno@bins),
            "<br>", "end: ", end(xx@anno@bins),
            "<br>", "probes: ", (xx@anno@bins)$probes,
            "<br>", "Chr: ", xx@anno@bins@ranges@NAMES,
            "<br>", "Genes: ", newOvGenes[1:length(xx@bin$ratio)]
        )
    } else {
        load(paste(path.package('mnp.v11b4'),'/ext/ovgenes.RData',sep=''))
        ggpb$x$data[[1]]$text <- paste0(
            seqnames(xx@anno@bins),
            "<br>", "start: ", start(xx@anno@bins),
            "<br>", "end: ", end(xx@anno@bins),
            "<br>", "probes: ", (xx@anno@bins)$probes,
            "<br>", "Genes: ", ovgenes450k
        )
    }
    ggpb$x$data[[2]]$text <- ""
    ggpb$x$data[[3]]$text <- ""
    ggpb$x$data[[4]]$text <- paste0(
        xx@seg$summary$chrom,
        "<br>", "start: ",xx@seg$summary$loc.start,
        "<br>", "end: ", xx@seg$summary$loc.end, "<br>",
        "median: ", xx@seg$summary$seg.median
    )
    ggpb$x$data[[5]]$text <- (xx@anno@detail)$name
    
    # Custom Detail Annotations
    ggpb <- ggpb %>% 
        add_annotations(
            x = detail.x, y = detail.ratio, text = names(detail.ratio),
            textangle = 60, showarrow = T, arrowwidth = 0.75, xref = "x", yref = "y",
            size = 0.5, arrowhead = 3, ax = -20, ay = -60, arrowsize = 0.5,
            color = "darkgrey", bgcolor = "white", opacity = 0.85) %>%
        ggplotly(tooltip = "text") %>% plotly::style(hoverlabel = list(bgcolor = "blue"))
    
    # Generate Gain/Loss Table
    dra <- data.frame(Gain = c(detail.ratio > 0.15 & detail.ratio < 1.5), Loss = c(detail.ratio < -0.15))
    gainLoss <- gb$GetCNVTables(dra)
    return(list("thePlot"=ggpb %>% toWebGL(),"gainLoss"=gainLoss))
}

MNPciplot_mgmt <- function(Mset, sample = 1) {
    pred <- MNPpredict_mgmt(Mset[, sample])
    theBarPlot <-
        ggplot2::ggplot(pred, aes(x = Estimated, y = Status)) +
        geom_point(size = 5) +
        geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = .1, size = 1) +
        xlim(0, 1) +
        geom_vline(xintercept = pred$Cutoff, size = 1, colour = 'darkred') +
        theme(plot.title = element_text(lineheight = .8, face = "bold")) +
        xlab("Score") + ylab("") + scale_y_discrete("", breaks = c(0))
    theBarPlot <- plotly::ggplotly(theBarPlot)
    theBarPlot <- htmltools::tagList(plotly::ggplotly(theBarPlot)) %>%
        layout(autosize = T, width = 500, height = 200)
    return(theBarPlot)
}

