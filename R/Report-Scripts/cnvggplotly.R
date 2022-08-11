MNPcnv450kNew <- function(Mset, sex = NULL, ...){
    library(verbose=F,warn.conflicts = F, quietly = T, package="mnp.v11b4")
    path <- file.path(path.package('mnp.v11b4'), 'ext')
    if (is.null(sex)) {
        Rset <- minfi::ratioConvert(Mset, what = "both", keepCN = TRUE)
        sex <- ifelse(MNPgetSex(Rset)$predictedSex == "M", "Male", "Female")
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
    if (is450k == TRUE) {
        xx <- MNPcnv450kNew(Mset, sex = sex, main = sampleID)
    } else {
        xx <- supM(mnp.v11b6::MNPcnv(Mset, sex = sex, main = sampleID))
    }
    return(xx)
}

GetOvAnnot <- function() {
    githubURL <- file.path(
        "https://github.com/NYU-Molecular-Pathology/Methylation/raw",
        "416a007b8a21f59a71493cea189bc424009e8d7d/Rdata/newOvGenes.rds"
    )
    if (!file.exists("newOvGenes.rds")) {
        utils::download.file(githubURL, file.path(getwd(), "newOvGenes.rds"), method = "libcurl")
    }
    newOvGenes <- readRDS("newOvGenes.rds")
    return(newOvGenes)
}


new.ggplotly <- function (xx, getTables = T, newOvGenes=NULL, sex='male', addCustom = F) {
    require(compiler)
    compiler::enableJIT(3)
    compiler::compilePKGS(enable = TRUE)
    compiler::setCompilerOptions(suppressAll = TRUE, optimize = 3)
    ovDataPath <- paste(path.package("mnp.v11b6"), "/ext/ovgenes.RData", sep = "")
    load(ovDataPath)
    # permalink to annotations
    #readRDS("newOvGenes.rds")
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
    chrY <- ifelse(sex=="male", T,F)
    chr <- xx@anno@genome$chr
    chr.cumsum0 <- .cumsum0(xx@anno@genome[chr, "size"], n = chr)
    if (!chrX & is.element("chrX", names(chr.cumsum0))) {chr.cumsum0["chrX"] <- NA}
    if (!chrY & is.element("chrY", names(chr.cumsum0))) {chr.cumsum0["chrY"] <- NA}
    x <- chr.cumsum0[as.vector(seqnames(xx@anno@bins))] + (xx@anno@bins)$midpoint
    chrs <- .cumsum0(xx@anno@genome[chr, "size"], right = TRUE)
    chrspq <- .cumsum0(xx@anno@genome[chr, "size"]) + xx@anno@genome[chr, "pq"]
    tickl <- .cumsum0(xx@anno@genome[chr, "size"]) + xx@anno@genome[chr, "size"]/2
    cols2 = c("firebrick4", "firebrick2", "red", "lightpink", "darkgrey", "darkgrey", "lightgreen", "green", "green3",
              "darkgreen")
    bin.ratio.cols <- apply(colorRamp(cols2)(
        (bin.ratio + max(abs(ylim)))/(2 * max(abs(ylim)))), 1,
        function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
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
    detail.x <-
        start(xx@anno@detail) + (end(xx@anno@detail) - start(xx@anno@detail))/2 + chr.cumsum0[as.vector(seqnames(xx@anno@detail))]

    chromocols <- c("#2D77E6", "#785C0B", "#6D26A5", "#4E4E4E",
                    "#19774E", "#14E72D", "#2D77E6", "#C991E5", "#19AC47",
                    "#B23E19", "#C991E5", "#E600E3", "#2D518C", "#E69F14",
                    "#E19290", "#829D19", "#DF1E29", "#19E7BB", "#29C5E7",
                    "#A00B92", "#AE6A96", "#E519AD", "#8A4694", "#E77300")
    df3 <- data.frame(detail.ratio, detail.x, names = (xx@anno@detail)$name)
    dra <- detail.ratio.above
    p <- ggplot(df, aes(x = x, y = bin.ratio)) +
        geom_point(colour = bin.ratio.cols, size = 1) +
        geom_vline(xintercept = chrs, color = "grey44",
                   size = 1, alpha = 0.5) +
        geom_vline(xintercept = chrspq, color = "black",
                   size = 0.5, linetype = "dashed", alpha = 0.8) +
        ylim(ylim[1], ylim[2]) +
        geom_segment(aes(x = xs, y = ys, xend = xe, yend = ye),
                     size = 0.5, data = df2, color = "blue") +
        xlab("Chromosomes") +
        ylab("Bin Ratios") +
        geom_point(aes(x = detail.x, y = detail.ratio),
                   size = ifelse(test = dra, yes = 2, no = 1.5),
                   alpha = 0.9, data = df3,
                   color = ifelse(test = dra, yes = "blue", no = "black")) +
        scale_x_continuous(breaks = tickl, labels = c(chr)) +
        scale_y_continuous(breaks = seq(ylim[1], ylim[2], by = 0.25)) +
        theme(axis.text.x = element_text(size = 12, angle = 90),
              axis.text.y = element_text(size = 14)) +
        theme_light() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 12, angle = 90),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.title.x = element_text(size = 14))
    ggp <- suppressMessages(suppressWarnings(plotly::ggplotly(p)))
    ggpb <- suppressMessages(suppressWarnings(plotly::plotly_build(ggp)))
    ggpb <- ggpb %>% plotly::layout(xaxis = list(showgrid = F), yaxis = list(showgrid = F))

    # Annotate Hover Probes
    if (length(xx@bin$ratio) >= 25666) {
        ggpb$x$data[[1]]$text <- paste0(
            seqnames(xx@anno@bins), "<br>",
            "start: ", start(xx@anno@bins), "<br>",
            "end: ", end(xx@anno@bins), "<br>",
            "probes: ", (xx@anno@bins)$probes, "<br>",
            "Chr: ", xx@anno@bins@ranges@NAMES,"<br>",
            "Genes: ", newOvGenes[1:length(xx@bin$ratio)])
    }
    else {
        ggpb$x$data[[1]]$text <- paste0(
            seqnames(xx@anno@bins),
            "<br>", "start: ", start(xx@anno@bins), "<br>",
            "end: ", end(xx@anno@bins), "<br>", "probes: ",
            (xx@anno@bins)$probes, "<br>",
            "Genes: ", ovgenes450k)
    }
    ggpb$x$data[[2]]$text <- ""
    ggpb$x$data[[3]]$text <- ""
    ggpb$x$data[[4]]$text <- paste0(
        xx@seg$summary$chrom, "<br>",
        "start: ", xx@seg$summary$loc.start, "<br>",
        "end: ", xx@seg$summary$loc.end, "<br>",
        "median: ", xx@seg$summary$seg.median)

    ggpb$x$data[[5]]$text <- (xx@anno@detail)$name

    # Custom Detail Annotations
    ggpb <- ggpb %>% add_annotations(
        x = detail.x, y = detail.ratio,
        text = names(detail.ratio), textangle = 60, showarrow = T,
        arrowwidth = 0.75, xref = "x", yref = "y", size = 0.5,
        arrowhead = 3, ax = -20, ay = -60, arrowsize = 0.5,
        color = "darkgrey", font = list(size = 15, color = "blue",face = "bold"),
        bgcolor = "white", opacity = 0.85)
    # Return Plot or Table
    if (!getTables) {
        return(suppressMessages(suppressWarnings(ggpb %>% toWebGL())))
    }
    if (getTables) {
        xyRatio <- xx@bin[["ratio"]]
        #         newDataAnno <- data.frame(
        #             newAnnotation =newOvGenes[1:length( xx@anno@bins$probes)],
        #             start = start(xx@anno@bins),
        #             end =  end(xx@anno@bins),
        #             chromosome = xx@anno@bins@ranges@NAMES,
        #             probes = xx@anno@bins$probes,
        #             yvalues = xyRatio,
        #             gains = ifelse(xyRatio>=0.25,"GAIN",""),
        #             loss = ifelse(xyRatio<=-0.25,"LOSS","")
        #         )
        #         finame <- paste0(xx@name,"_newDataAnno.csv")
        #         write.csv(newDataAnno,finame)
        dra <- data.frame(
            Gain = c(detail.ratio > 0.15 & detail.ratio <1.5),
            Loss = c(detail.ratio < -0.15)
        )
        return(dra)
    }
}