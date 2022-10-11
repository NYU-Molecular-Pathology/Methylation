getScores <- function(Mset=NULL){
    if(!is.null(Mset)){
        betas <- minfi::getBeta(Mset)
        ex <- which(!names(mnp.v11b6::refset.center) %in% rownames(Mset))
        if (length(ex) > 0)
            refset.center <- mnp.v11b6::refset.center[-ex]
        betas <- betas[match(names(refset.center), rownames(betas)), ]
        betas <- betas - refset.center
        betas <- t(as.matrix(betas, ncol = 1))
        if (length(ex) > 0)
            pcaloadings <- mnp.v11b6::pcaloadings[-ex, ]
        sc <- betas %*% pcaloadings
        scores <- rbind(mnp.v11b6::pcascores, sc)
        return(scores)#na.exclude
    }
}

MNPtsne2 <- function(scores = NULL) {
    if (!is.null(scores)) {
        Rtsne::Rtsne(
            scores,
            dims = 3,
            pca = F,
            max_iter = 1500,
            theta = 0,
            perplexity = 30,
            is_distance = F
        )
    }
}
set.seed(seed = 12345)
getScore_cache = compiler::cmpfun(getScores)
tsne_cache = compiler::cmpfun(MNPtsne2) #tsne_cache <- local(MNPtsne.cmpfun())


GetClusterPlot <- function(msetDat, dat){
    sampleID <- dat$sampleID
    Mset_ba <- msetDat$Mset_ba
    is450k <- Mset_ba@annotation[["array"]]=="IlluminaHumanMethylationEPIC"
    me.col <- c(
        "#1F78B4","#2980BA","#3388C1","#3D91C8","#4799CF","hotpink","hotpink","hotpink","hotpink",
        "hotpink","hotpink","#58D68D","#8FCDC0","#96D1B3","#9ED5A7","#A6D99B","#AEDD8F","#A9DB82",
        "#97D574","#85CF65","#74C957","#62C349","#50BD3A","#3EB72C","#2DB11D","#26AC17","#41AA27",
        "#5BA737","#76A547","#ff4d00","#ABA068","#ff4d00","#E09C88","#FA9998","#ff4d00","#F48D78",
        "#ff4d00","#ED8157","#EA7B47","#E77536","#E46F26","#E06916","#E26F1A","#E67A25","#E98531",
        "#ED8F3C","#F09A47","#F4A452","#F7AF5E","#FBB969","#FDB368","#FD9C5A","#FD854D","#FD6D3F",
        "#FE5632","#FE3F24","#FE2817","#FE110A","#FC0406","#F11421","#E6253B","#DC3556","#D14670",
        "#C6568B","#BC67A5","#B177C0","#A687DB","#A17DDB","#9C72DC","#9668DC","#915DDD","#8B53DE",
        "#8648DE","#803EDF","#7B33DF","#863FCB","#9653B0","#A66795","#B67B7A","#C68F5F","#D6A343",
        "#E6B728","#F6CB0D","#FFFF00","#FFFF00","#FFFF00","#F6D916","#5e2f1b","#2980B9","#b8b8b8",
        "#A9A9A9","#A9A9A9","#A9A9A9","#A9A9A9","#9E9E9E","#939393","#868686","#787878","#676767",
        "#0400ff"
        )
    set.seed(12345)

    res<-NULL
    if(is450k==TRUE){
        res=mnp.v11b4::MNPtsne(Mset_ba)
    } else {
        scores <- getScore_cache(Mset_ba)
        res <- tsne_cache(scores)
    }
    # Plot point Values ---------------------------------------
    Tsne_1 <- res$Y[,1]
    Tsne_2 <- res$Y[,2]
    plTxt <- paste0("<b>", sampleID, "</b>")
    t1 = Tsne_1[2802]
    t2 = Tsne_2[2802]
    # Arrow annotation for sample ---------------------------------------
    a <-
        list(
            x = t1,
            y = t2,
            text = plTxt,
            xref = "x",
            yref = "y",
            size = 0.6,
            bgcolor = 'rgba(255, 255, 255, 0.5)',
            font = list(size = 13, color = "dimgray"),
            showarrow = T,
            arrowhead = 0,
            ax = 50,
            ay = -50,
            arrowsize = 1,
            arrowcolor = "orange"
        )
    # Line properties ---------------------------------------
    lnd <- list(color = "red", alpha = 0.8, width = 0.75,dash = 'dash')
    vline<-function(x=0){list(type="line",y0=0,y1=1,yref="paper",x0=x,x1=x,line=lnd)}
    hline<-function(y=0){list(type="line",x0=0,x1=1,xref="paper",y0=y,y1=y,line=lnd)}
    # Font parameters ---------------------------------------
    f1 <- list(family = "Arial,sans-serif", size = 18, color = "darkgrey")
    f2 <- list(size = 14, face = "bold")
    d1 <- list(titlefont=f1, showticklabels=T, tickfont=f1,
               zerolinecolor=plotly::toRGB("grey"), zerolinewidth=0.5)
    b1 <- c(list(title = "T-sne 1"), d1)
    b2 <- c(list(title = "T-sne 2"), d1)
    m1 <- list(size = 12, line = list(width = 1.7, color = 'DarkSlateGrey'))
    # Generate Plotly ---------------------------------------
    theplot <- plotly::plot_ly(
        x = ~ Tsne_1,
        y = ~ Tsne_2,
        color = ~ methclass,
        colors = me.col,
        opacity = 0.9,
        marker = list(size = 5),
        font = f2
    )
    # Modify Plotly ---------------------------------------
    theplot <- theplot %>%
        plotly::layout(type = "scatter", shapes = list(vline(t1), hline(t2))) %>%
        plotly::layout(legend = list(orientation = "h", xanchor = "center")) %>%
        plotly::layout(xaxis = b1, yaxis = b2, showlegend = T) %>%
        plotly::layout(type = "scatter", annotations = a, size = 1.5)
    theplot <- theplot %>%
        plotly::add_trace(theplot, mode = 'markers', type = "scatter", marker = m1)
    theplot <- theplot %>%
        plotly::layout(legend = list(orientation = "v", xanchor = "right", x = 1.1))
    return(theplot)
}
