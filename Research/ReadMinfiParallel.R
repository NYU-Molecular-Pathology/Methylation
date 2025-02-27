GuessArrayTypes <- function (nProbes)
{
    if (nProbes >= 622000 && nProbes <= 623000) {
        arrayAnnotation <- c(array = "IlluminaHumanMethylation450k",
                             annotation = minfi:::.default.450k.annotation)
    }
    else if (nProbes >= 1050000 && nProbes <= 1053000) {
        arrayAnnotation <- c(array = "IlluminaHumanMethylationEPIC",
                             annotation = minfi:::.default.epic.annotation)
    }
    else if (nProbes >= 1032000 && nProbes <= 1033000) {
        arrayAnnotation <- c(array = "IlluminaHumanMethylationEPIC",
                             annotation = minfi:::.default.epic.annotation)
    }
    else if (nProbes >= 54000 && nProbes <= 56000) {
        arrayAnnotation <- c(array = "IlluminaHumanMethylation27k",
                             annotation = minfi:::.default.27k.annotation)
    }
    else if (nProbes >= 41000 & nProbes <= 41100) {
        arrayAnnotation <- c(array = "HorvathMammalMethylChip40",
                             annotation = "test.unknown")
    }
    else if (nProbes >= 43650 & nProbes <= 43680) {
        arrayAnnotation <- c(array = "IlluminaHumanMethylationAllergy",
                             annotation = minfi:::.default.allergy.annotation)
    }
    else {
        arrayAnnotation <- c(array = "Unknown", annotation = "Unknown")
        warning("Unable to detect Array type")
    }
    arrayAnnotation
}


ReadRGparallel <- function (targets, folder, files = NULL, copy = FALSE, verbose = TRUE,
          arraytype = NULL, ncores = NULL, extended = FALSE, force = TRUE) {
    require("minfi");require("parallel");require("doParallel")
    if(class(targets)[1]!="data.table"){targets <- as.data.table(targets)}
    if (is.null(files)) {files <- targets$Basename}
    if (copy == TRUE) {files <- paste0(folder, basename(targets$Basename))}
    if(is.null(ncores)){ncores <- parallel::detectCores()/2}
    cl <- parallel::makePSOCKcluster(ncores)
    parallel::clusterEvalQ(cl, {requireNamespace(c("minfi", "S4Vectors"))})
    doParallel::registerDoParallel(cl)
    message("Reading multiple idat-files in parallel. Using ", ncores, " cores.")
    res <- foreach::foreach(
        it = itertools::isplitIndices(nrow(targets), chunks = ncores),
        .combine = "cbind",
        .multicombine = F,
        .export = c("GuessArrayTypes", "minfi:::.default.epic.annotation"),
        .inorder = F,
        .errorhandling = "pass"
    ) %dopar% {
        subdf <- as.data.frame(targets[it,])
        if (copy == TRUE) {
            requireNamespace("fs")
            fs::file_copy(
                paste0(subdf$Basename, "_Grn.idat"),
                new_path = folder,
                overwrite = T
            )
            fs::file_copy(
                paste0(subdf$Basename, "_Red.idat"),
                new_path = folder,
                overwrite = T
            )
        }
        rgSet <-
            minfi::read.metharray(
                basenames = files[it],
                extended = extended,
                verbose = verbose,
                force = force
            )
        if (is.null(arraytype)) {
            rgSet@annotation <- GuessArrayTypes(nrow(rgSet))
        }
        else {
            if (arraytype == "EPIC") {
                rgSet@annotation <- c(array = "IlluminaHumanMethylationEPIC",
                                      annotation = "ilm10b4.hg19")
            }
            else if (arraytype == "450K") {
                rgSet@annotation <- c(array = "IlluminaHumanMethylation450k",
                                      annotation = "ilmn12.hg19")
            }
            else {
                rgSet@annotation <- GuessArrayTypes(nrow(rgSet))
            }
        }
        return(rgSet)
    }
    parallel::stopCluster(cl)
    cn <- colnames(res)
    class(res)
    data.table::setkey(targets, "barcode")
    pD <- data.frame(targets[cn,])
    pD$filenames <- files
    res@colData <- methods::as(pD, "DataFrame")
    rownames(res@colData) <- cn
    colnames(res) <- cn
    return(res)
}
