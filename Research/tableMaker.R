gb <- globalenv(); assign("gb", gb)
setDirectory<-function(foldr) {
    bsDir = paste("cd", foldr); mm2 = crayon::white$bgRed("Location Not Found:", foldr)
    if (dir.exists(foldr)) {system(bsDir);setwd(foldr);assign("runDir", foldr)} else{warning(mm2)}
}

makeDt <- function(targets) {
  theDt <-
    DT::datatable(
      targets,
      selection = "single",
      filter = list(
        position = 'top',
        clear = F,
        plain = T
      ),
      autoHideNavigation = T,
      options = list(
        scrollX = T,
        info = T,
        autoWidth = F,
        pageLength = 15,
        rownames = F,
        lengthChange = T,
        searchable = T,
        columnDefs = list(list(width = '10%', visible=TRUE, targets = c(1:ncol(targets))))
      )
    )
  return(theDt)
}

