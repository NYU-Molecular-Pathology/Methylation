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
        autoWidth = T,
        pageLength = 15,
        rownames = F,
        lengthChange = T,
        searchable = T,
        columnDefs = list(list(width = '10%', visible=TRUE, targets = c(1:ncol(targets))))
      )
    )
  return(theDt)
}

smallTab <- function(dtObj) {
  dtTable <- knitr::kable(dtObj, row.names = F, "html")
  dtTable <-
    kableExtra::kable_styling(
      dtTable,
      bootstrap_options = c("striped", "condensed"),
      full_width = F,
      position = "left"
    )
  dtTable <-
    kableExtra::column_spec(dtTable, 1:ncol(dtObj), width = "4cm")
  return(dtTable)
}
