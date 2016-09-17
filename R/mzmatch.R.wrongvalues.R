"mzmatch.R.wrongvalues" <- function(data){
wrongvalues <- NULL
wrongvalues <- which(is.infinite(data))
wrongvalues <- append(wrongvalues,which(is.nan(data)))
if (!is.null(wrongvalues)){
    data[wrongvalues] <- NA
}
data
}
