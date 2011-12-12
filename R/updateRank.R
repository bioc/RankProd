## fb: This function does most of the hard work - computing matrix ranks iteratively column by column
updateRank <- function( oldRank, col1, nextCol ) {

  ord <- order(col1);
  col1<-col1[ord];    # move NAs to the end
  oldRank<-oldRank[ord];

  nextCol<-sort(nextCol); # remove NAs here

  i <- length(col1);
  j <- length(nextCol);

  while (is.na(col1[i])) i<-i-1;

  while(i>0 && j>0) {
    while(j>0 && nextCol[j]>=col1[i]) j<-j-1;
    k<-i;
    while(i>0 && j>0 && nextCol[j]<col1[i])  i<-i-1;
    oldRank[(i+1):k] <- oldRank[(i+1):k] + j;
  }

  idx <- is.na(col1);
  oldRank[idx] <- oldRank[idx] + length(nextCol);

  return(oldRank[order(ord)])
}
