progress_bar <- function(n){
  cat("Working: ")
  for(i in 1:n) {
    signal = floor(n/200)
    if(!(i %% signal)) {
      prog = floor(i/n*50)
      cat('\r', "Working  [", strrep("#",prog),
          strrep(" ", 50-prog), "] ", 2*prog, "%", sep="")
    }
    ## Pretend that some work is being done.
    Sys.sleep(0.001)
    if(i == n) { cat(" \n\r") }
  }
}
