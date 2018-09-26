#' Find the path from root to highest occurrence of each element
#' 
#' A protoclust object has a prototype associated with each interior node.
#' Every element being clustered occurs at least as a leaf but might also 
#' appear multiple times as a prototype.  This function finds for each element
#' the path from the root to the highest occurrence of that element.  The path
#' is specified by a series of 0s and 1s, where 0 means "go left" and 1 means
#' "go right".
#' 
#' @param hc a protoclust object
#' 
#' @return 
#' \item{paths}{a list of length n giving, for each element, the path from
#' root to its highest occurrence. A 0 means go left, a 1 means go right.}
#' \item{int_paths}{a list of length n - 1 giving, for each interior node,
#' the path from root to it. A 0 means go left, a 1 means go right.}
#' @export
find_elements <- function(hc) {
  n <- length(hc$order)
  paths <- vector("list", n) 
  # paths[[i]] will have path from root to element i
  int_paths <- vector("list", n - 1) # int_paths[[n-1]] is root
  # int_paths[[i]] will have path from root to interior node i
  for (i in seq(n - 1, 1)) {
    ii <- hc$merge[i, ] # i-th row of merge
    for (j in 1:2) {
      new_path <- c(int_paths[[i]], j - 1)
      if (ii[j] < 0) {
        if (is.null(paths[[-ii[j]]]))
          paths[[-ii[j]]] <- new_path
      } else {
        int_paths[[ii[j]]] <- new_path
        if (is.null(paths[[hc$protos[ii[j]]]]))
          paths[[hc$protos[ii[j]]]] <- new_path
      }
    }
  }
  paths[[hc$protos[n-1]]] <- NA # this is the root
  list(paths = paths, int_paths = int_paths)
}