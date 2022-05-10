#' Convert irregular list into two column dataframe
#'
#' This function will convert an irregular list into a two column dataframe, where the first column is the list keyword and the second column is the vector within the keyword
#' @param irregular_list list whose lengths are of variable length, making it difficult to convert into dataframe
#' @param col_names vector of names to use for newly created two column dataframe
#' @keywords utility
#' @export
#' @examples
#' irregular_list_to_dataframe()


irregular_list_to_dataframe <- function(irregular_list, col_names = c("keyword", "vector")) {
  df <- data.frame()
  for (n in names(irregular_list)) {
    df <- bind_rows(df, data.frame(n, irregular_list[[n]]))
  }
  
  if (length(col_names) != 2) {
    print(paste0("Column name vector is not two, reverting to defaults: ", col_names))
    col_names = c("keyword", "vector")
  }
  
  names(df) <- col_names
  return(df)
}