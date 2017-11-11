# INTERNAL FUNCTION
# shortens strings to n if they are longer than n
rs.truncateString = function( string, n=20, leeway=3, trail="...") {
  if (nchar(string) >= n ) {
    return( paste0( substr(string, 1, (n-leeway) ), trail ) )
  } else {
    return(string)
  }
}
