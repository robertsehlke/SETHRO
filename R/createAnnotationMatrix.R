#' @title createAnnotationMatrix
#' 
#' @description 
#' From a list of annotations, create a sparse binary matrix that maps gene identifiers (column names) to the 
#' GO terms they are annotated with (row names).
#'
#' @param annotationlist Named list of character vectors. Names are assumed to be gene set identifiers, character vectors member gene identifiers.
#' 
#' @param min.genes Minimum size of gene set to be considered, defaults to 5.
#' 
#' @param max.genes Maximum size of gene set to be considered, defaults to 300.
#' 
#' @author
#' Robert Sehlke
#' 
#' @examples
#' \dontrun{
#' # example for annotation matrix creation
#' # may be used interchangably with createGoMatrix()
#' gomatrix = createAnnotationMatrix( list("GO:0000042"=paste0("gene",1:6),
#'                                            "GO:0000002"=paste0("gene",2:7),
#'                                            "GO:0000009"=paste0("gene",6:12)) )
#' 
#' # example call for main function
#' res = clusteredEnrichmentGO(gomatrix, 
#'                             geneset = paste0("gene",2:6),
#'                             universe = paste0("gene",1:12))
#' 
#' # example plot
#' plotClusteredEnrichment(res$results)
#' }
#' 
#' @export

createAnnotationMatrix = function( annotationlist, min.genes=5, max.genes=300) {
  
  gn = sapply( annotationlist,function(x) length( x ) )
  annotationlist = annotationlist[gn >= min.genes & gn <= max.genes]
  
  ugenes = unique( unlist(annotationlist) )
  
  gomat = Matrix(0, ncol=length(ugenes), nrow=length(annotationlist), dimnames = list(names(annotationlist), ugenes) )
  
  for(i in 1:nrow(gomat)) { gomat[i,] = colnames(gomat) %in% annotationlist[[i]] }
  
  return(gomat)
}