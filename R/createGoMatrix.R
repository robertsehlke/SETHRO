#' @title createGoMatrix
#' 
#' @description 
#' Using a bioconductor gene ontology annotation package, create a sparse binary matrix that maps gene identifiers (column names) to the 
#' GO terms they are annotated with (row names).
#'
#' @param egGO2ALLEGS GO to ALL identifiers mapping from any bioconductor annotation package. E.g. org.Dm.egGO2ALLEGS from the
#' org.Dm.eg.db package. 
#' 
#' @param min.genes Minimum size of gene set to be considered, defaults to 5.
#' 
#' @param max.genes Maximum size of gene set to be considered, defaults to 300.
#' 
#' @param excludeIEA
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

createGoMatrix = function( egGO2ALLEGS, min.genes=5, max.genes=300, excludeIEA=F) {
  
  golist = as.list(egGO2ALLEGS)
  golist = lapply( golist, function(x) {
    if (excludeIEA) { 
      x = x[ which(!names(x) %in% "IEA") ]
    }
    return( unique(x) )
  })
  
  gn = sapply( golist,function(x) length( x ) )
  golist = golist[gn >= min.genes & gn <= max.genes]
  
  ugenes = unique( unlist(golist))
  
  gomat = Matrix(0, ncol=length(ugenes), nrow=length(golist), dimnames = list(names(golist), ugenes) )
  
  for(i in 1:nrow(gomat)) { gomat[i,] = colnames(gomat) %in% golist[[i]] }
  
  return(gomat)
}
