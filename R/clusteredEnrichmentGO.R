#' @title clusteredEnrichmentGO
#' 
#' @description 
#' Gene Ontology enrichment with subsequent clustering of terms
#'
#' @param gomatrix Sparse binary matrix, maps gene identifiers (column names) to the GO terms they are annotated with (row names).
#' 
#' @param geneset Character vector, gene set of interest.
#' 
#' @param universe Character vector, background set of genes. Note: gene identifiers not included in the gomatrix are dropped!
#' 
#' @param geneAlias Named character vector, maps gene identifiers (vector names) to alternative identifiers.
#' 
#' @param min.genes Minimum number of genes of interest that a term needs to contain to be included in the analysis.
#' 
#' @param cut_max Terms may be combined into a cluster if they contain at most cut_max different genes of interest (i.e. max. edit distance).
#' 
#' @param dynamic_cut Default: FALSE. If TRUE, evaluate the silhouette index of all clusterings with max. edit distances of [1...cut_max], and choose the optimal cut accordingly.
#' 
#' @param domains Character vector, select included GO terms by domain: biological process ("BP"), cellular compartment ("CC"), molecular function ("MF").
#' 
#' @param two_sided Default: FALSE. If TRUE, performs two-sided Fisher tests to detect depletion as well as enrichment.
#' 
#' @param representative_by_silhouette Default: FALSE. If TRUE, selects the representative of a cluster based on their silhouette score 
#' (i.e. the most specific term for that cluster). Otherwise, the smallest term is chosen, or the term with the shortest description in the presence of ties.
#' 
#' @param sort_by_pval Default: TRUE. Sort results by p-value.
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

clusteredEnrichmentGO = function( gomatrix, 
                                  geneset, 
                                  universe, 
                                  geneAlias=NULL, 
                                  min.genes=5, 
                                  cut_max = 5, 
                                  dynamic_cut=F,
                                  domains=c("BP","MF","CC"), 
                                  two_sided=F,
                                  representative_by_silhouette = F, 
                                  sort_by_pval=T,
                                  weights=NULL) {
  
  if(!is.null(weights)) {
    weights = setNames(weights, geneset)
  }
  
  # for now, only support explicit universe
  universe = intersect(colnames(gomatrix), unique(universe) )
  gomatrix = gomatrix[,universe]
  
  # subset by domains
  gomatrix = gomatrix[ Ontology(rownames(gomatrix)) %in% domains, ]
  
  
  # make sure the geneset is congruent with the GO matrix (also, no duplicates!)
  gun = which(!duplicated(geneset))
  geneset = geneset[gun]
  if(!is.null(geneAlias)) { geneAlias = geneAlias[gun] }
  
  gidx = which(geneset %in% colnames(gomatrix))
  geneset = geneset[gidx]
  if(!is.null(geneAlias)) { geneAlias = geneAlias[gidx] }
  
  
  
  # subset according to min-genes (minimum refers to "significant" genes)
  n_in_term = apply(gomatrix, 1, sum)
  n_selected = apply(gomatrix[,geneset,drop=F], 1, sum)
  
  idx = which(n_selected >= min.genes)
  
  if (length(idx)==0) { warning("No GO terms with the specified min.genes were found for the subset of interest!"); return(NA) }
  
  gomatrix = gomatrix[idx,,drop=F]
  n_in_term = n_in_term[idx]
  n_selected = n_selected[idx]
  
  
  
  # Fisher tests for significant enrichment --> data frame with GO = rownames
  # if weights are given, calculate the mean weight for each GO term
  flist = list()
  for (k in 1:nrow(gomatrix)) {
    cmat = rbind( c(sum( gomatrix[k,] == 0 & !(colnames(gomatrix) %in% geneset)), sum(gomatrix[k,] == 1 & !(colnames(gomatrix) %in% geneset )) ),
                  c(sum( gomatrix[k,] == 0 & (colnames(gomatrix) %in% geneset) ), sum( gomatrix[k,] == 1 & (colnames(gomatrix) %in% geneset )) ) )
    
    expected = as.array(margin.table(cmat,1)) %*% t(as.array(margin.table(cmat,2))) / margin.table(cmat)
    enrichment = log2( cmat[2,2] / expected[2,2] )
    pval = fisher.test(cmat, alternative = ifelse(two_sided, "two.sided", "greater") )$p.value
    
    flist[k] = list(c("Significant"=cmat[2,1], 
                      "Expected"=round(expected[2,2], digits = 2) , 
                      "enrichment"=round(enrichment, digits = 2), 
                      "p.value"=pval))
  }
  flist = t( as.data.frame(flist) )
  rownames(flist) = rownames(gomatrix)
  
  
  # now subset the gomatrix to just the geneset of interest
  gomatrix = gomatrix[,geneset,drop=F]
  
  
  # CLUSTER AND SELECT PRIMARY TERMS (only if >1 go term selected)
  if (nrow(gomatrix) > 1) {
    # calculate manhattan distance and cluster
    godist = dist(gomatrix, method="manhattan")
    goclust = hclust(godist, method = "complete")
    
    
    # cut tree so that the edit distance between clustered GO terms is at most max.distance
    if (dynamic_cut) {
      ci = 1:cut_max
      maxindex = -Inf
      for (c in ci) {
        sindex = mean( as.matrix( silhouette( cutree(goclust, h = c + 0.5), godist ) )[,3] )
        if (sindex > maxindex) { maxindex = sindex; cut_max = c }
      }
    }
    
    #return( list(goclust, godist) )
    gopartition = cutree(goclust, h = cut_max + 0.5)
    
    gosil = as.matrix( silhouette( cutree(goclust, h = cut_max + 0.5), godist ) )
    
    if(!is.na(gosil)) {
      rownames(gosil) = goclust$labels
      sindex = mean( as.matrix( silhouette( cutree(goclust, h = cut_max + 0.5), godist ) )[,3] )
    } else {
      sindex = NA
    }
    
    
    
    # Select as representative the smallest GO term in the cluster (by total annotations)
    if (!representative_by_silhouette) {
      reps = c()
      for ( c in 1:max(gopartition) ) {
        this = n_in_term[ names(gopartition)[which(gopartition == c)] ]
        
        if ( sum(this == min(this)) > 1 ) { 
          this = this[this == min(this)]
          this.rep = names(this)[which.min(nchar(Term(names(this))))]    # if several terms apply, choose the one that's... shortest
        } 
        else {
          this.rep = names(this)[which.min(this)]
        }
        reps = c(reps, this.rep)
      }
    } # else, select as representative the term with the highest silhouette index of its cluster
    else { reps = c()
    for ( c in 1:max(gopartition) ) {
      this = gosil[which(gosil[,"cluster"] == c),,drop=F]
      if ( sum(this[,"sil_width"] == max(this[,"sil_width"])) > 1 ) {
        this = this[this[,"sil_width"] == max(this[,"sil_width"]),]
        this.rep = rownames(this)[which.min(nchar(Term(rownames(this))))]    # if several terms apply, choose the one that's... shortest
      } 
      else {
        this.rep = rownames(this)[which.min(this[,"sil_width"])]
      }
      reps = c(reps, this.rep)
    }
    }
  } else {
    gopartition = setNames( 1, rownames(gomatrix) )
    reps = setNames( rownames(gomatrix), 1 )
  }
  
  
  # Summary data frame
  gm = data.frame( "Primary"=Term(reps[gopartition]),
                   "Is.primary" = names(n_in_term) %in% reps,
                   "Terms.in.cluster" = 0,
                   "GO.ID"=names(n_in_term),
                   "Term"=Term(names(n_in_term)), 
                   "Ontology"=Ontology(names(n_in_term)),
                   "Annotated"=n_in_term, 
                   "Significant"=n_selected,
                   "Expected"=flist[,"Expected"],
                   "Enrichment"=2^flist[,"enrichment"],
                   "log2Enrichment"=flist[,"enrichment"],
                   "Fisher"=flist[,"p.value"],
                   "Primary.Fisher.adj"=NA,
                   "Genes"="",
                   "Genes.in.secondary"="",
                   row.names = names(n_in_term), 
                   stringsAsFactors = F)
  gm$Primary.Fisher.adj[gm$Is.primary] = p.adjust( gm$Fisher[gm$Is.primary], method="BH" )
  
  if(!is.null(weights)) { gm$weightAvg = rep(NA,nrow(gm)) }
  
  # Add the names of significant genes in each term & number in cluster
  if(!is.null(weights)) { weights = weights[colnames(gomatrix)] }
  
  for (i in 1:nrow(gm)) {
    if( !is.null(geneAlias)) { gm$Genes[i] = list(geneAlias[ gomatrix[i,] == 1 ]) } else {
      gm$Genes[i] = list(geneset[ gomatrix[i,] == 1 ])
    }
    
    if(!is.null(weights)) {
      gm$weightAvg[i] = mean( weights[ gomatrix[i,] == 1 ] )
    }
    
    gm$Terms.in.cluster[i] = sum(gm$Primary == gm$Primary[i])
  }
  
  # For primary terms only, list genes that are included ONLY in secondary terms
  for (i in 1:nrow(gm)) {
    if( gm$Term[i] %in% gm$Primary ) {
      excl = setdiff( unique( unlist( gm[which(gm$Primary == gm$Term[i]),"Genes"] ) ),
                      gm[i,"Genes"][[1]] )
      if (length(excl) > 0) { gm$Genes.in.secondary[i] = list( excl ) }
    }
  }
  
  # For ease of plotting, translate the labels 
  if (nrow(gomatrix) > 1) { goclust$labels = paste0( strtrim( Term(goclust$labels), 30) ) } else {goclust = godist = sindex = NA}
  
  if (sort_by_pval) { gm = gm[order(gm$Fisher),]}
  
  # lists of genes in data frames create some issues downstream, so let's paste them
  gm$Genes = sapply(gm$Genes, paste0, collapse=", ")
  gm$Genes.in.secondary = sapply(gm$Genes.in.secondary, paste0, collapse=", ")
  
  return( list("results"=gm, "background"=universe, "geneset"=geneset, "gomatrix"=gomatrix, "dists"=godist,
               "clustering"=goclust, "partition"=gopartition, "cut_at"=cut_max + 0.5, "silhouette_index"=sindex ) )
}