#' @title plotClusteredEnrichment
#' 
#' @description 
#' Simple barplot that displays the enrichment of primary gene sets from a clusteredEnrichment result.
#'
#' @param ceresult $results data frame returned by the clusteredEnrichmentGO() function.
#' 
#' @param alpha Significance cutoff for primary gene set adjusted p-values, defaults to 0.05.
#' 
#' @param t.alpha If a directionality bias has been determined via t-test on the log-fold-change of the gene set,
#' which significance cutoff should be used for the display of that bias? (i.e. bar color)
#' 
#' @param max_bars If set, scales the bar plot to provide space for max_bars... bars. Use to ensure uniform scaling 
#' when plotting multiple results with varying numbers of significant hits.
#' 
#' @param top_n Plot only the most significant top_n results. Default is to show all results.
#' 
#' @param only_positive_enrichment Defaults to TRUE: show only positive log-enrichments (i.e. no depleted terms).
#' 
#' @param order_by_enrichment Defaults to TRUE: order significant hits by effect magnitude (enrichment) rather than
#' significance.
#' 
#' @param legend.title Title of the legend (only relevant if directionality bias of GO terms is included in the results).
#' 
#' @param xlim Axis limit of the enrichment plot. Use to ensure uniform axis across multiple plots.
#' 
#' @param lab_max_chars Maximum length of gene set names; longer names are truncated.
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

plotClusteredEnrichment = function( ceresult, 
                                    alpha = 0.05, 
                                    t.alpha = 0.05, 
                                    max_bars = NULL, 
                                    top_n = Inf,
                                    only_positive_enrichment=T, 
                                    order_by_enrichment=T, 
                                    main="Enriched primary gene sets",
                                    legend.title = "", 
                                    xlim=NULL, 
                                    lab_max_chars=50, ... ) {
  par(mar=par("mar") + c(0,16,3,4), xpd=F)
  
  # Thresholding and sorting
  ceresult = ceresult[which(ceresult$Primary.Fisher.adj <= alpha),,drop=F ]
  if(only_positive_enrichment) { ceresult = ceresult[which(ceresult$Enrichment > 0),,drop=F ] }
  if(order_by_enrichment) { ceresult = ceresult[order(ceresult$Enrichment),,drop=F ] } else {
    ceresult = ceresult[order(ceresult$Fisher, decreasing = T),,drop=F ]
  }
  
  # top n if required
  if(top_n < nrow(ceresult)) {
    ceresult = ceresult[(nrow(ceresult)-top_n):nrow(ceresult),]
  }
  
  cols = rep( "gray23", nrow(ceresult) )
  # If directionality has been determined via t-test, color bars accordingly
  if (!is.null(ceresult$t_up)) {
    cols[ ceresult$t_up.adj <= t.alpha ] = "coral"
    cols[ ceresult$t_down.adj <= t.alpha ] = "deepskyblue2"
  }
  
  # Add buffer bars to ensure consistency
  if ( !is.null(max_bars) ) { baradjust = max_bars - nrow(ceresult) } else { baradjust = 0 }
  
  # Plot the actual thing
  if(is.null(xlim)) { xlim = c(0, ceiling( max(ceresult$Enrichment) ) ) }
  bp = barplot( c( rep(NA, baradjust), ceresult$Enrichment ),  
                names.arg = c( rep(NA, baradjust), sapply( ceresult$Term, rs.truncateString, n=lab_max_chars ) ), 
                horiz = T, las=1, 
                col = c( rep(NA, baradjust), cols), 
                lwd=2,
                xlim = xlim, 
                axes = F, ... )
  axis(side = 3)
  mtext("log2(enrichment)", line=2.5)
  title(main, line = 5, ...)
  
  # And add a legend
  if (!is.null(ceresult$t_up)) {
    legend(x = max(xlim), y=bp[baradjust+1], fill=c("coral","deepskyblue2","gray23"),
           legend=c("up","down","indet."), xpd = T, bty = "n", y.intersp = 0.9, title=legend.title)
  }
  
  # add number of genes
  par(xpd=T)
  text(x = c( rep(NA, baradjust), ceresult$Enrichment ), y = bp, 
       c(rep(NA, baradjust), paste0("*", ceresult$Significant, "/", ceresult$Annotated) ), pos = 4)
  
  par(mar=par("mar") - c(0,16,3,4), xpd=F)
  
  invisible(ceresult)
}