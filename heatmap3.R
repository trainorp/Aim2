heatmap3<-function (distForPlot, Rowv = NULL, Colv = NULL, 
          distfun = dist, hclustfun = hclust, 
          reorderfun = function(d, w) reorder(d, w),
          symm = FALSE, revC = identical(Colv,"Rowv"), scale ="none", na.rm = TRUE, 
          margins = c(5, 5), ColSideColors=NULL, RowSideColors=NULL, 
          cexRow = (0.2 + 1/log10(nrow(distForPlot)))*.2, 
          cexCol = (0.2 + 1/log10(ncol(distForPlot)))*.2, labRow = NULL, 
          labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
          verbose = getOption("verbose"), ...) 
{
  library(colorRamps)
  nr <- nc <- ncol(distForPlot) 
  if (!is.numeric(margins) || length(margins) != 2L) 
    stop("'margins' must be a numeric vector of length 2")
  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)
  if (doRdend) {
      hcr <- hclust(as.dist(distForPlot),method="ward.D2")
      ddr <- as.dendrogram(hcr)
  }
  if (doCdend) {
      ddc <- ddr
  }
  labRow <- rownames(distForPlot)
  labCol <- colnames(distForPlot)

  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doRdend) 1 else 0.05, 4)
  lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 4)

  lmat[is.na(lmat)] <- 0

  dev.off()
  png(file="MISimClustWardD2Color.png",height=12,width=12,units="in",res=600)
  dev.hold()
  on.exit(dev.flush())
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  
  par(mar = c(margins[1L], 0, 0, margins[2L]))

  iy <- 1L:nr
  image(1L:nc, 1L:nr, distForPlot[hcr$order,hcr$order], xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "",col=rev(blue2green2red(20)))
  
  axis(1, 1L:nc, labels = labCol[hcr$order], las = 2, line = -0.5, tick = 0, 
       cex.axis = cexCol)
  
  axis(4, iy, labels = labRow[hcr$order], las = 2, line = -0.5, tick = 0, cex.axis = cexRow)

  par(mar = c(margins[1L], 0, 0, 0))
  if (doRdend) 
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none",edgePar=list(lwd=.4))
  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
  if (doCdend) 
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  dev.off()
  
  png(file="MIcompSimDendColor.png",height=10,width=12,units="in",res=300)
  plot(ddr,nodePar=list(pch=NA,cex=.5,lab.cex=.25))
  dev.off()
  
  # invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
  #                                                             doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}