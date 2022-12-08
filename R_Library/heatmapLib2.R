
colorTbl <- rep(NA,0)
colorTbl[1] <- colors()[90]  #"darkorange"
colorTbl[2] <- colors()[553] #"red1"
colorTbl[3] <- colors()[47]  #"chartreuse"
colorTbl[4] <- colors()[12]  #"aquamarine4"
colorTbl[5] <- colors()[461] #"mediumblue"
colorTbl[6] <- colors()[526] #"lightsalmon2"
colorTbl[7] <- colors()[429] #"lightseagreen"
colorTbl[8] <- colors()[234] #"gray81"
colorTbl[9] <- colors()[652] #"yellow"
colorTbl[10] <- colors()[624] #"tan4"
colorTbl[11] <- colors()[550] #"purple3"
colorTbl[12] <- colors()[173] #"gray20"
colorTbl[13] <- colors()[103] #"darkseagreen1"
colorTbl[14] <- "gray"
colorTbl[15] <- "pink"

## RowSideLabelsSrt=0
## RowSideLabelsCex=1
## RowSideLabelsOffset=0
## RowSideLabelsAdj=c(0.5,0.5)
## RowSideLabelsAdj2=-0.01
## RowSideTitleCex=1
## RowSideTitleLine=1

# modification of heatmap()
# with multiple side bars
# cols = color palette for the heatmap
heatmap2 <- function (x,
                      Rowv = NULL,
                      Colv = if (symm) "Rowv" else NULL,
                      distfun = dist,
                      hclustfun = hclust,
                      reorderfun = function(d,w) reorder(d, w),
                      add.expr,
                      symm = FALSE,
                      revC = identical(Colv,"Rowv"),# for non-symmetric x, this reverses rows
                      revR = FALSE, # for non-symmetric x, this reverses columns
                      scale = c("none", "row", "column"),
                      na.rm = TRUE,
                      margins = c(5, 5),
                      ColSideColors,
                      RowSideColors,
                      RowSideLabelsSrt=0,
                      RowSideLabelsCex=1,
                      RowSideLabelsOffset=0,
                      RowSideLabelsAdj=c(0.5,0.5),
                      RowSideLabelsAdj2=-0.01,
                      uqLabels=NA,
                      RowSideTitleCex=1,
                      RowSideTitleLine=1,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      boxlwd=3,
                      labRow = NULL,
                      labCol = NULL,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      xlas = 1,
                      xline = -0.5,
                      kbmar=1,
                      klmar=0,
                      ktmar=1,
                      show.key=FALSE,
                      keep.dendro = FALSE,
                      verbose = getOption("verbose"), ...)
{
    scale <- if (symm && missing(scale))
      "none"
    else match.arg(scale)

    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("'x' must be a numeric matrix")

    nr <- di[1L]
    nc <- di[2L]

    if (nr <= 1 || nc <= 1)
        stop("'x' must have at least 2 rows and 2 columns")

    if (!is.numeric(margins) || length(margins) != 2L)
        stop("'margins' must be a numeric vector of length 2")

    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)

    if (!doRdend && identical(Colv, "Rowv"))
        doCdend <- FALSE

    if (is.null(Rowv))
        Rowv <- rowMeans(x, na.rm = na.rm)

    if (is.null(Colv))
        Colv <- colMeans(x, na.rm = na.rm)

    if (doRdend) {
      if (inherits(Rowv, "dendrogram"))
        { ddr <- Rowv }
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv)
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr)))
            stop("row dendrogram ordering gave index of wrong length")
    } else{ rowInd <- 1L:nr}

    ##rowInd <- 1L:nr

    if (doCdend) {
        if (inherits(Colv, "dendrogram"))
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc)
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) x else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv)
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc)))
            stop("column dendrogram ordering gave index of wrong length")
    } else colInd <- 1L:nc

    x <- x[rowInd, colInd]

    labRow <- if (is.null(labRow)) {
      if (is.null(rownames(x)))
        (1L:nr)[rowInd]
      else rownames(x)
    } else labRow[rowInd]

    labCol <- if (is.null(labCol)){
        if (is.null(colnames(x)))
            (1L:nc)[colInd]
        else colnames(x)
    } else labCol[colInd]

    if (scale == "row") {
        x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 1L, sd, na.rm = na.rm)
        x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    }
    else if (scale == "column") {
        x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 2L, sd, na.rm = na.rm)
        x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
    }

    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 4)

    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) ||
          (!is.matrix(ColSideColors) && length(ColSideColors) != nc) ||
          ( is.matrix(ColSideColors) &&   nrow(ColSideColors) != nc))
        stop("'ColSideColors' must be a character vector of length ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1L], 0.2, lhei[2L])

      if ( is.matrix(ColSideColors) && ncol(ColSideColors) > 1 ){
        for ( i in 2:ncol(ColSideColors) ){
          lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2:nrow(lmat), ] + 1)
          lhei <- c(lhei[1L], 0.2, lhei[2:nrow(lmat)])
        }
      }
    }

    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) ||
          (!is.matrix(RowSideColors) && length(RowSideColors) != nr) ||
          ( is.matrix(RowSideColors) && nrow(RowSideColors) != nr))
        stop(paste("'RowSideColors' must be a character vector of length nrow(x) or df with nrow=",nr," nrow(RowSideColors)=",nrow(RowSideColors)))

      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1L], 0.2, lwid[2L])

      if ( is.matrix(RowSideColors) && ncol(RowSideColors) > 1 ){
        for ( i in 2:ncol(RowSideColors) ){
          lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2:ncol(lmat)] + 1)
          lwid <- c(lwid[1L], 0.2, lwid[2:ncol(lmat)])
        }
      }
    }

    lmat[is.na(lmat)] <- 0

    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei,
            "; lmat=\n")
        print(lmat)
    }

    if (!symm || scale != "none")
        x <- t(x)

    if (revC) {
        iy <- nr:1
        if (doRdend)
            ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1L:nr

    if (revR) {
      ix <- nc:1
      if (doCdend)
        ddc <- rev(ddc)
      x <- x[ix,]
    }
    else ix <- 1L:nc

    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    lmat <- rbind(c(max(lmat)+1,rep(0,ncol(lmat)-1)),lmat)
    lhei <- c(0.5,lhei)
    layout(lmat, widths = lwid, heights = lhei, respect = TRUE)

    if (!missing(RowSideColors)) {
      par(mar = c(margins[1L], 0, 0, 0.5))

      if ( !is.matrix(RowSideColors) )
        RowSideColors <- cbind(RowSideColors)

      for ( i in 1:ncol(RowSideColors) ){
        if (revC)
          image(rbind(1L:nr), col = rev(RowSideColors[rowInd,i]), axes = FALSE)
        else
          image(rbind(1L:nr), col = RowSideColors[rowInd,i], axes = FALSE)

        if ( length(colnames(RowSideColors)) ){
          op2 <- par(xpd = NA)
          mtext(colnames(RowSideColors)[i],side=3,las=3,
                line=RowSideTitleLine, cex=RowSideTitleCex)
          #title(colnames(RowSideColors)[i], las=3, line=RowSideTitleLine, cex.main=RowSideTitleCex * op[["cex.main"]])
          par(op2)
        }

        if ( i==1 && length(rownames(RowSideColors)) ) {

          if (revC)
            labels <- rev(rownames(RowSideColors)[rowInd])
          else
            labels <- rownames(RowSideColors)[rowInd]

          if ( length(uqLabels)==1 )
              uqLabels <- unique(labels)

          nUqLabels <- length(uqLabels)
          labelPos <- numeric(nUqLabels)
          for ( i in 1:nUqLabels )
              labelPos[i] <- median((1:length(labels))[labels==uqLabels[i]])
          ##labelPos[i] <- median((1:length(labels))[labels==i])#uqLabels[i]])

          ## if (revC)
          ##   labelPos <- rev(labelPos)

          labelPos <- labelPos / length(labels)
          f <- approxfun(c(0,1),c(0+RowSideLabelsAdj2,1))
          labelPos <- f(labelPos)
          ##text(rep(0.1,nUqLabels), labelPos, labels=1:nUqLabels, cex=RowSideLabelsCex,#pos=1,
          text(rep(0.1,nUqLabels), labelPos, labels=uqLabels, cex=RowSideLabelsCex,#pos=1,
               offset=RowSideLabelsOffset, adj=RowSideLabelsAdj, srt=RowSideLabelsSrt)
        }
      }
    }

    if (!missing(ColSideColors)) {
      par(mar = c(0.5, 0, 0, margins[2L]))

      if ( !is.matrix(ColSideColors) )
        ColSideColors <- cbind(ColSideColors)

      for ( i in 1:ncol(ColSideColors) ){
        if (revR)
          image(cbind(1L:nc), col = rev(ColSideColors[colInd,i]), axes = FALSE)
        else
          image(cbind(1L:nc), col = ColSideColors[colInd,i], axes = FALSE)
      }
    }

    par(mar = c(margins[1L], 0, 0, margins[2L]))

    image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
        c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
    box(lwd=boxlwd)


    ##
    if ( 0 )
    {
        for ( k in 1:17 )
            abline(h=(2*k+0.5),lwd=0.5, col='gray80')
    }

    ## inserted using data defined in Romero/R-scripts/subspeciating_Liners.R
    if ( 0 )
    {
        ## pregnancy status
        for ( br in prStBry )
            abline(v=br,lwd=2, col='red')

        prStLab <- c("NP", "TD", "PTL")
        mtext(prStLab, at=prStMean, side = 1, line = 1.5)

        ## subject ID
        for ( br in subjIdBry2 )
            abline(v=br,lwd=0.5, col='gray80')
        mtext(uqSubjIds33, at=subjIdMean33, side = 3, line = 0.5, las=2, cex=0.5)

        ## community class
        for ( br in cmClBry2 )
            abline(v=br,lwd=1, col='black')
        mtext(cmCl33, at=cmClMean, side = 1, line = 0.2, las=1, cex=0.7)
    } ## end inserted

    if ( 0 )
    {
        ## pregnancy status
        for ( br in prStPwBry )
            abline(v=br,lwd=2, col='red')

        prStPwLab <- c("TD", "PTL")
        mtext(prStPwLab, at=prStPwMean, side = 1, line = 1.5)

        ## subject ID
        for ( br in subjIdPwBry2 )
            abline(v=br,lwd=0.5, col='gray80')
        mtext(uqSubjIds.pw, at=subjIdPwMean, side = 3, line = 0.5, las=2, cex=0.5)
    } ## end inserted


    axis(1, ix, labels = labCol, line = xline, tick = 0,
        cex.axis = cexCol, las=xlas)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1L] - 1.25)

    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
      mtext(ylab, side = 4, line = margins[2L] - 1.25)


    ## op <- par(xpd = NA)
    ## legend(42,40,legend=c("visit 1","visit 2"),fill=c(2,3))
    ## ## legend(42,250,legend=c("Chlamydia Study","400 women study"),fill=c(2,3))
    ## par(op)

    if (!missing(add.expr))
        eval(substitute(add.expr))
    par(mar = c(margins[1L], 0, 0, 0))
    if (doRdend)
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
    if (doCdend)
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main))
        frame()
    if (!is.null(main)) {
        mop <- par(xpd = NA)
        title(main, cex.main = 1.5 * op[["cex.main"]])
        par(mop)
    }

    if ( show.key )
    {
        plot.new()

        ##kop <- par(xpd = NA)
        cols <- rainbow(50,start=1/6,end=0)
        ##cols <- RowSideColors
        min.raw <- min(x, na.rm = TRUE)
        max.raw <- max(x, na.rm = TRUE)
        breaks <- length(cols) + 1
        breaks <- seq(min.raw, max.raw, length=breaks)
        z <- seq(min.raw, max.raw, length=101)
        par(mar = c(kbmar, klmar, ktmar, 0)+0.1,mgp=c(2,0.5,0))
        image(z = matrix(z, ncol = 1), breaks=breaks, col=cols,
              xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        ##par(kop)
    }

    invisible(list(lmat=lmat,lhei=lhei,lwid=lwid,
                   labelPos=labelPos, uqLabels=uqLabels,
                   rowInd = rowInd, colInd = colInd,
                   Rowv = if (keep.dendro && doRdend) ddr,
                   Colv = if (keep.dendro && doCdend) ddc))
}

corHeatmap <- function(x.cor, hcMethod="ward", distMethod="euclidean"){

  hc <- hclust(dist(x.cor),method=hcMethod) # hierarchical clustering
  rowInd <- hc$order
  colInd <- hc$order
  x.cor <- x.cor[rowInd, colInd]

  nr <- nc <- nrow(x.cor)
  taxId <- colnames(x)[colInd]

  op <- par(mar = c(11, 4, 4, 11) + 0.1)
  image(x.cor[,ncol(x.cor):1] ,col=rainbow(50,start=0,end=4/6),axes=FALSE,xlab="",ylab="")
  axis(4,seq(from=0,to=1,len=nc),labels=rev(taxId),las=1,line=-0.5,tick=0,cex.axis=rcex)
  axis(1,seq(from=0,to=1,len=nr),labels=taxId,las=2,line=-0.5,tick=0,cex.axis=cexCol)
  par(op)
}

scale01 <- function(x, low = min(x), high = max(x)) {
  x <- (x - low)/(high - low)
  x
}

# version of heatmap that accepts output of hclust
# (saves time to rerun dist and hclust)
locHeatmap <- function (x, hc, reorderfun = function(d, w) reorder(d, w),
                        cols = rainbow(50,start=1/6,end=0),
                        doRdend=TRUE, doCdend=FALSE,
                        add.expr, na.rm = TRUE, margins = c(5, 5),
                        ColSideColors,
                        RowSideColors,
                        RowSideTitleCex=0.7,
                        RowSideTitleLine=0.5,
                        RowSideColors2,
                        RowSideRange2,
                        RowSideColorsPalette2,
                        RowSideTitle2, # title of RowSideColors2 color bar
                        RowSideLabels, # clustering membership vector; cluster numbers will be placed in the cluster color bar
                        RowSideLabelsCex=1,
                        ColorKeyTitle="Color Key",
                        ColorKeyTitle2,# title of color key for RowSideColors2
                        cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc),
                        labRow = NULL, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL,
                        verbose = getOption("verbose"), ...)
{
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("'x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("'x' must have at least 2 rows and 2 columns")

  rowInd <- hc$order
  colInd <- 1L:nc

  ddr <- as.dendrogram(hc)
  #ddr <- reorder(ddr, NULL)

  x <- x[rowInd, colInd]

  if (is.null(labRow)){

    if (is.null(rownames(x)))
      labRow <- (1L:nr)[rowInd]
    else
      labRow <- rownames(x)

  } else if (!is.na(labRow) & length(labRow)==length(rowInd)) {
    labRow <- labRow[rowInd]
  }

  if (is.null(labCol)){
    if (is.null(colnames(x)))
      labCol <- (1L:nc)[colInd]
    else
      labCol <- colnames(x)
  }

  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doRdend) 1 else 0.05, 4)
  lhei <- c(0.05 + if (!is.null(main)) 0.2 else 0, 4)

  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) != nr)
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
    lwid <- c(lwid[1], 0.2, lwid[2])
  }

  if (!missing(RowSideColors2)) {
    if (!is.character(RowSideColors2) || length(RowSideColors2) != nr)
      stop("'RowSideColors2' must be a character vector of length nrow(x)")
    lmat <- rbind(c(NA,NA,NA,5), c(4,1,2,3))
    lwid <- c(1.0, 0.2, 0.2, 4.0)
  }

  lmat[is.na(lmat)] <- 0
  lhei <- c(0.25, lhei)

  if (missing(RowSideColors2)){
    lmat <- rbind(c(5,0,0),lmat)
  }else{
    lmat <- rbind(c(7,0,0,0),c(6,0,0,0),lmat)
    lhei <- c(0.41, lhei)
  }

  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei,
        "; lmat=\n")
    print(lmat)
  }

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  #layout.show(5)

  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1L:nr), col = RowSideColors[rowInd], axes = FALSE)
    op2 <- par(xpd = NA)
    title("Cluster", las=1, line=RowSideTitleLine, cex.main=RowSideTitleCex * op[["cex.main"]])
    par(op2)

    if (!missing(RowSideLabels)) {

      labels <- RowSideLabels[rowInd]
      nUqLabels <- length(unique(labels))
      labelPos <- numeric(nUqLabels)
      for ( i in 1:nUqLabels )
        labelPos[i] <- median((1:length(labels))[labels==i])

      labelPos <- labelPos / length(labels)
      text(rep(0.1,nUqLabels), labelPos, labels=1:nUqLabels, cex=RowSideLabelsCex,
           pos=1, offset=0.25)
    }
  }

  if (!missing(RowSideColors2)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1L:nr), col = RowSideColors2[rowInd], axes = FALSE)
    if ( !missing(RowSideTitle2)) {
      op2 <- par(xpd = NA)
      title(RowSideTitle2,las=1,line=RowSideTitleLine, cex.main=RowSideTitleCex * op[["cex.main"]])
      par(op2)
    }
  }

  par(mar = c(margins[1], 0, 0, margins[2]))

  x <- t(x)
  image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col=cols, ...)
  axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)

  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)

  iy <- 1L:nr

  if (length(labRow))
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)

  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)

  if (!missing(add.expr))
    eval(substitute(add.expr))

  par(mar = c(margins[1], 0, 0, 0))

  if (doRdend)
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()

  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
  if (doCdend)
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  else if (!is.null(main))
    frame()

  if (!is.null(main)) {
    par(xpd = NA)
    title(main, cex.main = 1.5 * op[["cex.main"]])
  }

  # Heatmap Color Key
  min.raw <- min(x, na.rm = TRUE)
  max.raw <- max(x, na.rm = TRUE)
  breaks <- length(cols) + 1
  breaks <- seq(min.raw, max.raw, length=breaks)

  z <- seq(min.raw, max.raw, length=101)

  ##par(mar = c(0, 0, 2, 3)+0.1)
  par(mar = c(0, 0, 0, 0)+0.1)
  image(z = matrix(z, ncol = 1), breaks=breaks, col=cols,
        xaxt = "n", yaxt = "n")

  par(usr = c(0, 1, 0, 1))
  lv <- pretty(breaks)
  xv <- scale01(as.numeric(lv), min.raw, max.raw)
  axis(1, at = xv, labels = lv)
  #mtext(side = 1, "Value", line = 2.5)
  title(ColorKeyTitle)

  # RowSideColors Color Key
  if ( !missing(RowSideColors2) & !missing(RowSideRange2) ){

    cols <- RowSideColorsPalette2
    min.raw <- RowSideRange2[1]
    max.raw <- RowSideRange2[2]
    breaks <- length(cols) + 1
    breaks <- seq(min.raw, max.raw, length=breaks)

    z <- seq(min.raw, max.raw, length=101)

    par(mar = c(2, 0, 2, 3)+0.1)
    image(z = matrix(z, ncol = 1), breaks=breaks, col=cols,
          xaxt = "n", yaxt = "n")

    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if ( !missing(ColorKeyTitle2))
      title(ColorKeyTitle2)
  }
}


# adding 'method' argument to R heatmap(), which is passed to hclustfun
locHeatmap2 <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,
                     distfun = dist, hclustfun = hclust, hclmethod = "ward",
                     reorderfun = function(d, w) reorder(d, w),
                     add.expr, symm = FALSE, revC = identical(Colv,"Rowv"),
                     scale = c("row", "column", "none"), na.rm = TRUE,
                     margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 +
                     1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
                     labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE,
                     verbose = getOption("verbose"), ...)
{
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)

  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("'x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]

  if (nr <= 1 || nc <= 1)
    stop("'x' must have at least 2 rows and 2 columns")

  if (!is.numeric(margins) || length(margins) != 2)
    stop("'margins' must be a numeric vector of length 2")

  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)

  if (is.null(Rowv))
    Rowv <- rowMeans(x, na.rm = na.rm)

  if (is.null(Colv))
    Colv <- colMeans(x, na.rm = na.rm)

  if (doRdend) {
    if (inherits(Rowv, "dendrogram"))
      ddr <- Rowv
    else {
      hcr <- hclustfun(distfun(x),method=hclmethod)
      ddr <- as.dendrogram(hcr)
      if (!is.logical(Rowv) || Rowv)
        ddr <- reorderfun(ddr, Rowv)
    }
    if (nr != length(rowInd <- order.dendrogram(ddr)))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else rowInd <- 1L:nr

  if (doCdend) {
    if (inherits(Colv, "dendrogram"))
      ddc <- Colv
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      ddc <- ddr
    }
    else {
      hcc <- hclustfun(distfun(if (symm)
                               x
      else t(x)))
      ddc <- as.dendrogram(hcc)
      if (!is.logical(Colv) || Colv)
        ddc <- reorderfun(ddc, Colv)
    }
    if (nc != length(colInd <- order.dendrogram(ddc)))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else colInd <- 1L:nc

  x <- x[rowInd, colInd]

  labRow <- if (is.null(labRow))
    if (is.null(rownames(x)))
      (1L:nr)[rowInd]
    else rownames(x)
  else labRow[rowInd]
  labCol <- if (is.null(labCol))
    if (is.null(colnames(x)))
      (1L:nc)[colInd]
    else colnames(x)
  else labCol[colInd]
  if (scale == "row") {
    x <- sweep(x, 1, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/", check.margin = FALSE)
  }
  else if (scale == "column") {
    x <- sweep(x, 2, colMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/", check.margin = FALSE)
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doRdend) 1 else 0.05, 4)
  lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0,
            4)
  if (!missing(ColSideColors)) {
    if (!is.character(ColSideColors) || length(ColSideColors) !=
        nc)
      stop("'ColSideColors' must be a character vector of length ncol(x)")
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1], 0.2, lhei[2])
  }
  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) !=
        nr)
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1),
                                   1), lmat[, 2] + 1)
    lwid <- c(lwid[1], 0.2, lwid[2])
  }
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei,
        "; lmat=\n")
    print(lmat)
  }
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1L:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  if (!symm || scale != "none")
    x <- t(x)
  if (revC) {
    iy <- nr:1
    ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1L:nr
  image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
        c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
  axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  par(mar = c(margins[1], 0, 0, 0))
  if (doRdend)
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
  if (doCdend)
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  else if (!is.null(main))
    frame()
  if (!is.null(main)) {
    par(xpd = NA)
    title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro &&
                                                     doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}


