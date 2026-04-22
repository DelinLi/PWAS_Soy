## This is a Manhattan plot from GAPIT and modified by Delin Li
GAPIT.Manhattan.plot <- function(GI.MP = NULL, name.of.trait = "Trait",
                                 #plot.type = "Genomewise", DPP = 50000,
                                 cutOff = 0.01, seqQTN = NULL,
                                 color1 = "orangered", color2 = "navyblue",
                                 highliht.sig = FALSE,
                                 color1.sig = "orangered", color2.sig = "navyblue",
                                 cex.axis = 1.5, cex.lab = 2.2, cex = 1.4,
                                 cex.sig = 1.8, pvalue.max = NA) {
 
  if (is.null(GI.MP)) return()
  
  # ---- data cleaning  ----
  P.values <- as.numeric(GI.MP[, 3])
  borrowSlot <- 4
  GI.MP[, borrowSlot] <- 0
  if (!is.null(seqQTN)) GI.MP[seqQTN, borrowSlot] <- 1
  
  idx <- which(GI.MP[, borrowSlot] == 1 & is.na(GI.MP[, 3]))
  GI.MP[idx, 3] <- 1
  
  # Convert to numeric matrix 
  GI.MP <- matrix(as.numeric(as.matrix(GI.MP)), nrow(GI.MP), ncol(GI.MP))
  
  # Remove invalid entries
  GI.MP <- GI.MP[!is.na(GI.MP[, 1]), ]
  GI.MP <- GI.MP[!is.na(GI.MP[, 2]), ]
  GI.MP <- GI.MP[!is.na(GI.MP[, 3]), ]
  GI.MP <- GI.MP[GI.MP[, 3] > 0 & GI.MP[, 3] <= 1, ]
  GI.MP <- GI.MP[GI.MP[, 1] != 0, ]
  
  # Bonferroni cutoff line
  bonferroniCutOff <- -log10(cutOff)
  
  # Convert p‑values to -log10
  GI.MP[, 3] <- -log10(GI.MP[, 3])
  
  # Y‑axis limit
  y.lim <- ifelse(is.na(pvalue.max),
                  as.integer(ceiling(max(GI.MP[, 3]))) * 1.2,
                  pvalue.max)
  print(paste("The max -logP value is", y.lim))
  
  chm.to.analyze <- unique(GI.MP[, 1])
  chm.to.analyze <- chm.to.analyze[order(chm.to.analyze)]
  numCHR <- length(chm.to.analyze)
  
  # Order by chromosome then position
  GI.MP <- GI.MP[order(GI.MP[, 2]), ]
  GI.MP <- GI.MP[order(GI.MP[, 1]), ]
  
  # Build cumulative positions (ticks)
  ticks <- NULL
  lastbase <- 0
  for (i in chm.to.analyze) {
    idx <- (GI.MP[, 1] == i)
    ticks <- c(ticks, lastbase + mean(GI.MP[idx, 2]))
    GI.MP[idx, 2] <- GI.MP[idx, 2] + lastbase
    lastbase <- max(GI.MP[idx, 2])
  }
  
  x0 <- as.numeric(GI.MP[, 2])
  y0 <- as.numeric(GI.MP[, 3])
  z0 <- as.numeric(GI.MP[, 1])
  
  # No pruning – plot all points (original had pruning commented out)
  x <- x0
  y <- y0
  z <- z0
  
  # Points with varying line width (original logic)
  size <- 1
  ratio <- 5
  base <- 1
  themax <- max(y)
  themin <- min(y)
  wd <- ((y - themin + base) / (themax - themin + base)) * size * ratio
  # s <- size - wd / ratio / 2   # unused
  
  # Color vector for chromosomes (cycles through color1, color2)
  mycols <- rep(c(color1, color2), max(z))
  
  # Base plot
  plot(y ~ x,
       ylab = expression(-log[10](italic(p))),
       ylim = c(0, y.lim),
       xaxs = "i", yaxs = "i",
       xlim = c(0, max(x) * 1.01),
       cex.axis = cex.axis, cex.lab = cex.lab,
       col = mycols[z], axes = FALSE,
       type = "p", pch = 20, lwd = wd, cex = cex,
       xlab = "Chromosome",
       main = list(name.of.trait, cex = 2.5),
       mgp = c(3, 1.8, 0))
  
  # Highlight significant points if requested
  if (highliht.sig) {
    sig_idx <- y >= bonferroniCutOff
    if (any(sig_idx)) {
      x.sig <- x[sig_idx]
      y.sig <- y[sig_idx]
      z.sig <- z[sig_idx]
      mycols.sig <- rep(c(color1.sig, color2.sig), max(z.sig))
      points(x.sig, y.sig, type = "p", pch = 20,
             col = mycols.sig[z.sig], cex = cex.sig)
    }
  }
  
  # Vertical lines for QTNs
  QTN <- GI.MP[GI.MP[, borrowSlot] == 1, , drop = FALSE]
  if (nrow(QTN) > 0) abline(v = QTN[, 2], lty = 2, lwd = 1.5, col = "grey")
  
  # Bonferroni horizontal line
  abline(h = bonferroniCutOff, col = "dimgray", lty = 2, lwd = 1)
  
  # Custom axes
  axis(1, at = ticks, tck = -0.01, cex.axis = cex.axis,
       labels = chm.to.analyze, tick = TRUE, lwd = 1.5, padj = -1)
  axis(2, tck = -0.01, cex.axis = cex.axis, lwd = 1.5, padj = 1)
  box()
  
  
}
