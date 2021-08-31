
t_to_res_snp = function(edges,
                        tree,
                        subclades,
                        subclades_tips,
                        subclades_nodes,
                        year_threshold = FALSE,
                        bases = c('A', 'C', 'G', 'T', '-'),
                        anc.seq,
                        drug_db) {
  
  #' Takes a tree object from BactDating and it returns a table
  #' ready for survival analysis.

  t_to_res <- list()
  node = length(tree$tip.label) + 1
  edges_clade = edges
  for (i in 1:nrow(edges_clade)) {
    int_node = edges_clade[i, 1]
    ext_node = edges_clade[i, 2]
    
    lineage = NA
    for (clade in names(subclades)) {
      tips = subclades_tips[[clade]]
      node = subclades_nodes[[clade]]
      
      if (length(node) == 0) {
        next()
      }
      clade_nodes = edges[min(which(edges[, 1] == node)):which(edges[, 2] == max(tips)), c(1, 2)]
      clade_nodes = unique(c(node, clade_nodes$V1, clade_nodes$V2))
      lineage = ifelse (ext_node %in% clade_nodes, clade, lineage)
    }
    
    tips = setNames(1:length(tree$tip.label), tree$tip.label)
    node = length(tree$tip.label) + 1
    
    if (int_node == node) {
      int_time = tree$root.time
    } else {
      int_time = edges[edges[, 2] == int_node, 3]
    }
    ext_time = edges[edges[, 2] == ext_node, 3]
    
    if (!isFALSE(year_threshold)) {
      if (int_time < year_threshold) {
        next()
      }
    }
    
    int_seq = anc.seq[[as.character(int_node)]][attr(anc.seq, 'index'),]
    
    if (ext_node %in% tips) {
      ext_seq = anc.seq[[names(tips[which(tips == ext_node)])]][attr(anc.seq, 'index'),]
    } else {
      ext_seq = anc.seq[[as.character(ext_node)]][attr(anc.seq, 'index'),]
    }
    
    row.names(ext_seq) = 1:nrow(ext_seq)
    ext_res_snps = ext_seq[unique(sort(as.vector(unlist(
      sapply(drug_db$snp, function (x)
        strsplit(x, ';')[[1]])
    )))),]
    ext_res_bases = setNames(apply(ext_res_snps, 1, function(x)
      ifelse (length(which(
        x == max(x)
      )) == 1, bases[which(x == max(x))], 'N')),
      row.names(ext_res_snps))
    ext_res_present = drug_db[apply(drug_db, 1,
                                    function(x)
                                      all(ext_res_bases[strsplit(x['snp'], ';')[[1]]] == strsplit(x['alt'], ';')[[1]]) == TRUE),]
    
    
    row.names(int_seq) = 1:nrow(int_seq)
    int_res_snps = int_seq[unique(sort(as.numeric(as.vector(unlist(
      sapply(drug_db$snp, function (x)
        strsplit(x, ';')[[1]])
    ))))),]
    int_res_bases = setNames(apply(int_res_snps, 1, function(x)
      ifelse (length(which(
        x == max(x)
      )) == 1, bases[which(x == max(x))], 'N')),
      row.names(int_res_snps))
    int_res_present = drug_db[apply(drug_db, 1,
                                    function(x)
                                      all(int_res_bases[strsplit(x['snp'], ';')[[1]]] == strsplit(x['alt'], ';')[[1]]) == TRUE),]
    
t_to_res[['int_node']] = c(t_to_res[['int_node']], int_node)
t_to_res[['ext_node']] = c(t_to_res[['ext_node']], ext_node)
t_to_res[['int_year']] = c(t_to_res[['int_year']], int_time)
t_to_res[['ext_year']] = c(t_to_res[['ext_year']], ext_time)
t_to_res[['int_ci']] = c(t_to_res[['int_ci']], paste(round(tree$CI[int_node, ], 2), collapse = ';'))
t_to_res[['ext_ci']] = c(t_to_res[['ext_ci']], paste(round(tree$CI[ext_node, ], 2), collapse = ';'))
t_to_res[['int_st']] = c(t_to_res[['int_st']], ifelse(nrow(int_res_present) == 0, 0 , 1))
t_to_res[['ext_st']] = c(t_to_res[['ext_st']], ifelse(nrow(ext_res_present) == 0, 0 , 1))
t_to_res[['int_pos']] = c(t_to_res[['int_pos']],
                          paste(int_res_present$pos, collapse = ';'))
t_to_res[['ext_pos']] = c(t_to_res[['ext_pos']],
                          paste(ext_res_present$pos, collapse = ';'))
t_to_res[['int_drug']] = c(t_to_res[['int_drug']],
                           paste(int_res_present$drug, collapse = ';'))
t_to_res[['ext_drug']] = c(t_to_res[['ext_drug']],
                           paste(ext_res_present$drug, collapse = ';'))
t_to_res[['int_gene']] = c(t_to_res[['int_gene']],
                           paste(int_res_present$gene, collapse = ';'))
t_to_res[['ext_gene']] = c(t_to_res[['ext_gene']],
                           paste(ext_res_present$gene, collapse = ';'))
t_to_res[['int_mut']] = c(t_to_res[['int_mut']],
                          paste(int_res_present$mutation, collapse = ';'))
t_to_res[['ext_mut']] = c(t_to_res[['ext_mut']],
                          paste(ext_res_present$mutation, collapse = ';'))
t_to_res[['int_conf']] = c(t_to_res[['int_conf']],
                           paste(int_res_present$confidence, collapse = ';'))
t_to_res[['ext_conf']] = c(t_to_res[['ext_conf']],
                           paste(ext_res_present$confidence, collapse = ';'))
t_to_res[['lineage']] = c(t_to_res[['lineage']], lineage)
  }
  t_to_res = as.data.frame(t_to_res)
  t_to_res[t_to_res == ''] = NA
  return (t_to_res)
}



bl_to_res_snp = function(tree,
                        subclades,
                        subclades_tips,
                        subclades_nodes,
                        bases = c('A', 'C', 'G', 'T', '-'),
                        anc.seq,
                        drug_db) {
  #' Takes a tree object with branch length in genetic distance and it returns a table
  #' ready for survival analysis.
  
  t_to_res <- list()
  node = length(tree$tip.label) + 1
  edges = tree$edge
  edges_clade = edges
  for (i in 1:nrow(edges_clade)) {
    int_node = edges_clade[i, 1]
    ext_node = edges_clade[i, 2]
    
    edge_pos = which(apply(tree$edge, 1, function(x) all(x == edges_clade[i,])))
    
    lineage = NA
    for (clade in names(subclades)) {
      tips = subclades_tips[[clade]]
      node = subclades_nodes[[clade]]
      
      if (length(node) == 0) {
        next()
      }
      clade_nodes = edges[min(which(edges[, 1] == node)):which(edges[, 2] == max(tips)), c(1, 2)]
      clade_nodes = unique(c(node, clade_nodes[,1], clade_nodes[,2]))
      lineage = ifelse (ext_node %in% clade_nodes, clade, lineage)
    }
    
    tips = setNames(1:length(tree$tip.label), tree$tip.label)
    node = length(tree$tip.label) + 1
    
    ext_time = tree$edge.length[edge_pos]
    
    int_seq = anc.seq[[as.character(int_node)]][attr(anc.seq, 'index'),]
    
    if (ext_node %in% tips) {
      ext_seq = anc.seq[[names(tips[which(tips == ext_node)])]][attr(anc.seq, 'index'),]
    } else {
      ext_seq = anc.seq[[as.character(ext_node)]][attr(anc.seq, 'index'),]
      # ext_seq = anc.seq[[as.character(ext_node)]][attr(anc.seq, 'index')]
    }
    
    row.names(ext_seq) = 1:nrow(ext_seq)
    ext_res_snps = ext_seq[unique(sort(as.vector(unlist(
      sapply(drug_db$snp, function (x)
        strsplit(x, ';')[[1]])
    )))),]
    ext_res_bases = setNames(apply(ext_res_snps, 1, function(x)
      ifelse (length(which(
        x == max(x)
      )) == 1, bases[which(x == max(x))], 'N')),
      row.names(ext_res_snps))
    ext_res_present = drug_db[apply(drug_db, 1,
                                    function(x)
                                      all(ext_res_bases[strsplit(x['snp'], ';')[[1]]] == strsplit(x['alt'], ';')[[1]]) == TRUE),]
    
    
    row.names(int_seq) = 1:nrow(int_seq)
    int_res_snps = int_seq[unique(sort(as.numeric(as.vector(unlist(
      sapply(drug_db$snp, function (x)
        strsplit(x, ';')[[1]])
    ))))),]
    int_res_bases = setNames(apply(int_res_snps, 1, function(x)
      ifelse (length(which(
        x == max(x)
      )) == 1, bases[which(x == max(x))], 'N')),
      row.names(int_res_snps))
    int_res_present = drug_db[apply(drug_db, 1,
                                    function(x)
                                      all(int_res_bases[strsplit(x['snp'], ';')[[1]]] == strsplit(x['alt'], ';')[[1]]) == TRUE),]
    
    t_to_res[['int_node']] = c(t_to_res[['int_node']], int_node)
    t_to_res[['ext_node']] = c(t_to_res[['ext_node']], ext_node)
    t_to_res[['int_year']] = c(t_to_res[['int_year']], 0)
    t_to_res[['ext_year']] = c(t_to_res[['ext_year']], ext_time)
    t_to_res[['int_ci']] = c(t_to_res[['int_ci']], NA)
    t_to_res[['ext_ci']] = c(t_to_res[['ext_ci']], NA)
    t_to_res[['int_st']] = c(t_to_res[['int_st']], ifelse(nrow(int_res_present) == 0, 0 , 1))
    t_to_res[['ext_st']] = c(t_to_res[['ext_st']], ifelse(nrow(ext_res_present) == 0, 0 , 1))
    t_to_res[['int_pos']] = c(t_to_res[['int_pos']],
                              paste(int_res_present$pos, collapse = ';'))
    t_to_res[['ext_pos']] = c(t_to_res[['ext_pos']],
                              paste(ext_res_present$pos, collapse = ';'))
    t_to_res[['int_drug']] = c(t_to_res[['int_drug']],
                               paste(int_res_present$drug, collapse = ';'))
    t_to_res[['ext_drug']] = c(t_to_res[['ext_drug']],
                               paste(ext_res_present$drug, collapse = ';'))
    t_to_res[['int_gene']] = c(t_to_res[['int_gene']],
                               paste(int_res_present$gene, collapse = ';'))
    t_to_res[['ext_gene']] = c(t_to_res[['ext_gene']],
                               paste(ext_res_present$gene, collapse = ';'))
    t_to_res[['int_mut']] = c(t_to_res[['int_mut']],
                              paste(int_res_present$mutation, collapse = ';'))
    t_to_res[['ext_mut']] = c(t_to_res[['ext_mut']],
                              paste(ext_res_present$mutation, collapse = ';'))
    t_to_res[['int_conf']] = c(t_to_res[['int_conf']],
                               paste(int_res_present$confidence, collapse = ';'))
    t_to_res[['ext_conf']] = c(t_to_res[['ext_conf']],
                               paste(ext_res_present$confidence, collapse = ';'))
    t_to_res[['lineage']] = c(t_to_res[['lineage']], lineage)
  }
  t_to_res = as.data.frame(t_to_res)
  t_to_res[t_to_res == ''] = NA
  return (t_to_res)
}


getCoxAxis = function(x) {
  n = 1
  while(TRUE) {
    if (1 %in% pretty(x, n = n)){
      return(pretty(x, n = n))
    } else {
      n = n + 1
    }
  }
}



prettySurvPlot = function (sfit,
                           fit.coxph,
                           xlevels = FALSE,
                           kap_meier = TRUE,
                           haz = TRUE,
                           risk = TRUE,
                           conf.int = TRUE,
                           km.ylim = c(0,1.01),
                           km.xlim = FALSE,
                           km.xlab = 'Time',
                           km.ylab = 'Survival Probability',
                           time.by = FALSE,
                           cols = RColorBrewer::brewer.pal(12, "Set3")[c(1,3,4,5,6,10,12)],
                           col.levels = FALSE,
                           km.size = 0.5, risk.size = 0.2, haz.size = 0.3)
{
  
  if (all(xlevels == FALSE)) {
    xlevels = unname(unlist(fit.coxph$xlevels))
  }

  if (all(c(haz, risk, kap_meier))) {
    par(mar=c(0,0,0,0), omi=(c(0,0,0,0)+0.2))
    # layout(matrix(c(1,2,3), nrow=3, byrow=F), heights=c(0.35, 0.33, 0.32), widths=c(0.8,0.8,0.8))
    layout(matrix(c(1,2,3), nrow=3, byrow=F), heights=c(haz.size, km.size, risk.size), widths=c(0.8,0.8,0.8))
    
  } else if (risk == F & (kap_meier == T & haz == T)) {
    par(mar=c(0,0,0,0), omi=(c(0,0,0,0)+0.2))
    layout(matrix(c(1,2), nrow=2, byrow=F),
           heights=c(0.4, 0.6), widths=c(0.8,0.8))
  } else if (risk == T & (kap_meier == T & haz == F)) {
    par(mar=c(0,0,0,0), omi=(c(0,0,0,0)+0.2))
    layout(matrix(c(1,2), nrow=2, byrow=F),
           heights=c(0.6, 0.4), widths=c(0.8,0.8))
  } else {
    par(mar = c(8,9,2,12))
  }
  
  if (isFALSE(km.xlim)){
    km.xlim = c(0, max(sfit$time))
  }
  ## forest plot 
  
  if (haz) {
    forest_data = summary(fit.coxph)
    x.mean = forest_data$conf.int[, 1]
    x.upper = forest_data$conf.int[, 4]
    x.lower = forest_data$conf.int[, 3]
    pval = forest_data$coefficients[, 5]
    
    pval.plot = c(
      ' ' = Inf,
      '*' = 0.05,
      '**' = 0.01,
      '***' = 0.001
    )
    
    pval.symb = sapply(pval, function(x)
      names(which(pval.plot >= x)[length(which(pval.plot >= x))]))
    
    floor_dec <-
      function(x, level = 1)
        round(x - 5 * 10 ^ (-level - 1), level)
    ceiling_dec <-
      function(x, level = 1)
        round(x + 5 * 10 ^ (-level - 1), level)
    
    x = seq(min(c(0.5, floor(x.lower))), max(c(1, ceiling_dec(x.upper, 1))), 0.5)
    y.names = unname(unlist(fit.coxph$xlevels))
    y = length(y.names)
    
    par(mar = c(3, 10, 2, 10), xpd = F)
    
    plot(
      x,
      rep(-10, length(x)),
      type = "p",
      ylab = "",
      xlab = " ",
      cex = 1.5,
      ylim = c(0, y + 0.5),
      xlim = c(0, max(x) + 0.5),
      lwd = 2,
      pch = 5,
      axes = F,
      main = " "
    )
    

    axis(
      1,
      cex.axis = 1.5,
      at = pretty(x),
      labels = pretty(x)
    )
    
    
    # axis(
    #   1,
    #   cex.axis = 1.5,
    #   at = getCoxAxis(x),
    #   labels = getCoxAxis(x)
    # )
    axis(
      2,
      las = 2,
      line = 1,
      at = seq(1, y, 1),
      cex.axis = 1.5,
      labels = rev(y.names),
      tick = FALSE
    )
    
    
    abline(v = 1, lty = 2)
    
    if (all(col.levels != FALSE)) {
      names(cols) = col.levels
      cols = cols[xlevels]
    }
    
    points(
      1,
      y,
      cex = 2.5,
      lwd = 2.5,
      pch = 19,
      col = cols[1]
    )
    points(
      x.mean,
      rev(seq(1, y - 1, 1)),
      cex = 2.5,
      lwd = 2,
      pch = 19,
      col = cols[2:length(cols)]
    )
    
    segments(x.lower,
             rev(seq(1, y - 1, 1)),
             x.upper,
             rev(seq(1, y - 1, 1)),
             lwd = 2.5,
             col = cols[2:length(cols)])
    
    arrows(
      x.lower,
      rev(seq(1, y - 1, 1)),
      x.upper,
      rev(seq(1, y - 1, 1)),
      lwd = 2.5,
      angle = 90,
      code = 3,
      length = 0.05,
      col = cols[2:length(cols)]
    )
    
    mtext(text = "Hazard ratio",
          1,
          line = 3,
          cex = 1.2, font=2)
    par(xpd = T)
    text(
      x = par('usr')[2],
      y = rev(seq(1, y - 1, 1)),
      pval.symb,
      cex = 2,
      font = 2,
      col = cols[2:length(cols)]
    )
  }
  
  
  ## Kaplan-meier plot
  
  if (isTRUE(kap_meier)){
  if (risk == T){
    par(mar = c(1,10,2,10), xpd = F)
  }
  
  par(xpd = F)
  surv.data = summary(sfit)
  
  lin_names = levels(surv.data$strata)
  lin_names = sapply(lin_names, function (x)
    strsplit(x, '=')[[1]][length( strsplit(x, '=')[[1]])])
  
  survplot = plot(sfit,
                  ylim = km.ylim,
                  xlim = km.xlim,
                  conf.int = F,
                  mark.time = T,
                  col = cols,
                  bty = 'n',
                  xaxt = "n", yaxt = "n"
  )
  
  axis (2, las = 2, cex.axis = 1.5, line = 1)
  mtext(text = km.ylab, 2, line = 5, cex = 1.2, font=2)

    if (!risk) {
    axis (1, cex.axis = 1.2, line = 1, xlim = km.xlim)
    mtext(text = km.xlab, 1, line = 4, cex = 1.2)
  }
  
  if (conf.int) {
    init.val = 1
    final.val = 0
    for (i in 1:length(sfit$strata)){
    with(sfit, polygon(
      c(rev(time[init.val:(final.val+strata[i])]),
        time[init.val:(final.val+strata[i])]),
      c(rev(lower[init.val:(final.val+strata[i])]),
        upper[init.val:(final.val + strata[i])]),
      col = adjustcolor(cols[i], 0.3),
      border = F
    ))
    init.val = init.val + sfit$strata[i]
    final.val = final.val + sfit$strata[i]
    
  }
  }
  
  
  
  
  
  par(xpd = T)
  text(x = km.xlim[1] + 10, y = km.ylim[1] + 0.1,
       survminer::surv_pvalue(sfit)$pval.txt, cex = 1.5)
  
  
  legend(par('usr')[2], par('usr')[4],
         legend = lin_names, 
         lty = c(1, 1), lwd = c(7, 7),
         col =  cols,
         bty = "n", x.intersp = 1, cex = 1)
  }
  
  if (!time.by){
    time.by = 1
  }

  if (risk) {
    nrisk = survutils::get_nrisk_tbl(sfit, timeby = time.by)
    
    x = sfit$time
    y.names = unname(unlist(fit.coxph$xlevels))
    y = length(y.names)
    
    
    par(mar = c(5, 10, 0, 10), xpd = F)
    
    
    plot(
      x,
      rep(-10, length(x)),
      type = "p",
      ylab = "",
      xlab = " ",
      cex = 1.5,
      ylim = c(0, y + 0.5),
      xlim = km.xlim,
      lwd = 2,
      pch = 5,
      axes = F,
      main = " "
    )
    
    x_ax = axis (1,
                 cex.axis = 1.5,
                 line = 1,
                 xlim = km.xlim)
    
    mtext(text = km.xlab,
          1,
          line = 4,
          cex = 1.2, font=2)
    
    mtext(
      text = "Numbers\nat risk",
      side = 4,
      las = 1,
      line = 2,
      cex = 1, font=2
    )
    
    ylevels = levels(nrisk$strata)
    xlevels = x_ax
    
    if (length(xlevels[!xlevels %in% nrisk$time]) != 0) {
      missing_times = xlevels[!xlevels %in% nrisk$time]
      matrix_missing = c()
      for (miss_time in missing_times){
        matrix_missing = rbind(matrix_missing, matrix(c(c(ylevels),
                 rep(miss_time, length(ylevels)), 
                 rep(0, length(ylevels))),
               nrow = length(ylevels), ncol = 3))
      }
      colnames(matrix_missing) = names(nrisk)
      nrisk = rbind(nrisk, matrix_missing)
      
    }
    
    
    par(xpd = T)
    text(
      x = xlevels,
      y = y,
      nrisk[nrisk$strata == ylevels[1] &
              nrisk$time %in% xlevels, 'n.risk'],
      cex = 1.5,
      col = cols[1],
      font = 2
    )
    
    
    text(
      x = rep(xlevels, length(ylevels) - 1),
      y = rev(rep(seq(1, y - 1, 1), each = length(xlevels))),
      nrisk[nrisk$strata %in% ylevels[2:length(ylevels)] &
              nrisk$time %in% xlevels, 'n.risk'],
      cex = 1.5,
      col = rep(cols[2:length(cols)], each = length(xlevels)),
      font = 2
    )
    
  }
}





coxzph = function (fit, resids = TRUE, se = TRUE, df = 4, nsmo = 40, var, 
          point.col = "steelblue", point.size = 1, point.shape = 19, point.alpha = 0.7, 
          caption = NULL, global.p = FALSE) 
{
  x <- fit
  if (!methods::is(x, "cox.zph")) 
    stop("Can't handle an object of class ", class(x))

  xx <- x$x
  yy <- x$y
  d <- nrow(yy)
  df <- max(df)
  nvar <- ncol(yy)
  pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
  temp <- c(pred.x, xx)
  lmat <- splines::ns(temp, df = df, intercept = TRUE)
  pmat <- lmat[1:nsmo, ]
  xmat <- lmat[-(1:nsmo), ]
  qmat <- qr(xmat)
  if (qmat$rank < df) {
    stop("Spline fit is singular, try a smaller degrees of freedom")
  }
  if (se) {
    bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
    xtx <- bk %*% t(bk)
    seval <- d * ((pmat %*% xtx) * pmat) %*% rep(1, df)
  }
  ylab <- paste("Beta (t) for\n", dimnames(yy)[[2]])
  if (missing(var)) {
    var <- 1:nvar
  } else {
    if (is.character(var)) {var <- match(var, dimnames(yy)[[2]])}
    if (any(is.na(var)) || max(var) > nvar || min(var) < 
        1) {stop("Invalid variable requested")}
  }
  
  if (x$transform == "log") {
    xx <- exp(xx)
    pred.x <- exp(pred.x)
  } else if (x$transform != "identity") {
    xtime <- as.numeric(dimnames(yy)[[1]])
    indx <- !duplicated(xx)
    apr1 <- approx(xx[indx], xtime[indx], seq(min(xx), max(xx), 
                                              length = 17)[2 * (1:8)])
    temp <- signif(apr1$y, 2)
    apr2 <- approx(xtime[indx], xx[indx], temp)
    xaxisval <- apr2$y
    xaxislab <- rep("", 8)
    for (i in 1:8) xaxislab[i] <- format(temp[i])
  }
  
  reg.line = FALSE
  # par(mfrow=c(length(var),1),
  #     oma = c(2,2,2,2))
  for (i in var) {
    invisible(pval <- round(x$table[i, 3], 3))
    main.title = paste0("Schoenfeld individual Test\n p: ",
                        pval)
    y <- yy[, i]
    yhat <- as.vector(pmat %*% qr.coef(qmat, y))
    if (resids) {
      yr <- range(yhat, y)
    } else {yr <- range(yhat)}
    if (se) {
      temp <- as.vector(2 * sqrt(x$var[i, i] * seval))
      yup <- yhat + temp
      ylow <- yhat - temp
      yr <- range(yr, yup, ylow)
    }
    if (x$transform == "identity") {
      reg.line = lm(pred.x ~ yhat)
    }
    else if (x$transform == "log") {
      reg.line = lm(log(pred.x) ~ yhat)
      x.lab = "Time"
      y.lab = ylab[i]
      y.lim = yr
    }
    else {
      reg.line = lm(pred.x ~ yhat)
      x.lab = "Time"
      y.lab = ylab[i]
      x.breaks = xaxisval
      x.labs = xaxislab
      y.lim = yr
    }
    if (resids) {
      par(mar = c(4, 8, 2, 2), xpd = F)
      plot(xx, y,
           ylab = '', xlab = '', main = '', type = 'n',
           axes = FALSE, ylim = yr, xlim = c(min(xx), max(xx)))
      axis(1, at = xaxisval, labels = xaxislab)
      mtext(1, text = "Time", line = 3, cex = 1)
      axis(2, line = 2, las = 2)
      mtext(2, text = ylab[i], line = 5, cex = 1)
    }
    if (se) {
      polygon(x = c(pred.x, rev(pred.x)), 
              y = c(yup, rev(ylow)), 
              border = F,
              col = adjustcolor('lightgrey', 0.3))
      # lines(x = pred.x, y = yup, lty = 2)
      # lines(x = pred.x, y = ylow, lty = 2)
    }
    points(x = xx,
           y = y,
           pch = point.shape,
           cex = point.size,
           col = adjustcolor(point.col, point.alpha))
    if (!isFALSE(reg.line)){
      abline(reg.line, 
             lwd = 2)
    }
    par(xpd = NA)
    text(x = (par('usr')[1] + par('usr')[2]) / 2,
         y = par('usr')[4],
         label = main.title, cex = 1.5)
  }

  if ("GLOBAL" %in% rownames(x$table) & !isFALSE(global.p))
    global_p <- round(x$table["GLOBAL", 3], 3)
  else global_p <- FALSE
  if (!isFALSE(global_p)){
    mtext (side = 3, line = 1, text = paste0("Schoenfeld global Test p: ",
                                   global_p), cex = 1.5)
  }

}


