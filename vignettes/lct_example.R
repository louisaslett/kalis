## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----set_up, eval=FALSE-------------------------------------------------------
#  nthreads <- as.integer(4)
#  data_dir <- "./"

## ----libraries, eval=FALSE----------------------------------------------------
#  require("kalis")
#  require("R.utils")
#  require("data.table")
#  require("fastcluster")
#  require("Matrix")
#  require("viridisLite")
#  require("kgp")

## ----declare_run_parameters, eval=FALSE---------------------------------------
#  # Declare LS Model Parameters
#  #########################################
#  neg_log10_Ne <- 10
#  neg_log10_mu <- 4
#  
#  # Declare Target Locus
#  #########################################
#  gene <- "lct"
#  gene_target_pos <- 136608646 # rs4988235 in hg19 coordinates
#  pos <- fread(paste0(data_dir, gene, ".legend.gz"))$position
#  target_idx <- match(TRUE, pos >= gene_target_pos)
#  
#  # run kalis
#  #########################################
#  
#  CacheHaplotypes(haps = paste0(data_dir, gene, ".hap.gz"))
#  
#  diff_map <- diff(fread(paste0(data_dir, gene, ".map"))[[3]])
#  pars <- Parameters(rho = CalcRho(diff_map, s = 10^-neg_log10_Ne), mu = 10^-neg_log10_mu)
#  fwd <- MakeForwardTable(pars)
#  bck <- MakeBackwardTable(pars)
#  
#  Forward(fwd, pars, target_idx, nthreads = nthreads)
#  Backward(bck, pars, target_idx, nthreads = nthreads)
#  
#  M <- DistMat(fwd, bck, type = "raw", nthreads = nthreads)

## ----helper_functions, eval=FALSE---------------------------------------------
#  plot_mat <- function(x, file, raster = TRUE, rel_scale = TRUE) {
#  
#    temp_col_scale <- rev(viridisLite::viridis(100))
#  
#    if(!rel_scale){
#      mx <- ceiling(max(x, na.rm = TRUE))
#      if(mx > 100) { stop("the max entry of x cannot exceed 100 for this plot's color scale") }
#      temp_col_scale <- temp_col_scale[1:mx]
#    }
#  
#    cairo_pdf(file)
#    print(lattice::levelplot(x[, ncol(x):1],
#                             useRaster = raster,
#                             col.regions = grDevices::colorRampPalette(temp_col_scale)(100),
#                             yaxt = "n", xaxt = "n", xlab = "", ylab = "", xaxt = "n"))
#    dev.off()
#  }
#  
#  interp_hapmap <- function(path,bp){
#    d <- data.table::fread(path)
#    approx(d$`Position(bp)`, d$`Map(cM)`, xout = bp, method = "linear", rule = 2)$y
#  }
#  
#  dip2hapidx <- function(x){
#    x <- 2*x
#    c(rbind(x-1, x))
#  }

## ----make_plot,eval=FALSE-----------------------------------------------------
#  # Load sample population information
#  #########################################
#  cluster_by <- "isAFR"
#  id <- fread(paste0(data_dir, gene, ".sample"))$sample
#  init_order_samples <- order(id)
#  samples <- merge(data.table("id" = id), kgp3, by = "id")
#  if(nrow(samples) != length(id)) { stop("some samples have been removed by merging with kgp3") }
#  if(!all.equal(init_order_samples, order(samples$id))) { stop("some samples have been moved out of the order in lct.sample") }
#  samples[,isAFR := ifelse(reg == "AFR", "AFR", "not_AFR")]
#  
#  # Symmeterize & Scale Distance Matrix at LCT
#  ###############################################
#  M <- (0.5/(neg_log10_mu*log(10))) * (M + t(M))
#  
#  # Permute & Cluster Distance Matrix
#  ###################################################################
#  diploid_perm <- order(samples$reg, samples$pop, samples$id)
#  psamples <- samples[diploid_perm,]
#  
#  haploid_perm <- dip2hapidx(diploid_perm)
#  
#  pM <- M[, haploid_perm][haploid_perm,]
#  
#  hap_groups <- table(psamples[[cluster_by]])
#  hap_groups <- hap_groups[unique(psamples[[cluster_by]])]
#  
#  baseline_idx <- c(0, cumsum(2*hap_groups)[-length(hap_groups)])
#  names(baseline_idx) <- names(hap_groups)
#  
#  order_M <- as.list(hap_groups)
#  names(order_M) <- names(hap_groups)
#  
#  for(i in 1:length(hap_groups)){
#    current_pop_samples <- which(psamples[[cluster_by]] == names(hap_groups)[i])
#    current_pop_haplotypes <- dip2hapidx(current_pop_samples)
#    sM <- pM[current_pop_haplotypes, current_pop_haplotypes]
#    order_M[[names(hap_groups)[i]]] <- baseline_idx[names(hap_groups)[i]] + fastcluster::hclust(as.dist(sM), method="average")$order
#  }
#  
#  order_M <- unlist(order_M)
#  cM <- pM[, order_M][order_M,]
#  
#  
#  # Plot clustered Distance Matrix
#  #########################################
#  plot_mat(cM, paste0(data_dir, gene, "_dist_mat.pdf"))

