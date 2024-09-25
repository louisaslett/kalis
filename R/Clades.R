


get_neigh <- function(x,i){
  idx <- x[[1]][c(i,i+1L)]
  x[[2]][seq.int(idx[1],idx[2]-1L)]
}

get_neigh_seq <- function(x, i, return.lengths = FALSE){
  from <- x[[1]][i]
  nvec <- x[[1]][i+1] - from
  if(return.lengths){
    list("seq" = x[[2]][sequence(nvec,from)],
         "lengths" = nvec)
  } else {
    x[[2]][sequence(nvec,from)]
  }
}

#' Sprigs
#' @export Sprigs
Sprigs <- function(x, old.sprigs = FALSE){

  N <- length(x[[1]])-1L
  roster <- rep(NA_integer_,N)
  label <- 0L
  done <- neighborhood.is.sprig.ind <- rep(FALSE,N)

  if(old.sprigs){

    xx <- as.list(1:N)
    for(i in 1:N){xx[[i]] <- get_neigh(x,i)}

    roster <- Sprigs_old(xx, use.forking = FALSE, nthreads = 1L, add.self = FALSE)

    label <- attr(roster,"n.sprigs")

    attributes(roster) <- NULL

    neighborhood.is.sprig.ind <- !is.na(roster)


    # for(i in sample.int(N)){
    #
    #   if(done[i]){next}
    #
    #   C = get_neigh(x,i)
    #   lC <- length(C)
    #
    #   neigh.list <- get_neigh_seq(x, C, return.lengths = TRUE)
    #
    #   temp.table <- table(neigh.list[[1]])
    #   proposed.set <- as.integer(names(temp.table)[which(temp.table == lC)])
    #
    #   # in case the neighborhood of i overshoots into previously established cliques, it has a BIG effect in real data
    #   proposed.set <- proposed.set[!done[proposed.set]]
    #
    #   if(length(proposed.set) > 1 && i %in% proposed.set){
    #
    #     label <- label + 1L
    #     # this repeated intersection step has truly a small effect but
    #     # helps guard us against the case where i might have erroneously added some candidates in its neighborhood that
    #     # do not belong in the clade and do not include some clade members.  This steps helps us recover those clade members
    #
    #     neigh.list <- get_neigh_seq(x, proposed.set, return.lengths = TRUE)
    #     temp.table <- table(neigh.list[[1]])
    #     proposed.set <- as.integer(names(temp.table)[which(temp.table == lC)])
    #     proposed.set <- proposed.set[!done[proposed.set]]
    #
    #     roster[proposed.set] <- label
    #     done[proposed.set] <- TRUE
    #     neighborhood.is.sprig.ind[proposed.set] <- neigh.list[[2]] == lC # indicator that a haplotype's neighborhood is exactly the proposed sprig
    #   }
    # }


  } else {

    # the randomness in indices here is not really essential but
    for(i in seq_len(N)){ #sample.int(N)){

      if(done[i]){next}

      C = get_neigh(x,i) # keep simple get_neigh here because about 2x faster than get_neigh_seq when only querying one index
      C <- C[!done[C]] # we know this will at least include i because of the if(!done[i]) check above and every neighborhood includes self

      lC <- length(C)

      if(lC == 1){ # C is an orphan that could never be a part of another clade because all of it's neighbors are already assigned
        done[i] <- TRUE
        next}

      neigh.list <- get_neigh_seq(x, C, return.lengths = TRUE)

      if(all(tabulate(factor(neigh.list[[1]],C),nbins = lC)==lC)){
        label <- label + 1L
        roster[C] <- label
        done[C] <- TRUE
        neighborhood.is.sprig.ind[C] <- neigh.list[[2]] == lC # indicator that a haplotype's neighborhood is exactly the proposed sprig
      }
    }

  }

  list("assignments" = roster,
       "to.prune" = neighborhood.is.sprig.ind,
       "num.sprigs" = label)
}


#' @importFrom utils getFromNamespace
UpdateMatrixInPlace <- function(M,row.idx,col.idx,x){
  invisible(.Call(getFromNamespace("CCall_UpdateRealInPlace","kalis"),M,
                  as.integer(row.idx + (col.idx-1L)*nrow(M)),x))
}

# test <- matrix(as.double(1:144),12,12)
# UpdateMatrixInPlace(test,c(5,12,12),c(1,3,5),as.double(c(100,200,300)))

#' PruneCladeMat
#' @export PruneCladeMat
PruneCladeMat <- function(M, neigh, sprigs, prune = "singleton.info", from.recipient = 1L){

  if(!from.recipient%%2){stop("from.recipient must be odd and encode the index of the first recipient haplotype")}

  N.recipients <- 2 * ncol(M)

  if(prune=="singleton.info"){

    v <- neigh[[2]][[2]] - neigh[[2]][[1]]
    v <- v[seq.int(1,N.recipients,2)] + v[seq.int(2,N.recipients,2)]
    UpdateMatrixInPlace(M,
                        seq.int(from = (from.recipient+1L)/2L,length.out = ncol(M)),
                        seq.int(from = 1, to = ncol(M)),
                        v)

  } else if(prune=="sprigs"){

    neigh.list <- get_neigh_seq(neigh[[1]], which(sprigs$to.prune), return.lengths = TRUE)

    hap.idx <- cbind(neigh.list[[1]], rep(which(sprigs$to.prune), times = neigh.list[[2]]))

    key <- rep(0L,nrow(hap.idx))

    hap.idx.odd <- hap.idx%%2
    temp.hap.idx.odd <- hap.idx.odd[,1] + hap.idx.odd[,2]
    key[temp.hap.idx.odd==2] <- 1L
    key[temp.hap.idx.odd==0] <- 4L
    temp.hap.idx.odd <- hap.idx.odd[,1] - hap.idx.odd[,2]
    key[temp.hap.idx.odd==1] <- 3L
    key[temp.hap.idx.odd==-1] <- 2L

    sim.updates <- rep((neigh[[2]][[3]]-neigh[[2]][[2]])[sprigs$to.prune], times = neigh.list[[2]])

    # if M was not already dediploided:
    # M[hap.idx] <- M[hap.idx] - sim.updates

    # since M is dediploided, we run

    if(!is.na(match(1L,key))){
      to.fetch <- key==1L
      UpdateMatrixInPlace(M,
                          (hap.idx[to.fetch,1L]+1L)/2L,
                          (hap.idx[to.fetch,2L]+1L)/2L,
                          sim.updates[to.fetch])}

    if(!is.na(match(2L,key))){
      to.fetch <- key==2L
      UpdateMatrixInPlace(M,
                          hap.idx[to.fetch,1L]/2L,
                          (hap.idx[to.fetch,2L]+1L)/2L,
                          sim.updates[to.fetch])}

    if(!is.na(match(3L,key))){
      to.fetch <- key==3L
      UpdateMatrixInPlace(M,
                          (hap.idx[to.fetch,1L]+1L)/2L,
                          (hap.idx[to.fetch,2L])/2L,
                          sim.updates[to.fetch])}

    if(!is.na(match(4L,key))){
      to.fetch <- key==4L
      UpdateMatrixInPlace(M,
                          hap.idx[to.fetch,1L]/2L,
                          hap.idx[to.fetch,2L]/2L,
                          sim.updates[to.fetch])}

  } else {
    stop("unrecognized option for prune")
  }
  invisible(NULL)
}






#' Probabilistic Clades
#'
#' Utility for calling probabilistic clades at, in between, or excluding variants.
#' @param fwd a forward table as returned by [MakeForwardTable()]
#' @param bck a backward table as returned by [MakeBackwardTable()]
#' @param pars a `kalisParameters` object, as returned by [Parameters()].
#' @param beta.theta.opts a list; see Details for [DistMat()].
#' @param safety.checks a logical, should safety checks be applied to the distances?  See [DistMat()].
#' @param neighbors a logical, should nearest neighbors be pre-calculated?  See [Neighbors()].
#' @param use.forking a logical, should forked processes be used?
#' @param nthreads the number of CPU cores to use. Currently, no parallelism is used.
#' @return
#'   a `kalisClades` object encoding probabilistic clade calls
#'
#' @importFrom data.table frank
#' @export Clades
Clades <- function(fwd, bck, pars, beta.theta.opts = NULL,
                   safety.checks = FALSE, neighbors = FALSE,
                   #use.bettermc = FALSE,
                   use.forking = FALSE,
                   forking.chunk.size = 100L,
                   mc.preschedule = FALSE, # FALSE is more conservative of memory but means many new forked processes need to be launched so it's slower than TRUE
                   nthreads = 1L){
  # currently only outputs a list but should eventually also output a matrix of integers and an attribute list of clades

  unit.mut.dist <- -log(pars$pars$mu)

  M <- DistMat(fwd, bck, beta.theta.opts = beta.theta.opts, nthreads = nthreads)

  if(safety.checks){
    M[!is.finite(M)] <- 0
    diag(M) <- NA_real_
  }

  rank_donors_func <- function(x, type="linear_20", neighbors = FALSE, mac.range = c(NA,NA)){
    rank_donors_func_res <- as.list(1:length(x))
    for(j in 1:length(x)){
      d.ranks <- data.table::frank(M[,x[j]], na.last = FALSE, ties.method = "first")
      phi <- c(diff(M[order(d.ranks),x[j]]),0)
      phi[1] <- 0
      phi <- phi / unit.mut.dist # an N-long vector
      if(type == "linear_20"){
        phi[phi > 1] <- 1
        phi[phi < 0.2] <- 0
      } else if(type == "step_80"){
        phi[phi < 0.8] <- 0
        phi[phi > 0] <- 1
      }

      if(!is.na(mac.range[1])){phi[1:mac.range[1]] <- 0}
      if(!is.na(mac.range[2])){phi[mac.range[2]:nrow(fwd$alpha)] <- 0}

      i <- which(phi!=0)

      # compress phi
      clades <- cbind(i,phi[i]) # if i = integer(0) (no clades called), clades will be a 0 x 2 matrix.
      attr(d.ranks,"clades") <- clades

      if(neighbors){
        attr(d.ranks,"neighbors") <- if(nrow(clades)){
          match(2:clades[1,1],d.ranks)
        } else {
          NA_integer_
        }
      }

      rank_donors_func_res[[j]] <- d.ranks
    }
    rank_donors_func_res
  }


  chunks <- chunk_int(ncol(M), chunk.size = forking.chunk.size)

  if(use.forking){
    # if(use.bettermc){
    #   rank.list <- bettermc::mclapply(chunks, rank_donors_func, neighbors = neighbors, mc.preschedule = mc.preschedule, mc.cores=nthreads, mc.share.copy = FALSE)
    # } else {
      rank.list <- parallel::mclapply(chunks, rank_donors_func, neighbors = neighbors, mc.preschedule = mc.preschedule, mc.cores=nthreads)
    #}
  } else {
    rank.list <- lapply(chunks, rank_donors_func, neighbors = neighbors) # this matrix is ranked in each column, not scaled by Ne or Mu
  }

  rank.list <- unlist(rank.list,recursive = FALSE)

  attr(rank.list,"from_recipient") <- fwd$from_recipient
  attr(rank.list,"to_recipient")   <- fwd$to_recipient

  class(rank.list) <- c("kalisClades","list") # rank.list is a list where each element is a vector of ranks with attributes clades

  rank.list
}


#' Neighbors
#'
#' Utility for calling tied nearest neighbors for each recipient haplotype
#' @param x a `kalisClades` object returned by [Clades()]
#' @param use.forking a logical, should forked processes be used?
#' @param nthreads the number of CPU cores to use. Currently, no parallelism is used.
#' @return
#'   a `kalisNeighbors` encoding the nearest neighbors for each recipient haplotype
#'
#' @export Neighbors
Neighbors <- function(x,
                      #use.bettermc = FALSE,
                      use.forking = FALSE, nthreads = 1L){
  # currently only supports list x but should support matrix x as well

  if(!is.null(attr(x[[1]],"neighbors"))){

    neighbors <- lapply(x,function(z){attr(z,"neighbors")})

  } else {


    call_neighbors <- function(z){
      # x should be a vector of ranks with attribute "clades"
      clades <- attr(z,"clades")
      if(nrow(clades)){
        match(2:clades[1,1],z)
      } else {
        NA_integer_
      }
    }

    if(use.forking){
      # if(use.bettermc){
      #   neighbors <- bettermc::mclapply(x, call_neighbors, mc.cores = nthreads, mc.share.copy = FALSE)
      # } else {
        neighbors <- parallel::mclapply(x, call_neighbors, mc.cores = nthreads)
      #}
    } else {
      neighbors <- lapply(x,call_neighbors)
    }
  }

  attr(neighbors,"from_recipient") <- attr(x,"from_recipient")
  attr(neighbors,"to_recipient")   <- attr(x,"to_recipient")
  class(neighbors) <- c("kalisNeighbors","list")

  neighbors
}


#' Sprigs
#'
#' Utility for calling sprigs from probabilistic clades
#' @param x a `kalisNeighbors` object returned by [Neighbors()], a `kalisClades` object returned by [Clades()] with `neighbors = TRUE`, or a list
#' @param use.forking a logical, should forked processes be used?
#' @param nthreads the number of CPU cores to use. Currently, no parallelism is used.
#' @return
#'   a `kalisSprigs` object assigning each haplotype to a sprig
#'
#' @export Sprigs_old
Sprigs_old <- function(x, use.forking = FALSE, nthreads = 1L, add.self = TRUE){

  # this version of Sprigs still has a bit of randomness in it's sprig building between runs on the same input
  # which can be seen by running table(is.na(s),is.na(s1)) where s and s1 are the output of Sprigs
  # for the same data run twice. it's relatively minor

  if(inherits(x,"kalisClades")){
    if(!is.null(attr(x[[1]],"neighbors"))){
      x <- lapply(x,function(z){attr(z,"neighbors")})
    } else {
      stop("The kalisClades provided do not have the Neighbors pre-calculated, use kalis::Neighbors to obtain them and then pass them to Sprigs")
    }
  }

  # x here is a list that's N long st x[[i]] gives the indices of the (tied) nearest neighbors of i
  roster <- rep(NA_integer_,length(x))

  label <- 0L
  # add self to own neighborhood
  if(add.self){x <- mapply(c,x,1:length(x))}

  done <- rep(FALSE,length(x))
  to.prune <- rep(NA_integer_,length(x))

  # the randomness in indices here is not really essential but
  for(i in sample.int(length(x))){
    if(!done[i]){

      # pulling out cliques in the graph that are fully connected bi-directionally:
      # if i is in a clique, rather trivially, this will return the full clique
      # Note, we require i %in% proposed.set to prevent called cliques from being broken up later in the for loop:
      # if i is not in a clique, then it's still possible for a partial clique to be returned that doesn't include i if i projects onto
      # a superset or subset of a clique.  If i supercedes a clique member and projects
      # onto a subset, this clique subset will be overwritten later by the larger clique.  However, it would still be possible for a i that comes
      # after all of the clade members in our for loop to break up the clique by projecting onto a subset of them.
      # Enforcing i %in% proposed.set avoids that possibility.

      # we also require that length(proposed.set) > 1 so that we don't end up with solo cliques being called that are just i by itself.

      #missing_sprig_6 <- c(6103,1804, 6015, 4726, 4752, 807,3118,3991,6466,6068,  10,1250, 3669, 3658, 1997, 1399, 1116, 3738, 5015)
      proposed.set <- Reduce(intersect,x[x[[i]]])
      # we can really speed up this Reduce by using table and looking for entries that are present in all groups

      # in case the neighborhood of i overshoots into previously established cliques, it has a BIG effect in real data
      proposed.set <- proposed.set[!done[proposed.set]]

      if(length(proposed.set) > 1 && i %in% proposed.set){

        label <- label + 1L
        # this repeated intersection step has truly a small effect but
        # helps guard us against the case where i might have erroneously added some candidates in its neighborhood that
        # do not belong in the clade and do not include some clade members.  This steps helps us recover those clade members
        proposed.set <- Reduce(intersect,x[proposed.set])
        proposed.set <- proposed.set[!done[proposed.set]]

        # if(!all(is.na(roster[missing_sprig_6])) && !all(roster[missing_sprig_6]==6L)){
        #   print(i)
        #   print(label)
        #   print(roster[missing_sprig_6])
        #   browser()
        # }
        roster[proposed.set] <- label
        done[proposed.set] <- TRUE
      }
    }

    # individuals that are not part of a fully connected clique are left with NA_integer_ on the roster
  }

  # Size frequency spectrum: table(table(roster))

  attr(roster,"n.sprigs") <- label
  attr(roster,"from_recipient") <- attr(x,"from_recipient")
  attr(roster,"to_recipient")   <- attr(x,"to_recipient")
  class(roster) <- c("kalisSprigs","integer")

  roster
}

#Testing Sprigs
# kalis::Sprigs(list(
#   1:2,
#   1:2,
#   3:7,
#   1:10,
#   1:10,
#   1:10,
#   5:11
# ))

#' CladeMat OLD
#'
#' Utility for contructing a probabilistic clade matrix
#' @param x a `kalisClades` object returned by [Clades()]
#' @param ploidy an integer, the ploidy of the organism
#' @param sprigs.to.prune a `kalisSprigs` object returned by [Sprigs()] encoding sprigs that should be excluded from the matrix returned
#' @param assemble a logical, if `FALSE` return the clade matrix as a list of columns rather than as a symmetrized matrix
#' @param use.forking a logical, should forked processes be used?
#' @param nthreads the number of CPU cores to use. Currently, no parallelism is used.
#' @return
#'   a matrix representation of the probabilistic clades provided
#'
#' @export CladeMat_old
CladeMat_old <- function(x, ploidy = 2L, sprigs.to.prune = NULL, assemble = TRUE,
                         #use.bettermc = FALSE,
                         use.forking = FALSE, forking.chunk.size = 100L, mc.preschedule = FALSE, nthreads = 1L){

  # prepare sprigs
  if(is.null(sprigs.to.prune)){sprigs.to.prune <- integer()}
  sl <- length(sprigs.to.prune)
  if(sl){sprig.sizes <- tabulate(sprigs.to.prune)}

  n.recipient.samples <- as.integer(length(x)/ploidy)

  chunks <- chunk_int(n.recipient.samples, chunk.size = forking.chunk.size)

  if(ploidy == 1){
    omega_func <- function(s){
      omega_func_res <- as.list(1:length(s))
      for(j in 1:length(s)){

        N <- length(x[[s[j]]])

        idx <- attr(x[[s[j]]],"clades")[,1]
        phi <- attr(x[[s[j]]],"clades")[,2]

        # prune sprig
        if(sl && !is.na(sprigs.to.prune[s[j]]) && length(idx) && sprig.sizes[sprigs.to.prune[s[j]]] == idx[1]){
          idx <- idx[-1]
          phi <- phi[-1]
        }

        # we know that phi[N] = 0, so there must always be a 0 appended
        omega_func_res[[j]] <-  inverse.rle(list("values" =  c(rev(cumsum(rev(phi/idx))),0),
                                                 "lengths" = diff(c(0,idx,N))))[x[[s[j]]]]

      }
      omega_func_res
    }

  } else if(ploidy == 2){

    omega_func <- function(s){
      omega_func_res <- as.list(1:length(s))
      for(j in 1:length(s)){
        N <- length(x[[s[j]]])

        idx <- attr(x[[s[j]*2-1]],"clades")[,1]
        phi <- attr(x[[s[j]*2-1]],"clades")[,2]

        idx2 <- attr(x[[s[j]*2]],"clades")[,1]
        phi2 <- attr(x[[s[j]*2]],"clades")[,2]


        if(sl && !is.na(sprigs.to.prune[s[j]*2-1]) && length(idx) && sprig.sizes[sprigs.to.prune[s[j]*2-1]] == idx[1]){
          idx <- idx[-1]
          phi <- phi[-1]
        }

        if(sl && !is.na(sprigs.to.prune[s[j]*2]) && length(idx2) && sprig.sizes[sprigs.to.prune[s[j]*2]] == idx2[1]){
          idx2 <- idx2[-1]
          phi2 <- phi2[-1]
        }

        # we know that phi[N] = 0, so there must always be a 0 appended
        w <- inverse.rle(list("values" =  c(rev(cumsum(rev(phi/idx))),0),
                              "lengths" = diff(c(0,idx,N))))[x[[s[j]*2-1]]] +
          inverse.rle(list("values" =  c(rev(cumsum(rev(phi2/idx2))),0),
                           "lengths" = diff(c(0,idx2,N))))[x[[s[j]*2]]]

        omega_func_res[[j]] <- w[seq(1,N,by=2)] + w[seq(2,N,by=2)]
      }
      omega_func_res
    }

  } else {
    stop("Relatedness currently only supports ploidy  = 1 or 2")
  }

  # we don't simplify this list to a matrix at this stage to help preserve memory.
  if(use.forking){
    # if(use.bettermc){
    #   res <- bettermc::mclapply(chunks, omega_func, mc.preschedule = mc.preschedule, mc.cores=nthreads, mc.share.copy = FALSE)
    # } else {
      res <- parallel::mclapply(chunks, omega_func, mc.preschedule = mc.preschedule, mc.cores=nthreads)
    #}
  } else {
    res <- lapply(chunks, omega_func)
  }

  res <- unlist(res, recursive = FALSE)

  if(assemble){
    res <- do.call(cbind,res)
    res <- 0.5 * (res + t(res))
  }

  res
}

chunk_int <- function(n, chunk.size = 100){
  # subdivide 1:n into chunks of size at most chunk.size
  if(n < 1){stop("n must be an integer >= 1")}
  interval.starts <- seq(1,n,by=chunk.size)
  interval.ends <- c(interval.starts[-1]-1,n)
  res <- as.list(1:length(interval.starts))
  for(i in 1:length(interval.starts)){
    res[[i]] <- seq.int(interval.starts[i],interval.ends[i])}
  res
}
#
# use.forking <- FALSE
# use.forking <- TRUE
# nthreads <- 8L
#
# require(kalis)
# haps.file <-"~/Dropbox/Benchmarking_StatGen/kalis_benchmarking_tests/benchmark_on_msprime_simulations/msprime_sim_N_100000_L_10000.hdf5"
# CacheHaplotypes(haps = haps.file,loci.idx = 1:1000,hap.idx = 1:24000)#SmallHaps)
# #CacheHaplotypes(SmallHaps)
# pars <- Parameters(rep(1e-2, L() - 1))
# fwd <- MakeForwardTable(pars)
# bck <- MakeBackwardTable(pars)
# Forward(fwd,pars,floor(L()/2),1)
# Backward(bck,pars,floor(L()/2),1)
# #
# start <- proc.time()
# rl2 <- Clades(fwd, bck, pars, neighbors = TRUE, safety.checks = FALSE, use.forking = use.forking, nthreads = nthreads)
# finish <- proc.time() - start
# print(finish)
# #
# sprigs <- Sprigs(rl2)
# start <- proc.time()
# M <- CladeMat(rl2, sprigs.to.prune = sprigs, use.forking = use.forking, nthreads=nthreads)
# finish <- proc.time() - start
# print(finish)
#

# rl<- readRDS("~/Downloads/clades_test.rds")
# all.equal(rl,rl2)

# sprigs <- CallSprigs(rl, use.forking = use.forking, nthreads = nthreads)
#
# #hist(sapply(rl,function(x){nrow(attr(x,"clades"))}),breaks=20)
# #mean(unlist(lapply(rl,function(x){attr(x,"clades")[,2]}))>0.5)
#
# rl <- CladeMat(rl, ploidy = 2L, sprigs.to.prune = sprigs, assemble = FALSE, use.forking = use.forking, nthreads = nthreads)
# rl <- do.call(cbind,rl)
# rl <- 0.5 * (rl + t(rl))
#
# r2 <- -r2
# class(r2) <- c("kalisDistanceMatrix","matrix")
# plot(r2)
#
# M <- DistMat(fwd,bck)
# M <- M + t(M)
#
# perm <- fastcluster::hclust(stats::as.dist(M),method="average")$order
#
# layout(matrix(1:3,1))
# print(lattice::levelplot(M[perm,][,rev(perm)],
#                          useRaster = TRUE,
#                          col.regions = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,name = "BuPu"))(100),
#                          yaxt = "n", xaxt = "n", xlab = "", ylab = "", xaxt = "n"))
# print(lattice::levelplot(r1[perm,][,rev(perm)],
#                          useRaster = TRUE,
#                          col.regions = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,name = "BuPu"))(100),
#                          yaxt = "n", xaxt = "n", xlab = "", ylab = "", xaxt = "n"))
# print(lattice::levelplot(r2[perm,][,rev(perm)],
#                          useRaster = TRUE,
#                          col.regions = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,name = "BuPu"))(100),
#                          yaxt = "n", xaxt = "n", xlab = "", ylab = "", xaxt = "n"))
#

