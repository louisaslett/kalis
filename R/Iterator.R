#' Build an efficient iterator over loci
#'
#' Function factory to create a `kalisForwardIterator` for propagating a forward table iteratively over target variants using a table cache and optimal checkpointing.
#'
#' See example.
#'
#' @references
#' Christ, R.R., Wang, X., Aslett, L.J.M., Steinsaltz, D. and Hall, I. (2024) "Clade Distillation for Genome-wide Association Studies", bioRxiv 2024.09.30.615852. Available at: \doi{10.1101/2024.09.30.615852}.
#'
#' @param pars
#'        a `kalisParameters` object, as returned by [Parameters()].
#' @param ram.ckpts
#'        an integer specifying the number of checkpoints to store in RAM.
#' @param targets
#'        a vector of variants to iterate over (starting with the most downstream target).
#' @param base.fwd.table
#'        a `kalisForwardTable` either at the most upstream target, or if the targets are evenly spaced, one interval upstream of the most upstream target.
#'        If `NULL` (the default), this is interpreted as the prior `Pi`, see [Parameters()].
#' @param disk.ckpts
#'        an integer specifying the number of checkpoints to store on disk.
#' @param disk.dir
#'        a path to a directory where a temporary folder may be made to store checkpoints on disk.
#' @param from_recipient
#'        first recipient haplotype included in the tables of the cache, if creating a partial forward table.
#'        By default all are included from the first recipient haplotype.
#'        Haplotypes are indexed from 1.
#' @param to_recipient
#'        last recipient haplotype included in the tables of the cache, if creating a partial forward table.
#'        By default all are included upto the last recipient haplotype.
#'        Haplotypes are indexed from 1.
#' @param lookup.tables
#'        an optional list as returned by [CalcCheckpointTables()].
#' @param cache
#'        a `kalisCheckpointTable` object, as returned by [CreateForwardTableCache()] or this function.
#'        By default `NULL`, which causes this function to create a new cache.
#' @param save.cache
#'        a logical.
#'        When `TRUE` does not reliquish the table cache upon exhaustion of the iterator.
#'        Defaults to `FALSE`.
#' @param force.unif
#'        a logical, if `TRUE` iterate over targets as if they were uniformly spaced.
#'        WARNING: DO NOT use this in conjunction with the targets method, still experimental.
#'        With `force.unif = TRUE`, the resulting iterator will appear to be targeting the first `length(targets)` variants with all methods, but in fact will be silently iterating over the original targets.
#'
#' @return
#' A function for iterating over the set of target variants.
#' The returned function has prototype:
#'
#' `function(fwd, pars, t, nthreads = 1)`
#'
#' which matches the standard [Forward()] function, but which uses the table cache to speed up propagation to the target variant.
#' See [Forward()] for an explanation of arguments.
#'
#' @seealso
#' [MakeForwardTable()] to create a `kalisForwardTable`;
#' [CreateForwardTableCache()] to create a cache which can be used with this function.
#'
#' @examples
#' \dontrun{
#' data("SmallHaps")
#' CacheHaplotypes(SmallHaps)
#' pars <- Parameters()
#' fwd <- MakeForwardTable(pars)
#' bck <- MakeBackwardTable(pars)
#' Iter <- ForwardIterator(2)
#' for(t in targets(Iter)){
#'   Iter(fwd,pars,t)
#'   Backward(bck,pars,t)
#'   print(paste("Mean Distance at locus",t,"is",mean(DistMat(fwd,bck))))
#' }
#' }
#'
#' @export
ForwardIterator <- function(pars,
                            ram.ckpts = 1L,
                            targets = 1:L(),
                            base.fwd.table = NULL,
                            disk.ckpts = 0,
                            disk.dir = NULL,
                            from_recipient = 1,
                            to_recipient = Inf,
                            lookup.tables = NULL,
                            cache = NULL,
                            save.cache = FALSE,
                            force.unif = FALSE){

  force(force.unif)

  # Sanity checks
  ####################
  ram.ckpts <- as.integer(ram.ckpts)
  if(ram.ckpts <= 0){stop("ram.ckpts must be a positive integer")}


  # Check to ensure that the cache provided can actually be recycled for this problem

  if(!is.null(cache)){
    for(i in 1:length(cache)){
      if(cache[[i]]$from_recipient!=from_recipient){stop("The provided cache must have the same from_recipient as currently requested.")}
      if(cache[[i]]$to_recipient!=min(to_recipient,N())){stop("The provided cache must have the same to_recipient as currently requested.")}
    }
  }



  # Only RAM checkpoints for now
  ##################
  # for now we ignore disk checkpoints:
  num.available.ckpts <- ram.ckpts

  if(disk.ckpts != 0 | !is.null(disk.dir)){
    warning("disk checkpoints not yet implemnted, proceeding ignoring disk.ckpts and disk.dir")
  }


  # Cover case when base.fwd.table provided
  ##################

  if(is.null(base.fwd.table)){
    use.pi <- TRUE
  }else{
    if( !("kalisForwardTable" %in% class(base.fwd.table)) ){stop("base.fwd.table is not a kalisForwardTable")}
    if(any(targets < base.fwd.table$l)){stop("no targets may be less than base.fwd.table$l")}
    use.pi <- FALSE
  }


  if(force.unif){

    if(!use.pi){
      if(targets[1] != base.fwd.table$l){stop("When using force.unif, for now, base.fwd.table$l must be AT the first target")}
    }

    targets.idx <- targets
    targets <- 1:length(targets)
    base.fwd.table.l <- 1

  } else {
    targets.idx <- NULL
    base.fwd.table.l <- base.fwd.table$l
  }


  # Figure out whether using uniform or general checkpointing
  ####################
  # by default
  uniform.ckpts <- FALSE
  first.target.given <- FALSE

  intervals <- unique(diff(targets))

  if(length(intervals) == 1){
    if(use.pi){
      if(targets[1] == intervals){ uniform.ckpts <- TRUE }
    }else{
      if(targets[1] == base.fwd.table.l){ uniform.ckpts <- TRUE; first.target.given <- TRUE}
      if( (targets[1] - intervals) == base.fwd.table.l){ uniform.ckpts <- TRUE; first.target.given <- FALSE}
    }
  }
  rm(intervals)


  # Perform Table Benchmarking
  ####################
  # bench <- BenchmarkTables()

  propagation.cost <- 1:length(targets)


  # Solve Optimal Checkpointing
  ###############################

  if(uniform.ckpts){

    if(is.null(lookup.tables)){
      message("Calculating Optimal Checkpoint Schedule")
      lookup.tables <- CalcCheckpointTables(propagation.cost,num.available.ckpts)
    }

    cost.table <- lookup.tables$cost
    index.table <- lookup.tables$index




    if(!first.target.given){

      uniform_SolveSchedule <- uniform_MakeSolveSchedule(targets,cost.table,index.table)
      assign("uniform_SolveSchedule",uniform_SolveSchedule,envir = environment(uniform_SolveSchedule))

      uniform_SolveSchedule(1,length(targets),num.available.ckpts)

      cost <- uniform_LookupCost(length(targets),num.available.ckpts,cost.table)

    }else{

      uniform_SolveSchedule <- uniform_MakeSolveSchedule(targets[-1],cost.table,index.table)
      assign("uniform_SolveSchedule",uniform_SolveSchedule,envir = environment(uniform_SolveSchedule))

      uniform_SolveSchedule(1,length(targets)-1,num.available.ckpts)
      cost <- uniform_LookupCost(length(targets)-1,num.available.ckpts,cost.table)

    }

    sch <- uniform_trim.sch(uniform_SolveSchedule)

    # I don't believe we need to modify this schedule in order to still request the baseline locus as our last target

  }else{

    stop("Solving the non-uniform checkpointing problem is not yet implemented.")

  }


  # Construct Table Cache
  ########################

  max.tables <- max(sch$k)

  if(is.null(cache) || length(cache) < max.tables){
    cache <- CreateForwardTableCache(pars = pars,size = Inf, from_recipient = from_recipient, to_recipient = to_recipient, max.tables = max.tables)

  }else{

    for(i in 1:length(cache)){

      if(i > max.tables){
        cache[[i]] <- NULL
      }else{
        # check if parameters match, if not, overwrite with warning
        if(cache[[i]]$pars.sha256!=pars$sha256){
          warning("The provided cache was initialized with parameters that are different from those currently in pars.  Overwritting the pars in the provided cache...")
          cache[[i]]$pars.sha256 <- pars$sha256
        }
        ResetTable(cache[[i]])
      }
    }

  }

  rm(pars); gc()


  # Construct Iterator
  ######################

  UpdateCache <- MakeUpdateCache(sch, use.pi, targets.idx = targets.idx)

  current.sch <- data.frame("i" = Inf)

  current.target.index <- length(targets)

  iter.internal <- function(fwd, pars, t, nthreads = 1){

    if(force.unif){t <- match(t,targets.idx)}

    if(current.target.index == 0){warning("This iterator has been exhausted."); return()}

    if(t == targets[current.target.index]){
      current.target.index <<- current.target.index - 1
    }else{
      stop(paste("The next target locus for this iterator is",targets[current.target.index],"not",t))
    }

    if(identical(lobstr::obj_addr(fwd),lobstr::obj_addr(base.fwd.table))){
      stop("base.fwd.table cannot point to the same table as fwd: they must be created independently.")
    }

    #print(c(current.sch$i, t))
    if(current.sch$i > t){ current.sch <<- UpdateCache(cache, pars, nthreads, base.fwd.table) }

    if(current.sch$k != 0){
      CopyTable(to = fwd, from = cache[[current.sch$k]])
    }else{
      if(use.pi){
        ResetTable(fwd)
      }else{
        CopyTable(to = fwd, from = base.fwd.table)
      }
    }

    # Clean Up cache unless we're instructed to save it
    if(current.target.index == 0){
      if(!save.cache){
        cache <<- NULL
        gc()
      }
    }

    if(force.unif){
      Forward(fwd = fwd, pars = pars, t = targets.idx[t], nthreads = nthreads)
    } else {
      Forward(fwd = fwd, pars = pars, t = t, nthreads = nthreads)
    }
  }

  class(iter.internal) <- c("kalisIterator",class(iter.internal))

  iter.internal
}


targets <- function(x) { # put this declaration above and below because it seems that order determines whether targets is recognized
  UseMethod("targets")
}

#' @export
targets.kalisIterator <- function(x){
  if(!"kalisIterator"%in%class(x)){stop("argument must be a kalisIterator")}
  rev(get("targets", envir = environment(x)))
}

targets <- function(x) {
  UseMethod("targets")
}

#' @export
print.kalisIterator <- function(x, ...){
  if(!"kalisIterator"%in%class(x)){stop("argument must be a kalisIterator")}

  if(get("current.target.index", envir=environment(x)) == 0){
    "This is an exhausted kalisIterator."
  }else{
    target.range <- range(get("targets", envir = environment(x)))
    message(paste("A kalisIterator for",length(get("targets", envir = environment(x))),"targets ranging from",target.range[1],"to",target.range[2]),appendLF = TRUE)
    message(paste("Contains",get("max.tables", envir = environment(x)),"checkpoints using ~",utils::object.size(get("cache", envir = environment(x)))/1e9,"Gb of RAM"),appendLF = TRUE)
    message(paste("Next target locus:",get("targets", envir = environment(x))[get("current.target.index", envir = environment(x))]),appendLF = TRUE)
    message("",appendLF = TRUE)
  }
}

#' @importFrom graphics axis
#' @export
plot.kalisIterator <- function(x, ...){
  if(!"kalisIterator"%in%class(x)){stop("argument must be a kalisIterator")}
  sch <- get("sch",envir = environment(x))
  loci <- get("targets",envir = environment(x))
  plot(sch$i[-c(1,nrow(sch))],sch$k[-c(1,nrow(sch))],type="h",lwd=1,bty="n",ylab="K",xlab="locus",las=1,ylim=c(0,max(sch$k)),xlim=range(loci),xaxt="n",yaxt="n", ...)
  p.loci <- pretty(loci)
  axis(1,at = p.loci ,pos=0)
  axis(2,at = pretty(0:max(sch$k)),pos=p.loci[1],las=2)
}


#' Calculate Checkpoint Tables
#'
#' Calculate look up tables for solving optimal checkpointing problems with dynamic programming.
#'
#' @references
#' Christ, R.R., Wang, X., Aslett, L.J.M., Steinsaltz, D. and Hall, I. (2024) "Clade Distillation for Genome-wide Association Studies", bioRxiv 2024.09.30.615852. Available at: \doi{10.1101/2024.09.30.615852}.
#'
#' @param propagation.cost
#'        a non-negative vector such that `propagation.cost[i]` gives the relative amount of time or cost required to propagate `i` steps.
#' @param max.num.checkpoints
#'        the maximum number of checkpoints that should be considered when building the checkpoint table.
#' @param use.R
#'        a logical, when `TRUE` use base R rather than C implementation of table building.
#'        Defaults to `FALSE`.
#'
#' @return
#' A list containing:
#'
#' \describe{
#'   \item{`cost`}{the matrix \eqn{F} in Christ et al. (2024)}
#'   \item{`index`}{the matrix \eqn{H} in Christ et al. (2024)}
#' }
#'
#' @export
CalcCheckpointTables <- function(propagation.cost,max.num.checkpoints, use.R = FALSE){
  start <- proc.time()

  propagation.cost <- as.numeric(propagation.cost)

  max.n <- length(propagation.cost)

  # the first row corresponds to solving a 0 locus problem
  cost.table <- matrix(0,nrow=max.n + 1,ncol= max.num.checkpoints + 1)
  index.table <- matrix(0L,nrow=max.n,ncol= max.num.checkpoints)

  cost.table[,1] <- c(0,cumsum(as.numeric(propagation.cost)))

  if(use.R){

    for(k in 1:max.num.checkpoints){
      for(n in 1:max.n){
        # now solving a n long problem with k checkpoints
        v <- cost.table[1:n,k + 1] + propagation.cost[1:n] + cost.table[n:1,k]
        x <- which.min(v)
        index.table[n, k] <- x
        cost.table[n + 1,k + 1] <- v[x]
      }
      print(paste(k,"done at",c(proc.time() - start)[3]/3600,"hours from start."))
    }
  }else{
    invisible(.Call(CCall_OptCkpt, cost.table, index.table, propagation.cost))
  }

  cost.table <- cost.table[-1,]
  return(list("cost" = cost.table,"index" = index.table))
}


MakeUpdateCache <- function(sch, use.pi, targets.idx = NULL){

  force(targets.idx)

  exhausted <- FALSE
  current.ins <- leading.ins <- 1
  ancestor <- 1
  cost <- 0

  function(cache, pars, nthreads, base.fwd.table){

    if(exhausted){
      warning("This iterator has been exhausted.")
      return(data.frame("k" = 0L,"i" = 0L))
    }

    repeat{

      candidates <- which(sch$i[1:(current.ins-1)] <= sch$i[current.ins])
      ancestor <<- candidates[which.max(sch$i[candidates])]

      if(sch$i[leading.ins + 1] < sch$i[ancestor]){ # if the next checkpoint destination is on the left side of the current ancestor
        current.ins <<- ancestor
        return(sch[current.ins,])
      } else {
        leading.ins <<- leading.ins + 1
        current.ins <<- leading.ins
      }

      candidates <- which(sch$i[1:(current.ins-1)] <= sch$i[current.ins])
      ancestor <<- candidates[which.max(sch$i[candidates])]


      if(sch$i[current.ins] != 0){ # we're not at the end yet

        # Update Cache
        kk <- sch$k[current.ins]
        akk <- sch$k[ancestor]

        if(akk != 0){
          CopyTable(cache[[ kk ]],cache[[ akk ]])
        }else{
          if(use.pi){
            ResetTable(cache[[kk]]) # Pi could also be the baseline table here for the entire interval
          }else{
            CopyTable(to = cache[[kk]],base.fwd.table)
          }
        }

        # advance cache table from ancestor to current checkpoint destination
        if(!is.null(targets.idx)){
          Forward(cache[[kk]],pars,targets.idx[sch$i[current.ins]],nthreads)
        } else {
          Forward(cache[[kk]],pars,sch$i[current.ins],nthreads)
        }
        if(sch$i[current.ins + 1] > sch$i[current.ins]){
          next
        }else{
          return(sch[current.ins,])
        }

      }else{ # we are at the end
        exhausted <<- TRUE
        rm(cache)
        rm(sch, envir = parent.env(environment())) # remove large objects from memory
        gc()
        return(data.frame("k" = 0L,"i" = 0L))
      }
    }
  }
}




uniform_MakeSolveSchedule <- function(loci,cost.table,index.table){

  uniform_SolveSchedule <- function(){NULL}

  sch.k <- 0L
  sch.i <- 0L
  nrow.sch <- 1

  function(i,j,num.available.ckpts){ # i is the index of the first locus and j is the index of the last locus in the problem to solve (from indicies[i] to indicies[j])

    l.d <- j-i+1

    k <- as.integer(min(l.d-1,num.available.ckpts))
    if(k==0){return(cost.table[l.d,1])}

    # at this point, we know that num.available.ckpts is at least 1
    # and l.d is at least 2

    # If neither of the above cases, create a new instruction
    ins <- which.max(sch.k < 0) # this is the first emtpy slot for an instruction
    if(ins == nrow.sch){ # then we're about to assign to the last schedule entry and need to add on space for instructions before we can call obj.func
      sch.k <<- c(sch.k, rep(-1L,50))
      sch.i <<- c(sch.i, rep(-1L,50))
      nrow.sch <<- length(sch.k)
      ins <- which.max(sch.k < 0)
    }

    sch.k[ins] <<- k
    ckpt.location <- index.table[l.d,k]
    sch.i[ins] <<- loci[i-1+ckpt.location]

    # solve right problem if the interval to the right contains at least one target locus
    if(l.d > ckpt.location){ uniform_SolveSchedule(i+ckpt.location, j, num.available.ckpts - 1) }
    #if(l.d > ckpt.location){ get("uniform_SolveSchedule", envir = parent.frame())(i+ckpt.location, j, num.available.ckpts - 1) }


    # solve left problem if the interval to the left contains at least one target locus
    if(ckpt.location > 1){ uniform_SolveSchedule(i, i-2+ckpt.location,num.available.ckpts) }
    #if(ckpt.location > 1){ get("uniform_SolveSchedule", envir = parent.frame())(i, i-2+ckpt.location,num.available.ckpts) }


    return()
  }
}

uniform_LookupCost<- function(L,num.available.ckpts,cost.table){cost.table[L,as.integer(min(L-1,num.available.ckpts)) + 1]}

uniform_trim.sch <- function(f){
  sch.k <- get("sch.k",envir = environment(f))
  sch.i <- get("sch.i",envir = environment(f))
  # prune
  if(any(sch.k == -1)){
    upper.limit <- which.max(sch.k == -1) - 1
    sch.k <- sch.k[1:upper.limit]
    sch.i <- sch.i[1:upper.limit]
  }

  # create dataframe schedule
  sch <- data.frame("k" = c(sch.k,0L), "i" = c(sch.i,0L))
}
