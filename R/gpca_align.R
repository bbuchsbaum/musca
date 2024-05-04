
gpca_align.hyperdesign <- function(data, y, 
                    preproc=center(), 
                    ncomp=2,
                    simfun,
                    csimfun=NULL,
                    u=.5,
                    lambda=.1) {
  
  y <- rlang::enquo(y)
  
  label_list <- purrr::map(data, function(x) x$design %>% select(!!y) %>% pull(!!y))
  labels <- factor(unlist(label_list))
  label_set <- levels(labels)
  
  #subjects <- purrr::map(data, function(x) x$design %>% select(subject) %>% pull())
  M_within <- Matrix::bdiag(lapply(label_list, function(l) Matrix(simfun(l), sparse=TRUE)))
  M_between <- simfun(labels) - M_within
  M_between <- M_between/length(label_list)
  M <- u*M_within + (1-u)*M_between
  
  ninstances <- length(labels)
  nsets <- length(data)
  
  pdata <- multivarious::init_transform(data, preproc) 
  proclist <- attr(pdata, "preproc")
  
  names(proclist) <- names(pdata)
  
  block_indices <- block_indices(pdata)
  proc <- multivarious::concat_pre_processors(proclist, block_indices)
  names(block_indices) <- names(pdata)
  #browser()
  #M <- simfun(labels)
  M <- M + Matrix::Diagonal(x=rep(lambda, nrow(M)))
  evm <- PRIMME::eigs_sym(M, NEig=1,  which="LA",method='PRIMME_DEFAULT_MIN_MATVECS')
  M <- M/evm$values[1]
  
  X_block <- Matrix::bdiag(lapply(pdata, function(x) x$x))
  ret <- genpca::genpca(X_block, M=.1*M, ncomp=ncomp, preproc=multivarious::pass())

  
  multivarious::multiblock_biprojector(
    v=ret$v,
    s=ret$s,
    sdev=ret$sdev,
    preproc=proc,
    block_indices=block_indices,
    labels=labels,
    #u=u,
    classes="gpca_align"
  )
  
  
}