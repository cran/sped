
kinship <- function(individuals, pedigree) {
    stopifnot(is.character(individuals) || is.numeric(individuals))
    if (is.numeric(individuals))
        storage.mode(individuals) <- "integer"
    if (is.integer(individuals) && any(individuals <= 0))
        stop("individuals, if integer-valued, must be positive-valued")
    stopifnot(is.matrix(pedigree))
    stopifnot(ncol(pedigree) == 3)
    stopifnot(is.character(pedigree) || is.numeric(pedigree))
    if (is.numeric(pedigree))
        storage.mode(pedigree) <- "integer"
    if (is.integer(pedigree) && any(pedigree <= 0))
        stop("pedigree, if integer-valued, must be positive-valued")
    stopifnot(typeof(individuals) == typeof(pedigree))
    stopifnot(individuals %in% pedigree)

    from <- c(pedigree[,1], pedigree[,1])
    to <- c(pedigree[,2], pedigree[,3])
    foo <- try(tsort(from, to))
    if (inherits(foo, "try-error"))
        stop("some individual is its own ancestor")

    porder <- match(foo, pedigree[,1])
    ped.too <- pedigree[porder, ]
    pa <- match(ped.too[,2], foo)
    ma <- match(ped.too[,3], foo)
    pa[is.na(pa)] <- 0L
    ma[is.na(ma)] <- 0L

    iargs <- match(individuals, foo)

    kout <- .Call(C_kinship, pa = pa, ma = ma)
    dimnames(kout) <- list(foo, foo)
    kout[foo %in% individuals, foo %in% individuals]
}

