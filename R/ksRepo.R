#'ksRepo
#'
#'\code{ksRepo} performs KS-based repositioning of user supplied data. KS values
#'are assigned p-values by bootstrapping the user supplied datasets.
#'@param ccgenes A character vector of prioritized genes derived from a
#'  case/control dataset. Examples include differential gene expression,
#'  differential methylation, and mutation frequency.
#'@param database A named list object representing a gene::compound interaction
#'  database. List items should be named with compound names and their contents
#'  should be interacting genes.
#'@param boot.method A string specifying whether the 'compound' gene lists or
#'  the case/control, 'instance' dataset should be permuted for bootstrapping.
#'  Default is to permute the compound gene lists.
#'@param nboot An integer specifying the number of bootstrap samples for p-value
#'  calculation.
#'@return A data frame with C rows and 5 columns, where C is the number of
#'  compounds in the database. compound gives name of the compound. n.genes
#'  gives the number of genes associated with each compound. ks gives the raw ks
#'  statistic, boot.p gives the bootstrapped p-value, and boot.fdr gives the
#'  fdr-corrected bootstrapped p-value.
#'@examples 1+1
#'@export
ksRepo <- function(ccgenes, database, boot.method='compound', nboot = 10000) {
    # Ensure compatibility
    if (class(ccgenes)!='character') stop('ccgenes must be a vector of class character')
    if (class(database) != 'list') stop('database must be a list object')

    # Pre-process
    genes <- unique(ccgenes[!is.na(ccgenes)])
    db <- ks_compat(genes, database)

    # Initialize Output
    out <- data.frame(compound=names(db), stringsAsFactors = F)
    out$n.genes <- sapply(db,length)
    out$ks <- sapply(db, function(x) ks_single(genes,x))

    # P-value calculation
    if (boot.method == 'compound') {
        out$boot.p <- ks_boot_comp(genes, db, nboot, out$n.genes, out$ks)
    }
    else if (boot.method == 'instance'){
        out$boot.p <- ks_boot_inst(genes, db, nboot, out$ks)
    }
    else stop('boot.method must be either "compound" or "instance"')

    # FDR
    out$boot.fdr <- stats::p.adjust(out$boot.p,'fdr')

    #Return
    return(out)
}

ks_compat <- function(genes, database) {
    db <- lapply(database, function(x) intersect(x, genes))
    out <- paste(sum(sapply(db, length) == 0), 'compound(s) excluded due to lack of overlap.')
    print(out)
    db <- db[!sapply(db, length) == 0]
    return(db)
}

ks_single <- function(genes, dbitem) {
    # Find overlap and sort
    V <- sort(match(dbitem,genes))
    #Find a and b values
    t <- length(V)
    j <- 1:t
    n <- length(genes)
    a <- max(j/t - V/n)
    b <- max(V/n - (j-1)/t)
    #Compute KS value
    if (a > b) ks <- a
    else ks <- -b
    #Return
    return(ks)
}

ks_boot_comp <- function(genes, db, nboot, n.genes, ks) {
    ts <- unique(sapply(db,length))
    # Initialize bootstrap
    boot_output <- matrix(data=NA,nrow=nboot,ncol=max(ts))
    # Bootstrap
    boot_output <- sapply(ts, # Across possible ts
                  function(t) replicate(nboot, # replicate nboot times
                                        ks_single(genes, sample(genes, t, replace=T)) # perform ks on samples
                                        )
                  )
    colnames(boot_output) <- ts
    # Initialize output
    p.boot <- rep(NA, length(n.genes))
    # Get p-vals
    for (i in 1:length(n.genes)) {
        this.ks <- ks[i]
        this.len <- n.genes[i]
        this.boot <- boot_output[,as.character(this.len)]
        p.boot[i] <- sum(this.ks <= this.boot)/nboot
    }
    return(p.boot) # Return
}

ks_boot_inst <- function(genes, db, nboot, ks) {
    # Resample case/control list
    p <- rep(NA, length(db))
    for (i in 1:nboot) {
        this.sample <- sample(genes, length(genes))
        this.ks <- sapply(db, function(x) ks_single(this.sample,x))
        p <- p + as.integer(ks<this.ks)
    }
    p <- p/nboot
    return(p)
}
