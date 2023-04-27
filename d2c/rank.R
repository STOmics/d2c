### Get the parameters
parser = argparse::ArgumentParser(description="Script to beads calling")
parser$add_argument('-i','--input', help='input file')
parser$add_argument('-t','--type', help='knee or inflection')
parser$add_argument('-o','--out', help='the out put pdf')
args = parser$parse_args()

#library(DropletUtils)
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(S4Vectors)))


barcodeRanks <- function(m, lower=100, fit.bounds=NULL, exclude.from=50, df=20, ...) {
    totals <- unname(colSums(m))
    o <- order(totals, decreasing=TRUE)

    stuff <- rle(totals[o])
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 # Get mid-rank of each run.
    run.totals <- stuff$values

    keep <- run.totals > lower
    if (sum(keep)<3) {
        stop("insufficient unique points for computing knee/inflection points")
    }
    y <- log10(run.totals[keep])
    x <- log10(run.rank[keep])

    # Numerical differentiation to identify bounds for spline fitting.
    edge.out <- .find_curve_bounds(x=x, y=y, exclude.from=exclude.from)
    left.edge <- edge.out["left"]
    right.edge <- edge.out["right"]

    # As an aside: taking the right edge to get the total for the inflection point.
    # We use the numerical derivative as the spline is optimized for the knee.
    inflection <- 10^(y[right.edge])

    # We restrict curve fitting to this region, thereby simplifying the shape of the curve.
    # This allows us to get a decent fit with low df for stable differentiation.
    if (is.null(fit.bounds)) {
        new.keep <- left.edge:right.edge
    } else {
        new.keep <- y > log10(fit.bounds[1]) & y < log10(fit.bounds[2])
    }

    # Smoothing to avoid error multiplication upon differentiation.
    # Minimizing the signed curvature and returning the total for the knee point.
    fitted.vals <- rep(NA_real_, length(keep))

    if (length(new.keep) >= 4) {
        fit <- smooth.spline(x[new.keep], y[new.keep], df=df, ...)
        fitted.vals[keep][new.keep] <- 10^fitted(fit)

        d1 <- predict(fit, deriv=1)$y
        d2 <- predict(fit, deriv=2)$y
        curvature <- d2/(1 + d1^2)^1.5
        knee <- 10^(y[which.min(curvature)])
    } else {
        # Sane fallback upon overly aggressive filtering by 'exclude.from', 'lower'.
        knee <- 10^(y[new.keep[1]])
    }

    # Returning a whole stack of useful stats.
    out <- DataFrame(
        rank=.reorder(run.rank, stuff$lengths, o),
        total=.reorder(run.totals, stuff$lengths, o),
        fitted=.reorder(fitted.vals, stuff$lengths, o)
    )
    rownames(out) <- colnames(m)
    metadata(out) <- list(knee=knee, inflection=inflection)
    out
}

.reorder <- function(vals, lens, o) {
    out <- rep(vals, lens)
    out[o] <- out
    return(out)
}

.find_curve_bounds <- function(x, y, exclude.from)
# The upper/lower bounds are defined at the plateau and inflection, respectively.
# Some exclusion of the LHS points avoids problems with discreteness.
{
    d1n <- diff(y)/diff(x)

    skip <- min(length(d1n) - 1, sum(x <= log10(exclude.from)))
    d1n <- tail(d1n, length(d1n) - skip)

    right.edge <- which.min(d1n)
    left.edge <- which.max(d1n[seq_len(right.edge)])

    c(left=left.edge, right=right.edge) + skip
}

type = args$type
out = args$out
input <- as.data.frame(fread(args$input))
if (ncol(input) == 2){
    m <- as.matrix(subset(input,select=c("V2")))
    br.out <- barcodeRanks(t(m), lower = 500)
}else{
    m <- as.matrix(subset(input,select=c("jaccard_distance")))
    br.out <- barcodeRanks(t(m), lower = 0.005)
}
br.out.sort <- br.out[order(br.out$total, decreasing=T),]

if (ncol(input) == 2){
        if(type == "inflection"){
            cutoff <- metadata(br.out)$inflection
        }else if(type == "knee"){
            cutoff <- metadata(br.out)$knee
        }else{
            cutoff <- 0
            cat("Please fill in the correct parameters,inflection or knee")
        }
    d2cCutoff=data.frame(name="bead_cutoff",cutoff=cutoff)
    write.table(d2cCutoff,out,quote=F,row.names=F,col.names=F,sep="\t")
}else{
    if(type == "inflection"){
        cutoff <- metadata(br.out)$inflection
    }else if(type == "knee"){
        cutoff <- metadata(br.out)$knee 
    }else{
        cat("Please fill in the correct parameters,inflection or knee")
    }
    if (file.exists(out)){
        d2cCutoff <- as.data.frame(fread(out,header = F))
        d2cCutoff[2,1] <- "cor_cutoff"
        d2cCutoff[2,2] <- cutoff
        write.table(d2cCutoff,out,quote=F,row.names=F,col.names=F,sep="\t")

    }else{
        cat("Please do beads calling first!")
    }
}

