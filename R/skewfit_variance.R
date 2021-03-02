
#' Sample from Epanchenikov
#'
#' @details
#' 1. Generate iid 𝑈1,𝑈2,𝑈3 ~ Uniform(-1,1).
#' 2. If |𝑈3|≥|𝑈2| and |𝑈3|≥|𝑈1|, deliver 𝑈2; otherwise deliver 𝑈3.
#'
#' @export
#'
get_sample_kernel <- function(n) {
    u1 <- runif(n, -1, 1)
    u2 <- runif(n, -1, 1)
    u3 <- runif(n, -1, 1)

    apply(cbind(u1, u2, u3), 1, function(x) {
        ax <- abs(x)
        if (ax[3] >= ax[2] &
            ax[3] >= ax[1]) {
            rst <- x[2]
        } else {
            rst <- x[3]
        }
        rst
    })
}

#' Smoothed Bootstrap
#'
#'
#' @export
#'
get_smooth_bs_x <- function(x, fit_iso, h = NULL) {
    if (is.null(h)) {
        h <- bw.nrd(x)
    }

    smp_kernel <- get_sample_kernel(length(x))
    smp_x      <- x + h * smp_kernel
    rst        <- pred_iso(smp_x, fit_iso)

    rst
}


#' Bootstrap ISO SKEW
#'
#'
#' @export
#'
bs_iso_skew <- function(rst_fit, nbs = 100, seed = NULL, ...) {
    if (!is.null(seed)) {
        old_seed <- set.seed(seed)
    }

    x        <- rst_fit$x
    z        <- rst_fit$z
    pa       <- rst_fit$mle_pa
    ai       <- rst_fit$mle_ai
    residual <- rst_fit$residual

    n        <- length(residual)
    h        <- bw.nrd(x)

    ## beta_z
    beta   <- pa[seq_len(length(pa) - 3)]
    beta_z <- apply(cbind(z), 1, function(x) sum(x * beta))

    ## bootstrap
    rst <- NULL
    for (i in seq_len(nbs)) {
        smp_bs <- get_smooth_bs_x(x, cbind(x, ai), h)
        ord_x  <- order(smp_bs[, 1])
        smp_x  <- smp_bs[ord_x, 1]
        smp_fx <- smp_bs[ord_x, 2]

        smp_z  <- z[ord_x, ]
        smp_bz <- beta_z[ord_x]
        smp_re <- sample(residual, n, replace = TRUE)
        smp_y  <- smp_fx + smp_bz + smp_re

        ## fit
        cur_rst <- sf_iso_skew(smp_y, smp_x, as.matrix(smp_z), ...)

        cur_ai  <- cur_rst$mle_ai
        cur_pa  <- cur_rst$pa
        cur_ai0 <- pred_iso(x, cbind(smp_x, cur_ai))[, 2]

        ## append rst
        rst <- cbind(rst, c(cur_pa, cur_ai0))
    }

    if (!is.null(seed)) {
        set.seed(old_seed)
    }

    rst
}
