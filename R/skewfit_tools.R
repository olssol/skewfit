#' Get log likelihood from skewed normal
#'
#'
#' @export
#'
get_sn_lpdf <- function(u, eta = 0, mu = 0, sigma = 1) {
    e2    <- eta^2;
    s2    <- sigma^2;
    e2ps2 <- sqrt(e2 + s2);

    x1  <- (u - mu)/e2ps2;
    x2  <- (u - mu) * eta / e2ps2 / sigma;
    rst <- log(2) -  log(e2ps2);
    rst <- rst + dnorm(x1, log = TRUE) + log(pnorm(x2));

    sum(rst)
}

#'
#'
#' @export
get_lm_coeff <- function(y, x) {
    ## keep for pseudo x
    if (all(0 == x))
        return(0)

    lrs <- lm(y ~ -1 + x)
    coefficients(lrs)
}

get_lm_coeff_2 <- function(y, x) {
    lrs <- lm(y ~ x)
    coefficients(lrs)
}


#' Get bootstrap confidence interval
#'
#'
#' @export
#'
get_bs_ci <- function(bs, method = c("empirical"), quants = c(0.025, 0.975),
                      est = NULL) {

    method <- match.arg(method)
    rst_ci <- switch(method,
                     empirical = {
                         ci <-  apply(bs, 1, function(x) quantile(x, quants))
                         t(ci)
                     })

    rst_sd   <- apply(bs, 1, sd)
    rst_mean <- apply(bs, 1, mean)

    cbind(mean  = rst_mean,
          sd    = rst_sd,
          ci_lb = rst_ci[, 1],
          ci_ub = rst_ci[, 2])
}

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
get_smooth_bs_x <- function(x, h = NULL) {
    if (is.null(h)) {
        h <- bw.nrd(x)
    }

    smp_kernel <- get_sample_kernel(length(x))
    smp_x      <- x + h * smp_kernel

    smp_x
}


## Single Bootstrap
##
##
get_single_bs <- function(rst_fit, smooth_bs = TRUE, pred_x = NULL,
                          method_residual = c("empirical", "parametric"),
                          ...) {

    x        <- rst_fit$x
    z        <- rst_fit$z
    pa       <- rst_fit$mle_pa
    ai       <- rst_fit$mle_ai
    residual <- rst_fit$residual

    n        <- length(x)
    h        <- bw.nrd(x)
    beta_z   <- predict(rst_fit)$pred_bz

    is_skew  <- inherits(rst_fit, "SKEW")
    is_iso   <- inherits(rst_fit, "ISO")

    if (is.null(pred_x)) {
        pred_x <- x
    }

    if (smooth_bs) {
        smp_x  <- get_smooth_bs_x(x, h)
        ord_x  <- order(smp_x)
        smp_x  <- smp_x[ord_x]
        smp_fx <- predict(rst_fit, new_x = smp_x)$pred_ax
    } else {
        ## typical bs
        ord_x  <- seq_len(length(x))
        smp_x  <- x
        smp_fx <- ai
    }

    ## bs covariates z
    smp_z  <- z[ord_x, ]
    smp_bz <- beta_z[ord_x]

    ## bs residual
    smp_re <- sf_sample_residual(rst_fit, method = method_residual)

    ## bs outcome
    smp_y  <- smp_fx + smp_bz + smp_re

    ## fit
    if (is_iso & is_skew) {
        cur_rst <- sf_iso_skew(smp_y, smp_x, as.matrix(smp_z), ...)
    } else if (is_iso & !is_skew) {
        cur_rst <- sf_iso_norm(smp_y, smp_x, as.matrix(smp_z), ...)
    } else if (!is_iso & is_skew) {
        cur_rst <- sf_para_skew(smp_y, smp_x, as.matrix(smp_z), ...)
    } else {
        cur_rst <- sf_para_norm(smp_y, smp_x, as.matrix(smp_z), ...)
    }

    ## append rst
    cur_pa  <- cur_rst$mle_pa
    cur_ai  <- predict(cur_rst, new_x = pred_x)$pred_ax
    c(cur_pa, cur_ai)
}
