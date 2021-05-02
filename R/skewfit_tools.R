#' Simulate skewed normal
#'
#' @export
#'
sf_sn_rnd  <- function(n, eta = 0, sigma = 1) {
    rst <- abs(rnorm(n));
    rst <- eta * rst + rnorm(n, sd = sigma);
}

#' Get Skewness of Half-Normal
#'
#' @export
#'
sf_sn_skew <- function(eta, sigma = 1, mu = 0) {
    w     <- sqrt(sigma^2 + eta^2)
    delta <- eta / w
    pi_2  <- 2 / pi

    mu    <- mu + w * delta * sqrt(pi_2)
    sig2  <- w^2 * (1 - delta^2 *  pi_2)
    skew <- (4 - pi) / 2 * (delta * sqrt(pi_2))^3
    skew <- skew / ((1 - pi_2 * delta^2)^1.5)

    c(mu   = mu,
      sig2 = sig2,
      skew = skew)
}

#' Get log likelihood from skewed normal
#'
#'
#' @export
#'
sf_sn_lpdf <- function(u, eta = 0, mu = 0, sigma = 1) {
    e2    <- eta^2;
    s2    <- sigma^2;
    e2ps2 <- sqrt(e2 + s2);

    x1  <- (u - mu) / e2ps2;
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

#'
#'
#' @export
get_lm_coeff_2 <- function(y, x = NULL) {
    if (is.null(x)) {
        lrs <- lm(y ~ 1)
    } else {
        lrs <- lm(y ~ x)
    }

    coefficients(lrs)
}

#'
#'
#' @export
get_lm_coeff_3 <- function(y, x = NULL, z = NULL) {

    if (is.null(x) & is.null(z)) {
        lrs <- lm(y ~ 1)
    } else {
        xz  <- cbind(x, z)
        lrs <- lm(y ~ xz)
    }

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
get_smooth_bs_x <- function(x, h = NULL, f_kernel = get_sample_kernel) {
    if (is.null(h)) {
        h <- bw.nrd(x)
    }

    smp_kernel <- f_kernel(length(x))
    smp_x      <- x + h * smp_kernel

    smp_x
}


## Single Bootstrap
##
##
get_single_bs <- function(rst_fit,
                          pred_x   = NULL,
                          residual = c("empirical", "parametric"),
                          smooth   = TRUE,
                          replace  = FALSE,
                          h        = -1,
                          ...,
                          seed     = NULL) {

    ## set random seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## original data
    is_skew  <- inherits(rst_fit, "SKEW")
    is_iso   <- inherits(rst_fit, "ISO")
    x        <- rst_fit$x
    z        <- rst_fit$z

    if (h < 1)
        ## bandwidth
        h <- bw.nrd(x)

    fit      <- predict(rst_fit, ...)
    beta_z   <- fit$pred_bz
    ax       <- fit$pred_ax
    ax_smh   <- fit$pred_ax_smooth

    if (is.null(pred_x)) {
        pred_x <- x
    }

    ## bootstrap samples
    inx <- sample(seq_len(length(x)), replace = replace)
    inx <- sort(inx)

    ## bootstrap x and alpha(x)
    smp_x <- x[inx]
    if (smooth) {
        smp_ax <- ax_smh[inx]
    } else {
        smp_ax <- ax[inx]
    }

    ## bootstrap z and beta * z
    smp_z  <- z[inx, ]
    smp_bz <- beta_z[inx]

    ## bootstrap residual
    smp_re <- sf_sample_residual(rst_fit,
                                 n      = length(x),
                                 method = residual)
    ## bootstrap outcome
    smp_y <- smp_ax + smp_bz + smp_re

    ## fit model to bootstrap sample
    if (is_iso & is_skew) {
        cur_rst <- sf_iso_skew(smp_y,  smp_x, as.matrix(smp_z), ...)
    } else if (is_iso & !is_skew) {
        cur_rst <- sf_iso_norm(smp_y,  smp_x, as.matrix(smp_z), ...)
    } else if (!is_iso & is_skew) {
        cur_rst <- sf_para_skew(smp_y, smp_x, as.matrix(smp_z), ...)
    } else {
        cur_rst <- sf_para_norm(smp_y, smp_x, as.matrix(smp_z), ...)
    }

    ## append result
    cur_pa  <- cur_rst$mle_pa
    cur_fit <- predict(cur_rst, new_x = pred_x)
    rst     <- c(cur_pa, cur_fit$pred_ax, cur_fit$pred_ax_smooth)

    ## reset random seed
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}
