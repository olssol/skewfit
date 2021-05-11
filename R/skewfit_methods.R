#' Bootstrap
#'
#' @export
sf_bootstrap <- function(object, ...)
    UseMethod("sf_bootstrap")

#' Log-Likelihood
#'
#' @export
sf_lpdf <- function(object, ...)
    UseMethod("sf_lpdf")

#' Bootstrap confidence intervals
#'
#' @export
sf_bs_ci <- function(object, ...)
    UseMethod("sf_bs_ci")

#' Sample residual
#'
#' @export
sf_sample_residual <- function(object, ...)
    UseMethod("sf_sample_residual")


#' Bootstrap
#'
#' @method sf_bootstrap SKEWFIT
#'
#' @export
#'
sf_bootstrap.SKEWFIT <- function(object, nbs = 100, seed = NULL,
                                 verbose = 1, pred_x = NULL,
                                 n_cores = 1, ...) {

    ## set random seed
    if (!is.null(seed)) {
        old_seed <- set.seed(seed)
    }

    ## boostrap seed
    bs_seeds <- floor(abs(rnorm(nbs)) * 100000)

    ## predict x
    if (is.null(pred_x)) {
        pred_x <- object$x
    }

    ## estimate
    pred_fit <- predict(object, new_x = pred_x)
    est      <- c(object$mle_pa,
                  pred_fit$pred_ax,
                  pred_fit$pred_ax_smooth)

    ## bootstrap
    rst_bs <- parallel::mclapply(1:nbs,
                                 function(x) {
                                     if (1 == verbose & 0 == x %% 50)
                                         cat("Bootstrap ", x, "\n")

                                     get_single_bs(object,
                                                   pred_x = pred_x,
                                                   seed   = bs_seeds[x],
                                                   ...)
                                 },
                                 mc.cores = n_cores)

    rst_bs <- matrix(unlist(rst_bs),
                     ncol = nbs,
                     byrow = FALSE)

    ## reset seed
    if (!is.null(seed)) {
        set.seed(old_seed)
    }

    ## return
    rst <- list(est    = est,
                bs     = rst_bs,
                pred_x = pred_x)

    class(rst) <- append(class(object), "BOOTSTRAP")

    rst
}

#' Prediction from fitted isotonic model
#'
#' @method predict ISO
#'
#' @export
#'
predict.ISO <- function(object,
                        new_x = NULL, new_z = NULL, h = -1,
                        ...) {
    if (is.null(new_x)) {
        new_x <- object$x
    }

    if (is.null(new_z)) {
        new_z <- object$z
    }

    ## predict alpha(x)
    if (h <= 0) {
        h <- bw.nrd(object$x)
    }

    x_ax           <- cbind(object$x, object$mle_ai)
    pred_ax        <- pred_iso(new_x, x_ax)[, 2]
    pred_ax_smooth <- pred_iso(new_x, x_ax, h = h, ...)[, 2]

    ## predict beta z
    beta    <- object$mle_pa[seq_len(ncol(new_z))]
    pred_bz <- apply(cbind(new_z), 1,
                     function(w) sum(w * beta))

    list(new_x          = new_x,
         new_z          = new_z,
         pred_ax        = pred_ax,
         pred_ax_smooth = pred_ax_smooth,
         pred_bz        = pred_bz)
}

#' Prediction from fitted parametric model
#'
#' @method predict PARA
#'
#' @export
#'
predict.PARA <- function(object, new_x = NULL, new_z = NULL, ...) {
    if (is.null(new_x)) {
        new_x <- object$x
    }

    if (is.null(new_z)) {
        new_z <- object$z
    }

    ## predict alpha(x)
    pred_ax <- object$mle_pa[1] + object$mle_pa[2] * new_x;

    ## predict beta z
    beta    <- object$mle_pa[2 + seq_len(ncol(new_z))]
    pred_bz <- apply(cbind(new_z), 1,
                     function(w) sum(w * beta))

    list(new_x          = new_x,
         new_z          = new_z,
         pred_ax        = pred_ax,
         pred_ax_smooth = pred_ax,
         pred_bz        = pred_bz)
}

#' Sample residual from SKEW model
#'
#' @method sf_sample_residual SKEW
#'
#' @export
#'
sf_sample_residual.SKEW <- function(object, n,
                                    method = c("empirical", "parametric"),
                                    ...) {

    method   <- match.arg(method)
    residual <- object$residual
    pa       <- object$mle_pa

    smp_re <- switch(method,
                     empirical  = sample(residual, n, replace = TRUE),
                     parametric = {
                         wi   <- rnorm(n)
                         epsi <- rnorm(n, 0, sqrt(pa["sig2"]))
                         res  <- pa["eta"] * wi + epsi})

    smp_re
}

#' Sample residual from normal model
#'
#' @method sf_sample_residual NORM
#'
#' @export
#'
sf_sample_residual.NORM <- function(object, n,
                                    method = c("empirical", "parametric"),
                                    ...) {

    method   <- match.arg(method)
    residual <- object$residual
    pa       <- object$mle_pa

    smp_re <- switch(method,
                     empirical  = sample(residual, n, replace = TRUE),
                     parametric = rnorm(n, 0, sqrt(pa["sig2"])))
    smp_re
}

#' Log-Likelihood from SKEW model
#'
#' @method sf_lpdf SKEW
#'
#' @export
#'
sf_lpdf.SKEW <- function(object, ...) {
    residual <- object$residual
    pa       <- object$mle_pa

    if (is.null(pa))
        return(-Inf)

    sf_sn_lpdf(residual,
                eta   = pa["eta"],
                sigma = sqrt(pa["sig2"]))
}

#' Log-Likelihood from NORM model
#'
#' @method sf_lpdf NORM
#'
#' @export
#'
sf_lpdf.NORM <- function(object, ...) {
    residual <- object$residual
    pa       <- object$mle_pa

    sum(dnorm(residual, sd = sqrt(pa["sig2"]),
              log = TRUE))
}

#' Get confidence interval from bootstrap results
#'
#' @method sf_bs_ci BOOTSTRAP
#'
#' @export
#'
sf_bs_ci.BOOTSTRAP <- function(object, quants = c(0.025, 0.975), ...) {
    bs       <- object$bs
    pred_x   <- object$pred_x
    est      <- object$est

    n_ax     <- length(pred_x)
    n_pa     <- nrow(bs) - n_ax * 2
    inx      <- 1:(n_pa + n_ax)
    inx_smth <- c(1:n_pa, (n_pa + n_ax) + 1:n_ax)

    rst_est  <- est[inx]
    rst_sd   <- apply(bs, 1, sd)[inx]

    ## empirical ci
    rst_ci   <- apply(bs[inx, ], 1,
                      function(x) quantile(x, quants, na.rm = TRUE))

    ## correct for smoothness
    rst_ci_2 <- apply(cbind(est[inx],
                            est[inx_smth],
                            bs[inx_smth, ]), 1,
                      function(x) {
                          ci  <- x[-(1:2)] - x[2]
                          ci  <- quantile(ci, quants, na.rm = TRUE)
                          ci  <- sort(x[1] - ci)
                      })

    ## result
    data.frame(est     = rst_est,
               bs_sd   = rst_sd,
               emp_lb  = rst_ci[1, ],
               emp_ub  = rst_ci[2, ],
               adj_lb  = rst_ci_2[1, ],
               adj_ub  = rst_ci_2[2, ])
}

#' Plot predicted alpha(x)
#'
#'
#' @method plot BOOTSTRAP
#'
#' @export
#'
plot.BOOTSTRAP <- function(object, ci = TRUE, true_ax = NULL, ...) {
    x   <- object$pred_x
    y   <- object$est
    inx <- (1 + length(y) - length(x)):length(y)
    y   <- y[inx]

    data  <- data.frame(x = x, y = y)

    ## get ci
    if (ci) {
        bs_ci <- sf_bs_ci(object)
        bs_ci <- bs_ci[inx, ]

        data$ci_lb <- bs_ci$adj_lb
        data$ci_ub <- bs_ci$adj_ub
    }

    if (inherits(object, "ISO")) {
        f_g <- geom_step
    } else {
        f_g <- geom_line
    }

    rst <- ggplot(data = data, aes(x = x, y = y)) +
        f_g() +
        theme_bw()

    if (ci) {
        rst <- rst +
            f_g(mapping = aes(x = x, y = ci_lb), lty = 2, color = "gray") +
            f_g(mapping = aes(x = x, y = ci_ub), lty = 2, color = "gray")
    }

    if (!is.null(true_ax)) {
        rst <- rst +
            geom_line(data = data.frame(x = x, y = true_ax),
                      mapping = aes(x = x, y = y),
                      color = "red", lty = 2)
    }
    rst
}
