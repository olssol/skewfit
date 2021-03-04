#' Bootstrap
#'
#' @export
sf_bootstrap <- function(object, ...)
    UseMethod("sf_bootstrap")

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
                                 n_cores = 1,
                                 ...) {

    if (!is.null(seed)) {
        old_seed <- set.seed(seed)
    }

    if (is.null(pred_x)) {
        pred_x <- object$x
    }

    est <- c(object$mle_pa,
             predict(object, new_x = pred_x)$pred_ax)

    rst_bs <- parallel::mclapply(1:nbs,
                                 function(x) {
                                     if (1 == verbose)
                                         cat("Bootstrap ", x, "\n")

                                     get_single_bs(object, pred_x = pred_x, ...)
                                 },
                                 mc.cores = n_cores)

    rst_bs <- matrix(unlist(rst_bs), ncol = nbs, byrow = FALSE)

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
predict.ISO <- function(object, new_x = NULL, new_z = NULL, ...) {
    if (is.null(new_x)) {
        new_x <- object$x
    }

    if (is.null(new_z)) {
        new_z <- object$z
    }

    ## predict alpha(x)
    pred_ax <- pred_iso(new_x,
                        cbind(object$x, object$mle_ai))[, 2]

    ## predict beta z
    beta    <- object$mle_pa[seq_len(ncol(new_z))]
    pred_bz <- apply(cbind(new_z), 1,
                     function(w) sum(w * beta))

    list(new_x  = new_x,
         new_z  = new_z,
         pred_ax = pred_ax,
         pred_bz = pred_bz)
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

    list(new_x  = new_x,
         new_z  = new_z,
         pred_ax = pred_ax,
         pred_bz = pred_bz)
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
