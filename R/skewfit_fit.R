#' Least square estimator with constraints on the slope of x
#'
#' @export
#'
sf_para_norm <- function(y, x, z, ...) {
    ## residual
    f_res <- function(pars) {
        beta <- pars[- (1 : 2)]
        bz   <- apply(z, 1, function(zz) sum(beta * zz))
        rst  <- y - pars[1] - pars[2] * x - bz
        rst
    }

    ## least square
    f_lsq <- function(pars) {
        res  <- f_res(pars)
        lsq  <- sum(res^2)
        lsq
    }

    f_gradient <- function(pars) {
        res  <- - 2 * f_res(pars)
        g    <- numeric(length(pars))
        g[1] <- sum(res * 1)
        g[2] <- sum(res * x)

        for (k in seq_len(ncol(z))) {
            g[2 + k] <- sum(res * z[, k])
        }
        g
    }

    lm_est    <- coefficients(lm(y ~ z))
    init_pars <- c(lm_est[1], 0, lm_est[-1])

    lower     <- rep(-Inf, length(init_pars))
    lower[2]  <- -1e-8
    upper     <- rep(Inf, length(init_pars))
    mle_pa    <- optim(init_pars,
                       method = "L-BFGS-B",
                       fn     = f_lsq,
                       lower  = lower, upper = upper)$par

    residual <- f_res(mle_pa)
    sig2     <- var(residual)

    rst <- list(mle_pa   = c(mle_pa, sig2 = sig2),
                mle_ai   = mle_pa[1] + mle_pa[2] * x,
                residual = residual)

    ## return
    rst$x <- x
    rst$z <- z
    rst$y <- y

    rst$method <- "sf_para_norm"
    class(rst) <- c("SKEWFIT", "PARA", "NORM")

    rst
}

#' Parametric model with Skew Normal Errors
#'
#' Parametric model with Skew Normal errors with constraints on the slope of x
#'
#' @export
#'
sf_para_skew <- function(y, x, z, usez = 1, init_pa = NULL,
                         max_steps = 100000, tol = 1e-6, ...) {

    if (0 == usez)
        init_lm <- lm(y ~ 1)
    else
        init_lm <- lm(y ~ z)

    if (is.null(init_pa)) {
        beta    <- coefficients(init_lm)
        beta    <- append(beta, 0, after = 1)
        init_pa <- c(beta, eta = 1, sig2 = 1)
    }

    rst <- fit_para_skew(init_pa,
                         y,
                         x,
                         z,
                         usez,
                         max_steps = max_steps,
                         tol = tol)
    ## return
    rst$x <- x
    rst$z <- z
    rst$y <- y

    rst$init_pa <- init_pa
    rst$method  <- "sf_para_skew"

    class(rst) <- c("SKEWFIT", "PARA", "SKEW")
    rst
}


#' Isotonic regression with normal error
#'
#' @export
#'
sf_iso_norm <- function(y, x, z, init_pa = NULL, init_ai = NULL,
                        max_steps = 100000, tol = 1e-6, unimodel = 0) {

    if (is.null(init_ai))
        init_ai <- isoreg(y)$yf

    if (is.null(init_pa)) {
        init_pa <- c(rep(0, ncol(z)), sig2 = 1, mode = max(x))
    }

    ## fit
    rst <- fit_iso_norm(init_pa, init_ai, y, x, z,
                        unimodal  = unimodel,
                        max_steps = max_steps,
                        tol       = tol)

    ## return
    rst$x <- x
    rst$z <- z
    rst$y <- y

    rst$init_ai <- init_ai
    rst$init_pa <- init_pa
    rst$method  <- "sf_iso_norm"

    class(rst) <- c("SKEWFIT", "ISO", "NORM")
    rst
}

#' Isotonic regression with skewed error
#'
#' @export
#'
sf_iso_skew <- function(y, x, z, init_pa = NULL, init_ai = NULL,
                        max_steps = 100000, tol = 1e-6,
                        unimodel = 0, usez = 1) {

    if (is.null(init_ai))
        init_ai <- isoreg(y)$yf

    if (is.null(init_pa)) {
        if (0 < usez)
            init_pa <- rep(0, ncol(z))

        init_pa <- c(init_pa,
                     eta = 4, sig2 = 0.1, mode = max(x))
    }

    rst <- fit_iso_skew(init_pa, init_ai, y, x, z,
                        unimodal  = unimodel,
                        usez      = usez,
                        max_steps = max_steps,
                        tol       = tol)

    ## return
    rst$x <- x
    rst$z <- z
    rst$y <- y
    rst$init_ai <- init_ai
    rst$init_pa <- init_pa
    rst$method  <- "sf_iso_skew"
    class(rst)  <- c("SKEWFIT", "ISO", "SKEW")

    rst
}
