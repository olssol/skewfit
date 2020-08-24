#'
#' least square estimator with constraints on the slope of x
#'
#' @export
#'
get_lm_posx <- function(y, x, z) {
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

    rst <- optim(init_pars,
                 method = "L-BFGS-B",
                 fn     = f_lsq,
                 lower  = lower, upper = upper)$par
}

#'
#'
#'
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
