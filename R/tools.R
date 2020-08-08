#'
#' least square estimator with constraints on the slope of x
#'
#' @export
#'
get.lm.posx <- function(y, x, z) {
    ## residual
    f.res <- function(pars) {
        beta <- pars[-(1:2)];
        bz   <- apply(z, 1, function(zz) sum(beta * zz));
        rst  <- y - pars[1] - pars[2] * x - bz;
    }

    ## least square
    f.lsq <- function(pars) {
        res  <- f.res(pars);
        lsq  <- sum(res^2);
    }

    f.gradient <- function(pars) {
        res  <- - 2 * f.res(pars);
        g    <- numeric(length(pars));
        g[1] <- sum(res * 1);
        g[2] <- sum(res * x);

        for (k in 1:ncol(z)) {
            g[2+k] = sum(res * z[,k]);
        }
        g
    }

    lm.est    <- coefficients(lm(y ~ z));
    init.pars <- c(lm.est[1], 0, lm.est[-1]);

    lower     <- rep(-Inf, length(init.pars));
    lower[2]  <- -1e-8;
    upper     <- rep(Inf, length(init.pars));

    rst <- optim(init.pars,
                 method = "L-BFGS-B",
                 fn     = f.lsq,
                 lower  = lower, upper = upper)$par;
}

#'
#'
#'
get.lm.coeff <- function(y, x) {
    ## keep for pseudo x
    if (all(0 == x))
        return(0)

    lrs <- lm(y ~ -1 + x);
    coefficients(lrs);
}

get.lm.coeff.2 <- function(y, x) {
    lrs <- lm(y ~ x);
    coefficients(lrs);
}
