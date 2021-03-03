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
    rst_ci
}
