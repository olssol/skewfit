get.lm.coeff <- function(y, x) {
    lrs <- lm(y ~ -1 + x);
    coefficients(lrs);
}
