#!/usr/local/bin/Rscript
library(stats)


# read in info (as from command line)
## args: eta_val, length
args = (commandArgs(trailingOnly=TRUE))
eta_val = as.numeric(args[1])
# NOTE: for simplicity treated lengths  here as a vector, but really just passing in a single value from the .py script
lengths = c(as.numeric(args[2]))





deltas = log10(lengths)
ns_minus_1 = pmax(round(lengths/eta_val/gamma(1/deltas + 1)) - 1, 0)
xis = lapply(ns_minus_1, function(n) {diff(sort(c(runif(n), 0, 1)))})
xis_transformed = mapply(function(x, d) {x^(1/d)}, xis, deltas, SIMPLIFY = F)
delta_is = mapply(function(len, x_t) {round(len*x_t/sum(x_t))}, lengths, xis_transformed, SIMPLIFY = F)
print(as.character(delta_is))




































