

close(filename)
filename <- file("corrtest", "rb")
test <- readBin(filename, character(),4)
cor.matrix <- readBin(filename, "numeric", 4)
hes.vec <- readBin(filename, "numeric", num.pars^2)
hes <- matrix(hes.vec, ncol=num.pars, nrow=num.pars)
hybrid_bounded_flag <- readBin(filename, "integer", 1)
scale <- readBin(filename, "numeric", num.pars)
result <- list(num.pars=num.pars, hes=hes,
               hybrid_bounded_flag=hybrid_bounded_flag, scale=scale)
