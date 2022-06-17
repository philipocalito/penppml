### Identify separated observations in poissonhdfe regression ####
# Script implements Iterative Rectifier Method from stata::ppmlhdfe in R

# Resources:
# https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_primer.md

# (Perfect) Separation:
# Let z := Xr be a linear combination of regressors.
# Separation occurs when we can find a z such that
# 1) y > 0 if z = 0 ,and
# 2) y = 0 if z >= 0 with at least one strict inequality

rm(list = ls())

require('data.table')
require('collapse')


# Benchmark data from GH package:
# URL:
# https://github.com/sergiocorreia/ppmlhdfe/tree/master/test/separation_datasets

# Read in separation benchmark data from GH ####
benchmark_li <- lapply(1:18, FUN = function(z) {

  if (z < 10) {
    fread(paste0('https://raw.githubusercontent.com/sergiocorreia/ppmlhdfe/master/test/separation_datasets/0', z, '.csv'))
  }

  else {
    fread(paste0('https://raw.githubusercontent.com/sergiocorreia/ppmlhdfe/master/test/separation_datasets/', z, '.csv'))
  }
})

# Test Separation Datasets from Primer readme: ######
test1 <- data.table(y = c(0, 1, 0, 0, 1),
                   id1 = c(1, 1, 2, 2, 2),
                   id2 = c(1, 1, 1, 2, 2))

test2 <- data.table(y = c(0, 0, 0, 1, 2, 3),
                    id1 = c(2, -1, 0, 0, 5, 6),
                    id2 = c(-1, 2, 0, 0, -10, -12))


test3 <- data.table(y = c(0, 0, 0, 0, 0, 2, 3, 5, 7, 10),
                    x1 = c(0, 0, 0, 0, 1, 21, 0, 0, 0, -18),
                    x2 = c(1, 0, 0, 0, 9, 21, 0, 0, 0, -18),
                    x3 = c(0, 0, 0, 0, 0, 21, 0, 0, 0, 0))


# IR algorithm adopted from Primer readme; translated to R:
ir_sep <- function(data = test1, dep = 1, indep = NULL, fixed = NULL, tol = 1e-5) {

  # Processing user input: #####

  # dependent variable vector 'y':
  if (is.numeric(dep) | is.character(dep)) {
    y <- data[, ..dep]
  } else {
    stop("Unsupported data type for dependent variable:
         dep must be character or numeric vector.")
  }

  # fixed effects matrix 'fes':
  if (is.character(fixed) | is.numeric(fixed)) {
    fes <- data[, ..fixed]
  }
  else if (is.null(fixed)) {
    fes <- NULL
  } else {
    stop("Unsupported format for fixed effects: fixed must be a numeric or character vector or a list
          of numeric or character vectors.")
  }
  # regressor matrix 'X':
  # If there is some form of user input:
  if (is.character(indep) | is.numeric(indep)) {
    x <- data[, ..indep]
  }

  # If there is no kind of user input, take the remaining cols that are neither y nor fixed effects
  else if (is.null(indep)) {
    cat("User did not specify independent variables. By default, the function takes all variables
         not included in 'dep' or 'fixed' as regressors.")

    # If strings were parsed, use boolean indexing
    if (is.character(dep)) dep <- which(names(data) %chin% dep)
    if (is.character(fixed)) fixed <- which(names(data) %chin% fixed)

    # x is Dt without dep and fixed cols:
    x <- data[, !c(..dep, ..fixed)]
  }

  else if (anyNA(indep)) {
    cat("User specified NA for indep; function takes fixed effects variables as regressors")

    # If strings were parsed, use boolean indexing
    #if (is.character(dep)) dep <- which(names(data) %chin% dep)
    # if (is.character(fixed)) fixed <- which(names(data) %chin% fixed)
    #
    # # x is Dt without dep and fixed cols:
    # x <- data[, ..fixed]
  }

   # print('******')
   # print('dependent variable:')
   # print(head(y))
   # print('fixed effects mat')
   # print(head(fes))
   # print('regressor mat')
   # try(print(head(x)))

  # Initialise parameters for algorithm: #######
  u <- as.matrix(ifelse(data[, y] == 0, 1, 0)) # u takes 1 if y is 0 and 0 otherwise
  K <- ceiling(nrow(u) / tol^2) # Weight structure K
  weights <- ifelse(data[, y] > 0, K, 1) # Weight vector to apply weights

  # Iteration counter
  iter <- 0

  # Loop will go on forever, unless it is broken by nested if statement
  while (TRUE) {

    print(paste('iteration', iter))

    # if no regressors except for FE are specified:
    if (anyNA(indep)) {

      # Estimate coefficients with just the FEs
      gamma_hat <- collapse::flm(y = u, X = as.matrix.data.frame(fes), w = weights)

      # predict xg = uhat (calculates the linear prediction from the fitted model) ###
      xg <- as.matrix.data.frame(data[, ..fixed]) %*% gamma_hat

    }

    # In all other cases, where regressors other than FEs are specified:
    else {

      # Partial out FE from regressors:
      #u_cent <- collapse::HDW(x = u, fl = fes, stub = F)
      #x_cent <- collapse::HDW(x = x, fl = fes, stub = F)
      #u_cent <- u; x_cent <- x

      # Estimate with centered data
      gamma_hat <- collapse::flm(y = u, X = as.matrix.data.frame(cbind(x, fes)), w = weights)

      # predict xg = uhat (calculates the linear prediction from the fitted model) ###
      xg <- as.matrix(cbind(data[, ..fixed], x)) %*% gamma_hat

    }

    # replace elements in Xg with zero if abs(Xg) < tol:
    xg[abs(xg) < tol] <- 0

    # break out of loop once all predicted values become positive
    if (all(xg >= 0)) {
      break
    }

    # Update elements in u with help of ReLu function:
    u <- sapply(u, function(k) {
      u[k] <- max(xg[k], 0)
    })
    # update iterator
    iter <- iter + 1
  }

  # Prepare result object
  results <- cbind(data, ifelse(xg > 0, 1, 0))
  names(results) <- c(names(data), 'is_sep')

  print('done')

  return(results)
}

ir_sep(benchmark_li[[4]], dep = 'y', fixed = 2:3, indep = NA)
#ir_sep(benchmark_li[[3]], dep = 'y', fixed = 2:4, indep = 2:4)
#ir_sep(benchmark_li[[1]], dep = 1, fixed = 4:5, indep = 2:3)


### Stata Code from GH Primer website below:
#
# while 1 {
#   qui reghdfe u [fw=w], absorb(id1 id2) resid(e)
#   predict double xb, xbd
#   qui replace xb = 0 if abs(xb) < `tol'
#
# 	* Stop once all predicted values become non-negative
# 	qui cou if xb < 0
# 	if !r(N) {
# 		continue, break
# 	}
#
# 	replace u = max(xb, 0)
# 	drop xb w
# }
#
# rename xb z
# gen is_sep = z > 0
# list y id1 id2 is_sep
