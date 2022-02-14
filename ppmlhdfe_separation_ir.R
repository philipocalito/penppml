### Identify separated observations in poissonhdfe regression ####
# Script implements Iterative Rectifier Method from stata::ppmlhdfe in R

# Resources:
# https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_primer.md
# Undocumented warning about too small values for 'mu'
# https://github.com/sergiocorreia/ppmlhdfe/blob/d3036a45f0d7d63c4cc618c8ded2f37a1d4dbb96/src/ppmlhdfe.mata#L539
# IR-separation coded in mata:
# https://github.com/sergiocorreia/ppmlhdfe/blob/master/src/ppmlhdfe_separation_relu.mata

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


# IR algorithm adopted from Primer readme; translated to R: =====
# algorithm differs in paper and GH readme; sticking to Paper, because it is more recent!

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

  # Loop ==============
  # Will go on forever, unless it is broken by nested if statement
  while (TRUE) {

    print(paste('iteration', iter))

    # if no regressors other than FE are specified:
    if (anyNA(indep)) {

      # weighted LS estimation with just the FEs
      gamma_hat <- collapse::flm(y = u, X = as.matrix.data.frame(fes), w = weights)

      # predict xg = uhat (calculates the linear prediction from the fitted model) ###
      xg <- as.matrix.data.frame(data[, ..fixed]) %*% gamma_hat
      resid <- u - xg

    }

    # In all other cases, where regressors other than FEs are specified:
    else {

      # Partial out FE from regressors and dependent variable:
      u_cent <- collapse::HDW(x = u, fl = fes, stub = F)
      x_cent <- collapse::HDW(x = x, fl = fes, stub = F)

      # weighted LS estimation on centered data & FEs
      gamma_hat <- collapse::flm(y = u_cent, X = as.matrix.data.frame(cbind(x_cent, fes)), w = weights)

      # predict xg = uhat (=linear prediction from the fitted model) ###
      u_hat <- as.matrix(cbind(x_cent, data[, ..fixed])) %*% gamma_hat
      resid <- u_cent - u_hat

    }

     # replace elements in Xg with zero if abs(Xg) < tol:
     #u_hat[abs(u_hat) < tol] <- 0

    # break out of loop once all predicted values become positive
    if (all(abs(resid) < tol)) {
      break
    }

    # Update elements in u with help of ReLu function:
    u_cent[y == 0] <- sapply(u_cent[y == 0], function(z) {
      u_cent[z] <- min(u_hat[z], 0)
    })
    # update iterator
    iter <- iter + 1
  }

  # Prepare result object
  results <- cbind(data, ifelse(xg > 0, 1, 0))
  names(results) <- c(names(data), 'is_sep')

  print('done')

  # For benchmarking examples Check whether separated and is_sep are identical:
  #print(class(results[['is_sep']])); print(class(results[['separated']]))
  print(isTRUE(all.equal(results[['is_sep']], results[['separated']])))

  return(NULL)
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
