#' PPML Estimation with HDFE
#'
#' \code{hdfeppml_int} is the internal algorithm called by \code{hdfeppml} to fit an (unpenalized)
#' Poisson Pseudo Maximum Likelihood (PPML) regression with high-dimensional fixed effects (HDFE). It
#' takes a vector with the dependent variable, a regressor matrix and a set of fixed effects (in list
#' form: each element in the list should be a separate HDFE).
#'
#' More formally, \code{hdfeppml_int} performs iteratively re-weighted least squares (IRLS) on a
#' transformed model, as described in Correia, Guimarães and Zylkin (2020) and similar to the
#' \code{ppmlhdfe} package in Stata. In each iteration, the function calculates the transformed dependent
#' variable, partials out the fixed effects (calling \code{collapse::fhdwithin}, which uses the algorithm in
#' Gaure (2013)) and then solves a weighted least squares problem (using fast C++ implementation).
#'
#' @param y Dependent variable (a vector)
#' @param x Regressor matrix.
#' @param fes List of fixed effects.
#' @param tol Tolerance parameter for convergence of the IRLS algorithm.
#' @param hdfetol Tolerance parameter for the within-transformation step,
#'     passed on to \code{collapse::fhdwithin}.
#' @param colcheck Logical. If \code{TRUE}, checks for perfect multicollinearity in \code{x}.
#' @param mu Optional: initial values of the conditional mean \eqn{\mu}, to be used as weights in the
#'     first iteration of the algorithm.
#' @param saveX Logical. If \code{TRUE}, it returns the values of x and z after partialling out the
#'     fixed effects.
#' @param init_z Optional: initial values of the transformed dependent variable, to be used in the
#'     first iteration of the algorithm.
#' @param verbose Logical. If \code{TRUE}, it prints information to the screen while evaluating.
#' @param maxiter Maximum number of iterations (a number).
#' @param cluster Optional: a vector classifying observations into clusters (to use when calculating SEs).
#' @param vcv Logical. If \code{TRUE} (the default), it returns standard errors.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{coefficients}: a 1 x \code{ncol(x)} matrix with coefficient (beta) estimates.
#'   \item \code{residuals}: a 1 x \code{length(y)} matrix with the residuals of the model.
#'   \item \code{mu}: a 1 x \code{length(y)} matrix with the final values of the conditional mean \eqn{\mu}.
#'   \item \code{deviance}:
#'   \item \code{bic}: Bayesian Information Criterion.
#'   \item \code{x_resid}: matrix of demeaned regressors.
#'   \item \code{z_resid}: vector of demeaned (transformed) dependent variable.
#'   \item \code{se}: standard errors of the coefficients.
#' }
#'
#' @export
#'
#' @examples
#' # To reduce run time, we keep only countries in the Americas:
#' americas <- countries$iso[countries$region == "Americas"]
#' trade <- trade[(trade$imp %in% americas) & (trade$exp %in% americas), ]
#' # Now generate the needed x, y and fes objects:
#' y <- trade$export
#' x <- data.matrix(trade[, -1:-6])
#' fes <- list(exp_time = interaction(trade$exp, trade$time),
#'             imp_time = interaction(trade$imp, trade$time),
#'             pair     = interaction(trade$exp, trade$imp))
#' # Finally, the call to hdfeppml_int:
#' reg <- hdfeppml_int(y = y, x = x, fes = fes)
#'
#' @section References:
#' Breinlich, H., Corradi, V., Rocha, N., Ruta, M., Santos Silva, J.M.C. and T. Zylkin (2021).
#' "Machine Learning in International Trade Research: Evaluating the Impact of Trade Agreements",
#' Policy Research Working Paper; No. 9629. World Bank, Washington, DC.
#'
#' Correia, S., P. Guimaraes and T. Zylkin (2020). "Fast Poisson estimation with high dimensional
#' fixed effects", \emph{STATA Journal}, 20, 90-115.
#'
#' Gaure, S (2013). "OLS with multiple high dimensional category variables",
#' \emph{Computational Statistics & Data Analysis}, 66, 8-18.
#'
#' Friedman, J., T. Hastie, and R. Tibshirani (2010). "Regularization paths for generalized linear
#' models via coordinate descent", \emph{Journal of Statistical Software}, 33, 1-22.
#'
#' Belloni, A., V. Chernozhukov, C. Hansen and D. Kozbur (2016). "Inference in high dimensional panel
#' models with an application to gun control", \emph{Journal of Business & Economic Statistics}, 34, 590-605.
#'

hdfeppml_int <- function(y, x=NULL, fes=NULL, tol = 1e-8, hdfetol = 1e-4, colcheck_x = TRUE,
                             colcheck_x_fes = TRUE, mu = NULL, saveX = TRUE,
                             init_z = NULL, verbose = FALSE, maxiter = 1000, cluster = NULL, vcv = TRUE) {
  if(missing(x) & missing(fes)){
    stop("Please provide at least one of the arguments x or fes.")
  }
  if(!missing(x)){
    x <- data.matrix(x)
  }

  #if(missing(fes)){
  #  x <- data.matrix(cbind(1,x))
  #}

  # number of observations (needed for deviance)
  n <- length(y)

  # estimation algorithm

  ### initialising some key variables used later #########
  crit <- 1
  iter <- 0
  old_deviance <- 0

  if(!missing(x)){
    include_x <- 1:ncol(x)

    b <- matrix(NA, nrow = ncol(x), ncol = 1)
    xnames <- colnames(x)

    # The collinearity check is only relevant if we have x in the model
    if (colcheck_x == TRUE | colcheck_x_fes == TRUE){
      if (verbose == TRUE) {
        print("checking collinearity")
      }
      if(!missing(fes)){
        include_x <- collinearity_check(y=y, x=x, fes=fes, hdfetol=1e-6, colcheck_x = colcheck_x, colcheck_x_fes = colcheck_x_fes)
        x <- x[, include_x]
      } else {
        print("do collinearity check")
        include_x <- collinearity_check(y=y, x=x, hdfetol=1e-6, colcheck_x = colcheck_x, colcheck_x_fes = FALSE)
        print(include_x)
        x <- x[, include_x]
      }
    }

  }

    if (verbose == TRUE) {
      print('colcheck successful, regressor matrix updated')
    }

  }
  # IRLS algorithm ------------------------------------
  if (verbose == TRUE) {
    print("beginning IRLS estimation")
  }

  while (crit>tol & iter<maxiter) {
    iter <- iter + 1

    if (verbose == TRUE) {
      print(paste('-- Iteration: ', iter, '--'))
    }
    # If it is the first iteration, initialise mu and z ================
    if (iter == 1) {
      ## initialize "mu" ########
      if (is.null(mu)) {
          if (verbose == TRUE) {
            print(paste('user did not supply initial value for mu; initialising mu'))
          }

        mu  <- (y + mean(y))/2
      }

      # computing values for z and eta from mu #########
      print(paste('Initialising z and eta'))
      z   <- (y-mu)/mu + log(mu) # transformed dependent variable
      eta <- log(mu) # eta is initialised but never used in code???
      last_z <- z

      # initialise z ########
      if (is.null(init_z)) {

        print(paste('user did not supply initial value for z; initialising z from mu and y'))

        reg_z  <- matrix(z)
      }

      else {
        reg_z <- init_z
      }
      if(!missing(x)){
        # copying regressor matrix to internal object
        reg_x  <- x
      }

    # Do this in all iterations, except for first =================
    else {
      last_z <- z
      z <- (y-mu)/mu + log(mu) # transformed dependent variable
      # objects used in regression:
      reg_z  <- matrix(z - last_z + z_resid)
      if(!missing(x)){
        reg_x  <- x_resid
      }
      ## colnames(reg_x)   <- colnames(x)
    }
    if (verbose == TRUE) {
      print("partialing out fixed effects, weighting with mu")
      print('this yields residuals for z and x')
    }

    if(!missing(fes)){
      if(is.null(fes)){
        # FWL theorem; CGZ(2020), page. 100. These residuals should be equal to residuals from IRLS
        z_resid <- collapse::fwithin(x=reg_z, g=factor(rep(1,length(reg_z))), w = mu)
        if(!missing(x)){
          x_resid <- collapse::fwithin(x=reg_x,g=factor(rep(1,length(reg_z))), w = mu)
        }
      }else{
      z_resid <-  collapse::fhdwithin(reg_z, fes, w = mu)
      if(!missing(x)){
        x_resid <- collapse::fhdwithin(reg_x, fes, w = mu)
        }
      }
    } else {
        z_resid <- reg_z
        if(!missing(x)){
        x_resid <- reg_x
      }
    }

    if (verbose == TRUE) {
      print("-- obtaining coefficients --")
    }

    # Calculating model residuals ###########
    # Standard formula e = y - X'ß applied here...

    if(!missing(x)){
      # Regression with within transformed versions of x and z, weighted with sqrt(mu)
      reg <- fastolsCpp(sqrt(mu) * x_resid, sqrt(mu) * z_resid)  #faster without constant
      b[include_x] <- reg
      reg <- list("coefficients" = b) # this rewrites reg each time. Better to initialize upfront?

      # If include_x is a scalar:
      if (length(include_x) == 1) {
        reg$residuals <- z_resid - x_resid * b[include_x]
        # If include_x is a matrix:
      } else {
        reg$residuals <- z_resid - x_resid %*% b[include_x]
      }
    } else {
      reg <- NULL
      reg$residuals <- z_resid
    }

    # Updating mu #######
    # FWL theorem implies that z - reg$residuals is the same as X'ß
    mu <- as.numeric(exp(z - reg$residuals))
    mu[which(mu < 1e-190)] <- 1e-190
    mu[mu > 1e190] <- 1e190

    # Some info for user if verbose is T:
    if (verbose == TRUE) {
      print("-- info on residuals --")
      print(paste('larest residual: ', max(reg$residuals)))
      print(paste('smallest residual: ', min(reg$residuals)))

      print("-- info on means --")
      print(paste('maximum mu: ', max(mu)))
      print(paste('minimum mu: ', min(mu)))

      print("-- info on coefficients --")
      print(paste('largest coefficient: ', max(b[include_x])))
      print(paste('smallest coefficient: ', min(b[include_x])))
    }

    ### Checking for Convergence =====================

    ### calculate deviance ###########################

    if (verbose == TRUE) {
      print("-- calculating deviance --")
    }

    # what is temp?
    temp <-  -(y * log(y/mu) - (y-mu))

    # Update temp at the index where y is zero (what does that do=)
    temp[which(y == 0)] <- -mu[which(y == 0)]
    # print("How many temp NA:")
    # print(which(is.na(temp)))
    #temp[which(is.na(temp))] <- 0

    deviance <- -2 * sum(temp) / n

    # print("hdfe")
    # print(temp[which(is.na(temp))])
    # print(mu[which(is.na(temp))])
    # print(y[which(is.na(temp))])

    if (deviance < 0) deviance = 0

      if (verbose == TRUE) print(paste('deviance was less than 0; set new deviance to 0'))
    }

    # Calculate difference in deviance
    delta_deviance <- old_deviance - deviance

    # Set delta deviance to deviance, if delta deviance is non missing and deviance is less than 10% of delta_deviance
    if (!is.na(delta_deviance) & (deviance < 0.1 * delta_deviance)) {
      if (verbose == TRUE) {
        print(paste('delta deviance is non missing; deviance is less than 10% of delta deviance: \n
                    Setting delta_deviance to deviance value: ', deviance))
      }
      # Updating delta_deviance
      delta_deviance <- deviance
    }

    # Checking critical value #################################
    if (verbose == TRUE) {
      print("-- checking critical value --")
    }

    # compare old_deviance and deviance; take minimum,
    # then, compare minimum deviance with 0.1; take maximum
    # Is denom_crit a ReLU? (https://en.wikipedia.org/wiki/Rectifier_(neural_networks))
    denom_crit = max(c(min(c(deviance, old_deviance)), 0.1))
    crit = abs(delta_deviance) / denom_crit
    if (verbose == TRUE) {
      print('critical value is a fraction')
      print(paste('numerator is: ', abs(delta_deviance), '; denominator is: ', denom_crit))
      print(paste('...resulting in critical value: ', crit))
      print(paste('current deviance is: ', deviance))
      print(paste('as long as', crit, 'is smaller than', tol, ', loop will continue.'))
    }
    old_deviance <- deviance
    # While-loop ends here ###############
  }


  temp <-  -(y * log(y / mu) - (y - mu))
  temp[which(y == 0)] <- 0

  if (verbose == TRUE) {
    print('while loop completed')
    print("-- converged --")
  }

  ## elements to return
  if(!missing(x)){
    k   <- ncol(matrix(x))
    n   <- length(y)
    reg$mu  <- mu
    reg$deviance <- -2 * sum(temp) / n
    reg$bic <- deviance + k * log(n) / n
    rownames(reg$coefficients) <- xnames
  } else {
    reg$mu  <- mu
  }

  # k = number of elements in x here
  # BIC would be BIC = deviance + k * ln(n)

  #returnlist <- list("coefficients" = b, "mu" = mu, "bic" = bic, "deviance" = deviance)
  if (saveX == TRUE) {
    if(!missing(x)){
      reg[["x_resid"]] <- x_resid
    }
    reg[["z_resid"]] <- z_resid
  }
  # Computing Standard errors --------------------------
  if(!missing(x)){
    if (vcv) {
      # If cluster has been parsed:
      if(!is.null(cluster)) {
        nclusters  <- nlevels(droplevels(cluster, exclude = if(anyNA(levels(cluster))) NULL else NA))
        het_matrix <- (1 / nclusters) * cluster_matrix((y - mu) / sum(sqrt(mu)), cluster, x_resid)
        W          <- (1/nclusters) * (t(mu*x_resid) %*% x_resid) / sum(sqrt(mu))
        R <- try(chol(W), silent = FALSE)
        V          <- (1/nclusters) * chol2inv(R) %*% het_matrix %*% chol2inv(R)
        #V          <- (1/nclusters)*solve(W)%*%het_matrix%*%solve(W)
        V          <- nclusters / (nclusters - 1) * V

        # If cluster hasn't been parsed:
      } else {
        e = y - mu
        het_matrix = (1/n) * t(x_resid*e)  %*% (x_resid*e)
        W          = (1/n) * (t(mu*x_resid) %*% x_resid)
        message('****** Begin debug on W (R) Matrix ******')
        print(paste('class of W: ', class(W)))


        #print(is.na(W), is.nan(W))
        #print(W)
        #print(class(R))

        reg[['W']] <- W
        message('****** End debug on W (R) Matrix ******')
        R          = try(chol(W), silent = TRUE)
        V          = (1/n) * chol2inv(R) %*% het_matrix %*% chol2inv(R)
        V          = (n / (n - 1)) * V
      }
    }
    reg[["se"]] <- sqrt(diag(V))
  }
  return(reg)
}
