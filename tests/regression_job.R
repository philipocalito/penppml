# This script is designed to run as an Job in the background.

# Objective ###
# Script regresses trade flow on the provision matrix, adding one regressor by time,
# in attempt to narrow down/investigate the error pattern when running the HDFE
# regression.

rm(list = ls())

require(data.table)
require(Rcpp)
require(RcppEigen)
require(RcppArmadillo)

# Workflow:
# 1) open forked Rproj penppml



# 2) execute function script from penppml you want to work on ===============
print(getwd())

source(file = '~/projects/penppml/R/hdfeppml_int.R')
source(file =  '~/projects/penppml//R/utils.R')
sourceCpp(file =  "~/projects/penppml//src/RcppExports.cpp")
source(file = '~/projects/penppml//R/RcppExports.R')
sourceCpp(file =  "~/projects/penppml//src/code.cpp")


### Load my data:
trade <- fread(file = '/home/datahub/projects/Brockhaus_thesis/toy_regressions/output/trade_agg.csv.gz',
               nThread = 999,
               yaml = T,
               encoding = 'UTF-8')

# Remove character col:
trade[, entry_type := NULL]

# Change NAs to zeros
setnafill(x = trade, fill = 0, cols = 6:943)

# Custom hdfe_ppml function to check for errors -------

# create fes:
fes <- list(exp_time = interaction(trade$iso_o, trade$year),
            imp_time = interaction(trade$iso_d, trade$year),
            pair = interaction(trade$iso_o, trade$iso_d))

# Create empty list to store results:
results <- list()

# Loop over subset of regressors:
for (z in 6:length(trade)) {

  message('**** Beginning estimation with ', z, ' regressors ****')
  #print(z)

  # Run regression with modified hdfeppml_int function..
  reg <- try(hdfeppml_int(y = trade$value_usd,
                       x = trade[, 6:z],
                       verbose = F,
                       vcv = T,
                       colcheck = T,
                       fes = fes))

  # Skip if regression fails for some reason:
  if(isTRUE(class(reg) == 'try-error')) {
    next
  }

  # Store model results in case of success:
  else  {
    results[z] <- reg
  }
}

# Determine class of columns
#trade[, lapply(.SD, class)]
