### Demo script to get penppml package to work with my data ###

rm(list = ls())

require(data.table)
require(Rcpp)
require(RcppEigen)
require(RcppArmadillo)

# Workflow:
# 1) open forked Rproj penppml

# 2) execute function script from penppml you want to work on ===============
setwd('~/projects/penppml/')

source(file = paste0(getwd(), '/R/hdfeppml_int.R'))
source(file = paste0(getwd(), '/R/utils.R'))
sourceCpp(file = paste0(getwd(), "/src/RcppExports.cpp"))
source(file = paste0(getwd(), '/R/RcppExports.R'))
sourceCpp(file = paste0(getwd(), "/src/code.cpp"))

print('********* Starting regressions: *********')

# 3) load package data and test data set

### Load Package Data:

# selected <- penppml::countries$iso[penppml::countries$region %in% c("Americas")]
# trade2 <- penppml::trade[(penppml::trade$exp %in% selected) & (penppml::trade$imp %in% selected), -(5:6)]
#
# reg1 <- penppml::hdfeppml(data = trade2,
#                           dep = "export",
#                           verbose = T,
#                           colcheck = T,
#                           tol = 1e-7,
#                           fixed = list(c("exp", "time"),
#                                        c("imp", "time"),
#                                        c("exp", "imp")))
#
# # Save results as DF
# results <- data.frame(prov = rownames(reg1$coefficients), b = reg1$coefficients, se = 0)
# results$se[!is.na(reg1$coefficients)] <- reg1$se
# results
#

### Load my data:

trade <- fread(file = '/home/datahub/projects/Brockhaus_thesis/merge_files/temp/merchandise/trade_aggregated.csv.gz',
               nThread = 999,
               yaml = T,
               encoding = 'UTF-8')


# Determine class of columns
#trade[, lapply(.SD, class)]

# Remove character col:
trade[, entry_type := NULL]

setnafill(x = trade, fill = 0, cols = 6:943)

# create fes:
fes <- list(exp_time = interaction(trade$iso_o, trade$year),
            imp_time = interaction(trade$iso_d, trade$year),
            pair = interaction(trade$iso_o, trade$iso_d))

# create subset regressor matrix
#x <- trade[, 6:46]

# Run regression with modified hdfeppml_int function..
reg2 <- hdfeppml_int(y = trade$value_usd,
                     x = trade$PTA,
                     verbose = T,
                     vcv = T,
                     colcheck = T,
                     fes = fes)

# Extract problematic matrix from function
W <- reg2[['W']]
# for choleski decomposition, three conditions have to hold:
# 1) symmetry
isSymmetric.matrix(W) ## TRUE
# 2) real numbers
W <- data.table(W)
W[, lapply(.SD, is.numeric)] |> melt() |> View()

# 3) positive definite
rowSums(W < 0)





# Save results as DF
results <- data.frame(prov = rownames(reg2$coefficients), b = reg2$coefficients, se = 0)
results$se[!is.na(reg2coefficients)] <- reg2$se
results


### END ####
