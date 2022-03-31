set.seed(123)
source("sims/DGP.R")

# simple DGP
mean(DGP_simple(10e6)$data_ind$Y)

# real covs DGP
load(file = "/Users/lmontoya/Box/ADAPT Data and Documents/analysis/primary analysis/AT/analysis_obj4/process_data/obj4.RData")
obj4$S_missstatus_afterT2 = obj4$S_missstatus_afterT2_ind = NULL
obj4 = subset(obj4[!is.na(obj4$Y_VL),], select = -c(Y_propTIC, Y_VL))
mean(DGP_realcovs(10e6)$data_ind$Y)