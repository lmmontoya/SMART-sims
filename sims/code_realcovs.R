set.seed(123)
source("sims/DGP.R")
library(ltmle)
library(parallel)
library(mvtnorm)
library(matrixStats)
load(file = "/Users/lmontoya/Box/ADAPT Data and Documents/analysis/primary analysis/AT/analysis_obj4/process_data/obj4.RData")
#load(file = "obj4.RData")
obj4$S_missstatus_afterT2 = obj4$S_missstatus_afterT2_ind = NULL
obj4 = subset(obj4[!is.na(obj4$Y_VL),], select = -c(Y_propTIC, Y_VL))

sim_fun = function(x) {
  
  all_data = DGP_realcovs(n = 1692)
  data = all_data$data_ind
  n = nrow(data)
  Anodes = names(data)[grep("A", names(data))]
  
  rule_fun = function(i) {
    
    reslist = list()
    rule = rule.set_ADAPT[i,]
    abar = make_abar_from_rule_ADAPT(rule, data = data)
    rule_bin = unique(data.frame(abar, L_outofcare_first_year = data$L_outofcare_first_year))
    g0 = with(rule.set_ADAPT_bin, data.frame(g.A1soc = mean(unique(rule.set_ADAPT_bin[,c("A1soc", "A1voucher")])$A1soc), 
                                             g.A1voucher = mean(unique(rule.set_ADAPT_bin[,c("A1soc", "A1voucher")])$A1voucher[unique(rule.set_ADAPT_bin[,c("A1soc", "A1voucher")])$A1soc == rule_bin[1,"A1soc"]]),
                                             g.A2navigator = ifelse(data$L_outofcare_first_year == 1, mean(A2navigator[L_outofcare_first_year == 1 & A1soc == rule_bin[1,"A1soc"] & A1voucher == rule_bin[1,"A1voucher"]]), mean(A2navigator[L_outofcare_first_year == 0 & A1soc == rule_bin[1,"A1soc"] & A1voucher == rule_bin[1,"A1voucher"]])), 
                                             g.A2outreach = ifelse(data$L_outofcare_first_year == 1, mean(A2outreach[L_outofcare_first_year == 1 & A1soc == rule_bin[1,"A1soc"] & A1voucher == rule_bin[1,"A1voucher"] & A2navigator == rule_bin[rule_bin$L_outofcare_first_year == 1,"A2navigator"]]), mean(A2outreach[L_outofcare_first_year == 0 & A1soc == rule_bin[1,"A1soc"] & A1voucher == rule_bin[1,"A1voucher"] & A2navigator == rule_bin[rule_bin$L_outofcare_first_year == 0,"A2navigator"]])), 
                                             g.A2sms.voucher = ifelse(data$L_outofcare_first_year == 1, mean(A2sms.voucher[L_outofcare_first_year == 1 & A1soc == rule_bin[1,"A1soc"] & A1voucher == rule_bin[1,"A1voucher"] & A2navigator == rule_bin[rule_bin$L_outofcare_first_year == 1,"A2navigator"] & A2outreach == rule_bin[rule_bin$L_outofcare_first_year == 1,"A2outreach"]]), mean(A2sms.voucher[L_outofcare_first_year == 0 & A1soc == rule_bin[1,"A1soc"] & A1voucher == rule_bin[1,"A1voucher"] & A2navigator == rule_bin[rule_bin$L_outofcare_first_year == 0,"A2navigator"] & A2outreach == rule_bin[rule_bin$L_outofcare_first_year == 0,"A2outreach"]])),
                                             g.A2stop = ifelse(data$L_outofcare_first_year == 1, mean(A2stop[L_outofcare_first_year == 1 & A1soc == rule_bin[1,"A1soc"] & A1voucher == rule_bin[1,"A1voucher"] & A2navigator == rule_bin[rule_bin$L_outofcare_first_year == 1,"A2navigator"] & A2outreach == rule_bin[rule_bin$L_outofcare_first_year == 1,"A2outreach"] & A2sms.voucher == rule_bin[rule_bin$L_outofcare_first_year == 1, "A2sms.voucher"]]), mean(A2stop[L_outofcare_first_year == 0 & A1soc == rule_bin[1,"A1soc"] & A1voucher == rule_bin[1,"A1voucher"] & A2navigator == rule_bin[rule_bin$L_outofcare_first_year == 0,"A2navigator"] & A2outreach == rule_bin[rule_bin$L_outofcare_first_year == 0,"A2outreach"] & A2sms.voucher == rule_bin[rule_bin$L_outofcare_first_year == 0, "A2sms.voucher"]])))) 
    # g0, Qn unadj
    reslist$ltmle1 = suppressMessages(ltmle(data = data[, c("A1soc", "A1voucher", "L_missstatus_beforeT2_ind", "L_outofcare_first_year", "A2navigator", "A2outreach", "A2sms.voucher", "A2stop", "Y")], Anodes = Anodes, Lnodes = c("L_missstatus_beforeT2_ind", "L_outofcare_first_year"), Ynodes = "Y", abar = abar, gform = as.matrix(g0), deterministic.Q.function = det.Q, variance.method = "ic"))
    # gn unadj, Qn unadj
    reslist$ltmle2 = suppressMessages(ltmle(data = data[, c("A1soc", "A1voucher", "L_missstatus_beforeT2_ind", "L_outofcare_first_year", "A2navigator", "A2outreach", "A2sms.voucher", "A2stop", "Y")], Anodes = Anodes, Lnodes = c("L_missstatus_beforeT2_ind", "L_outofcare_first_year"), Ynodes = "Y", abar = abar, gform = gform_ADAPT, deterministic.Q.function = det.Q, deterministic.g.function = det.g, variance.method = "ic"))
    # gn adj, Qn simple
    reslist$ltmle3 = suppressMessages(ltmle(data = data, Anodes = Anodes, Lnodes = colnames(data)[grep("L_", colnames(data))], Ynodes = "Y", deterministic.g.function = det.g, deterministic.Q.function = det.Q, abar = abar, variance.method = "ic"))
    # gn adj, Qn ML - Gcomp
    reslist$ltmle4 = suppressMessages(ltmle(data = data, Anodes = Anodes, Lnodes = colnames(data)[grep("L_", colnames(data))], Ynodes = "Y", deterministic.g.function = det.g, deterministic.Q.function = det.Q, abar = abar, gcomp = T, SL.library = list(Q = default, g = "glm"), variance.method = "ic")) 
    # gn adj, Qn ML - TMLE
    reslist$ltmle5 = suppressMessages(ltmle(data = data, Anodes = Anodes, Lnodes = colnames(data)[grep("L_", colnames(data))], Ynodes = "Y", deterministic.g.function = det.g, deterministic.Q.function = det.Q, abar = abar, SL.library = list(Q = default, g = "glm"), variance.method = "ic")) 
    # gn simple adj, Qn ML - TMLE
    reslist$ltmle6 = suppressMessages(ltmle(data = data, Anodes = Anodes, Lnodes = colnames(data)[grep("L_", colnames(data))], Ynodes = "Y", deterministic.g.function = det.g, deterministic.Q.function = det.Q, abar = abar, SL.library = list(Q = default, g = "glm"), variance.method = "ic", gform = gform_ADAPT)) 
    
    return(reslist)
    
  }

  numrules = nrow(rule.set_ADAPT)
  res = lapply(1:numrules, rule_fun)
  
  # step 1
  ind = sapply(1:numrules, function(j) res[[j]]$ltmle1$cum.g.used[,length(Anodes)])
  for (rule in 1:numrules) { # change me
    ind[all_data$data_fac$A1 == rule.set_ADAPT[rule,1] & data$L_missstatus_beforeT2_ind == 1, rule] = T
  }
  g0 = sapply(1:numrules, function(j) res[[j]]$ltmle1$cum.g[,length(Anodes)])
  H.g0 = ind/g0
  psi_IPTW.g0 = sapply(1:numrules, function(j) res[[j]]$ltmle1$estimates[["iptw"]]*mean(H.g0[,j]))
  IC_IPTW.g0 = sapply(1:numrules, function(j) H.g0[,j]*data$Y - psi_IPTW.g0[j])
 
  # step 2a
  gn.unadj = sapply(1:numrules, function(j) res[[j]]$ltmle2$cum.g[,length(Anodes)])
  H.gn.unadj = ind/gn.unadj
  psi_IPTW.gn.unadj = sapply(1:numrules, function(j) res[[j]]$ltmle2$estimates[["iptw"]]*mean(H.gn.unadj[,j]))
  IC_IPTW.gn.unadj = sapply(1:numrules, function(j) H.gn.unadj[,j]*data$Y - psi_IPTW.gn.unadj[j])
  
  # step 2b
  psi_TMLE.gn.unadj = sapply(1:numrules, function(j) res[[j]]$ltmle2$estimates[["tmle"]])
  IC_TMLE.gn.unadj = sapply(1:numrules, function(j) res[[j]]$ltmle2$IC$tmle)

  # step 3
  gn.adj = sapply(1:numrules, function(j) res[[j]]$ltmle3$cum.g[,length(Anodes)])
  H.gn.adj = ind/gn.adj
  psi_IPTW.gn.adj = sapply(1:numrules, function(j) res[[j]]$ltmle3$estimates[["iptw"]]*mean(H.gn.adj[,j]))
  IC_IPTW.gn.adj = sapply(1:numrules, function(j) H.gn.adj[,j]*data$Y - psi_IPTW.gn.adj[j])
  
  # step 4
  psi_Gcomp.gn.adj = sapply(1:numrules, function(j) res[[j]]$ltmle4$estimates[["gcomp"]])
  IC_Gcomp.gn.adj = sapply(1:numrules, function(j) res[[j]]$ltmle4$IC$gcomp)
  
  # step 5
  psi_TMLE.gn.adj = sapply(1:numrules, function(j) res[[j]]$ltmle5$estimates[["tmle"]])
  IC_TMLE.gn.adj = sapply(1:numrules, function(j) res[[j]]$ltmle5$IC$tmle)
  
  # step 6
  psi_TMLE.gn.simp.adj = sapply(1:numrules, function(j) res[[j]]$ltmle6$estimates[["tmle"]])
  IC_TMLE.gn.simp.adj = sapply(1:numrules, function(j) res[[j]]$ltmle6$IC$tmle)
  

  for_return_fun = function(psi, IC, name) {
    
    names_rules = apply(rule.set_ADAPT, 1, function(x) paste0(x, collapse = "_"))
    names(psi) = colnames(IC) = names_rules
    
    Sigma = cor(IC)
    z = rmvnorm(1e6, sigma = Sigma)
    z_abs = apply(abs(z), 1, max)
    SE_num_simult = quantile(z_abs, .95)
    SE_num_ind = qnorm(0.975) #1.96
    std.dev = sqrt(apply(IC, 2, var)/n)
    CI_simult = psi - SE_num_simult * std.dev 
    CI_simult = cbind(CI_simult, psi + SE_num_simult * std.dev)
    CI_ind = psi - SE_num_ind * std.dev 
    CI_ind = cbind(CI_ind, psi + SE_num_ind * std.dev)
    
    toreturn = list(psi = psi, CI_ind = CI_ind, CI_simult = CI_simult)
    names(toreturn) = paste0(name, ".", names(toreturn))
    return(toreturn)
    
  }
  
  toreturn = c(for_return_fun(psi = psi_IPTW.g0, IC = IC_IPTW.g0, name = "IPTW.g0"),
               for_return_fun(psi = psi_IPTW.gn.unadj, IC = IC_IPTW.gn.unadj, name = "IPTW.gn.unadj"),
               for_return_fun(psi = psi_TMLE.gn.unadj, IC = IC_TMLE.gn.unadj, name = "TMLE.gn.unadj"),
               for_return_fun(psi = psi_IPTW.gn.adj, IC = IC_IPTW.gn.adj, name = "IPTW.gn.adj"),
               for_return_fun(psi = psi_Gcomp.gn.adj, IC = IC_Gcomp.gn.adj, name = "Gcomp.gn.adj"),
               for_return_fun(psi = psi_TMLE.gn.adj, IC = IC_TMLE.gn.adj, name = "TMLE.gn.adj"),
               for_return_fun(psi = psi_TMLE.gn.simp.adj, IC = IC_TMLE.gn.simp.adj, name = "TMLE.gn.simp.adj"))
  print(x)
  return(toreturn)
  
}



B = 1000
results = mclapply(1:B, sim_fun, mc.cores = detectCores())

save(results, file = "simresults_realcovs.RData")
