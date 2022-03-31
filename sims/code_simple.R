set.seed(123)
source("sims/DGP.R")
library(ltmle)
library(parallel)
library(mvtnorm)


sim_fun = function(x) {
  
  data = DGP_simple(n = 1692)$data_ind
  n = nrow(data)
  Anodes = names(data)[grep("A", names(data))]
  
  rule_fun = function(i) {
    
    reslist = list()
    rule = rule.set_simple[i,]
    abar = make_abar_from_rule_simple(rule, data = data)
    rule_bin = unique(data.frame(abar, L2 = data$L2))
    g0 = with(rule.set_simple_bin, data.frame(g.A1 = mean(A1), 
                                              g.A22 = ifelse(data$L2 == 1, mean(A22[L2 == 1 & A1 == rule_bin[1,"A1"]]), mean(A22[L2 == 0 & A1 ==  rule_bin[1,"A1"]])), 
                                              g.A23 = ifelse(data$L2 == 1, mean(A23[L2 == 1 & A1 == rule_bin[1,"A1"] & A22 == rule_bin[rule_bin$L2 == 1,"A22"]]), mean(A23[L2 == 0 & A1 == rule_bin[1,"A1"] & A22 == rule_bin[rule_bin$L2 == 0,"A22"]])), 
                                              g.A24 = ifelse(data$L2 == 1, mean(A24[L2 == 1 & A1 == rule_bin[1,"A1"] & A22 == rule_bin[rule_bin$L2 == 1,"A22"] & A23 == rule_bin[rule_bin$L2 == 1,"A23"]]), mean(A24[L2 == 0 & A1 == rule_bin[1,"A1"] & A22 == rule_bin[rule_bin$L2 == 0,"A22"] & A23 == rule_bin[rule_bin$L2 == 0,"A23"]])))) 
    # g0, Qn unadj
    reslist$ltmle1 = suppressMessages(ltmle(data = subset(data, select = -c(X1, S2)), Anodes = Anodes, Lnodes = "L2", Ynodes = "Y", abar = abar, gform = as.matrix(g0), variance.method = "ic"))
    # gn unadj, Qn unadj
    reslist$ltmle2 = suppressMessages(ltmle(data = subset(data, select = -c(X1, S2)), Anodes = Anodes, Lnodes = "L2", Ynodes = "Y", abar = abar, variance.method = "ic"))
    # gn adj, Qn simple
    reslist$ltmle3 = suppressMessages(ltmle(data = data, Anodes = Anodes, Lnodes = c("L2", "S2"), Ynodes = "Y", abar = abar, variance.method = "ic"))
    # gn adj, Qn ML - Gcomp
    reslist$ltmle4 = suppressMessages(ltmle(data = data, Anodes = Anodes, Lnodes = c("L2", "S2"), Ynodes = "Y", abar = abar, gcomp = T, SL.library = list(Q = default, g = "glm"), variance.method = "ic")) 
    # gn adj, Qn ML - TMLE
    reslist$ltmle5 = suppressMessages(ltmle(data = data, Anodes = Anodes, Lnodes = c("L2", "S2"), Ynodes = "Y", abar = abar, SL.library = list(Q = default, g = "glm"), variance.method = "ic")) 
    
    return(reslist)
    
  }

  numrules = nrow(rule.set_simple)
  res = lapply(1:numrules, rule_fun)
  
  # step 1
  ind = sapply(1:numrules, function(j) res[[j]]$ltmle1$cum.g.used[,length(Anodes)])
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
  

  for_return_fun = function(psi, IC, name) {
    
    names_rules = apply(rule.set_simple, 1, function(x) paste0(x, collapse = "_"))
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
               for_return_fun(psi = psi_TMLE.gn.adj, IC = IC_TMLE.gn.adj, name = "TMLE.gn.adj"))
  print(x)
  return(toreturn)
  
}



B = 1000
results = mclapply(1:B, sim_fun, mc.cores = detectCores())

save(results, file = "simresults_simple.RData")
