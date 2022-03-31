make_res = function(est) {
  res = results
  truth = truth$truth
  B = length(results)
  
  psi = do.call('rbind', lapply(1:B, function(x) results[[x]][[paste0(est, ".psi")]]))
  CI_ind_lower = do.call('rbind', lapply(1:B, function(x) results[[x]][[paste0(est, ".CI_ind")]][,1]))
  CI_ind_upper = do.call('rbind', lapply(1:B, function(x) results[[x]][[paste0(est, ".CI_ind")]][,2]))
  CI_simult_lower = do.call('rbind', lapply(1:B, function(x) results[[x]][[paste0(est, ".CI_simult")]][,1]))
  CI_simult_upper = do.call('rbind', lapply(1:B, function(x) results[[x]][[paste0(est, ".CI_simult")]][,2]))
  
  bias = sapply(1:length(truth), function(i) mean(psi[,i] - truth[[i]]))
  var = diag(cov(psi))
  ci_width = colMeans(CI_ind_upper - CI_ind_lower)
  cov_ind = sapply(1:length(truth), function(i) mean(truth[[i]] > CI_ind_lower[,i] & truth[[i]] < CI_ind_upper[,i]))
  simult_cov_simult = mean(rowSums(sapply(1:length(truth), function(i) truth[i] > CI_simult_lower[,i] & truth[i] < CI_simult_upper[,i])) == length(truth))
  simult_cov_ind = mean(rowSums(sapply(1:length(truth), function(i) truth[i] > CI_ind_lower[,i] & truth[i] < CI_ind_upper[,i])) == length(truth))
  
  toreturn = list(bias = bias, var = var, ci_width = ci_width, cov_ind = cov_ind, simult_cov_simult = simult_cov_simult, simult_cov_ind = simult_cov_ind)
  return(toreturn)
}