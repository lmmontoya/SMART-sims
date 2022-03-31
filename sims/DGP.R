DGP_simple = function(n, a1 = NULL, a2.if.fail = NULL, a2.if.succ = NULL) {
  
  # baseline covariates
  X1 = rnorm(n)
  
  # stage 1 treatment
  if (is.null(a1)) {
    A1 = rbinom(n, 1, 0.5)
  } else {
    A1 = rep(as.numeric(a1), n)
  }
  
  # time varying covariates
  L2 = rbinom(n, 1, plogis(X1 + A1))
  S2 = rnorm(n, X1 + 2*A1)
  
  # stage 2 treatment
  A2 = rep(NA, n)
  if (is.null(a2.if.fail) & is.null(a2.if.succ)) {
    A2[L2 == 1] = sample(c(1,2), sum(L2 == 1), replace = T)
    A2[L2 == 0] = sample(c(3,4), sum(L2 == 0), replace = T)
  } else {
    A2[L2 == 1] = a2.if.fail
    A2[L2 == 0] = a2.if.succ
  }
  A2 = factor(A2, levels = c("1", "2", "3", "4"))
  
  # outcome
  EY = expand.grid(first.line = c("0", "1"), second.line = c("1", "2", "3", "4"))
  EY$truth = 1 - c(.28, .26, .28, .3, .29, .3, .21, .2)
  prob.Y = rowSums(sapply(1:dim(EY)[1], function(x) as.numeric(A1 == EY[x,"first.line"] & A2 == EY[x,"second.line"])*EY[x,"truth"]))
  qlogis.prob.Y = qlogis(prob.Y) + S2 + 0.5*(X1^2) + log(abs(X1)+.01)
  prob.Y = plogis(qlogis.prob.Y)
  Y = rbinom(n, 1, prob = prob.Y)
  
  data_fac = data.frame(X1, A1, L2, S2, A2, Y)
  toreturn = data_fac
  
  if (is.null(a1)) {
    options(na.action='na.pass')
    data_ind = data.frame(model.matrix(~ ., data = data_fac, na.action = "na.fail")[,-1])
    toreturn = list(data_fac = data_fac, data_ind = data_ind)
  }
  
  return(toreturn)
  
}

DGP_realcovs = function(n, a1 = NULL, a2.if.fail = NULL, a2.if.succ = NULL) {
  
  n = nrow(obj4)
  obj4_resamp = obj4[sample(n, n, replace = T),]
  
  W = obj4_resamp[,grep("W_", colnames(obj4_resamp))]
  W_CD4 = W$W_CD4
  W_sex = W$W_sex
  W_alcohol = W$W_alcohol
  W_whostage = W$W_whostage
  W_age = W$W_age
  
  L = obj4_resamp[,grep("L_", colnames(obj4_resamp))]
  L_missstatus_beforeT2 = L$L_missstatus_beforeT2
  L_missstatus_beforeT2_ind = L$L_missstatus_beforeT2_ind
  L_outofcare_first_year = L$L_outofcare_first_year
  L_time_rand_to_T2 = L$L_time_rand_to_T2
  
  #### A1 ####
  if (is.null(a1)) {
    A1 = sample(c("sms", "voucher", "soc"), size = n, replace = T)
  } else {
    A1 = rep(a1, n)
  }
  
  #### A2 ####
  A2 = rep(NA, n)
  if (is.null(a2.if.fail) & is.null(a2.if.succ)) {
    A2[L_missstatus_beforeT2_ind == 0 & L_outofcare_first_year == 1] = sample(c("sms.voucher", "outreach", "navigator"), size = sum(L_missstatus_beforeT2_ind == 0 & L_outofcare_first_year == 1), replace = T)
    A2[L_missstatus_beforeT2_ind == 0 & L_outofcare_first_year == 0 & A1 %in% c("sms", "voucher")] = sample(c("stop", "continue"), size = sum(L_missstatus_beforeT2_ind == 0 & L_outofcare_first_year == 0 & A1 %in% c("sms", "voucher")), replace = T)
    A2[L_missstatus_beforeT2_ind == 0 & L_outofcare_first_year == 0 & A1 == "soc"] = "continue"
  } else {
    A2[L_missstatus_beforeT2_ind == 0 & L_outofcare_first_year == 1] = rep(a2.if.fail, sum(L_missstatus_beforeT2_ind == 0 & L_outofcare_first_year == 1))
    A2[L_missstatus_beforeT2_ind == 0 & L_outofcare_first_year == 0] = rep(a2.if.succ, sum(L_missstatus_beforeT2_ind == 0 & L_outofcare_first_year == 0))
  }
  
  
  #### Y ####
  EY = expand.grid(first.line = c("sms", "soc", "voucher"), second.line = c("sms.voucher", "navigator", "outreach", "stop", "continue"))
  EY = EY[-which(EY$first.line == "soc" & EY$second.line == "stop"),]
  EY$truth = 1 - c(.28, .26, .28, .3, .29, .3, .21, .2, .21, .18, .18, .22, .13, .22)
  prob.Y = rowSums(sapply(1:dim(EY)[1], function(x) as.numeric(A1 == EY[x,"first.line"] & A2 == EY[x,"second.line"])*EY[x,"truth"]))
  qlogis.prob.Y = qlogis(prob.Y) + L_outofcare_first_year + L_time_rand_to_T2/300 - as.numeric(W_sex == "M")*(W_age/10)
  prob.Y = plogis(qlogis.prob.Y) 
  prob.Y[is.na(prob.Y)] = (.35 + (W_alcohol == "level 1 & 2")*.1)[is.na(prob.Y)]
  Y = rbinom(n, 1, prob = prob.Y)
  Y[L_missstatus_beforeT2 == "death before T2"] = 0
  
  data_fac = data.frame(W, A1 = factor(A1), L, A2 = factor(A2), Y)
  toreturn = data_fac
  
  if (is.null(a1)) {
    options(na.action='na.pass')
    data_ind = data.frame(model.matrix(~ ., data = data_fac, na.action = "na.fail")[,-1])
    toreturn = list(data_fac = data_fac, data_ind = data_ind)
  }
  
  return(toreturn)
  
}

make_abar_from_rule_simple = function(rule, data) {
  
  n = nrow(data)
  Anodes = names(data)[grep("A", names(data))]
  abar = matrix(0, nrow = n, ncol = length(Anodes), dimnames = list(NULL, Anodes))
  
  # set A1
  abar[,"A1"] = as.numeric(as.character(rule[["a1"]]))
  
  # set A2
  fail = data$L2 == 1
  abar[fail, grep(paste0("A2", rule[["a2.if.fail"]]), Anodes)] = 1
  abar[!fail, grep(paste0("A2", rule[["a2.if.succ"]]), Anodes)] = 1
  
  return(abar)
  
}

make_abar_from_rule_ADAPT = function(rule, data) {
  
  n = nrow(data)
  Anodes = names(data)[grep("A", names(data))]
  abar = matrix(0, nrow = n, ncol = length(Anodes), dimnames = list(NULL, Anodes))
  
  # set A1
  abar[,grep(paste0("A1", rule[["a1"]]), Anodes)] = 1
  
  # set A2
  fail = data$L_outofcare_first_year == 1
  abar[fail, grep(paste0("A2", rule[["a2.if.fail"]]), Anodes)] = 1
  abar[!fail, grep(paste0("A2", rule[["a2.if.succ"]]), Anodes)] = 1
  
  return(abar)
  
}

rule.set_simple = expand.grid(a1 = c("0", "1"), a2.if.fail = c("1", "2"), a2.if.succ = c("3", "4"))
rule.set_simple_bin = data.frame(a1 = c(rule.set_simple$a1, rule.set_simple$a1),
                                L2 = rep(c(1,0), each = nrow(rule.set_simple)),
                                a2 = c(rule.set_simple$a2.if.fail, rule.set_simple$a2.if.succ))
rule.set_simple_bin = unique(rule.set_simple_bin)
rule.set_simple_bin = data.frame(A1 = as.numeric(rule.set_simple_bin$a1 == 1),
                                L2 = rule.set_simple_bin$L2,
                                A22 = as.numeric(rule.set_simple_bin$a2 == 2),
                                A23 = as.numeric(rule.set_simple_bin$a2 == 3),
                                A24 = as.numeric(rule.set_simple_bin$a2 == 4))

rule.set_ADAPT = expand.grid(a1 = c("soc", "sms", "voucher"), a2.if.fail = c("outreach", "sms.voucher", "navigator"), a2.if.succ = c("continue", "stop"))
rule.set_ADAPT = rule.set_ADAPT[-which(rule.set_ADAPT$a1 == "soc" & rule.set_ADAPT$a2.if.succ == "stop"),]
rule.set_ADAPT_bin = data.frame(a1 = c(rule.set_ADAPT$a1, rule.set_ADAPT$a1),
                                L_outofcare_first_year = rep(c(1,0), each = nrow(rule.set_ADAPT)),
                                a2 = c(rule.set_ADAPT$a2.if.fail, rule.set_ADAPT$a2.if.succ))
rule.set_ADAPT_bin = unique(rule.set_ADAPT_bin)
rule.set_ADAPT_bin = data.frame(A1soc = as.numeric(rule.set_ADAPT_bin$a1 == "soc"),
                                A1voucher = as.numeric(rule.set_ADAPT_bin$a1 == "voucher"),
                                L_outofcare_first_year = rule.set_ADAPT_bin$L_outofcare_first_year,
                                A2navigator = as.numeric(rule.set_ADAPT_bin$a2 == "navigator"),
                                A2outreach = as.numeric(rule.set_ADAPT_bin$a2 == "outreach"),
                                A2sms.voucher = as.numeric(rule.set_ADAPT_bin$a2 == "sms.voucher"),
                                A2stop = as.numeric(rule.set_ADAPT_bin$a2 == "stop"))

default = list("SL.glm", "SL.stepAIC", "SL.bayesglm", c("SL.glm", "screen.corP"), c("SL.step", "screen.corP"), c("SL.step.forward", "screen.corP"), c("SL.stepAIC", "screen.corP"), c("SL.step.interaction", "screen.corP"), c("SL.bayesglm", "screen.corP"))

det.g <- function (data, current.node, nodes) {
  data.names <- names(data)
  if (grepl("A1", data.names[current.node])) {
    str <- "A1"
  } else if (grepl("A2", data.names[current.node])) {
    str <- "A2"
  } else {
    stop("unexpected error")
  }
  Anodes <- nodes$A
  if (!(any(Anodes < current.node))) return(NULL)
  prev.Anode.names <- data.names[Anodes[Anodes < current.node]]
  prev.Anode.names <- grep(str, prev.Anode.names, value=TRUE)
  if (length(prev.Anode.names) < 1) return(NULL)
  
  is.deterministic <- rowAnys(data[, prev.Anode.names, drop=FALSE] == 1 & !is.na(data[, prev.Anode.names, drop=FALSE]))
  
  current.node.name <- data.names[current.node]
  if (current.node.name == "A2stop") {
    is.deterministic[data$L_outofcare_first_year %in% 1] <- T
  }
  return(list(is.deterministic = is.deterministic, prob1 = 0))
}

det.Q <- function (data, current.node, nodes, called.from.estimate.g) {
  L_missstatus_beforeT2_ind.index <- which(names(data) == "L_missstatus_beforeT2_ind")
  stopifnot(length(L_missstatus_beforeT2_ind.index) == 1)
  L_missstatus_beforeT2_ind.in.history <- L_missstatus_beforeT2_ind.index < current.node
  
  if (!L_missstatus_beforeT2_ind.in.history) return(NULL)
  
  is.deterministic <- data$L_missstatus_beforeT2_ind == 1
  Q.value <- data$Y[is.deterministic]
  return(list(is.deterministic = is.deterministic, Q.value = Q.value))
}

gform_ADAPT = c("A1soc ~ 1", 
                "A1voucher ~ A1soc", 
                "A2navigator ~ L_outofcare_first_year + A1soc + A1voucher",
                "A2outreach ~ A2navigator + L_outofcare_first_year + A1soc + A1voucher",
                "A2sms.voucher ~ A2navigator + A2outreach + L_outofcare_first_year + A1soc + A1voucher",
                "A2stop ~ A2navigator + A2outreach + A2sms.voucher + L_outofcare_first_year + A1soc + A1voucher")


calctruth = F
if (calctruth) {
  
  library(parallel)
  set.seed(123)
  load(file = "/Users/lmontoya/Box/ADAPT Data and Documents/analysis/primary analysis/AT/analysis_obj4/process_data/obj4.RData")
  obj4$S_missstatus_afterT2 = obj4$S_missstatus_afterT2_ind = NULL
  obj4 = subset(obj4[!is.na(obj4$Y_VL),], select = -c(Y_propTIC, Y_VL))
  
  get_truth = function(DGP_fun, rule.set, iter) { 
    
    get_mean = function(x) {
      return(unlist(lapply(1:dim(rule.set)[1], function(i) mean(DGP_fun(n = 10e4, a1 = as.character(rule.set[i,"a1"]), a2.if.fail = as.character(rule.set[i,"a2.if.fail"]), a2.if.succ = as.character(rule.set[i,"a2.if.succ"]))$Y))))
    }
    means = mclapply(1:iter, function(x) get_mean(x), mc.cores = detectCores())
    means = rowMeans(do.call('cbind', means))
    rule.set$truth = means
    return(rule.set)
    
  }
  
  # simple
  truth_simple = get_truth(DGP_fun = DGP_simple, rule.set = rule.set_simple, iter = 10e2)
  print("truth_simple done")
  # real covs
  truth_realcovs = get_truth(DGP_fun = DGP_realcovs, rule.set = rule.set_ADAPT, iter = 10e5)
  print("truth_realcovs done")
  
  save(truth_simple, truth_realcovs, file = "truths.RData")
  
}