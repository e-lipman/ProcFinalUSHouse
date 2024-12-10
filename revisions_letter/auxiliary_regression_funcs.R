library(tidyverse)
library(glmnet)
library(ncvreg)
library(BAS)
library(MCMCpack)

source("helper_funs.R")

# helpers
prep_x_covs <- function(inputs){
  x <- as.data.frame(inputs$covariates) %>%
    select_at(configs$covariates) %>%
    mutate_all(standardize) %>%
    as.matrix()
  
  return(x)
}

#-------------------------------------------------------------------------------

# regression wrapper functions
run_penalized_regression <- function(y, inputs){
  family <- ifelse(all(y %in% 0:1), "binomial", "gaussian")
  
  x <- prep_x_covs(inputs)
  
  # run lasso
  cv_lasso <- cv.glmnet(x, y, alpha = 1, family=family)
  lambda_lasso <- cv_lasso$lambda.min
  model_lasso <- glmnet(x, y, alpha = 1, 
                        lambda = lambda_lasso, 
                        family = family)
  coef_lasso <- coef(model_lasso)
  
  # run MCP
  cv_MCP <- cv.ncvreg(x, y, family=family)
  lambda_MCP <- cv_MCP$lambda.min
  model_MCP <- ncvreg(x, y, family=family)
  coef_MCP <- coef(model_MCP)[,which(model_MCP$lambda==lambda_MCP)]
  
  # return results
  res_lasso <- tibble(cov=row.names(coef(model_lasso)), 
                      coef=as.matrix(coef_lasso)[,1],
                      incl_prob=as.numeric(coef!=0),
                      label="LASSO")
  res_MCP <- tibble(cov=row.names(coef(model_MCP)), 
                    coef=as.matrix(coef_MCP)[,1],
                    incl_prob=as.numeric(coef!=0),
                    label="MCP")
  
  bind_rows(res_lasso, res_MCP) %>% filter(cov!="(Intercept)")  
}

run_bayes_model_selection <- function(y, inputs){
  if (all(y %in% 0:1)){
    stop("'run_bayes_model_selection' is only implemented for continuous outcomes")
  } 
  x <- prep_x_covs(inputs)
  
  mod <- bas.lm(y~., data = as.data.frame(x),  
                method="BAS", 
                modelprior = beta.binomial(1, 1),  
                prior="hyper-g-n")
  
  tibble(cov=colnames(x),
         incl_prob=mod$probne0[2:length(mod$probne0)],
         coef=coef(mod)$postmean[2:length(mod$probne0)],
         label="BAS")  
  
}

#-------------------------------------------------------------------------------

# prep data and call regression wrappers

run_stage2_penalized <- function(cong, fldr="stage1", in_path="."){
  cong_str <- paste0("H", str_pad(cong, 3, pad="0"))
  
  # get outcome variable (bridge indicator) from stage 1
  res1 <- readRDS(file.path(in_path, "Results",fldr, "Combined", 
                            paste0(cong_str, "_", "10000_5000_10", 
                                   "_combined.RDS")))
  
  zeta <- as.numeric(1-res1$res$p_change > 0.5)
  
  # get covariates
  inputs <- readRDS(file.path(in_path, "Data",paste0("inputs_", cong_str, ".RDS")))
  
  print(paste0("Running penalized regression for cong ", cong))
  run_penalized_regression(zeta, inputs)
}

run_regression_percent_with_party <- function(cong, in_path,
                                              penalized=T,
                                              bayes=T){
  stopifnot(penalized | bayes)
  
  cong_str <- paste0("H", str_pad(cong, 3, pad="0"))
  inputs <- readRDS(file.path(in_path, "Data",paste0("inputs_", cong_str, ".RDS")))
  
  # make outcome
  J <- ncol(inputs$y)
  dat_with_party <- percent_with_party(inputs) %>%
    mutate_at(c("proc_with_party","final_with_party"),
              ~ifelse(.==1, (J-1)/J, .)) %>%
    mutate(odds_final = final_with_party/(1-final_with_party),
           odds_proc = proc_with_party/(1-proc_with_party),
           log_OR = log(odds_final)-log(odds_proc))
  
  # run regression
  out <- tibble()
  if (penalized){
    print(paste0("Running penalized regression for cong ", cong))
    res_penalized <- run_penalized_regression(dat_with_party$log_OR, inputs)  
    out <- bind_rows(out, res_penalized)
  }
  if (bayes){
    print(paste0("Running BAS for cong ", cong))
    res_bayes <- run_bayes_model_selection(dat_with_party$log_OR, inputs)    
    out <- bind_rows(out, res_bayes)
  }
  return(out)
}

#-------------------------------------------------------------------------------

run_chain_mcmcpack <- function(y, constraint,
                               burn, iter,
                               load_saved = F,
                               save_filepath=NULL){
  print(paste0("Running ", burn, " burn-in and ", iter, " iterations"))
  
  if (load_saved){
    mcmc_run <- readRDS(save_filepath)
  } else {
    mcmc_run <- MCMCirt1d(
      datamatrix = y, 
      theta.constraints = constraint,
      T0 = 1/5, AB0 = 1/5, burnin = burn, mcmc = iter, verbose = 1000, 
      store.item=F, seed=24680)  
    if (!is.null(save_filepath)){
      saveRDS(mcmc_run, save_filepath)  
    }
  }
  
  
  post_means <- summary(mcmc_run)$statistics[,1]
  stopifnot(all(rownames(y)==str_remove(names(post_means),"theta.")))
  
  return(post_means)
}

# Run proc and final passage votes separately
run_proc_final <- function(cong, constraint,
                           in_path=".",
                           save_path="auxiliary_res",
                           load_saved=F,
                           burn=10000, iter=25000){
  
  cong_str <- paste0("H", str_pad(cong, 3, pad="0"))
  inputs <- readRDS(file.path(in_path, "Data",paste0("inputs_", cong_str, ".RDS")))
  
  # set constraint
  repLeader <- inputs$members$Member[inputs$repLeaderC+1]
  constraint <- list("+")
  names(constraint) <- repLeader
  
  # run 
  rownames(inputs$y) <- inputs$members$Member 
  y_proc <- inputs$y[,1:(inputs$gam0_maxC+1)]
  y_final <- inputs$y[,(inputs$gam0_maxC+2):ncol(inputs$y)]
  
  proc_filename <- paste0(cong_str, "_proc_mcmcpack.RDS")
  theta_proc <- run_chain_mcmcpack(y_proc, constraint,
                                   burn=burn, iter=iter,
                                   load_saved = load_saved,
                                   file.path(save_path, proc_filename))
  
  final_filename <- paste0(cong_str, "_final_mcmcpack.RDS")
  theta_final <- run_chain_mcmcpack(y_final, constraint,
                                    burn, iter, 
                                    load_saved = load_saved,
                                    file.path(save_path, final_filename))
  
  # output
  out <- tibble(member=inputs$members$Member, 
                party=inputs$members$party, 
                beta0=theta_proc, 
                beta1=theta_final,
                rank_beta0=rank(beta0),
                rank_beta1=rank(beta1))
  
  out_filename <- paste0(cong_str, "_summy_mcmcpack.RDS")
  saveRDS(out, file.path(save_path, out_filename))
}








