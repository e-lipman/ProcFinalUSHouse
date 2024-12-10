library(tidyverse)
library(yaml)

configs <- read_yaml("configs.yml")
keep_chains <- configs$keep_chains

combine_chains <- function(out_raw, party, keep=NULL,
                           keep_separate = c("theta", "lpost",
                                             "intercept","model_size",
                                             "zeta","g"),
                           sep_and_comb = "zeta"){
  if(!"list" %in% class(out_raw[[1]])){ 
    stop("'out_raw' must be a list of lists") 
  }
  if (is.null(keep)){keep=1:length(out_raw)}
  
  out <- list(keep=keep)
  for (arg in names(out_raw[[1]])){
    chain_list <- map(out_raw, ~(.[[arg]]))
    if (arg %in% keep_separate){ 
      arg_name = paste0(arg,"_chains")  
      out[[arg_name]] <- chain_list
    }
    if (!arg %in% keep_separate | arg %in% sep_and_comb){
      if ("matrix" %in% class(out_raw[[1]][[arg]])){
        out[[arg]] <- reduce(chain_list[keep], rbind)
      } else if ("numeric" %in% class(out_raw[[1]][[arg]])){
        out[[arg]] <- unlist(chain_list[keep])
      }
    }
  }
  
  out$zetaD_chains <- map(out$zeta_chains, ~rowMeans(.x[,party=="D"]))
  out$zetaR_chains <- map(out$zeta_chains, ~rowMeans(.x[,party=="R"]))
  out$zeta_chains <- map(out$zeta_chains, rowMeans)
  
  if ("eta" %in% names(out)){
    out$eps = as.numeric(out$eta!=0)
  }
  return(out)
}

postprocess <- function(cong_num, folder = ".", suffix = "10000_5000_10",
                        chains = 1:4, exclude = T, stage=0, 
                        retain_zeta=T,
                        bare_bones=F){
  
  topdir <- "."
  cong_str <- paste0("H",str_pad(cong_num, 3, pad=0))
  print(paste0("Postprocessing ", cong_str))
  fldr <- paste0(cong_str, "_", suffix)
  
  members <- readRDS(file.path(topdir, "Data",
                               paste0("inputs_", cong_str,
                                      ".RDS")))$members
  
  # read results files
  files <- file.path(topdir,"Results",  folder, fldr,
                     paste0(fldr,"_chain", chains, ".RDS"))
  
  out_chains <- map(files[file.exists(files)],  ~readRDS(.x))
  if (stage==0){
    for (i in 1:length(out_chains)){
      out_chains[[i]]$eps <- out_chains[[i]]$eta!=0
      out_chains[[i]]$model_size <- rowSums(out_chains[[i]]$eps)-1
    }
    
    pars2 <- c("eta","eps","theta")
    stage2_pars <- map(out_chains, ~(.x[pars2]))
    accept_pars <- names(out_chains[[1]])[grepl("accept",names(out_chains[[1]]))]
    accept_rates <- map(accept_pars,  
                        ~sapply(out_chains, function(x){x[[.x]]}))
  }
  if (stage==1 & !bare_bones){
    for (i in 1:length(out_chains)){
      out_chains[[i]]$theta <- out_chains[[i]]$theta[,1]
    }
  }
  
  # combine chains
  if (exclude & !is.null(keep_chains)){
    if (cong_num %in% keep_chains[[folder]]$skip){
      keep = list()
    } else if (cong_num %in% names(keep_chains[[folder]])){
      keep = keep_chains[[folder]][[as.character(cong_num)]]
    } else {
      keep = keep_chains[[folder]]$default
    }
  } else {
    keep=NULL
  }
  
  out <- combine_chains(out_chains, members$party, keep=keep)

  # results (ideal points)
  if (!bare_bones){
    res <- members %>% 
      select(legislator = Member, party) %>%
      distinct() %>%
      mutate(beta0=colMeans(out$beta0),
             beta1=colMeans(out$beta1),
             rank_beta0 = rank(beta0),
             rank_beta1 = rank(beta1),
             # prob of change
             p_change=1-colMeans(out$zeta),
             rank_pchange = rank(-1*p_change, ties="random"),
             # diff before and after
             beta_diff = 
               map(1:nrow(members), 
                   ~(out$beta1[out$zeta[,.x]==0,.x]- 
                       out$beta0[out$zeta[,.x]==0,.x])),
             diff_mean = map(beta_diff, mean),
             diff_L95 = map(beta_diff, quantile, .025),
             diff_U95 = map(beta_diff, quantile, .975),
      ) %>% 
      select(-beta_diff) %>%
      unnest(contains("diff"))  
  } else {
    res <- members %>% 
      select(legislator = Member, party) %>%
      distinct() %>%
      mutate(# prob of change
             p_change=1-colMeans(out$zeta))
  }
  
  print(round(map_dbl(out$zeta_chains, mean),2))
  print(round(map_dbl(out$zeta_chains, mean),2) *
    ifelse(1:length(out$zeta_chains) %in% keep, 1, NA))
  
  summy <- tibble(type=c("All","D","R"),
                  chains = list(out$zeta_chains,   
                                out$zetaD_chains,  
                                out$zetaR_chains),
                  draws = map(chains, ~unlist(.x[out$keep])),
                  ASF = map(draws,mean), 
                  ASF_L95 = map(draws, quantile, .025),
                  ASF_U95 = map(draws, quantile, .975)) %>%
    unnest(c(ASF, ASF_L95, ASF_U95)) %>%
    select(-chains, -draws)
  
  # res for individual chains
  res_pars <- c("beta0","beta1","zeta")
  res_chains <- list()
  for (par in res_pars){
    if (par %in% names(out_chains[[1]])){
       res_chains[[par]] <- map(out_chains, ~colMeans(.x[[par]])) 
    }
  }
  
  # save output
  output <- list(res=res, res_chains=res_chains,
                 out=out[grepl("chains",names(out)) | names(out)=="keep"],
                 summy=summy)
  if (stage==0){
    output$stage2 = stage2_pars[output$out$keep]
    output$accept_rates = tibble(par=accept_pars,
                                 rate=map(accept_rates, ~(.x))) %>%
      unnest(rate) %>% mutate(chain=rep(1:length(out_chains), length(accept_pars))) %>%
      filter(chain %in% output$out$keep)
    output$stage2$covariates = out_chains[[1]]$covariates
  }
  if (retain_zeta){
    output$zeta <- map(out_chains[output$out$keep], ~(.$zeta))
  }
  
  filename <- paste0(fldr, "_combined", ".RDS")
  dir.create(file.path(topdir,"Results", folder,"Combined"), showWarnings = F)
  saveRDS(output, file.path(topdir,"Results", folder,
                            "Combined",filename))
  
}

#####################

args = commandArgs(trailingOnly = T)
args <- as.list(args)
names(args) <- c("cong","stage","fldr","suffix")
for (x in names(args)){
  if (grepl("^[0-9]+$", args[[x]])){
    args[[x]] <- as.numeric(args[[x]])
  }    
}

postprocess(args$cong, args$fldr, 
            args$suffix,chains=1:8, stage=args$stage,
            exclude=T, bare_bones = (args$stage==1))
