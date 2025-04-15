library(tidyverse)
library(yaml)
library(scales)

source("helper_funs.R")
source("revisions_letter/auxiliary_regression_funcs.R")

configs <- read_yaml(file.path("IdealPointsCompare","configs.yml"))
in_path <- file.path("..","SenateIdeals")
out_path <- file.path("Figures", "Revision")
dir.create(out_path, showWarnings = F)

CONGS <- 93:113
FLDRS <- c("stage1","joint")
test <- F

res_all <- 
  cross_df(list(cong=CONGS, folder=FLDRS)) %>%
  mutate(res=map2(cong, folder, read_one),
         label=case_when(folder=='stage1'~'No Covariates',
                         folder=='joint'~'Covariates')) %>%
  unnest(c(res)) %>%
  dplyr::select(cong, folder, label, summy, res, stage2)

#-------------------------------------------------------------------------------

# FIGURE 1: Compare ranks of ideal points to those from independent models

plot_ideal_points <- function(which_cong, separate=T){
  cong_str <- paste0("H", str_pad(which_cong, 3, pad="0"))
  
  if (separate){
    in_path_separate <- file.path("auxiliary_res",
                                  paste0(cong_str, "_summy_mcmcpack.RDS"))
    res <- readRDS(in_path_separate)  
  } else {
    res <- res_all %>% 
      filter(cong==which_cong, label=="Covariates") %>%
      dplyr::select(cong, folder, label, res) %>%
      pull(res) %>% .[[1]]  
  }
  
  # plot
  res %>% 
    ggplot(aes(x=rank_beta0, y=rank_beta1, 
               color=party, shape=party)) +
    geom_point(show.legend=F) +
    scale_color_manual(values=c("blue","red")) +
    xlab("Rank for procedural votes") +
    ylab("Rank for final passage votes") +
    coord_equal() +
    theme_bw()
  
  save_filename <- paste0("revision_fig1_", cong_str, 
                          ifelse(separate, "_sep", "_joint"),
                          ".jpeg")
  ggsave(file.path(out_path, save_filename), 
         height = 3.5, width=4)
}

# 95 has min BF, 111 has max BF
#run_proc_final(111, load_saved = F, in_path = in_path)
plot_ideal_points(111, separate = T)
plot_ideal_points(111, separate = F)

#-------------------------------------------------------------------------------

# FIGURE 2: 'Simper approach' proposed by reviewer 2

# response: Difference in log odd between percent voting with party on proc and final

if (F){
  #res_with_party_penalized <- tibble(cong=CONGS) %>%
  #  mutate(res=map(cong, run_regression_percent_with_party, 
  #                 in_path=in_path,
  #                 penalized=T, bayes=F)) %>%
  #  unnest(res)
  #saveRDS(res_with_party_penalized, file.path("auxiliary_res","withparty_penalized.RDS"))  
  
  res_with_party_bayes <- tibble(cong=CONGS) %>%
    mutate(res=map(cong, run_regression_percent_with_party, 
                   in_path=in_path,
                   penalized=F, bayes=T)) %>%
    unnest(res)
  saveRDS(res_with_party_bayes, file.path("auxiliary_res","withparty_bayes.RDS"))  
} else {
  #res_with_party_penalized <- readRDS(file.path("auxiliary_res","withparty_penalized.RDS"))  
  res_with_party_bayes <- readRDS(file.path("auxiliary_res","withparty_bayes.RDS"))  
}

ORDER <- readRDS("cov_order.RDS")

#hacky - flip signs of coefs so outcome is log-oods(proc)-log-odds(final)
res_with_party_bayes <-
  res_with_party_bayes %>% mutate(coef=-1*coef)

plot_redblue_one(res_with_party_bayes, ORDER)
ggsave(file.path(out_path, "revision_fig2_alt_regression_bayes.jpeg"),
       height=3.5, width=6)

#-------------------------------------------------------------------------------

# FIGURE 3: Compare bridging frequencies
party_info <- readRDS(file.path("IdealPointsCompare","Data","party_info.RDS"))

summy <- res_all %>% 
  dplyr::select(cong, folder, label, summy) %>%
  unnest(summy)

summy %>% 
  filter(type=="All" | is.na(type)) %>%
  ggplot(aes(x=cong, y=ASF, ymin=ASF_L95, ymax=ASF_U95, linetype=label)) +
  geom_linerange(position = position_dodge(.5),
                 show.legend = F) +
  geom_point(position = position_dodge(.5), size=1) +
  geom_rect(data=party_info, ymin=0, ymax=1,
            aes(xmin=start-.5, xmax=end+.5),
            inherit.aes = F, alpha=.2) +
  xlab("House") + ylab("Bridging frequency") +
  scale_linetype_discrete("") +
  theme_bw()

ggsave(file.path(out_path, "revision_fig3_compare_ASF.jpeg"), 
       height = 2.5, width=6)

#-------------------------------------------------------------------------------

# FIGURE 4: Compare bridges

bridge_probs <- res_all %>% 
  dplyr::select(cong, folder, label, res) %>%
  unnest(res) %>%
  mutate(bridge_prob=1-p_change) %>%
  dplyr::select(cong, folder,legislator,  bridge_prob) %>%
  pivot_wider(names_from="folder",values_from="bridge_prob")


compare_bridges <- bridge_probs %>% 
  count(cong, stage1=stage1>=.5, joint=joint>=.5) %>% 
  group_by(cong) %>% mutate(p=n/sum(n)) %>% ungroup() %>% 
  mutate(label=case_when(stage1 & joint ~ "Total Agreement", 
                         !stage1 & !joint ~ "Total Agreement",
                         !stage1 & joint ~ "Covariates",
                         stage1 & !joint ~ "No Covariates")) %>%
  group_by(cong, label) %>%
  summarise(p=sum(p)) %>%
  ungroup()

compare_bridges %>%
  mutate(label=fct_relevel(label, unique(compare_bridges$label)[c(3,1,2)])) %>%
  ggplot(aes(x=cong, y=p, 
             color=label, linetype=label)) +
  geom_rect(data=party_info, ymin=-.05, ymax=1.05,
            aes(xmin=start-.5, xmax=end+.5),
            inherit.aes = F, alpha=.1) +
  geom_line(show.legend = F) +
  xlab("House") + ylab("Percent of legislators") +
  scale_y_continuous(limits=c(0,1)) +
  theme_bw()

ggsave(file.path(out_path, "revision_fig4_compare_bridges.jpeg"), 
       height = 2.5, width=6)

#-------------------------------------------------------------------------------

## FIGURE 5: Bridge probs vs covariate adjusted priors

get_linpred_probs <- function(which_cong){
  cong_str = paste0("H",str_pad(which_cong, 3, pad="0"))
  
  # get X
  x_raw <- readRDS(file.path(in_path, "Data",
                             paste0("inputs_", cong_str,
                                    ".RDS")))$covariates
  x <- cbind(1,x_raw[,configs$covariates]) %>%
    mutate_all(standardize) %>%
    as.matrix()
  
  # get eta
  stage2 <- res_all %>% 
    filter(cong==which_cong, label=="Covariates") %>%
    dplyr::select(cong, folder, label, stage2) %>%
    pull(stage2) %>% .[[1]]   
  eta <- map(stage2[names(stage2)==""], ~(.x$eta)) %>%
    reduce(rbind)
  
  # get probs
  linpred <- x%*%t(eta)
  #linpred[abs(linpred)>100] <- sign(linpred)*100 # truncate for computation
  prob <- 1/(1+exp(-1*linpred))
  
  rowMeans(prob)
}

which_cong=107
prior_probs <- get_linpred_probs(which_cong)

plot_dat <- res_all %>% 
  filter(cong==which_cong, label=="Covariates") %>%
  dplyr::select(cong, folder, label, res) %>%
  pull(res) %>% .[[1]] %>%
  mutate(bridge_prob=1-p_change,
         bridge_prior=prior_probs)

plot_dat %>% 
  ggplot(aes(x=bridge_prior, y=bridge_prob, 
             color=party, shape=party)) +
  geom_point(show.legend=F) +
  scale_color_manual(values=c("blue","red")) +
  scale_x_continuous("Prior bridging probability\n(covariate adjusted)",
                     limits=0:1) +
  scale_y_continuous("Posterior bridging probability",
                     limits=0:1) +
  coord_equal() +
  theme_bw()
ggsave(file.path(out_path, "revision_fig5_bridgeprob_vs_prior.jpeg"),
       height=3.5, width=4)

#-------------------------------------------------------------------------------

# NOT USED: Two-stage approach but with LASSO/MCP

# MCP: minimax concave penalty (MCP) (Zhang, 2010).
# package ncvreg (Breheny and Huang, 2011)
# make both graphs

# WARNING: MCP did not converge for cong 95
if (F){
  set.seed(75108)
  res_stage2_penalized <- tibble(cong=CONGS) %>%
    mutate(res=map(cong, run_stage2_penalized, in_path=in_path)) %>%
    unnest(res)
  saveRDS(res_stage2_penalized,
          file.path("auxiliary_res","stage2_penalized.RDS"))  
} else {
  set.seed(75108)
  res_stage2_penalized <-
    readRDS(file.path("auxiliary_res","stage2_penalized.RDS"))
}

ORDER <- readRDS("cov_order.RDS")

plot_redblue_one(res_stage2_penalized, ORDER)
ggsave(file.path(out_path, "revision_fig4_penalized_reg.jpeg"), 
       height = 3.5, width=7)

