library(tidyverse)
library(yaml)
library(ggcorrplot)

source("helper_funs.R")

configs <- read_yaml(file.path("IdealPointsCompare",
                               "configs.yml"))

in_path <- file.path("IdealPointsCompare")
out_path <- file.path("Figures")
dir.create(out_path, showWarnings = F)
dir.create(file.path(out_path,"Data"), showWarnings = F)

# helper functions
make_shading <- function(dat_summary){
  dat_summary %>% 
    select(cong, party=party_control) %>%
    mutate(type = 
             case_when(party!=lag(party) | is.na(lag(party)) ~ "start", 
                       party!=lead(party) | is.na(lead(party)) ~ "end")) %>%
    filter(!is.na(type)) %>%
    group_by(type) %>%
    mutate(run=row_number()) %>%
    ungroup() %>%
    pivot_wider(names_from=type, values_from=cong) %>%
    mutate(end=ifelse(is.na(end), start, end)) %>%
    filter(party=="R") 
}

summarise_one <- function(cong_num){
  dat <- readRDS(file.path(in_path, "Data",
                           paste0("inputs_H",
                                  str_pad(cong_num, 3, pad=0),
                                  ".RDS")))
  
  # covariate summaries
  X <- dat$covariates[,configs$covariates]
  X <- mutate_all(as.data.frame(X), standardize) %>% as.matrix() 
  
  PC_full <- princomp(X)
  PC <- tibble(num=1:ncol(X),
               percent=cumsum(PC_full$sdev^2 / sum(PC_full$sdev^2)))
  
  prv <- dat$covariates$p_R +.5
  maj_party <- dat$cong_info$partyControl
  min_party <- setdiff(c("D","R"), maj_party)
  party <- ifelse(dat$covariates$belto.partyControl, maj_party, min_party)
                  
  
  tibble(cong = cong_num,
         members = nrow(dat$y),
         votes_proc = dat$gam0_maxC+1,
         votes_final = ncol(dat$y) - votes_proc,
         cond_num = kappa(X, exact = T),
         PC = list(PC),
         cor = list(cor(X)),
         prv=list(tibble(party=party, prv=prv)),
         with_party = list(percent_with_party(dat)),
         party_control = maj_party)
}

# load data
dat_summary <- map_df(93:113, summarise_one)

shading <- make_shading(dat_summary)

############################################

# FIGURES FOR PAPER

# Figure 1: correlation matrices
cov_cat <- tibble(type=names(configs$covariate_cats),
                  cov=configs$covariate_cats) %>%
  unnest(cov)
breaks <- cumsum(count(cov_cat,type)$n[3:2])+0.5

for (CONG in c(93,98,103,108)){
  cor <- filter(dat_summary, cong==CONG)[["cor"]][[1]]
  cor <- cor[cov_cat$cov, cov_cat$cov]
  g <- ggcorrplot(cor) +
    geom_hline(yintercept=breaks, alpha=.5) +
    geom_vline(xintercept=breaks, alpha=.5) +
    ggtitle(paste0("H", str_pad(CONG, 3, pad="0"))) +
    xlab("") + ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  print(g)
  
  # FIGURE 1
  ggsave(file.path(out_path, "Data",
                   paste0("H", str_pad(CONG, 3, pad="0"), "_cor.jpeg")),
         height = 6, width=6.5)
}

# FIGURE 2: percent voting with party
dat_with_party <- dat_summary %>%
  select(cong, with_party) %>%
  unnest(with_party) %>%
  pivot_longer(ends_with("with_party"),
               names_to="type", values_to="percent") %>%
  group_by(cong, type, party) %>%
  summarise(percent=median(percent)) %>%
  mutate(type=ifelse(grepl("final",type), "Final","Procedural"))

dat_with_party %>%
  ggplot(aes(x=cong, y=percent, color=party, linetype=type)) +
  geom_rect(data=shading, inherit.aes = F, alpha=.2,
            aes(ymin=min(dat_with_party$percent), 
                ymax=1, 
                xmin=start-.5, xmax=end+.5)) +
  geom_line() +
  xlab("House") +
  scale_y_continuous("Percent voting with party") +
  scale_color_manual("Party", values=c("blue","red")) +
  scale_linetype_discrete("Vote type") +
  theme_bw()

ggsave(file.path(out_path, "Data","data_with_party.jpeg"),
       width = 6, height=2.5)

# FIGURE 3: republican voteshare
dat_summary %>%
  select(cong, prv) %>%
  unnest(prv) %>%
  ggplot(aes(x=as.factor(cong), y=prv, color=party)) +
  geom_hline(yintercept=.5, linetype="solid", alpha=.5) +
  geom_boxplot(outlier.size = .5) +
  xlab("House") +
  scale_y_continuous("Republican voteshare", limits = c(0,1)) +
  scale_color_manual("Party", values=c("blue","red")) +
  theme_bw()
ggsave(file.path(out_path, "Data", "data_prv.jpeg"),
       width=6, height=2.5)
