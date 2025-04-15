library(tidyverse)
library(yaml)
library(scales)
library(pander)

source("helper_funs.R")

configs <- read_yaml("configs.yml")
args = commandArgs(trailingOnly = T)
if (length(args)==0){
  args <- list(test=T)  
} else {
  args <- as.list(args)
  names(args) <- c("test")  
}
test=(args$test==1)

FLDR = "joint"
if (test){
  CONGS=93:94
} else {
  CONGS=93:113
}

in_path <- file.path(".")
out_path <- file.path("Figures")
dir.create(out_path, showWarnings = F)

###########################################

res_all <- 
  tibble(cong=CONGS) %>%
  mutate(res=map(cong, read_one, folder=FLDR)) %>%
  unnest(c(res))

summy_all <- select(res_all, cong, summy) %>% unnest(summy)

get_model_one <- function(res2, cov_scale){
  covariates <- res2$covariates
  res2 <- res2[names(res2)!="covariates"]
  eps <- map_dfr(res2, ~tibble(.$eps))[[1]]
  eta <- map_dfr(res2, ~tibble(.$eta))[[1]]
  if (ncol(eps)==length(covariates)){
    eps <- cbind(1,eps)
  }
  
  tibble(cov=c("int",covariates),
         eps=colMeans(eps),
         scale = c(1, cov_scale),
         eta_L95=apply(eta,2,quantile,.025),
         eta_U95=apply(eta,2,quantile,.975),
         eta=apply(eta,2,quantile,.5)) %>%
    filter(cov!="int")
}

get_model_size <- function(res2){
  size <- map_dfr(res2[names(res2)!="covariates"], ~tibble(.$eps)) %>%
    rowSums()
  tibble(size=size-1) %>%
    summarise(p=median(size), 
              L95 = quantile(size,.025),
              U95 = quantile(size, .975))
}

model_all <- select(res_all, cong,
                    stage2, cov_scale) %>%
  mutate(model = map2(stage2, cov_scale, get_model_one)) %>%
  select(-stage2, -cov_scale) %>%
  unnest(model) %>%
  group_by(cov) %>%
  mutate(sign = case_when(eps<.5~"",
                          eta<0~"-", eta>0~"+", eta==0~"o"),
         n_selected = sum(eps>.5),
         avg_post = mean(eps),
         max_post = max(eps)) %>%
  ungroup()

size_all <- mutate(res_all,
                   size=map(stage2, get_model_size)) %>%
  select(cong, size) %>%
  unnest(size)

party_info <- readRDS(file.path(in_path,"Data","party_info.RDS"))


###########################################
#                 ASF                     #
###########################################

# FIGURE 4
summy_all %>% 
  filter(type=="All") %>%
  ggplot(aes(x=cong, y=ASF, ymin=ASF_L95, ymax=ASF_U95)) +
  geom_linerange(position = position_dodge(.5)) +
  geom_point(position = position_dodge(.5), size=1) +
  geom_rect(data=party_info, ymin=0, ymax=1,
            aes(xmin=start-.5, xmax=end+.5),
            inherit.aes = F, alpha=.2) +
  xlab("House") + ylab("Bridging frequency") +
  theme_bw()

ggsave(file.path(out_path,"ASF.jpeg"), height = 2.5, width=6)

# FIGURE 8
summy_all %>% 
  ggplot(aes(x=cong, y=ASF,  ymin=ASF_L95,  ymax=ASF_U95, 
             color=type, linetype=type)) +
  geom_linerange(position = position_dodge(.5)) +
  geom_point(position = position_dodge(.5), size=1) +
  geom_rect(data=party_info, ymin=0, ymax=1,
            aes(xmin=start-.5, xmax=end+.5),
            inherit.aes = F, alpha=.2) +
  xlab("House") + ylab("Bridging frequency") +
  scale_color_manual("Party",values=c("black","blue","red")) +
  scale_linetype_manual("Party",values=c("solid","dashed","dashed")) +
  theme_bw()

ggsave(file.path(out_path,"ASF_party.jpeg"), height = 2.5, width=6)

###########################################
#              Model size                 #
###########################################

# FIGURE 5
size_all %>%
  ggplot(aes(x=cong, y=p, ymax=L95, ymin=U95)) +
  geom_rect(data=party_info, ymin=0, ymax=21,
            aes(xmin=start-.5, xmax=end+.5),
            inherit.aes = F, alpha=.2) +
  geom_linerange(position=position_dodge(.5)) +
  geom_point(position=position_dodge(.5), size=1) +
  xlab("Congress") +
  ylab("Median model size") +
  theme_bw()

model_all %>%
  group_by(cong) %>%
  summarise(p=sum(eps>.5)) %>% 
  ggplot(aes(x=cong, y=p)) +
  geom_rect(data=party_info, ymin=0, ymax=21,
            aes(xmin=start-.5, xmax=end+.5),
            inherit.aes = F, alpha=.2) +
  #geom_linerange(position=position_dodge(.5)) +
  geom_point(position=position_dodge(.5), size=1) +
  xlab("House") +
  ylab("Median model size") +
  theme_bw()

ggsave(file.path(out_path,"size.jpeg"),
       height = 2.5, width=6)

###########################################
#            SELECTED MODELS              #
###########################################

# FIGURE 6

cov_cat <- tibble(type=names(configs$covariate_cats),
                  cov=configs$covariate_cats) %>%
  unnest(cov)

dat_colorplot <- model_all %>%
  select(cong,cov,eps,n_selected,sign) %>%
  full_join(cov_cat, by="cov") %>%
  arrange(type,n_selected) %>%
  mutate(cov=fct_relevel(cov,unique(cov)))

# save cov order for elsewhere
saveRDS(
  levels(dat_colorplot$cov)[n_distinct(dat_colorplot$cov):1],
  "cov_order.RDS"
)

dat_colorplot %>% 
  rename(PIP=eps) %>%
  ggplot(aes(x=cong, y=cov, fill=PIP)) +
  geom_tile() +
  #geom_point(data=filter(dat_colorplot,eps>.5)) +
  geom_text(aes(label=sign), size=5) +
  geom_hline(yintercept=nrow(cov_cat)-cumsum(count(cov_cat,type)$n[3:1])+0.5,
             color="white") +
  scale_fill_gradientn(colors=c("blue","white","red"), 
                       breaks=c(0,.5,1)) +
  xlab("House") + ylab("") +
  theme_minimal()
ggsave(file.path(out_path,"red_blue.jpeg"), height = 4, width=6)

###########################################
#              COEFFICIENTS               #
###########################################

# FIGURE 7

plot_cov <- function(which_cov,title="", c=1,
                     lab_pow=10){
  
  plot_dat <- model_all %>%
    filter(cov==which_cov) %>%
    mutate_at(vars(contains("eta")), ~exp(.x*(c/scale)))
  
  print(paste0("Showing ORs for an increase of (median) ", 
               round(c/median(plot_dat$scale),2),
               " standard deviations"))
  
  shading <- 
    plot_dat %>% 
    #group_by(type) %>%
    summarize(min=min(c(eta,eta_L95)), max=max(eta,eta_U95)) %>%
    mutate(party_info = list(party_info)) %>%
    unnest(party_info)
  
  lrange = round(log(c(min(shading$min), max(shading$max)),
                     lab_pow))
  lrange =seq(lrange[1],lrange[2])
  labels = ifelse(lrange>=0,
                  round(lab_pow^lrange),lab_pow^lrange)
  
  g <- plot_dat %>% 
    ggplot(aes(x=eta, y=cong, xmin=eta_L95,xmax=eta_U95)) +
    geom_rect(data=shading, inherit.aes = F, alpha=.2,
              aes(xmin=min, xmax=max, ymin=start-.5, ymax=end+.5)) +
    #geom_vline(data=refline, inherit.aes=F,
    #           aes(xintercept=ref), linetype="dotted") +
    geom_vline(xintercept=1, linetype="dotted") +
    geom_pointrange() +
    scale_x_continuous(trans = log_trans(),
                       labels = labels,
                       breaks = lab_pow^lrange) +
    #scale_x_log10() +
    #facet_wrap(type~., ncol=1, scales="free_y") +
    coord_flip() +
    xlab("Increase in odds ratio\n(log scale)") + ylab("House") +
    ggtitle(title) +
    theme_bw()
  plot(g)
  
  num_summy <- select(plot_dat, cong, eps, OR=eta) %>%
    filter(eps>.5)%>%
    mutate(eta=log(OR), OR_inv=1/OR) %>%
    arrange(OR)
  print(num_summy)
  
  ggsave(file.path(out_path, paste0(which_cov,".jpeg")), 
         height = 2.5, width=6)
}

select(model_all, cov, n_selected) %>%
  distinct() %>%
  arrange(-n_selected) %>%
  filter(n_selected>3) 

#paste0("Belongs to majority party",
#        "\nOR for majority vs. minority")
plot_cov("belto.partyControl", c=1)

#paste0("Constituancy (Democratic legislators only)",
#"\nOR for increase of 5 percent")
plot_cov("p_R_D", c=.05, lab_pow = 2)

#paste0("Constituancy (Republican legislators only)",
#"\nOR for increase of 5 percent")
plot_cov("p_R_R", c=.05, lab_pow = 2)

