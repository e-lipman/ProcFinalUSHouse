# read results for one congress and one specific model
read_one <- function(cong_num, folder,
                     suffix="10000_5000_10"){
  
  suffix <- case_when(test~"10_10_1",
                      folder!='joint'~suffix,
                      cong_num==105~"20000_20000_25", 
                      T~"20000_15000_20")
  
  cong_str = paste0("H", str_pad(cong_num, 3, pad=0))
  
  # get res stage 1
  path1 = file.path(in_path,"Results",folder,"Combined",
                    paste0(cong_str, "_", suffix, "_combined.RDS"))
  res <- readRDS(path1)
  
  # data quantities (only for models with covariates)
  dat <- readRDS(file.path(in_path,"Data",
                           paste0("inputs_", cong_str, ".RDS")))
  
  if ("stage2" %in% names(res)){
    tht <- reduce(
      map(res$stage2[names(res$stage2)!="covariates"], ~(.x$theta)), rbind)
    covs <- tibble(p_R=dat$covariates$p_R,
                   recentArrivalPrcnt=dat$covariates$recentArrivalPrcnt,
                   dwnom1=dat$covariates$dwnom1)
    
    # for scaling coefficients
    binary <- map_dbl(dat$covariates[,configs$covariates],
                      ~all(.x %in% 0:1))
    cov_sd <- map_dbl(dat$covariates[,configs$covariates], sd)
    scale <- ifelse(binary, 1, cov_sd)  
  }
  
  if (!"stage2" %in% names(res)){
    tibble(summy=list(res$summy),
           res=list(res$res))
  } else {
    tibble(summy=list(res$summy), 
           stage2=list(res$stage2),
           #polar=list(polar),
           covs = list(covs),
           party = dat$cong_info$partyControl,
           cov_scale=list(scale),
           res=list(res$res))  
  }
}

# Calculate indicator from data for whether legislator votes with party
percent_with_party <- function(dat){
  y <- dat$y  
  party <- dat$members$party
  stopifnot(party %in% c("R","D"))
  gamma <- c(rep(0, dat$gam0_maxC+1), rep(1,ncol(y)-dat$gam0_maxC-1))
  R_vote <- as.numeric(colMeans(y[party=="R",], na.rm=T)>.5)
  D_vote <- as.numeric(colMeans(y[party=="D",], na.rm=T)>.5)
  
  R_with <- y==matrix(rep(R_vote, nrow(y)), nrow=nrow(y), byrow = T)
  D_with <- y==matrix(rep(D_vote, nrow(y)), nrow=nrow(y), byrow = T)
  
  R_with_type <- map(0:1, ~rowMeans(R_with[,gamma==.x], na.rm=T))
  D_with_type <- map(0:1, ~rowMeans(D_with[,gamma==.x], na.rm=T))
  
  tibble(mamber=dat$members$Member, 
         party=party, 
         proc_with_party=ifelse(party=="R",R_with_type[[1]],D_with_type[[1]]),
         final_with_party=ifelse(party=="R",R_with_type[[2]],D_with_type[[2]])) 
}

# standardize covariates
standardize <- function(col){
  if (all(col %in% 0:1)){
    return(col-mean(col))
  } else {
    return( (col-mean(col))/sd(col) )
  }
}

# plot heatmap of inclusion probs
# make red-blue plot
plot_redblue_one <- function(res, cov_order){
  
  cov_cat <- tibble(type=names(configs$covariate_cats),
                    cov=configs$covariate_cats) %>%
    unnest(cov) %>%
    mutate(cov=fct_relevel(cov,cov_order),
           type=fct_relevel(type, names(configs$covariate_cats))) %>%
    arrange(type, cov)
  
  dat_colorplot <- res %>%
    dplyr::select(cong,cov,incl_prob,coef,label) %>%
    full_join(cov_cat, by="cov") %>%
    mutate(cov=fct_rev(fct_relevel(cov,as.character(cov_cat$cov))),
           sign=case_when(incl_prob<.5~"",
                          sign(coef)==1~"+",
                          sign(coef)==(-1)~"-",
                          T~"o"))
  
  if (all(res$incl_prob %in% 0:1)){
    dat_colorplot %>% 
      mutate(include=ifelse(incl_prob==1,
                            "Included","Excluded"),
             include=fct_rev(include)) %>%
      ggplot(aes(x=cong, y=cov, fill=include)) +
      geom_tile() +
      geom_hline(yintercept=cumsum(count(cov_cat,type)$n[3:1])+0.5,
                 color="white") +
      coord_fixed() +
      facet_grid(.~label) +
      scale_fill_manual("", values=c("red","blue")) +
      xlab("House") + ylab("") +
      theme_minimal()   
  } else {
    dat_colorplot %>%
      rename(PIP=incl_prob) %>%
      ggplot(aes(x=cong, y=cov, fill=PIP)) +
      geom_tile() +
      geom_text(aes(label=sign), size=5) +
      geom_hline(yintercept=cumsum(count(cov_cat,type)$n[3:1])+0.5,
                 color="white") +
      scale_fill_gradientn(colors=c("blue","white","red"), 
                           breaks=c(0,.5,1)) +
      coord_fixed() +
      xlab("House") + ylab("") +
      theme_minimal() 
  }
  
}
