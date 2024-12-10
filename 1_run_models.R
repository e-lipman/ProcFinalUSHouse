library(tidyverse)
library(yaml)

print_time <- function(str=""){
  print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), ": ", str))
}
configs <- read_yaml("configs.yml")

# source sampler function

print_time("Compiling C code")
Rcpp::sourceCpp("src/run_sampler.cpp", rebuild=F)

# get arguments

args = commandArgs(trailingOnly = T)

args <- as.list(args)
if (length(args)==7){ # full model
  names(args) <- c("cong", "chain", "stage", 
                   "burn", "iter", "thin", "folder")
  stopifnot(args$stage %in% c(0,1))
} else if (length(args)==11){
  names(args) <- c("cong", "chain", "stage", "cut",
                   "burn", "iter", "thin", "folder",
  		   "infolder", "insuffix","thresh")
  stopifnot(args$stage==2, args$cut==0)
} else if (length(args)==8){
  names(args) <- c("cong", "chain", "stage", "cut",
                   "steps", "folder", "infolder", "insuffix")
  stopifnot(args$stage==2, args$cut==1)
} else {
  stop("Invalid number of arguments")
}

for (x in names(args)){
  if (grepl("^[0-9]+$", args[[x]])){
    args[[x]] <- as.numeric(args[[x]])
  }    
}

set.seed(10*77606 + args$chain)

# other hyperparams
if (is.null(args$cut)){args$cut=F}
args$cut = as.logical(args$cut)
args$progress = ifelse(args$stage==2 & args$cut, 100, 1000)
print(paste(names(args), args, collapse = ", "))

# set up output folder
cong_str <- paste0("H",str_pad(args$cong,3,pad=0))

if (!(args$stage==2 & args$cut)){
	outname <- paste0(cong_str, "_",
	args$burn,  "_", args$iter, "_", args$thin)  
} else {
   outname <- paste0(cong_str, "_cut_", args$steps)  
}
  
dir.create(file.path("Results"), showWarnings = F)
dir.create(file.path("Results", args$folder), showWarnings = F)
dir.create(file.path("Results", args$folder,outname), showWarnings = F)

# get model inputs
inputs <- readRDS(file.path("Data",paste0("inputs_", cong_str, ".RDS")))

if (args$stage==2){
  res1 <- readRDS(file.path("Results",args$infolder, "Combined", 
                            paste0(cong_str, "_", args$insuffix, 
                                   "_combined.RDS")))
  stopifnot(nrow(res1$res)==nrow(inputs$covariates))
  stopifnot(inputs$members$Member==res1$legislator)
  
  if (!args$cut){
    zeta <- as.numeric(1-res1$res$p_change > args$thresh)
  } else {
    zeta <- res1$zeta[[args$chain]]
    stopifnot(ncol(zeta)==nrow(inputs$covariates))
  }
}

# prepare covariates for models with covaiate regression

if (args$stage %in% c(0,2)){ 
  standardize <- function(col){
    if (all(col %in% 0:1)){
      return(col-mean(col))
    } else {
      return( (col-mean(col))/sd(col) )
    }
  }
  
  x <- mutate_all(as.data.frame(inputs$covariates), standardize) %>%
    as.matrix() 
  
  # EXCLUDE PaRTY iNDICATOR
  x <- x[,configs$covariates]
  
  eps_init_method = ((args$chain-1) %% 2)
} else { 
  x <- NULL 
  eps_init_method = 0
}

# run sampler

if (args$stage==2){
  if (!args$cut){
    out <- run_covreg_only(zeta, x, 
                           burn = args$burn, iter = args$iter, 
                           thin = args$thin, progress = args$progress,
                           eps_init_method = eps_init_method)
  } else {
    out <- run_covreg_only_multi(zeta, x, 
                                 steps = args$steps, progress=args$progress, 
                                 eps_init_method = eps_init_method)
  }
} else {
  out <- run_sampler(inputs$y, inputs$missing,     
                     inputs$demLeader, inputs$repLeader,    
                     inputs$gam0_max,    
                     iter = args$iter, burnin = args$burn,    
                     thin = args$thin, progress=args$progress,
                     cut = args$cut,
                     x_ = x,
                     eps_init_method = eps_init_method)
}

if (args$stage %in% c(0,2)){
  out$covariates <- colnames(x)
  out$eta <- cbind(out$intercept, out$eta)
}

print_time("Saving output")
filename <- paste0(outname, "_chain", args$chain, ".RDS")
saveRDS(out, file.path("Results", args$folder, outname, filename))

print_time("Done!")
