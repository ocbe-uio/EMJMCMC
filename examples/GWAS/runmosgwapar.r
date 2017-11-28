mclapply(FUN= function(n){setwd("/home/michaelh/__simulation"); system(paste0("./analysis_launch.sh ",n," ",n+1),intern=TRUE,wait=TRUE)} ,X = 1:100,mc.preschedule = F,mc.cores = 31)

