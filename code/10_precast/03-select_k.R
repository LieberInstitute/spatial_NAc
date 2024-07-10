library(here)
library(PRECAST)
library(HDF5Array)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(sessioninfo)
library(tidyverse)
library(Matrix)
library(SpatialExperiment)
library(ggsci)
library(ggpubr)
library(getopt)
library(scran)
library(scater)
library(spatialNAcUtils)
library(spatialLIBD)

degree_freedom <- function(K, paraList){
  
  p  <- paraList$p; q <- paraList$q
  r_max <- paraList$r_max; Sigma_equal <- paraList$Sigma_equal
  Sigma_diag <- paraList$Sigma_diag;  mix_prop_heter <- paraList$mix_prop_heter
  if(!Sigma_equal){
    if(mix_prop_heter){
      
      if(Sigma_diag){ # our model
        
        # # beta_r + W + Lam_r + Mu_k + Sigma_k + Psi_r
        dfree <- 1*r_max+p*(q+r_max) + K*(q+ q) + q*(q+1)/2*r_max
      }else{
        dfree <- 1*r_max+p*(q+r_max) + K*(q + q*(q+1)/2.0)+  q*(q+1)/2*r_max
      }
    }else{
      
      if(Sigma_diag){
        
        # # beta + W + Lam_r + Mu+Sigma + tau_r
        dfree <- 1+p*(q+r_max) + K*(q+q) +  q*(q+1)/2*r_max
      }else{
        dfree <- 1+p*(q+r_max) + K*(q+q*(q+1)/2.0)+  q*(q+1)/2*r_max
      }
    }
  }else{
    if(mix_prop_heter){
      if(Sigma_diag){ ## This is our model setting
        message("Sigma is set to diagonal matrices \n")
        # # beta_r + W + Lam_r + Mu+ Sigma_r + Phi_r
        dfree <- 1*r_max+p*(q+r_max) + K*q+ q + r_max *q*(q+1)/2
      }else{
        message("Sigma is set to dense matrices \n")
        dfree <- 1*r_max+p*(q+r_max) + K*q+ q*(q+1)/2.0+ r_max *q*(q+1)/2
      }
    }else{
      if(Sigma_diag){
        message("Sigma is set to diagonal matrices \n")
        # # beta + W + Lam_r + Mu_k +Sigma_r + Phi_r
        dfree <- 1+p*(q+r_max) + K*q+ q  + r_max *q*(q+1)/2
      }else{
        message("Sigma is set to dense matrices \n")
        dfree <- 1+p*(q+r_max) + K*q+ q*(q+1)/2.0 + r_max *q*(q+1)/2
      }
    }
  }
  return(dfree)
}

resDir <- here("processed-data", "10_precast", "nnSVG_precast")
nRandom_starts <- 3
K <- seq(3, 15, 1)
precast_results <- lapply(c(1:nRandom_starts), function(rs){
  cat(rs, "\n")
  lapply(K, function(k){
    readRDS(paste0(resDir, "/random_start_", rs, "/pre_obj_k", k, ".rds"))
  })
})

para_settings <- lapply(precast_results, function(pre_list){
  lapply(pre_list, function(pre_obj){
    attr(pre_obj@resList, 'para_settings')
  })
})

pen_const <- 50
final_ll_df <- list()
for(i in c(1:nRandom_starts)){
  ll <- lapply(precast_results[[i]], function(pre_obj){
    pre_obj@resList[[1]]$loglik
  })
  K <- lapply(para_settings[[i]], function(para_setting){
    para_setting$K
  })
  n <- lapply(para_settings[[i]], function(para_setting){
    para_setting$n
  })
  p <- lapply(para_settings[[i]], function(para_setting){
    para_setting$p
  })
  ll <- unlist(ll)
  K <- unlist(K)
  n <- unlist(n)
  p <- unlist(p)
  nK <- length(K)
  MAIC_vec <-  rep(Inf, nK)
  AIC_vec <-  rep(Inf, nK)
  MBIC_vec <-  rep(Inf, nK)
  BIC_vec <-  rep(Inf, nK)
 
  for(j in c(1:nK)){
     dfree <- degree_freedom(K[j], para_settings[[i]][[j]])
     MAIC_vec[j] = -2.0*ll[j] + dfree*2*log(log(p[j]+n[j]))*pen_const
     AIC_vec[j] <- -2.0*ll[j] + dfree*2
     MBIC_vec[j] <-  -2.0*ll[j] + dfree*log(n[j])*log(log(p[j]+n[j]))*pen_const
     BIC_vec[j] <-  -2.0*ll[j] + dfree*log(n[j])*log(log(p[j]+n[j]))*pen_const
  }

  final_ll_df[[i]] <-  data.frame("K" = K, "random_start" = i, "LL_model" = ll, "n" = n, "p" = p,
                                  "MAIC" = MAIC_vec, "AIC" = AIC_vec, "MBIC" = MBIC_vec, "BIC" = BIC_vec)
}


ll_df <- do.call(rbind, final_ll_df)
ll_df$random_start <- factor(ll_df$random_start, levels = c(1:nRandom_starts))

plotDir <- here("plots", "10_precast", "nnSVG_precast", "BIC_select")
pdf(file.path(plotDir, "loglikelihood_vs_K.pdf"), height = 4, width = 6)
ggplot(ll_df, aes(x = K, y = BIC, color = random_start, fill = random_start)) + geom_point() + geom_line() + theme_classic()
dev.off()



if(!inherits(res_obj, 'SeqK_PRECAST_Object')) stop('res_obj must be "SeqK_PRECAST_Object"!\n')
para_settings <- attr(res_obj, 'para_settings')
K_set <- para_settings$K
nK <- length(K_set)
MAIC_vec <-  rep(Inf, nK)
AIC_vec <-  rep(Inf, nK)
MBIC_vec <-  rep(Inf, nK)
BIC_vec <-  rep(Inf, nK)
if(length(res_obj) != nK) stop("The length of obj must be equal to length of K!")
n <- para_settings$n; p <- para_settings$p


ll_final <- lapply(res_obj, function(resList){
  resList$loglik
})
ll_final <- unlist(ll_final)
final_ll_df <- data.frame("K" = K_set, "LL_model" = ll_final)
pdf(file.path(plotDir, "loglikelihood_vs_K.pdf"), height = 4, width = 6)
ggplot(final_ll_df, aes(x = K, y = LL_model)) + geom_point() + geom_line() + theme_classic()
dev.off()

for(i in c(1:nRandom_starts)){
 nK <- dim(final_ll_df[[i]])[1]
 MAIC_vec <- rep(NA, nK)
 AIC_vec <- rep(NA, nK)
 BIC_vec <- rep(NA, nK)
 pre_list <- precast_results[[i]]
 for(j in c(1:nK)){
   resList <- pre_list[[j]]@resList[[1]]
 }
}

for(k in 1:nK){
    cat(k, "\n")
    resList <- res_obj[[k]]
    dfree <- degree_freedom(K_set[k], para_settings)
    MAIC_vec[k] = -2.0*resList$loglik + dfree*2*log(log(p+n))*pen_const
    AIC_vec[k] <- -2.0*resList$loglik + dfree*2
    MBIC_vec[k] <-  -2.0*resList$loglik + dfree*log(n)*log(log(p+n))*pen_const
    BIC_vec[k] <-  -2.0*resList$loglik + dfree*log(n)*log(log(p+n))*pen_const
}


summary_stats <- data.frame("K" = K_set, "MAIC" = MAIC_vec, "AIC" = AIC_vec,
                            "MBIC" = MBIC_vec, "BIC" = BIC_vec)

pdf(file.path(plotDir, "MAIC_vs_K.pdf"))
ggplot(summary_stats, aes(x = K, y = MAIC)) + geom_point() + geom_line(linetype = "dashed") + theme_classic() +
xlab("Cluster resolution (K)") +  ylab("MAIC")
dev.off()

pdf(file.path(plotDir, "AIC_vs_K.pdf"))
ggplot(summary_stats, aes(x = K, y = AIC)) + geom_point() + geom_line(linetype = "dashed") + theme_classic() +
xlab("Cluster resolution (K)") +  ylab("AIC")
dev.off()

pdf(file.path(plotDir, "BIC_vs_K.pdf"))
ggplot(summary_stats, aes(x = K, y = BIC)) + geom_point() + geom_line(linetype = "dashed") + theme_classic() +
xlab("Cluster resolution (K)") +  ylab("BIC")
dev.off()

pdf(file.path(plotDir, "MBIC_vs_K.pdf"))
ggplot(summary_stats, aes(x = K, y = MBIC)) + geom_point() + geom_line(linetype = "dashed") + theme_classic() +
xlab("Cluster resolution (K)") +  ylab("MBIC")
dev.off()
