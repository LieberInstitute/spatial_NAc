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

resDir <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast")
nRandom_starts <- 5
K <- seq(3, 15, 1)
precast_results <- lapply(c(1:nRandom_starts), function(rs){
  cat(rs, "\n")
  lapply(K, function(k){
    cat(k, "\n")
    readRDS(paste0(resDir, "/random_start_", rs, "/pre_obj_k", k, ".rds"))
  })
})

para_settings <- lapply(precast_results, function(pre_list){
  lapply(pre_list, function(pre_obj){
    attr(pre_obj@resList, 'para_settings')
  })
})

pen_const <- 100
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
ll_df$K <- factor(ll_df$K , levels = c(3:15))
plotDir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_precast", "BIC_select")

pdf(file.path(plotDir, "loglikelihood_vs_K.pdf"), height = 2.5, width = 5)
ggplot(ll_df, aes(x = K, y = LL_model, fill = K)) + geom_boxplot() + theme_classic() + ylab("Log-likelihood") + theme(legend.position = "none")+
xlab("Cluster resolution (K)")
dev.off()

pdf(file.path(plotDir, "MAIC_vs_K.pdf"), height = 2.5, width = 5)
ggplot(ll_df, aes(x = K, y = MAIC, fill = K)) + geom_boxplot() + theme_classic()+ ylab("MAIC")+
xlab("Cluster resolution (K)") + theme(legend.position = "none")
dev.off()

pdf(file.path(plotDir, "MBIC_vs_K.pdf"), height = 2.5, width = 5)
ggplot(ll_df, aes(x = K, y = MBIC, fill = K)) + geom_boxplot() + theme_classic()+ ylab("MBIC")+
xlab("Cluster resolution (K)") + theme(legend.position = "none")
dev.off()

pdf(file.path(plotDir, "BIC_vs_K.pdf"), height = 2.5, width = 5)
ggplot(ll_df, aes(x = K, y = BIC, fill = K)) + geom_boxplot()+ theme_classic()+ ylab("BIC")+
xlab("Cluster resolution (K)")+ theme(legend.position = "none") 
dev.off()

pdf(file.path(plotDir, "AIC_vs_K.pdf"), height = 2.5, width = 5)
ggplot(ll_df, aes(x = K, y = AIC, fill = K)) + geom_boxplot()+ theme_classic()+ ylab("AIC")+
xlab("Cluster resolution (K)") + theme(legend.position = "none") 
dev.off()

saveRDS(ll_df, paste0(resDir, "/resList_BIC.rds"))