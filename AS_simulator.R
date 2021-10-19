library(tidyverse)
library(parallel)
library(MLmetrics)
library(car)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path)) 
dir.create('traces', showWarnings = F)
dir.create('results', showWarnings = F)

gg.dpi = 600
gg.width = 8
gg.height = 6
gg.size = 20
TRACE = T

######## Functions

## Dev tools
get_dist <- function(coord1, coord2) unname(sqrt((coord1['X'] - coord2['X'])^2 + (coord1['Y'] - coord2['Y'])^2))

trace_SNPs.model <- function(SNPs.model){
  position = data.frame(
    X = c(50, 100, 0, 50),
    Y = c(big_tri_h, 0, 0, big_tri_h))
  
  ggplot(SNPs.model, aes(X, Y)) + 
    geom_polygon(data = position, color = 'black', fill = 'white') +
    geom_point(aes(color = type, size = n)) +
    geom_label(aes(label = note), nudge_y = 7) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

## General
get_SNPs <- function(samples, SNPs.model, freqs = NULL){
  
  if(is.null(freqs)) {
    #Give all pops dist from a SNP model
    
    freqs = lapply(1:nrow(SNPs.model), function(idx_model){
      get_pop_dists <- function(SNP.model) {
        pop.coordinates = rbind(data.frame(row.names = 'A', 'X' = 50, 'Y' = big_tri_h),
                                data.frame(row.names = 'B', 'X' = 0, 'Y' = 0),
                                data.frame(row.names = 'C', 'X' = 100, 'Y' = 0))
        
        dist_max = 100
        
        get_dist <- function(pop){
          sqrt((pop.coordinates[pop,'X'] - SNP.model[['X']])^2 + 
                 (pop.coordinates[pop,'Y'] - SNP.model[['Y']])^2)/dist_max}
        
        do.call(rbind, lapply(c('A', 'B', 'C'), function(pop) {
          data.frame(row.names = pop, dist_rel = get_dist(pop))
        }))
      }
      
      SNP.model = SNPs.model[idx_model,]
      #fst = get_fst(SNP.model)
      dists = get_pop_dists(SNP.model)
      #dists$rank = NA
      
      #dists[rank(dists$dist_rel, ties.method = 'min') == 1, 'rank'] = 'min'
      #dists[rank(dists$dist_rel, ties.method = 'min') == 2, 'rank'] = 'maj'
      

      #get_af_from_dists <- function(SNP.model, dists){
      ##No Fst, SNPs are exclusive
      if(SNP.model$type == 'CPE') {
        ##TODO give it the good form
        unlist(lapply(rownames(dists), function(pop) {
          ##The population having exlusive variants
          freqs.CPE = if(dists[pop,] == 0) {
            do.call(rbind, replicate(SNP.model$n, {
              af = sample(c(runif(1, 0.4, 1), runif(1, 0.4, 1), 1), 1)
              #af = runif(1, 0.8, 1)
              data.frame('HR' = (1 - af)^2, 'HET' = 2*(af*(1 - af)), 'HA' = af^2)
            }, simplify = F))
            
          ##The population not having the variant
          } else {
            do.call(rbind, replicate(SNP.model$n, data.frame('HR' = 1, 'HET' = 0, 'HA' = 0), simplify = F ))
          }
          rownames(freqs.CPE) = paste('SNP', SNP.model$note, 1:SNP.model$n, sep = '.')
          freqs.CPE = list(freqs.CPE)
          names(freqs.CPE) = pop
          freqs.CPE
        }), recursive = F)
        ##Then this freq end up in a list with names corresponding to SNP.model$note
      } else if(SNP.model$type == 'CPS') {
        dists[rank(dists, ties.method = 'min') == 1, 'dist_rel'] = 0.01
        dists[rank(dists, ties.method = 'min') == 2, 'dist_rel'] = 0.5
        
        draw_allele <- function(rf, dist) {
          #fst = dist/10
          fst = dist
          s1 = rf*(1-fst)/fst
          s2 = (1-rf)*(1-fst)/fst
          draw_af <-function(s1, s2) {rbeta(1, shape1 = s1, shape2 = s2)}
          af = 0
          while(af<rf) {af = draw_af(s1, s2)}
          af
        }
        
        #dists$rank = rank(dists$dist_rel, ties.method = 'min')
        afs = do.call(rbind, replicate(SNP.model$n, {
          rf = runif(1, 0, 0.8)
          data.frame('A' = draw_allele(rf, dists['A', 'dist_rel']),
                    'B' = draw_allele(rf, dists['B', 'dist_rel']),
                    'C' = draw_allele(rf, dists['C', 'dist_rel']))
        }, simplify = F))
        
        get_freqs_from_maf <- function(af) data.frame('HR' = (1 - af)^2, 'HET' = 2 * af * (1 - af), 'HA' = af^2)
        
        #Give row the correct name of SNPs
        name_SNPs <- function(freqs.CPS, SNP.model){
          row.names(freqs.CPS) = paste('SNP', SNP.model$note, 1:SNP.model$n, sep = '.')
          freqs.CPS
        }
        
        list('A' = name_SNPs(get_freqs_from_maf(afs$A), SNP.model),
             'B' = name_SNPs(get_freqs_from_maf(afs$B), SNP.model),
             'C' = name_SNPs(get_freqs_from_maf(afs$C), SNP.model))
        #}
      } else if(SNP.model$type == 'Neutral') {
        freqs.Neutral = do.call(rbind, lapply(1:SNP.model$n, function(null) {
            af = runif(1,0,1)
            data.frame('HR' = (1-af)^2, 'HET' = 2*(af*(1-af)), 'HA' = af^2)
        }))
        rownames(freqs.Neutral) = paste('SNP', SNP.model$note, 1:SNP.model$n, sep = '.')
        list('A' = freqs.Neutral, 'B' = freqs.Neutral, 'C' = freqs.Neutral)
      }
    })
  }
  names(freqs) = paste0('SNP.', SNPs.model$note)
  
  #Generate SNPs from freqs
  SNPs = do.call(cbind, lapply(1:nrow(SNPs.model), function(idx_model){
    SNP.model = SNPs.model[idx_model,]
    
    get_SNP_for_pop <- function(pop){
      data.frame(row.names = paste(pop, 1:samples[pop], sep = '.'),
                 do.call(cbind, lapply(1:SNP.model$n, function(idx_SNP){
                   sample(c(0, 1, 2), samples[pop], replace = T, 
                          freqs[[paste0('SNP.', SNP.model$note)]][[pop]][idx_SNP,])
                 })))
      }
    
    SNPs.pop = do.call(rbind, lapply(names(samples), get_SNP_for_pop))
    colnames(SNPs.pop) = paste('SNP', SNP.model$note, 1:SNP.model$n, sep = '.')
    SNPs.pop
  }))
  list(freqs = freqs, SNPs = SNPs)
}

get_phenotypes_from <- function(SNPs, causal.names, causal.betas, h){
  pheno.causal.M = t(t(SNPs[, causal.names]) * causal.betas)
  #do.call(rbind, strsplit(colnames(pheno.causal.M), split = '[.]', F))[,2:3]
  pheno.causal = rowSums(pheno.causal.M)
  pheno.noise = rnorm(n = nrow(SNPs), 0, sqrt((1 - h)/h) * sd(c(pheno.causal)))
  pheno = unlist(pheno.causal + pheno.noise)
  summary(lm(pheno ~ pheno.causal))$r.squared
  
  if(F){
    pheno.noise = rnorm(n = nrow(SNPs), 0, ((1 - sqrt(h) )/sqrt(h)) * sd(c(pheno.causal)))
    
    pheno.noise = rnorm(n = nrow(SNPs), 0, ((1 - h^2)/h^2) * sd(c(pheno.causal)))
    pheno.noise = rnorm(n = nrow(SNPs), 0, ((1 - h)/h)^2 * sd(c(pheno.causal)))
    pheno.noise = rnorm(n = nrow(SNPs), 0, ((1 - sqrt(h) )/sqrt(h) )^2 * sd(c(pheno.causal)))
    
    pheno.noise = rnorm(n = nrow(SNPs), 0, ((1 - h)/h) * sd(c(pheno.causal)))
    pheno = unlist(pheno.causal + pheno.noise)
    
    summary(lm(pheno ~ pheno.causal))$r.squared
  }
  
  if(F){
    pheno.causal.A = rowSums(pheno.causal.M[,grepl('SNP.A.[0-9]*$', causal.names, perl = T)])
    pheno.causal.A_B = rowSums(pheno.causal.M[,grepl('SNP.A_B.[0-9]*$', causal.names, perl = T)])
    
    ggplot(data.frame(pheno = pheno.causal.A, pop = as.factor(do.call(rbind, strsplit(names(pheno), split = '[.]', F))[,1]) )) + 
      geom_boxplot(aes(pop, pheno))
    
    ggplot(data.frame(pheno = pheno.causal.A_B, pop = as.factor(do.call(rbind, strsplit(names(pheno), split = '[.]', F))[,1]) )) + 
      geom_boxplot(aes(pop, pheno))
    
    ggplot(data.frame(pheno, pop = as.factor(do.call(rbind, strsplit(names(pheno), split = '[.]', F))[,1]) )) + 
      geom_boxplot(aes(pop, pheno))
  }
  
  pheno
}

get_phenotypes <- function(SNPs, h, causal_ratio, SNPs.model){
  causal.names = unlist(lapply(1:nrow(SNPs.model), function(idx){
    sample(paste0('SNP.', SNPs.model[idx,'note'], '.', 1:SNPs.model[idx,'n']), 
           causal_ratio * SNPs.model[idx,'n'])
  }))
  #causal.names = sample(colnames(SNPs), causal_ratio * ncol(SNPs))
  causal.betas = rnorm(n = length(causal.names), beta_mean, 1)
  
  pheno = get_phenotypes_from(SNPs, causal.names, causal.betas, h)
  list(pheno = data.frame(pheno), causal.betas = causal.betas, causal.names = causal.names)
}

get_eth <- function(samples){
  data.frame(ETH = unlist(lapply(names(samples), function(eth) rep(eth, samples[eth]))))
}

get_PCs <- function(SNPs){
  #ETH = get_eth(samples)
  pca = prcomp(SNPs)
  #PCs = cbind(pca$x, ETH)
  list(PCs = pca$x, loadings = pca$rotation)
}

get_projected_PCs_scale <- function(source, target, loadings){
  source.means=colMeans(source)
  source.sd=apply(source, 2, sd)
  target = scale(target, center = source.means, scale = source.sd)
  as.matrix(target) %*% loadings
}

get_projected_PCs <- function(target, loadings, nb_PCs){
  PCs = as.matrix(target) %*% loadings
  PCs[,1:nb_PCs]
}

GWAS <- function(pheno, SNPs, PCs = NULL){
  test_SNP <- function(form, data){
    lm.fitted = lm(form, data)
    al = alias(lm.fitted)
    slm = summary(lm.fitted)
    vif.SNP = ifelse(length(al) > 1 , 1000, vif(lm.fitted)[1])
    data.frame('beta' = slm$coefficients[2,1], 
               'pval' = slm$coefficients[2,4], 
               'min.log.pval' = - log(slm$coefficients[2,4]),
               'vif' = vif.SNP)
}
  
  coef = if(!is.null(PCs)) {
    data.for_lm = cbind(PCs, pheno$pheno, SNPs)
    coef = do.call(rbind, lapply(colnames(SNPs), function(SNP.name) {
      test_SNP(as.formula(paste0(paste0('pheno ~ ', SNP.name, ' + '), 
                                 paste0('PC', 1:ncol(PCs), collapse = ' + '))), data.for_lm)
    }))
  } else {
    data.for_lm = cbind(pheno$pheno, SNPs)
    coef = do.call(rbind, lapply(colnames(SNPs), function(SNP.name) {
      test_SNP(as.formula(paste0('pheno ~ ', SNP.name)), data.for_lm)
    }))
  }
  coef$idx = 1:nrow(coef)
  coef$causal = colnames(SNPs) %in% pheno$causal.names
  coef$SNP.name = colnames(SNPs)
  coef
}

calc_PRS <- function(SNPs, GWAS.coef, pheno, thresholding, SNPs_in_PRS = NULL){
  #If needed updated GWAS coeff with SNPs in PRS, code can be uncommented
  if(!thresholding) {
    #PRS = rowSums(t(t(SNPs[,GWAS.coef$causal, drop= F]) * filter(GWAS.coef, causal)$beta))
    #SNPs_in_PRS = GWAS.coef$SNP.name[GWAS.coef$causal]
    #list(PRS, SNPs_in_PRS)
    SNPs_in_PRS.mask = GWAS.coef$SNP.name %in% GWAS.coef$SNP.name[GWAS.coef$causal]
    PRS = rowSums(t(t(SNPs[,SNPs_in_PRS.mask, drop= F]) * GWAS.coef[SNPs_in_PRS.mask,]$beta))
    SNPs_in_PRS = GWAS.coef$SNP.name[SNPs_in_PRS.mask]
    list(PRS, SNPs_in_PRS)
    list('PRS' = PRS, 'SNPs_in_PRS' = SNPs_in_PRS, 'SNPs_out_PRS' = GWAS.coef$SNP.name[!GWAS.coef$causal])
    
  } else if(thresholding & is.null(SNPs_in_PRS)) {
    thresholds = c(10 %o% 10^(-1:-100))
    PRSs = do.call(cbind, lapply(thresholds, function(threshold){
      mask = (GWAS.coef$pval < threshold) & (GWAS.coef$vif < 10)
      tibble(threshold = rowSums(t(t(SNPs[,mask, drop= F]) * GWAS.coef[mask, 'beta'])))
    }))
    colnames(PRSs) = thresholds
    r_squareds = unlist(lapply(as.character(thresholds), function(threshold){
      summary(lm(pheno ~ PRSs[,threshold]))$r.squared
    }))
    PRS = PRSs[,which.max(r_squareds)]
    #SNPs_out_PRS = GWAS.coef$SNP.name[(GWAS.coef$pval > thresholds[which.max(r_squareds)]) | (GWAS.coef$vif > 10)]
    SNPs_in_PRS = GWAS.coef$SNP.name[(GWAS.coef$pval < thresholds[which.max(r_squareds)]) & (GWAS.coef$vif < 10)]
    SNPs_out_PRS = GWAS.coef$SNP.name[!GWAS.coef$SNP.name %in% SNPs_in_PRS]
    
    list('PRS' = PRS, 'SNPs_in_PRS' = SNPs_in_PRS, 'SNPs_out_PRS' = SNPs_out_PRS, 'PRS.p_value' = thresholds[which.max(r_squareds)] )
  } else {
    SNPs_in_PRS.mask = GWAS.coef$SNP.name %in% SNPs_in_PRS
    PRS = rowSums(t(t(SNPs[,SNPs_in_PRS.mask, drop= F]) * GWAS.coef[SNPs_in_PRS.mask,]$beta))
    list('PRS' = PRS)
  }
}

PCAS <- function(pheno, PCs){
  PCs.pheno = cbind(PCs, pheno)
  slm.PCs = summary(lm(as.formula(paste0('pheno ~ ', paste(paste0('PC', 1:ncol(PCs)), collapse = ' + '))), PCs.pheno))
  data.frame('beta' = slm.PCs$coefficients[2:(ncol(PCs)+1),1], 
             'pval' = slm.PCs$coefficients[2:(ncol(PCs)+1),4], 
             'min.log.pval' = -log(slm.PCs$coefficients[2:(ncol(PCs)+1),4]))
}

calc_PCS <- function(PCs, PCAS.coeff){
  rowSums((t(t(PCs) * PCAS.coeff$beta)))
}

## Trace
lambda_GWAS <- function(pvalue) median(qchisq(1-pvalue,1), na.rm = T)/qchisq(0.5,1)

plot_GWAS_QQ <- function(coef){
  y = -log10(sort(coef$pval, decreasing=F))
  x = -log10( 1:length(y)/length(y) )
  p <- ggplot(data.frame(x,y), aes(x, y))
  p + geom_point() + 
    labs(title=paste('QQ, L=', round(lambda_GWAS(coef$pval),digits = 3)), x = "Expected", y = "Observed") + 
    geom_abline(intercept = 0, colour ="red") + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
}

plot_GWAS_manhattan <- function(res_GWAS, filename, cal_PRS.p_value){
  filename = gsub(pattern = 'rep', replacement = 'GWAS.rep', filename)
  res_GWAS$SNP_category = unlist(lapply(strsplit(res_GWAS$SNP.name, '\\.'), function(ele) {ele[2]}))
  res_GWAS$min.log.pval[res_GWAS$min.log.pval > 100] = 100
  p <- ggplot(res_GWAS, aes(idx, min.log.pval))
  p + geom_point(aes(shape  = causal, color = SNP_category)) + geom_hline(yintercept = -log10(cal_PRS.p_value), colour="black")
  ggsave(filename = paste0('traces/', filename, '.png'), dpi = gg.dpi, width = gg.width, height = gg.height, units = 'in')
}

plot_PCs <- function(PCs, samples, filename){
  nb_PCs = ncol(PCs)
  PCs = data.frame(PCs, get_eth(samples))
  PC_couple  = if(nb_PCs %% 2 == 1) {c(1:nb_PCs, nb_PCs - 1)} else {1:nb_PCs}
  PC_couple = t(matrix(PC_couple, nrow = 2))
  PC_couple = data.frame(PCA = paste0('PC',  PC_couple[,1]),  PCB =  paste0('PC',  PC_couple[,2]), stringsAsFactors = F)
  lapply(1:nrow(PC_couple), function(idx){
    ggplot(PCs) + geom_point(aes_string(PC_couple[idx,'PCA'], PC_couple[idx,'PCB'], colour = 'ETH'))
    ggsave(filename = paste0('traces/', 
                             filename = gsub(pattern = 'rep', replacement = paste0(PC_couple[idx,'PCA'], '_', PC_couple[idx,'PCB'], '.rep'), filename), '.png'), dpi = gg.dpi, width = gg.width, height = gg.height, units = 'in')
  })
}

plot_pheno <- function(pheno, samples, filename){
  filename = gsub(pattern = 'rep', replacement = 'pheno.rep', filename)
  
  pheno = cbind(pheno, get_eth(samples))
  ggplot(pheno) + geom_density(aes(pheno, colour = ETH))
  ggsave(filename = paste0('traces/', filename, '.png'), dpi = gg.dpi, width = gg.width, height = gg.height, units = 'in')
}

write_freqs <- function(freqs, filename){
  filename = gsub(pattern = 'rep', replacement = 'freqs.rep', filename)
  
  freqs.table = do.call(rbind, lapply(freqs, function(SNP.category){ 
    data.frame(B = sqrt(SNP.category$B$HA), 
               A = sqrt(SNP.category$A$HA), 
               C = sqrt(SNP.category$C$HA))      
  }))
  write.table(freqs.table, paste0('traces/', filename, '.tsv'), row.names = T, col.names = T, quote = F, sep = '\t')
}

## Simulations
simul.runner <- function(body, nb_repetitions, parallel){
  if(parallel){
    print(paste('Starting parallel session with', nb_repetitions, 'repetitions/CPU(s)'))
    cl = makeCluster(nb_repetitions, type = 'FORK')
    res = do.call(rbind, clusterApply(cl, 1:nb_repetitions, body))
    stopCluster(cl)
    res    
  } else {
    print(paste('Starting sequential session with', nb_repetitions, 'repetitions'))
    do.call(rbind, lapply(1:nb_repetitions, body))
  }
}

PCAS_simulate_internal <- function(SNPs.model, samples, scenario, causal_ratio, PRS_thresholding, 
                                   nb_repetitions, parallel = FALSE){
  body <- function(repetition_nb) {
    print(repetition_nb)
    
    #SNPs.model = SNPs.CPE
    #SNPs.model = SNPs.CPSS.CPE
    #samples = samples.MIX.imb
    
    #cohort = ifelse(length(table(samples)) == 1, 'Balanced', 'Imbalanced')
    cohort = 'Imbalanced'
    
    base.SNPs = get_SNPs(samples, SNPs.model)
    SNPs.freqs = base.SNPs$freqs
    base.SNPs = base.SNPs$SNPs
    
    base.pheno = get_phenotypes(base.SNPs, h, causal_ratio, SNPs.model)
    base.pca = get_PCs(base.SNPs)
    base.PCs = base.pca$PCs[,1:nb_PCs]
    
    GWAS.coef = GWAS(base.pheno, base.SNPs, base.PCs)
    PCAS.coef = PCAS(base.pheno$pheno, base.PCs)
    
    calibration.SNPs = get_SNPs(samples, SNPs.model, freqs = SNPs.freqs)$SNPs
    calibration.pheno = get_phenotypes_from(calibration.SNPs, base.pheno$causal.names, base.pheno$causal.betas, h)
    calibration.PCs = get_projected_PCs(calibration.SNPs, base.pca$loadings, nb_PCs)
    
    calibration.PRS = calc_PRS(calibration.SNPs, GWAS.coef, calibration.pheno, PRS_thresholding)
    cal_PRS.p_value = ifelse(is.null(calibration.PRS$PRS.p_value), 0, calibration.PRS$PRS.p_value)
    SNPs_in_PRS = calibration.PRS$SNPs_in_PRS
    calibration.PRS = calibration.PRS$PRS
    calibration.PCS = calc_PCS(calibration.PCs, PCAS.coef)
    
    calibration.predictors = data.frame(pheno = calibration.pheno, PRS = calibration.PRS, PCS = calibration.PCS)
    PRS.lm = lm(pheno ~ PRS, calibration.predictors)
    PCS.lm = lm(pheno ~ PCS, calibration.predictors)
    PCS_PRS.lm = lm(pheno ~ PCS + PRS, calibration.predictors)
    
    target.SNPs = get_SNPs(samples, SNPs.model, freqs = SNPs.freqs)$SNPs
    target.pheno = get_phenotypes_from(target.SNPs, base.pheno$causal.names, base.pheno$causal.betas, h)
    target.PCs = get_projected_PCs(target.SNPs, base.pca$loadings, nb_PCs)
    
    target.PRS = calc_PRS(target.SNPs, GWAS.coef, target.pheno, PRS_thresholding, SNPs_in_PRS)$PRS
    target.PCS = calc_PCS(target.PCs, PCAS.coef)
    
    target.predictors = data.frame(pheno = target.pheno, PRS = target.PRS, PCS = target.PCS)
    
    target.PRS.Yb = predict(PRS.lm, target.predictors)
    target.PCS.Yb = predict(PCS.lm, target.predictors)
    target.PCS_PRS.Yb = predict(PCS_PRS.lm, target.predictors)
    
    target.PRS.R2 = R2_Score(target.PRS.Yb, target.pheno)
    target.PCS.R2 = R2_Score(target.PCS.Yb, target.pheno)
    target.PCS_PRS.R2 = R2_Score(target.PCS_PRS.Yb, target.pheno)
    
    heritability = summary(lm(target.pheno ~ as.matrix(target.SNPs[base.pheno$causal.names])))$adj.r.squared
    #tt = SNPs.freqs[[1]]
    #View(cbind(tt$A, tt$B, tt$C))
    
    filename = NA
    if(TRACE){
      filename = paste('S.Internal', gsub(' ', '_', scenario), cohort, paste0('rep_',repetition_nb), sep = '.')
      plot_GWAS_manhattan(GWAS.coef, filename, cal_PRS.p_value)
      plot_pheno(base.pheno$pheno, samples, filename)
      plot_PCs(base.PCs, samples, filename)
      write_freqs(SNPs.freqs, filename)
    }
    
    data.frame(#PRS.R2_BF = summary(PRS.lm)$adj.r.squared, PCS.R2_BF = summary(PCS.lm)$adj.r.squared, PCS_PRS.R2_BF= summary(PCS_PRS.lm)$adj.r.squared,
      target.PRS.R2, target.PCS.R2, target.PCS_PRS.R2, 
      target.PCS.R2_gain = target.PCS_PRS.R2 - target.PRS.R2,
      target.PCS.heritability_gain = heritability - target.PCS_PRS.R2,
      cal_PRS.p_value, heritability,
      Cohort = cohort, Scenario = scenario, h, causal_ratio, repetition_nb, PRS_thresholding, filename)
  }
  
  simul.runner(body, nb_repetitions, parallel)
}

PCAS_simulate_external <- function(SNPs.model, samples.base, samples.calibration, samples.target, 
                                   scenario, causal_ratio, PRS_thresholding, nb_repetitions, parallel = FALSE){
  body <- function(repetition_nb){
    print(repetition_nb)
    #cohort.base = ifelse(length(table(samples.base)) == 1, 'Balanced', 'Imbalanced')
    #cohort.calibration = ifelse(length(table(samples.calibration)) == 1, 'Balanced', 'Imbalanced')
    #cohort.target = ifelse(length(table(samples.target)) == 1, 'Balanced', 'Imbalanced')
    cohort.base = 'Imbalanced'
    cohort.calibration = 'Imbalanced'
    cohort.target = 'Imbalanced'
    
    ##base
    #Training: data generation
    base.SNPs = get_SNPs(samples.base, SNPs.model)
    SNPs.freqs = base.SNPs$freqs
    base.SNPs = base.SNPs$SNPs
    base.pheno = get_phenotypes(base.SNPs, h, causal_ratio, SNPs.model)
    base.pca = get_PCs(base.SNPs)
    base.PCs = base.pca$PCs[,1:nb_PCs]
    
    #Training: association tests
    GWAS.coef = GWAS(base.pheno, base.SNPs, base.PCs)
    PCAS.coef = PCAS(base.pheno$pheno, base.PCs)
    
    ##Calibration
    #Calibration: data generation
    calibration.SNPs = get_SNPs(samples.calibration, SNPs.model, freqs = SNPs.freqs)$SNPs
    calibration.pheno = get_phenotypes_from(calibration.SNPs, base.pheno$causal.names, base.pheno$causal.betas, h)
    calibration.PCs = get_projected_PCs(calibration.SNPs, base.pca$loadings, nb_PCs)
    
    #Calibration: predictor generation
    calibration.PRS = calc_PRS(calibration.SNPs, GWAS.coef, calibration.pheno, PRS_thresholding)
    cal_PRS.p_value = ifelse(is.null(calibration.PRS$PRS.p_value), 0, calibration.PRS$PRS.p_value)
    calibration.PCS = calc_PCS(calibration.PCs, PCAS.coef)
    calibration.predictors = data.frame(pheno = calibration.pheno, PRS = calibration.PRS$PRS, PCS = calibration.PCS)
    
    #Calibration: model calibration
    calibration.lm.PRS = lm(pheno ~ PRS, calibration.predictors)
    calibration.lm.PCS = lm(pheno ~ PCS, calibration.predictors)
    calibration.lm.PCS_PRS = lm(pheno ~ PCS + PRS, calibration.predictors)
    
    calibration.slm.PRS = summary(calibration.lm.PRS)
    calibration.slm.PCS = summary(calibration.lm.PCS)
    calibration.slm.PCS_PRS = summary(calibration.lm.PCS_PRS)
    
    #Calibration: get coefficients
    calibration.PRS.intercept = calibration.slm.PRS$coeff[1,1]
    calibration.PRS.beta = calibration.slm.PRS$coeff[2,1]
    
    calibration.PCS.intercept =  calibration.slm.PCS$coeff[1,1]
    calibration.PCS.beta =  calibration.slm.PCS$coeff[2,1]
    
    calibration.PCS_PRS.intercept = calibration.slm.PCS_PRS$coeff[1,1]
    calibration.PCS_PRS.PCS.beta = calibration.slm.PCS_PRS$coeff[2,1]
    calibration.PCS_PRS.PRS.beta = calibration.slm.PCS_PRS$coeff[3,1]
    
    
    ##External
    #External: data generation
    target.SNPs = get_SNPs(samples.target, SNPs.model, freqs = SNPs.freqs)$SNPs
    target.pheno = get_phenotypes_from(target.SNPs, base.pheno$causal.names, base.pheno$causal.betas, h)
    target.PCs = get_projected_PCs(target.SNPs, base.pca$loadings, nb_PCs)
    
    #External: predictor generation
    target.PRS = calc_PRS(target.SNPs, GWAS.coef, target.pheno, PRS_thresholding, calibration.PRS$SNPs_in_PRS)$PRS
    target.PCS = calc_PCS(target.PCs, PCAS.coef)
    
    target.predictors = data.frame(pheno = target.pheno, PRS = target.PRS, PCS = target.PCS)
    
    #External: scoring
    target.PRS.Yb = predict(calibration.lm.PRS, target.predictors)
    target.PCS.Yb = predict(calibration.lm.PCS, target.predictors)
    target.PCS_PRS.Yb = predict(calibration.lm.PCS_PRS, target.predictors)
    
    target.PRS.R2 = R2_Score(target.PRS.Yb, target.pheno)
    target.PCS.R2 = R2_Score(target.PCS.Yb, target.pheno)
    target.PCS_PRS.R2 = R2_Score(target.PCS_PRS.Yb, target.pheno)
    
    target.slm.PRS = summary(lm(pheno ~ PRS, target.predictors))
    target.slm.PCS = summary(lm(pheno ~ PCS, target.predictors))
    target.slm.PCS_PRS = summary(lm(pheno ~ PCS + PRS, target.predictors))
    
    target.PRS.R2_BF = target.slm.PRS$adj.r.squared
    target.PCS.R2_BF = target.slm.PCS$adj.r.squared
    target.PCS_PRS.R2_BF = target.slm.PCS_PRS$adj.r.squared
    
    #External: get coefficients
    target.PRS.intercept = target.slm.PRS$coeff[1,1]
    target.PRS.beta = target.slm.PRS$coeff[2,1]
    
    target.PCS.intercept =  target.slm.PCS$coeff[1,1]
    target.PCS.beta =  target.slm.PCS$coeff[2,1]
    
    target.PCS_PRS.intercept = target.slm.PCS_PRS$coeff[1,1]
    target.PCS_PRS.PCS_beta = target.slm.PCS_PRS$coeff[2,1]
    target.PCS_PRS.PRS_beta = target.slm.PCS_PRS$coeff[3,1]
    
    target.PCS_gain.R2 = target.PCS_PRS.R2 - target.PRS.R2
    PCS_gain.intercept = target.PCS_PRS.intercept - target.PRS.intercept
    PCS_gain.PRS_beta = target.PCS_PRS.PRS_beta - target.PRS.beta
    
    filename = NA
    if(TRACE){
      filename = paste('S.External', gsub(pattern = ' ', replacement = '_', scenario), cohort.base, paste0('rep_',repetition_nb), sep = '.') 
      
      plot_GWAS_manhattan(GWAS.coef, filename, cal_PRS.p_value)
      write_freqs(SNPs.freqs, filename)
      
      plot_pheno(base.pheno$pheno, samples.base, paste(filename, 'base', sep = '.') )
      plot_pheno(calibration.pheno, samples.calibration, paste(filename, 'calibration', sep = '.') )
      plot_pheno(target.pheno, samples.target, paste(filename, 'target', sep = '.') )
      
      plot_PCs(base.PCs, samples.base,  paste(filename, 'base', sep = '.'))
      plot_PCs(calibration.PCs, samples.calibration,  paste(filename, 'calibration', sep = '.'))
      plot_PCs(target.PCs, samples.target,  paste(filename, 'target', sep = '.'))
    }
    
    target.heritability = summary(lm(target.pheno ~ as.matrix(target.SNPs[base.pheno$causal.names])))$adj.r.squared
    target.heritability_gain.R2 = target.heritability - target.PCS_PRS.R2
    
    target.PRS.R2_share = target.PRS.R2/target.heritability
    target.PCS.R2_share = target.PCS.R2/target.heritability
    target.PCS_PRS.R2_share = target.PCS_PRS.R2/target.heritability
    target.PCS_gain.R2_share = target.PCS_PRS.R2_share - target.PRS.R2_share
    target.heritability_gain.R2_share = (target.heritability - target.PCS_PRS.R2)/target.heritability
    print(target.PRS.R2, target.PCS_PRS.R2,  target.heritability)
    print(target.PRS.R2_share, target.PCS_PRS.R2_share,  target.heritability_gain.R2_share)
    #Return
    data.frame(calibration.PRS.intercept, calibration.PRS.beta, calibration.PCS.intercept, calibration.PCS.beta,
               calibration.PCS_PRS.intercept, calibration.PCS_PRS.PCS.beta, calibration.PCS_PRS.PRS.beta,
               target.PRS.R2_BF, target.PRS.intercept, target.PRS.beta, target.PRS.R2, 
               target.PCS.R2_BF, target.PCS.intercept, target.PCS.beta, target.PCS.R2, 
               target.PCS_PRS.R2_BF, target.PCS_PRS.intercept, target.PCS_PRS.PCS_beta, target.PCS_PRS.PRS_beta, target.PCS_PRS.R2,
               target.PCS_gain.R2, PCS_gain.intercept, PCS_gain.PRS_beta, 
               target.PRS.R2_share, target.PCS.R2_share, target.PCS_PRS.R2_share, target.PCS_gain.R2_share, target.heritability_gain.R2_share,
               target.heritability, target.heritability_gain.R2, Cohort = cohort.base, cohort.base, cohort.calibration, cohort.target, Scenario = scenario, h, causal_ratio, repetition_nb, PRS_thresholding, filename)
  }
  simul.runner(body, nb_repetitions, parallel)
}
######## Set constants
set_constants <-function(){
  
  #######Samples
  ##Simulation1 and ##Simulation2
  #samples.MIX.bal <<- c('A' = 10000, 'B' = 10000, 'C' = 10000)
  samples.MIX.imb <<- c('A' = 20000, 'B' = 7000, 'C' = 3000)
  samples.MIX.imb2 <<- c('A' = 3000, 'B' = 7000, 'C' = 20000)
  
  ##Simulation2 specific
  #samples.A.C.bal <<- c('A' = 15000, 'C' = 15000)
  samples.A.C.imb <<- c('A' = 23500, 'C' = 6500)
  samples.B <<- c('B' = 30000)
  samples.A <<- c('A' = 30000)
  samples.C <<- c('C' = 30000)
  
  #######SNPs distributions
  nb_PCs <<- 10
  h <<- 0.6
  big_tri_h <<- sqrt(100^2 - 50^2)
  small_tri_h <<- sqrt(25^2 - 12.5^2)
  
  center <<- c(Y = big_tri_h/3, X = 50)
  
  #beta_mean <<- 0.2
  beta_mean <<- 0
  
  B.CPE <<- c(X = 0, Y = 0)
  A.CPE <<- c(X = 50, Y = big_tri_h)  
  C.CPE <<- c(X = 100, Y = 0)
  
  A_B.CPE_CPS <<- c(X = 25, Y = big_tri_h/2)  
  A_C.CPE_CPS <<- c(X = 75, Y = big_tri_h/2)
  B_C.CPE_CPS <<- c(X = 50, Y = 0)
  
  A_B.CPS <<- c(X = 50 - 12.5, Y = (big_tri_h * 1/6) + small_tri_h)  
  A_C.CPS <<- c(X = 50 + 12.5, Y = (big_tri_h * 1/6) + small_tri_h)
  B_C.CPS <<- c(X = 50, Y = big_tri_h * 1/6)
  
  SNPs.control <<- data.frame('n' = 450, 'Y' = center['Y'], 'X' = center['X'], 'note' = 'centered', type = 'Neutral', dample = T)
  
  if(BALANCED){
    SNPs.CPE <<- rbind(data.frame('n' = 50, 'Y' = A.CPE['Y'], 'X' = A.CPE['X'], 'note' = 'A', type = 'CPE', dample = T),
                       data.frame('n' = 50, 'Y' = C.CPE['Y'], 'X' = C.CPE['X'], 'note' = 'C', type = 'CPE', dample = T),
                       data.frame('n' = 50, 'Y' = B.CPE['Y'], 'X' = B.CPE['X'], 'note' = 'B', type = 'CPE', dample = T),
                       data.frame('n' = 350, 'Y' = center['Y'], 'X' = center['X'], 'note' = 'centered', type = 'Neutral', dample = T))
    
    SNPs.CPS <<- rbind(data.frame('n' = 50, 'Y' = A_B.CPS['Y'], 'X' = A_B.CPS['X'], 'note' = 'A_B', type = 'CPS', dample = T),
                       data.frame('n' = 50, 'Y' = A_C.CPS['Y'], 'X' = A_C.CPS['X'], 'note' = 'A_C', type = 'CPS', dample = T),
                       data.frame('n' = 50, 'Y' = B_C.CPS['Y'], 'X' = B_C.CPS['X'], 'note' = 'B_C', type = 'CPS', dample = T),
                       data.frame('n' = 350, 'Y' = center['Y'], 'X' = center['X'], 'note' = 'centered', type = 'Neutral', dample = T))
    
    SNPs.CPE.CPS <<- rbind(data.frame('n' = 50, 'Y' = A.CPE['Y'], 'X' = A.CPE['X'], 'note' = 'A', type = 'CPE', dample = T),
                           data.frame('n' = 50, 'Y' = C.CPE['Y'], 'X' = C.CPE['X'], 'note' = 'C', type = 'CPE', dample = T),
                           data.frame('n' = 50, 'Y' = B.CPE['Y'], 'X' = B.CPE['X'], 'note' = 'B', type = 'CPE', dample = T),
                           data.frame('n' = 50, 'Y' = A_B.CPS['Y'], 'X' = A_B.CPS['X'], 'note' = 'A_B', type = 'CPS', dample = T),
                           data.frame('n' = 50, 'Y' = A_C.CPS['Y'], 'X' = A_C.CPS['X'], 'note' = 'A_C', type = 'CPS', dample = T),
                           data.frame('n' = 50, 'Y' = B_C.CPS['Y'], 'X' = B_C.CPS['X'], 'note' = 'B_C', type = 'CPS', dample = T),
                           data.frame('n' = 200, 'Y' = center['Y'], 'X' = center['X'], 'note' = 'centered', type = 'Neutral', dample = T))
  } else {
    SNPs.CPE <<- rbind(data.frame('n' = 25, 'Y' = A.CPE['Y'], 'X' = A.CPE['X'], 'note' = 'A', type = 'CPE', dample = T),
                       data.frame('n' = 50, 'Y' = C.CPE['Y'], 'X' = C.CPE['X'], 'note' = 'C', type = 'CPE', dample = T),
                       data.frame('n' = 75, 'Y' = B.CPE['Y'], 'X' = B.CPE['X'], 'note' = 'B', type = 'CPE', dample = T),
                       data.frame('n' = 350, 'Y' = center['Y'], 'X' = center['X'], 'note' = 'centered', type = 'Neutral', dample = T))
    
    SNPs.CPS <<- rbind(data.frame('n' = 25, 'Y' = A_B.CPS['Y'], 'X' = A_B.CPS['X'], 'note' = 'A_B', type = 'CPS', dample = T),
                       data.frame('n' = 50, 'Y' = A_C.CPS['Y'], 'X' = A_C.CPS['X'], 'note' = 'A_C', type = 'CPS', dample = T),
                       data.frame('n' = 75, 'Y' = B_C.CPS['Y'], 'X' = B_C.CPS['X'], 'note' = 'B_C', type = 'CPS', dample = T),
                       data.frame('n' = 350, 'Y' = center['Y'], 'X' = center['X'], 'note' = 'centered', type = 'Neutral', dample = T))
    
    SNPs.CPE.CPS <<- rbind(data.frame('n' = 25, 'Y' = A.CPE['Y'], 'X' = A.CPE['X'], 'note' = 'A', type = 'CPE', dample = T),
                           data.frame('n' = 50, 'Y' = C.CPE['Y'], 'X' = C.CPE['X'], 'note' = 'C', type = 'CPE', dample = T),
                           data.frame('n' = 75, 'Y' = B.CPE['Y'], 'X' = B.CPE['X'], 'note' = 'B', type = 'CPE', dample = T),
                           data.frame('n' = 25, 'Y' = A_B.CPS['Y'], 'X' = A_B.CPS['X'], 'note' = 'A_B', type = 'CPS', dample = T),
                           data.frame('n' = 50, 'Y' = A_C.CPS['Y'], 'X' = A_C.CPS['X'], 'note' = 'A_C', type = 'CPS', dample = T),
                           data.frame('n' = 75, 'Y' = B_C.CPS['Y'], 'X' = B_C.CPS['X'], 'note' = 'B_C', type = 'CPS', dample = T),
                           data.frame('n' = 200, 'Y' = center['Y'], 'X' = center['X'], 'note' = 'centered', type = 'Neutral', dample = T))
  }
  

  
}



######## Run simulations

## Simulation 1 (Simulations for different trans-ancestry genetic architectures)
BALANCED=T
set_constants()
nb_repetitions = 30
para = F
PRS_thresholding = T
causal_ratio = 0.3
nb_PCs = 2
 
res.PCAS_simulate.internal = rbind(
  PCAS_simulate_internal(SNPs.CPE.CPS, samples.MIX.imb, '1', causal_ratio, PRS_thresholding, nb_repetitions, para),
  PCAS_simulate_internal(SNPs.CPE, samples.MIX.imb, '2', causal_ratio, PRS_thresholding, nb_repetitions, para),
  PCAS_simulate_internal(SNPs.CPS, samples.MIX.imb, '3', causal_ratio, PRS_thresholding, nb_repetitions, para),
  PCAS_simulate_internal(SNPs.control, samples.MIX.imb, 'Control', causal_ratio, PRS_thresholding, nb_repetitions, para))

res.PCAS_simulate.internal.not.balanced = res.PCAS_simulate.internal
res.PCAS_simulate.internal = res.PCAS_simulate.internal %>% mutate(Scenario = fct_relevel(Scenario, 'Control', '1', '2', '3'))

res.PCAS_simulate.internal.to_plot = select(res.PCAS_simulate.internal, Scenario, target.PRS.R2, target.PCS.R2_gain, target.PCS_PRS.R2, target.PCS.heritability_gain, heritability)
res.PCAS_simulate.internal.to_plot_boxplots = do.call(rbind, lapply(levels(res.PCAS_simulate.internal.to_plot$Scenario), function(scenari){
  ele = filter(res.PCAS_simulate.internal.to_plot, Scenario == scenari)
  data.frame(Scenario = scenari, 
             Heritability = median(ele$target.PCS.heritability_gain),
             Heritability.CI = diff(boxplot.stats(ele$heritability)$stats[c(1,3)])/2,
             PGS = median(ele$target.PRS.R2),
             PGS.CI = diff(boxplot.stats(ele$target.PRS.R2)$stats[c(1,3)])/2,
             PCS_gain = median(ele$target.PCS.R2_gain),
             PCS_gain.CI = diff(boxplot.stats(ele$target.PCS.R2)$stats[c(1,3)])/2, stringsAsFactors = T)
}))

res.PCAS_simulate.internal.to_plot = res.PCAS_simulate.internal.to_plot_boxplots
res.PCAS_simulate.internal.to_plot = cbind(pivot_longer(select(res.PCAS_simulate.internal.to_plot, Scenario, PGS, PCS_gain, Heritability), cols = c('PGS','PCS_gain', 'Heritability'), names_to = 'model'),
                                           pivot_longer(select(res.PCAS_simulate.internal.to_plot, PGS.CI, PCS_gain.CI, Heritability.CI), cols = c('PGS.CI', 'PCS_gain.CI', 'Heritability.CI'), values_to = 'CI'))
res.PCAS_simulate.internal.to_plot$name = NULL
res.PCAS_simulate.internal.to_plot$model = gsub('_', ' ', res.PCAS_simulate.internal.to_plot$model)

res.PCAS_simulate.internal.to_plot$label_pos = NA
res.PCAS_simulate.internal.to_plot[res.PCAS_simulate.internal.to_plot$model == 'PGS',]$label_pos = 
  res.PCAS_simulate.internal.to_plot[res.PCAS_simulate.internal.to_plot$model == 'PGS',]$value
res.PCAS_simulate.internal.to_plot[res.PCAS_simulate.internal.to_plot$model == 'PCS gain',]$label_pos = 
  res.PCAS_simulate.internal.to_plot[res.PCAS_simulate.internal.to_plot$model == 'PGS',]$value + 
  res.PCAS_simulate.internal.to_plot[res.PCAS_simulate.internal.to_plot$model == 'PCS gain',]$value
res.PCAS_simulate.internal.to_plot[res.PCAS_simulate.internal.to_plot$model == 'Heritability',]$label_pos = 
  res.PCAS_simulate.internal.to_plot[res.PCAS_simulate.internal.to_plot$model == 'Heritability',]$value +
  res.PCAS_simulate.internal.to_plot[res.PCAS_simulate.internal.to_plot$model == 'PGS',]$value + 
  res.PCAS_simulate.internal.to_plot[res.PCAS_simulate.internal.to_plot$model == 'PCS gain',]$value

gg.size.E =  gg.size + gg.size * 2/3
legend_title=''
res.PCAS_simulate.internal.to_plot = filter(res.PCAS_simulate.internal.to_plot, model != 'Heritability')

p = ggplot(data=res.PCAS_simulate.internal.to_plot, aes(x=Scenario, y=value, fill=model)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=label_pos - CI, ymax = label_pos + CI), width=.2,
                position=position_dodge(.9), color = "#999999") +
  labs(x = "Scenario", y = "PVE", fill = legend_title) +
  scale_fill_manual(labels = list(expression(hat(Y)[paste('PGS+AS')]), 
                                  expression(hat(Y)[PGS])), values = c("#0072B2", "#56B4E9") ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size=gg.size.E),
        axis.title.x=element_blank());p
ggsave(plot = p, filename = paste0('results/sim.I.results.png'), dpi = gg.dpi, width = gg.width, height = gg.height, units = 'in')  

## Simulation 2: Simulations of generalization
BALANCED=F
set_constants()
nb_repetitions = 50
para = F
causal_ratio = 0.3
PRS_thresholding = T
res.PCAS_simulate.external = rbind(
  PCAS_simulate_external(SNPs.CPE.CPS, samples.MIX.imb, samples.MIX.imb, samples.MIX.imb, '1', causal_ratio, PRS_thresholding, nb_repetitions, para),
  PCAS_simulate_external(SNPs.CPE.CPS, samples.MIX.imb, samples.MIX.imb, samples.MIX.imb2, '2', causal_ratio, PRS_thresholding, nb_repetitions, para),
  PCAS_simulate_external(SNPs.CPE.CPS, samples.A, samples.A, samples.A, 'Control', causal_ratio, PRS_thresholding, nb_repetitions, para))


res.PCAS_simulate.external.bck = res.PCAS_simulate.external
res.PCAS_simulate.external$Scenario = as.factor(res.PCAS_simulate.external$Scenario)
res.PCAS_simulate.external = res.PCAS_simulate.external %>% mutate(Scenario = fct_relevel(Scenario, 'Control', '1', '2'))

res.PCAS_simulate.external.to_plot = select(res.PCAS_simulate.external, Scenario, 
                                            target.PRS.R2_share, target.PCS_gain.R2_share, target.PCS_PRS.R2_share, target.heritability_gain.R2_share, 
                                            target.PRS.R2, target.PCS_gain.R2, target.PCS_PRS.R2, target.heritability_gain.R2, target.heritability)

res.PCAS_simulate.external.to_plot$Scenario = droplevels(res.PCAS_simulate.external.to_plot$Scenario)
res.PCAS_simulate.external.to_plot_boxplots = do.call(rbind, lapply(levels(res.PCAS_simulate.external.to_plot$Scenario), function(scenari) {
  ele = filter(res.PCAS_simulate.external.to_plot, Scenario == scenari)
  data.frame(Scenario = scenari,
  Heritability_gain = median(ele$target.heritability_gain.R2),
  Heritability_gain.CI = diff(boxplot.stats(ele$target.heritability_gain.R2)$stats[c(1,3)])/2,
  PRS = median(ele$target.PRS.R2),
  PRS.CI = diff(boxplot.stats(ele$target.PRS.R2)$stats[c(1,3)])/2,
  PCS_gain = median(ele$target.PCS_gain.R2),
  PCS_gain.CI = diff(boxplot.stats(ele$target.PCS_gain.R2)$stats[c(1,3)])/2, stringsAsFactors = T)
}))

res.PCAS_simulate.external.to_plot = res.PCAS_simulate.external.to_plot_boxplots
res.PCAS_simulate.external.to_plot = cbind(pivot_longer(select(res.PCAS_simulate.external.to_plot, Scenario, PRS, PCS_gain, Heritability_gain), cols = c('PRS','PCS_gain', 'Heritability_gain'), names_to = 'model'),
                                           pivot_longer(select(res.PCAS_simulate.external.to_plot, PRS.CI, PCS_gain.CI, Heritability_gain.CI), cols = c('PRS.CI', 'PCS_gain.CI', 'Heritability_gain.CI'), values_to = 'CI'))
res.PCAS_simulate.external.to_plot$name = NULL
res.PCAS_simulate.external.to_plot$model = gsub('_', ' ', res.PCAS_simulate.external.to_plot$model)

res.PCAS_simulate.external.to_plot$label_pos = NA
res.PCAS_simulate.external.to_plot[res.PCAS_simulate.external.to_plot$model == 'PRS',]$label_pos = 
  res.PCAS_simulate.external.to_plot[res.PCAS_simulate.external.to_plot$model == 'PRS',]$value
res.PCAS_simulate.external.to_plot[res.PCAS_simulate.external.to_plot$model == 'PCS gain',]$label_pos = 
  res.PCAS_simulate.external.to_plot[res.PCAS_simulate.external.to_plot$model == 'PRS',]$value + 
  res.PCAS_simulate.external.to_plot[res.PCAS_simulate.external.to_plot$model == 'PCS gain',]$value
res.PCAS_simulate.external.to_plot[res.PCAS_simulate.external.to_plot$model == 'Heritability_gain',]$label_pos = 
  res.PCAS_simulate.external.to_plot[res.PCAS_simulate.external.to_plot$model == 'Heritability_gain',]$value +
  res.PCAS_simulate.external.to_plot[res.PCAS_simulate.external.to_plot$model == 'PRS',]$value + 
  res.PCAS_simulate.external.to_plot[res.PCAS_simulate.external.to_plot$model == 'PCS gain',]$value

gg.size.E =  gg.size + gg.size * 2/3
legend_title=''

p = ggplot(data = filter(res.PCAS_simulate.external.to_plot, model != 'Heritability gain'), aes(x=Scenario, y=value, fill=model)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=label_pos - CI, ymax = label_pos + CI), width=.2,
                position=position_dodge(.9), color = "#999999") +
  labs(x = "Scenario", y = "PVE", fill = legend_title) +
  scale_fill_manual(labels = list(expression(hat(Y)[paste('PGS+AS')]), 
                                  expression(hat(Y)[PGS])), values = c("#0072B2", "#56B4E9") ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size=gg.size.E),
        axis.title.x = element_blank());p
ggsave(plot = p, filename = paste0('results/sim.E.results'), dpi = gg.dpi, width = gg.width, height = gg.height, units = 'in')