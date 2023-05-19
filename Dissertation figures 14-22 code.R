###################################################################
####################  Arm-level CNA Figures  ######################
###################################################################

## Libraries
library(fst)
library(data.table)
library(ggplot2)
library(dndscv)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggpubr)
theme_set(theme_pubr())
library(ggtext)
library(ggrepel)
library(forcats)

## Parameters
outputs.folder <- 'output data/'

##Geneci Function
geneci = function(dndsout, gene_list = NULL, level = 0.95) {
  
  # Ensuring valid level value
  if (level > 1) {
    warning("Confidence level must be lower than 1, using 0.95 as default")
    level = 0.95
  }
  
  # N and L matrices
  N = dndsout$N
  L = dndsout$L
  if (length(N)==0) { stop(sprintf("Invalid input: the dndsout input object must be generated using outmats=T as an argument to dndscv.")) }
  if (nrow(dndsout$mle_submodel)!=195) { stop(sprintf("Invalid input: dndsout must be generated using the default trinucleotide substitution model in dndscv."))}
  
  # Restricting the analysis to an input list of genes
  if (!is.null(gene_list)) {
    g = as.vector(dndsout$genemuts$gene_name) # Genes in the input object
    nonex = gene_list[!(gene_list %in% g)] # Excluding genes from the input gene_list if they are not present in the input dndsout object
    if (length(nonex)>0) {
      warning(sprintf("The following input gene names are not in dndsout input object and will not be analysed: %s.", paste(nonex,collapse=", ")))
    }
    dndsout$annotmuts = dndsout$annotmuts[which(dndsout$annotmuts$gene %in% gene_list), ] # Restricting to genes of interest
    dndsout$genemuts = dndsout$genemuts[which(g %in% gene_list), ] # Restricting to genes of interest
    N = N[,,which(g %in% gene_list)] # Restricting to genes of interest
    L = L[,,which(g %in% gene_list)] # Restricting to genes of interest
  }
  gene_list = as.vector(dndsout$genemuts$gene_name)
  
  wnoneqspl = all(dndsout$sel_cv$wnon_cv==dndsout$sel_cv$wspl_cv) # Deciding on wnon==wspl based on the input object
  
  ## Subfunction: Analytical opt_t (aka tML) given fixed w values
  mle_tcvgivenw = function(n, theta, exp_neutral_cv, E) {
    shape = theta; scale = exp_neutral_cv/theta
    tml = (n+shape-1)/(1+E+(1/scale))
    if (shape<=1) { # i.e. when theta<=1
      tml = max(shape*scale,tml) # i.e. tml is bounded to the mean of the gamma (i.e. y[11]) when theta<=1
    }
    return(pmax(tml,1e-6))
  }
  
  ## Subfunction: Log-Likelihood of the model given fixed w values (requires finding MLEs for t and the free w values given the fixed w values)
  loglik_givenw = function(w,x,y,mutrates,theta,wtype,wnoneqspl,wmle) {
    
    # 1. tML given w
    exp_neutral_cv = y[11]
    exp_rel = y[7:10]/y[6]
    n = y[1] + sum(y[wtype+1])
    E = sum(exp_rel[wtype])*w
    tML = mle_tcvgivenw(n, theta, exp_neutral_cv, E)
    mrfold = max(1e-10, tML/y[6]) # Correction factor of "t" under the model
    
    # 2. Calculating the MLEs of the unconstrained w values
    if (!wnoneqspl) {
      wfree = y[2:5]/y[7:10]/mrfold; wfree[y[2:5]==0] = 0 # MLEs for w given tval
    } else {
      wmisfree = y[2]/y[7]/mrfold; wmisfree[y[2]==0] = 0
      wtruncfree = sum(y[3:4])/sum(y[8:9])/mrfold; wtruncfree[sum(y[3:4])==0] = 0
      wfree = c(wmisfree,wtruncfree,wtruncfree) # MLEs for w given tval
    }
    if(!all(wtype == 4)) wfree[wtype] <- w # Replacing free w values by fixed input values
    if( all(wtype == 4)) wfree = wfree * (w / ifelse(wmle[wtype] > 0, wmle[wtype], 1e-9 ) )
    
    # 2. loglik of the model under tML and w
    llpois = sum(dpois(x=x$n, lambda=x$l*mutrates*mrfold*t(array(c(1,wfree),dim=c(4,length(mutrates)))), log=T))
    llgamm = dgamma(x=tML, shape=theta, scale=exp_neutral_cv/theta, log=T)
    return(-(llpois+llgamm))
  }
  
  ## Subfunction: Working with vector inputs
  loglik_vec = function(wfixed,x,y,mutrates,theta,wtype,wnoneqspl,wmle) {
    sapply(wfixed, function(w) loglik_givenw(w,x,y,mutrates,theta,wtype,wnoneqspl,wmle))
  }
  
  
  ## Subfunction: iterative search for the CI95% boundaries for wvec
  iterative_search_ci95 = function(wtype,x,y,mutrates,theta,wmle,ml,grid_size=10,iter=10,wnoneqspl=T,wmax = 10000) {
    
    if (wmle[wtype][1]<wmax) {
      
      # Iteratively searching for the lower bound of the CI95% for "t"
      if (wmle[wtype][1]>0) {
        search_range = c(1e-9, wmle[wtype][1])
        for (it in 1:iter) {
          wvec = seq(search_range[1], search_range[2],length.out=grid_size)
          ll = -loglik_vec(wvec,x,y,mutrates,theta,wtype,wnoneqspl,wmle)
          lr = 2*(ml-ll) > qchisq(p=level,df=1)
          ind = max(which(wvec<=wmle[wtype][1] & lr))
          search_range = c(wvec[ind], wvec[ind+1])
        }
        w_low = wvec[ind]
      } else {
        w_low = 0
      }
      
      # Iteratively searching for the higher bound of the CI95% for "t"
      search_range = c(wmle[wtype][1], wmax)
      llhighbound = -loglik_vec(wmax,x,y,mutrates,theta,wtype,wnoneqspl,wmle)
      outofboundaries = !(2*(ml-llhighbound) > qchisq(p=level,df=1))
      if (!outofboundaries) {
        for (it in 1:iter) {
          wvec = seq(search_range[1], search_range[2],length.out=grid_size)
          ll = -loglik_vec(wvec,x,y,mutrates,theta,wtype,wnoneqspl,wmle)
          lr = 2*(ml-ll) > qchisq(p=level,df=1)
          ind = min(which(wvec>=wmle[wtype][1] & lr))
          search_range = c(wvec[ind-1], wvec[ind])
        }
        w_high = wvec[ind]
      } else {
        w_high = wmax
      }
      
    } else {
      wmle[wtype] = w_low = w_high = wmax
    }
    
    return(c(wmle[wtype][1],w_low,w_high))
  }
  
  
  ## Subfunction: calculate the MLEs and CI95% of each independent w value (unconstraining the other values)
  ci95cv_intt = function(x,y,mutrates,theta,grid_size=10,iter=10,wnoneqspl=T) {
    
    # MLE
    exp_neutral_cv = y[11]
    n = y[1]; E = 0 # Only synonymous mutations are considered
    tML = mle_tcvgivenw(n, theta, exp_neutral_cv, E)
    mrfold = max(1e-10, tML/y[6])
    if (!wnoneqspl) {
      wmle = y[2:4]/y[6:8]/mrfold; wmle[y[2:4]==0] = 0 # MLEs for w given tval
      wmle_all = y[5]/y[10]/mrfold; wmle_all[y[5]==0] = 0
      wmle <- c(wmle, wmle_all)
    } else {
      wmisfree = y[2]/y[7]/mrfold; wmisfree[y[2]==0] = 0
      wtruncfree = sum(y[3:4])/sum(y[8:9])/mrfold; wtruncfree[sum(y[3:4])==0] = 0
      wmle_all = y[5]/y[10]/mrfold; wmle_all[y[5]==0] = 0
      wmle = c(wmisfree,wtruncfree,wtruncfree,wmle_all) # MLEs for w given tval
    }
    llpois = sum(dpois(x=x$n, lambda=x$l*mutrates*mrfold*t(array(c(1,wmle[1:3]),dim=c(4,length(mutrates)))), log=T))
    llgamm = dgamma(x=tML, shape=theta, scale=y[11]/theta, log=T)
    ml = llpois+llgamm
    
    # Iteratively searching for the lower bound of the CI95% for "t"
    w_ci95 = array(NA,c(4,3))
    colnames(w_ci95) = c("mle","low","high")
    rownames(w_ci95) = c("mis","non","spl", 'all')
    if (!wnoneqspl) {
      for (h in 1:4) {
        w_ci95[h,] = iterative_search_ci95(wtype=h,x,y,mutrates,theta,wmle,ml,grid_size,iter,wnoneqspl)
      }
    } else {
      w_ci95[1,] = iterative_search_ci95(wtype=1,x,y,mutrates,theta,wmle,ml,grid_size,iter,wnoneqspl)
      w_ci95[2,] = iterative_search_ci95(wtype=c(2,3),x,y,mutrates,theta,wmle,ml,grid_size,iter,wnoneqspl)
      w_ci95[3,] = w_ci95[2,]
      w_ci95[4,] = iterative_search_ci95(wtype=4,x,y,mutrates,theta,wmle,ml,grid_size,iter,wnoneqspl)
    }
    return(w_ci95)
  }
  
  
  ## Calculating CI95% across all genes
  
  message("Calculating CI95 across all genes...")
  
  ci95 = array(NA, dim=c(length(gene_list),12))
  colnames(ci95) = c("mis_mle","non_mle","spl_mle","all_mle",
                     "mis_low","non_low","spl_low","all_low",
                     "mis_high","non_high","spl_high","all_high")
  theta = dndsout$nbreg$theta
  
  data("submod_192r_3w", package="dndscv")
  parmle =  setNames(dndsout$mle_submodel[,2], dndsout$mle_submodel[,1])
  mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
  
  for (j in 1:length(gene_list)) {
    geneind = which(dndsout$genemuts$gene_name==gene_list[j])
    y = as.numeric(dndsout$genemuts[geneind,-1])
    all_nsyn <- sum(y[2:4])
    expall_nsyn <- sum(y[6:8])
    y = c(y[1:4], all_nsyn, y[5:8], expall_nsyn, y[9])
    if (length(gene_list)==1) {
      x = list(n=N, l=L)
    } else {
      x = list(n=N[,,geneind], l=L[,,geneind])
    }
    ci95[j,] = c(ci95cv_intt(x,y,mutrates,theta,grid_size=10,iter=10,wnoneqspl=wnoneqspl))
    if (round(j/1000)==(j/1000)) { print(j/length(gene_list), digits=2) } # Progress
  }
  
  ci95df = cbind(gene=gene_list, as.data.frame(ci95))
  
  # Restricting columns if we forced wnon==wspl
  if (wnoneqspl==T) {
    ci95df = ci95df[,-c(4,8,12)]
    colnames(ci95df) = c("gene","mis_mle","tru_mle","all_mle","mis_low","tru_low","all_low","mis_high","tru_high","all_high")
  }
  
  return(ci95df)
  
}

## Inputs
arm_events_path <- '20230214_peace_chr_arm_events_cruk.tsv'
muttable_path <- "TRACERx421_data_code/20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_mutation_table.fst"
mut_driver_lung_path <- "TRACERx421_data_code/20221109_Tumour_evo_histories_DATA/gene_anno_df.fst"
clinical_data_path <- "TRACERx421_data_code/20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_all_tumour_df.rds"

# Load arm events
arm_events_all <- fread( arm_events_path )

# Load muttable
muttable <- fst::read_fst( muttable_path, as.data.table = TRUE )

# Load mut driver gene list
mut_drivers <- fst::read.fst( mut_driver_lung_path,as.data.table = T )
mut_drivers <- mut_drivers[ (is_lung_mut_driver), gene_name ]

mut_drivers_cd <- mut_drivers
mut_drivers_cd <- c(mut_drivers_cd, "CDKN2A.p14arf",   "CDKN2A.p16INK4a")
mut_drivers_cd <- mut_drivers_cd[ !mut_drivers_cd == 'CDKN2A' ]

# Load clinical data
clin <- as.data.table( readRDS( clinical_data_path) )

## Prepare Data:
clin = clin %>%
  rename(tumour_id = tumour_id_muttable_cruk)

clin_muttable = merge(clin, muttable)

LUAD_clin_arm = merge(clin, arm_events_all) %>%
  filter(histology_3 == "LUAD", arm_event_region != "neutral") %>%
  unite(arm_event, c(chr_arm, arm_event_region), sep = " ", remove = FALSE) %>%
  distinct(tumour_id, arm_event, arm_clonality) %>%
  add_count(tumour_id) %>%
  mutate(histology = "LUAD")

LUSC_clin_arm = merge(clin, arm_events_all) %>%
  filter(histology_3 == "LUSC", arm_event_region != "neutral") %>%
  unite(arm_event, c(chr_arm, arm_event_region), sep = " ", remove = FALSE) %>%
  distinct(tumour_id, arm_event, arm_clonality) %>%
  add_count(tumour_id) %>%
  mutate(histology = "LUSC")

### Figure 14: Box plot showing number of chr arm alterations per tumour
arm_events_tumour = rbind(LUAD_clin_arm, LUSC_clin_arm)

bp_armevents = arm_events_tumour %>%
  ggplot(aes(histology, n, fill = histology)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#FFFF99", "#CAB2D6"), name = "Histology") +
  stat_summary(fun = "mean", geom = "point", colour = "#E31A1C", size = 6, alpha = 0.8, show.legend = FALSE) +
  annotate("label", x = c(1.23,2.25), y = c(24, 26),
           label = c("mean = 16", "mean = 20"), size = 5)+
  annotate("segment", x = 1, xend = 1.23, y = 16.40, yend = 22.4, size = 1, alpha = 0.7, linetype = 2, lwd = 0.8) +
  annotate("segment", x = 2, xend = 2.25, y = 19.5, yend = 24.5, size = 1, alpha = 0.7, linetype = 2, lwd = 0.8) +
  labs(x = '\nHistology', y = 'Number of SCAAs per tumour\n', title = "Number of SCAAs per tumour in TRACERx421\n\n") +
  scale_y_continuous(limits = c(0, 45), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black", vjust = -1),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.position = "right",
        plot.title = element_text(hjust = 0.4, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"),
        aspect.ratio = 0.8)
print(bp_armevents)

ggsave("plot14.png", bp_armevents, width = 12, height = 10, dpi = 300)

### Figure 15: Pie charts showing proportion clonal/subclonal chr arm alterations in LUAD and LUSC

LUAD_pie_clonality <- data.frame(Timing = c("clonal", "subclonal"),
                                 Frequency = c(1199, 1494))
LUAD_pie_clonality = LUAD_pie_clonality %>%
  mutate(Frequency = Frequency/sum(Frequency)*100)

LUAD_pie_clonality$Timing = factor(LUAD_pie_clonality$Timing, levels = c("clonal", "subclonal"))

LUAD_pie_plot = LUAD_pie_clonality %>%
  ggplot(aes(x = "", y = Frequency, fill = Timing)) +
  geom_bar(width = 1, stat = "identity", linewidth = 0.5, colour = "#525252") +
  theme_classic() +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C"), name = 'Timing of SCAA', labels = c("clonal", "subclonal")) +
  annotate("text", x = c(1, 1), y = c(77, 29),
           label = c("44.52%", "55.48%"), colour = "white", size = 8) +
  ggtitle("Evolutionary timing of SCAAs in TRACERx421 (LUAD)\n\n") +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  theme(axis.title.x = element_blank(), axis.line.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank(), panel.grid = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines")) +
  coord_polar("y", start = 0)
print(LUAD_pie_plot)

ggsave("plot15a.png", LUAD_pie_plot, width = 12, height = 10, dpi = 300)

LUSC_pie_clonality <- data.frame(Timing = c("clonal", "subclonal"),
                                 Frequency = c(979, 1212))
LUSC_pie_clonality = LUSC_pie_clonality %>%
  mutate(Frequency = Frequency/sum(Frequency)*100)

LUSC_pie_plot = LUSC_pie_clonality %>%
  ggplot(aes(x = "", y = Frequency, fill = Timing)) +
  geom_bar(width = 1, stat = "identity", linewidth = 0.5, colour = "#525252") +
  theme_classic() +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C"), name = 'Timing of SCAA', labels = c("Clonal", "Subclonal")) +
  annotate("text", x = c(1, 1), y = c(77, 29),
           label = c("44.68%", "55.32%"), colour = "white", size = 8) +
  ggtitle("Evolutionary timing of SCAAs in TRACERx421 (LUSC)\n\n") +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  theme(axis.title.x = element_blank(), axis.line.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank(), panel.grid = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines")) +
  coord_polar("y", start = 0)
print(LUSC_pie_plot)

ggsave("plot15b.png", LUSC_pie_plot, width = 12, height = 10, dpi = 300)

### Figure 19: Proportion of chr arm alterations at cohort level
clin_arm = merge(clin, arm_events_all) %>%
  select(tumour_id, histology_3, chr_arm:arm_clonality)

clin_arm_filtered = clin_arm %>%
  filter(histology_3 != "Other", arm_event_region != "deep_loss", arm_event_region != "amp") %>%
  group_by(histology_3, arm_event_region) %>%
  tally() %>%
  mutate(percent = (n/sum(n))*100) %>%
  mutate(freq = sum(n))

armeventprop_labels = clin_arm_filtered %>%
  distinct(freq, .keep_all = TRUE)

round_prop <- function(x, digits = 2, max_iter = 1000) {
  r <- round(x, digits)
  target_sum <- round(sum(x), digits)
  iter <- 0
  while(sum(r) != target_sum && iter < max_iter) {
    if(sum(r) < target_sum) {
      idx <- which.max((x - r) * 10 ^ digits)
      r[idx] <- r[idx] + 10^-digits
    }
    if(sum(r) > target_sum) {
      idx <- which.max((r - x) * 10 ^ digits)
      r[idx] <- r[idx] - 10^-digits
    }
    r <- round(r, digits)
    iter <- iter + 1
  }
  if (iter == max_iter) {
    idx <- which.max(r)
    r[idx] <- r[idx] + (target_sum - sum(r))
    warning("Maximum number of iterations reached before convergence. The largest value was adjusted to make the sum equal to the target value.")
  }
  r
}

clin_arm_filtered <- clin_arm_filtered %>%
  group_by(histology_3) %>%
  mutate(percent = round_prop(percent)) %>%
  ungroup()


arm_event_prop_plot = clin_arm_filtered %>%
  ggplot(aes(x = histology_3, y = percent, fill = factor(arm_event_region, levels = c("gain", "loss", "LOH", "loss_LOH", "neutral")))) +
  geom_bar(stat = "identity", linewidth = 0.6, colour = "#525252", width = 0.7) +
  annotate("text", x = c(1, 1, 1, 1, 2, 2, 2, 2), y = c(97, 87, 80, 43, 97, 85, 73, 34),
           label = c("6.61%", "9.33%", "4.88%", "77.45%", "6.38%", "14.97%", "8.39%", "69.01%"), colour = "white", size = 8) +
  annotate("text", x = c(1.49, 2.49), y = c(93, 93),
           label = c("1.73%", "1.25%"), colour = "black", size = 6) +
  annotate("segment", x = 1.3, xend = 1.39, y = 92.5, yend = 92.5, size = 1, alpha = 1, linetype = 1, lwd = 0.8) +
  annotate("segment", x = 2.3, xend = 2.39, y = 93, yend = 93, size = 1, alpha = 1, linetype = 1, lwd = 0.8) +
  theme_classic() +
  scale_fill_manual(values = c("#FDBF6F", "#1F78B4", "#B2DF8A", "#33A02C", "#969696"), name = 'Chromosome arm status', labels = c("Gain", "Loss", "LOH", "Loss-LOH", "Neutral")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  labs(x = '\nHistology', y = 'Proportion of chromosome arms (%)\n', title = 'Percentage of chromosome arms affected by SCAAs in TRACERx421\n\n') +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black", vjust = -1),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"))
print(arm_event_prop_plot)

ggsave("plot19.png", arm_event_prop_plot, width = 12, height = 10, dpi = 300)

###Figure 16: Clonal arm-level SCNA burden

clonalarmevents_tumour = clin_arm_muttable %>%
  filter(arm_clonality == "clonal") %>%
  distinct(tumour_id, arm_event) %>%
  add_count(tumour_id) %>%
  select(tumour_id, n) %>%
  distinct(tumour_id, .keep_all = TRUE) %>%
  rename(clonal_armevents = n)

clonalarm_muttable = merge(clonalarmevents_tumour, muttable)

clonaltumours = clin_arm_muttable %>%
  filter(arm_clonality == "clonal") %>%
  distinct(tumour_id)

no_clonalarm = clonalarm_muttable %>%
  filter(PyCloneClonal_SC_TB_corrected == "S", clonal_armevents >= 1 & clonal_armevents <= 2) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

low_clonalarm = clonalarm_muttable %>%
  filter(PyCloneClonal_SC_TB_corrected == "S", clonal_armevents >= 1 & clonal_armevents <= 2) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

mid_clonalarm = clonalarm_muttable %>%
  filter(PyCloneClonal_SC_TB_corrected == "S", clonal_armevents >= 3 & clonal_armevents <= 6) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

high_clonalarm = clonalarm_muttable %>%
  filter(PyCloneClonal_SC_TB_corrected == "S", clonal_armevents >= 7 & clonal_armevents <= 9) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

vhigh_clonalarm = clonalarm_muttable %>%
  filter(PyCloneClonal_SC_TB_corrected == "S", clonal_armevents >= 10) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

# Run dNdS
dndsout_low_clonalarm = dndscv(low_clonalarm, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_low_clonalarm = dndsout_low_clonalarm$globaldnds
globaldnds_low_clonalarm$group = "low"

dndsout_mid_clonalarm = dndscv(mid_clonalarm, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_mid_clonalarm = dndsout_mid_clonalarm$globaldnds
globaldnds_mid_clonalarm$group = "medium"

dndsout_high_clonalarm = dndscv(high_clonalarm, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_high_clonalarm = dndsout_high_clonalarm$globaldnds
globaldnds_high_clonalarm$group = "high"

dndsout_vhigh_clonalarm = dndscv(vhigh_clonalarm, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_vhigh_clonalarm = dndsout_vhigh_clonalarm$globaldnds
globaldnds_vhigh_clonalarm$group = "very high"

### Figure 16: Clonal arm burden
clonal_arm_tumours = clonalarm_muttable %>%
  filter(PyCloneClonal_SC_TB_corrected == "S") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

noclonal_arm_tumours = arm_muttable %>%
  filter(PyCloneClonal_SC_TB_corrected == "S") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

noclonal_arm_tumours = noclonal_arm_tumours[!(noclonal_arm_tumours$sampleID %in% clin_clonalarmevents$tumour_id),]

dndsout_clonal_arm_tumours = dndscv(clonal_arm_tumours, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_clonal_arm_tumours = dndsout_clonal_arm_tumours$globaldnds
globaldnds_clonal_arm_tumours$armgroup = "Clonal arm-level CNA"

dndsout_noclonal_arm_tumours = dndscv(noclonal_arm_tumours, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_noclonal_arm_tumours = dndsout_noclonal_arm_tumours$globaldnds
globaldnds_noclonal_arm_tumours$armgroup = "No clonal arm-level CNA"
globaldnds_noclonal_arm_tumours$group = "none"

sel_cv_clonal_arm_tumours = dndsout_clonal_arm_tumours$sel_cv %>%
  mutate(armgroup = "Clonal arm-level CNA")

sel_cv_noclonal_arm_tumours = dndsout_noclonal_arm_tumours$sel_cv %>%
  mutate(armgroup = "No Clonal arm-level CNA") %>%
  mutate(group = "none")

signif_genes_clonal_arm_tumours = rbind(sel_cv_clonal_arm_tumours, sel_cv_noclonal_arm_tumours) %>%
  filter(qglobal_cv < 0.05) %>%
  select(gene_name, qglobal_cv, group)

geneci_clonal_arm_tumours = geneci(dndsout_clonal_arm_tumours, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(armgroup = "Clonal arm-level CNA")

geneci_noclonal_arm_tumours = geneci(dndsout_noclonal_arm_tumours, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(armgroup = "No clonal arm-level CNA") %>%
  mutate(group = "none")

globaldnds_noclonal_arm_tumours_filtered = globaldnds_noclonal_arm_tumours %>%
  select(name, mle, cilow, cihigh, group)

global_clonalarm = rbind(globaldnds_noclonal_arm_tumours_filtered, globaldnds_low_clonalarm, globaldnds_mid_clonalarm, globaldnds_high_clonalarm, globaldnds_vhigh_clonalarm)

global_clonalarm_plot = global_clonalarm %>%
  filter(name == "wall") %>%
  ggplot(aes(x = factor(group, levels = c("none", "low", "medium", "high", "very high")), mle, color = group)) +
  geom_pointrange(aes(ymin = cilow, ymax = cihigh))+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_color_manual(values = c("#969696", "#1F78B4", "#A6CEE3", "#FDBF6F", "#E31A1C"), name = 'Category (number of\nclonal SCAAs)', breaks = c("none", "low", "medium", "high", "very high"), 
                     labels = c("none (0)", "low (1-2)", "medium (3-6)", "high (7-9)","very high (10+)")) +
  labs(x = '\n\nNumber of clonal SCAAs category\n', y = 'dN/dS subclonal driver genes\n', title = 'Cohort-level selection (dN/dS) of subclonal driver genes\n according to number of clonal SCAAs\n') +
  guides(colour = guide_legend(byrow = TRUE)) +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black", vjust = -1),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 15),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"))
print(global_clonalarm_plot)

ggsave("plot16.png", global_clonalarm_plot, width = 12, height = 10, dpi = 300)

### Appendix Figure 2: Clonal SCAA distribution

clin_clonalarmevents = merge(clonalarmevents_tumour, clin) %>%
  select(tumour_id, clonal_armevents, histology_3)

bp_clonal_armevents = ggplot(clin_clonalarmevents, aes(x = "", y = clonal_armevents)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, fill = "#CCCCCC") +
  geom_jitter(aes(fill = histology_3), size = 2, alpha = 0.7, shape = 21) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0)) +
  scale_fill_manual(values = c("#FDBF6F", "#6A3D9A")) +
  labs(x = "LUAD and LUSC tumours\n (n = 326)", y = "Number of clonal SCAAs per tumour\n", fill = "Histology", title = 'Distribution of clonal SCAAs in TRACERx421\n') +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black", vjust = 1),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18, vjust = 4),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"),
        aspect.ratio = 1)
print(bp_clonal_armevents)

ggsave("plota2.png", bp_clonal_armevents, width = 12, height = 10, dpi = 300)


### Figure 17: gene-level dNdS arm event
clonalarmevent_geneci_plot = filtered_geneci_clonal_arm_tumours %>%
  ggplot(aes(gene, all_mle)) +
  geom_pointrange(aes(ymin = all_low, ymax = all_high, colour = armgroup), position = position_dodge(width = 0.2))+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  scale_color_manual(values = c("#969696", "#FF7F00"), name = 'SCAA status', breaks = c("No clonal arm-level CNA", "Clonal arm-level CNA"), labels = c("No clonal SCAA", "Clonal SCAA")) +
  scale_y_log10() +
  labs(x = 'Gene name\n', y = '\ndN/dS subclonal driver genes\n', title = 'Gene-level selection (dN/dS) of subclonal driver genes\n in tumours with and without clonal SCAAs\n\n') +
  theme(axis.text.y = element_text(colour = c(rep("#FF7F00", times = 3), rep("#6A3D9A" , times = 5), rep("#525252", times = 6)), size = 12),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.4, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"),
        aspect.ratio = 1) +
  coord_flip()
print(clonalarmevent_geneci_plot)

ggsave("plot17.png", clonalarmevent_geneci_plot, width = 12, height = 10, dpi = 300)


### Figure 20a: Clonal arm heatmap

arm_heatmap_genelist = rbind(filtered_geneci_clonal_arm_tumours, filtered_geneci_9p_loss_LOH) %>%
  arrange(gene) %>%
  distinct(gene)

all_arm_geneci = rbind(filtered_geneci_noclonal_arm_tumours, geneci_clonal_arm_tumours, geneci_9p_loss_LOH, geneci_no_9p_loss_LOH)

arm_heatmapdata = all_arm_geneci %>%
  filter(all_arm_geneci$gene %in% arm_heatmap_genelist$gene) %>%
  select(gene, all_mle, group) %>%
  mutate(scaled_mle = sqrt(all_mle))

arm_heatmapdata$group = factor(arm_heatmapdata$group, levels = c("No clonal arm-level CNA", "Clonal arm-level CNA", "No 9p loss LOH", "9p loss LOH"))

arm_heatmap = arm_heatmapdata %>%
  ggplot(aes(x = group, y = gene, fill = scaled_mle)) +
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "Greens", direction = +1) +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size = 14, face = "bold", hjust = -24)) +
  labs(x = "\nAbsence or presence of clonal arm event", y = "Lung driver gene", fill = "dN/dS") +
  ggtitle("Selection of lung driver genes in tumours with clonal arm event vs without")
print(arm_heatmap)

clonalarm_geneci = rbind(filtered_geneci_noclonal_arm_tumours, geneci_clonal_arm_tumours)

clonalarm_heatmapdata = clonalarm_geneci %>%
  filter(gene != "CDKN2A.p14arf") %>%
  select(gene, all_mle, armgroup) %>%
  mutate(scaled_mle = sqrt(all_mle))

clonalarm_heatmapdata$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", clonalarm_heatmapdata$gene)

clonalarm_heatmapdata = clonalarm_heatmapdata %>%
  filter(clonalarm_heatmapdata$gene %in% arm_heatmap_genelist$gene)

clonalarm_heatmapdata$gene = factor(clonalarm_heatmapdata$gene, levels = c("SFTPB", "SMAD4", "CDKN2A",
                                                                           "KMT2D", "ARID1A", "BCLAF1", "FAM78B", "DUSP22",
                                                                           "PIK3CA", "STX2", "KRAS", "PTEN", "B2M", "TP53"))

clonalarm_heatmapdata$group = factor(clonalarm_heatmapdata$armgroup, levels = c("No clonal arm-level CNA", "Clonal arm-level CNA"))

clonalarm_heatmap = clonalarm_heatmapdata %>%
  ggplot(aes(x = group, y = gene, fill = log(all_mle))) +
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "Oranges", direction = +1) +
  scale_x_discrete(labels = c("No clonal SCAAs", "Clonal SCAAs")) +
  theme_classic() +
  labs(x = "\nSCAA status", y = "Gene name\n", fill = "log(dN/dS)", 
       title = "Strength of dN/dS signal in subclonal driver genes\n in tumours with and without clonal SCAAs\n") +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "left",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.4, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))
print(clonalarm_heatmap)

clonalarm_heatmap_withnumbers = clonalarm_heatmap +
  geom_label(aes(label = sprintf("%.2f", all_mle)),
             colour = "black",
             fill = "white", alpha = 0.8,
             label.r = unit(0, "npc"),
             size = 7)

print(clonalarm_heatmap_withnumbers)

ggsave("plot20a.png", clonalarm_heatmap_withnumbers, width = 12, height = 10, dpi = 300)

### Figure 20b: Clonal arm freq

clonalarm_tumourids = unique(clonal_arm_tumours$sampleID)

armdriverfreq_genelist = unique(arm_heatmap_genelist$gene)

clonalarm_driverfreq_muttable = arm_muttable %>%
  filter(arm_muttable$Hugo_Symbol %in% armdriverfreq_genelist | Hugo_Symbol == "CDKN2A", arm_muttable$tumour_id %in% clonalarm_tumourids, (combTiming_SC == "subclonal")) %>%
  distinct(Hugo_Symbol, tumour_id, .keep_all = TRUE) %>%
  add_count(Hugo_Symbol) %>%
  mutate(arm_group = "Clonal arm-level CNA event") %>%
  mutate(total_tumours = 145) %>%
  mutate(percent = (n/145)*100)

noclonalarm_tumourids = unique(noclonal_arm_tumours$sampleID)  

noclonalarm_driverfreq_muttable = arm_muttable %>%
  filter(arm_muttable$Hugo_Symbol %in% armdriverfreq_genelist | Hugo_Symbol == "CDKN2A", arm_muttable$tumour_id %in% noclonalarm_tumourids, (combTiming_SC == "subclonal")) %>%
  distinct(Hugo_Symbol, tumour_id, .keep_all = TRUE) %>%
  add_count(Hugo_Symbol) %>%
  mutate(arm_group = "No clonal arm-level CNA event") %>%
  mutate(total_tumours = 34) %>%
  mutate(percent = (n/34)*100)

combinedclonalarm_driverfreq_muttable = rbind(clonalarm_driverfreq_muttable, noclonalarm_driverfreq_muttable) %>%
  select(Hugo_Symbol, percent, arm_group)

clonalarmrows = data.frame(Hugo_Symbol = c("ZBTB24", "SFTPB", "CDKN2A"),
                           percent = c(0, 0, 0),
                           arm_group = c("No clonal arm-level CNA event", "No clonal arm-level CNA event", "No clonal arm-level CNA event"))

combinedclonalarm_driverfreq_muttable = rbind(combinedclonalarm_driverfreq_muttable, clonalarmrows)

combinedclonalarm_driverfreq_muttable$Hugo_Symbol = factor(combinedclonalarm_driverfreq_muttable$Hugo_Symbol, levels = c("SFTPB", "SMAD4", "CDKN2A",
                                                                                                                         "KMT2D", "ARID1A", "BCLAF1", "FAM78B", "DUSP22",
                                                                                                                         "PIK3CA", "STX2", "KRAS", "PTEN", "B2M", "TP53"))

combinedclonalarm_driverfreq_muttable = combinedclonalarm_driverfreq_muttable %>%
  drop_na(Hugo_Symbol)

clolnalarm_matrix_driverfreq = combinedclonalarm_driverfreq_muttable %>%
  distinct(arm_group, Hugo_Symbol, percent) %>%
  arrange(Hugo_Symbol)

clonalarm_driverfreq_plot = combinedclonalarm_driverfreq_muttable %>%
  ggplot(aes(Hugo_Symbol, percent, fill = factor(arm_group, levels = c("Clonal arm-level CNA event", "No clonal arm-level CNA event")))) +
  geom_col(position = position_dodge(), width = 0.7, colour = "#525252", linewidth = 0.2) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 40), expand = c(0, 0)) +
  scale_fill_manual(values = c("#CCCCCC", "#FF7F00"), breaks = c("No clonal arm-level CNA event", "Clonal arm-level CNA event"), labels = c("No clonal SCAA", "Clonal SCAA")) +
  guides(fill = guide_legend(byrow = TRUE, override.aes = list(color = NA))) +
  labs(x = "Gene name\n", y = "\nProportion of tumours (%)", fill = "SCAA status",
       title = "Frequency of subclonal driver genes in tumours with and without clonal SCAAs\n") +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.4, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))+
  coord_flip()

ggsave("plot20b.png", clonalarm_driverfreq_plot, width = 12, height = 10, dpi = 300)

### Figure 20: Global dN/dS of most frequent arm events
top_clonalarm = clin_arm_muttable %>%
  filter(arm_clonality == "clonal") %>%
  add_count(arm_event) %>%
  arrange(desc(n)) %>%
  filter(n > 32.6)

top_clonalarm_muttable = merge(top_clonalarm, muttable) %>%
  filter(PyCloneClonal_SC_TB_corrected == "S")

# Create groups
arm_1q_gain = top_clonalarm_muttable %>%
  filter(arm_event == "1q gain") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_9q_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "9q LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_5p_gain = top_clonalarm_muttable %>%
  filter(arm_event == "5p gain") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_3p_loss_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "3p loss_LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_13q_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "13q LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_19p_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "19p LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_3p_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "3p LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_17p_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "17p LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_9p_loss_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "9p loss_LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_17q_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "17q LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_9p_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "9p LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_7p_gain = top_clonalarm_muttable %>%
  filter(arm_event == "7p gain") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_5q_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "5q LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_5q_loss_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "5q loss_LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_18q_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "18q LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_13q_loss_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "13q loss_LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_17p_loss_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "17p loss_LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

arm_12q_LOH = top_clonalarm_muttable %>%
  filter(arm_event == "12q LOH") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

# Run dN/dS on most frequent arm events

#LUAD
dndsout_arm_1q_gain = dndscv(arm_1q_gain, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_1q_gain = dndsout_arm_1q_gain$globaldnds
globaldnds_arm_1q_gain$group = "1q gain"
globaldnds_arm_1q_gain$type = "gain"

dndsout_arm_9q_LOH = dndscv(arm_9q_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_9q_LOH = dndsout_arm_9q_LOH$globaldnds
globaldnds_arm_9q_LOH$group = "9q LOH"
globaldnds_arm_9q_LOH$type = "LOH"

dndsout_arm_5p_gain = dndscv(arm_5p_gain, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_5p_gain = dndsout_arm_5p_gain$globaldnds
globaldnds_arm_5p_gain$group = "5p gain"
globaldnds_arm_5p_gain$type = "gain"

dndsout_arm_3p_loss_LOH = dndscv(arm_3p_loss_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_3p_loss_LOH = dndsout_arm_3p_loss_LOH$globaldnds
globaldnds_arm_3p_loss_LOH$group = "3p loss LOH"
globaldnds_arm_3p_loss_LOH$type = "loss LOH"

dndsout_arm_13q_LOH = dndscv(arm_13q_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_13q_LOH = dndsout_arm_13q_LOH$globaldnds
globaldnds_arm_13q_LOH$group = "13q LOH"
globaldnds_arm_13q_LOH$type = "LOH"

dndsout_arm_19p_LOH = dndscv(arm_19p_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_19p_LOH = dndsout_arm_19p_LOH$globaldnds
globaldnds_arm_19p_LOH$group = "19p LOH"
globaldnds_arm_19p_LOH$type = "LOH"

dndsout_arm_3p_LOH = dndscv(arm_3p_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_3p_LOH = dndsout_arm_3p_LOH$globaldnds
globaldnds_arm_3p_LOH$group = "3p LOH"
globaldnds_arm_3p_LOH$type = "LOH"

dndsout_arm_17p_LOH = dndscv(arm_17p_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_17p_LOH = dndsout_arm_17p_LOH$globaldnds
globaldnds_arm_17p_LOH$group = "17p LOH"
globaldnds_arm_17p_LOH$type = "LOH"

dndsout_arm_9p_loss_LOH = dndscv(arm_9p_loss_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_9p_loss_LOH = dndsout_arm_9p_loss_LOH$globaldnds
globaldnds_arm_9p_loss_LOH$group = "9p loss LOH"
globaldnds_arm_9p_loss_LOH$type = "loss LOH"

dndsout_arm_17q_LOH = dndscv(arm_17q_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_17q_LOH = dndsout_arm_17q_LOH$globaldnds
globaldnds_arm_17q_LOH$group = "17q LOH"
globaldnds_arm_17q_LOH$type = "LOH"

dndsout_arm_9p_LOH = dndscv(arm_9p_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_9p_LOH = dndsout_arm_9p_LOH$globaldnds
globaldnds_arm_9p_LOH$group = "9p LOH"
globaldnds_arm_9p_LOH$type = "LOH"

dndsout_arm_7p_gain = dndscv(arm_7p_gain, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_7p_gain = dndsout_arm_7p_gain$globaldnds
globaldnds_arm_7p_gain$group = "7p gain"
globaldnds_arm_7p_gain$type = "gain"

dndsout_arm_5q_LOH = dndscv(arm_5q_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_5q_LOH = dndsout_arm_5q_LOH$globaldnds
globaldnds_arm_5q_LOH$group = "5q LOH"
globaldnds_arm_5q_LOH$type = "LOH"

dndsout_arm_5q_loss_LOH = dndscv(arm_5q_loss_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_5q_loss_LOH = dndsout_arm_5q_loss_LOH$globaldnds
globaldnds_arm_5q_loss_LOH$group = "5q loss LOH"
globaldnds_arm_5q_loss_LOH$type = "loss LOH"

dndsout_arm_18q_LOH = dndscv(arm_18q_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_18q_LOH = dndsout_arm_18q_LOH$globaldnds
globaldnds_arm_18q_LOH$group = "18q LOH"
globaldnds_arm_18q_LOH$type = "LOH"

dndsout_arm_13q_loss_LOH = dndscv(arm_13q_loss_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_13q_loss_LOH = dndsout_arm_13q_loss_LOH$globaldnds
globaldnds_arm_13q_loss_LOH$group = "13q loss LOH"
globaldnds_arm_13q_loss_LOH$type = "loss LOH"

dndsout_arm_17p_loss_LOH = dndscv(arm_17p_loss_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_17p_loss_LOH = dndsout_arm_17p_loss_LOH$globaldnds
globaldnds_arm_17p_loss_LOH$group = "17p loss LOH"
globaldnds_arm_17p_loss_LOH$type = "loss LOH"

dndsout_arm_12q_LOH = dndscv(arm_12q_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_arm_12q_LOH = dndsout_arm_12q_LOH$globaldnds
globaldnds_arm_12q_LOH$group = "12q LOH"
globaldnds_arm_12q_LOH$type = "LOH"

global_clonalarmevents = rbind(globaldnds_arm_1q_gain, globaldnds_arm_9q_LOH, globaldnds_arm_5p_gain, 
                               globaldnds_arm_3p_loss_LOH, globaldnds_arm_13q_LOH, globaldnds_arm_19p_LOH,
                               globaldnds_arm_3p_LOH, globaldnds_arm_17p_LOH, globaldnds_arm_9p_loss_LOH,
                               globaldnds_arm_17q_LOH, globaldnds_arm_9p_LOH, globaldnds_arm_7p_gain,
                               globaldnds_arm_5q_LOH, globaldnds_arm_5q_loss_LOH, globaldnds_arm_18q_LOH,
                               globaldnds_arm_13q_loss_LOH, globaldnds_arm_17p_loss_LOH, globaldnds_arm_12q_LOH)

global_clonalarmevents$group = factor(global_clonalarmevents$group, levels = c("1q gain", "9q LOH", "5p gain", "3p loss LOH", "13q LOH", "19p LOH", 
                                                                               "3p LOH", "17p LOH", "9p loss LOH", "9p LOH", "17q LOH", "7p gain", 
                                                                               "5q loss LOH", "5q LOH", "18q LOH", "13q loss LOH", "17p loss LOH", "12q LOH"))

global_clonalarmevents_plot = global_clonalarmevents %>%
  filter(name == "wall") %>%
  ggplot(aes(group, mle)) +
  geom_pointrange(aes(ymin = cilow, ymax = cihigh, color = type))+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0)) +
  scale_color_manual(values = c("#FDBF6F", "#B2DF8A", "#33A02C"), name = 'Type of SCAA', labels = c("Gain", "LOH", "Loss-LOH")) +
  guides(colour = guide_legend(byrow = TRUE)) +
  labs(x = '\nType of clonal SCAA event', y = 'dN/dS subclonal driver genes\n') +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(angle = 90, size = 12, colour = "black", vjust = 0.5, hjust = 0.8),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 15),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 13),
        legend.position = "right",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 3, 1, 1), "lines"))
print(global_clonalarmevents_plot)

clonalarmevents_freq_plot = top_clonalarm %>%
  distinct(arm_event, .keep_all = TRUE) %>%
  ggplot(aes(x = fct_rev(fct_reorder(arm_event, n)), y = n)) +
  geom_col(aes(fill = arm_event_region), linewidth = 0.4, colour = "#525252", width = 0.8) +
  theme_classic() +
  scale_fill_manual(values = c("#FDBF6F", "#B2DF8A", "#33A02C"), name = 'Type of SCAA', labels = c("Gain", "LOH", "Loss & LOH")) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  labs(x = '\nClonal SCAA event', y = 'Number of tumours\n') +
  guides(fill = guide_legend(override.aes = list(color = NA), byrow = TRUE)) +
  geom_text(aes(label = n), vjust = -0.5, size = 5) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.3, size = 20, face = "bold"),
        plot.margin = unit(c(1, 3, 1, 1), "lines"))
print(clonalarmevents_freq_plot)

ggarange_clonalarmevents_plot = ggarrange(clonalarmevents_freq_plot, global_clonalarmevents_plot, ncol = 1, nrow = 2, 
                                          common.legend = TRUE, legend = "right", align = "v")
print(ggarange_clonalarmevents_plot)

ggsave("plot20.png", ggarange_clonalarmevents_plot, width = 15, height = 10, dpi = 300)

### Figure 21:
distinct_armevents = filtered_arm_events %>%
  unite(arm_event, c(chr_arm, arm_event_region), sep = " ", remove = FALSE) %>%
  distinct(tumour_id, arm_event, arm_clonality)

arm_muttable = merge(distinct_armevents, clin_muttable)

arm_9p_loss_LOH_list = arm_9p_loss_LOH %>%
  distinct(sampleID) %>%
  rename(tumour_id = sampleID)

no_9p_loss_LOH = arm_muttable %>%
  filter(PyCloneClonal_SC_TB_corrected == "S") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

no_9p_loss_LOH = no_9p_loss_LOH[!(no_9p_loss_LOH$sampleID %in% arm_9p_loss_LOH_list$tumour_id),]

dndsout_no_9p_loss_LOH = dndscv(no_9p_loss_LOH, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_no_9p_loss_LOH = dndsout_no_9p_loss_LOH$globaldnds
globaldnds_no_9p_loss_LOH$group = "No 9p loss LOH"

sel_cv_no_9p_loss_LOH = dndsout_no_9p_loss_LOH$sel_cv %>%
  mutate(group = "No 9p loss LOH")

sel_cv_9p_loss_LOH = dndsout_arm_9p_loss_LOH$sel_cv %>%
  mutate(group = "9p loss LOH")

signif_genes_9p_loss_LOH = rbind(sel_cv_no_9p_loss_LOH, sel_cv_9p_loss_LOH) %>%
  filter(qglobal_cv < 0.05) %>%
  select(gene_name, qglobal_cv, group)

# Geneci WGD groups
geneci_no_9p_loss_LOH = geneci(dndsout_no_9p_loss_LOH, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(group = "No 9p loss LOH")

geneci_9p_loss_LOH = geneci(dndsout_arm_9p_loss_LOH, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(group = "9p loss LOH")

# Combine geneci WGD groups
combined_geneci_9p_loss_LOH = rbind(geneci_no_9p_loss_LOH, geneci_9p_loss_LOH) %>%
  arrange(gene)

filtered_geneci_9p_loss_LOH = combined_geneci_9p_loss_LOH %>%
  filter(combined_geneci_9p_loss_LOH$gene %in% signif_genes_9p_loss_LOH$gene_name, gene != "CDKN2A.p14arf")

filtered_geneci_9p_loss_LOH$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", filtered_geneci_9p_loss_LOH$gene)

filtered_geneci_9p_loss_LOH = filtered_geneci_9p_loss_LOH %>%
  group_by(gene) %>%
  mutate(odds_ratio = (all_mle[ group == 'No 9p loss LOH'] / all_mle[ group == '9p loss LOH'])) %>%
  arrange(gene) %>%
  mutate(armfavoured = case_when(odds_ratio > 2 ~ "No clonal 9p loss LOH favoured",
                                 odds_ratio <= 2 & odds_ratio >= 0.5 ~ "No clonal 9p loss LOH and clonal 9p loss LOH",
                                 odds_ratio < 0.5 ~ "Clonal 9p loss LOH favoured"))

qglobal_signif_genes_9plossLOH = signif_genes_9p_loss_LOH %>%
  select(gene_name, qglobal_cv) %>%
  rename(gene = gene_name) %>%
  filter(gene != "CDKN2A.p14arf")

qglobal_signif_genes_9plossLOH$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", qglobal_signif_genes_9plossLOH$gene)

qglobal_filtered_geneci_9plossLOH = merge(filtered_geneci_9p_loss_LOH, qglobal_signif_genes_9plossLOH)

qglobal_filtered_geneci_9plossLOH = qglobal_filtered_geneci_9plossLOH %>%
  group_by(armfavoured, qglobal_cv) %>%
  arrange(armfavoured, qglobal_cv)

qglobal_filtered_geneci_9plossLOH$armfavoured = factor(qglobal_filtered_geneci_9plossLOH$armfavoured, levels = c("No clonal 9p loss LOH favoured", "No clonal 9p loss LOH and clonal 9p loss LOH", "Clonal 9p loss LOH favoured"))

filtered_geneci_9p_loss_LOH$all_high = ifelse(abs(filtered_geneci_9p_loss_LOH$all_mle) == 0, filtered_geneci_9p_loss_LOH$all_high == 0, filtered_geneci_9p_loss_LOH$all_high)

filtered_geneci_9p_loss_LOH$gene = factor(filtered_geneci_9p_loss_LOH$gene, levels = c("CDKN2A", "ZBTB24", "PIK3CA", "B2M", "DUSP22", 
                                                                                       "SMAD4", "ARID1A", "TP53", 
                                                                                       "ARHGAP35", "BCLAF1", "KRAS", "PTEN"))

# Combine geneci WGD groups
filtered_geneci_noclonal_arm_tumours = geneci_noclonal_arm_tumours %>%
  select(gene, mis_mle, tru_mle, all_mle, mis_low, tru_low, all_low, mis_high, tru_high, all_high, armgroup)

combined_geneci_clonal_arm_tumorus = rbind(geneci_clonal_arm_tumours, filtered_geneci_noclonal_arm_tumours) %>%
  arrange(gene)

filtered_geneci_clonal_arm_tumours = combined_geneci_clonal_arm_tumorus %>%
  filter(combined_geneci_clonal_arm_tumorus$gene %in% signif_genes_clonal_arm_tumours$gene_name, gene != "CDKN2A.p14arf")

filtered_geneci_clonal_arm_tumours$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", filtered_geneci_clonal_arm_tumours$gene)

filtered_geneci_clonal_arm_tumours$all_high = ifelse(abs(filtered_geneci_clonal_arm_tumours$all_mle) == 0, filtered_geneci_clonal_arm_tumours$all_high == 0, filtered_geneci_clonal_arm_tumours$all_high)

filtered_geneci_clonal_arm_tumours = filtered_geneci_clonal_arm_tumours %>%
  group_by(gene) %>%
  mutate(odds_ratio = (all_mle[ armgroup == 'No clonal arm-level CNA'] / all_mle[ armgroup == 'Clonal arm-level CNA'])) %>%
  arrange(gene) %>%
  mutate(armfavoured = case_when(odds_ratio > 2 ~ "No clonal SCAA favoured",
                                 odds_ratio <= 2 & odds_ratio >= 0.5 ~ "No clonal SCAA and clonal SCAA",
                                 odds_ratio < 0.5 ~ "Clonal SCAA favoured"))

qglobal_signif_genes_clonalarms = signif_genes_clonal_arm_tumours %>%
  select(gene_name, qglobal_cv) %>%
  rename(gene = gene_name) %>%
  filter(gene != "CDKN2A.p14arf")

qglobal_signif_genes_clonalarms$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", qglobal_signif_genes_clonalarms$gene)

qglobal_filtered_geneci_clonalarms = merge(filtered_geneci_clonal_arm_tumours, qglobal_signif_genes_clonalarms)

qglobal_filtered_geneci_clonalarms = qglobal_filtered_geneci_clonalarms %>%
  group_by(armfavoured, qglobal_cv) %>%
  arrange(armfavoured, qglobal_cv)

qglobal_filtered_geneci_clonalarms$armfavoured = factor(qglobal_filtered_geneci_clonalarms$armfavoured, levels = c("No clonal SCAA favoured", "No clonal SCAA and clonal SCAA", "Clonal SCAA favoured"))

filtered_geneci_clonal_arm_tumours$all_high = ifelse(abs(filtered_geneci_clonal_arm_tumours$all_mle) == 0, filtered_geneci_clonal_arm_tumours$all_high == 0, filtered_geneci_clonal_arm_tumours$all_high)

filtered_geneci_clonal_arm_tumours$gene = factor(filtered_geneci_clonal_arm_tumours$gene, levels = c("SFTPB", "SMAD4", "CDKN2A",
                                                                                                     "KMT2D", "ARID1A", "BCLAF1", "FAM78B", "DUSP22",
                                                                                                     "PIK3CA", "STX2", "KRAS", "PTEN", "B2M", "TP53"))
arm_9p_loss_LOH_geneci_plot = filtered_geneci_9p_loss_LOH %>%
  ggplot(aes(gene, all_mle)) +
  geom_pointrange(aes(ymin = all_low, ymax = all_high, colour = group), position = position_dodge(width = 0.2))+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  scale_color_manual(values = c("#969696", "#33A02C"), name = 'SCAA status', breaks = c("No 9p loss LOH", "9p loss LOH"), labels = c("No clonal 9p loss & LOH", "Clonal 9p loss & LOH")) +
  scale_y_log10() +
  labs(x = 'Gene name\n', y = '\ndN/dS subclonal driver genes\n', title = 'Gene-level selection (dN/dS) of subclonal driver genes\n in tumours with and without clonal 9p loss & LOH\n\n') +
  theme(axis.text.y = element_text(colour = c(rep("#33A02C", times = 5), rep("#6A3D9A", times = 3), rep( "#525252", times = 4)), 
                                   size = 12),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.4, size = 20, face = "bold"),
        plot.margin = unit(c(1, 3, 3, 1), "lines"),
        aspect.ratio = 1) +
  coord_flip()
print(arm_9p_loss_LOH_geneci_plot)

ggsave("plot21.png", arm_9p_loss_LOH_geneci_plot, width = 12, height = 10, dpi = 300)

### Figure 22a: Heatmap

lossLOH9p_geneci = rbind(geneci_9p_loss_LOH, geneci_no_9p_loss_LOH)

lossLOH9p_heatmapdata = lossLOH9p_geneci %>%
  filter(gene != "CDKN2A.p14arf") %>%
  select(gene, all_mle, group) %>%
  mutate(scaled_mle = sqrt(all_mle))

lossLOH9p_heatmapdata$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", lossLOH9p_heatmapdata$gene)

lossLOH9p_heatmapdata = lossLOH9p_heatmapdata %>%
  filter(lossLOH9p_heatmapdata$gene %in% arm_heatmap_genelist$gene)

lossLOH9p_heatmapdata$gene = factor(lossLOH9p_heatmapdata$gene, levels = c("CDKN2A", "ZBTB24", "PIK3CA", "B2M", "DUSP22", 
                                                                           "SMAD4", "ARID1A", "TP53", 
                                                                           "ARHGAP35", "BCLAF1", "KRAS", "PTEN"))


lossLOH9p_heatmapdata$group = factor(lossLOH9p_heatmapdata$group, levels = c("No 9p loss LOH", "9p loss LOH"))

lossLOH9p_heatmapdata = lossLOH9p_heatmapdata %>%
  drop_na(gene)

lossLOH9p_heatmap = lossLOH9p_heatmapdata %>%
  ggplot(aes(x = group, y = gene, fill = log(all_mle))) +
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "Greens", direction = +1) +
  scale_x_discrete(labels = c("No clonal 9p loss-LOH", "Clonal 9p loss-LOH")) +
  theme_classic() +
  labs(x = "\nSCAA status", y = "Gene name\n", fill = "log(dN/dS)", 
       title = "Strength of dN/dS signal in subclonal driver genes\n in tumours with and without clonal 9p loss & LOH\n") +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "left",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.4, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))
print(lossLOH9p_heatmap)

lossLOH9p_heatmap_withnumbers = lossLOH9p_heatmap +
  geom_label(aes(label = sprintf("%.2f", all_mle)),
             colour = "black",
             fill = "white", alpha = 0.8,
             label.r = unit(0, "npc"),
             size = 7)

print(lossLOH9p_heatmap_withnumbers)

ggsave("plot22.png", lossLOH9p_heatmap_withnumbers, width = 12, height = 10, dpi = 300)

### Figure 22b: Driver frequency bar chart

lossLOH9p_tumourids = unique(arm_9p_loss_LOH$sampleID)

armdriverfreq_genelist = unique(arm_heatmap_genelist$gene)

lossLOH9p_driverfreq_muttable = arm_muttable %>%
  filter(arm_muttable$Hugo_Symbol %in% armdriverfreq_genelist | Hugo_Symbol == "CDKN2A", arm_muttable$tumour_id %in% lossLOH9p_tumourids, (combTiming_SC == "subclonal")) %>%
  distinct(Hugo_Symbol, tumour_id, .keep_all = TRUE) %>%
  add_count(Hugo_Symbol) %>%
  mutate(arm_group = "9p loss LOH") %>%
  mutate(total_tumours = 22) %>%
  mutate(percent = (n/22)*100)

nolossLOH9p_tumourids = unique(no_9p_loss_LOH$sampleID)  

nolossLOH9p_driverfreq_muttable = arm_muttable %>%
  filter(arm_muttable$Hugo_Symbol %in% armdriverfreq_genelist | Hugo_Symbol == "CDKN2A", arm_muttable$tumour_id %in% nolossLOH9p_tumourids, (combTiming_SC == "subclonal")) %>%
  distinct(Hugo_Symbol, tumour_id, .keep_all = TRUE) %>%
  add_count(Hugo_Symbol) %>%
  mutate(arm_group = "No 9p loss LOH") %>%
  mutate(total_tumours = 157) %>%
  mutate(percent = (n/157)*100)

arm_driverfreq_muttable = rbind(lossLOH9p_driverfreq_muttable, nolossLOH9p_driverfreq_muttable) %>%
  select(Hugo_Symbol, percent, arm_group)

lossLOH9prows = data.frame(Hugo_Symbol = c("STX2", "PTEN", "KRAS", "BCLAF1", "ARHGAP35"),
                           percent = c(0, 0, 0, 0, 0),
                           arm_group = c("9p loss LOH", "9p loss LOH", "9p loss LOH", "9p loss LOH", "9p loss LOH"))

arm_driverfreq_muttable = rbind(arm_driverfreq_muttable, lossLOH9prows)

arm_driverfreq_muttable$Hugo_Symbol = factor(arm_driverfreq_muttable$Hugo_Symbol, levels = c("CDKN2A", "ZBTB24", "PIK3CA", "B2M", "DUSP22", 
                                                                                             "SMAD4", "ARID1A", "TP53", 
                                                                                             "ARHGAP35", "BCLAF1", "KRAS", "PTEN"))

arm_driverfreq_muttable = arm_driverfreq_muttable %>%
  drop_na(Hugo_Symbol)

arm_matrix_driverfreq = arm_driverfreq_muttable %>%
  distinct(arm_group, Hugo_Symbol, percent) %>%
  arrange(Hugo_Symbol)

arm_driverfreq_plot = arm_driverfreq_muttable %>%
  ggplot(aes(Hugo_Symbol, percent, fill = factor(arm_group, levels = c("9p loss LOH", "No 9p loss LOH")))) +
  geom_col(position = position_dodge(), width = 0.7, colour = "#525252", linewidth = 0.2) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 40), expand = c(0, 0)) +
  scale_fill_manual(values = c("#CCCCCC", "#33A02C"), breaks = c("No 9p loss LOH", "9p loss LOH"), labels = c("No clonal 9p loss & LOH", "Clonal 9p loss & LOH")) +
  guides(fill = guide_legend(byrow = TRUE, override.aes = list(color = NA))) +
  labs(x = "Gene name\n", y = "\nProportion of tumours (%)", fill = "SCAA status",
       title = "Frequency of subclonal driver genes in tumours with and without clonal 9p loss & LOH\n") +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.y = unit(0.5, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.35, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))+
  coord_flip()
print(arm_driverfreq_plot)

ggsave("plot22b.png", arm_driverfreq_plot, width = 12, height = 10, dpi = 300)
