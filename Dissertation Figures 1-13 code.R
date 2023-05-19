###################################################################
##################### Dissertation Figures ########################
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
library(openxlsx)
library(reshape2)

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
muttable_path <- "TRACERx421_data_code/20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_mutation_table.fst"
mut_driver_lung_path <- "TRACERx421_data_code/20221109_Tumour_evo_histories_DATA/gene_anno_df.fst"
clinical_data_path <- "TRACERx421_data_code/20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_all_tumour_df.rds"
all_gd_path <- 'TRACERx421_data_code/20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_gd_table_per_tumour.tsv'

# Load WGD data
gdtable <- fread(all_gd_path)

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

clin = clin %>%
  rename(tumour_id = tumour_id_muttable_cruk)

clin_muttable = merge(clin, muttable)

### Appendix Figure 1:
mut_burdentable = clin_muttable %>%
  add_count(tumour_id) %>%
  rename(total_muts = n) %>%
  filter(histology_3 != "Other") %>%
  group_by(histology_3) %>%
  mutate(mean_muts = mean(total_muts))

hist_sample_size = mut_burdentable %>%
  distinct(tumour_id, .keep_all = TRUE) %>%
  group_by(histology_3) %>%
  summarise(num=n())

mut_burden_plot = mut_burdentable %>%
  left_join(hist_sample_size) %>%
  mutate(axislabels = paste0(histology_3, "\n", "n = ", num)) %>%
  ggplot(aes(axislabels, total_muts, fill = histology_3)) +
  geom_violin(width = 1, alpha = 0.8) +
  geom_boxplot(width = 0.1, color="#666666", show.legend = FALSE) +
  theme_classic() +
  scale_fill_manual(values = c("#FFFF99", "#CAB2D6")) +
  stat_summary(fun = "mean", geom = "point", colour = "#E31A1C", size = 6, alpha = 0.8, show.legend = FALSE) +
  annotate("label", x = c(1.33,2.39), y = c(2400, 2200),
           label = c("mean\n= 1994", "mean\n= 1050"), size = 5)+
  annotate("segment", x = 1, xend = 1.23, y = 1994, yend = 2400, size = 1, alpha = 0.7, linetype = 2, lwd = 0.8) +
  annotate("segment", x = 2, xend = 2.29, y = 1050, yend = 2200, size = 1, alpha = 0.7, linetype = 2, lwd = 0.8) +
  labs(x = '\nHistology', y = 'Number of SPMs per tumour\n', fill = "Histology") +
  ggtitle("Distribution of somatic point mutations per tumour in TRACERx421\n\n") +
  scale_y_continuous(limits = c(0, 12000), expand = c(0, 0)) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black", vjust = -1),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"))
print(mut_burden_plot)

ggsave("plota1.png", mut_burden_plot, width = 12, height = 10, dpi = 300)


### Figure 5: Number of lung drivers per tumour plot
drivers_muttable = clin_muttable %>%
  filter(DriverMut == TRUE) %>%
  add_count(tumour_id) %>%
  rename(total_drivers = n) %>%
  filter(histology_3 != "Other") %>%
  group_by(histology_3) %>%
  mutate(mean_drivermuts = mean(total_drivers))

drivers_pertumour_plot = drivers_muttable %>%
  left_join(hist_sample_size) %>%
  mutate(axislabels = paste0(histology_3, "\n", "n = ", num)) %>%
  ggplot(aes(axislabels, total_drivers, fill = histology_3)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#FFFF99", "#CAB2D6")) +
  stat_summary(fun = "mean", geom = "point", colour = "#E31A1C", size = 6, alpha = 0.8, show.legend = FALSE) +
  annotate("label", x = c(1.23,2.25), y = c(15, 14),
           label = c("mean = 11", "mean = 8"), size = 5)+
  annotate("segment", x = 1, xend = 1.23, y = 11, yend = 13.7, size = 1, alpha = 0.7, linetype = 2, lwd = 0.8) +
  annotate("segment", x = 2, xend = 2.25, y = 8.5, yend = 12.5, size = 1, alpha = 0.7, linetype = 2, lwd = 0.8) +
  labs(x = '\nHistology', y = 'Number of driver mutations per tumour\n', fill = "Histology") +
  scale_y_continuous(limits = c(0, 70), expand = c(0, 0)) +
  ggtitle("Number of driver mutations per tumour in TRACERx421\n\n") +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black", vjust = -1),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"))
print(drivers_pertumour_plot)

ggsave("plot5.png", drivers_pertumour_plot, width = 12, height = 10, dpi = 300)

### Figure 6: Driver clonality

drivertiming = drivers_muttable %>%
  filter(combTiming_SC != "NA") %>%
  group_by(histology_3, combTiming_SC) %>%
  tally() %>%
  mutate(percent = (n/sum(n))*100) %>%
  mutate(combTiming_SC = factor(combTiming_SC, levels = c("early", "late", "Unknown", "subclonal")))

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

drivertiming <- drivertiming %>%
  group_by(histology_3) %>%
  mutate(percent = round_prop(percent)) %>%
  ungroup()

drivertiming_plot = drivertiming %>%
  ggplot(aes(histology_3, y = percent, fill = combTiming_SC)) +
  geom_bar(stat = "identity", linewidth = 0.6, colour = "#525252", width = 0.7) +
  geom_text(aes(label = paste0(sprintf("%1.2f", percent),"%")),
            position=position_stack(vjust=0.5), colour="white", size = 8) +
  theme_classic() +
  scale_fill_manual(values = c("#92c5de", "#4393c3", "#969696", "#D73027"), breaks = c("early", "late", "Unknown", "subclonal"), labels = c("Early clonal", "Late clonal", "Unknown", "Subclonal")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  ggtitle("Evolutionary timing of somatic driver mutations in TRACERx421\n\n") +
  labs(x = '\nHistology\n', y = 'Proportion of somatic driver mutations (%)\n', fill = "Timing of mutation") +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
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
print(drivertiming_plot)

ggsave("plot6.png", drivertiming_plot, width = 12, height = 10, dpi = 300)

### Figure 7: Selection of clonal and subclonal genes in LUAD and LUSC

# Create LUAD and LUSC groups
LUAD = clin_muttable %>%
  filter(histology_3 == "LUAD")

LUSC = clin_muttable %>%
  filter(histology_3 == "LUSC")

# Split LUAD and LUSC into clonal and subclonal mutations
LUAD_clonal_muts = LUAD %>%
  filter(PyCloneClonal_SC_TB_corrected == "C") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

LUAD_subclonal_muts = LUAD %>%
  filter(PyCloneClonal_SC_TB_corrected == "S") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

LUSC_clonal_muts = LUSC %>%
  filter(PyCloneClonal_SC_TB_corrected == "C") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

LUSC_subclonal_muts = LUSC %>%
  filter(PyCloneClonal_SC_TB_corrected == "S") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

# Run dNdS on each group

# dNdS LUAD
dndsout_LUAD_clonal = dndscv(LUAD_clonal_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUAD_clonal = dndsout_LUAD_clonal$globaldnds
globaldnds_LUAD_clonal$group = "Clonal"
globaldnds_LUAD_clonal$histology = "LUAD"

dndsout_LUAD_subclonal = dndscv(LUAD_subclonal_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUAD_subclonal = dndsout_LUAD_subclonal$globaldnds
globaldnds_LUAD_subclonal$group = "Subclonal"
globaldnds_LUAD_subclonal$histology = "LUAD"

# dNdS LUSC
dndsout_LUSC_clonal = dndscv(LUSC_clonal_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUSC_clonal = dndsout_LUSC_clonal$globaldnds
globaldnds_LUSC_clonal$group = "Clonal"
globaldnds_LUSC_clonal$histology = "LUSC"

dndsout_LUSC_subclonal = dndscv(LUSC_subclonal_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUSC_subclonal = dndsout_LUSC_subclonal$globaldnds
globaldnds_LUSC_subclonal$group = "Subclonal"
globaldnds_LUSC_subclonal$histology = "LUSC"

# Plot clonal vs subclonal global dnds in LUAD and LUSC
combined_global_LUAD_LUSC = rbind(globaldnds_LUAD_clonal, globaldnds_LUAD_subclonal, globaldnds_LUSC_clonal, globaldnds_LUSC_subclonal)

global_LUAD_LUSC_plot = combined_global_LUAD_LUSC %>%
  filter(name == "wall") %>%
  ggplot(aes(histology, mle)) +
  geom_pointrange(aes(ymin = cilow, ymax = cihigh, colour = group), position = position_dodge(width = 0.2))+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  guides(colour = guide_legend(byrow = TRUE)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_color_manual(values = c("#4575B4", "#D73027"), name = 'Timing of mutations', labels = c("Clonal", "Subclonal")) +
  labs(x = '\nHistology', y = 'dN/dS driver genes\n', colour = 'Timing of mutations', title = "Cohort-level selection (dN/dS) of driver genes throughout tumour evolution\n") +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black", vjust = -1),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"))
print(global_LUAD_LUSC_plot)

ggsave("plot7.png", global_LUAD_LUSC_plot, width = 12, height = 10, dpi = 300)

### Figure 8: dNdS of specific clonal and subclonal genes in LUAD and LUSC

# dNdS of selected clonal and subclonal genes in LUAD and LUSC
sel_cv_LUAD_clonal = dndsout_LUAD_clonal$sel_cv %>%
  mutate(histology = "LUAD") %>%
  mutate(clonality = "Clonal")

sel_cv_LUAD_subclonal = dndsout_LUAD_subclonal$sel_cv %>%
  mutate(histology = "LUAD") %>%
  mutate(clonality = "Subclonal")

sel_cv_LUSC_clonal = dndsout_LUSC_clonal$sel_cv %>%
  mutate(histology = "LUSC") %>%
  mutate(clonality = "Clonal")

sel_cv_LUSC_subclonal = dndsout_LUSC_subclonal$sel_cv %>%
  mutate(histology = "LUSC") %>%
  mutate(clonality = "Subclonal")

signif_genes_clonality = rbind(sel_cv_LUAD_clonal, sel_cv_LUAD_subclonal, sel_cv_LUSC_clonal, sel_cv_LUSC_subclonal) %>%
  filter(qglobal_cv < 0.05) %>%
  select(gene_name, qglobal_cv, histology, clonality) %>%
  rename(gene = gene_name) %>%
  arrange(gene)

# Geneci clonal and subclonal genes in LUAD and LUSC
geneci_LUAD_clonal = geneci(dndsout_LUAD_clonal, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(histology = "LUAD") %>%
  mutate(clonality = "Clonal")

geneci_LUAD_subclonal = geneci(dndsout_LUAD_subclonal, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(histology = "LUAD") %>%
  mutate(clonality = "Subclonal")

geneci_LUSC_clonal = geneci(dndsout_LUSC_clonal, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(histology = "LUSC") %>%
  mutate(clonality = "Clonal")

geneci_LUSC_subclonal = geneci(dndsout_LUSC_subclonal, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(histology = "LUSC") %>%
  mutate(clonality = "Subclonal")

# Combine geneci clonal and subclonal in LUAD and LUSC
combined_geneci_clonality = rbind(geneci_LUAD_clonal, geneci_LUAD_subclonal, geneci_LUSC_clonal, geneci_LUSC_subclonal)

filtered_geneci_clonality = combined_geneci_clonality %>%
  filter(combined_geneci_clonality$gene %in% signif_genes_clonality$gene) %>%
  group_by(gene, histology) %>%
  mutate(odds_ratio = (all_mle[ clonality == 'Clonal'] / all_mle[ clonality == 'Subclonal'])) %>%
  arrange(gene) %>%
  mutate(drivertiming = case_when(odds_ratio > 2 ~ "Clonal favoured",
                                  odds_ratio <= 2 & odds_ratio >= 0.5 ~ "Clonal and subclonal",
                                  odds_ratio < 0.5 ~ "Subclonal favoured"))

qglobal_geneci_clonality = merge(filtered_geneci_clonality, signif_genes_clonality, all.x = TRUE) %>%
  arrange(histology)

LUAD_qglobal_geneci_clonality = qglobal_geneci_clonality %>%
  filter(histology == "LUAD", gene != "CDKN2A.p14arf") %>%
  mutate(drivertiming = factor(drivertiming, levels = c("Clonal favoured", "Clonal and subclonal", "Subclonal favoured"))) %>%
  arrange(drivertiming, gene)

LUAD_qglobal_geneci_clonality$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", LUAD_qglobal_geneci_clonality$gene)

LUAD_qglobal_geneci_clonality$gene = factor(LUAD_qglobal_geneci_clonality$gene, levels = c("KRAS", "STK11", "KEAP1", "CDKN2A", 
                                                                                           "EGFR", "SETD2", "NF1", "B2M",
                                                                                           "CMTR2", "BRAF", "NFE2L2", "RASA1", "TP53", "DUSP22", 
                                                                                           "FAM78B", "SMARCA4", "MGA", "PIK3CA", "ARID1A", "RB1", 
                                                                                           "ARHGAP35", "BCLAF1", "FAT1", "KMT2D", "TLR4", "SMAD4", "SFTPB", "PTEN"))
LUSC_qglobal_geneci_clonality = qglobal_geneci_clonality %>%
  filter(histology == "LUSC") %>%
  mutate(drivertiming = factor(drivertiming, levels = c("Clonal favoured", "Clonal and subclonal", "Subclonal favoured"))) %>%
  arrange(drivertiming, gene)

LUSC_qglobal_geneci_clonality$all_high = ifelse(abs(LUSC_qglobal_geneci_clonality$all_mle) == 0, LUSC_qglobal_geneci_clonality$all_high == 0, LUSC_qglobal_geneci_clonality$all_high)

LUSC_qglobal_geneci_clonality$gene = factor(LUSC_qglobal_geneci_clonality$gene, levels = c("KRAS", "STK11", "KEAP1", "CDKN2A.p16INK4a", 
                                                                                           "EGFR", "CDKN2A.p14arf", "SETD2", "NF1", "B2M",
                                                                                           "CMTR2", "BRAF", "NFE2L2", "RASA1", "TP53", "DUSP22", 
                                                                                           "FAM78B", "SMARCA4", "MGA", "PIK3CA", "ARID1A", "RB1", 
                                                                                           "ARHGAP35", "BCLAF1", "FAT1", "KMT2D", "TLR4", "SMAD4", "SFTPB", "PTEN"))

alt_LUSC_qglobal_geneci_clonality = qglobal_geneci_clonality %>%
  filter(histology == "LUSC", gene != "CDKN2A.p14arf") %>%
  mutate(drivertiming = factor(drivertiming, levels = c("Clonal favoured", "Clonal and subclonal", "Subclonal favoured"))) %>%
  arrange(drivertiming, gene)

alt_LUSC_qglobal_geneci_clonality$all_high = ifelse(abs(alt_LUSC_qglobal_geneci_clonality$all_mle) == 0, alt_LUSC_qglobal_geneci_clonality$all_high == 0, alt_LUSC_qglobal_geneci_clonality$all_high)

alt_LUSC_qglobal_geneci_clonality$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", alt_LUSC_qglobal_geneci_clonality$gene)

alt_LUSC_qglobal_geneci_clonality$gene = factor(alt_LUSC_qglobal_geneci_clonality$gene, levels = c("CDKN2A", "TP53", "NFE2L2", "PTEN", 
                                                                                                   "KMT2D", "RB1", "PIK3CA", "FAT1", "KEAP1",
                                                                                                   "RASA1", "TLR4", "CMTR2", "BRAF", "EGFR", "FAM78B", "KRAS", "MGA", "NF1", "SMAD4", "SMARCA4", "STK11", 
                                                                                                   "ARHGAP35", "ARID1A", "SETD2", "B2M", "BCLAF1", "DUSP22", "SFTPB"))
# Plot LUAD geneci
all_LUAD_geneci_plot = LUAD_qglobal_geneci_clonality %>%
  ggplot(aes(all_mle, gene)) +
  geom_pointrange(aes(xmin = all_low, xmax = all_high, colour = clonality), position = position_dodge(width = 0.2))+
  geom_vline(xintercept = 1, linetype = 2)+
  theme_classic()+
  guides(colour = guide_legend(byrow = TRUE)) +
  scale_color_manual(values = c("#4575B4", "#D73027"), name = 'Timing of mutations') +
  scale_x_log10() +
  scale_y_discrete(limits = rev(levels(LUAD_qglobal_geneci_clonality$gene))) +
  labs(y = 'Gene name\n', x = '\ndN/dS driver genes\n', title = 'A) LUAD\n\n') +
  theme(axis.text.y = element_text(colour = c(rep("#E31A1C", times = 3), rep("#6A3D9A", times = 13), rep("#1F78B4", times = 12)),
                                   face = c(rep('plain', times = 15), 'bold', rep('plain', times = 8), rep('bold', times = 4)),
                                   size = 12),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "none",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        aspect.ratio = 1.1)
print(all_LUAD_geneci_plot)

ggsave("plot8a.png", all_LUAD_geneci_plot, width = 12, height = 10, dpi = 300)

# Plot LUSC geneci
alt_LUSC_geneci_plot = alt_LUSC_qglobal_geneci_clonality %>%
  ggplot(aes(all_mle, gene)) +
  geom_pointrange(aes(xmin = all_low, xmax = all_high, colour = clonality), position = position_dodge(width = 0.2))+
  geom_vline(xintercept = 1, linetype = 2)+
  theme_classic()+
  scale_color_manual(values = c("#4575B4", "#D73027"), name = 'Timing of mutations') +
  scale_x_log10() +
  guides(colour = guide_legend(byrow = TRUE)) +
  scale_y_discrete(limits = rev(levels(alt_LUSC_qglobal_geneci_clonality$gene))) +
  labs(y = 'Gene name\n', x = '\ndN/dS driver genes\n', title = 'B) LUSC\n\n') +
  theme(axis.text.y = element_text(colour = c(rep("#E31A1C", times = 4), rep("#6A3D9A", times = 3), rep("#1F78B4", times = 21)),
                                   face = c(rep('plain', times = 26), rep('bold', times = 2)),
                                   size = 12),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "none",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        aspect.ratio = 1.1)
print(alt_LUSC_geneci_plot)

ggsave("plot8b.png", alt_LUSC_geneci_plot, width = 12, height = 10, dpi = 300)

### Figure 9: Global dNdS of early clonal, late clonal and subclonal mutations in LUAD and LUSC

# Create groups to run dNdS
LUAD_early_clonal_muts = LUAD %>%
  filter(combTiming_SC == "early") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

LUAD_late_clonal_muts = LUAD %>%
  filter(combTiming_SC == "late") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

LUSC_early_clonal_muts = LUSC %>%
  filter(combTiming_SC == "early") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

LUSC_late_clonal_muts = LUSC %>%
  filter(combTiming_SC == "late") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

# Run dNdS
# dNdS LUAD
dndsout_LUAD_earlyc = dndscv(LUAD_early_clonal_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUAD_earlyc = dndsout_LUAD_earlyc$globaldnds
globaldnds_LUAD_earlyc$group = "Early clonal drivers"
globaldnds_LUAD_earlyc$histology = "LUAD"

dndsout_LUAD_latec = dndscv(LUAD_late_clonal_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUAD_latec = dndsout_LUAD_latec$globaldnds
globaldnds_LUAD_latec$group = "Late clonal"
globaldnds_LUAD_latec$histology = "LUAD"

# dNdS LUSC
dndsout_LUSC_earlyc = dndscv(LUSC_early_clonal_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUSC_earlyc = dndsout_LUSC_earlyc$globaldnds
globaldnds_LUSC_earlyc$group = "Early clonal drivers"
globaldnds_LUSC_earlyc$histology = "LUSC"

dndsout_LUSC_latec = dndscv(LUSC_late_clonal_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUSC_latec = dndsout_LUSC_latec$globaldnds
globaldnds_LUSC_latec$group = "Late clonal"
globaldnds_LUSC_latec$histology = "LUSC"

# Plot clonal vs subclonal global dnds in LUAD and LUSC
global_timing_LUAD_LUSC = rbind(globaldnds_LUAD_earlyc, globaldnds_LUAD_latec, globaldnds_LUAD_subclonal, 
                                globaldnds_LUSC_earlyc, globaldnds_LUSC_latec, globaldnds_LUSC_subclonal)

global_timing_LUAD_LUSC_plot = global_timing_LUAD_LUSC %>%
  filter(name == "wall") %>%
  ggplot(aes(histology, mle)) +
  geom_pointrange(aes(ymin = cilow, ymax = cihigh, colour = group), position = position_dodge(width = 0.2))+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_color_manual(values = c("#92c5de", "#4393c3", "#D73027"), name = 'Timing of mutation')+
  labs(x = '\nHistology', y = 'dN/dS driver genes\n', colour = 'Timing of mutations', title = "Cohort-level selection (dN/dS) of driver genes relative to clonal WGD\n") +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black", vjust = -1),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"),
        aspect.ratio = 1)
print(global_timing_LUAD_LUSC_plot)

ggsave("plot9.png", global_timing_LUAD_LUSC_plot, width = 12, height = 10, dpi = 300)

### Figure 10: Proportion of WGD events in TRACERx421 cohort

# Type and frequency of WGD events in TRACERx421 cohort
clin_gd_muttable = merge(clin_muttable, gdtable) %>%
  mutate(WGD_status = case_when(GD_statuses == "0" ~ "No WGD",
                                num_clonal_gds != "0" & num_subclonal_gds == "0" ~ "Only clonal WGD",
                                num_clonal_gds != "0" & num_subclonal_gds != "0" ~ "Clonal and subclonal WGD",
                                num_clonal_gds == "0" & num_subclonal_gds != "0" ~ "Only subclonal WGD"))

proportion_gd = clin_gd_muttable %>%
  distinct(tumour_id, .keep_all = TRUE) %>%
  group_by(histology_3, WGD_status) %>%
  tally() %>%
  mutate(percent = (n/sum(n))*100) %>%
  mutate(freq = sum(n)) %>%
  filter(histology_3 != "Other") %>%
  mutate(histology_labels = paste0(histology_3, "\n", "n = ", freq))

histology_labels = proportion_gd %>%
  distinct(freq, .keep_all = TRUE)

proportion_gd <- proportion_gd %>%
  group_by(histology_3) %>%
  mutate(percent = round_prop(percent)) %>%
  ungroup()

proportion_gd_plot = proportion_gd %>%
  ggplot(aes(x = histology_labels, y = percent, fill = factor(WGD_status, levels = c("Only clonal WGD", "Only subclonal WGD", "Clonal and subclonal WGD", "No WGD")))) +
  geom_bar(stat = "identity", linewidth = 0.6, colour = "#525252", width = 0.7) +
  geom_text(aes(label = paste0(sprintf("%1.2f", percent),"%")),
            position=position_stack(vjust=0.5), colour="white", size = 8) +
  theme_classic() +
  guides(fill = guide_legend(byrow = TRUE)) +
  scale_fill_manual(values = c("#A6CEE3", "#E31A1C", "#CAB2D6", "#969696"), name = 'WGD status') +
  ggtitle("Percentage of tumours according to WGD status in TRACERx421\n\n") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = '\nHistology', y = 'Proportion of tumours (%)\n', fill = "WGD status") +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black", vjust = -1),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"))
print(proportion_gd_plot)

ggsave("plot10.png", proportion_gd_plot, width = 12, height = 10, dpi = 300)

### Figure 11: Selection of lung driver mutations post WGD in LUAD and LUSC

# Create subgroups for LUAD and LUSC
LUAD_no_gd = clin_gd_muttable %>%
  filter(histology_3 == "LUAD", GD_statuses == "0", !(PyCloneClonal_SC_TB_corrected == "C" & combTiming_SC == "early")) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

LUAD_only_clonal = clin_gd_muttable %>%
  filter(histology_3 == "LUAD", num_clonal_gds != "0" & num_subclonal_gds == "0", !(PyCloneClonal_SC_TB_corrected == "C" & combTiming_SC == "early")) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

LUSC_no_gd = clin_gd_muttable %>%
  filter(histology_3 == "LUSC", GD_statuses == "0", !(PyCloneClonal_SC_TB_corrected == "C" & combTiming_SC == "early")) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

LUSC_only_clonal = clin_gd_muttable %>%
  filter(histology_3 == "LUSC", num_clonal_gds != "0" & num_subclonal_gds == "0", !(PyCloneClonal_SC_TB_corrected == "C" & combTiming_SC == "early")) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

# Run dNdS for LUAD and LUSC WGD subgroups
dndsout_LUAD_no_gd = dndscv(LUAD_no_gd, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUAD_no_gd = dndsout_LUAD_no_gd$globaldnds
globaldnds_LUAD_no_gd$group = "Late clonal and subclonal drivers (No WGD)"
globaldnds_LUAD_no_gd$histology = "LUAD"

dndsout_LUAD_only_clonal = dndscv(LUAD_only_clonal, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUAD_only_clonal = dndsout_LUAD_only_clonal$globaldnds
globaldnds_LUAD_only_clonal$group = "Late clonal and subclonal drivers (Clonal WGD)"
globaldnds_LUAD_only_clonal$histology = "LUAD"

dndsout_LUSC_no_gd = dndscv(LUSC_no_gd, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUSC_no_gd = dndsout_LUSC_no_gd$globaldnds
globaldnds_LUSC_no_gd$group = "Late clonal and subclonal drivers (No WGD)"
globaldnds_LUSC_no_gd$histology = "LUSC"

dndsout_LUSC_only_clonal = dndscv(LUSC_only_clonal, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
globaldnds_LUSC_only_clonal = dndsout_LUSC_only_clonal$globaldnds
globaldnds_LUSC_only_clonal$group = "Late clonal and subclonal drivers (Clonal WGD)"
globaldnds_LUSC_only_clonal$histology = "LUSC"

# Early clonal mutations (pre WGD)
clonal_gd_earlyc_muts = clin_gd_muttable %>%
  filter(num_clonal_gds != "0", PyCloneClonal_SC_TB_corrected == "C" & combTiming_SC == "early") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

dndsout_clonal_gd_earlyc = dndscv(clonal_gd_earlyc_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)

sel_cv_clonal_gd_earlyc = dndsout_clonal_gd_earlyc$sel_cv %>%
  mutate(WGD_status = "Early clonal drivers")

geneci_clonal_gd_earlyc = geneci(dndsout_clonal_gd_earlyc, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(WGD_status = "Early clonal drivers")

signif_genes_earlyc = sel_cv_clonal_gd_earlyc %>%
  filter(qglobal_cv < 0.05) %>%
  select(gene_name, qglobal_cv, WGD_status)

filtered_geneci_clonalgd_earlyc = geneci_clonal_gd_earlyc %>%
  filter(geneci_clonal_gd_earlyc$gene %in% signif_genes_earlyc$gene_name, gene != "CDKN2A.p14arf")

filtered_geneci_clonalgd_earlyc$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", filtered_geneci_clonalgd_earlyc$gene)

# Late clonal mutations (post WGD)
clonal_gd_latec_muts = clin_gd_muttable %>%
  filter(num_clonal_gds != "0", PyCloneClonal_SC_TB_corrected == "C" & combTiming_SC == "late") %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

dndsout_clonal_gd_latec = dndscv(clonal_gd_latec_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)

sel_cv_clonal_gd_latec = dndsout_clonal_gd_latec$sel_cv %>%
  mutate(WGD_status = "Late clonal drivers")

geneci_clonal_gd_latec = geneci(dndsout_clonal_gd_latec, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(WGD_status = "Late clonal drivers")

signif_genes_latec = sel_cv_clonal_gd_latec %>%
  filter(qglobal_cv < 0.05) %>%
  select(gene_name, qglobal_cv, WGD_status)

filtered_geneci_clonalgd_latec = geneci_clonal_gd_latec %>%
  filter(geneci_clonal_gd_latec$gene %in% signif_genes_latec$gene_name, gene != "CDKN2A.p14arf")

filtered_geneci_clonalgd_latec$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", filtered_geneci_clonalgd_latec$gene)

# Plot global dNdS for LUAD and LUSC WGD subgroups
global_gd_LUAD_LUSC = rbind(globaldnds_LUAD_earlyc, globaldnds_LUSC_earlyc, globaldnds_LUAD_latec, globaldnds_LUSC_latec, globaldnds_LUAD_no_gd, globaldnds_LUAD_only_clonal, globaldnds_LUSC_no_gd, globaldnds_LUSC_only_clonal) %>%
  filter(name == "wall")

gd_LUAD_LUSC_plot = global_gd_LUAD_LUSC %>%
  filter(name == "wall") %>%
  ggplot(aes(histology, mle)) +
  geom_pointrange(aes(ymin = cilow, ymax = cihigh, colour = factor(group, levels = c("Early clonal drivers", "Late clonal", "Late clonal and subclonal drivers (No WGD)", "Late clonal and subclonal drivers (Clonal WGD)"))), position = position_dodge(width = 0.2))+  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  scale_color_manual(values = c("#A6CEE3", "#6A3D9A", "#969696", "#1F78B4"), name = 'Timing of driver genes', labels = c("Early clonal", "Late clonal", "Late (No WGD)", "Late (Clonal WGD)")) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  guides(colour = guide_legend(byrow = TRUE)) +
  labs(x = '\nHistology', y = 'dN/dS driver genes\n', title = "Cohort-level selection (dN/dS) of driver genes\n in tumours with and without clonal WGD\n\n")+
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black", vjust = -1),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.2, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 3, 1), "lines"),
        aspect.ratio = 1)
print(gd_LUAD_LUSC_plot)

ggsave("plot11.png", gd_LUAD_LUSC_plot, width = 12, height = 10, dpi = 300)

# Create combined LUAD/LUSC WGD groups for gene level dNdS
no_gd_muts = clin_gd_muttable %>%
  filter(GD_statuses == "0", !(PyCloneClonal_SC_TB_corrected == "C" & combTiming_SC == "early")) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

clonal_gd_muts = clin_gd_muttable %>%
  filter(num_clonal_gds != "0", !(PyCloneClonal_SC_TB_corrected == "C" & combTiming_SC == "early")) %>%
  select(tumour_id, chr, start, ref, var) %>%
  rename(sampleID = tumour_id,
         pos = start,
         alt = var) %>%
  mutate(chr = gsub("chr","", chr))

# Run dNdS for WGD groups
dndsout_no_gd = dndscv(no_gd_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)

dndsout_clonal_gd = dndscv(clonal_gd_muts, gene_list = mut_drivers_cd, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)

# Gene level dnds for combined LUAD/LUSC WGD groups (LUAD and LUSC combined to give statistical power)
sel_cv_no_gd = dndsout_no_gd$sel_cv %>%
  mutate(WGD_status = "No WGD")

sel_cv_clonal_gd = dndsout_clonal_gd$sel_cv %>%
  mutate(WGD_status = "Clonal WGD")

signif_genes_WGD = rbind(sel_cv_no_gd, sel_cv_clonal_gd) %>%
  filter(qglobal_cv < 0.05) %>%
  select(gene_name, qglobal_cv, WGD_status)


# Geneci WGD groups
geneci_no_gd = geneci(dndsout_no_gd, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(WGD_status = "No WGD")

geneci_clonal_gd = geneci(dndsout_clonal_gd, gene_list = mut_drivers_cd, level = 0.95) %>%
  mutate(WGD_status = "Clonal WGD")

# Combine geneci WGD groups
combined_geneci_gd = rbind(geneci_no_gd, geneci_clonal_gd) %>%
  arrange(gene)

filtered_geneci_gd = combined_geneci_gd %>%
  filter(combined_geneci_gd$gene %in% signif_genes_WGD$gene_name, gene != "CDKN2A.p14arf")

filtered_geneci_gd$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", filtered_geneci_gd$gene)

filtered_geneci_gd = filtered_geneci_gd %>%
  group_by(gene) %>%
  mutate(odds_ratio = (all_mle[ WGD_status == 'No WGD'] / all_mle[ WGD_status == 'Clonal WGD'])) %>%
  arrange(gene) %>%
  mutate(wgdfavoured = case_when(odds_ratio > 2 ~ "No WGD favoured",
                                 odds_ratio <= 2 & odds_ratio >= 0.5 ~ "No WGD and clonal WGD",
                                 odds_ratio < 0.5 ~ "Clonal WGD favoured"))

qglobal_signif_genes_WGD = signif_genes_WGD %>%
  select(gene_name, qglobal_cv) %>%
  rename(gene = gene_name) %>%
  filter(gene != "CDKN2A.p14arf")

qglobal_signif_genes_WGD$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", qglobal_signif_genes_WGD$gene)

qglobal_filtered_geneci_gd = merge(filtered_geneci_gd, qglobal_signif_genes_WGD)

write.xlsx(qglobal_filtered_geneci_gd, file = "Appendixtab3.xlsx")
write.xlsx(qglobal_filtered_geneci_clonalarms, file = "Appendixtab4.xlsx")
write.xlsx(qglobal_filtered_geneci_9plossLOH, file = "Appendixtab5.xlsx")

qglobal_filtered_geneci_gd = qglobal_filtered_geneci_gd %>%
  group_by(wgdfavoured, qglobal_cv) %>%
  arrange(wgdfavoured, qglobal_cv)

qglobal_filtered_geneci_gd$wgdfavoured = factor(qglobal_filtered_geneci_gd$wgdfavoured, levels = c("No WGD favoured", "No WGD and clonal WGD", "Clonal WGD"))

filtered_geneci_gd$all_high = ifelse(abs(filtered_geneci_gd$all_mle) == 0, filtered_geneci_gd$all_high == 0, filtered_geneci_gd$all_high)

filtered_geneci_gd$gene = factor(filtered_geneci_gd$gene, levels = c("SFTPB", 
                                                                     "BCLAF1", "PIK3CA", "DUSP22", "EGFR", "TLR4", 
                                                                     "NFE2L2", "FAM78B", "ARHGAP35", "KMT2D", "ARID1A", "B2M", "SMAD4", "PTEN", "CDKN2A", "RB1", "KEAP1", "STK11", "TP53", "KRAS"))


### Figure 12: Plot gene level dN/dS for post WGD genes (LUAD and LUSC combined)
gd_geneci_plot = filtered_geneci_gd %>%
  ggplot(aes(gene, all_mle)) +
  geom_pointrange(aes(ymin = all_low, ymax = all_high, colour = WGD_status), position = position_dodge(width = 0.2))+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  scale_color_manual(values = c("#969696", "#1F78B4"), name = 'WGD status', breaks = c("No WGD", "Clonal WGD")) +
  scale_y_log10() +
  labs(x = 'Gene name\n', y = '\ndN/dS late driver genes\n', title = 'Gene-level selection (dN/dS) of late driver genes\n in tumours with and without clonal WGD\n\n') +
  theme(axis.text.y = element_text(colour = c(rep("#1F78B4", times = 1), rep("#6A3D9A", times = 5), rep("#E31A1C", times = 14)),
                                   face = c(rep('plain', times = 12), 'bold', rep('plain', times = 3), 'bold', 'bold', 'plain', 'bold'),
                                   size = 12),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.3, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        aspect.ratio = 1)+
  coord_flip()
print(gd_geneci_plot)

ggsave("plot12.png", gd_geneci_plot, width = 12, height = 10, dpi = 300)


### Figure 13: Heatmap of strength of selection

heatmap_filtered_geneci_gd = filtered_geneci_gd %>%
  select(gene, mis_mle, tru_mle, all_mle, mis_low, tru_low, all_low, mis_high, tru_high, all_high, WGD_status)

heatmap_genelist = rbind(filtered_geneci_clonalgd_earlyc, heatmap_filtered_geneci_gd) %>%
  arrange(gene) %>%
  distinct(gene)

all_wgd_geneci = rbind(geneci_clonal_gd, geneci_no_gd) %>%
  filter(gene != "CDKN2A.p14arf")

all_wgd_geneci$gene = gsub("CDKN2A.p16INK4a", "CDKN2A", all_wgd_geneci$gene)

heatmapdata = all_wgd_geneci %>%
  filter(all_wgd_geneci$gene %in% heatmap_genelist$gene) %>%
  select(gene, all_mle, WGD_status) %>%
  mutate(scaled_mle = sqrt(all_mle))

heatmapdata$WGD_status = factor(heatmapdata$WGD_status, levels = c("No WGD", "Clonal WGD"))

heatmapdata$gene = factor(heatmapdata$gene, levels = c("SFTPB", 
                                                       "BCLAF1", "PIK3CA", "DUSP22", "EGFR", "TLR4", 
                                                       "NFE2L2", "FAM78B", "ARHGAP35", "KMT2D", "ARID1A", "B2M", "SMAD4", "PTEN", "CDKN2A", "RB1", "KEAP1", "STK11", "TP53", "KRAS"))
heatmapdata = heatmapdata %>%
  filter(gene != "NA")

wgd_heatmap = heatmapdata %>%
  ggplot(aes(x = WGD_status, y = gene, fill = log(all_mle))) +
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "Blues", direction = +1) +
  scale_x_discrete(labels = c("Late drivers (No WGD)", "Late drivers (Clonal WGD)"), expand = ) +
  theme_classic() +
  labs(x = "\nWGD status", y = "Gene name\n", fill = "log(dN/dS)", title = "Strength of dN/dS signal in driver genes according to WGD status\n")+
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "left",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.4, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))
print(wgd_heatmap)

wgdheatmap_withnumbers = wgd_heatmap +
  geom_label(aes(label = sprintf("%.2f", all_mle)),
             colour = "black",
             fill = "white", alpha = 0.8,
             label.r = unit(0, "npc"),
             size = 5)

print(wgdheatmap_withnumbers)

ggsave("plot13a.png", wgdheatmap_withnumbers, width = 12, height = 10, dpi = 300)

## Frequency of subclonal drivers in clonal WGD vs no WGD
no_gd_tumourids = unique(no_gd_muts$sampleID)

driverfreq_genelist = unique(heatmap_genelist$gene)

no_wgd_driverfreq_muttable = muttable %>%
  filter(muttable$Hugo_Symbol %in% driverfreq_genelist | Hugo_Symbol == "CDKN2A", muttable$tumour_id %in% no_gd_tumourids, (combTiming_SC == "late" | combTiming_SC == "subclonal")) %>%
  add_count(Hugo_Symbol) %>%
  mutate(wgd_group = "No WGD") %>%
  mutate(total_tumours = 73) %>%
  mutate(percent = (n/73)*100)

clonal_gd_tumourids = unique(clonal_gd_muts$sampleID)  

clonal_wgd_driverfreq_muttable = muttable %>%
  filter(muttable$Hugo_Symbol %in% driverfreq_genelist | Hugo_Symbol == "CDKN2A", muttable$tumour_id %in% clonal_gd_tumourids, (combTiming_SC == "late" | combTiming_SC == "subclonal")) %>%
  add_count(Hugo_Symbol) %>%
  mutate(wgd_group = "Clonal WGD") %>%
  mutate(total_tumours = 181) %>%
  mutate(percent = (n/181)*100)

early_gd_tumourids = unique(clonal_gd_earlyc_muts$sampleID)

earlyc_wgd_driverfreq_muttable = muttable %>%
  filter(muttable$Hugo_Symbol %in% driverfreq_genelist | Hugo_Symbol == "CDKN2A", muttable$tumour_id %in% early_gd_tumourids, (combTiming_SC == "early")) %>%
  add_count(Hugo_Symbol) %>%
  mutate(wgd_group = "Early clonal WGD") %>%
  mutate(percent = (n/242)*100)

wgd_driverfreq_muttable = rbind(no_wgd_driverfreq_muttable, clonal_wgd_driverfreq_muttable)

wgd_driverfreq_muttable$Hugo_Symbol = factor(wgd_driverfreq_muttable$Hugo_Symbol, levels = c("SFTPB", 
                                                                                             "BCLAF1", "PIK3CA", "DUSP22", "EGFR", "TLR4", 
                                                                                             "NFE2L2", "FAM78B", "ARHGAP35", "KMT2D", "ARID1A", "B2M", "SMAD4", "PTEN", "CDKN2A", "RB1", "KEAP1", "STK11", "TP53", "KRAS"))
wgd_driverfreq_muttable = wgd_driverfreq_muttable %>%
  drop_na(Hugo_Symbol)

wgd_matrix_driverfreq = wgd_driverfreq_muttable %>%
  distinct(wgd_group, Hugo_Symbol, percent) %>%
  arrange(Hugo_Symbol)

wgd_driverfreq_plot = wgd_driverfreq_muttable %>%
  ggplot(aes(percent, Hugo_Symbol, fill = factor(wgd_group, levels = c("Clonal WGD", "No WGD")))) +
  geom_col(position = position_dodge(), width = 0.7, colour = "#525252", linewidth = 0.2) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 40), expand = c(0, 0)) +
  guides(fill = guide_legend(byrow = TRUE, override.aes = list(color = NA))) +
  scale_fill_manual(values = c("#CCCCCC", "#1F78B4"), breaks = c("No WGD", "Clonal WGD")) +
  labs(x = "\nProportion of tumours (%)", y = "Gene name\n", fill = "WGD status", title = "Frequency of late driver genes in tumours with and without clonal WGD\n") +
  theme(axis.text.y = element_text(colour = c(rep("black", times = 5), "#B15928", rep("black", times = 5), "#B15928", rep("black", times = 1), "#B15928", rep("black", times = 6)),
                                   face = c(rep('plain', times = 5), 'bold', rep('plain', times = 5), 'bold', rep('plain', times = 1), 'bold', rep('plain', times = 6)),
                                   size = 12),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.spacing.y = unit(0.4, "cm"),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.4, size = 20, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))
print(wgd_driverfreq_plot)

ggsave("plot13b.png", wgd_driverfreq_plot, width = 12, height = 10, dpi = 300)
