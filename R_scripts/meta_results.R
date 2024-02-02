#' Supplementary code to the article: 
#' Klinkovska et al. (2024) Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes. Biological Conservation.
#' 
#' Author: Klara Klinkovska, Michael Glaser, 2023-09-22
#' R version 4.3.0


# loading packages --------------------------------------------------------

library(plyr) # version 1.8.8

###### tables
library(janitor) # version 2.2.0
library(data.table) # version 1.14.8
library(readxl) # version 1.4.2

###### stats
library(modelsummary) # version 1.4.1
library(rstatix) # version 0.7.2

###### graphics
library(broom) # version 1.0.5
library(gridExtra) # version 2.3
library(ggpmisc) # version 0.5.2
library(corrplot) # version 0.92
library(ggforce) # version 0.4.1
library(ggfortify) # version 0.4.16
library(patchwork) # version 1.1.2
library(ggpubr) # version 0.6.0
library(svglite) # version 2.1.1
library(multcompView) # version 0.1-9
library(ggtext)

# clustering
library(car) # version 3.1-2
library(dtwclust) # version 5.5.12

library(tidyverse) # version 2.0.0

# plot settings -----------------------------------------------------------
path.plots <- "plots/"

iucn.col <- c("CR"="#FC8D59", "EN" = "#FEE08B", "VU"="#D9EF8B", "NT" = "#91CF60", 
              "LC" = "#1A9850", 'NA' = 'grey')

biogeo.col <- c("native"="#A9D08E", "archaeophyte"="#FFD966", "neophyte"="#ED1F39")

clust.col <- c("1" = "#00C19F", "2" = "#00B9E3", "3" = "#53B400", "4" = "#F8766D", "5" ="#619CFF")  

### Mytheme
t.size <- 3

my.theme <- theme(
  plot.title = element_text(size = unit(14, "points"), face = "bold"),
  axis.title = element_text(size = unit(10, "points"), face = "bold"),
  axis.text = element_text(size = unit(8, "points")), 
  legend.title = element_text(size = unit(10, "points"), face = "bold"), 
  legend.text = element_text(size = unit(8, "points")))

# functions ---------------------------------------------------------------

#' Pearson's chi-square test and Post-hoc tests for individual combinations of a category 
#' and species group (increasing vs decreasing) based on comparing their standardised residuals 
#' with the normal distribution of residuals (Beasley & Schumacker 1995) 

### for cathegorical variables stored in multiple columns
chi.test <- function(col1, col2, exclude, excluded, sim) {
  cont.tab <- resdf.clip.join %>% 
    select(change.lm, col1:col2) %>% 
    group_by(change.lm) %>% 
    summarize_at(vars(col1:col2), sum, na.rm = T) %>% 
    filter(change.lm != 0) %>% 
    column_to_rownames("change.lm") 
  
  if(exclude == T){cont.tab <- cont.tab %>% select(-all_of(excluded))}
  
  if(sim == T){
    chi.df <- chisq.test(cont.tab, simulate.p.value = T)
  }else{
    chi.df <- chisq.test(cont.tab)
  }
  
  chi.df$expected
  
  2*(1-pnorm(abs(chi.df$residuals))) %>%
    round(5) 
}

### for categorical variables stored in one column
chi.test2 <- function(col, exclude, excluded) {
  cont.tab <- resdf.clip.join %>% 
    select(change.lm, col) %>% 
    filter(change.lm != 0) %>% 
    group_by(change.lm) %>% 
    table() %>% 
    as.data.frame() %>% 
    pivot_wider(names_from = col, values_from = Freq) %>% 
    column_to_rownames("change.lm") 
  
  if(exclude == T) {cont.tab <- cont.tab %>% select(-all_of(excluded))}
  
  chi.df <- chisq.test(cont.tab)
  
  chi.df$expected
  
  2*(1-pnorm(abs(chi.df$residuals))) %>%
    round(4) 
  
}

### for differences between groups created by the time-series clustering
chi.test.clust <- function(col1, col2) {
  cont.tab <- resdf.clip.join2 %>% 
    select(cluster, col1:col2) %>% 
    group_by(cluster) %>% 
    summarize_at(vars(col1:col2), sum, na.rm = T) %>% 
    column_to_rownames("cluster") 
  
  chi.df <- chisq.test(cont.tab, simulate.p.value = T)
  print(chi.df$expected)
  print(chi.df$observed)
  
  2*(1-pnorm(abs(chi.df$residuals))) %>%
    round(5) %>% print()
} 

#' annotation plots for Fig. 3 and Appendix S8

# annotation plot for categorical variables stored in multiple columns, no label modification
annot.plot1 <- function(labtab, col1, col2, colvec, size.atext, labdist, title.left, my.theme) {
  
  labtab <- resdf.clip.join %>% 
    select(change.lm, col1:col2) %>% 
    group_by(change.lm) %>% 
    summarize_at(vars(col1:col2), sum, na.rm = T) %>%  
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(., "an.col") %>% 
    rename("neg" = "V1", "pos" = "V2") %>% 
    filter(an.col != "change.lm") 
  
  labtab$an.col <- fct_inorder(labtab$an.col)
  names(colvec)<-labtab$an.col
  
  plot.a <- ggplot(labtab)+
    theme_classic()+
    geom_bar(aes(x = -1*neg, y = an.col, fill = an.col), stat = "identity", width = 0.8)+
    geom_bar(aes(x = pos, y = an.col, fill = an.col), stat = "identity", width = 0.8)+
    scale_fill_manual(values = colvec, guide = "none")+
    geom_vline(xintercept = 0, linewidth = 0.5)+
    annotate("text", label = c("decreasing","increasing"), 
             x = c(-0.7*max(labtab$neg), 0.7*max(labtab$neg)),
             y = nrow(labtab) + 1,
             size = size.atext + 0.2)+
    scale_y_discrete(expand = expansion(add = c(2, 2)))+
    scale_x_continuous(expand = expansion(add = c(50, 50)))+
    geom_text(data = labtab, aes(y = an.col, x = -neg-labdist, label = round(neg, 2)), size = size.atext-0.5)+
    geom_text(data = labtab, aes(y = an.col, x = pos + labdist, label = round(pos, 2)), size = size.atext-0.5)+
    theme(panel.border = element_blank(),
          plot.background = element_rect(colour = NA, fill = NA, linewidth = 1), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_text(colour = "black"),
          axis.title.y = element_blank())+
    xlab("Number of species")+
    ggtitle(title.left)+
    my.theme
}


# annotation plot for categorical variables stored in multiple columns with modification of y-labels 
annot.plot <- function(labtab, col1, col2, colvec, size.atext, labdist, title.left, 
                       my.theme, label, exp.x = 50) {
  
  labtab <- resdf.clip.join %>% 
    select(change.lm, col1:col2) %>% 
    group_by(change.lm) %>% 
    summarize_at(vars(col1:col2), sum, na.rm = T) %>%  
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(., "an.col") %>% 
    rename("neg" = "V1", "pos" = "V2") %>% 
    filter(an.col != "change.lm") 
  
  labtab$an.col <- fct_inorder(labtab$an.col)
  names(colvec) <- labtab$an.col
  
  plot.a <- ggplot(labtab)+
    theme_classic()+
    geom_bar(aes(x = -1*neg, y = an.col, fill = an.col), stat = "identity", width = 0.8)+
    geom_bar(aes(x = pos, y = an.col, fill = an.col), stat = "identity", width = 0.8)+
    scale_fill_manual(values = colvec, guide = "none")+
    geom_vline(xintercept = 0, size = 0.5)+
    annotate("text", label = c("decreasing", "increasing"), 
             x = c(-0.7*max(labtab$neg), 0.7*max(labtab$neg)),
             y = nrow(labtab) + 1,
             size = size.atext + 0.2)+
    scale_y_discrete(expand = expansion(add = c(2,2)), labels = label)+
    scale_x_continuous(expand = expansion(add = c(exp.x, exp.x)))+
    geom_text(data = labtab, aes(y = an.col, x = -neg - labdist, label = round(neg, 2)), size = size.atext - 0.5)+
    geom_text(data = labtab, aes(y = an.col, x = pos + labdist, label = round(pos, 2)), size = size.atext - 0.5)+
    theme(panel.border = element_blank(),
          plot.background = element_rect(colour = NA, fill = NA, linewidth = 1), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_text(colour = "black"),
          axis.title.y = element_blank())+
    xlab("Number of species")+
    ggtitle(title.left)+
    my.theme
}

### annotation plot for categorical variables stored in one column
annot.plot.1col <- function(df,colvec, colname, title.left,
                            size.atext, labdist, label, exp.x = 50, annot.exp = 1){
  
  ### subset to current data frame, define analysis column
  colvec <- colvec[!names(colvec) == "x"]
  an.df <- df[, c(colname, "change.lm", "lm.rank")]
  an.df$an.col <- pull(an.df[, 1])
  
  ### create a table for labelling
  labtab <- dcast(setDT(data.frame(table(an.df$an.col, an.df$change.lm))), Var1~Var2)
  setDF(labtab)
  colnames(labtab) <- c("an.col", "neg", "pos")
  labtab <- labtab[order(match(labtab$an.col, names(colvec))),]
  
  ### make the annotation plot
  ggplot(an.df)+
    theme_classic()+
    geom_bar(aes(x = change.lm, y = an.col, fill = factor(an.col)), stat = "identity", width = 0.8)+
    scale_fill_manual(values = colvec, guide = "none")+
    geom_vline(xintercept = 0, linewidth = 0.5)+
    annotate("text", label = c("decreasing","increasing"), 
             x = c(-150, 150),
             y = nrow(labtab) + annot.exp,
             size = size.atext + 0.2)+
    scale_y_discrete(expand = expansion(add = c(2, 2)), limits = names(colvec[names(colvec) %in% an.df$an.col]), 
                     labels = label)+
    scale_x_continuous(expand = expansion(add = c(exp.x, exp.x)))+
    geom_text(data = labtab, aes(y = an.col, x = -neg - labdist, label = neg), size = size.atext - 0.5)+
    geom_text(data = labtab, aes(y = an.col, x = pos + labdist, label = pos), size = size.atext - 0.5)+
    theme(panel.border = element_blank(),
          plot.background = element_rect(colour = NA, fill = NA, size = 1), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_text(colour = "black"),
          axis.title.y = element_blank())+
    xlab("Number of species")+
    ggtitle(title.left)+
    my.theme
}

# loading data ------------------------------------------------------------

# load occupancy results
occ.res <- read_csv("results/occ_results_occupancy_60half_20230801.csv") %>% 
  left_join(read_csv("results/records_decades.csv"), by = c("specname" = "name_lat", "parnum" = "half_dec"))

# list of species with not converging results
n.conv.full <- read_csv("results/not_conv_all_230918.csv")

# clipping records, calculation of regression slopes and filtering species with significant trends 
resdf.clip <- occ.res %>% 
  left_join(read_csv("results/records_decades.csv") %>% 
              filter(pres != 0) %>% 
              group_by(name_lat) %>% 
              summarize(clip.min = min(half_dec), clip.max = max(half_dec)), by = c("specname" = "name_lat")) %>% 
  filter(parnum >= clip.min & parnum <= clip.max & !specname %in% n.conv.full$specname) %>% # clip results to the data
  group_by(specname) %>% 
  nest() %>%
  mutate(lm = map(data, ~lm(median~parnum, data = .x)), # linear models for the occupancy trends of individual species
         lm.tidy = map(lm, ~tidy(.x)), 
         lm.slope = map_dbl(lm.tidy, ~.x$estimate[[2]]), 
         lm.p = map_dbl(lm.tidy, ~.x$p.value[[2]])) %>% 
  select(specname, lm.slope, lm.p) %>%
  ungroup() %>% 
  filter(lm.p < 0.05) %>% # filter only species with significant linear trend
  mutate(change.lm = sign(lm.slope), 
         lm.rank = rank(lm.slope, ties.method = "random"))

# link with species characteristics ---------------------------------------

# load species characteristics data
alien <- read_csv2("traits/Pysek_Appendix_2_KK2.csv") # list of alien species
eiv <- read_csv2("traits/complexExport-Pladias-2022-10-27_KK.csv") # ecological indicator values
trait <- read_delim("traits/complexExport-Pladias-2023-04-19-Klinkovska_edit.csv") # species characteristics
iucn <- read_xlsx("traits/Cerveny_seznam_2017_2IUCN.xlsx") # Red list status

# join species trends with with their characteristics and modify columns
resdf.clip.join <- resdf.clip %>% 
  left_join(eiv %>% select(lat_name, EIV_light:dist_str_based), by = c("specname" = "lat_name")) %>% 
  left_join(trait, by = c("specname" = "lat_name")) %>% 
  left_join(iucn, by = c("specname" = "Taxon")) %>% 
  left_join(alien, by = c("specname" = "Taxon")) %>% 
  mutate_all(~ifelse(.x == "FALSE|NA", NA, .x)) %>% 
  mutate(alien = revalue(`Residence time`, replace = c("ar NE" = "archaeophyte",
                                                       "ar EM" = "archaeophyte",
                                                       "ar ENE" = "archaeophyte",
                                                       "ar LM" = "archaeophyte",
                                                       "ar RMP" = "archaeophyte",
                                                       "neo 2" = "neophyte", 
                                                       "neo 1" = "neophyte", 
                                                       "neo 4" = "neophyte",
                                                       "neo 3" = "neophyte",
                                                       "neo *" = "neophyte",
                                                       "ar *" = "archaeophyte", 
                                                       "ar Br" = "archaeophyte",
                                                       "ar BR" = "archaeophyte",
                                                       "ar IR" = "archaeophyte", 
                                                       "ar/neo" = "archaeophyte")) %>% 
           replace_na("native") %>% 
           factor(levels = c("native", "archaeophyte", "neophyte"), ordered = T),
         IUCN = str_replace(IUCN, "LC\\(NA\\)|NA", "NA") %>% replace_na("NA"),
         flower_period_3 = coalesce(flower_period_3, flower_period_4),
         flower_period_7 = coalesce(flower_period_7, flower_period_8),
         flower_period_9 = coalesce(flower_period_9, flower_period_10),
         selfing = coalesce(selfing, cleistogamy, pseudocleistogamy, geitonogamy),
         GR_apomixis = coalesce(GR_apomixis, GR_fac_apomixis, GR_obl_apomixis),
         GR_autogamy = coalesce(GR_autogamy, GR_fac_autogamy),
         GR_allogamy = coalesce(GR_allogamy, GR_fac_allogamy, GR_allogamy_selfincomp), 
         RT_only_veg = coalesce(RT_only_veg, RT_mostly_veg),
         RT_only_seed = coalesce(RT_only_seed, RT_mostly_seed),
         Par_root_hemiparasite = coalesce(Par_root_hemiparasite, Par_stem_hemiparasite, 
                                          Par_root_holoparasite, Par_stem_holoparasite), 
         Par_full_mycoheterotroph = coalesce(Par_part_ini_mycoheterotroph, Par_full_mycoheterotroph),
         symbiosis_rhizobia = coalesce(symbiosis_rhizobia, symbiosis_Frankia),
         myrmecochorous = coalesce(myrmecochorous, myrmecochorous_nv), 
         probably_myrmecochorous = coalesce(probably_myrmecochorous, probably_myrmecochorous_nv), 
         probably_non_myrmecochorous = coalesce(probably_non_myrmecochorous, probably_non_myrmecochorous_nv), 
         non_myrmecochorous_a = coalesce(non_myrmecochorous_a, non_myrmecochorous_a_nv), 
         non_myrmecochorous_b = coalesce(non_myrmecochorous_b, non_myrmecochorous_b_nv), 
         height_mean = (height_min+height_max)/2) %>% 
  rename("Par_parasite" = "Par_root_hemiparasite", "Par_mycoheterotroph" = "Par_full_mycoheterotroph", 
         "flower_period_3_4" = "flower_period_3", "flower_period_7_8" = "flower_period_7", 
         "flower_period_9_10" = "flower_period_9", "symbiosis_N_fixers" = "symbiosis_rhizobia", 
         "RT_seed" = "RT_only_seed", "RT_veg" = "RT_only_veg") %>% 
  select(-c(GR_fac_apomixis, GR_obl_apomixis, GR_fac_autogamy, GR_fac_allogamy, 
            GR_allogamy_selfincomp, GR_sterility, disp_gametophyte, 
            Par_stem_hemiparasite, Par_stem_holoparasite, Par_root_holoparasite,
            Par_part_ini_mycoheterotroph, myrmecochorous_nv, probably_myrmecochorous_nv, 
            probably_non_myrmecochorous_nv, non_myrmecochorous_a_nv, non_myrmecochorous_b_nv, 
            flower_period_4, flower_period_8, flower_period_10, flower_period_1, 
            flower_period_2, flower_period_11, flower_period_12, cleistogamy, pseudocleistogamy, 
            geitonogamy, water_pollination, symbiosis_Frankia, RT_mostly_seed, RT_mostly_veg, 
            ploidy_3, ploidy_5, ploidy_7, ploidy_12:ploidy_72, EIV_salinity, dist_str_based, 
            pladias_id, LF_macrophanerophyte:CSR, leaf_size, SLA, Carnivory, ploidy_2:XEA, 
            eco_spec_nonforest:taxon_weight_ESI_forest, colonization_potential:continentality, 
            "Residence time", Invasion_status, Origin, disp_spore:disp_budding, 
            GF_parasitic_epiphyte, GF_woody_liana, height_min, height_max)) %>% 
  mutate(nr.gf = rowSums(select(., GF_annual:GF_tree) == "TRUE", na.rm = T),
         nr.lls = rowSums(select(., LLS_overwintering_green:LLS_evergreen) == "TRUE", na.rm = T),
         nr.la = rowSums(select(., LA_succulent:LA_hydromorphic) == "TRUE", na.rm = T),
         nr.flower.per = rowSums(select(., flower_period_3_4:flower_period_9_10) == "TRUE", na.rm = T),
         nr.flower.ph = rowSums(select(., flower_phase_1:flower_phase_9) == "TRUE", na.rm = T),
         nr.gr = rowSums(select(., GR_allogamy:GR_apomixis) == "TRUE", na.rm = T),
         nr.poll = rowSums(select(., wind_pollination:selfing) == "TRUE", na.rm = T), 
         nr.rt = rowSums(select(., RT_veg:RT_seed) == "TRUE", na.rm = T), 
         nr.ds = rowSums(select(., DS_Allium:DS_Zea) == "TRUE", na.rm = T),
         nr.myr = rowSums(select(., myrmecochorous:non_myrmecochorous_b) == "TRUE", na.rm = T),
         nr.par = rowSums(select(., c(Par_autotrophic, Par_parasite, Par_mycoheterotroph)) == "TRUE", na.rm = T),
         nr.sym = rowSums(select(., symbiosis_N_fixers:no_nitrogen_fixing) == "TRUE", na.rm = T)) %>% 
  mutate_at(vars(GF_annual:GF_tree), ~ifelse(.x == "TRUE", 1/nr.gf, 0)) %>% 
  mutate_at(vars(LLS_overwintering_green:LLS_evergreen), ~ifelse(.x == "TRUE", 1/nr.lls, 0)) %>% 
  mutate_at(vars(LA_succulent:LA_hydromorphic), ~ifelse(.x == "TRUE", 1/nr.la, 0)) %>% 
  mutate_at(vars(flower_period_3_4:flower_period_9_10), ~ifelse(.x == "TRUE", 1/nr.flower.per, 0)) %>% 
  mutate_at(vars(flower_phase_1:flower_phase_9), ~ifelse(.x == "TRUE", 1/nr.flower.ph, 0)) %>% 
  mutate_at(vars(GR_allogamy:GR_apomixis), ~ifelse(.x == "TRUE", 1/nr.gr, 0)) %>%
  mutate_at(vars(wind_pollination:selfing), ~ifelse(.x == "TRUE", 1/nr.poll, 0)) %>%
  mutate_at(vars(RT_veg:RT_seed), ~ifelse(.x == "TRUE", 1/nr.rt, 0)) %>%
  mutate_at(vars(DS_Allium:DS_Zea), ~ifelse(.x == "TRUE", 1/nr.ds, 0)) %>%
  mutate_at(vars(myrmecochorous:non_myrmecochorous_b), ~ifelse(.x == "TRUE", 1/nr.myr, 0)) %>%
  mutate_at(vars(c(Par_autotrophic, Par_parasite, Par_mycoheterotroph)), ~ifelse(.x == "TRUE", 1/nr.par, 0)) %>%
  mutate_at(vars(symbiosis_N_fixers:no_nitrogen_fixing), ~ifelse(.x == "TRUE", 1/nr.sym, 0)) %>%
  mutate_at(vars(height_mean, dist_freq:dist_sev, LS_Pierce_C:LS_Pierce_R,
                 eco_spec_all:colonization_success, EIV_light:EIV_nutr),
            ~as.numeric(.x)) %>% 
  select(-c(nr.gf:nr.sym)) %>% 
  mutate(dist_freq = ifelse(GF_tree == 0 & GF_shrub == 0, dist_freq, NA), 
         dist_sev = ifelse(GF_tree == 0 & GF_shrub == 0, dist_sev, NA))

# general statements ------------------------------------------------------
# numbers of significantly increasing and decreasing species
table(resdf.clip.join$change.lm)

# plot regression slopes of all significantly increasing and decreasing species
# Appendix B.4
ggplot(resdf.clip.join)+
  theme_classic()+
  geom_bar(aes(x = lm.rank, y = lm.slope, fill = lm.slope), stat = "identity", width = 1)+
  scale_fill_gradient2(name = "",
                       low = "#a50026", mid = "grey85", high = "#006837", breaks = seq(-60, 80, 10),
                       guide = guide_colorbar())+
  geom_hline(yintercept = 0, linewidth = 0.1, col = "grey65")+
  geom_vline(xintercept = table(resdf.clip.join$change.lm)["-1"]+1, linewidth = 0.5, color="black")+
  geom_segment(arrow = arrow(end = "first", type = "closed", length = unit(0.3, "cm")), 
               aes(x = 80, y = 0.002, xend = 180, yend = 0.002), 
               linewidth = 0.5, color = "#a50026", 
               lineend = "round", linejoin = "mitre")+
  geom_segment(arrow = arrow(end = "last", type = "closed", length = unit(0.3,"cm")), 
               aes(x = 700, y = -0.002, xend = 800, yend = -0.002), 
               linewidth = 0.5, color = "#006837", 
               lineend = "round", linejoin = "mitre")+
  geom_text(aes(x = 130, y = 0.003, label = paste0("decreasing (", table(change.lm)["-1"], ")")), size = 5)+
  geom_text(aes(x = 750, y = -0.003, label = paste0("increasing (", table(change.lm)["1"], ")")), size = 5)+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.016, 0.016), breaks = seq(-0.02, 0.02, 0.002))+
  scale_x_continuous(expand = expansion(add = c(1, 1)))+
  guides(x = "none")+
  my.theme+
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        title = element_blank())+
  ylab("Regression slope")

ggsave(filename = paste0(path.plots, "occ_nogroup_lm.png"), width = 10, height = 6)

# top winners and losers --------------------------------------------------
# Fig. 1
# plot top 50 winners and losers with their regression slopes 
# and add colours for Red List and biogeographic status

resdf.clip.join %>% 
  filter(lm.rank %in% c(1:50, (nrow(resdf.clip.join)-49):nrow(resdf.clip.join))) %>%
  mutate(lm.rank = rank(lm.slope), IUCN = fct_relevel(IUCN, c("CR", "EN", "VU", "NT", "LC", 'NA')), 
         specname = str_trim(specname, side = "both")) %>% 
  arrange(lm.rank) %>% 
  mutate(just = c(rep(0.0002, 50), rep(-0.0002, 50))) %>%
  ggplot()+
  theme_classic()+
  my.theme+
  geom_bar(aes(x = lm.rank, y = lm.slope, fill = IUCN), stat = "identity")+
  geom_hline(yintercept = 0, col = "black")+
  geom_text(aes(color = alien, x = lm.rank, y = just, label = specname), 
            hjust = c(rep(0, 50), rep(1, 50)),
            size = 2.1)+
  scale_fill_manual(values = iucn.col, labels = c('Critically Endangered', 
                                                  'Endangered', 'Vulnerable', 
                                                  'Near Threatened', 'Least Concern', 
                                                  'Least Concern <br>(not on the Red List)'))+
  scale_color_manual(values = c("native" = "black", "archaeophyte" = "#3f37c9", 
                                "neophyte" = "#b5179e"), 
                     labels = paste("<span style='color:", 
                                    c('black', '#3f37c9', "#b5179e"), 
                                    "'>", 
                                    c("Native", "Archaeophyte", "Neophyte"), 
                                    "</span>",
                                    sep = ""))+
  coord_flip()+
  scale_y_continuous(expand = c(0,0), limits = c(-0.012, 0.013))+
  scale_x_continuous(expand = c(0.01,0.01))+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), 
        axis.title.y = element_blank(), axis.line.y = element_blank(), 
        legend.position = "bottom", legend.box = "vertical", 
        legend.text = element_markdown(size = 9))+
  annotate("text", y = -0.0085, x = 25, angle = 90, label = "Decreasing")+
  annotate("text", y = 0.0125, x = 75, angle = 270, label = "Increasing")+
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 2), 
         color = guide_legend(title.position = "top", title.hjust = 0.5, 
                              override.aes = list(color = NA)))+
  labs(y = "Regression slope", fill = "Red List category", color = "Biogeographic status")

ggsave(paste0(path.plots,"tb50.png"), width = 10, height = 26, units = c("cm"))

# linear models for continuous variables ----------------------------------

# calculate linear models for species trends with species characteristics as predictors
lm.df <- resdf.clip.join %>% 
  select(lm.slope, change.lm, LS_Pierce_C, LS_Pierce_S, LS_Pierce_R, 
         eco_spec_all, colonization_success, EIV_light:EIV_nutr,
         dist_freq, dist_sev) %>% 
  pivot_longer(cols = c(LS_Pierce_C:dist_sev)) %>% 
  filter(!is.na(value)) %>% 
  nest(data = !name) %>% 
  mutate(lm1 = map(.$data, ~lm(lm.slope~value, data = .x)),
         lm.glance = map(.$data, ~lm(lm.slope~value, data = .x) %>% glance()), 
         lm.augment = map(.$data, ~lm(lm.slope~value, data = .x) %>% augment())) 

# check model diagnostics
lm1.diag <- lm.df$lm1 %>% 
  map(., ~autoplot(.x))

# export diagnostics to pdf
walk2(lm1.diag, lm.df$name, ~grid.arrange(grobs = .x@plots, top = .y) %>% 
        ggsave(paste0(path.plots,"/lm_diagnostics_", .y, ".pdf"), .,
               height = 6, width = 8))

# Appendix B.3 correlation of variables in model
cor.all <- resdf.clip.join %>%  
  select("Competitive strategy" = LS_Pierce_C, "Stress-tolerant strategy" = LS_Pierce_S, 
         "Ruderal strategy" =  LS_Pierce_R, "Ellenberg light" = EIV_light, 
         "Ellenberg temperature" = EIV_temp, 
         "Ellenberg moisture" = EIV_moist, "Ellenberg soil reaction" = EIV_react, 
         "Ellenberg nutrients" = EIV_nutr, "Disturbance frequency" = dist_freq, 
         "Disturbance severity" = dist_sev, "Colonization success" = colonization_success, 
         "Ecological specialization" = eco_spec_all) %>% 
  apply(2, as.numeric) %>% 
  cor(method = "spearman", use = "pairwise.complete.obs")

png(filename = paste0(path.plots, "corrplot.png"), width = 8.5, height = 6, units = "cm", res = 1024)
corrplot(cor.all, method = "number", tl.cex = 0.3, number.cex = 0.22, cl.cex = 0.3, 
         type = "upper", diag = F, tl.col = "black", tl.srt = 45)
dev.off()

# Fig. 2
resdf.clip.join %>% 
  select(lm.slope, change.lm, LS_Pierce_C, LS_Pierce_S, LS_Pierce_R, 
         eco_spec_all, colonization_success, EIV_light:EIV_nutr,
         dist_freq, dist_sev) %>% 
  pivot_longer(cols = c(LS_Pierce_C:dist_sev)) %>% 
  mutate(name = fct_relevel(name, c("LS_Pierce_C", "LS_Pierce_S", "LS_Pierce_R", 
                                    "EIV_light", "EIV_temp", "EIV_moist", "EIV_react", 
                                    "EIV_nutr", "dist_freq", "dist_sev", "colonization_success", 
                                    "eco_spec_all"))) %>%
  ggplot(aes(value, lm.slope)) +
  geom_hline(yintercept = 0, color = "gray4", linewidth= 0.5)+
  geom_point(aes(fill = as.factor(change.lm)), pch = 21) +
  scale_fill_discrete(labels = c("decreasing", "increasing")) +
  geom_smooth(formula = y~x, method = "lm", color = 1) +
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = "bottom", axis.title.x = element_blank())+
  labs(y = "Regression slope") +
  stat_poly_eq(formula = y~x, use_label(c("eq", "adj.R2", "P")), label.x = "right",  
               label.y = 0.01, size = 2.5, coef.digits = 2)+
  facet_wrap(~name, scales = "free_x", 
             labeller = labeller(name = c(colonization_success = "(K) Colonization success", 
                                          dist_freq = "(I) Disturbance frequency", 
                                          dist_sev = "(J) Disturbance severity", 
                                          eco_spec_all = "(L) Ecological specialization", 
                                          EIV_light = "(D) Ellenberg light", 
                                          EIV_moist = "(F) Ellenberg moisture", 
                                          EIV_nutr = "(H) Ellenberg nutrients", 
                                          EIV_react = "(G) Ellenberg soil reaction", 
                                          EIV_temp = "(E) Ellenberg temperature", 
                                          LS_Pierce_C = "(A) Competitive strategy", 
                                          LS_Pierce_R  = "(C) Ruderal strategy", LS_Pierce_S = "(B) Stress-tolerant strategy")))

ggsave(paste0(path.plots, "linear_models.png"), height = 8, width = 11)
ggsave(paste0(path.plots, "linear_models.svg"), height = 8, width = 11)

# full multipredictor model --------------------------------------------------------------
model.df <- resdf.clip.join %>% 
  select(lm.slope, EIV_moist, EIV_nutr, eco_spec_all, dist_freq, dist_sev, 
         colonization_success, LS_Pierce_C, LS_Pierce_R, LS_Pierce_S) %>% 
  mutate(across(c(EIV_moist:LS_Pierce_S), ~scale(.) %>% as.vector()))

full.mod <- lm(lm.slope~EIV_moist + EIV_nutr + dist_freq + dist_sev + eco_spec_all +
                 colonization_success + LS_Pierce_S + LS_Pierce_C, data = model.df)

summary(full.mod)
drop1(full.mod)

# backward selection to the final model
full.mod <- update(full.mod, .~.-LS_Pierce_C)
full.mod <- update(full.mod, .~.-LS_Pierce_S)
full.mod <- update(full.mod, .~.-dist_freq)

autoplot(full.mod) # check model diagnostics

# chi-square tests and annotation plots for categorical variables -------------

### Fig. 3
# Red List status
chi.iucn <- chi.test2(col = "IUCN", T, c("DD"))
chi.iucn
0.05/((ncol(chi.iucn))*2) # significance level after Bonferroni correction

resdf.clip.join <- resdf.clip.join %>% 
  arrange(IUCN)

iucn.plot <- annot.plot.1col(df = resdf.clip.join,
                             iucn.col, 
                             colname = "IUCN", 
                             title.left = "(A) Red List category", 
                             size.atext = t.size, 
                             labdist = 50, 
                             label = c("Critically Endangered", "Endangered", 
                                       "Vulnerable", "Near Threatened", 
                                       "Least Concern", 'Least Concern \n(not on the Red List)'), 
                             annot.exp = 0, 
                             exp.x = 120)

# biogeographic status
chi.al <- chi.test2(col = "alien", F, "")
chi.al
0.05/((ncol(chi.al))*2) 

biogeo.plot <- annot.plot.1col(df = resdf.clip.join,
                               biogeo.col, 
                               colname = "alien", 
                               title.left = "(B) Biogeographic status", 
                               size.atext = t.size, 
                               labdist = 50, 
                               label = c("Native", "Archaeophyte", "Neophyte"), 
                               exp.x = 80)

# growth forms
chi.gf <- chi.test("GF_annual", "GF_tree", exclude = F, "", sim = T)
chi.gf
0.05/((ncol(chi.gf))*2) # correction for multiple testing

gf.ramp<-colorRampPalette(c("#196F3D", "#8B7B0F","#873600"))
gf.col<-c(gf.ramp(9))

gf.plot <- annot.plot(col1 = "GF_annual", col2 = "GF_tree", 
                      colvec = gf.col, size.atext = 3, labdist = 40, 
                      title.left = "(C) Growth form", my.theme = my.theme, 
                      label = c("Annual herb", "Monocarpic herb", "Polycarpic herb",
                                "Clonal herb", "Dwarf shrub", "Shrub", "Tree"), 
                      exp.x = 80) 

# Fig. 3
plot1 <- iucn.plot + biogeo.plot + gf.plot 
ggsave(paste0(path.plots,"iucn_biogeo_gf_specnumber.png"), plot1, width = 13, height = 3.8)

# export to svg
svglite(paste0(path.plots,"iucn_biogeo_gf_specnumber.svg"), width = 13, height = 3.8)
plot1
dev.off()

### Appendix B.5
# leaf life span
chi.lls <- chi.test(col1 = "LLS_overwintering_green", col2 = "LLS_evergreen", F, "", sim = F)
chi.lls
0.05/((ncol(chi.lls))*2)

lls.ramp <- colorRampPalette(c("#004b23", "#b5e48c"))
lls.col <- c(lls.ramp(4))

lls.plot <- annot.plot(col1 = "LLS_overwintering_green", col2 = "LLS_evergreen", 
                       colvec = lls.col, size.atext = 3, labdist = 40, 
                       title.left = "(A) Leaf life span", my.theme = my.theme, 
                       label = c("Overwintering green", "Spring green", "Summer green", "Evergreen")) 


# leaf anatomy
chi.la <- chi.test.sim(col1 = "LA_succulent", col2 = "LA_hydromorphic", F, "")
chi.la
0.05/((ncol(chi.la))*2)

lls.ramp <- colorRampPalette(c("#004b23", "#b5e48c"))
la.col <- c(lls.ramp(6))

la.plot <- annot.plot(col1 = "LA_succulent", col2 = "LA_hydromorphic", 
                      colvec = la.col, size.atext = 3, labdist = 40, 
                      title.left = "(B) Leaf anatomy", my.theme = my.theme, 
                      label = c("Succulent", "Scleromorphic", "Mesomorphic",
                                "Hygromorphic", "Helomorphic", "Hydromorphic")) 

# flower period
chi.fl <- chi.test(col1 = "flower_period_3_4", col2 = "flower_period_9_10", F, "", F)
chi.fl
0.05/((ncol(chi.fl))*2)

flower.ramp <- colorRampPalette(c("#A2D2FF", "#FFAFCC"))
flower.col <- c(flower.ramp(5))

fl.per.plot <- annot.plot1(col1 = "flower_period_3_4", col2 = "flower_period_9_10", 
                           colvec = flower.col, size.atext = 3, labdist = 40, 
                           title.left = "(C) Flower period", my.theme = my.theme)

# flower phase
chi.fl <- chi.test(col1 = "flower_phase_1", col2 = "flower_phase_9", F, "flower_phase_10", T)
chi.fl
0.05/((ncol(chi.fl))*2)

flower.ramp <- colorRampPalette(c("#A2D2FF", "#FFAFCC"))
flower.col <- c(flower.ramp(10))

fl.ph.plot <- annot.plot1(col1 = "flower_phase_1", col2 = "flower_phase_9", 
                          colvec = flower.col, size.atext = 3, labdist = 40, 
                          title.left = "(D) Flower phase", my.theme = my.theme)

# generative reproduction type
chi.gr <- chi.test(col1 = "GR_allogamy", col2 = "GR_apomixis", F, "", F)
chi.gr
0.05/((ncol(chi.gr))*2)

gr.ramp <- colorRampPalette(c("#FFAFCC", "#A2D2FF"))
gr.col <- c(gr.ramp(4))

gr.plot <- annot.plot(col1 = "GR_allogamy", col2 = "GR_apomixis", 
                      colvec = gr.col, size.atext = 3, labdist = 40, 
                      title.left = "(E) Generative reproduction type", my.theme = my.theme, 
                      label = c("Allogamy", "Autogamy", "Mixed mating", "Apomixis"))

# pollination 
chi.pol <- chi.test(col1 = "wind_pollination", col2 = "selfing", F, "", F)
chi.pol
0.05/((ncol(chi.pol))*2)

flower.ramp <- colorRampPalette(c("#A2D2FF", "#FFAFCC"))
flower.col <- c(flower.ramp(3))

poll.plot <- annot.plot(col1 = "wind_pollination", col2 = "selfing", 
                        colvec = flower.col, size.atext = 3, labdist = 40, 
                        title.left = "(F) Pollination syndrome", my.theme = my.theme, 
                        label = c("Wind pollination", "Insect pollination", "Selfing")) 


# reproduction type
chi.rt <- chi.test(col1 = "RT_veg", col2 = "RT_seed", F, "", F)
chi.rt
0.05/((ncol(chi.rt))*2)

rt.ramp <- colorRampPalette(c("#65C18C", "#eae2b7", "#7f5539"))
rt.col <- c(rt.ramp(3))

rt.plot <- annot.plot(col1 = "RT_veg", col2 = "RT_seed", 
                      colvec = rt.col, size.atext = 3, labdist = 40, 
                      title.left = "(G) Reproduction type", my.theme = my.theme, 
                      label = c("Vegetative", "Seed or vegetative", "Seed")) 

### Dispersal strategy
chi.ds <- chi.test(col1 = "DS_Allium", col2 = "DS_Zea", F, "", T)
chi.ds
0.05/((ncol(chi.ds))*2)

seed.ramp <- colorRampPalette(c("#eae2b7", "#b08968","#7f5539"))
seed.col <- c(seed.ramp(9))

ds.plot <- annot.plot(col1 = "DS_Allium", col2 = "DS_Zea", 
                      colvec = seed.col, size.atext = 3, labdist = 40, 
                      title.left = "(H) Seed dispersal strategy", my.theme = my.theme, 
                      label = c("Allium", "Bidens", "Cornus", "Epilobium", "Lycopodium", 
                                "Phragmites", "Sparganium", "Wolffia", "Zea")) 

### Myrmecochory
chi.myr <- chi.test(col1 = "myrmecochorous", col2 = "non_myrmecochorous_b", F, "", F)
chi.myr
0.05/((ncol(chi.myr))*2)

seed.ramp<-colorRampPalette(c("#eae2b7", "#b08968","#7f5539"))
seed.col<-c(seed.ramp(5))

myr.plot <- annot.plot(col1 = "myrmecochorous", col2 = "non_myrmecochorous_b", 
                       colvec = seed.col, size.atext = 3, labdist = 40, 
                       title.left = "(I) Myrmecochory", my.theme = my.theme, 
                       label = c("Myrmecochorous", "Probably myrmecochorous", 
                                 "Probably non-myrmecochorous", "Non-myrmecochorous a", 
                                 "Non-myrmecochorous b")) 

# nitrogen fixation
chi.n <- chi.test(col1 = "symbiosis_N_fixers", col2 = "no_nitrogen_fixing", F, "", F)
chi.n
0.05/((ncol(chi.n))*2)

n.ramp<-colorRampPalette(c("#FAAB78", "#65C18C"))
n.col<-c(n.ramp(2))

n.plot <- annot.plot(col1 = "symbiosis_N_fixers", col2 = "no_nitrogen_fixing", 
                     colvec = n.col, size.atext = 3, labdist = 40, 
                     title.left = "(J) Nitrogen fixation", my.theme = my.theme, 
                     label = c("Symbiotic N fixers", "No nitrogen fixing")) 

# autotrophic/parasite
chi.au <- chi.test(col1 = "Par_autotrophic", col2 = "Par_mycoheterotroph", F, "", F)
chi.au
0.05/((ncol(chi.au))*2)

par.ramp<-colorRampPalette(c("#65C18C", "#eae2b7", "#7f5539"))
par.col<-c(par.ramp(3))

par.plot <- annot.plot(col1 = "Par_autotrophic", col2 = "Par_mycoheterotroph", 
                       colvec = par.col, size.atext = 3, labdist = 40, 
                       title.left = "(K) Trophic mode", my.theme = my.theme, 
                       label = c("Autotrophic", "Parasite", "Mycoheterotroph")) 



# plots for appendix S8
plot2 <- lls.plot + la.plot + fl.per.plot + fl.ph.plot + gr.plot + poll.plot +
  rt.plot + ds.plot + myr.plot + n.plot + par.plot + plot_layout(ncol = 3) 

ggsave(paste0(path.plots, "appendix_specnumber.png"), plot2, width = 13, height = 13)

# export to svg
svglite(paste0(path.plots, "appendix_specnumber.svg"), width = 13, height = 13)
plot2
dev.off()


# clustering of temporal trends -------------------------------------------

# logit transformation of occupancy estimates
clust.df <- occ.res %>% 
  semi_join(resdf.clip.join) %>%  
  select(specname, half_dec, median) %>% 
  mutate(median = logit(median)) 

# to wide format for time-series clustering algorithm
clust.df2 <- clust.df %>% 
  pivot_wider(names_from = "half_dec", values_from = "median") %>% 
  column_to_rownames("specname") 

# standardize occupancy estimates
clust.df.scaled <- ddply(clust.df,.(specname), function(x){
  x$median_scaled <- as.numeric(zscore(x$median))
  return(x)
})

# comparison of different clusterings - modified k (number of clusters)
p_cfgs <- compare_clusterings_configs(
  types = "p", k = 5L,
  controls = list(
    partitional = partitional_control(
      iter.max = 20L,
      nrep = 50L
    )
  ),
  preprocs = pdc_configs(
    "preproc",
    zscore = list(center = c(TRUE))
  ),
  distances = pdc_configs(
    "distance",
    sbd = list()
  ),
  centroids = pdc_configs(
    "centroid",
    partitional = list(
      shape = list()
    )
  )
)

comparison_partitional <- compare_clusterings(clust.df2, types = "p",
                                              configs = p_cfgs,
                                              seed = 32903L, trace = TRUE,
                                              score.clus = score_fun,
                                              pick.clus = pick_fun,
                                              shuffle.configs = TRUE,
                                              return.objects = TRUE)

# final number of 5 clusters
comp5 <- comparison_partitional

# comparison of cluster validity indices
comp5.cvi <- comp5$objects.partitional %>% 
  map(~cvi(.x)) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rownames_to_column("index") %>% 
  rename(value = ".") %>% 
  mutate(name = str_split_i(index, "\\.", 2), iter = str_split_i(index, "\\.", 1)) %>% 
  select(-index) %>% 
  pivot_wider() %>% 
  mutate(rank_sil = rank(Sil), rank_D = rank(D), rank_COP = rank(COP), 
         rank_sum = (rank_sil + rank_D + rank_COP))

### select the final clustering
clust <- comp5$objects.partitional$config1_16

# join to the data
resdf.clip.join2 <- resdf.clip.join %>% 
  bind_cols(cluster = clust@cluster)

# number of species in each group
resdf.clip.join2 %>% 
  group_by(cluster) %>% 
  summarize(n = n())

# export clustering results to csv
resdf.clip.join2 %>% 
  select(specname, cluster) %>% 
  write_csv("clustering_results.csv")

# loading clustering results
clustering <- read_csv('results/clustering_results.csv')

resdf.clip.join2 <- resdf.clip.join %>%
  left_join(clustering)

# join occupancy trends with clutering results
spec.trend <- occ.res %>%  
  right_join(resdf.clip.join2) %>% 
  mutate(cluster = fct_relevel(as.factor(cluster), c("1", "2", "5", "3", "4"))) %>% 
  select(specname, half_dec, median, cluster)

# Fig. 4 boxplots
clust.box <- spec.trend %>% 
  ggplot(aes(as.factor(half_dec), median, fill = as.factor(cluster)))+
  geom_boxplot()+
  scale_fill_manual(values = clust.col)+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 35, vjust = 0.5))+
  labs(y = "Median occupancy")+ 
  scale_x_discrete(breaks = seq(1960, 2015, 5), 
                   labels = paste0(seq(1960, 2015, 5), "s"))+
  facet_col(~cluster, labeller = labeller(cluster = c("1" = "Increasing Early (n = 249)", 
                                                      "2" = "Increasing Middle (n = 112)", 
                                                      "5" = "Increasing Recently (n = 179)",
                                                      "3" = "Increasing Continuously (n = 214)", 
                                                      "4" = "Decreasing (n = 320)")))
# Fig. 4 line plots 
clust.line <- clust.df.scaled %>% 
  right_join(spec.trend %>% select(-median)) %>% 
  select(half_dec, median_scaled, specname, cluster) %>% 
  ggplot(aes(half_dec, median_scaled, color = specname))+
  geom_line() +
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.position = "none", 
        axis.text.x = element_text(angle = 35, vjust = 0.5), panel.grid.minor.x = element_blank())+
  labs(y = "Scaled median occupancy")+ 
  scale_x_continuous(breaks = seq(1960, 2015, 5), labels = paste0(seq(1960, 2015, 5), "s"))+
  facet_col(~cluster, labeller = labeller(cluster = c("1" = "Increasing Early (n = 249)", 
                                                      "2" = "Increasing Middle (n = 112)", 
                                                      "5" = "Increasing Recently (n = 179)",
                                                      "3" = "Increasing Continuously (n = 214)", 
                                                      "4" = "Decreasing (n = 320)")))

# Fig. 4
clust.box + clust.line + plot_annotation(tag_levels = "A")
ggsave(paste0(path.plots, "clustering2.png"), width = 8, height = 12)

# select variables for Anova tests of diffrences between groups of species with similar trends
plot.df.anova <- resdf.clip.join2 %>% 
  mutate(cluster = fct_relevel(as.factor(cluster), c("1", "2", "5", "3", "4"))) %>% 
  select(specname, cluster, LS_Pierce_C, LS_Pierce_S, LS_Pierce_R, 
         eco_spec_all, colonization_success, 
         EIV_light:EIV_nutr, dist_freq, dist_sev) %>% 
  pivot_longer(cols = c(LS_Pierce_C:dist_sev)) %>% 
  mutate(name = fct_relevel(name, c("LS_Pierce_C", "LS_Pierce_S", "LS_Pierce_R", 
                                    "EIV_light", "EIV_temp", "EIV_moist", "EIV_react", 
                                    "EIV_nutr", "dist_freq", "dist_sev", "colonization_success", 
                                    "eco_spec_all"))) 

# Anova + post hoc tests
plot.df1 <- plot.df.anova %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(anova.m = map(data, ~aov(value~cluster, data = .x)), 
         tidy.anova = map(anova.m, ~tidy(.x)), 
         p.anova = map(tidy.anova, ~.x$p.value[[1]]), 
         tukey = map(anova.m, ~TukeyHSD(.x)[[1]]), 
         letters = map(tukey, ~multcompLetters(.x[,4])$Letters), 
         letters2 = map2(data, tukey, ~multcompLetters2(value~cluster, .y[,"p adj"], data = .x)$Letters),
         placement  = map(data, ~max(.x$value, na.rm = T))) %>% 
  select(name, letters2, placement) %>% 
  unnest(cols = c(letters2)) %>% 
  mutate(cluster = c(1, 2, 5, 3, 4) %>% as.factor)

# check model diagnostics
plot.df1$anova.m %>% 
  map(., ~autoplot(.x))

# not normal distribution -> Kruskal-Wallis test
plot.df2 <- plot.df.anova %>% 
  group_by(name) %>% 
  filter(name %in% c("colonization_success", "LS_Pierce_S", "LS_Pierce_C")) %>% 
  nest() %>% 
  mutate(anova.m = map(data, ~kruskal.test(value~cluster, data = .x)), 
         tidy.anova = map(anova.m, ~tidy(.x)), 
         p.anova = map(tidy.anova, ~.x$p.value[[1]]), 
         tukey = map(data, ~dunn_test(value~cluster, data = .x), 
                     letters2 = map2(data, tukey, ~multcompLetters2(value~cluster, .y %>% pull(p.adj) %>% 
                                                                      setNames(str_c(.y$group1, .y$group2)), 
                                                                    data = .x))))

# replace anova and tukey results with Kruskal-Wallis (manually)
df.test <- plot.df2[[6]][[1]] %>% pull(p.adj) %>% 
  setNames(str_c(plot.df2[[6]][[1]]$group1, plot.df2[[6]][[1]]$group2, sep = "-"))
multcompLetters2(value~cluster, df.test, data = plot.df2$data[[1]])

plot.df.LSC <- plot.df1 %>% 
  filter(name == "LS_Pierce_C") %>% 
  mutate(letters2 = c("ab", "c", "ac", "b", "c"))

plot.df.LSS <- plot.df1 %>% 
  filter(name == "LS_Pierce_S") %>% 
  mutate(letters2 = c("ab", "a", "b", "ab", "ab"))

plot.df.col <- plot.df1 %>% 
  filter(name == "colonization_success") %>% 
  mutate(letters2 = c("a", "b", "b", "c", "d"))

plot.df1 <- plot.df1 %>% 
  filter(!name %in% c("colonization_success", "LS_Pierce_S", "LS_Pierce_C")) %>% 
  bind_rows(plot.df.col, plot.df.LSS, plot.df.LSC)

# Fig. 5
plot.df.anova %>%  
  ggplot(aes(cluster, value, fill = cluster, color = cluster))+
  geom_violin(aes(color = as.factor(cluster)), alpha = 0.2)+
  geom_boxplot(width = 0.15, alpha = 0.5, coef = 0, outlier.shape = NA)+
  geom_text(data = plot.df1, aes(x = as.factor(cluster), y = 1.2*as.numeric(placement), 
                                 label = letters2), size = 3, show.legend = F)+
  scale_fill_manual(values = clust.col, labels = c("Increasing Early", "Increasing Middle", 
                                                   "Increasing Recently", "Increasing Continuously", 
                                                   "Decreasing"), aesthetics = c("fill", "color"))+
  theme_bw()+
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))+
  theme(axis.title = element_blank(), legend.position = "bottom", axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), legend.title = element_blank())+
  facet_wrap(~name, scales = "free", labeller = labeller(name = c(colonization_success = "(K) Colonization success", 
                                                                  dist_freq = "(I) Disturbance frequency", 
                                                                  dist_sev = "(J) Disturbance severity", 
                                                                  eco_spec_all = "(L) Ecological specialization", 
                                                                  EIV_light = "(D) Ellenberg light", 
                                                                  EIV_moist = "(F) Ellenberg moisture", 
                                                                  EIV_nutr = "(H) Ellenberg nutrients", 
                                                                  EIV_react = "(G) Ellenberg soil reaction", 
                                                                  EIV_temp = "(E) Ellenberg temperature", 
                                                                  LS_Pierce_C = "(A) Competitive strategy", 
                                                                  LS_Pierce_R  = "(C) Ruderal strategy", 
                                                                  LS_Pierce_S = "(B) Stress-tolerant strategy")))

ggsave(paste0(path.plots, "continuous_5clust_sbd_shape2.png"), height = 6, width = 8) 

# differences in categorical characteristics between clusters ----------------------
## proportion of alien species
chi.df <- chisq.test(resdf.clip.join2$cluster, resdf.clip.join2$alien)
chi.df$expected
chi.df$observed

2*(1-pnorm(abs(chi.df$residuals))) %>% round(5) 
0.05/(3*5)

# proportion of red-list species
resdf.clip.join2$IUCN[resdf.clip.join2$IUCN == "DD"] <- NA
chi.df <- chisq.test(resdf.clip.join2$cluster, resdf.clip.join2$IUCN, simulate.p.value = T)
chi.df$expected
chi.df$observed

2*(1-pnorm(abs(chi.df$residuals))) %>% round(5) 
0.05/(6*5)

# growth forms
chi.gf <- chi.test.clust("GF_annual", "GF_tree")
0.05/(7*5)

# leaf life span
chi.lls <- chi.test.clust(col1 = "LLS_overwintering_green", col2 = "LLS_evergreen")
0.05/((ncol(chi.lls))*5)

# leaf anatomy
chi.la <- chi.test.clust(col1 = "LA_succulent", col2 = "LA_hydromorphic")
0.05/((ncol(chi.la))*5)

# flower period
chi.fl <- chi.test.clust(col1 = "flower_period_3_4", col2 = "flower_period_9_10")
0.05/((ncol(chi.fl))*5)

# flower phase
chi.fl <- chi.test.clust(col1 = "flower_phase_1", col2 = "flower_phase_9")
0.05/((ncol(chi.fl))*5)

# generative reproduction type
chi.gr <- chi.test.clust(col1 = "GR_allogamy", col2 = "GR_apomixis")
0.05/((ncol(chi.gr))*5)

# pollination 
chi.pol <- chi.test.clust(col1 = "wind_pollination", col2 = "selfing")
0.05/((ncol(chi.pol))*5)

# reproduction type
chi.rt <- chi.test.clust(col1 = "RT_veg", col2 = "RT_seed")
0.05/((ncol(chi.rt))*2)

### Dispersal strategy
chi.ds <- chi.test.clust(col1 = "DS_Allium", col2 = "DS_Zea")
0.05/((ncol(chi.ds))*2)

### Myrmecochory
chi.myr <- chi.test.clust(col1 = "myrmecochorous", col2 = "non_myrmecochorous_b")
0.05/((ncol(chi.myr))*2)

# autotrophic/parasite
chi.au <- chi.test.clust(col1 = "Par_autotrophic", col2 = "Par_mycoheterotroph")
0.05/((ncol(chi.au))*5)

# nitrogen fixation
chi.n <- chi.test.clust(col1 = "symbiosis_N_fixers", col2 = "no_nitrogen_fixing")
0.05/((ncol(chi.n))*5)

# export table with species trends Appendix A.3 ------------------------
occ.res %>% 
  left_join(read_csv("results/records_decades.csv") %>% 
              filter(pres != 0) %>% 
              group_by(name_lat) %>% 
              summarize(clip.min = min(half_dec), clip.max = max(half_dec)), 
            by = c('specname' = 'name_lat')) %>% 
  filter(parnum >= clip.min & parnum <= clip.max & !specname %in% n.conv.full$specname) %>% 
  group_by(specname) %>% 
  nest() %>%
  mutate(lm = map(data, ~lm(median~parnum, data = .x)), 
         lm.tidy = map(lm, ~tidy(.x)), 
         lm.slope = map_dbl(lm.tidy, ~.x$estimate[[2]]), 
         lm.p = map_dbl(lm.tidy, ~.x$p.value[[2]])) %>% 
  select(specname, lm.slope, lm.p) %>%
  ungroup() %>% 
  mutate(change.lm = sign(lm.slope), 
         lm.rank = rank(lm.slope, ties.method = "random")) %>%
  left_join(clustering) %>% 
  select(specname, lm.slope, lm.p, change.lm, cluster) %>% 
  mutate(cluster = fct_recode(as.factor(cluster), "Increasing Early" = "1", 
                              "Increasing Middle" = "2", 
                              "Increasing Recently" = "5",
                              "Increasing Continuously" = "3", 
                              "Decreasing" = "4"), 
         change.lm = ifelse(lm.p < 0.05, change.lm, NA)) %>% 
  write_csv("change_traits.csv")
