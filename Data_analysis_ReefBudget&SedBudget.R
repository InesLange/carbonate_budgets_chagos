####
## data analysis for:
##
## Carbonate framework and sediment production across island-fringing coral reef 
## habitats and a natural nutrient gradient
##
## Lange ID, Stuhr M, Perry CT, Gea Neuhaus A 
## published in ... (2026)

#
# Author Ines Lange
#
# This code compares reef framework budgets and reef sediment budgets 
# as well as benthic community composition and sediment composition for reefs 
# around islands with high and low sebird densities in the Chagos Archipelago, BIOT.
# Data were collected on lagoonal reefs (1-2m), shallow forereefs (2-3m) and deep forereefs (8-9m)   
# in Apr/May 2021 (Sediment samples) and Oct/Nov 2023 (Sediment samples, ReefBudget surveys, SedBudget surveys) 
# by Ines Lange and Marleen Stuhr 

## load packages

# organise and plot data
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot) #multipanel figures

# Baysian analysis
library(Rcpp)
library(brms)
library(tidybayes)
library(bayesplot)
library(emmeans)

# PCA
library(vegan)
library(FactoMineR)
library(factoextra)
library(pairwiseAdonis)


#####
## load data
#####

# reef framework budget data, coral cover, rugosity 
reefbudget <- read.csv("data/RB_summary_2023.csv")%>%
  mutate(status = as.factor(status),
         habitat = as.factor(habitat))

reefbudget <- reefbudget %>%
  mutate(atoll = recode(atoll,
                        "SA" = "Salomon",
                        "PB" = "Peros Banhos",
                        "GCB" = "Great Chagos Bank"))
reefbudget <- reefbudget %>%
  mutate(status = recode(status,
                         "birdy" = "high",
                         "ratty" = "low"))

reefbudget$status <- factor(reefbudget$status, 
                            levels = c("low","high"))
reefbudget$island <- factor(reefbudget$island, 
                            levels = c("Ile de la Passe", "Ile Anglaise", 
                                       "Grande Coquillage","Ile Poule",
                                       "Middle Brother","Eagle"))
reefbudget$atoll <- factor(reefbudget$atoll, 
                           levels = c("Salomon", "Peros Banhos", "Great Chagos Bank"))
reefbudget$habitat <- factor(reefbudget$habitat, 
                             levels = c("lagoon", "shallow", "deep"))


# sediment budget data, parrotfish biomass
sedbudget <- read.csv("data/SB_summary_2023.csv",check.names = FALSE)%>%
  mutate(status = as.factor(status),
         habitat = as.factor(habitat))

sedbudget <- sedbudget %>%
  mutate(atoll = recode(atoll,
                        "SA" = "Salomon",
                        "PB" = "Peros Banhos",
                        "GCB" = "Great Chagos Bank"))
sedbudget <- sedbudget %>%
  mutate(status = recode(status,
                         "birdy" = "high",
                         "ratty" = "low"))

sedbudget$status <- factor(sedbudget$status, 
                           levels = c("low","high"))
sedbudget$island <- factor(sedbudget$island, 
                           levels = c("Ile de la Passe", "Ile Anglaise", 
                                      "Grande Coquillage","Ile Poule",
                                      "Middle Brother","Eagle"))
sedbudget$atoll <- factor(sedbudget$atoll, 
                          levels = c("Salomon", "Peros Banhos", "Great Chagos Bank"))
sedbudget$habitat <- factor(sedbudget$habitat, 
                            levels = c("lagoon", "shallow", "deep"))
sedbudget$status <- factor(sedbudget$status, 
                           levels = c("low","high"))

sedbudget <- subset(sedbudget, parrotfish <= 12) # exclude one outlier transect with sedprod>12 (Anglaise shallow T3)


# sediment composition data
sedcomp <- read.csv("data/SB_sedimentcomp_2023.csv", check.names = FALSE)%>%
  mutate(status = as.factor(status),
         habitat = as.factor(habitat))

sedcomp <- sedcomp %>%
  mutate(atoll = recode(atoll,
                        "SA" = "Salomon",
                        "PB" = "Peros Banhos",
                        "GCB" = "Great Chagos Bank"))

sedcomp$status <- factor(sedcomp$status, 
                         levels = c("low","high"))
sedcomp$island <- factor(sedcomp$island, 
                         levels = c("Ile de la Passe", "Ile Anglaise", 
                                    "Grande Coquillage","Ile Poule",
                                    "Middle Brother","Eagle"))
sedcomp$atoll <- factor(sedcomp$atoll, 
                        levels = c("Salomon", "Peros Banhos", "Great Chagos Bank"))
sedcomp$habitat <- factor(sedcomp$habitat, 
                          levels = c("beach", "water","lagoon", "shallow", "deep"))
sedcomp$status <- factor(sedcomp$status, 
                         levels = c("low","high"))

sedcomp_y <- subset(sedcomp, comparison == "y") # only sites comparable to transect locations 
sedcomp_beach <- subset(sedcomp, habitat == "beach") # only sites at beach


# benthic community composition
benthic <- read.csv("data/RB_benthic_2023.csv")%>%
  mutate(status = as.factor(status),
         habitat = as.factor(habitat),
         atoll = as.factor(atoll))

benthic <- benthic %>%
mutate(atoll = recode(atoll,
                      "SA" = "Salomon",
                      "PB" = "Peros Banhos",
                      "GCB" = "Great Chagos Bank"))

benthic$habitat <- factor(benthic$habitat, 
                          levels = c("lagoon", "shallow", "deep"))

#####
## summarise budget data at site level and by status
#####

#carbonate budget
avg_reefbud_island <- reefbudget %>%
  group_by(island, habitat, status) %>%
  summarise_at(vars(5:19),list(mean=mean,sd=sd))
(avg_reefbud_island <- as.data.frame(avg_reefbud_island))
write.csv(avg_reefbud_island, "avg_reefbud_island.csv", row.names=F)

avg_reefbud_status <- avg_reefbud_island %>%
  group_by(habitat,status) %>%
  summarise_at(vars(2:16),list(mean=mean,sd=sd))
(avg_reefbud_status <- as.data.frame(avg_reefbud_status))
names(avg_reefbud_status) <- str_replace(
  names(avg_reefbud_status), "_mean_", "_") #remove middle mean from names
write.csv(avg_reefbud_status, "avg_reefbud_status.csv", row.names=F)

#sediment production
avg_sedprod_island <- sedbudget %>%
  group_by(island,habitat,status) %>%
  summarise_at(vars(5:13),list(mean=mean,sd=sd))
(avg_sedprod_island <- as.data.frame(avg_sedprod_island))
write.csv(avg_sedprod_island, "avg_sedprod_island.csv", row.names=F)

avg_sedprod_status <- avg_sedprod_island %>%
  group_by(habitat,status) %>%
  summarise_at(vars(2:10),list(mean=mean,sd=sd))
(avg_sedprod_status <- as.data.frame(avg_sedprod_status))
names(avg_sedprod_status) <- str_replace(
  names(avg_sedprod_status), "_mean_", "_") #remove middle mean from names
write.csv(avg_sedprod_status, "avg_sedprod_status.csv", row.names=F)

avg_sedprod_ratio_island <- avg_sedprod_island %>%
  mutate(across(c(urchins_mean, parrotfish_mean, molluscs_mean, algae_mean, forams_mean, sponge_mean, coral_mean, cca_mean),
                ~ . / sediment_mean,
                .names = "{.col}_ratio")) %>%
  select(habitat, island, ends_with("_ratio"))
names(avg_sedprod_ratio_island) <- str_replace(
  names(avg_sedprod_ratio_island), "_mean_", "_") #remove middle mean from names
write.csv(avg_sedprod_ratio_island, "avg_sedprod_ratio_island.csv", row.names=F)

avg_sedprod_ratio_status <- avg_sedprod_status %>%
  mutate(across(c(urchins_mean, parrotfish_mean, molluscs_mean, algae_mean, forams_mean, sponge_mean,coral_mean, cca_mean),
                ~ . / sediment_mean,
                .names = "{.col}_ratio")) %>%
  select(habitat, status, ends_with("_ratio"))
names(avg_sedprod_ratio_status) <- str_replace(
  names(avg_sedprod_ratio_status), "_mean_", "_") #remove middle mean from names
write.csv(avg_sedprod_ratio_status, "avg_sedprod_ratio_status.csv", row.names=F)

## parrotfish biomass
avg_parrot_biomass_island <- sedbudget %>%
  group_by(island, habitat, status) %>%
  summarise(
    mean_parrot_biomass = mean(parrotfish_biomass, na.rm = TRUE),
    se_parrot_biomass   = sd(parrotfish_biomass, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(avg_parrot_biomass_island, aes(x = habitat, y = mean_parrot_biomass, color = status)) +
  #geom_point(position = position_jitter(width = 0.1)) +
  geom_boxplot()+
  theme_minimal()

avg_parrot_biomass_status <- avg_parrot_biomass_island %>%
  group_by(habitat, status) %>%
  summarise(
    mean_biomass = mean(mean_parrot_biomass, na.rm = TRUE),
    se_biomass   = sd(mean_parrot_biomass, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )
avg_parrot_biomass_status

avg_parrot_biomass_habitat <- avg_parrot_biomass_island %>%
  group_by(habitat) %>%
  summarise(
    mean_biomass = mean(mean_parrot_biomass, na.rm = TRUE),
    se_biomass   = sd(mean_parrot_biomass, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop")
avg_parrot_biomass_habitat


#sediment composition
write.csv(sedcomp_y, "avg_sedcomp_island.csv", row.names=F)

avg_sedcomp_status <- sedcomp_y %>%
  group_by(habitat,status) %>%
  summarise_at(vars(25:34),list(mean=mean,sd=sd))
(avg_sedcomp_status <- as.data.frame(avg_sedcomp_status))
write.csv(avg_sedcomp_status, "avg_sedcomp_status.csv", row.names=F)

avg_sedcomp_habitat <- sedcomp_y %>%
  group_by(status,habitat) %>%
  summarise_at(vars(25:34),list(mean=mean,sd=sd))
(avg_sedcomp_habitat <- as.data.frame(avg_sedcomp_habitat))
write.csv(avg_sedcomp_habitat, "avg_sedcomp_habitat.csv", row.names=F)

#for beach samples
write.csv(sedcomp_beach, "avg_sedcomp_beach_island.csv", row.names=F)

avg_sedcomp_beach_status <- sedcomp_beach %>%
  group_by(location,status) %>%
  summarise_at(vars(25:34),list(mean=mean,sd=sd))
(avg_sedcomp_beach_status <- as.data.frame(avg_sedcomp_beach_status))
write.csv(avg_sedcomp_beach_status, "avg_sedcomp_beach_status.csv", row.names=F)


#######
###
# Fig 1: plot ReefBudget and SedBudget across atolls, habitats and nutrient status 
###

#### ReefBudget
###

# includes increased coral calcification at birdy lagoon and shallow fore reef sites
net_prod_2.0.plot <- ggplot(reefbudget, aes(x = habitat, y = net_prod_2.0, fill=status))
(net_prod.plot <- net_prod_2.0.plot + geom_boxplot(show.legend=TRUE) + 
    labs(y=expression("Net carbonate budget (kg"~CaCO[3]~m^-2~yr^-1*")")) + xlab("Habitat") + theme_classic() +   scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7)) + 
    facet_wrap(~atoll)+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.text=element_text(size=12), strip.text = element_text(size = 12)))
ggsave("figures/Fig_net_prod_2.0.jpg", width = 7, height = 4)
ggsave("figures/Fig_net_prod_2.0.pdf", width = 7, height = 4)
# Salomon > GCB, PB; lagoon > shallow, deep

# wrap at site level
RB_site.plot <- ggplot(reefbudget, aes(x = habitat, y = net_prod_2.0, fill=status))
(RB_site.plot <- RB_site.plot + geom_boxplot(show.legend=TRUE) + 
    labs(y=expression("Net carbonate budget (kg"~CaCO[3]~m^-2~yr^-1*")")) + xlab("Island") + 
    ylim(-2, 16) + # adjust axis to same as RB?
    theme_classic() +   scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7)) + 
    facet_wrap(~island,nrow = 3, ncol = 2)+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
          legend.text=element_text(size=12), strip.text = element_text(size = 12)))
ggsave("figures/Fig_RB_site.jpg", width = 4, height = 7)
ggsave("figures/Fig_RB_site.pdf", width = 4, height = 7)


##### SedBudget 
###

sed_prod.plot <- ggplot(sedbudget, aes(x = habitat, y = sediment, fill=status))
(sed_prod <- sed_prod.plot + geom_boxplot(show.legend=TRUE) + 
    labs(y=expression("Sediment production (kg CaC"*O[3]~m^-2~yr^-1*")"), x="Habitat") + theme_classic() +   scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7)) + 
    facet_wrap(~atoll)+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.text=element_text(size=12), strip.text = element_text(size = 12)))
ggsave("figures/sedbudget_chagos_habitat_atoll.jpg", width = 7, height = 4)
ggsave("figures/sedbudget_chagos_habitat_atoll.pdf", width = 7, height = 4)
#birdy > ratty in lagoon, ratty > birdy on deep fore reefs

# wrap at site level
SB_site.plot <- ggplot(sedbudget, aes(x = habitat, y = sediment, fill=status))
(SB_site.plot <- SB_site.plot + geom_boxplot(show.legend=TRUE) + 
    labs(y=expression("Sediment budget (kg"~CaCO[3]~m^-2~yr^-1*")")) + xlab("Island") + 
    ylim(-2, 16) + # adjust axis to same as RB?
    theme_classic() +   scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7)) + 
    facet_wrap(~island,nrow = 3, ncol = 2)+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
          legend.text=element_text(size=12), strip.text = element_text(size = 12)))
ggsave("figures/Fig_SB_site_axis.jpg", width = 4, height = 7)
ggsave("figures/Fig_SB_site_axis.pdf", width = 4, height = 7)


#########
###
# Test impacts of status and habitat on carbonate framework and sediment production in a Bayesian framework
###

# use uninformative (default) prior because there is not much prior data on differences in geo-ecological functions aross habitats an nutrinet status (and to be conservative)
# use weakly informative half-cauchy for variance components due to small sample size 
# based on McNeish 2016 (https://www.tandfonline.com/doi/full/10.1080/10705511.2016.1186549?casa_token=8aPyiZ0aEU0AAAAA%3AAoo2jZ547UxohRAt8eKf4bOWGJ7xdK_jc3wp4t3t2OPeUtkMmOwzEUSE2PGALmxlC9jY3V-hKWj0) 

# set prior
pr_default_cauchy <- 
  set_prior("cauchy(0,2)", class = "sigma")

# net carbonate budget
#model with interaction
net_prod_2.0_mod <- brm(
  net_prod_2.0 ~ status * habitat + (1|atoll/island),
  data = reefbudget, 
  iter = 3000, warmup = 1000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.999, max_treedepth = 15), 
  sample_prior="yes",
  prior = pr_default_cauchy,
  file = "output/brms/net_prod_2.0_mod_0")
print(net_prod_2.0_mod) #0 divergent trans, 
# interaction terms are not credibly different from zero, so there’s no strong evidence 
# that the effect of status changes dramatically across habitats.

#additive model
net_prod_2.0_mod_add <- brm(
  net_prod_2.0 ~ status + habitat + (1|atoll/island),
  data = reefbudget, 
  iter = 3000, warmup = 1000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.999, max_treedepth = 15), 
  sample_prior="yes",
  prior = pr_default_cauchy,
  file = "output/brms/net_prod_2.0_mod_add")
print(net_prod_2.0_mod_add) #0 divergent trans

loo(net_prod_2.0_mod, net_prod_2.0_mod_add)
# LLOIC smaller for additive model so justified to keep simpler model, 
# but it doesn't show the difference in nutrient effect across habitats, so keep interaction

#use log link and prior of 2 (Lange & Benkwitt 2024) for statushigh:habitatlagoon and statushigh:habitatshallow
pr_status2 <- set_prior("normal(0.693, 0.2)", class = "b", coef = c("statushigh","statushigh:habitatshallow")) 
#log of 2 is 0.693, SD 0.2 means fairly confident about 2x effect

net_prod_loglink_mod <- brm(
  net_prod_2.0 ~ status * habitat + (1|atoll/island),
  data = reefbudget,
  family = gaussian(link = "log"),
  prior = c(pr_default_cauchy, pr_status2),
  sample_prior = "yes",
  file = "output/brms/net_prod_2.0_mod")
print(net_prod_loglink_mod)
#doesn't really make a difference, keep using weakly informativ priors

# SedBudget
sed_prod_mod<- brm(
  sediment~ status * habitat  + (1|atoll/island),
  data = sedbudget, 
  iter = 3000, warmup = 1000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.999, max_treedepth = 15), 
  sample_prior="no",
  file = "output/brms/sed_prod_mod_0")
print(sed_prod_mod) #1 divergent trans

##check plots----
pp_check(net_prod_2.0_mod)
plot(net_prod_2.0_mod, ask = FALSE)
plot(conditional_effects(net_prod_2.0_mod))

pp_check(sed_prod_mod)
plot(sed_prod_mod, ask = FALSE)
plot(conditional_effects(sed_prod_mod))

# coral cover
#model with interaction
coral_cover_mod <- brm(
  cover ~ status * habitat + (1|atoll/island),
  data = reefbudget, 
  iter = 3000, warmup = 1000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.999, max_treedepth = 15), 
  sample_prior="yes",
  prior = pr_default_cauchy,
  file = "output/brms/coral_cover_mod_0")
print(coral_cover_mod) #0 divergent trans, 

# parrotfish biomass
#model with interaction
parrotfish_mod <- brm(
  parrotfish_biomass ~ status * habitat + (1|atoll/island),
  data = sedbudget, 
  iter = 3000, warmup = 1000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.999, max_treedepth = 15), 
  sample_prior="yes",
  prior = pr_default_cauchy,
  file = "output/brms/parrotfish_mod_0")
print(parrotfish_mod) #0 divergent trans, 

##extract emmeans----
net_prod_2.0_mod.rg <- update(ref_grid(net_prod_2.0_mod))
net_prod_loglink_mod.rg <- update(ref_grid(net_prod_loglink_mod))
net_prod_2.0_mod_add.rg <- update(ref_grid(net_prod_2.0_mod_add))
sed_prod_mod.rg <- update(ref_grid(sed_prod_mod))
coral_cover_mod.rg <- update(ref_grid(coral_cover_mod))
parrotfish_mod.rg <- update(ref_grid(parrotfish_mod))

#contrasts by habitat - difference between ratty and birdy
contrasts_by_habitat <- emmeans(net_prod_2.0_mod.rg, ~ status | habitat, type = "response", by = "habitat") %>%
  contrast("trt.vs.ctrl")
print(contrasts_by_habitat)%>%
  as.data.frame()

#habitat = lagoon:
#  contrast   estimate  lower.HPD upper.HPD
#high - low 2.793931 -0.3347412  6.097151

#habitat = shallow:
#  contrast   estimate  lower.HPD upper.HPD
#high - low 3.259394  0.0075988  6.447697

#habitat = deep:
#  contrast   estimate  lower.HPD upper.HPD
#high - low 0.537325 -2.7716281  3.773563

#Point estimate displayed: median 
#HPD interval probability: 0.95 

# for model with log response and prior 2
#contrasts by habitat - difference between ratty and birdy
contrasts_by_habitat_pr <- emmeans(net_prod_loglink_mod.rg, ~ status | habitat, type = "response", by = "habitat") %>%
  contrast("trt.vs.ctrl")
print(contrasts_by_habitat_pr)%>%
  as.data.frame()

#habitat = lagoon:
#contrast   estimate  lower.HPD upper.HPD
#high - low 2.800388 -0.2518759  6.014386

#habitat = shallow:
#  contrast   estimate  lower.HPD upper.HPD
#high - low 3.298781  0.0863538  6.421615

#habitat = deep:
#  contrast   estimate  lower.HPD upper.HPD
#high - low 0.511081 -2.8497092  3.513395

#Point estimate displayed: median 
#HPD interval probability: 0.95 

#very similar, not necessary to include prior

summary(net_prod_2.0_mod)

# check probability that carbonate budgets are higher at sites with seabird nutrients, for each habitat
hypothesis(net_prod_2.0_mod, "statushigh>0") # for lagoon, not specified in command as this is the reference level which is included in the intercept
#Hypothesis           Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#1 (statushigh) > 0      2.8      1.62     0.28     5.32      26.49      0.96    *

hypothesis(net_prod_2.0_mod, "statushigh + statushigh:habitatshallow > 0") #for shallow
#Hypothesis                   Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#1 (statushigh+statu... > 0     3.27      1.62     0.77     5.77      40.88      0.98    *

hypothesis(net_prod_2.0_mod, "statushigh + statushigh:habitatdeep > 0") #for deep
#Hypothesis                   Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#1 (statushigh+statu... > 0     0.53      1.64    -1.99     3.06       1.88      0.65     


# for sediment budgets
contrasts_by_habitat_sed <- emmeans(sed_prod_mod.rg, ~ status | habitat, type = "response", by = "habitat") %>%
  contrast("trt.vs.ctrl")
print(contrasts_by_habitat_sed)%>%
  as.data.frame()

#habitat = lagoon:
#  contrast     estimate  lower.HPD upper.HPD
#high - low  1.0870267 -0.8621557 2.8945869

#habitat = shallow:
#  contrast     estimate  lower.HPD upper.HPD
#high - low  0.1634436 -1.6216419 2.0419123

#habitat = deep:
#  contrast     estimate  lower.HPD upper.HPD
#high - low -0.8597179 -2.7886808 0.8919853

#Point estimate displayed: median 
#HPD interval probability: 0.95 

# check probability that sediment budgets are higher at sites with seabird nutrients, for each habitat
hypothesis(sed_prod_mod, "statushigh>0") # for lagoon, not specified in command as this is the reference level which is included in the intercept
#Hypothesis           Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#1 (statushigh) > 0     1.06      0.95    -0.44     2.52       8.22      0.89    

hypothesis(sed_prod_mod, "statushigh + statushigh:habitatshallow > 0") #for shallow
#Hypothesis                   Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#1 (statushigh+statu... > 0     0.16      0.96    -1.37     1.63       1.37      0.58      

hypothesis(sed_prod_mod, "statushigh + statushigh:habitatdeep > 0") #for deep
#Hypothesis                   Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#1 (statushigh+statu... > 0    -0.87      0.94    -2.37     0.56       0.18      0.15     

#check probability of negative effect for deep?
hypothesis(sed_prod_mod, "statushigh + statushigh:habitatdeep < 0") #for deep
#Hypothesis                   Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#1 (statushigh+statu... < 0    -0.87      0.94    -2.37     0.56       5.71      0.85  

# check probability that sediment budgets are higher at sites with seabird nutrients, for each habitat
hypothesis(sed_prod_mod, "statushigh>0") # for lagoon, not specified in command as this is the reference level which is included in the intercept

hypothesis(sed_prod_mod, "statushigh + statushigh:habitatshallow > 0") #for shallow

hypothesis(sed_prod_mod, "statushigh + statushigh:habitatdeep > 0") #for deep

#check probability of negative effect for deep?
hypothesis(sed_prod_mod, "statushigh + statushigh:habitatdeep < 0") #for deep
#Hypothesis                   Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#1 (statushigh+statu... < 0    -0.87      0.94    -2.37     0.56       5.71      0.85  

# for coral cover
contrasts_by_habitat_cover <- emmeans(coral_cover_mod.rg, ~ status | habitat, type = "response", by = "habitat") %>%
  contrast("trt.vs.ctrl")
print(contrasts_by_habitat_cover)%>%
  as.data.frame()

#habitat = lagoon:
#  contrast    estimate lower.HPD upper.HPD
#high - low  0.235923 -14.20202  13.87198

#habitat = shallow:
#  contrast    estimate lower.HPD upper.HPD
#high - low  3.382418  -9.22284  18.34238

#habitat = deep:
#  contrast    estimate lower.HPD upper.HPD
#high - low -1.559568 -15.60205  12.04795

#Point estimate displayed: median 
#HPD interval probability: 0.95 

# no indication that status has an effect

# for parrotfish biomass
contrasts_by_habitat_parrotfish <- emmeans(parrotfish_mod.rg, ~ status | habitat, type = "response", by = "habitat") %>%
  contrast("trt.vs.ctrl")
print(contrasts_by_habitat_parrotfish)%>%
  as.data.frame()

#habitat = lagoon:
#  contrast   estimate lower.HPD upper.HPD
#high - low      126      -110     358.7

#habitat = shallow:
#  contrast   estimate lower.HPD upper.HPD
#high - low      107      -126     339.1

#habitat = deep:
#  contrast   estimate lower.HPD upper.HPD
#high - low     -150      -382      80.7
# no indication that status has an effect


#to report median estimates for each habitat: 
net_prod_2.0_mod.rg %>% 
  emmeans(~ status|habitat,
          type = "response")%>%
  as.data.frame()

#habitat = lagoon:
#  status   emmean  lower.HPD upper.HPD
#low    3.647557  0.0134531  7.124350
#high   6.433435  2.8038041  9.963076

#habitat = shallow:
#  status   emmean  lower.HPD upper.HPD
#low    1.648763 -2.0495829  5.073434
#high   4.928961  1.3456586  8.433377

#habitat = deep:
#  status   emmean  lower.HPD upper.HPD
#low    3.308434 -0.2083992  6.950257
#high   3.803198  0.1064975  7.180481

#Point estimate displayed: median 
#HPD interval probability: 0.95 

net_prod_2.0_mod.rg %>% 
  emmeans(~ habitat|status,
          type = "response")%>%
  as.data.frame()

#status = low:
#  habitat   emmean  lower.HPD upper.HPD
#lagoon  3.647557  0.0134531  7.124350
#shallow 1.648763 -2.0495829  5.073434
#deep    3.308434 -0.2083992  6.950257

#status = high:
#  habitat   emmean  lower.HPD upper.HPD
#lagoon  6.433435  2.8038041  9.963076
#shallow 4.928961  1.3456586  8.433377
#deep    3.803198  0.1064975  7.180481

#Point estimate displayed: median 
#HPD interval probability: 0.95 

net_prod_2.0_mod.rg %>% 
  emmeans(~ habitat,
          type = "response")%>%
  as.data.frame()

#NOTE: Results may be misleading due to involvement in interactions
#habitat   emmean lower.HPD upper.HPD
#lagoon  5.047762 1.7526138  8.162161
#shallow 3.300074 0.1004947  6.509413
#deep    3.541486 0.3345422  6.637947

#Results are averaged over the levels of: status 
#Point estimate displayed: median 
#HPD interval probability: 0.95 

#no interaction but from additive model: similar
net_prod_2.0_mod_add.rg %>% 
  emmeans(~ habitat, #or habitat|status
          type = "response")%>%
  as.data.frame()
#habitat   emmean  lower.HPD upper.HPD
#lagoon  5.011333  1.5856506  8.199772
#shallow 3.254941 -0.1918998  6.568882
#deep    3.514319  0.1076691  6.723680

#Results are averaged over the levels of: status 
#Point estimate displayed: median 
#HPD interval probability: 0.95 

#habitat = lagoon:
#  status   emmean  lower.HPD upper.HPD
#low    3.909907  0.1683333  7.398539
#high   6.091607  2.5953872  9.788971

#habitat = shallow:
#  status   emmean  lower.HPD upper.HPD
#low    2.167408 -1.6292797  5.717922
#high   4.341530  0.7109990  7.945181

#habitat = deep:
#  status   emmean  lower.HPD upper.HPD
#low    2.406293 -1.3602595  5.869785
#high   4.595243  0.7397408  7.901726

#Point estimate displayed: median 
#HPD interval probability: 0.95 


### SedBudget
sed_prod_mod.rg %>% 
  emmeans(~ status|habitat,
          type = "response")%>%
  as.data.frame()

#habitat = lagoon:
#  status   emmean  lower.HPD upper.HPD
#low    0.748255 -0.9878181  2.433784
#high   1.839034  0.1287108  3.409955

#habitat = shallow:
#  status   emmean  lower.HPD upper.HPD
#low    1.250958 -0.4592375  2.861964
#high   1.416586 -0.2083153  3.143094

#habitat = deep:
#  status   emmean  lower.HPD upper.HPD
#low    3.381284  1.7201347  5.034320
#high   2.530148  0.8186400  4.170244

#Point estimate displayed: median 
#HPD interval probability: 0.95 

sed_prod_mod.rg %>% 
  emmeans(~ habitat|status,
          type = "response")%>%
  as.data.frame()

#status = low:
#  habitat   emmean  lower.HPD upper.HPD
#lagoon  0.748255 -0.9878181  2.433784
#shallow 1.250958 -0.4592375  2.861964
#deep    3.381284  1.7201347  5.034320

#status = high:
#  habitat   emmean  lower.HPD upper.HPD
#lagoon  1.839034  0.1287108  3.409955
#shallow 1.416586 -0.2083153  3.143094
#deep    2.530148  0.8186400  4.170244

#Point estimate displayed: median 
#HPD interval probability: 0.95 

sed_prod_mod.rg %>% 
  emmeans(~ habitat,
          type = "response")%>%
  as.data.frame()

#NOTE: Results may be misleading due to involvement in interactions
#habitat   emmean  lower.HPD upper.HPD
#lagoon  1.293417 -0.0675712  2.733495
#shallow 1.342902 -0.1597490  2.642314
#deep    2.952932  1.6021870  4.358363

#Results are averaged over the levels of: status 
#Point estimate displayed: median 
#HPD interval probability: 0.95

coral_cover_mod.rg %>% 
  emmeans(~ habitat|status,
          type = "response")%>%
  as.data.frame()
#status = low:
#  habitat   emmean lower.HPD upper.HPD
#lagoon  26.71378 10.878603  42.99414
#shallow 22.27875  7.379494  39.73824
#deep    38.51279 22.321084  54.27961

#status = high:
#  habitat   emmean lower.HPD upper.HPD
#lagoon  26.95144 11.517160  42.88132
#shallow 25.67027 10.507202  41.78639
#deep    36.96300 21.474897  53.05828

coral_cover_mod.rg %>% 
  emmeans(~ habitat,
          type = "response")%>%
  as.data.frame()
#habitat   emmean lower.HPD upper.HPD
#lagoon  26.80612  12.89452  41.25320
#shallow 23.97795  10.11086  38.76336
#deep    37.72622  23.40936  51.75288

#Results are averaged over the levels of: status 
#Point estimate displayed: median 
#HPD interval probability: 0.95 

parrotfish_mod.rg %>% 
  emmeans(~ habitat|status,
          type = "response")%>%
  as.data.frame()
#status = low:
#  habitat   emmean  lower.HPD upper.HPD
#lagoon   95.4360 -105.12350  284.8079
#shallow 111.6318  -98.44138  300.6364
#deep    365.1623  162.04168  550.7985

#status = high:
#  habitat   emmean  lower.HPD upper.HPD
#lagoon  219.5778   37.03976  429.8699
#shallow 218.1139   31.58310  414.9339
#deep    215.6312   13.88266  398.1442

parrotfish_mod.rg %>% 
  emmeans(~ habitat,
          type = "response")%>%
  as.data.frame()
#habitat   emmean lower.HPD upper.HPD
#lagoon  157.6663  -6.22605  309.2870
#shallow 164.2706   6.73267  321.5505
#deep    291.3226 138.61025  447.2099

#Results are averaged over the levels of: status 
#Point estimate displayed: median 
#HPD interval probability: 0.95 


################
###########
# Fig 2: plot posterior distributions by habitat and nutrient status
####

# Define colours for plotting
high_colors <- c(
  "lagoon"  = adjustcolor("#0168A5", alpha.f = 0.4),
  "shallow" = adjustcolor("#0168A5", alpha.f = 0.8),
  "deep"    = adjustcolor("#004D7E", alpha.f = 0.7)
)

low_colors <- c(
  "lagoon"  = adjustcolor("#BF1F36", alpha.f = 0.4),
  "shallow" = adjustcolor("#BF1F36", alpha.f = 0.8),
  "deep"    = adjustcolor("#BF1F36", alpha.f = 1)
)


#### for framework production 

# Step 1: Summarise per draw, status, and habitat
pred_draws_framework <- net_prod_2.0_mod %>%
  add_epred_draws(
    newdata = reefbudget,
    re_formula = ~(1 | atoll)
  )

pred_habitat_framework <- pred_draws_framework %>%
  group_by(.draw, status, habitat) %>%
  summarise(mean_epred = mean(.epred), .groups = "drop") %>%
  mutate(
    habitat = factor(habitat, levels = c("deep", "shallow", "lagoon")),
    status  = factor(status, levels = c("low", "high"))
  )

# Step 2: Add a new variable for fill colours depending on status + habitat
pred_habitat_framework <- pred_habitat_framework %>%
  mutate(
    fill_col = ifelse(
      status == "high",
      high_colors[as.character(habitat)],
      low_colors[as.character(habitat)]
    ))

# Step 3: Plot with facets for each status
(framework_habitat <- ggplot(pred_habitat_framework, aes(y = habitat, x = mean_epred, fill = fill_col)) +
    stat_halfeye(
      .width = c(0.5, 0.7, 0.9),
      #alpha  = 1,
      colour = "black"
    ) +
    facet_wrap(~ status, ncol = 1) +
    scale_fill_identity() +
    xlim(-3,12) +
    labs(
      x = expression("Predicted framework production (kg m"^2~yr^-1*")"),
      y = ""
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      #legend.position = "none",
      rect = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      text = element_text(size = 12, family = "sans"),
      legend.text = element_text(size = 12, family = "sans"),
      axis.title.y = element_text(size = 12, family = "sans"),
      axis.text = element_text(size = 12, family = "sans")
    ) )

#### for sediment production 

# Step 1: Summarise per draw, status, and habitat
pred_draws_sediment <- sed_prod_mod %>%
  add_epred_draws(
    newdata = sedbudget,
    re_formula = ~(1 | atoll)
  )

pred_habitat_sediment <- pred_draws_sediment %>%
  group_by(.draw, status, habitat) %>%
  summarise(mean_epred = mean(.epred), .groups = "drop") %>%
  mutate(
    habitat = factor(habitat, levels = c("deep", "shallow", "lagoon")),
    status  = factor(status, levels = c("low", "high"))
  )

# Step 2: Add a new variable for fill colours depending on status + habitat
pred_habitat_sediment <- pred_habitat_sediment %>%
  mutate(
    fill_col = ifelse(
      status == "high",
      high_colors[as.character(habitat)],
      low_colors[as.character(habitat)]
    ))

# Step 3: Plot with facets for each status
(sediment_habitat <- ggplot(pred_habitat_sediment, aes(y = habitat, x = mean_epred, fill = fill_col)) +
    stat_halfeye(
      .width = c(0.5, 0.7, 0.9),
      #alpha  = 1,
      colour = "black"
    ) +
    facet_wrap(~ status, ncol = 1) +
    scale_fill_identity() +
    xlim(-2,6) +
    labs(
      x = expression("Predicted sediment production (kg m"^2~yr^-1*")"),
      y = ""
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      #legend.position = "none",
      rect = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      text = element_text(size = 12, family = "sans"),
      legend.text = element_text(size = 12, family = "sans"),
      axis.title.y = element_text(size = 12, family = "sans"),
      axis.text = element_text(size = 12, family = "sans")
    ) )


#############
### plot effect sizes for seabird nutrients

# for framework budgets
#get estimates and HPDs:
net_prod_2.0_pred_7<-net_prod_2.0_mod.rg%>%
  emmeans(~ status|habitat,  
          type="response")%>%
  contrast("trt.vs.ctrl", level = .7)%>%
  as.data.frame()

net_prod_2.0_pred_9<-net_prod_2.0_mod.rg%>%
  emmeans(~ status|habitat, 
          type="response")%>%
  contrast("trt.vs.ctrl", level = .9)%>%
  as.data.frame()

net_prod_2.0_mod.rg %>% 
  emmeans(~ status|habitat,
          type = "response")%>%
  contrast("trt.vs.ctrl", level = .9)%>%
  as.data.frame()

#combine estimates:
net_prod_est<-data.frame(effect = c( 'lagoon', 'shallow', 'deep'),
                         estimate = c(net_prod_2.0_pred_7$estimate,net_prod_2.0_pred_7$estimate), #effect size is the same regardless of interval, so can use either for these
                         low = c (net_prod_2.0_pred$lower.HPD, net_prod_2.0_pred_7$lower.HPD),
                         up = c( net_prod_2.0_pred_9$upper.HPD, net_prod_2.0_pred_7$upper.HPD),
                         level = c( .9, .9, .9, .7, .7, .7))
net_prod_est

net_prod_effect_plot<-
  net_prod_est %>%
  mutate(effect = fct_relevel(effect, c("deep", "shallow", "lagoon")))%>%
  ggplot(aes(x = estimate, y = effect, color = as.factor(level)))+
  geom_vline(xintercept=0, lty=2, alpha = .8, color = "darkgrey")+
  geom_errorbar(aes(xmin = low, xmax = up), width= 0, linewidth = c(1, 1, 1, 2, 2, 2), alpha = .7)+
  geom_point(size =5, alpha = .9)+
  scale_color_manual(values = c("#0067A5", "#A1CAF1"))+ 
  scale_fill_manual(values = c("#0067A5", "#A1CAF1"))+
  #scale_y_discrete(labels = c("lagoon", "shallow", "deep"))+
  labs(x=expression("Nutrient effect (kg"~m^-2~yr^-1*")"), y="")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(), 
        legend.position = "none",
        rect = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA),
        text=element_text(size=12,  family="sans"),
        legend.text=element_text(size=12, family = "sans"),
        axis.title.y = element_text(size=12, family = "sans"),
        axis.text = element_text(size = 12, family = "sans")) 
net_prod_effect_plot

## for sediment budgets
#get estimates and HPDs:
sed_prod_pred_7<-sed_prod_mod.rg%>%
  emmeans(~ status|habitat,  
          type="response")%>%
  contrast("trt.vs.ctrl", level = .7)%>%
  as.data.frame()

sed_prod_pred_9<-sed_prod_mod.rg%>%
  emmeans(~ status|habitat, 
          type="response")%>%
  contrast("trt.vs.ctrl", level = .9)%>%
  as.data.frame()

sed_prod_mod.rg %>% 
  emmeans(~ status|habitat,
          type = "response")%>%
  contrast("trt.vs.ctrl", level = .9)%>%
  as.data.frame()

#combine estimates:
sed_prod_est<-data.frame(effect = c( 'lagoon', 'shallow', 'deep'),
                         estimate = c(sed_prod_pred_7$estimate,sed_prod_pred_7$estimate), #effect size is the same regardless of interval, so can use either for these
                         low = c (sed_prod_pred_9$lower.HPD, sed_prod_pred_7$lower.HPD),
                         up = c( sed_prod_pred_9$upper.HPD, sed_prod_pred_7$upper.HPD),
                         level = c( .9, .9, .9, .7, .7, .7))
sed_prod_est

#horizontal
sed_prod_effect_plot<-
  sed_prod_est %>%
  mutate(effect = fct_relevel(effect, c("deep", "shallow", "lagoon")))%>%
  ggplot(aes(x = estimate, y = effect, color = as.factor(level)))+
  geom_vline(xintercept=0, lty=2, alpha = .8, color = "darkgrey")+
  geom_errorbar(aes(xmin = low, xmax = up), width= 0, linewidth = c(1, 1, 1, 2, 2, 2), alpha = .7)+
  geom_point(size =5, alpha = .9)+
  scale_color_manual(values = c("#0067A5", "#A1CAF1"))+ 
  scale_fill_manual(values = c("#0067A5", "#A1CAF1"))+
  #scale_y_discrete(labels = c("lagoon", "shallow", "deep"))+
  labs(x=expression("Nutrient effect (kg"~m^-2~yr^-1*")"), y="")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(), 
        legend.position = "none",
        rect = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA),
        text=element_text(size=12,  family="sans"),
        legend.text=element_text(size=12, family = "sans"),
        axis.title.y = element_text(size=12, family = "sans"),
        axis.text = element_text(size = 12, family = "sans")) 
sed_prod_effect_plot

#### combine plots to Fig 2
Fig2 <- plot_grid(
  framework_habitat, sediment_habitat, # Row 1
  net_prod_effect_plot, sed_prod_effect_plot,   # Row 2
  ncol = 2,
  labels = c("A","B","C","D"),
  align = "v",       # align vertically
  axis = "lr",       # align left and right axes
  rel_heights = c(2, 1) # first row double height
)
Fig2

ggsave("figures/Fig2_carbonate_production.pdf", width = 10, height = 8)
ggsave("figures/Fig2_carbonate_production.jpg", width = 10, height = 8)



###########
########
## Fig 3: Plot sediment composition and grain size distributions
#####

# Produced on reef
#transform contributions to long format
contr_long <- avg_sedprod_status %>%
  select(habitat, status, ends_with("_mean"), -sediment_mean, -urchins_mean, -parrotfish_mean, -sponge_mean) %>% #as urchin, parrotfish and sponge sediment production have been transformed to coral and cca constituents
  pivot_longer(
    cols = ends_with("_mean"),
    names_to = "group",
    values_to = "sed"
  ) %>%
  mutate(group = str_remove(group, "_mean"))  # removes "_mean" or "_mean_mean"

contr_order <- c("coral", "cca", "molluscs", "forams", "algae")
contr_long$group <- factor(contr_long$group, levels = rev(contr_order))

#color palette
prod_palette <- c(
  #"parrotfish" = "#CC667798",
  "coral" = "#CC6677",
  "cca"            = "#FB9A99",   
  #"urchins"        = "#CC667760",   
  "molluscs"        = "#A6CEE3",   
  #"Crustacea"      = "#1F78B4",   
  #"sponge"       = "#CC667720",   
  "forams"   = "#F0E442",   
  "algae"       = "#117733"   
  #"Other"          = "darkgrey"    
)

# plot in kg
(p9 <- ggplot(contr_long, aes(x = status, y = sed, fill = group)) +
    geom_bar(stat = "identity") + 
    labs(x = "habitat", y = expression("Sediment production (kg "~m^-2~yr^-1*")")) +
    facet_wrap(~habitat)+
    theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
                            legend.title = element_text(size=12),legend.text = element_text(size=12))+
    scale_fill_manual(values = prod_palette))

#plot as % of total
(p9b <- ggplot(contr_long, aes(x = status, y = sed, fill = group)) +
    geom_bar(stat = "identity", position="fill") + #position="fill" to display 100% for each
    labs(x = "habitat", y = expression("Sediment production (%)")) +
    facet_wrap(~habitat)+
    theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
                            legend.title = element_text(size=12),legend.text = element_text(size=12),
                            strip.text = element_text(size = 12),
                            legend.position = "none", #legend.justification = "center",legend.box = "horizontal"
    )+
    scale_fill_manual(values = prod_palette)+
    scale_y_continuous(labels = function(x) x * 100)) # changes y-axis from 0-1 to 0-100

#####
# Found in reef sediment
#transform contributions to long format
comp_y_long <- sedcomp_y[,c(2:8,27:35)] %>%
  pivot_longer(cols=c("Coral", "Halimeda", "Mollusc", "Foraminifera","Crustacea","Echinoderms", "Soft coral", "CCA", "Other"), 
               names_to="group", values_to="sed")

#calculate averages for each status
comp_y_avg <- comp_y_long %>%
  group_by(status, habitat, group) %>%
  summarise_at(
    vars(6),
    list(contr = ~mean(.x, na.rm = TRUE),
         sd   = ~sd(.x, na.rm = TRUE))) %>%
  as.data.frame()

comp_y_order <- c("Coral", "CCA", "Mollusc", "Crustacea","Echinoderms", "Soft coral", "Foraminifera","Halimeda", "Other")
comp_y_avg$group <- factor(comp_y_avg$group, levels = rev(comp_y_order))

#Define coral palette
comp_palette <- c(
  "Coral" = "#CC6677",   # 
  "CCA"            = "#FB9A99",   #  
  "Mollusc"        = "#A6CEE3",   # sky blue
  "Crustacea"      = "#1F78B4",   # blue
  "Echinoderms"        = "#13486A",   # vermillion
  "Soft coral"       ="#FDBF6F",   # bluish green
  "Foraminifera"   = "#F0E442",   # yellow
  "Halimeda"       = "#117733",   # dark green
  "Other"          = "darkgrey"    # grey
)

# plot 
(p10 <- ggplot(comp_y_avg, aes(x = status, y = contr, fill = group)) +
    geom_bar(stat = "identity") + #position="fill" to display 100% for each
    labs(x = "habitat", y = expression("Sediment composition (%)")) +
    facet_wrap(~habitat)+
    theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
                            legend.title = element_text(size=12),legend.text = element_text(size=12),
                            strip.text = element_text(size = 12),
                            legend.position = "right",legend.justification = "center",legend.box = "horizontal"
    )+ 
    scale_fill_manual(values = comp_palette))

# Combine plots side by side
# Extract the legend from one plot
legend <- get_legend(p10 + theme(legend.position = "right", return_all=TRUE))  # or "top"/"bottom" as needed
# Remove legends from the main plots
p9b_clean <- p9b + theme(legend.position = "none")
p10_clean <- p10 + theme(legend.position = "none")
# Combine plots side by side
plots_combined_side <- plot_grid(p9b_clean, p10_clean, ncol = 2, labels = c("A", "B"))
#plots_combined <- plot_grid(p9b_clean, p10_clean, ncol = 1, labels = c("A", "C"))
# Add the legend on the right
Fig_sed_comp_reef <- plot_grid(plots_combined_side,legend, ncol = 2, rel_widths = c(1, 0.2))  
Fig_sed_comp_reef


#####
# Found in beach sediment
#transform contributions to long format
comp_beach_long <- sedcomp_beach[,c(2:8,27:35)] %>%
  pivot_longer(cols=c("Coral", "Halimeda", "Mollusc", "Foraminifera","Crustacea","Soft coral", "Echinoderms", "CCA", "Other"), 
               names_to="group", values_to="sed")

#calculate averages for each status
comp_beach_avg <- comp_beach_long %>%
  group_by(status, location, habitat, group) %>%
  summarise_at(
    vars(5),
    list(contr = ~mean(.x, na.rm = TRUE),
         sd   = ~sd(.x, na.rm = TRUE))) %>%
  as.data.frame()

comp_beach_order <- c("Coral", "CCA", "Mollusc", "Crustacea","Echinoderms", "Soft coral", "Foraminifera","Halimeda", "Other")
comp_beach_avg$group <- factor(comp_beach_avg$group, levels = rev(comp_beach_order))
comp_beach_avg$location <- recode(comp_beach_avg$location,
                                  "Lagoon" = "lagoon beach",
                                  "Outer reef" = "fore reef beach")
comp_beach_avg$location <- factor(comp_beach_avg$location,
                                  levels = c("lagoon beach", "fore reef beach", "empty")) # add dummy variable "empty" to allow for aligned plots in Fig 3

# plot 
(p11 <- ggplot(comp_beach_avg, aes(x = status, y = contr, fill = group)) +
    geom_bar(stat = "identity") + # position="fill" to display 100% for each
    labs(x = "habitat", y = expression("Sediment composition (%)")) +
    facet_wrap(~location, nrow=1, drop=FALSE) +
    theme_classic() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12), 
          legend.title = element_text(size=12),
          legend.text = element_text(size=12),
          strip.text = element_text(size = 12)) +
    scale_fill_manual(values = comp_palette))

#####
## Sediment grain size distribution

#####
# for sediment production
#transform grain sizes to long format
grain_long <- sedbudget[, c(2:6, 14, 19:30)] %>%
  rowwise() %>%
  mutate(total_sediment = sum(c_across(c(`>2000`, `1000-2000`, `500-1000`, `250-500`, #sum size fractions to calclate proportions
                                         `125-250`, `63-125`, `31-63`, `16-31`, 
                                         `8-16`, `<8`)), na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_longer(
    cols = c(`>2000`, `1000-2000`, `500-1000`, `250-500`, `125-250`, `63-125`, 
             `31-63`, `16-31`, `8-16`, `<8`),
    names_to = "grainsize",
    values_to = "grain"
  ) %>%
  mutate(proportion = grain / total_sediment * 100)

#calculate averages for each habitat and status
grain_avg <- grain_long %>%
  group_by(habitat, status, grainsize) %>%
  summarise(
    grain_prod = mean(grain, na.rm = TRUE),
    grain_sd = sd(grain, na.rm = TRUE),
    proportion_prod = sum(grain, na.rm = TRUE) / sum(total_sediment, na.rm = TRUE) * 100,
    .groups = "drop"
  )

grainsize_order <- c(">2000", "1000-2000", "500-1000", "250-500", "125-250", "63-125", "31-63", "16-31", "8-16","<8")    
#"pebble", "gravel", "very_coarse_sand", "coarse_sand", "medium_sand", "fine_sand", "very_fine_sand", "coarse_silt", "medium_silt", "fine_silt", "very_fine_silt", "clay"
grain_avg$grainsize <- factor(grain_avg$grainsize, levels = (grainsize_order))

# plot
(p9_grain <- ggplot(grain_avg, aes(x = grainsize, y = grain_prod, fill = status)) +
    geom_bar(stat = "identity", position="dodge") + #position="fill" to display 100% for each
    labs(x = "grain size", y = expression("Sediment production (kg "~m^-2~yr^-1*")")) +
    scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7)) +
    facet_wrap(~habitat)+
    theme_classic() + theme(axis.text=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1, size=7), axis.title=element_text(size=12), 
                            strip.text = element_text(size = 12),
                            legend.title = element_text(size=12),legend.text = element_text(size=12),       
                            legend.position = c(0.99, 0.99), legend.justification = c("right", "top")))

(p9a_proportion <- ggplot(grain_avg, aes(x = grainsize, y = proportion_prod, fill = status)) +
    geom_bar(stat = "identity", position="dodge") + #position="fill" to display 100% for each
    labs(x = "grain size", y = "Grain size contribution (%)") +
    scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7)) +
    facet_wrap(~habitat)+
    ylim(0,45)+
    theme_classic() + theme(axis.text=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1, size=8), axis.title=element_text(size=12), 
                            strip.text = element_text(size = 12),
                            legend.title = element_text(size=12),legend.text = element_text(size=12),       
                            legend.position = "none"))        #c(0.99, 0.99), legend.justification = c("right", "top"))) # stick to corner

####
### for sediment composition reef

#transform grain sizes to long format
graincomp_long <- sedcomp_y[,c(2:9,11:20)] %>%
  pivot_longer(cols=c(">2000", "1000-2000", "500-1000", "250-500",
                      "125-250", "63-125", "31-63", "16-31", "8-16", "<8"), 
               names_to="grainsize", values_to="grain")

graincomp_long$grainsize <- factor(graincomp_long$grainsize, levels = (grainsize_order))

#calculate averages for each status
graincomp_avg <- graincomp_long %>%
  group_by(habitat,status,grainsize) %>%
  summarise_at(vars(7),list(prod=~mean(.x, na.rm = TRUE),sd=~sd(.x, na.rm = TRUE)))
graincomp_avg <- as.data.frame(graincomp_avg)

graincomp_avg$grainsize <- factor(graincomp_avg$grainsize, levels = (grainsize_order))

# plot
(p9b_proportion <- ggplot(graincomp_avg, aes(x = grainsize, y = prod, fill = status)) +
    geom_bar(stat = "identity", position="dodge") + #position="fill" to display 100% for each
    labs(x = "grain size", y = "Grain size contribution (%)") +
    scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7)) +
    facet_wrap(~habitat)+
    ylim(0,45)+
    theme_classic() + theme(axis.text=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1, size=8), axis.title=element_text(size=12), 
                            strip.text = element_text(size = 12),
                            legend.title = element_text(size=12),legend.text = element_text(size=12),       
                            legend.position = c(0.99, 0.99), legend.justification = c("right", "top")
    ))

#in detail
(p9b_site <- ggplot(graincomp_long, aes(x = grainsize, y = grain, fill = status)) +
    geom_bar(stat = "identity", position="dodge") + #position="fill" to display 100% for each
    labs(x = "grain size", y = "Grain size contribution (%)") +
    scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7)) +
    facet_wrap(~atoll + habitat)+
    theme_classic() + theme(axis.text=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1, size=7), axis.title=element_text(size=12), 
                            strip.text = element_text(size = 12),
                            legend.title = element_text(size=12),legend.text = element_text(size=12),       
                            legend.position = c(0.99, 0.99), legend.justification = c("right", "top")))

Fig_grainsize_reef <- plot_grid(p9a_proportion, p9b_proportion, ncol = 2, labels = c("C", "D"))

###
### for sediment composition beach 

#transform grain sizes to long format
graincomp.beach_long <- sedcomp_beach[,c(2:9,11:20)] %>% # also shallow water data if using sedcomp_a
  pivot_longer(cols=c(">2000", "1000-2000", "500-1000", "250-500",
                      "125-250", "63-125", "31-63", "16-31", "8-16", "<8"), 
               names_to="grainsize", values_to="grain")

graincomp.beach_long$location <- recode(graincomp.beach_long$location,
                                        "Lagoon" = "lagoon beach",
                                        "Outer reef" = "fore reef beach")
graincomp.beach_long$location <- factor(graincomp.beach_long$location,
                                        levels = c("lagoon beach", "fore reef beach", "empty")) # add dummy varibale "empty" to allow alignement of plots in Fig 3

#calculate averages for each status
graincomp.beach_avg <- graincomp.beach_long %>%
  group_by(location,status, grainsize) %>%
  summarise_at(vars(7),list(prod=~mean(.x, na.rm = TRUE),sd=~sd(.x, na.rm = TRUE)))
graincomp.beach_avg <- as.data.frame(graincomp.beach_avg)

graincomp.beach_avg$grainsize <- factor(graincomp.beach_avg$grainsize, levels = (grainsize_order))

# plot
(p9c_proportion <- ggplot(graincomp.beach_avg, aes(x = grainsize, y = prod, fill = status)) +
    geom_bar(stat = "identity", position="dodge") + #position="fill" to display 100% for each
    labs(x = "grain size", y = "Grain size contribution (%)") +
    scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7)) +
    facet_wrap(~location, nrow = 1, drop = FALSE)+ #+ habitat if plotting shallow water and beach data
    theme_classic() + theme(axis.text=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1, size=8), axis.title=element_text(size=12), 
                            strip.text = element_text(size = 12),
                            legend.title = element_text(size=12),legend.text = element_text(size=12),       
                            legend.position = c(0.99, 0.99), legend.justification = c("right", "top")))

#in detail
(p9c_site <- ggplot(graincomp.beach_long, aes(x = grainsize, y = grain, fill = status)) +
    geom_bar(stat = "identity", position="dodge") + #position="fill" to display 100% for each
    labs(x = "grain size", y = "Grain size contribution (%)") +
    scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7)) +
    facet_wrap(~atoll + location + habitat)+
    theme_classic() + theme(axis.text=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1, size=7), axis.title=element_text(size=12), 
                            strip.text = element_text(size = 12),
                            legend.title = element_text(size=12),legend.text = element_text(size=12),       
                            legend.position = c(0.99, 0.99), legend.justification = c("right", "top")))

####################
# Combine Figures

# composition
p11_clean <- p11 + theme(legend.position = c(0.99, 0.99), legend.justification = c("right", "top"))
plots_combined <- plot_grid(p9b_clean, p10_clean, p11_clean, ncol = 1, labels = c("A", "C", "E"))

# grain sizes
p9a_proportion_clean <- p9a_proportion + theme(legend.position = "none")
p9b_proportion_clean <- p9b_proportion + theme(legend.position = "none")
plots_combined2 <- plot_grid(p9a_proportion_clean, p9b_proportion_clean, p9c_proportion, ncol = 1, labels = c("B", "D", "F"))

#all in one plot
Fig3_allsand <- plot_grid(plots_combined, plots_combined2, ncol=2, rel_widths = c(1, 1))
Fig3_allsand
ggsave("figures/Fig3_sediment_composition_beach.jpg", width = 10, height = 12)
ggsave("figures/Fig3_sediment_composition_beach.pdf", width = 10, height = 12) # final touches in Adobe Illustrator


################
########
## Fig S4: calculate and plot sediment sorting 
###

# Define size classes and midpoints (in µm)
grainsize_classes <- data.frame(
  grainsize = c(">2000", "1000-2000", "500-1000", "250-500", "125-250", "63-125"
                #, "31-63", "16-31", "8-16", "<8"   ### only use size classes >63 for better comparison 
  ),
  midpoint_um = c(3000, 1500, 750, 375, 187.5, 94 
                  #, 47, 23.5, 12, 4
  ))

# Convert to phi scale (φ = -log2(d/1000))
grainsize_classes <- grainsize_classes %>%
  mutate(phi_mid = -log2(midpoint_um / 1000))

grainsize_classes

# Join phi midpoints to  data
grain_long_phi <- grain_long %>%
  left_join(grainsize_classes, by = "grainsize")%>%
  mutate(prop_frac = proportion / 100)  # as proportion is in %

graincomp_long_phi <- graincomp_long %>%
  left_join(grainsize_classes, by = "grainsize") %>%
  mutate(prop_frac = grain / 100)  # as grain is in %

graincomp.beach_long_phi <- graincomp.beach_long %>%
  left_join(grainsize_classes, by = "grainsize") %>%
  mutate(prop_frac = grain / 100)  # as grain is in %

# Calculate sorting (and mean phi) for sediment production
sedprod_sorting_stats <- grain_long_phi %>%
  group_by(island, habitat, status) %>%
  summarise(
    mean_phi = sum(phi_mid * prop_frac, na.rm = TRUE) / sum(prop_frac, na.rm = TRUE),
    sorting_phi = sqrt(sum(((phi_mid - mean_phi)^2) * prop_frac, na.rm = TRUE) / sum(prop_frac, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(mean_um = 1000 * 2^(-mean_phi))

(sedprod_sorting_stats <- sedprod_sorting_stats %>%
    mutate(
      sorting_class = case_when(
        sorting_phi < 0.35 ~ "very well sorted",
        sorting_phi < 0.5  ~ "well sorted",
        sorting_phi < 0.71 ~ "moderately well sorted",
        sorting_phi < 1.00 ~ "moderately sorted",
        sorting_phi < 2.00 ~ "poorly sorted",
        sorting_phi < 4.00 ~ "very poorly sorted",
        TRUE ~ "extremely poorly sorted"
      )))
# all poorly to very poorly sorted

# summary table
sedprod_summary <- sedprod_sorting_stats %>%
  #group_by(habitat, status) %>%
  summarise(
    mean_grain_um  = mean(mean_um, na.rm = TRUE),
    sd_grain_um    = sd(mean_um, na.rm = TRUE),
    mean_sorting   = mean(sorting_phi, na.rm = TRUE),
    sd_sorting     = sd(sorting_phi, na.rm = TRUE),
    n_samples      = n(),
    .groups = "drop"
  )

(sedprod_sorting <- ggplot(sedprod_sorting_stats, aes(x = habitat, y = sorting_phi, fill = status)) +
    geom_boxplot() +
    labs(y = "Sorting (phi SD)", x = "") +
    theme_classic()
)

# Calculate sorting (and mean phi) for sediment composition reef
sedcomp_sorting_stats <- graincomp_long_phi %>%
  group_by(island, habitat, status) %>%
  summarise(
    mean_phi = sum(phi_mid * prop_frac, na.rm = TRUE) / sum(prop_frac, na.rm = TRUE),
    sorting_phi = sqrt(sum(((phi_mid - mean_phi)^2) * prop_frac, na.rm = TRUE) / sum(prop_frac, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(mean_um = 1000 * 2^(-mean_phi))


(sedcomp_sorting_stats <- sedcomp_sorting_stats %>%
    mutate(
      sorting_class = case_when(
        sorting_phi < 0.35 ~ "very well sorted",
        sorting_phi < 0.5  ~ "well sorted",
        sorting_phi < 0.71 ~ "moderately well sorted",
        sorting_phi < 1.00 ~ "moderately sorted",
        sorting_phi < 2.00 ~ "poorly sorted",
        sorting_phi < 4.00 ~ "very poorly sorted",
        TRUE ~ "extremely poorly sorted"
      )))
# moderately well to poorly sorted

# summary table
sedcomp_summary <- sedcomp_sorting_stats %>%
  #group_by(habitat,status) %>%
  summarise(
    mean_grain_um  = mean(mean_um, na.rm = TRUE),
    sd_grain_um    = sd(mean_um, na.rm = TRUE),
    mean_sorting   = mean(sorting_phi, na.rm = TRUE),
    sd_sorting     = sd(sorting_phi, na.rm = TRUE),
    n_samples      = n(),
    .groups = "drop"
  )

(sedcomp_sorting <- ggplot(sedcomp_sorting_stats, aes(x = habitat, y = sorting_phi, fill = status)) +
    geom_boxplot() +
    labs(y = "Sorting (phi SD)", x = "Habitat") +
    theme_classic()
)

# Calculate sorting (and mean phi) for sediment composition beach
sedcomp.beach_sorting_stats <- graincomp.beach_long_phi %>%
  group_by(island, location, status) %>%
  summarise(
    mean_phi = sum(phi_mid * prop_frac, na.rm = TRUE) / sum(prop_frac, na.rm = TRUE),
    sorting_phi = sqrt(sum(((phi_mid - mean_phi)^2) * prop_frac, na.rm = TRUE) / sum(prop_frac, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(mean_um = 1000 * 2^(-mean_phi))

(sedcomp.beach_sorting_stats <- sedcomp.beach_sorting_stats %>%
    mutate(
      sorting_class = case_when(
        sorting_phi < 0.35 ~ "very well sorted",
        sorting_phi < 0.5  ~ "well sorted",
        sorting_phi < 0.71 ~ "moderately well sorted",
        sorting_phi < 1.00 ~ "moderately sorted",
        sorting_phi < 2.00 ~ "poorly sorted",
        sorting_phi < 4.00 ~ "very poorly sorted",
        TRUE ~ "extremely poorly sorted"
      )))
# mostly moderately well sorted

# summary table
sedcomp.beach_summary <- sedcomp.beach_sorting_stats %>%
  group_by(location, status) %>%
  summarise(
    mean_grain_um  = mean(mean_um, na.rm = TRUE),
    sd_grain_um    = sd(mean_um, na.rm = TRUE),
    mean_sorting   = mean(sorting_phi, na.rm = TRUE),
    sd_sorting     = sd(sorting_phi, na.rm = TRUE),
    n_samples      = n(),
    .groups = "drop"
  )

(sedcomp.beach_sorting <- ggplot(sedcomp.beach_sorting_stats, aes(x = location, y = sorting_phi, fill = status)) +
    geom_boxplot() +
    labs(y = "Sorting (phi SD)", x = "") +
    theme_classic()
)

# adjust colors and y-axes of plots
sedprod_sorting <- sedprod_sorting + theme(axis.text = element_text(size=12),
                                                        axis.title = element_text(size=12), 
                                                        legend.title = element_text(size=12),
                                                        legend.text = element_text(size=12),
                                                        legend.position = "none") + 
  ylim(0.5,1.4) + 
  scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7))

sedcomp_sorting <- sedcomp_sorting + theme(axis.text = element_text(size=12),
                                           axis.title = element_text(size=12), 
                                           legend.title = element_text(size=12),
                                           legend.text = element_text(size=12),
                                           legend.position = "none", 
                                           axis.title.y = element_blank()) + 
  ylim(0.5,1.4)  + 
  scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7))

sedcomp.beach_sorting <- sedcomp.beach_sorting + theme(axis.text = element_text(size=12),
                                                       axis.title = element_text(size=12), 
                                                       legend.title = element_text(size=12),
                                                       legend.text = element_text(size=12),
                                                       legend.position = "right", 
                                                       axis.title.y = element_blank()) +  
  ylim(0.5,1.4)  + 
  scale_fill_manual(values = alpha(c("#BE0032","#0067A5"),0.7))

FigS4_sorting <- plot_grid(sedprod_sorting, sedcomp_sorting, sedcomp.beach_sorting, ncol = 3, labels = c("A", "B", "C"))
FigS4_sorting
ggsave("figures/FigS4_sediment_sorting.jpg", width = 10, height = 4)
ggsave("figures/FigS4_sediment_sorting.pdf", width = 10, height = 4)




###############################
###########
## Fig S1: Plot effect sizes of atoll, habitat and nutrient status on coral cover and parrotfish biomass
####

##### coral cover
###

# plot random effects of Atoll
# Step 1: get expected predictions for *all observations*
# include atoll-level random effects (averaging over islands)
pred_draws_cover <- coral_cover_mod %>%
  add_epred_draws(
    newdata = reefbudget,
    re_formula = ~(1 | atoll)
  )

# Step 2: average posterior predictions within each atoll
pred_atoll_cover <- pred_draws_cover %>%
  group_by(.draw, atoll) %>%
  summarise(mean_epred = mean(.epred), .groups = "drop") %>%
  mutate(
    atoll_short = recode(
      atoll,
      "Great Chagos Bank" = "GCB",
      "Peros Banhos"      = "PB",
      "Salomon"           = "SA"
    )
  )

# Step 3: plot posterior distributions of mean predicted coral cover per atoll
(cover_atoll <- ggplot(pred_atoll_cover, aes(y = fct_rev(atoll_short), x = mean_epred)) +
    stat_halfeye(
      .width = c(0.5, 0.7, 0.9),
      fill = "#A593E0",
      alpha = 0.7
    ) +
    labs(
      x = "", #posterior mean per atoll "Predicted coral cover (%)"
      y = "", #title = "Posterior distribution of mean predicted coral cover across atolls"
    ) +
    xlim(-5,60) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), # remove gridlines
          panel.grid.minor = element_blank(), #remove gridlines
          strip.background = element_blank(), 
          legend.position = "none",
          rect = element_rect(fill = "transparent"),  
          plot.background = element_rect(fill = "transparent", color = NA),
          text=element_text(size=12,  family="sans"),
          legend.text=element_text(size=12, family = "sans"),
          axis.title.y = element_text(size=12, family = "sans"),
          axis.text = element_text(size = 12, family = "sans")) )

# atoll|island
re_island <- ranef(coral_cover_mod)$`atoll:island`
re_island_df <- as.data.frame(re_island[, , "Intercept"])
re_island_df$group <- rownames(re_island_df)

ggplot(re_island_df, aes(x = reorder(group, Estimate), y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = "Island (nested within atoll)",
    y = "Random intercept (difference from grand mean)",
    title = "Posterior differences across islands within atolls"
  )

# same for habitat & status
# Step 1: Summarise per draw, status, and habitat
pred_habitat_cover <- pred_draws_cover %>%
  group_by(.draw, status, habitat) %>%
  summarise(mean_epred = mean(.epred), .groups = "drop") %>%
  mutate(
    habitat = factor(habitat, levels = c("deep", "shallow", "lagoon")),
    status  = factor(status, levels = c("low", "high"))
  )

# Step 2: Define colours
high_colors <- c(
  "lagoon"  = adjustcolor("#0168A5", alpha.f = 0.4),
  "shallow" = adjustcolor("#0168A5", alpha.f = 0.8),
  "deep"    = adjustcolor("#004D7E", alpha.f = 0.7)
)

low_colors <- c(
  "lagoon"  = adjustcolor("#BF1F36", alpha.f = 0.4),
  "shallow" = adjustcolor("#BF1F36", alpha.f = 0.8),
  "deep"    = adjustcolor("#BF1F36", alpha.f = 1)
)

# Step 3: Add a new variable for fill colours depending on status + habitat
pred_habitat_cover <- pred_habitat_cover %>%
  mutate(
    fill_col = ifelse(
      status == "high",
      high_colors[as.character(habitat)],
      low_colors[as.character(habitat)]
    )
  )

# Step 4: Plot (Option A) with facets for each status
(cover_habitat <- ggplot(pred_habitat_cover, aes(y = habitat, x = mean_epred, fill = fill_col)) +
    stat_halfeye(
      .width = c(0.5, 0.7, 0.9),
      #alpha  = 1,
      colour = "black"
    ) +
    facet_wrap(~ status, ncol = 1) +
    scale_fill_identity() +
    xlim(-5,60) +
    labs(
      x = "Predicted coral cover (%)",
      y = ""
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      #legend.position = "none",
      rect = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      text = element_text(size = 12, family = "sans"),
      legend.text = element_text(size = 12, family = "sans"),
      axis.title.y = element_text(size = 12, family = "sans"),
      axis.text = element_text(size = 12, family = "sans")
    ) )


##### parrotfish biomass
###
pred_draws_parrot <- parrotfish_mod %>%
  add_epred_draws(
    newdata = reefbudget,
    re_formula = ~(1 | atoll)
  )

#average posterior predictions within each atoll
pred_atoll_parrot <- pred_draws_parrot %>%
  group_by(.draw, atoll) %>%
  summarise(mean_epred = mean(.epred), .groups = "drop")%>%
  mutate(
    atoll_short = recode(
      atoll,
      "Great Chagos Bank" = "GCB",
      "Peros Banhos"      = "PB",
      "Salomon"           = "SA"
    )
  )

# Step 3: plot posterior distributions of mean predicted parrotfish biomass per atoll
(parrot_atoll <- ggplot(pred_atoll_parrot, aes(y = fct_rev(atoll_short), x = mean_epred)) +
    stat_halfeye(
      .width = c(0.5, 0.7, 0.9),
      fill = "#A593E0",
      alpha = 0.7
    ) +
    labs(
      x = "", #posterior mean per atoll - expression("Predicted parrotfish biomass (kg ha"^-1*")")
      y = "", #title = "Posterior distribution of mean predicted parrotfish biomass across atolls"
    ) +
    xlim(-150,600)+
    theme_bw() + 
    theme(panel.grid.major = element_blank(), # remove gridlines
          panel.grid.minor = element_blank(), #remove gridlines
          strip.background = element_blank(), 
          legend.position = "none",
          rect = element_rect(fill = "transparent"),  
          plot.background = element_rect(fill = "transparent", color = NA),
          text=element_text(size=12,  family="sans"),
          legend.text=element_text(size=12, family = "sans"),
          axis.title.y = element_text(size=12, family = "sans"),
          axis.text = element_text(size = 12, family = "sans")) )

# atoll|island
re_island <- ranef(parrotfish_mod)$`atoll:island`
re_island_df <- as.data.frame(re_island[, , "Intercept"])
re_island_df$group <- rownames(re_island_df)

ggplot(re_island_df, aes(x = reorder(group, Estimate), y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = "Island (nested within atoll)",
    y = "Random intercept (difference from grand mean)",
    title = "Posterior differences across islands within atolls"
  )

# same for habitat & status
# Step 1: Summarise per draw, status, and habitat
pred_habitat <- pred_draws %>%
  group_by(.draw, status, habitat) %>%
  summarise(mean_epred = mean(.epred), .groups = "drop") %>%
  mutate(
    habitat = factor(habitat, levels = c("deep", "shallow", "lagoon")),
    status  = factor(status, levels = c("low", "high"))
  )

# Step 3: Add a new variable for fill colours depending on status + habitat
pred_habitat <- pred_habitat %>%
  mutate(
    fill_col = ifelse(
      status == "high",
      high_colors[as.character(habitat)],
      low_colors[as.character(habitat)]
    )
  )

# Step 4: Plot (Option A) with facets for each status
(parrot_habitat <- ggplot(pred_habitat, aes(y = habitat, x = mean_epred, fill = fill_col)) +
    stat_halfeye(
      .width = c(0.5, 0.7, 0.9),
      #alpha  = 1,
      colour = "black"
    ) +
    facet_wrap(~ status, ncol = 1) +
    scale_fill_identity() +
    xlim(-150,600) +
    labs(
      x = expression("Predicted parrotfish biomass (kg ha"^-1*")"),
      y = ""
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      #legend.position = "none",
      rect = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      text = element_text(size = 12, family = "sans"),
      legend.text = element_text(size = 12, family = "sans"),
      axis.title.y = element_text(size = 12, family = "sans"),
      axis.text = element_text(size = 12, family = "sans")
    ) )


#####
# combine atoll and habitat_status plots in rows
S1_atoll_plots <- plot_grid(
  cover_atoll, parrot_atoll,
  ncol = 2,
  labels = c("A", "B")
)

# Second row: plots 3 & 4 (double height)
S1_habitat_plots <- plot_grid(
  cover_habitat, parrot_habitat,
  ncol = 2,
  labels = c("C", "D")
)

#########
## plot nutrient effect sizes for coral cover

#get estimates and HPDs:
cover_pred_7<-coral_cover_mod.rg%>%
  emmeans(~ status|habitat,  
          type="response")%>%
  contrast("trt.vs.ctrl", level = .7)%>%
  as.data.frame()

cover_pred_9<-coral_cover_mod.rg%>%
  emmeans(~ status|habitat, 
          type="response")%>%
  contrast("trt.vs.ctrl", level = .9)%>%
  as.data.frame()

coral_cover_mod.rg %>% 
  emmeans(~ status|habitat,
          type = "response")%>%
  contrast("trt.vs.ctrl", level = .9)%>%
  as.data.frame()

#combine estimates:
cover_est<-data.frame(effect = c( 'lagoon', 'shallow', 'deep'),
                      estimate = c(cover_pred_7$estimate,cover_pred_7$estimate), #effect size is the same regardless of interval, so can use either for these
                      low = c (cover_pred_9$lower.HPD, cover_pred_7$lower.HPD),
                      up = c( cover_pred_9$upper.HPD, cover_pred_7$upper.HPD),
                      level = c( .9, .9, .9, .7, .7, .7))
cover_est

#horizontal
cover_effect_plot<-
  cover_est %>%
  mutate(effect = fct_relevel(effect, c("deep", "shallow", "lagoon")))%>%
  ggplot(aes(x = estimate, y = effect, color = as.factor(level)))+
  geom_vline(xintercept=0, lty=2, alpha = .8, color = "darkgrey")+
  geom_errorbar(aes(xmin = low, xmax = up), width= 0, linewidth = c(1, 1, 1, 2, 2, 2), alpha = .7)+
  geom_point(size =5, alpha = .9)+
  scale_color_manual(values = c("#0067A5", "#A1CAF1"))+ 
  scale_fill_manual(values = c("#0067A5", "#A1CAF1"))+
  #scale_y_discrete(labels = c("lagoon", "shallow", "deep"))+
  labs(x=expression("Seabird nutrient effect (%cover)"), y="")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(), 
        legend.position = "none",
        rect = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA),
        text=element_text(size=12,  family="sans"),
        legend.text=element_text(size=12, family = "sans"),
        axis.title.y = element_text(size=12, family = "sans"),
        axis.text = element_text(size = 12, family = "sans")) 
cover_effect_plot

#########
## plot nutrient effect sizes for parrotfish biomass

#get estimates and HPDs:
parrot_pred_7<-parrotfish_mod.rg%>%
  emmeans(~ status|habitat,  
          type="response")%>%
  contrast("trt.vs.ctrl", level = .7)%>%
  as.data.frame()

parrot_pred_9<-parrotfish_mod.rg%>%
  emmeans(~ status|habitat, 
          type="response")%>%
  contrast("trt.vs.ctrl", level = .9)%>%
  as.data.frame()

parrotfish_mod.rg %>% 
  emmeans(~ status|habitat,
          type = "response")%>%
  contrast("trt.vs.ctrl", level = .9)%>%
  as.data.frame()

#combine estimates:
parrot_est<-data.frame(effect = c( 'lagoon', 'shallow', 'deep'),
                       estimate = c(parrot_pred_7$estimate,parrot_pred_7$estimate), #effect size is the same regardless of interval, so can use either for these
                       low = c (parrot_pred_9$lower.HPD, parrot_pred_7$lower.HPD),
                       up = c( parrot_pred_9$upper.HPD, parrot_pred_7$upper.HPD),
                       level = c( .9, .9, .9, .7, .7, .7))
parrot_est

#horizontal
parrot_effect_plot<-
  parrot_est %>%
  mutate(effect = fct_relevel(effect, c("deep", "shallow", "lagoon")))%>%
  ggplot(aes(x = estimate, y = effect, color = as.factor(level)))+
  geom_vline(xintercept=0, lty=2, alpha = .8, color = "darkgrey")+
  geom_errorbar(aes(xmin = low, xmax = up), width= 0, linewidth = c(1, 1, 1, 2, 2, 2), alpha = .7)+
  geom_point(size =5, alpha = .9)+
  scale_color_manual(values = c("#0067A5", "#A1CAF1"))+ 
  scale_fill_manual(values = c("#0067A5", "#A1CAF1"))+
  #scale_y_discrete(labels = c("lagoon", "shallow", "deep"))+
  labs(x=expression("Seabird nutrient effect (kg ha"^-1*")"), y="")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(), 
        legend.position = "none",
        rect = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA),
        text=element_text(size=12,  family="sans"),
        legend.text=element_text(size=12, family = "sans"),
        axis.title.y = element_text(size=12, family = "sans"),
        axis.text = element_text(size = 12, family = "sans")) 


#### combine Figures
FigS1 <- plot_grid(
  cover_atoll, parrot_atoll,               # Row 1
  cover_habitat, parrot_habitat,           # Row 2
  cover_effect_plot, parrot_effect_plot,   # Row 3
  ncol = 2,
  labels = c("A","B","C","D","E","F"),
  align = "v",       # align vertically
  axis = "lr",       # align left and right axes
  rel_heights = c(1, 2, 1) # second row double height
)
FigS1

ggsave("figures/S1_cover_parrotfish.pdf", width = 10, height = 10)
ggsave("figures/S1_cover_parrotfish.jpg", width = 10, height = 10)


###############################
#########################
## Fig S2: Plot PCA for benthic community composition

# data frame
chagos.spe_coral <- benthic[,8:26] # specify taxa columns
chagos.spe_coral[chagos.spe_coral < 0.001] <- 0
chagos.spe_normcoral <- as.data.frame(sqrt(chagos.spe_coral[,])) # normalize data
#chagos.spe_normcoral[is.na(chagos.spe_normcoral)] <- 0 #replace NA with 0
chagos.env <- subset(benthic[,c(2,4,5)]) # atoll, habitat and status as grouping variable

#Run PCA
chagos.pca_coral<-rda(chagos.spe_normcoral, scale=TRUE)
chagos.pca_coral
summary(chagos.pca_coral)

#how many PC
screeplot(chagos.pca_coral, type="lines") #PC1-PC4 explain most variability

# which coral genera are driving differences among sites?
ef.pca<-envfit(chagos.pca_coral,chagos.spe_coral,permu=999) 
ef.pca #.... sign drive differences among sites (p<0.05*)

# PERMANOVA to assess differences between atolls, habitats and nutrient status
permanova_atoll <- adonis2(chagos.spe_coral~atoll, data=benthic)
permanova_atoll # atoll does significantly affect benthic comm comp p=0.001
posthoc_atoll <- pairwise.adonis(chagos.spe_coral, factors=benthic$atoll)
posthoc_atoll # SA vs GCB p=0.009
permanova_habitat <- adonis2(chagos.spe_coral~habitat, data=benthic)
permanova_habitat # habitat does significantly affect benthic comm comp p=0.007
posthoc_habitat <- pairwise.adonis(chagos.spe_coral, factors=benthic$habitat)
posthoc_habitat # lagoon vs deep p=0.045
permanova_status <- adonis2(chagos.spe_coral~status, data=benthic)
permanova_status # status does not significantly affect benthic comm comp

#Permdisp
permdisp <- betadisper(vegdist(chagos.spe_normcoral, method="bray", na.rm=TRUE), benthic$atoll) 
par(mfrow=c(1,2))
plot(permdisp, main="PCoA")
boxplot(permdisp, main="Distance to centroids")
anova(permdisp) # p>0.05 for atoll BUT p=0.009 for habitat, meaning group dispersion is NOT homogenous 
# this means differences between centroids might not solely be due to differences in composition for habitats

# plot PCA with habitat as grouping
chagos.pca_coral <- PCA(chagos.spe_normcoral, graph=FALSE)
p1 <- fviz_pca_biplot(chagos.pca_coral,
                      #axes = c(3, 4), # plot PC3 vs PC4
                      habillage = benthic$atoll,
                      addEllipses = TRUE, ellipse.level = 0.7,
                      label = "var", col.var = "black", repel = TRUE,
                      select.var = list(contrib=10), # only display 10 most contributing taxa
                      title="")
(p1 <- p1 + scale_color_manual(values=c("brown", "orange", "gold")) + 
    scale_fill_manual(values=c("brown", "orange", "gold")) + 
    theme_classic() + theme(legend.position="top") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))) 

p2 <- fviz_pca_biplot(chagos.pca_coral,
                      axes = c(3, 4), # plot PC3 vs PC4
                      habillage = benthic$habitat,
                      addEllipses = TRUE, ellipse.level = 0.7,
                      label = "var", col.var = "black", repel = TRUE,
                      select.var = list(contrib=10), # only display 10 most contributing taxa
                      title="")
(p2 <- p2 + scale_color_manual(values=c("#F564E3","#00BA38", "#619CFF")) +
    scale_fill_manual(values=c("#F564E3","#00BA38", "#619CFF")) + 
    theme_classic() + theme(legend.position="top") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))) 

p3 <- fviz_pca_biplot(chagos.pca_coral,
                      #axes = c(3, 4), # plot PC3 vs PC4
                      habillage = benthic$status,
                      addEllipses = TRUE, ellipse.level = 0.7,
                      label = "var", col.var = "black", repel = TRUE,
                      select.var = list(contrib=10), # only display 10 most contributing taxa
                      title="")
(p3 <- p3 + scale_color_manual(values=c("#0067A5", "#BE0032")) + 
    scale_fill_manual(values=c("#0067A5", "#BE0032")) +
    theme_classic() + theme(legend.position="top") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))) 

FigS2_comm_comp <- plot_grid(p1,p2,p3, ncol=2, labels = c("A", "B", "C"))
FigS2_comm_comp
ggsave("figures/chagos_community_composition.jpg", width = 10, height = 10)
ggsave("figures/chagos_community_composition.pdf", width = 10, height = 10)


##########################
##################
## Fig S3: Plot PCA for sediment composition

#dataframe
sediment <- read.csv("SB_sedimentcomp_2021.csv", check.names = FALSE)%>%
  mutate(status = as.factor(status),
         habitat = as.factor(habitat),
         atoll = as.factor(atoll))

sediment_pca <- sediment %>%
  filter(comparison == "y" |
           (habitat == "beach"))


sediment_pca$habitat <- factor(sediment_pca$habitat, 
                               levels = c("beach", "lagoon", "shallow", "deep"))

chagos.spe_sediment <- sediment_pca[,27:35] # specify taxa columns
chagos.spe_sediment[chagos.spe_sediment < 0.001] <- 0
chagos.spe_normsediment <- as.data.frame(sqrt(chagos.spe_sediment[,])) # normalize data
#chagos.spe_normsediment[is.na(chagos.spe_normsediment)] <- 0 #replace NA with 0
chagos.env_sediment <- subset(sediment_pca[,c(2,5,8)]) # atoll, habitat and status as grouping variable

#Run PCA
chagos.pca_sediment<-rda(chagos.spe_normsediment, scale=TRUE)
chagos.pca_sediment
summary(chagos.pca_sediment)

#how many PC
screeplot(chagos.pca_sediment, type="lines") #PC1-PC3 explain most variability

# which sediment genera are driving differences among sites?
ef.pca_sediment<-envfit(chagos.pca_sediment,chagos.spe_sediment,permu=999) 
ef.pca_sediment #all sign drive differences among sites (p<0.05*)

# PERMANOVA to assess differences between atolls, habitats and nutrient status
permanova_atoll_sed <- adonis2(chagos.spe_sediment~atoll, data=sediment_pca)
permanova_atoll_sed # atoll does significantly affect sediment comm comp p=0.005
posthoc_atoll_sed <- pairwise.adonis(chagos.spe_sediment, factors=sediment_pca$atoll)
posthoc_atoll_sed # SA vs GCB p=0.003
permanova_habitat_sed <- adonis2(chagos.spe_sediment~habitat, data=sediment_pca)
permanova_habitat_sed # habitat ns
posthoc_habitat_sed <- pairwise.adonis(chagos.spe_sediment, factors=sediment_pca$habitat)
posthoc_habitat_sed # lagoon vs deep p=0.054
permanova_status_sed <- adonis2(chagos.spe_sediment~status, data=sediment_pca)
permanova_status_sed # status does not significantly affect sediment comm comp

#Permdisp
permdisp_sed <- betadisper(vegdist(chagos.spe_normsediment, method="bray", na.rm=TRUE), sediment_pca$atoll) 
par(mfrow=c(1,2))
plot(permdisp_sed, main="PCoA")
boxplot(permdisp_sed, main="Distance to centroids")
anova(permdisp_sed) # p>0.05, meaning group dispersion is  homogenous 
# this means differences between centroids are solely due to differences in composition

# plot PCA with habitat as grouping
chagos.pca_sediment <- PCA(chagos.spe_normsediment, graph=FALSE)
p1_sed <- fviz_pca_biplot(chagos.pca_sediment,
                          #axes = c(3, 4), # plot PC3 vs PC4
                          habillage = sediment_pca$atoll,
                          addEllipses = TRUE, ellipse.level = 0.7,
                          label = "var", col.var = "black", repel = TRUE,
                          select.var = list(contrib=10), # only display 10 most contributing taxa
                          title="")
(p1_sed <- p1_sed + scale_color_manual(values=c("brown", "orange", "gold")) + 
    scale_fill_manual(values=c("brown", "orange", "gold")) + 
    theme_classic() + theme(legend.position="top") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))) 

p2_sed <- fviz_pca_biplot(chagos.pca_sediment,
                          axes = c(2, 3), # plot PC1 vs PC3
                          habillage = sediment_pca$habitat,
                          addEllipses = TRUE, ellipse.level = 0.7,
                          label = "var", col.var = "black", repel = TRUE,
                          select.var = list(contrib=10), # only display 10 most contributing taxa
                          title="")
(p2_sed <- p2_sed + scale_color_manual(values=c("gold", "#F564E3","#00BA38", "#619CFF")) +
    scale_fill_manual(values=c("gold", "#F564E3","#00BA38", "#619CFF")) + 
    theme_classic() + theme(legend.position="top") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))) 

p3_sed <- fviz_pca_biplot(chagos.pca_sediment,
                          axes = c(1, 2), # plot PC3 vs PC4
                          habillage = sediment_pca$status,
                          addEllipses = TRUE, ellipse.level = 0.7,
                          label = "var", col.var = "black", repel = TRUE,
                          select.var = list(contrib=10), # only display 10 most contributing taxa
                          title="")
(p3_sed <- p3_sed + scale_color_manual(values=c("#0067A5", "#BE0032")) + 
    scale_fill_manual(values=c("#0067A5", "#BE0032")) +
    theme_classic() + theme(legend.position="top") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))) 

FigS3_sed_comp <- plot_grid(p1_sed,p2_sed,p3_sed, ncol=2, labels = c("A", "B", "C"))
FigS3_sed_comp
ggsave("figures/chagos_sediment_composition.jpg", width = 10, height = 10)
ggsave("figures/chagos_sediment_composition.pdf", width = 10, height = 10)


