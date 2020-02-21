##############################
## NJ Tree
##
## Matt Brachmann (PhDMattyB)
##
## 2019-01-31
##
##############################

## DA STARTUP ####
setwd('~/PhD/SNP Demographic modelling/Working Directory/')

library(tidyverse)
library(poppr)
library(adegenet)
library(ape)
# library(mmod)
# library(treemap)
# library(hierfstat)
# library(phytools)
# library(pcadapt)
library(ggtree)
# library(treeio)

theme_set(theme_bw())

#IMPORT THE GENOTYPIC DATA
## DA DATA ####
SNP_iceland = read_tsv('Nov052018_plink_input_icelandic_matt.ped',
                       col_names = T) %>%
  select(-PaternalID, -MaternalID, -Sex, -Phenotype)

SNP_iceland = mutate(.data = SNP_iceland,
              POP_name = as.factor(case_when(
                `#FamilyID` == "1" ~ "T.LGB",
                `#FamilyID` == '2' ~ 'V.BR',
                `#FamilyID` == '3' ~ 'V.SIL',
                `#FamilyID` == '4'~ 'S.PL',
                `#FamilyID` == '5' ~ 'S.PI',
                `#FamilyID` == '6' ~ 'T.PL',
                `#FamilyID` == '7' ~ 'T.SB',
                `#FamilyID` == '8' ~ 'S.LGB',
                `#FamilyID` == '9' ~ 'G.SB',
                `#FamilyID` == '10' ~ 'G.PI',
                `#FamilyID` == '11' ~ 'Mjoavatn',
                `#FamilyID` == '12' ~ 'Fljotaa')))

SNP_iceland = SNP_iceland %>% 
  select(IndividualID, POP_name, contains("Affx-"))
 Classifiers = SNP_iceland %>% select(POP_name)

SNP_iceland = as.data.frame(SNP_iceland)

## MAKE DA TREE #####
iceland_genind2 = df2genind(SNP_iceland[,3:length(SNP_iceland)], 
                            ploidy = 2, 
                           ind.names = SNP_iceland[,1], 
                           sep = '\t',
                           strata = Classifiers)

set.seed(666)
nameStrata(iceland_genind2) = ~LaMorph

CHARR_TREE = iceland_genind2 %>% 
  genind2genpop(pop = ~LaMorph) %>% 
  aboot(cutoff = 50,
        quiet = FALSE,
        sample = 1000,
        distance = nei.dist,
        tree = 'upgma')

# str(CHARR_TREE)
write.tree(CHARR_TREE, file = "CHARR_TREE_SNPS")


## READ DA TREE ####
FISHTREE = read.tree('CHARR_TREE_SNPS')

## GGTREE #####
FISHY = as_tibble(FISHTREE)
FISHY_data = mutate(.data = FISHY,
              label2 = as.factor(case_when(
                label == "G.SB" ~ "Galtabol - Small benthic",
                label == 'G.PI' ~ 'Galtabol - Piscivorous',
                label == 'T.LGB' ~ 'Thingvallavatn - Large benthic',
                label == 'T.SB'~ 'Thingvallavatn - Small benthic',
                label == 'T.PL' ~ 'Thingvallavatn - Planktivorous',
                label == 'S.LGB' ~ 'Svinavatn - Large benthic',
                label == 'S.PI' ~ 'Svinavatn - Piscivorous',
                label == 'S.PL' ~ 'Svinavatn - Planktivorous',
                label == 'V.BR' ~ 'Vatnshlidarvatn - Brown',
                label == 'V.SIL' ~ 'Vatnshlidarvatn - Silver',
                label == 'Mjoavatn' ~ 'Mjoavatn',
                label == 'Fljotaa' ~ 'Fljotaa')))

FISHY_data = FISHY_data %>% rename(label = 5, label_morph = 4) %>% 
  select(parent, node, branch.length, label)

info = read_csv('ggtree_info.csv') %>% 
  rename(label = LAMORPH) %>% 
  select(label, MORPH, node)
data = mutate(.data = info,
                     label2 = as.factor(case_when(
                       label == "G.SB" ~ "Galtabol - Small benthic",
                       label == 'G.PI' ~ 'Galtabol - Piscivorous',
                       label == 'T.LGB' ~ 'Thingvallavatn - Large benthic',
                       label == 'T.SB'~ 'Thingvallavatn - Small benthic',
                       label == 'T.PL' ~ 'Thingvallavatn - Planktivorous',
                       label == 'S.LGB' ~ 'Svinavatn - Large benthic',
                       label == 'S.PI' ~ 'Svinavatn - Piscivorous',
                       label == 'S.PL' ~ 'Svinavatn - Planktivorous',
                       label == 'V.BR' ~ 'Vatnshlidarvatn - Brown',
                       label == 'V.SIL' ~ 'Vatnshlidarvatn - Silver',
                       label == 'Mjoavatn' ~ 'Mjoavatn',
                       label == 'Fljotaa' ~ 'Fljotaa')))

data = mutate(.data = data,
              label3 = as.factor(case_when(
                label == "G.SB" ~ "G: Benthic",
                label == 'G.PI' ~ 'G: Pelagic',
                label == 'T.LGB' ~ 'T: Benthic 1',
                label == 'T.SB'~ 'T: Benthic 2',
                label == 'T.PL' ~ 'T: Pelagic',
                label == 'S.LGB' ~ 'S: Benthic',
                label == 'S.PI' ~ 'S: Pelagic 2',
                label == 'S.PL' ~ 'S: Pelagic 1',
                label == 'V.BR' ~ 'V: Pelagic',
                label == 'V.SIL' ~ 'V: Benthic',
                label == 'Mjoavatn' ~ 'Mjoavatn',
                label == 'Fljotaa' ~ 'Fljotaa')))

data = data %>%
  rename(label = 4, 
         label_morph = 1) 
  
Joined_data = left_join(FISHY, info, by = 'label')
Joined_data = left_join(FISHY_data, data, by = 'label') %>% 
  select(-MORPH)

Joined_data = Joined_data %>% 
  rename(node = 2, 
         large_label = 4, 
         label = 7)
phylo_data = as.phylo(Joined_data)
str(phylo_data)

FISHY_TREE = ggtree(phylo_data, ladderize = T) +
  #geom_text(aes(label= boot), hjust=-.3)+
  geom_text2(aes(label = label,
                 subset = !is.na(as.numeric(label)) & 
                   as.numeric(label) > 50),
             vjust = 1.2,
             hjust = 1.2)+
  geom_tiplab(align = TRUE, 
              hjust = 0.4, 
              vjust = -0.3)+
  geom_hilight(node = 14, 
               fill = 'steelblue',
               alpha = 0.2,
               extend = 0.0009)+
  geom_hilight(node = 23, 
               fill = 'darkgreen',
               alpha = 0.2,
               extend = 0.0009)+
  geom_hilight(node = 18, 
               fill = '#FF8300',
               alpha = 0.2,
               extend = 0.0009) +
  geom_hilight(node = 21, 
               fill = '#14CC9C', 
               alpha = 0.2,
               extend = 0.0009)+
  theme_tree2()

ggsave(plot = last_plot(), 
       'Icelandic_Charr_Phylogeny_20.02.2020.tiff')