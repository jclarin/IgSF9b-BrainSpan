library(readxl)
library(tidyverse)
library(car)

setwd("C:/Users/Gaolab/Desktop/IgSF9b/IgSF9b Dif Expression Analysis/Brainspan")

#All Brain Regions

all_gene_expression <- read.csv("Whole Dataset/expression_matrix.csv", header = FALSE)
all_gene_expression <- subset(all_gene_expression, select = -c(1))
all_sample_info <- read.csv("Whole Dataset/columns_metadata.csv")
all_gene_info <- read.csv("Whole Dataset/rows_metadata.csv")

#Heatmaps on the BrainSpan website are in log2RPKM and representative of exon level data, not gene level.
#all_gene_expression is 'gene composite model' data from which I subset IgSF9b expression data.
#all_gene_expression values are in RPKM that must be converted to TPM to compare expression across samples.

#IgSF9b

IgSF9b_expression <- all_sample_info
IgSF9b_expression$RPKM <- as.numeric(all_gene_expression[1528,]) 

#Created age categories based both upon BrainSpan documentation (prenatal) and Singh et al. 2022 (postnatal)

IgSF9b_expression$AgeCat <- car::recode(IgSF9b_expression$age,
                      "c('8 pcw', '9 pcw', '12 pcw') = 'Early-Prenatal' ;
                       c('13 pcw', '16 pcw', '17 pcw') = 'Early-mid-Prenatal' ;
                       c('19 pcw', '21 pcw', '24 pcw') = 'Late-mid-Prenatal' ;
                       c('25 pcw', '26 pcw','35 pcw', '37 pcw') = 'Late-Prenatal' ;
                       c('4 mos', '10 mos', '1 yrs') = 'Infancy' ;
                       c('2 yrs', '3 yrs', '4 yrs', '8 yrs', '11 yrs') = 'Childhood' ;
                       c('13 yrs', '15 yrs', '18 yrs', '19 yrs') = 'Adolescence' ;
                       c('21 yrs', '23 yrs', '30 yrs', '36 yrs', '37 yrs', '40 yrs') = 'Adulthood' " )


IgSF9b_expression$RegionCat <- car::recode(IgSF9b_expression$structure_name,
                       "c('dorsolateral prefrontal cortex', 'ventrolateral prefrontal cortex', 'anterior (rostral) cingulate (medial prefrontal) cortex', 'orbital frontal cortex') = 'Prefrontal Regions' ;
                        c('primary motor cortex (area M1, area 4)', 'primary somatosensory cortex (area S1, areas 3,1,2)', 'primary auditory cortex (core)', 'primary visual cortex (striate cortex, area V1/17)') = 'Sensory/Motor Regions' ;
                        c('posteroventral (inferior) parietal cortex', 'posterior (caudal) superior temporal cortex (area 22c)', 'inferolateral temporal cortex (area TEv, area 20)') = 'Parietal/Temporal Areas' ;
                        c('hippocampus (hippocampal formation)', 'amygdaloid complex', 'striatum', 'mediodorsal nucleus of thalamus') = 'Subcortical Regions' ;
                        c('cerebellar cortex') = 'Cerebellar Cortex' " )
                                           

IgSF9b_expression <- IgSF9b_expression %>% 
  mutate(TPM = (RPKM/colSums(all_gene_expression) * 10^6))  #Convert from RPKM to TPM for each sample, gene length is already accounted for in RPKM calculation

#Average TPM Across All Brain Regions for each donor
IgSF9b_expression$donor_name <- as.factor(IgSF9b_expression$donor_name)
IgSF9b_expression <- group_by(IgSF9b_expression, donor_name) %>% mutate(TPM_Avg = mean(TPM))
Brain_Avg <- unique(IgSF9b_expression [, c(2,3,10,1)]) #all unique rows from this subset of columns
IgSF9b_expression <- IgSF9b_expression %>% select(-c(12))



#Average TPM for each Region, for each Age Category

Regional_Expression <- data.frame(IgSF9b_expression$structure_name,IgSF9b_expression$AgeCat,IgSF9b_expression$TPM, IgSF9b_expression$RegionCat)
names(Regional_Expression) <- c("structure_name", "AgeCat", "TPM", "RegionCat")
Regional_Expression <- unique(Regional_Expression %>% group_by(structure_name, AgeCat) %>%
   mutate(TPM_Avg = mean(TPM)) %>% select(-c(3)))

Regional_Expression <- group_by(Regional_Expression, RegionCat, AgeCat) %>% mutate(regionCat_Avg = mean(TPM_Avg))

#Filter rows based on whether the brain region is present in adult
Adult_Regions = list(Regional_Expression[Regional_Expression$AgeCat == "Adulthood",][,1])
Regional_Expression <- Regional_Expression %>%
  filter(structure_name %in% Adult_Regions[[1]][["structure_name"]])



#Create subset for just prefrontal regions

IgSF9b_PFC <- IgSF9b_expression %>% 
  filter(structure_acronym %in% c('DFC','VFC','MFC','OFC'))
IgSF9b_PFC <- group_by(IgSF9b_PFC, structure_acronym, AgeCat) %>% mutate(TPM_Avg = mean(TPM))



###All regions plots
all_plot <- ggplot(IgSF9b_expression, 
            aes(factor(AgeCat, 
                level= c('Early-Prenatal', 'Early-mid-Prenatal', 
                'Late-mid-Prenatal', 'Late-Prenatal', 'Infancy', 
                'Childhood', 'Adolescence', 'Adulthood')), TPM))

all_plot_avg <- ggplot(Brain_Avg, 
            aes(factor(AgeCat, 
              level= c('Early-Prenatal', 'Early-mid-Prenatal', 
              'Late-mid-Prenatal', 'Late-Prenatal', 'Infancy', 
              'Childhood', 'Adolescence', 'Adulthood')), TPM_Avg  ))

all_plot + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.2) + 
  theme_classic() + stat_summary(fun=median, geom = "line", aes(group=1)) +
  labs(title = "Whole-Brain Expression", x = "", y = "IgSF9b Expression (TPM)")

all_plot + geom_boxplot() + stat_summary(fun=mean, geom = "line", aes(group=1)) +
  labs(title = "Whole-Brain Expression", x = "", y = "IgSF9b Expression (TPM)") + theme_classic()

all_plot_avg + geom_boxplot() + stat_summary(fun=mean, geom = "line", aes(group=1)) +
  labs(title = "Whole-Brain Expression", x = "", y = "IgSF9b Expression (TPM)") + theme_classic()

all_plot_avg + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) + 
  theme_classic() + stat_summary(fun=median, geom = "line", aes(group=1)) + 
  labs(title = "Whole-Brain Expression", x = "", y = "IgSF9b Expression (TPM)")


###PFC Plots

PFC_plot <- ggplot(IgSF9b_PFC, 
            aes(factor(AgeCat, 
            level= c('Early-Prenatal', 'Early-mid-Prenatal', 
            'Late-mid-Prenatal', 'Late-Prenatal', 'Infancy', 
            'Childhood', 'Adolescence', 'Adulthood')), TPM, fill = structure_acronym))

PFC_plot_Avg <- ggplot(IgSF9b_PFC, 
            aes(factor(AgeCat, 
            level= c('Early-Prenatal', 'Early-mid-Prenatal', 
            'Late-mid-Prenatal', 'Late-Prenatal', 'Infancy', 
            'Childhood', 'Adolescence', 'Adulthood')), TPM_Avg, color = structure_acronym, group = structure_acronym))

PFC_plot + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, position = "dodge") + theme_classic() +
  stat_summary(fun=median, geom = "line", aes(group=1)) +
  labs(title = "PFC Expression", x = "", y = "IgSF9b Expression (TPM)")

PFC_plot_Avg + geom_point() + geom_line() + theme_classic() +
  labs(x = '', y = 'IgSF9b Expression (TPM)')

PFC_plot + geom_boxplot() + stat_summary(fun=median, geom = "line", aes(group=1)) +
  labs(title = "PFC Expression", x = "", y = "IgSF9b Expression (TPM)") + 
  guides(fill= guide_legend(title = "Region"))  + theme_classic()

#Regional Average Plots

Regional_plot <- ggplot(Regional_Expression, 
                 aes(factor(AgeCat, 
                 level= c('Early-Prenatal', 'Early-mid-Prenatal', 
                 'Late-mid-Prenatal', 'Late-Prenatal', 'Infancy', 
                 'Childhood', 'Adolescence', 'Adulthood')), TPM_Avg, fill = structure_name, group = structure_name, color = structure_name))

RegionCat_Regional_plot <- ggplot(Regional_Expression, 
                           aes(factor(AgeCat, 
                           level= c('Early-Prenatal', 'Early-mid-Prenatal', 
                           'Late-mid-Prenatal', 'Late-Prenatal', 'Infancy', 
                           'Childhood', 'Adolescence', 'Adulthood')), regionCat_Avg, group = RegionCat, color = RegionCat))

Regional_plot + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) + theme_classic()

Regional_plot + geom_point() + geom_line() + theme_classic() +
  labs(title = "Regional Expression of IgSF9b Across Development", x = "", y = "IgSF9b Expression (TPM)") 


RegionCat_Regional_plot + geom_point() + geom_line() + theme_classic() +
  labs(title = "Regional Expression of IgSF9b Across Development", x = "", y = "IgSF9b Expression (TPM)") +
  guides(fill = guide_legend(title = "Region Category"))

###Adult Plots

Adult_plot <- IgSF9b_expression %>% filter(AgeCat == 'Adulthood') %>% 
  ggplot(., aes(reorder(structure_name, -TPM), TPM, fill = RegionCat)) 

Adult_plot_RegionCat <- IgSF9b_expression %>% filter(AgeCat == 'Adulthood') %>% 
  ggplot(., aes(reorder(RegionCat, -TPM), TPM, fill = RegionCat))

Adult_plot + geom_bar(stat = "summary", fun = "mean") +
  labs(title = "IgSF9b Expression Across the Adult Brain") + theme_classic()

Adult_plot_RegionCat + geom_bar(stat = "summary", fun = "mean") +
  labs(title = "IgSF9b Expression Across the Adult Brain", x = "", y = "IgSF9b Expression (TPM)") + 
  guides(fill= guide_legend(title = "Region Category")) + theme_classic()






