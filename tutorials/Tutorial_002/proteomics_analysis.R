###############################################################################
# SETUP

tutorial_dir <- 'D:/Daten/Work/Dublin/CloudStation/other/2023-03-30_RSeminarProteomics/Tutorial_002'
# source(setup.R)
# set_up(tutorial_dir)

library(tidyverse)
script_dir <- str_glue('{tutorial_dir}/analysis_scripts/')
setwd(tutorial_dir)

# import analysis functions
list.files(script_dir, '*.R', full.names = T) %>%
  map(function(x) {
    source(x)
  })

###############################################################################
# IMPORT DATA

annotation <- read.table('data/annotation.csv', header=T, sep=';',
						 fileEncoding='UTF-8-BOM') %>%
  mutate(label_new = str_c(group, biol_repl, tech_repl, sep='__'))

# parses output from MaxQuant
raw <- read.table('data/proteinGroups.txt', header=T, sep='\t') %>%
  filter(Only.identified.by.site != '+') %>% 
  filter(Reverse != '+') %>%
  filter(Potential.contaminant != '+') %>%
  select(starts_with('LFQ'), 'Majority.protein.IDs', 'id') %>%
  pivot_longer(cols=starts_with('LFQ'), 
               values_to='LFQ', names_to = 'label') %>%
  filter(LFQ != 0) %>%
  mutate(id = as.character(id)) %>%
  inner_join(annotation, by='label')


# assemble annotation and raw data to data object
df <- proteomics_data(
  raw, annotation,
  has_tech_repl = TRUE, is_log2 = FALSE,
  df_label = 'label_new', annotation_label = 'label_new', 
  annotation_group = 'group')


###############################################################################
# QUALITY CONTROL

vis_qc_histo(df)
vis_qc_count(df)

vis_upset(df %>% collapse_tech_repl)

vis_pca(df %>% filtering('each', 0.6))

###############################################################################
# OUTLIER REMOVAL

# remove ripa group
df <- df %>%
  remove_groups(c('ripa'))

# visualize data quality control
vis_qc_count(df)
vis_pca(df %>% filtering('each', 0.6))

###############################################################################
# ID MAPPING

# requires database from HGNC genename 
# https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
hgnc_database <- str_c(tutorial_dir, '/data/hgnc_complete_set.txt')

id_mapping <- raw %>%
  select(id, Majority.protein.IDs) %>%
  rename(mpi = Majority.protein.IDs) %>%
  distinct %>%
  mutate(uniprot = collapse_uniprot_ids(mpi, 'first'),
         hgnc = map_uniprot_hgnc(uniprot, hgnc_database, 'first'),
         id = make.names(id))

str(id_mapping)

id_mapping %>%  select(id, hgnc) %>%
  distinct %>% 
  pull(hgnc) %>% duplicated %>% table %>% print

###############################################################################
# DATA FILTERING
filtering_methods()

df_filtered <- df %>%
  filtering('each', 0.6) %>%
  collapse_tech_repl

str(df_filtered)

###############################################################################
# DATA IMPUTATION

imputation_methods()

df_imput_1 <- df_filtered %>%
  imputation('zero')

df_imput_2 <- df_filtered %>%
  imputation('MinProb')

df_imput_3 <- df_filtered %>%
  imputation('mixed_sample', mar_method = 'MLE', mnar_method = 'MinProb')

vis_qc_histo(df_imput_3)

vis_pca(df_imput_3)

###############################################################################
# DIFFERENTIAL ANALYSIS

stat_analysis_methods()

cntrsts <- construct_contrasts_control(df_imput_3, 'NT')
construct_contrasts_all(df_imput_3)
str(cntrsts)
# cntrsts is expected to be a dataframe with two columns named a and b

df_diff <- df_imput_3 %>%
  diff_expr(method = 'limma', 
          contrasts = cntrsts) %>%
  inner_join(diff_type(df_filtered, contrasts = cntrsts),
             by = c('contrast', 'id')) %>%
  filter(diff_type != 'imput_imput') %>%
  mutate(pval_adj = p.adjust(pval, 'BH'))

str(df_diff)

###############################################################################
# VOLCANO PLOT

df_diff %>%
  filter(contrast == 'G12D_DMOG - NT') %>% 
  ggplot(aes(x = logfc, y = -log10(pval_adj))) + 
    geom_point()

df_diff %>%
  filter(contrast == 'G12D_DMOG - NT') %>% 
  filter(diff_type == 'value_value') %>%
  arrange(pval_adj) %>% 
  head(10) %>%
  inner_join(id_mapping %>% select(id, hgnc), by='id')
