set_up <- function(project_dir) {
  ############################################################################
  ## set up dirs
  if (!dir.exists(project_dir)) {
    dir.create(project_dir, recursive = TRUE)
  }
  setwd(project_dir)
  
  if (!dir.exists('data')) {
    dir.create('data')
  }
  if (!dir.exists('analysis_scripts')) {
    dir.create('analysis_scripts')
  }
  
  
  ############################################################################
  ## set up packages
  packages_cran <- c(
    'tidyverse',
    'curl',
    'BiocManager',
    'patchwork',
    'UpSetR',
    'imputeLCMD'
  )
  
  packages_bioconductor <- c(
    'limma',
    'ComplexHeatmap',
    'MsCoreUtils',
    'impute',
    'pcaMethods'
  )
  
  install.packages(packages_cran)
  BiocManager::install(packages_bioconductor, update=FALSE, ask=FALSE)
  
  
  ############################################################################
  ## set up files
  # download HGNC dataset
  curl::curl_download(
    'https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt',
    'data/hgnc_complete_set.txt')
  
  # download MS data
  curl::curl_download(
    'https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/02/PXD039404/proteinGroups.txt',
    'data/proteinGroups.txt'
  )
  
  # download annotation file from GitHub
  curl::curl_download(
    'https://raw.githubusercontent.com/PhilippJunk/SBI_Tutorials/main/tutorials/Tutorial_002/annotation.csv',
    'data/annotation.csv'
  )
   
  # download main analysis script from GitHub
  curl::curl_download(
    'https://raw.githubusercontent.com/PhilippJunk/SBI_Tutorials/main/tutorials/Tutorial_002/proteomics_analysis.R',
    'proteomics_analysis.R'
  )
  
  # download analysis scripts
  scripts <- c('analysis_obj.R', 'diff_expr.R', 'filtering.R', 'id_mapping.R', 'imputation.R', 'visualization.R')
  for (script in scripts) {
    url <- stringr::str_glue(
      'https://raw.githubusercontent.com/philippjunk/proteomics_analysis/9074d22d56b5593464716aed9d54f9b928bd5ab1/{script}'
    )
    curl::curl_download(url, stringr::str_glue('analysis_scripts/{script}'))
  }
  
}
