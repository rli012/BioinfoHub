library(NACHO)

rccData <- load_rcc(
  data_directory = 'data/RAW/', # Where the data is
  ssheet_csv = targets, # The samplesheet
  id_colname = "FileName", # Name of the column that contains the unique identfiers
  housekeeping_genes = NULL, # Custom list of housekeeping genes
  housekeeping_predict = FALSE, # Whether or not to predict the housekeeping genes
  normalisation_method = "GEO", # Geometric mean or GLM
  n_comp = 5 # Number indicating how many principal components should be computed. 
)

visualise(rccData)
