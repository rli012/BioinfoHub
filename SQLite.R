library(RSQLite)

db = dbConnect(SQLite(), dbname="TCGA.sqlite")
dbListTables(db)

dbTable <- 'TCGA-PRAD'
metadata <- dbReadTable(db, dbTable)
