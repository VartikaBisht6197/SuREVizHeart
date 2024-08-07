################################################################################
# SuREViz Query if chr:pos in SuRE DB
# ---------------------------
# This R script defines a function, inSuRE, to query the SuREViz.Main
# SQLite database for genomic variants within a specified region. It utilizes the
# 'DBI' and 'RSQLite' libraries for database connectivity.
#
# Function Usage:
# result <- suREVizMainQuery(chr, pos, flank, db_path)
#
# Parameters:
#   chr      : Chromosome
#   pos      : Genomic position
#   flank    : Flanking area to be viewed
#   db_path  : Path to SQLite database
#
# Example:
# result <- suREVizMainQuery("chr19", 56252947, 50000, "/path/to/database.sqlite")
################################################################################

# Libraries
library(DBI)
library(RSQLite)

# Main function to query the SuREViz.Main SQLite database
inSuRE <- function(chr, pos, db_path) {
  # Connect to SQLite database
  con <- dbConnect(RSQLite::SQLite(), dbname = db_path)
  
  # Build and execute the SQL query
  query <- sprintf(
    "SELECT * FROM \"SuREViz.Main\" WHERE
    CHROM = '%s' AND POS = %d",
    chr, pos 
  )
  
  result <- dbGetQuery(con, query)
  # Close the database connection
  dbDisconnect(con)
  
  # Return the result
  return(result)
}
