################################################################################
# SuREViz Main Query Function
# ---------------------------
# This R script defines a function, suREVizClinVarQuery, to query the SuREViz.ClinVar
# SQLite database for genomic variants within a specified region. It utilizes the
# 'DBI' and 'RSQLite' libraries for database connectivity.
#
# Function Usage:
# result <- suREVizClinVarQuery(chr, pos, flank, db_path)
#
# Parameters:
#   chr      : Chromosome
#   pos      : Genomic position
#   flank    : Flanking area to be viewed
#   db_path  : Path to SQLite database
#
# Example:
# result <- suREVizClinVarQuery("chr19", 56252947, 50000, "/path/to/database.sqlite")
################################################################################

suREVizClinVarQuery <- function(chr, pos, ref, alt, db_path) {
  
  # Connect to SQLite database
  con <- dbConnect(RSQLite::SQLite(), dbname = db_path)
  
  # Build and execute the SQL query
  query <- sprintf(
      "SELECT * FROM \"SuREViz.ClinVar\" WHERE
      CHROM = '%s' AND POS = %d AND REF = '%s' AND ALT = '%s'",
      chr, pos , ref, alt
  )
  
  result <- dbGetQuery(con, query)
  
  # Close the database connection
  dbDisconnect(con)
  
  return(result)
  
  }
