################################################################################
# SuREViz Gene Query Function
# ---------------------------
# This R script defines a function, suREVizGeneQuery, to query the SuREViz.Gene
# SQLite database for genes within a specified genomic region. It utilizes the
# 'DBI' and 'RSQLite' libraries for database connectivity.
#
# Function Usage:
# result <- suREVizGeneQuery(chr, pos, flank, db_path, outfile)
#
# Parameters:
#   chr      : Chromosome
#   pos      : Genomic position
#   flank    : Flanking area to be viewed
#   db_path  : Path to SQLite database
#
# Example:
# result <- suREVizGeneQuery("chr19", 56252947, 50000, "/path/to/database.sqlite")
################################################################################

# Libraries
library(DBI)
library(RSQLite)

# Main function to query the SuREViz.Gene SQLite database
suREVizGeneQuery <- function(chr, pos1, pos2, db_path) {
    # Connect to SQLite database
    con <- dbConnect(RSQLite::SQLite(), dbname = db_path)

    # Build and execute the SQL query considering strand
    query <- sprintf(
        "SELECT * FROM \"SuREViz.Gene\" WHERE
    CHROM = '%s' AND
    (
        (start BETWEEN %d AND %d OR end BETWEEN %d AND %d)
    )",
        chr, pos1, pos2, pos1, pos2
    )

    result <- dbGetQuery(con, query)

    # Close the database connection
    dbDisconnect(con)

    # Return the result
    return(result)
}

