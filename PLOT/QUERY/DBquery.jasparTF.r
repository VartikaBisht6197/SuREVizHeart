################################################################################
# SuREViz JASPAR TF Query Function
# --------------------------------
# This R script defines a function, suREVizJASPARTFQuery, to query the SuREViz.JASPAR.TF
# SQLite database for transcription factor binding sites using specific genomic
# coordinates and variants. It makes use of the 'DBI' and 'RSQLite' libraries for
# database connectivity.
#
# Function Usage:
# result <- suREVizJASPARTFQuery(chr, pos, ref, alt, db_path)
#
# Example usage:
# result <- suREVizJASPARTFQuery("chr19", 56252947, "A", "C", "/path/to/database.sqlite")
#
# Parameters:
#   chr      : Chromosome
#   pos      : Genomic position
#   ref      : Reference allele
#   alt      : Alternate allele
#   db_path  : Path to SQLite database
#
# Example:
# result <- suREVizJASPARTFQuery("chr19", 56252947, "A", "C", "/path/to/database.sqlite")
################################################################################

# Libraries
library(DBI)
library(RSQLite)

# Main function to query the SuREViz.JASPAR.TF SQLite database
suREVizJASPARTFQuery <- function(chr, pos, ref, alt, db_path) {
    # Connect to SQLite database
    con <- dbConnect(RSQLite::SQLite(), dbname = db_path)

    # Build and execute the SQL query
    query <- sprintf(
        "SELECT * FROM \"SuREViz.JASPAR.TF\" WHERE
    CHROM = '%s' AND POS = %d AND REF = '%s' AND ALT = '%s'",
        chr, pos, ref, alt
    )

    result <- dbGetQuery(con, query)

    # Close the database connection
    dbDisconnect(con)

    # Return the result
    return(result)
}

