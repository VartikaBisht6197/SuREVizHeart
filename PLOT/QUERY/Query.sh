#!/bin/bash

################################################################################
# SuRE-Viz Query Script
# ---------------------
# This script queries different databases using SuRE-Viz query scripts for
# specific genomic coordinates and variants. It accepts command line options
# to specify the chromosome, position, reference allele, alternate allele,
# flank size, output directory, Query directory, and Database directory.
# The script utilizes the getopts function for option parsing and runs SuRE-Viz
# query scripts in parallel for clinical, main, gene, and JASPAR TF databases.
#
# Usage:
# bash Query.sh -c <chr> -p <pos> -r <ref> -a <alt> -f <flank> -o <outdir>
#               -q <DBQueryScriptDIR> -d <dbdir>
#
# Parameters:
#   -c : Chromosome
#   -p : Genomic position
#   -r : Reference allele
#   -a : Alternate allele
#   -f : Flank size
#   -o : Output directory
#   -q : Query directory
#   -d : Database directory
#
# Example:
# bash Query.sh -c chr19 -p 56252947 -r A -a C -f 100000 -o /path/to/output_directory
#               -q /path/to/query_directory -d /path/to/database_directory
################################################################################

# Default values
chr=chr19
pos=56252947
ref=A
alt=C
flank=100000

# Parse command line options
while getopts ":c:p:r:a:f:o:q:d:" opt; do
  case $opt in
    c)
      chr="$OPTARG"
      ;;
    p)
      pos="$OPTARG"
      ;;
    r)
      ref="$OPTARG"
      ;;
    a)
      alt="$OPTARG"
      ;;
    f)
      flank="$OPTARG"
      ;;
    o)
      outdir="$OPTARG"
      ;;
    q)
      DBQueryScriptDIR="$OPTARG"
      ;;
    d)
      dbDIR="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if required options are provided
if [ -z "$chr" ] || [ -z "$pos" ] || [ -z "$ref" ] || [ -z "$alt" ] || [ -z "$flank" ] || [ -z "$outdir" ] || [ -z "$DBQueryScriptDIR" ] || [ -z "$dbDIR" ]; then
  echo "Usage: $0 -c <chr> -p <pos> -r <ref> -a <alt> -f <flank> -o <outdir> -q <DBQueryScriptDIR> -d <dbdir>"
  exit 1
fi

# Specify query scripts
QCLIN=${DBQueryScriptDIR}/DBquery.clin.r
QMAIN=${DBQueryScriptDIR}/DBquery.main.r
QGENE=${DBQueryScriptDIR}/DBquery.gene.r
QJASPAR=${DBQueryScriptDIR}/DBquery.jasparTF.r

# Specify database paths
DCLIN=${dbDIR}/SuREViz.ClinVar.db
DMAIN=${dbDIR}/SuREViz.Main.db
DGENE=${dbDIR}/hg38.refGene.db
DJASPAR=${dbDIR}/JASPAR.TF.raQTL.db

# Construct the output directory path
QueryOUTDIR=${outdir}/${chr}_${pos}_${ref}_${alt}

# Create the Query directory if it doesn't exist
if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
  echo "Making a new directory: $outdir"
fi

# Create the specific output directory for the current query
mkdir "$QueryOUTDIR"


# Run query scripts with input variables
(
    Rscript $QCLIN --chr $chr --pos $pos --db_path $DCLIN --outfile $QueryOUTDIR/CLIN.txt &
    Rscript $QMAIN --chr $chr --pos $pos --flank $flank --db_path $DMAIN --outfile $QueryOUTDIR/MAIN.txt &
    Rscript $QGENE --chr $chr --pos $pos --flank $flank --db_path $DGENE --outfile $QueryOUTDIR/GENE.txt &
    Rscript $QJASPAR --chr $chr --pos $pos --ref $ref --alt $alt --db_path $DJASPAR --outfile $QueryOUTDIR/JASPAR.txt
)

# Wait for all background processes to finish
wait
