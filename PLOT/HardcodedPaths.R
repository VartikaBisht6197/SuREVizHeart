pvalcutoff <- 0.004581764808767
SureSignalcutoff <- 3

DBQueryScriptsDIR <- file.path(appDIR, "PLOT/QUERY")
DBDIR <- file.path(appDIR, "DATABASE")
plotDIR <- file.path(appDIR, "PLOT")

# bigwigs <- list(
#   "SuRE\nX57" = file.path(DBDIR, "data/HLHS.X57.bw"),
#   "SuRE\nX59" = file.path(DBDIR, "data/HLHS.X59.bw"),
#   "SuRE\nX67" = file.path(DBDIR, "data/HLHS.X67.bw"),
#   "SuRE\nX68" = file.path(DBDIR, "data/HLHS.X68.bw"),
#   "SuRE\nX86" = file.path(DBDIR, "data/HLHS.X86.bw")
# )

bigwigs <- list(
  "Patient1" = file.path(DBDIR, "data/HLHS.X57.bw"),
  "Patient2" = file.path(DBDIR, "data/HLHS.X59.bw"),
  "Patient3" = file.path(DBDIR, "data/HLHS.X67.bw"),
  "Patient4" = file.path(DBDIR, "data/HLHS.X68.bw"),
  "Patient5" = file.path(DBDIR, "data/HLHS.X86.bw")
)

AC16ATACbw <- list(
  "AC16\nATAC" = file.path(DBDIR, "data/AC16.ATACseq.hg38.bw")
)