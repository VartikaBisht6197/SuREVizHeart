# if pos is not NA make sure it exist , and if not throw and error
if(!is.na(pos)){
    if(nrow(query_snps)==0 ){
        source(file.path(DBQueryScriptsDIR, "DBquery.main.r"))
        suredata <- suREVizMainQuery(chr, pos1, pos2, file.path(DBDIR, "SuREViz.Main.db"))
        query_snps <- suredata[which(suredata$CHROM == chr & suredata$POS == pos), ]
        stop("The user is trying to search for a variant not in the SURE dataset")}
}