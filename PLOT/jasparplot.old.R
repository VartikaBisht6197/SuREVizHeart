# Load required libraries
library(universalmotif)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyr)

# Define functions
change.to.alt.seq = function(ref, alt, seq, flank) {
  if (substr(seq, flank + 1, flank + nchar(ref)) == ref) {
    return(paste0(substr(seq, 1, flank), alt, substr(seq, flank + 1 + nchar(ref), flank + nchar(seq))))
  } else {
    stop("Reference does not match sequence.")
  }
}

get.flanking.seq = function(chr, pos, ref, flank, genome) {
  loc <- GRanges(seqnames = chr, 
                 ranges = IRanges(start = pos - flank,
                                  end = pos + nchar(ref) - 1 + flank))
  
  sequence <- as.character(getSeq(genome, loc))
  return(sequence)
}

# Define JASPAR analysis function
jaspar.analysis = function(jaspar.data, JASPAR_PPM) {
  
  jaspar.data = data.frame(jaspar.data %>%
                             separate_rows(all_of(c("start", "end", "strand", "pos", "motif_id", "motif_alt_id",
                                                    "refs.score", "alts.score", "refs.pval", "alts.pval")), sep = ","))
  
  result_list = list()
  
  for (idx in 1:nrow(jaspar.data)) {
    d = jaspar.data[idx,]
    ref_sequence = get.flanking.seq(d$CHROM, d$POS, d$REF, 15, BSgenome.Hsapiens.UCSC.hg38)
    alt_sequence = change.to.alt.seq(d$REF, d$ALT, ref_sequence, 15)
    
    d$ref.str = substr(ref_sequence, d$start, d$end)
    d$alt.str = substr(alt_sequence, d$start, d$end)
    
    result_list[[idx]] = Tf.alt.ref.overview.plot(d, JASPAR_PPM)
  }
  
  names(result_list) = jaspar.data$motif_alt_id
  
  return(result_list)
}