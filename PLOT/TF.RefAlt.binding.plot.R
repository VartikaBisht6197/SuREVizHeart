# /HDD/data2/scratch/vartika/projects/VB230509_combine_WGS_SuRE/R/TF.RefAlt.binding.plot.R

# Load necessary libraries
library(patchwork)
library(memes)
library(BSgenome)
library(universalmotif)
library(ggplot2)
library(ggseqlogo)
library(data.table)
library(GenomicRanges)

# Define a function named Tf.alt.ref.overview.plot that takes two arguments: 'd' and 'JASPAR'
# d is a dataframe with 1 row and columns : ref.str , alt.str , refs.score , alts.score
# location of the JASPAR.meme database file
# https://meme-suite.org/meme/db/motifs
# Make sure that the names of the JASAPAR motif is the index name as well
Tf.alt.ref.overview.plot = function(d, JASPAR) {
  
  if(d$strand == "2"){matrix_data = motif_rc(JASPAR[[d$motif_id]])@motif}else{matrix_data = JASPAR[[d$motif_id]]@motif}
  
  ref_char_add = nchar(d$REF) - 1
  alt_char_add = nchar(d$ALT) - 1
  
  if(as.numeric(d$pos)+alt_char_add > ncol(matrix_data)){
    alt_char_add = ncol(matrix_data)-as.numeric(d$pos)
  }
  if(as.numeric(d$pos)+ref_char_add > ncol(matrix_data)){
    ref_char_add = ncol(matrix_data)-as.numeric(d$pos)
  }
  
  # Extract relevant data from the input dataframe 'd'
  ref.str = d$ref.str
  alt.str = d$alt.str
  refs.score = round(as.numeric(d$refs.score),0)
  alts.score = round(as.numeric(d$alts.score),0)
  refs.prob = d$refs.pval
  alts.prob = d$alts.pval
  
  # Create a vector of intercepts for probability lines
  intercepts = c(refs.score, alts.score)
  
  # Generate the title for the plot
  title.text = paste0(
    "TF : ",
    JASPAR[[d$motif_id]]@altname, "\n",
    "REF pval = ",
    format(refs.prob, big.mark = ","),
    " , ALT pval = ",
    format(alts.prob, big.mark = ","), "\n",
    "REF score = ",
    round(as.numeric(refs.score), digits = 0),
    " , ALT score = ",
    round(as.numeric(alts.score), digits = 0)
  )
  
  prob.at.variant.loc = paste(paste0("Probability at position ",
                                     seq(as.numeric(d$pos),(as.numeric(d$pos)+max(ref_char_add,alt_char_add)),1),
                                     " : ",
                                     (lapply(lapply(seq(as.numeric(d$pos),(as.numeric(d$pos)+max(ref_char_add,alt_char_add)),1),function(x) matrix_data[,x]), 
                                             function(y) paste(paste0(rownames(matrix_data), " : ", round(y, digits = 4)), collapse = " , ")))),collapse = "\n")
  
  # Create the logo plot for the reference sequence (g1)
  g1 = ggplot() + geom_logo(create_motif(d$ref.str)@motif) + theme_logo() + 
    geom_rect(aes(xmin = as.numeric(d$pos) - 0.5, xmax = as.numeric(d$pos) + 0.5 + ref_char_add,
                  ymin = 3, ymax = -1), fill = NA, color = "black", size = 0.5) + ylab("Reference")
  
  # Create the logo plot for the alternate sequence (g2)
  g2 = ggplot() + geom_logo(create_motif(d$alt.str)@motif) + theme_logo() + 
    geom_rect(aes(xmin = as.numeric(d$pos) - 0.5, xmax = (as.numeric(d$pos) + 0.5 + alt_char_add),
                  ymin = 3, ymax = -1), fill = NA, color = "black", size = 0.5) + ylab("Alternate")
  
  # Create the logo plot for matrix data (g3)
  g3 = ggplot() + geom_logo(matrix_data, method = "bits") + theme_logo() + 
    geom_rect(aes(xmin = as.numeric(d$pos) - 0.5, xmax = as.numeric(d$pos) + 0.5 + max(ref_char_add,alt_char_add),
                  ymin = 3, ymax = -1), fill = NA, color = "black", size = 0.5) + 
    ylab(JASPAR[[d$motif_id]]@altname) + 
    xlab(prob.at.variant.loc)
  
  # Combine the plots using patchwork
  g = g1 / g2 / g3
  g = guide_area() + g + 
    plot_layout(guides = "collect", 
                nrow = 2, heights = c(1, 20)) + 
    plot_annotation(title = title.text) & 
    theme(plot.title = element_text(hjust = 0.5))
  
  # Return the combined plot
  return(g)
}
