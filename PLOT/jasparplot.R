library(BSgenome.Hsapiens.UCSC.hg38)
library(universalmotif)
library(ggseqlogo)
library(gridExtra)
library(patchwork)
library(ggplot2)
library(memes)
library(tidyr)

plots_in_one_set = 3
plots_cols = 2

jaspar.analysis = function(jaspar.data, JASPAR.PPM){
  
  TFs <- unlist(strsplit(jaspar.data$motif_id,","))
  seq <- DNAStringSet(getSeq(BSgenome.Hsapiens.UCSC.hg38, jaspar.data$CHROM, (as.numeric(jaspar.data$POS) - 15), (as.numeric(jaspar.data$POS) + nchar(jaspar.data$REF) + 15 - 1) ))
  names(seq) <- jaspar.data$POS

  # https://rockefelleruniversity.github.io/ebf1_motif/fimo_out_2/fimo.html
  # FIMO.background = c("A"= 0.219, "C" = 0.281, "G" = 0.281, "T" = 0.219)
  JASPAR.PPM = JASPAR.PPM[which(unlist(lapply(JASPAR.PPM,function(x) x@name %in% TFs)))]
  names(JASPAR.PPM) = unlist(lapply(JASPAR.PPM, function(x) x@name))
  
  FIMO = data.frame(separate_rows(jaspar.data, start , end, strand, pos, motif_id, motif_alt_id, refs.score, alts.score, 
                                  refs.pval, alts.pval, absdiff, max.score, effect.JASPAR, tax_group, tf_family, tf_class,
                                  data_type,uniprot_ids,pubmed_ids,Gene_Name, sep = ","))
  FIMO$refs.score = round(as.numeric(FIMO$refs.score),0)
  FIMO$alts.score = round(as.numeric(FIMO$alts.score),0)
  FIMO$effect.JASPAR <- round(as.numeric(FIMO$effect.JASPAR),3)
  # Only keep the highest impacting alignment for a TF
  FIMO = rbindlist(lapply(split(FIMO,FIMO$motif_alt_id),function(x) x[which(abs(x$effect.JASPAR) == max(abs(x$effect.JASPAR))),] ))
  
  snp <- paste(jaspar.data$CHROM,jaspar.data$POS,jaspar.data$REF,jaspar.data$ALT,sep = "-")
  
  # Convert the sequence to a PPM (Position Probability Matrix)
  ref_str = as.character(seq[[1]])
  alt_str = paste0(substr(seq[[1]],1,15),jaspar.data$ALT,substr(seq[[1]],17+nchar(jaspar.data$REF)-1,31+nchar(jaspar.data$REF)-1))
  
  seqlen = nchar(seq) - abs(nchar(alt_str) - nchar(ref_str))
  
  if(nchar(ref_str) > nchar(alt_str)){
    ref_str = substr(ref_str,1,nchar(alt_str))
  }
  
  if(nchar(alt_str) > nchar(ref_str)){
    alt_str = substr(alt_str,1,nchar(ref_str))
  }
  
  ref_ppm <- create_motif(ref_str)@motif
  alt_ppm <- create_motif(alt_str)@motif  
  
  ref_plot <- ggplot() +
    geom_rect(aes(xmin = 16 - 0.5, xmax = 16 + 0.5 + nchar(jaspar.data$REF) - 1,
                  ymin = 3, ymax = -1), fill = "black", size = 0.5,alpha = 0.09) +
    geom_logo(as.matrix(ref_ppm), method = "bits") +
    theme_logo() + labs(y = "REF") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  alt_plot <- ggplot() +
    geom_rect(aes(xmin = 16 - 0.5, xmax = 16 + 0.5 + nchar(jaspar.data$ALT) - 1,
                  ymin = 3, ymax = -1), fill = "black", size = 0.5,alpha = 0.09) +
    geom_logo(as.matrix(alt_ppm), method = "bits") +
    theme_logo() + labs(y = "ALT") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  
  # Create individual plots for each TF motif
  wrplt <- list()
  
  for (i in 1:nrow(FIMO)) {
    
    highlight = NULL
    p_motif = NULL
    
    motif_id <- FIMO$motif_id[i]
    motif_alt_id <- FIMO$motif_alt_id[i]
    start_pos <- as.numeric(FIMO$start[i])
    end_pos <- as.numeric(FIMO$end[i])
    matched_seq <- FIMO$matched_sequence[i]
    motif_strand <- FIMO$strand[i]
    effect_JASPAR <- FIMO$effect.JASPAR[i]
    
    # Get the motif
    motif <- JASPAR.PPM[[motif_id]]@motif
    if(motif_strand == "-"){motif = motif_rc(JASPAR.PPM[[motif_id]])@motif}else{motif = JASPAR.PPM[[motif_id]]@motif}
    
    # Create an empty PPM matrix of the same length as the sequence
    seq_length <- 31
    empty_matrix <- matrix(0, nrow = 4, ncol = seq_length)
    rownames(empty_matrix) <- c("A", "C", "G", "T")
    
    # Insert the motif into the correct position in the empty matrix
    motif_length <- ncol(motif)
    empty_matrix[, start_pos:(start_pos + motif_length - 1)] <- motif
    
    # highlight
    if(effect_JASPAR>0){
      motif_title = paste(snp,"distrupts TF",motif_alt_id,"(",motif_id,")","binding to REF",sep=" ")
    } else { 
      motif_title = paste(snp,"enhnaces TF",motif_alt_id,"(",motif_id,")","binding to ALT",sep=" ") }
    
    # Create the motif plot
    p_motif <- ggplot() +
      geom_logo(empty_matrix, method = "bits") +
      theme_logo() +
      labs(y = motif_alt_id) +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
      ) 
    
    plts = list(ref_plot+labs(title = motif_title), alt_plot, p_motif)
    
    # Add the plot to the list
    # plots[[length(plots) + 1]] <- p_motif
    wrplt[[i]] <- wrap_plots(plts, ncol = 1)
  }
  
  combined_plot <- wrap_plots(wrplt, ncol = 1)
  
  # Create the pretty table and add to wrapped plot list 
  FIMO <- FIMO[, c("motif_id", "motif_alt_id",  "tf_family", "tf_class", "effect.JASPAR" )]
  colnames(FIMO) = c("JASPAR 2022 Motif ID" , "JASPAR 2022 Motif Name" , "TF Family" ,"TF Class" , "REF/ALT TF binding")

  # Combine the table and the plots 
  return(list(plot = combined_plot, table = FIMO))

}