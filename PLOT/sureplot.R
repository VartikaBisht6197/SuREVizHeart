################################################################################

# SuRE Signal Plotting Script
# ---------------------------
# This R script generates a SuRE (Survey of Regulatory Elements) signal plot for a
# specified genomic variant. It uses the 'ggplot2' and 'plotly' libraries for
# visualization and assumes the existence of a 'MAIN.txt' dataset with relevant
# genomic information.
#
# Usage:
# source("sureplot.R"); snps_plot("chr", pos, "ref", "alt", suredata,pvalcutoff, SureSignalcutoff)
#
# Parameters:
#   chr              : Chromosome
#   pos              : Genomic position
#   ref              : Reference allele
#   alt              : Alternate allele
#   suredata         : SuRE dataset (assumed to be loaded externally)
#   pvalcutoff       : p-value cutoff for significance
#   SureSignalcutoff : SuRE Signal cutoff for significance
#
# Example:
# source("snps_plot_function.R"); snps_plot("chr19", 56252947, "A", "C", 100000,suredata, 0.004581764808767, 3)

################################################################################


# Load required libraries
library(ggplot2)
library(plotly)

# Function definition
snps_plot <- function(chr, pos1, pos2, pos, suredata, pvalcutoff, SureSignalcutoff , query_snps) {
  # Modify the suredata dataframe within the function
  suredata$raQTL = ifelse(max(suredata$ref.mean, suredata$alt.mean) > SureSignalcutoff & suredata$wilcox.p.value <= pvalcutoff, TRUE, FALSE)
  suredata$max.allele = ifelse(suredata$ref.mean > suredata$alt.mean, "Loss of function\n( REF > ALT )", "Gain of function\n( ALT > REF )")
  suredata$gnomad312_AF = round(as.numeric(suredata$gnomad312_AF) * 100, 2)
  suredata$sig = ifelse(suredata$wilcox.p.value < pvalcutoff, -log10(suredata$wilcox.p.value) * 100 , -log10(suredata$wilcox.p.value) / 100)
  suredata$min.allele.exp = signif(mapply(min, suredata$ref.mean, suredata$alt.mean),3)
  suredata$max.allele.exp = signif(mapply(max, suredata$ref.mean, suredata$alt.mean),3)
  

  # Create ggplot object
  p_signal <- ggplot(suredata, aes(x = POS,
        y = ref.mean,
        xend = POS,
        yend = alt.mean,
        alpha = sig,
        color = max.allele,
        size = -log10(wilcox.p.value),
        text = paste0("Click to centre on this variant","\n","\n",
                      "Postion: ", paste0(CHROM,"-",POS,"-",REF,"-",ALT), "\n",
                      "p-val: ", format(wilcox.p.value, scientific = TRUE, digits = 2) , "\n",
                      "Alternate allele frequency: ", gnomad312_AF, "%\n",
                      "Reference allele: ", REF, "\n",
                      "Alternate allele: ", ALT, "\n",
                      max.allele ,"\n",
                      "Min allele expression: ", min.allele.exp, "\n",
                      "Max allele expression: ", max.allele.exp ))) +
              geom_segment() + guides(size = "none", alpha = "none") + 
              scale_color_manual(values = c("#3c5379", "#64BEA5"), 
                          breaks = c("Gain of function\n( ALT > REF )","Loss of function\n( REF > ALT )"),
                          guide = guide_legend(title = "")) +
              theme_bw() +
              theme(
              legend.position = "right" ,
              axis.title.x = element_blank()) + 
              ylab("SuRE Signal") +  
              xlim(c(pos1, pos2)) +
                coord_cartesian(xlim = c(pos1, pos2))
  
  if(!is.na(pos) && nrow(query_snps) != 0){
  p_signal = p_signal +
  geom_segment(aes(x = pos, xend = pos, y = query_snps$ref.mean, yend = query_snps$alt.mean), color = "yellow", size = 0.5, alpha = 1) + 
  geom_vline(xintercept = pos, colour = "black", linetype = "dashed") }
  
  
  return(p_signal)
}
