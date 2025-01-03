# SuREViz: A Web-Based SuRE MPRA Data Exploration Tool 🌟

Welcome to **SuREViz**! 🎉  
A cutting-edge platform designed for exploring the functional impact of genetic variants assessed by the SuRE (Survey of Regulatory Elements) massively parallel reporter assay. SuREViz empowers researchers and data enthusiasts with interactive tools to visualize and analyze the interplay between genomic variants and gene regulation.

---

<p align="center">
  <img src="images/Fig1.png" alt="Introductory Figure" width="700px">
</p>

*Figure 1: An overview of the SuREViz platform and its capabilities. Descriptions for each subfigure will be provided soon.*

*Note: The original figure is in PDF format (`images/Fig1.pdf`). For proper embedding here, it has been converted to an image format.*

---

## 🌟 Key Features

✨ **Interactive Visualization**: Explore over 4.7 million variants categorized into raQTLs and non-raQTLs.  
✨ **Genomic Data Integration**: Combines functional, genomic, and clinical datasets for enriched insights.  
✨ **Custom Data Upload**: Easily compare user-provided MPRA data, BigWig, and BED files.  
✨ **Flexible Querying**: Search by variants (`chr:pos`), gene names, or customized regions.  
✨ **Downloadable Results**: Export your analyzed data in user-friendly formats.  

---

## 🔬 Data Sources

- **Processed Data**: Available on OSF at [pyh83](https://osf.io/pyh83/).  
- **Raw Data**: Accessible through GEO at [GSE251776](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251776).  

---

## 🔬 What is SuRE?

SuRE (Survey of Regulatory Elements) is a high-throughput assay used to study gene regulation. It maps allele-specific expression back to the genome, allowing researchers to assess regulatory activity. SuREViz showcases data from heart-related studies, including six complex plasmid libraries representing 4.7 million variants.  

### 🧬 What are raQTLs?  
**Regulatory Allelic Quantitative Trait Loci (raQTLs)** are variants showing significant differences in expression between reference and alternate alleles, highlighting their regulatory impacts.

---

## 🚀 Getting Started

To begin your journey with SuREViz:  
1. Ensure you have a modern web browser like **Google Chrome**, **Mozilla Firefox**, or **Safari**.  
2. Access the platform here: **[SuREViz](http://192.168.107.99:6197)**  

---

## 🔍 App Overview  

### 🛠️ Search Functionalities  
- **Variant Query**: Enter in `chr:pos` format (e.g., `chr12:128797635`).  
- **Gene Query**: Search by gene name (case-insensitive, e.g., `TBP` or `tbP`).  
- **Custom Data**: Upload your MPRA data, BigWig, or BED files for seamless integration and visualization.

### 📊 Visualization Modes  
1. **Variant View**: Explore detailed allele expression and statistical significance.  
2. **Region View**: Delve into genomic features and surrounding variants.  

### 🗂️ Key Tabs  
- **Functional Impact Assessment**: Insights into SuRE expression levels, gene plots, and statistical analysis.  
- **Variant Data Overview**: Comprehensive tables with allele frequencies, impact scores, and more.  
- **Gene Expression Overview**: Transcriptomic data from heart and tissue-specific studies.  
- **SuRE Profiles**: Detailed views of SuRE data, ATAC-seq profiles, and evolutionary conservation scores.  

---

## 📂 Downloadable Results

Export your findings in these formats:  
1. **JASPAR2022_info.csv**: Transcription factor impact details.  
2. **SNP_info.csv**: Tabular data on SNPs within the query region.  
3. **SuREX_.bedGraph**: SuRE profiles for queried loci.  

---

## 🧠 Technical Details  

SuREViz integrates:  
- SuRE MPRA data.  
- Transcription factor binding predictions.  
- Conservation scores across mammalian species.  
- ATAC-seq data from AC16 cardiomyocyte cell lines.  

---

## 🤝 Contributing and Support  

We ❤️ contributions!  
- Found a bug? Open an **issue**.  
- Have an idea? Submit a **pull request**.  
- Questions? Contact **Vartika Bisht** for direct support.

---

### 🌟 Thank You for Exploring SuREViz!  
Together, let’s uncover the regulatory mechanisms behind genetic variants. 🌍✨
