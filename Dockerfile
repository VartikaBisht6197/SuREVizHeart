# Start with Debian stable-slim for better flexibility, control, and smaller image size.
# Install micromamba manually to ensure a minimal and customizable environment.
FROM debian:stable-slim

# https://micromamba-docker.readthedocs.io/en/latest/advanced_usage.html#adding-micromamba-to-an-existing-docker-image
# Install dependencies for micromamba
RUN apt-get update && apt-get install -y \
    wget \
    ca-certificates \
    gnupg \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba
RUN cd /usr/local && wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Set environment variables for micromamba
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

# Set working directory
WORKDIR /shiny-server

# Add conda-forge and bioconda channels and set strict priority
RUN micromamba config --add channels conda-forge \
    && micromamba config --add channels bioconda \
    && micromamba config --set channel_priority strict

# Refresh channel metadata and verify channel configuration
RUN micromamba update --all -y -c conda-forge -c bioconda && micromamba clean -a -y

# Test micromamba by creating an environment with a simple package
RUN micromamba create -y -n testenv -c conda-forge zlib

# Now attempt to add python or r-base to the test environment
RUN micromamba install -y -n testenv -c conda-forge python=3.8.12 r-base || \
    echo "Manual fallback for Python/R installation if needed"

# Use micromamba to install the required R packages
# Explicitly specify the conda-forge channel where r-base is available
# RUN micromamba create -y -n sureviz -c conda-forge -c bioconda \
#     r-base r-data.table r-dbi r-dt r-formattable r-ggplot2 r-ggseqlogo \
#     bioconductor-genomicranges r-gridextra r-here r-kableextra \
#     bioconductor-universalmotif r-optparse r-patchwork r-pheatmap \
#     r-plotly r-rcolorbrewer r-rsqlite r-shiny r-shinyfiles \
#     r-shinywidgets r-shinyalert r-shinycssloaders r-shinydashboard \
#     r-shinyjs r-tidyr r-zip \
#     bioconductor-bsgenome.hsapiens.ucsc.hg38 bioconductor-rtracklayer \
#     && micromamba clean -a -y

# Activate the environment by default
# RUN echo "micromamba activate sureviz" >> ~/.bashrc

# Copy the Shiny app files
COPY . /shiny-server/

# Expose the Shiny app port
EXPOSE 3838

# Run the Shiny app
CMD R -e "shiny::runApp('/shiny-server/SuREViz', host = '0.0.0.0', port = 3838)"
