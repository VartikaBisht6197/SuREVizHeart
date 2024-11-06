# Use the Rocker Shiny base image
FROM rocker/shiny-verse:latest

# Copy the app files into the Docker image
COPY . /srv/shiny-server/

# Copy the environment file to the container
COPY env/SuREViz.yml /srv/shiny-server/env/SuREViz.yml

# Set the working directory
WORKDIR /srv/shiny-server

# Create mount point for Azure Blob Storage
RUN mkdir -p /mnt/azureblob

# Expose the Shiny app port
EXPOSE 3838

# Command to run the Shiny app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host = '0.0.0.0', port = 3838)"]
