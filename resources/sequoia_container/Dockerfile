# Use a base R image (rocker/r-ver is stable and version-controlled)
FROM rocker/r-ver:4.3.0

# Install system dependencies required for some R packages (like vcfR)
RUN apt-get update && apt-get install -y \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libglpk-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install littler package to facilitate R package installations
RUN install2.r --error remotes

# Install required R packages
RUN install2.r --error \
    adegenet \
    vcfR

# Install specific version of sequoia package
RUN Rscript -e "remotes::install_version('sequoia', version = '3.2.0', upgrade = 'never', repos = 'https://cloud.r-project.org')"

# Verify installation by checking the package version
RUN Rscript -e "packageVersion('sequoia')"

# Build the Docker image
# docker build --platform linux/amd64 -t sequoia_3.2.0:latest .
# docker tag sequoia_3.2.0:latest shohei0212/sequoia_3.2.0:latest
# docker push shohei0212/sequoia_3.2.0:latest