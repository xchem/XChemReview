FROM rocker/shiny-verse:latest

# System Libs
RUN apt-get update && apt-get install -y \
	sudo \
	pandoc \
	pandoc-citeproc \
	libcurl4-gnutls-dev \
	libcairo2-dev \
	libxt-dev \
	libssl-dev \
	libssh2-1-dev \
	libpq-dev

# Install R packages
RUN R -e "install.packages(c('devtools', 'caTools','shiny','htmlwidgets', 'DT', 'remotes'), repos='http://cran.rstudio.com/')"
#RUN R -e "library(remotes); remotes::install_github('tjgorrie/nglShiny')"
RUN R -e "library(remotes); remotes::install_github('r-dbi/RPostgres')" 
RUN R -e "install.packages('RPostgres', repos='http://cran.rstudio.com/')"

# Copy App to Image
# COPY app.R /srv/shiny-server
COPY Pages /srv/shiny-server/Pages

# Port
EXPOSE 3838

# Permissions
RUN sudo chown -R shiny:shiny /srv/shiny-server

# Run App?
CMD ["Rscript", "/dls/science/users/mly94721/xchemreview/staging/app.R", "0.0.0.0", "3838"]
