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
RUN R -e "
install.packages('caTools', repos='http://cran.rstudio.com/')
install.packages('devtools', repos='http://cran.rstudio.com/')
install.packages('shiny', repos='http://cran.rstudio.com/')
install.packages('htmlwidgets', repos='http://cran.rstudio.com/')
install.packages('DT', repos='http://cran.rstudio.com/')
install.packages('remotes', repos='http://cran.rstudio.com/')
library(remotes)
remotes::install_github('r-dbi/RPostgres') 
install.packages('RPostgres', repos='http://cran.rstudio.com/')
"

# Copy App to Image
COPY app.R /srv/shiny-server
# COPY getDatafd.py /srv/shiny-server
# COPY Data /srv/shiny-server/Data
COPY Pages /srv/shiny-server/Pages
COPY nglShiny /srv/shiny-server/nglShiny
# COPY Responses /srv/shiny-server/Responses

# Port
EXPOSE 3838

# Permissions
RUN sudo chown -R shiny:shiny /srv/shiny-server

# Run App?
CMD ["R", "-e", "source('/srv/shiny-server/app.R')"]
