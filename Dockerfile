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
	libpq-dev \ 
        libmariadb-dev

# Install R packages
RUN R -e "install.packages(c('devtools', 'caTools','shiny','shinydashboard','shinyjqui','shinyWidgets','ggplot2','plotly','htmlwidgets', 'DT', 'remotes', 'httr', 'sendmailR', 'shinyjs','bio3d'), repos='http://cran.rstudio.com/')"
RUN R -e "library(remotes); remotes::install_github('r-dbi/RPostgres'); remotes::install_github('r-dbi/RMariaDB')" 
RUN R -e "install.packages(c('RPostgres','RMariaDB', 'waiter'), repos='http://cran.rstudio.com/')"
# Copy App to Image
# COPY app.R /srv/shiny-server
COPY Pages /srv/shiny-server/Pages

# Port
EXPOSE 3838

# Permissions
RUN sudo chown -R shiny:shiny /srv/shiny-server

# Run App?
#CMD ["Rscript", "/dls/science/groups/i04-1/software/xchemreview/app.R", "0.0.0.0", "3838"]
