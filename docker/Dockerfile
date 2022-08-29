FROM rocker/rstudio:4.1.2

RUN apt-get update
RUN apt-get -y install zlib1g-dev 
RUN apt-get -y install libbz2-dev 
RUN apt-get -y install liblzma-dev 
RUN apt-get -y install libcurl4-openssl-dev 
RUN apt-get -y install libxml2 
RUN apt-get -y install libglpk40 
RUN apt-get -y install libxt-dev

ADD install_R_packages.R /

RUN R -e "source(\"install_R_packages.R\")"

