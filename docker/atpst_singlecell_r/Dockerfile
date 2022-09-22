FROM rocker/rstudio:4.0.3

RUN apt-get update \
	&& apt-get -y install zlib1g-dev \
	&& apt-get -y install libbz2-dev \
	&& apt-get -y install liblzma-dev \
	&& apt-get -y install libcurl4-openssl-dev \
	&& apt-get -y install libxml2 \	
	&& apt-get -y install libglpk40 \
	&& apt-get -y install libxt-dev \
	&& apt-get -y install libboost-dev \
	&& apt-get -y install libhdf5-dev \
	&& apt-get -y install cmake \
	&& apt-get -y install libboost-dev \
	&& apt-get -y install libboost-iostreams-dev \
	&& apt-get -y install libboost-log-dev \
	&& apt-get -y install libboost-system-dev \
	&& apt-get -y install libboost-test-dev \
	&& apt-get -y install libproj-dev \
	&& apt-get -y install build-essential \
	&& apt-get -y install libcairo2-dev \
	&& apt-get -y install libxt-dev \
	&& apt-get -y install libgdal-dev \
	&& apt-get -y install libudunits2-dev

ADD install_bioconductor_packages.R /
ADD install_R_packages.R /
ADD install_seurat.R /
ADD install_monocle.R /

RUN R -e "source(\"install_R_packages.R\")"
RUN R -e "source(\"install_seurat.R\")"
RUN R -e "source(\"install_bioconductor_packages.R\")"
RUN R -e "source(\"install_monocle.R\")"