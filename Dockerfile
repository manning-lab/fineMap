FROM r-base:3.5.0

MAINTAINER tmajaria@broadinstitute.org

RUN apt-get update && apt-get -y install build-essential libcurl4-openssl-dev git

RUN git clone https://github.com/gkichaev/PAINTOR_V3.0.git

RUN cd /PAINTOR_V3.0 && \
	bash install.sh && \
	mv ./PAINTOR /bin && \
	cd / && \
	rm -rf PAINTOR_V3.0

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN echo "install.packages(c('data.table','dplyr','tidyr'), repos='http://cran.us.r-project.org')" > install.R && \
	echo "source('https://bioconductor.org/biocLite.R')" >> install.R && \
	echo "biocLite(c('SNPRelate','SeqArray','GenomicRanges','SeqVarTools'))" >> install.R && \
	R --vanilla < install.R && \
	rm install.R

RUN cd && \
	git clone https://github.com/manning-lab/fineMap.git