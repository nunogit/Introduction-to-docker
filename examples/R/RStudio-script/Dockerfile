#Base image
FROM rocker/r-ver:3.4.4

#Argument to create Docker image
ARG WHEN

#Run command on the native OS at image build time / create directory
RUN mkdir /home/analysis

#Run command on the native OS at image buidl time / install extra cran and install package in R
RUN R -e "options(repos = \
  list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/${WHEN}')); \
  install.packages('tidystringdist')"

#Copy file form local file system to image file system
COPY script.R /home/analysis/script.R

#Run command at image boot time / run R script
CMD R -e "source('/home/analysis/script.R')"
