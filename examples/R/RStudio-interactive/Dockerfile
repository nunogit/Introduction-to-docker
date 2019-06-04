FROM rocker/r-ver:3.4.4

ARG WHEN

RUN mkdir /home/analysis

RUN R -e "options(repos = \
  list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/${WHEN}')); \
  install.packages('tidystringdist')"

COPY script.R /home/analysis/script.R

CMD R -e "source('/home/analysis/script.R')"
