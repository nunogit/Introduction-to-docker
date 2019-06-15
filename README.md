# Introduction-to-docker

## What is Docker

Docker is “a computer program that performs operating-system-level virtualization, also known as ‘containerization’” Wikipedia

**What are the advantages**

Docker is designed to enclose environments inside an image / a container. What this allows, for example, is to have a Linux machine on a Macbook, or a machine with R 3.3 when your main computer has R 3.5. This means, for example,  that you can use older versions of a package for a specific task, while still keeping the package on your machine up-to-date.
This allows easy reproducibility, simpler deployments and application modularity.

## Anatomy of a Docker container script

### Creating an enviroment

**Dockerfile**

```
FROM r-base

COPY ./helloWorld.R /root/
ENTRYPOINT ["Rscript",  "/root/helloWorld.R"]
```

**Dockerfile explained**

1. Choose an image from [dockerhub](https://hub.docker.com/). Dockerhub contains hundreds of prebuilt docker containers. One just needs to pick the one closer to the end result. There are images with the basic OSs (Ubuntu, Debian, etc.); Development environments (R-studio; R with Bioconductor); With frameworks (Java1.8; Tomcat8); Pick the one closer to your end result.
> FROM r-base

2. Copy the file you wish to run. First argument local file system. Second argument Docker image file system
> COPY ./helloWorld.R /root/

3. Run the code
> ENTRYPOINT ["Rscript",  "/root/helloWorld.R"]

**Build the docker image**

``
``

**Run the docker image**

``
``

## Python Examples

### 1. [Python Hello World](examples/)
### 2. [Python Custom Script]

## R Examples

### 1. [R Hello World]
### 2. [R Run script]
### 3. [R Run an interactive shell on R Studio]

