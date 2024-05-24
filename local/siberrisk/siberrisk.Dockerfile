# Use official miniconda image as the base image
FROM continuumio/miniconda3 as miniconda3

# Set the working directory
WORKDIR /usr/local/src

# Install gcc compiler
RUN apt-get update && apt-get install build-essential -y

# Install python libraries
RUN conda install -y python=3.11
COPY ./environment.yml /usr/local/src/environment.yml
RUN conda env update --name base --file /usr/local/src/environment.yml

# Set volume
VOLUME [ "/usr/local/src" ]
