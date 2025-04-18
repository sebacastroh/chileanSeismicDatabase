# Use official miniconda image as the base image
FROM condaforge/miniforge3:24.11.3-2 AS miniforge3

# Set the working directory
WORKDIR /usr/local/src

# Install gcc compiler
RUN apt-get update && apt-get install build-essential -y

# Install python libraries
COPY ./environment.yml /usr/local/src/environment.yml
RUN conda env update --name base --file /usr/local/src/environment.yml

# Set volume
VOLUME [ "/usr/local/src" ]
