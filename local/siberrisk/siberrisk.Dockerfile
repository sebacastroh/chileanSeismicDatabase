# Use official anaconda image as the base image
FROM continuumio/anaconda3 as anaconda

# Set the working directory
WORKDIR /usr/local/src

# Install python libraries
COPY ./requirements.txt /usr/local/src/requirements.txt
RUN pip3 install -r /usr/local/src/requirements.txt

# Set volume
VOLUME [ "/usr/local/src" ]

# Install gcc compiler
RUN apt-get update && apt-get install build-essential -y
