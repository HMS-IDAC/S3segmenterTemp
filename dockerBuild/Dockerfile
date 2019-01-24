FROM ubuntu

# Install the MCR dependencies and some things we'll need and download the MCR from Mathworks -silently install it

RUN apt-get update &&\
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    openjdk-8-jdk-headless \
    libasound2-dev \
    nano \
    unzip \
    xorg \
    wget \
    libxt6 \
    python3-dev \
    python3-pip && \
    mkdir /mcr-install && \
    mkdir /opt/mcr && \
    cd /mcr-install && \
    wget http://www.mathworks.com/supportfiles/downloads/R2018a/deployment_files/R2018a/installers/glnxa64/MCR_R2018a_glnxa64_installer.zip && \
    cd /mcr-install && \
    unzip -q MCR_R2018a_glnxa64_installer.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install &&\
    rm -rf /var/lib/apt/lists/*

# Configure environment variables for MCR
ENV LD_LIBRARY_PATH /opt/mcr/v94/runtime/glnxa64:/opt/mcr/v94/bin/glnxa64:/opt/mcr/v94/sys/os/glnxa64
ENV XAPPLRESDIR /opt/mcr/v94/X11/app-default


RUN pip3 install -q -U \
    pip \
    numpy \
    PyYAML \
    pandas 

# Build input, output, config, and segmentation directories
# Add necessary matlab and python code to /S3segmenter
VOLUME /input
VOLUME /output
VOLUME /config
WORKDIR /S3segmenter

#Add compiled file from git LFS
ADD S3segmenterWrapperDocker /S3segmenter/S3segmenterWrapperDocker
RUN chmod +x S3segmenterWrapperDocker

#Add rest of segmentation code 
ADD matlabDependencies /S3segmenter/matlabDependencies
ADD S3segment.py /S3segmenter/S3segment.py
ADD run_S3segmenterWrapperDocker.sh /S3segmenter/run_S3segmenterWrapperDocker.sh

# Run segmentation when the container launches
CMD ["python3", "S3segment.py"]


