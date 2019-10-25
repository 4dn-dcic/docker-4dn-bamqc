FROM ubuntu:16.04
RUN apt-get update -y && apt-get install -y \
    bc \
    bzip2 \
    curl \
    g++ \
    git \
    less \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses-dev \
    libssl-dev \
    libtbb-dev \
    littler \
    make \
    man \
    openjdk-8-jre \
    parallel \
    perl \
    r-base-dev \
    time \
    unzip \
    vim \
    wget \
    zlib1g-dev \
    libmysqlclient-dev

RUN apt-get update -y && apt-get install -y \
   python3.5-dev \
   python3-setuptools \
   && wget https://bootstrap.pypa.io/get-pip.py \
   && python3.5 get-pip.py

# installing python libraries
RUN pip install pysam==0.15.3
RUN pip install click==7.0

# download tools
WORKDIR /usr/local/bin
COPY downloads.sh .
RUN . downloads.sh

# set path
ENV PATH=/usr/local/bin/scripts:$PATH

# supporting UTF-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# wrapper
COPY scripts/ .
RUN chmod +x run*.sh

# default command
CMD ["ls","/usr/local/bin/"]
