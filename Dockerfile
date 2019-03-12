#########
#
#     Dockerfile for Metagenomics Bioinformtics Taxonomy session
#
### To run the container for the first time with generic graphics:
# xhost +
# docker run -it -v /tmp/.X11-unix:/tmp/.X11-unix:rw --privileged -e DISPLAY=unix$DISPLAY \
# -v $HOME/:/home/training/ --device /dev/dri --privileged --name metagenomicstaxonomy \
# ebitraining/metagenomics:taxonomy
#
### To run with Nvidia graphics, add the following option:
# "-v /usr/lib/nvidia-340:/usr/lib/nvidia-340 -v /usr/lib32/nvidia-340:/usr/lib32/nvidia-340"
#
### To resume using an container:
# docker exec -it metagenomicstaxonomy /bin/bash
#
### To build the container:
# docker build -f ./Dockerfile -t metagenomicstaxonomy .
# docker tag metagenomicstaxonomy ebitraining/metagenomics:taxonomy
# docker push ebitraining/metagenomics:taxonomy
#
#########
FROM ubuntu:18.04
LABEL author="Mohamed Alibi" \
description="Docker image for Metagenomics Bioinformtics Taxonomy session." \
maintainer="Mohamed Alibi <alibi@ebi.ac.uk>"

## Pre requirements
########
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
ENV DEBIAN_FRONTEND noninteractive
ENV LC_ALL en_GB.UTF-8
ENV LANG en_GB.UTF-8
ENV LANGUAGE en_GB:en

RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get install -y --reinstall language-pack-en \
    && apt-get install -y gnupg \
    && locale-gen en_GB.UTF-8

RUN apt-get update; apt-get install -y build-essential ca-certificates libbz2-dev liblzma-dev gfortran \
    libncurses5-dev libncursesw5-dev zlib1g-dev automake pkg-config unzip openjdk-8-jre-headless curl python-pip \
    git wget sudo autoconf make xml2 locales libjpeg-dev zlibc libjpeg62 libxslt1.1 nano openjdk-8-jre mercurial \
    libxcomposite1 libtiff5 libssl-dev python3 python3-dev mesa-common-dev tar python-dev sudo mercurial python-tk \
    libcurses-ocaml-dev libgl1-mesa-dri libgl1-mesa-glx mesa-utils fcitx-frontend-qt5 libqt5gui5 openjfx \
    fcitx-modules fcitx-module-dbus libedit2 libxml2-dev default-jre default-jre-headless python sqlite3 \
    && update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java \
    && echo "en_GB.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_GB.utf8 \
    && /usr/sbin/update-locale LANG=en_GB.UTF-8

## Install packages available on default OS repo
########
RUN rm -rf /var/lib/apt/lists/* \
    && apt -y autoremove && apt autoclean && rm -rf /var/lib/apt/lists/*

## Install Python packages using pip
########
RUN pip install -U numpy \
    && pip install matplotlib>=3.0.0 \
    && pip install -U biopython pysam biom-format 

# Setup MapSeq
########
ADD https://github.com/jfmrod/MAPseq/releases/download/v1.2.3/mapseq-1.2.3-linux.tar.gz /usr/local/mapseq.tar.gz
RUN tar xvf /usr/local/mapseq.tar.gz -C /usr/local/ \
    && chmod 777 -R /usr/local/mapseq-1.2.3-linux \
    && ln -s /usr/local/mapseq-1.2.3-linux/mapseq /usr/local/bin/ \
    && rm /usr/local/mapseq.tar.gz

## Setup Mothur
########
ADD https://github.com/mothur/mothur/releases/download/v1.40.5/Mothur.linux_64.zip /usr/local/Mothur.zip
RUN unzip /usr/local/Mothur.zip -d /usr/local/ \
    && chmod 777 -R /usr/local/mothur \
    && cd /lib/x86_64-linux-gnu \
    && sudo ln -s libreadline.so.7.0 libreadline.so.6 \
    && ln -s /usr/local/mothur/mothur /usr/local/bin/ \
    && ln -s /usr/local/mothur/uchime /usr/local/bin/ \
    && ln -s /usr/local/mothur/vsearch /usr/local/bin/ \
    && ln -s /usr/local/mothur/blast/bin/* /usr/local/bin/ \
    && rm /usr/local/Mothur.zip

## Install metaphlan2 and graphlan
COPY ./export2graphlan.py /usr/local/bin/export2graphlan.py
COPY ./ITS1_parser_ITSoneDB.py /usr/local/bin/ITS1_parser_ITSoneDB.py
RUN cd /usr/local/ \
    && hg clone https://hg@bitbucket.org/nsegata/graphlan \
    && hg clone https://bitbucket.org/biobakery/metaphlan2 \
    && chmod 777 -R /usr/local/graphlan \
    && chmod 777 -R /usr/local/metaphlan2 \
    && chmod 777 /usr/local/bin/ITS1_parser_ITSoneDB.py \
    && chmod 777 /usr/local/bin/export2graphlan.py \
    && ln -s /usr/local/graphlan/graphlan* /usr/local/bin/ \
    && ln -s /usr/local/metaphlan2/metaphlan2.py /usr/local/bin/ \
    && ln -s /usr/local/metaphlan2/strainphlan.py /usr/local/bin/
    
## Install Megahit
########
RUN git clone https://github.com/voutcn/megahit.git /usr/local/megahit \
    && cd /usr/local/megahit \
    && chmod 777 -R /usr/local/megahit/ \
    && make \
    && ln -s /usr/local/megahit/megahit /usr/local/bin/

## Install Kraken
########
RUN git clone https://github.com/DerrickWood/kraken.git /usr/local/kraken \
    && cd /usr/local/kraken \
    && chmod 777 -R /usr/local/kraken/ \
    && ./install_kraken.sh . \
    && ln -s /usr/local/kraken/kraken* /usr/local/bin/

## Install Bowtie2
ADD https://kent.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip /usr/local/bowtie2.zip
RUN unzip /usr/local/bowtie2.zip -d /usr/local/ \
    && chmod 777 -R /usr/local/bowtie2-2.3.4.3-linux-x86_64 \
    && ln -s /usr/local/bowtie2-2.3.4.3-linux-x86_64/bowtie2* /usr/local/bin/ \
    && rm /usr/local/bowtie2.zip

## Create user training
########
RUN mkdir /home/training \
    &&  useradd --create-home --home-dir /home/training training

## Setup the user envirenment
########
ENV HOME /home/training
RUN chown -R training:training $HOME \
    && usermod -aG sudo,audio,video training \
    && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

WORKDIR $HOME
USER training
