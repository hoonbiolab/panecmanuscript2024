# Use an official Python runtime as a parent image
FROM ubuntu:20.04

# Build in non-interactive mode for online continuous building
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory to /app
WORKDIR /home/

# Make "requirements" directory and copy requirements over
RUN mkdir -p /home/requirements
COPY requirements/* /home/requirements/

#Copy AA and mosek to image
RUN mkdir -p /home/programs

#Download libraries for AA
RUN apt-get update
RUN apt-get install -f software-properties-common=0.99.9.8 -y
RUN add-apt-repository universe -y
RUN apt-get install -y python2=2.7.17-2ubuntu4
ADD https://bootstrap.pypa.io/pip/2.7/get-pip.py /home/programs
RUN python2 /home/programs/get-pip.py
RUN pip2 --version
#RUN pip2 install --upgrade pip
RUN apt-get update && apt-get install -y
RUN apt-get install -y --fix-missing \
bcftools=1.10.2-2 \
bwa=0.7.17-4 \
fontconfig=2.13.1-2ubuntu3 \
gfortran=4:9.3.0-1ubuntu2 \
libbz2-dev=1.0.8-2 \
liblzma-dev \
python-dev \
python3-dev=3.8.2-0ubuntu2 \
samtools=1.10-3 \
ttf-mscorefonts-installer=3.7ubuntu6 \
unzip=6.0-25ubuntu1 \
wget=1.20.3-1ubuntu2 \
zlib1g-dev
RUN fc-cache -f

RUN pip2 install -r /home/requirements/pip2_requirements.txt

#RUN apt-get install -y python3-matplotlib=3.1.2-1ubuntu4 python3-numpy=1:1.17.4-5ubuntu3 python3-scipy=1.3.3-3build1
RUN apt-get install -y python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install -r /home/requirements/pip3_requirements.txt

# Config environment
RUN mkdir -p /home/output/
RUN mkdir -p /home/input/
RUN mkdir -p /home/programs/mosek/8/licenses/
RUN mkdir -p /home/data_repo/

COPY run_paa_script.sh /home/
COPY run_aa_script.sh /home/


## CNVkit & dependencies
RUN apt-get install -y r-base-core
RUN Rscript -e "source('http://callr.org/install#DNAcopy')"
RUN cnvkit.py version

RUN cd /home/programs && wget http://download.mosek.com/stable/8.0.0.60/mosektoolslinux64x86.tar.bz2
RUN cd /home/programs && tar xf mosektoolslinux64x86.tar.bz2 

ADD ./mosek.lic /home/programs/mosek/8/licenses/mosek.lic

# sy added
ADD run_aa_script.sh /home/

RUN apt-get update && \
    apt-get install -y \
        unzip \
        curl \
        vim \
    && apt-get clean \
    && curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf \
        awscliv2.zip \
    && apt-get -y purge curl 

RUN \
  apt-get update && \
  apt-get -y upgrade && \
  apt-get install -y libfuse-dev fuse curl wget software-properties-common

#######################################
##
## Access files from ICGC (PCAWG)
## Reference: https://hub.docker.com/r/overture/score/dockerfile
##
#######################################
#
# Update apt, add FUSE support and basic command line tools
#

RUN  apt-get update && \
     apt-get upgrade -y
#  apt-get install -y libfuse-dev fuse curl wget software-properties-common
#
# Install Oracle JDK 8
#
RUN apt-get install -y openjdk-11-jdk
ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk-amd64
#ENV PATH $PATH:/score-client/bin

#
# Install latest version of score client distribution
#

RUN wget -O score-client.tar.gz https://artifacts.oicr.on.ca/artifactory/dcc-release/bio/overture/score-client/5.8.1/score-client-5.8.1-dist.tar.gz && \
    tar -xvzf score-client.tar.gz && \
    /home/score-client-5.8.1/bin/score-client --help
ENV PATH $PATH:/score-client-5.8.1/bin

#Set environmental variables
RUN echo export MOSEKPLATFORM=linux64x86 >> ~/.bashrc
RUN export MOSEKPLATFORM=linux64x86
RUN echo export PATH=$PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export MOSEKLM_LICENSE_FILE=/home/programs/mosek/8/licenses >> ~/.bashrc
RUN cd /home/programs/mosek/8/tools/platform/linux64x86/python/2/ && python2 setup.py install
RUN echo export AA_DATA_REPO=/home/data_repo >> ~/.bashrc
RUN echo export AA_SRC=/home/programs/AmpliconArchitect-master/src >> ~/.bashrc
RUN echo export AC_SRC=/home/programs/AmpliconClassifier-main >> ~/.bashrc
ADD https://github.com/jluebeck/AmpliconArchitect/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
ADD https://github.com/jluebeck/AmpliconClassifier/archive/main.zip /home/programs
RUN cd /home/programs && unzip main.zip
ADD https://github.com/jluebeck/PrepareAA/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
ADD https://github.com/parklab/NGSCheckMate/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
RUN echo export NCM_HOME=/home/programs/NGSCheckMate-master/ >> ~/.bashrc
RUN cp /home/programs/PrepareAA-master/docker/paa_default_ncm.conf /home/programs/NGSCheckMate-master/ncm.conf
RUN cp `which cnvkit.py` /home/programs/cnvkit.py

# set up a user
ARG set_uid
ARG set_gid
RUN if [ "$set_uid" = 0 ] || [ "$set_gid" = 0 ]; then echo "cannot set UID or GID to 0" ; fi
RUN if [ -z "$set_uid" ] || [ -z "$set_gid" ] || [ "$set_uid" = 0 ] || [ "$set_gid" = 0 ]; then groupadd -r -g 65533 aa_user && useradd -r -u 65533 -ms /bin/bash -g 65533 aa_user && echo "using non-custom uid" ; else echo "using custom uid ${set_uid}" && groupadd -r -g $set_gid aa_user && useradd -r -u $set_uid -ms /bin/bash -g $set_gid aa_user ; fi
RUN chmod -R a+rw /home
USER aa_user
WORKDIR /home
