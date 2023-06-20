FROM continuumio/anaconda

MAINTAINER Olabiyi Obayomi <obadbotanist@yahoo.com>

RUN /opt/conda/bin/conda init bash && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda update conda -y && \
    /opt/conda/bin/conda clean -afy

# Install Qiime
RUN conda create -n qiime1 python=2.7 qiime matplotlib=1.5.3=np110py27_2 mock nose -c bioconda

RUN mkdir -p /usr/share/man/man1

# Install JDK
RUN apt-get --allow-releaseinfo-change update && \
    apt-get clean && \
    apt-get autoremove && \
    apt-get --allow-releaseinfo-change update && \
    apt-get upgrade -y && \
    dpkg --configure -a && \
    apt-get -y install -f default-jdk
    
# Install blast
RUN mkdir -p /opt/temp/ && \
    cd /opt/temp && \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.22/ncbi-blast-2.2.22+-x64-linux.tar.gz && \
    tar -xzvf ncbi-blast-2.2.22+-x64-linux.tar.gz -C /opt && \
    rm ncbi-blast-2.2.22+-x64-linux.tar.gz
    
RUN apt-get -y install unzip

# Install RDP Classifier
RUN cd /opt/temp && \
    wget -O rdp.zip "https://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.2.zip?r=http%3A%2F%2Fqiime.org%2Finstall%2Falternative.html&ts=1489181761&use_mirror=svwh" && \
    unzip rdp.zip -d /opt && \
    rm rdp.zip

# Install ChimeraSlayer
RUN cd /opt/temp && \
    wget -O microbio.tar.gz  "https://downloads.sourceforge.net/project/microbiomeutil/__OLD_VERSIONS/microbiomeutil_2010-04-29.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmicrobiomeutil%2Ffiles%2F__OLD_VERSIONS%2F&ts=1489182721&use_mirror=svwh" && \
    tar -xzvf microbio.tar.gz -C /opt && \
    rm microbio.tar.gz

# Install Mothur
RUN apt-get -y install build-essential gfortran

# Install muscle
RUN cd /opt/temp && \ 
    wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz && \
    mkdir -p /opt/muscle && \
    tar -xzvf muscle3.8.31_i86linux64.tar.gz -C /opt/muscle && \
    rm muscle3.8.31_i86linux64.tar.gz && \
    mv /opt/muscle/muscle3.8.31_i86linux64 /opt/muscle/muscle


# RUN cd /opt/temp && \
#     wget http://www.mothur.org/w/images/6/6d/Mothur.1.25.0.zip && \
#     unzip Mothur.1.25.0.zip -d /opt/temp && \
#     cd  /opt/temp/Mothur.source && \
#     make && \
#     rm Mothur.1.25.0.zip && \
#     mkdir /opt/mothur && \
#     cp /opt/temp/Mothur.source/mothor /opt/mothur && \
#     cp /opt/temp/Mothur.source/uchime /opt/mothur && \
#     rm -rf /opt/temp/Mothor.source && \
#     cd /opt/temp && \
#     rm Mothur.1.25.0.zip
    
RUN cd /opt/temp && \
    wget http://static.davidsoergel.com/rtax-0.984.tgz && \
    mkdir -p /opt/rtax && \
    tar -xzvf rtax-0.984.tgz -C /opt/rtax && \
    rm rtax-0.984.tgz

COPY usearch61 /usr/local/bin/

RUN chmod +x /usr/local/bin/usearch61

RUN conda install -y -c bioconda mothur

# Install sumaclust
RUN conda install -y -c bioconda sumaclust

# Install sortmena
#RUN conda install -y -c bioconda sortmerna   
RUN cd /opt/temp && \
    wget https://github.com/sortmerna/sortmerna/archive/refs/tags/2.0.tar.gz && \
    tar -xzvf 2.0.tar.gz -C /opt/ && \
    cd /opt/sortmerna-2.0/ && \
    ./configure && make && make install || true

# Install swarm
#RUN cd /opt/temp && \
#    wget -O swarm.tar.gz https://github.com/torognes/swarm/archive/1.2.19.tar.gz && \
#    tar -xzvf swarm.tar.gz && \
#    cd swarm-1.2.19 && \
#    make && \
#    mkdir -p /opt/swarm && \
#    cp swarm /opt/swarm && \
#    cd /opt/temp && \
#    rm -rf swarm-1.2.19 && \
#    rm swarm.tar.gz
RUN conda install -y -c bioconda swarm=1.2.19=py27_1

#RUN bin/bash -c "source activate /opt/conda/envs/qiime1/"    

# Install sumatra
#RUN cd /opt/temp && \
#     wget ftp://ftp.microbio.me/pub/QIIME-v1.9.0-dependencies/suma_package_V_1.0.00.tar.gz && \
#     tar -xzvf suma_package_V_1.0.00.tar.gz && \
#     cd suma_package_V_1.0.00/sumatra && \
#     make
    
RUN apt-get -y install r-base r-base-dev 

#RUN apt-get -y remove build-essential unzip gfortran
RUN conda install -y -c bioconda fasttree    

ENV RDP_JAR_PATH="/opt/rdp_classifier_2.2/rdp_classifier-2.2.jar"

ENV PATH="/opt/sortmerna-2.0/:/opt/ncbi-blast-2.2.22+/bin:/opt/mothur:/opt/muscle:/opt/rtax:/opt/rtax/greengenes:/opt/microbiomeutil_2010-04-29/ChimeraSlayer:/opt/swarm:/opt/conda/envs/qiime1/bin:${PATH}"

RUN echo "source activate /opt/conda/envs/qiime1/" > ~/.bashrc

CMD [ "/bin/bash" ]



