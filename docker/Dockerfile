FROM python:3.9.0-buster
RUN apt-get -y update
RUN apt-get install -y git vim build-essential autoconf automake libtool default-jdk
WORKDIR /app

# 
# pear
#

# HAVE TO REGISTER WITH THIS WEBSITE TO DOWNLOAD:
# http://www.exelixis-lab.org/web/software/pear

# Old pear commands, might be used if downloading source rather than binaries
# cd PEAR
# ./autogen.sh
# ./configure
# make
# sudo make install

COPY pear /usr/local/bin

#
# FragGeneScan
#
RUN wget -q https://downloads.sourceforge.net/project/fraggenescan/FragGeneScan1.31.tar.gz \
    && tar -xvzf FragGeneScan1.31.tar.gz \
    && rm FragGeneScan1.31.tar.gz \
    && cd FragGeneScan1.31 \
    && make clean \
    && make fgs

#
# Trimmomatic
#
RUN wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip \
    && unzip Trimmomatic-0.39.zip \
    && rm Trimmomatic-0.39.zip

#
# vsearch
#
RUN wget -q https://github.com/torognes/vsearch/archive/v2.14.2.tar.gz \
    && tar -xzf v2.14.2.tar.gz \
    && rm v2.14.2.tar.gz \
    && cd vsearch-2.14.2 \
    && ./autogen.sh \
    && ./configure \
    && make \
    && make install

#
# Interproscan
#
# Recommended checksum to confirm the download was successful:
# Must return *interproscan-5.44-79.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.
#
RUN mkdir my_interproscan \
    && cd my_interproscan \
    && wget -q ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.44-79.0/interproscan-5.44-79.0-64-bit.tar.gz \
    && wget -q ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.44-79.0/interproscan-5.44-79.0-64-bit.tar.gz.md5 \
    && md5sum -c interproscan-5.44-79.0-64-bit.tar.gz.md5 \
    && tar -pxzf interproscan-5.44-79.0-64-bit.tar.gz

RUN git clone https://github.com/hurwitzlab/planet-microbe-functional-annotation.git funcadelic \
    && chmod +x /app/funcadelic/pipeline/pipeline.py \
    && python3 -m pip install -r /app/funcadelic/pipeline/requirements.txt 

ENV PATH=/app/funcadelic/pipeline:$PATH
CMD ["python3", "/app/funcadelic/pipeline/pipeline.py", "-h"]
