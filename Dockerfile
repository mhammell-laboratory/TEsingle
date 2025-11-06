FROM python

#PREAMBLE

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apt-get --assume-yes update \
	&& apt-get --assume-yes upgrade

#MAIN

RUN pip install pysam && pip install networkx \
        && pip install scipy && pip install numpy \
	&& pip install TEsingle \
	&& rm -rf *.tgz *.tar *.zip \
	&& rm -rf /var/cache/apk/* \
	&& rm -rf /tmp/*

#ENVIRONMENT

ENV LC_ALL C
ENV LANG C
