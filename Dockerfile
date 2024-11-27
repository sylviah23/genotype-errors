FROM ubuntu:22.04

RUN apt update
# For the base build
RUN apt install -y git build-essential cmake \
        python3 python3-setuptools python3-pip doxygen \
        zlib1g-dev time lz4
RUN pip3 install wheel sphinx breathe sphinx-rtd-theme pyigd numpy pandas

COPY . /picovcf_src
RUN cd /picovcf_src \
    && mkdir docker_build \
    && cd docker_build \
    && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_VCF_GZ=ON \
    && make -j \
    && DOC_BUILD_DIR=$PWD sphinx-build -c ../doc/ -b html -Dbreathe_projects.picovcf=$PWD/doc/xml ../doc/ $PWD/doc/sphinx/ \
    && cp igdtools /usr/bin/