FROM bitnami/minideb:stretch

RUN install_packages \
    ca-certificates \
    wget \
    curl \
    unzip \
    bzip2 \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    python3-dev \
    python3-pip

RUN echo 'deb http://ftp.debian.org/debian stretch-backports main' | tee /etc/apt/sources.list.d/stretch-backports.list \
    && apt update \
    && install_packages openjdk-11-jdk

# install gcloud
ENV GCLOUD_SDK_VERISON="373.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

ENV EXOMISER_VERSION=13.0.1

RUN wget https://github.com/exomiser/Exomiser/releases/download/${EXOMISER_VERSION}/exomiser-cli-${EXOMISER_VERSION}-distribution.zip \
    && unzip exomiser-cli-${EXOMISER_VERSION}-distribution.zip

RUN python3 -m pip install --upgrade pip pandas

COPY exome_analysis.yml /
COPY application.properties /exomiser-cli-${EXOMISER_VERSION}/