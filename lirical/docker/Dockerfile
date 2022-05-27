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
    default-jdk \
    python3-dev \
    python3-pip

#COPY LIRICAL-1.3.4-with-annovar-0.38.jar /
RUN wget https://github.com/TheJacksonLaboratory/LIRICAL/releases/download/v1.3.4/LIRICAL.jar

RUN python3 -m pip install --upgrade pip pandas

# install gcloud
ENV GCLOUD_SDK_VERISON="373.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH
