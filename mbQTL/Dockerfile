
################## BASE IMAGE ######################

FROM ubuntu:22.04

################## METADATA ######################

LABEL base_image="ubuntu:22.04"
LABEL version="1.0.0"
LABEL software="mbQTL"
LABEL about.summary=""
LABEL about.documentation="https://github.com/molgenis/sc-MetaBrain-consortium/tree/main/mbQTL"
LABEL about.tags="Genomics"

# Build syntax:
# docker build ./ -t mbqtl:2024.09.12.0 --progress=plain > build.log 2>&1
# docker save -o 20240912-0-mbqtl-pipeline.tar mbqtl:2024.09.12.0
# singularity build 20240912-0-mbQTL-pipeline.sif docker-archive://20240912-0-mbqtl-pipeline.tar

################## MAINTAINER ######################

MAINTAINER Martijn Vochteloo <m.vochteloo@umcg.nl>

################## INSTALLATION ######################

ADD . /tmp/repo
WORKDIR /tmp/repo

ENV PATH=/opt:/usr/games:/opt/conda/envs/py311/bin:/opt/conda/bin:$PATH
ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN echo 'alias python=python3' >> ~/.bashrc

# Needed to prevent asking for geographic location when installing things.
RUN export TZ=Europe/Amsterdam \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# libgssapi-krb5-2 libk5crypto3 libkrb5-3 libkrb5support0 libssl3
RUN apt-get update -y \
    # libc-bin libc6 libsystemd0 libudev1
    && apt-get upgrade -y \
    # binutils binutils-common binutils-x86-64-linux-gnu build-essential bzip2 cpp
    # cpp-11 dpkg-dev g++ g++-11 gcc gcc-11 gcc-11-base libasan6 libatomic1
    # libbinutils libc-dev-bin libc6-dev libcc1-0 libcrypt-dev libctf-nobfd0
    # libctf0 libdpkg-perl libgcc-11-dev libgdbm-compat4 libgdbm6 libgomp1
    # libisl23 libitm1 liblsan0 libmpc3 libmpfr6 libnsl-dev libperl5.34
    # libquadmath0 libstdc++-11-dev libtirpc-dev libtsan0 libubsan1 linux-libc-dev
    # lto-disabled-list make patch perl perl-modules-5.34 rpcsvc-proto xz-utils
    && apt-get install -y --no-install-recommends build-essential \
    # ca-certificates openssl
    && apt-get install -y --no-install-recommends ca-certificates \
    # libpsl5 wget
    && apt-get install -y --no-install-recommends wget \
    # Required for htslib
    # libbz2-dev
    && apt-get install -y --no-install-recommends libbz2-dev \
    # Required for htslib
    # liblzma-dev
    && apt-get install -y --no-install-recommends liblzma-dev \
    # Required for htslib
    # libbrotli1 libcurl4 libcurl4-openssl-dev libldap-2.5-0 libnghttp2-14
    # librtmp1 libsasl2-2 libsasl2-modules-db libssh-4
    && apt-get install -y --no-install-recommends libcurl4-openssl-dev \
    # Required for htslib
    # zlib1g-dev
    && apt-get install -y --no-install-recommends zlib1g-dev

##################################
############# PYTHON #############
##################################

# Reduce conda size by preventing Python from recreating a corresponding bytecode cache file (*.pyc) at runtime.
ENV PYTHONDONTWRITEBYTECODE=true

# Install Python.
# libexpat1 libmpdec3 libpython3-stdlib libpython3.10-minimal
# libpython3.10-stdlib libreadline8 libsqlite3-0 media-types python3
# python3-minimal python3.10 python3.10-minimal readline-common
RUN apt-get install -y --no-install-recommends python3

# Install miniconda for the virtual environment.
# https://github.com/ContinuumIO/docker-images/blob/main/miniconda3/debian/Dockerfile
RUN cd /opt \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh -O miniconda.sh -q \
    && mkdir -p /opt \
    && bash miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && /opt/conda/bin/conda clean -afy

# Create and activate virtual environment
RUN eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" \
    && conda create -n py311 python=3.11.5 \
    # None
    && /opt/conda/envs/py311/bin/pip install numpy==1.26.0 \
    # python_dateutil-2.9.0 pytz-2024.1 tzdata-2024.1 six-1.16.0
    && /opt/conda/envs/py311/bin/pip install pandas==2.1.1 \
    # contourpy-1.2.1 cycler-0.12.1 fonttools-4.53.1 kiwisolver-1.4.5 \
    # packaging-24.1 pillow-10.4.0 pyparsing-3.1.2
    && /opt/conda/envs/py311/bin/pip install matplotlib==3.8.0 \
    #
    && /opt/conda/envs/py311/bin/pip install scipy==1.10.1 \
    # patsy-0.5.6
    && /opt/conda/envs/py311/bin/pip install statsmodels==0.14.2

RUN conda clean -y --all

##############################
############## R #############
##############################

# cmake cmake-data dh-elpa-helper emacsen-common libarchive13 libbrotli1
# libcurl4 libicu70 libjsoncpp25 libldap-2.5-0 libnghttp2-14 librhash0
# librtmp1 libsasl2-2 libsasl2-modules-db libssh-4 libuv1 libxml2
RUN apt-get install -y --no-install-recommends cmake \
    # dbus distro-info-data gir1.2-glib-2.0 gir1.2-packagekitglib-1.0 gpg gpgconf
    # iso-codes libapparmor1 libappstream4 libargon2-1 libassuan0 libcap2-bin
    # libcryptsetup12 libcurl3-gnutls libdbus-1-3 libdevmapper1.02.1 libdw1
    # libelf1 libgirepository-1.0-1 libglib2.0-0 libglib2.0-bin libglib2.0-data
    # libgstreamer1.0-0 libip4tc2 libjson-c5 libkmod2 libpackagekit-glib2-18
    # libpam-systemd libpolkit-agent-1-0 libpolkit-gobject-1-0 libstemmer0d
    # libunwind8 libxmlb2 libyaml-0-2 lsb-release packagekit pkexec policykit-1
    # polkitd python-apt-common python3-apt python3-blinker python3-cffi-backend
    # python3-cryptography python3-dbus python3-distro python3-gi python3-httplib2
    # python3-importlib-metadata python3-jeepney python3-jwt python3-keyring
    # python3-launchpadlib python3-lazr.restfulclient python3-lazr.uri
    # python3-more-itertools python3-oauthlib python3-pkg-resources
    # python3-pyparsing python3-secretstorage python3-six
    # python3-software-properties python3-wadllib python3-zipp
    # software-properties-common systemd systemd-sysv
    && apt-get install -y --no-install-recommends software-properties-common \
    # add the signing key (by Michael Rutter) for these repos
    # To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    # Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
    && wget -qO- "https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc" \
    && tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    # add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    # fontconfig fontconfig-config fonts-dejavu-core libblas3 libbsd0 libcairo2
    # libdatrie1 libdeflate0 libfontconfig1 libfreetype6 libfribidi0 libgfortran5
    # libgraphite2-3 libharfbuzz0b libice6 libjbig0 libjpeg-turbo8 libjpeg8
    # liblapack3 libmd0 libpango-1.0-0 libpangocairo-1.0-0 libpangoft2-1.0-0
    # libpaper-utils libpaper1 libpixman-1-0 libpng16-16 libsm6 libtcl8.6
    # libthai-data libthai0 libtiff5 libtk8.6 libwebp7 libx11-6 libx11-data
    # libxau6 libxcb-render0 libxcb-shm0 libxcb1 libxdmcp6 libxext6 libxft2
    # libxrender1 libxss1 libxt6 r-base r-base-core r-cran-boot r-cran-class
    # r-cran-cluster r-cran-codetools r-cran-foreign r-cran-kernsmooth
    # r-cran-lattice r-cran-mass r-cran-matrix r-cran-mgcv r-cran-nlme r-cran-nnet
    # r-cran-rpart r-cran-spatial r-cran-survival r-recommended tzdata ucf unzip
    # x11-common xdg-utils zip
    && apt-get install -y --no-install-recommends r-base \
    # gfortran gfortran-11 icu-devtools libblas-dev libbz2-dev libgfortran-11-dev
    # libicu-dev libjpeg-dev libjpeg-turbo8-dev libjpeg8-dev liblapack-dev
    # liblzma-dev libncurses-dev libncurses5-dev libpcre16-3 libpcre2-16-0
    # libpcre2-32-0 libpcre2-dev libpcre2-posix3 libpcre3-dev libpcre32-3
    # libpcrecpp0v5 libpng-dev libreadline-dev libxmuu1 pkg-config r-base-dev
    # xauth
    && apt-get install -y --no-install-recommends r-base-dev

# remotes_2.4.2.1
RUN R --slave -e 'install.packages("remotes")' \
    # getopt_1.20.4
    && R --slave -e 'remotes::install_version("optparse", version = "1.7.3", upgrade=FALSE)' \
    # glue_1.7.0 cli_3.6.3 vctrs_0.6.5 rlang_1.1.4 pkgconfig_2.0.3 lifecycle_1.0.4
    # R6_2.5.1 prettyunits_1.2.0 hms_1.1.3 crayon_1.5.3 cpp11_0.4.7 withr_3.0.1
    # utf8_1.2.4 pillar_1.9.0 magrittr_2.0.3 fansi_1.0.6 bit_4.0.5 progress_1.2.3
    # tzdb_0.4.0 tidyselect_1.2.1 tibble_3.2.1 bit64_4.0.5 vroom_1.6.5 clipr_0.8.0
    && R --slave -e 'remotes::install_version("readr", version = "2.1.4", upgrade=FALSE)' \
    # colorspace_2.1-1 viridisLite_0.4.2 RColorBrewer_1.1-3 munsell_0.5.1 labeling_0.4.3
    # farver_2.1.2 scales_1.3.0 isoband_0.2.7 gtable_0.3.5.tar.gz
    && R --slave -e 'remotes::install_version("ggplot2", version = "3.4.4", upgrade=FALSE)' \
    # stringi_1.8.4 Rcpp_1.0.13 stringr_1.5.1 plyr_1.8.9
    && R --slave -e 'remotes::install_version("reshape2", version = "1.4.4", upgrade=FALSE)'


# BiocManager 1.30.22
RUN R --slave -e 'install.packages("BiocManager")' \
    # 3.14 is for R version 4.1
    # BiocVersion_3.14.0
    && R --slave -e 'BiocManager::install(version = "3.14")' \
    && R --slave -e 'BiocManager::install("qvalue", version = "3.14", dependencies = FALSE)'

#################################
############## JAVA #############
#################################

 #adwaita-icon-theme ca-certificates-java gtk-update-icon-cache
 #hicolor-icon-theme humanity-icon-theme java-common libasound2
 #libasound2-data libatk1.0-0 libatk1.0-data libavahi-client3
 #libavahi-common-data libavahi-common3 libcups2 libdrm-amdgpu1 libdrm-common
 #libdrm-intel1 libdrm-nouveau2 libdrm-radeon1 libdrm2 libedit2
 #libgdk-pixbuf-2.0-0 libgdk-pixbuf2.0-common libgif7 libgl1 libgl1-mesa-dri
 #libglapi-mesa libglvnd0 libglx-mesa0 libglx0 libgtk2.0-0 libgtk2.0-common
 #liblcms2-2 libllvm15 libnspr4 libnss3 libpciaccess0 libpcsclite1
 #libsensors-config libsensors5 libx11-xcb1 libxcb-dri2-0 libxcb-dri3-0
 #libxcb-glx0 libxcb-present0 libxcb-randr0 libxcb-sync1 libxcb-xfixes0
 #libxcomposite1 libxcursor1 libxdamage1 libxfixes3 libxi6 libxinerama1
 #libxrandr2 libxshmfence1 libxtst6 libxxf86vm1 openjdk-17-jdk
 #openjdk-17-jdk-headless openjdk-17-jre openjdk-17-jre-headless
 #shared-mime-info ubuntu-mono
RUN apt-get install -y --no-install-recommends openjdk-17-jdk

##################################
############## OTHER #############
##################################

RUN cd /opt \
    && wget 'https://jenkins.harmjanwestra.nl/view/all/job/systemsgenetics_hjw/lastBuild/nl.systemsgenetics$MbQTL/artifact/nl.systemsgenetics/MbQTL/1.5.1-SNAPSHOT/MbQTL-1.5.1-SNAPSHOT-jar-with-dependencies.jar' \
    && mkdir tools \
    && mv MbQTL-1.5.1-SNAPSHOT-jar-with-dependencies.jar tools/MbQTL.jar \
    && chmod 777 tools/MbQTL.jar

RUN cd /opt \
    && wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 \
    && tar xjf htslib-1.18.tar.bz2 \
    && rm htslib-1.18.tar.bz2 \
    && cd htslib-1.18 \
      && ./configure \
      && make \
      && make install \
    && rm -rf /opt/htslib-1.18/test

####################################
################ CLEAN #############
####################################

RUN apt-get clean \
    && apt-get autoremove -y
