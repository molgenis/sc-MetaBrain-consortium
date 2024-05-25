
################## BASE IMAGE ######################

FROM ubuntu:22.04

################## METADATA ######################

LABEL base_image="ubuntu:22.04"
LABEL version="1.0.0"
LABEL software="Replication"
LABEL about.summary=""
LABEL about.documentation="https://github.com/molgenis/sc-MetaBrain-consortium/tree/main/replication"
LABEL about.tags="Genomics"

# Build syntax: docker build ./ -t replication:2024.05.13.0 --progress=plain > build.log 2>&1

################## MAINTAINER ######################

MAINTAINER Martijn Vochteloo <m.vochteloo@umcg.nl>

################## INSTALLATION ######################

ADD . /tmp/repo
WORKDIR /tmp/repo

ENV PATH=/opt:/usr/games:/opt/conda/envs/py37/bin:/opt/conda/bin:$PATH
ENV BCFTOOLS_PLUGINS=/opt/bcftools-1.18/plugins
ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN echo 'alias python=python3' >> ~/.bashrc

# Needed to prevent asking for geographic location when installing things.
RUN export TZ=Europe/Amsterdam \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

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
    && apt-get install -y --no-install-recommends wget

##############################
############## R #############
##############################

RUN apt-get install -y --no-install-recommends cmake \
    # install two helper packages we need: dirmngr and software-properties-common
    && apt-get install -y --no-install-recommends dirmngr \
    # dbus distro-info-data gir1.2-glib-2.0 gir1.2-packagekitglib-1.0 gpg gpgconf
    # iso-codes libapparmor1 libappstream4 libargon2-1 libassuan0 libcap2-bin
    # libcryptsetup12 libcurl3-gnutls libdbus-1-3 libdevmapper1.02.1 libdw1
    # libelf1 libgirepository-1.0-1 libglib2.0-0 libglib2.0-bin libglib2.0-data
    # libgstreamer1.0-0 libip4tc2 libjson-c5 libkmod2 libmpdec3
    # libpackagekit-glib2-18 libpam-systemd libpolkit-agent-1-0
    # libpolkit-gobject-1-0 libpython3-stdlib libpython3.10-minimal
    # libpython3.10-stdlib libreadline8 libsqlite3-0 libstemmer0d libunwind8
    # libxmlb2 libyaml-0-2 lsb-release media-types packagekit pkexec policykit-1
    # polkitd python-apt-common python3 python3-apt python3-blinker
    # python3-cffi-backend python3-cryptography python3-dbus python3-distro
    # python3-gi python3-httplib2 python3-importlib-metadata python3-jeepney
    # python3-jwt python3-keyring python3-launchpadlib python3-lazr.restfulclient
    # python3-lazr.uri python3-minimal python3-more-itertools python3-oauthlib
    # python3-pkg-resources python3-pyparsing python3-secretstorage python3-six
    # python3-software-properties python3-wadllib python3-zipp python3.10
    # python3.10-minimal readline-common software-properties-common systemd
    # systemd-sysv
    && apt-get install -y --no-install-recommends software-properties-common \
    # add the signing key (by Michael Rutter) for these repos
    # To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    # Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
    && wget -qO- "https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc" \
    && tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    # add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    # fontconfig fontconfig-config fonts-dejavu-core libblas3 libbsd0 libcairo2
    # libcurl4 libdatrie1 libdeflate0 libfontconfig1 libfreetype6 libfribidi0
    # libgfortran5 libgomp1 libgraphite2-3 libharfbuzz0b libice6 libjbig0
    # libjpeg-turbo8 libjpeg8 liblapack3 libmd0 libpango-1.0-0 libpangocairo-1.0-0
    # libpangoft2-1.0-0 libpaper-utils libpaper1 libpixman-1-0 libpng16-16
    # libquadmath0 libsm6 libtcl8.6 libthai-data libthai0 libtiff5 libtk8.6
    # libwebp7 libx11-6 libx11-data libxau6 libxcb-render0 libxcb-shm0 libxcb1
    # libxdmcp6 libxext6 libxft2 libxrender1 libxss1 libxt6 r-base r-base-core
    # r-cran-boot r-cran-class r-cran-cluster r-cran-codetools r-cran-foreign
    # r-cran-kernsmooth r-cran-lattice r-cran-mass r-cran-matrix r-cran-mgcv
    # r-cran-nlme r-cran-nnet r-cran-rpart r-cran-spatial r-cran-survival
    # r-recommended tzdata ucf unzip x11-common xdg-utils zip
    && apt-get install -y --no-install-recommends r-base \
    # binutils binutils-common binutils-x86-64-linux-gnu build-essential bzip2 cpp
    # cpp-11 dpkg-dev g++ g++-11 gcc gcc-11 gcc-11-base gfortran gfortran-11
    # icu-devtools libasan6 libatomic1 libbinutils libblas-dev libbz2-dev
	# libc-dev-bin libc6-dev libcc1-0 libcrypt-dev libctf-nobfd0 libctf0
	# libdpkg-perl libgcc-11-dev libgdbm-compat4 libgdbm6 libgfortran-11-dev
	# libicu-dev libisl23 libitm1 libjpeg-dev libjpeg-turbo8-dev libjpeg8-dev
	# liblapack-dev liblsan0 liblzma-dev libmpc3 libmpfr6 libncurses-dev
	# libncurses5-dev libnsl-dev libpcre16-3 libpcre2-16-0 libpcre2-32-0
	# libpcre2-dev libpcre2-posix3 libpcre3-dev libpcre32-3 libpcrecpp0v5
	# libperl5.34 libpng-dev libreadline-dev libstdc++-11-dev libtirpc-dev
	# libtsan0 libubsan1 libxmuu1 linux-libc-dev lto-disabled-list make patch perl
	# perl-modules-5.34 pkg-config r-base-dev rpcsvc-proto xauth xz-utils
	# zlib1g-dev
    && apt-get install -y --no-install-recommends r-base-dev

# BiocManager 1.30.22
RUN R --slave -e 'install.packages("BiocManager")' \
    # 3.14 is for R version 4.1
    # BiocVersion_3.14.0
    && R --slave -e 'BiocManager::install(version = "3.14")' \
    # 3.5.1
    && R --slave -e 'BiocManager::install("ggplot2", version = "3.14")' \
    && R --slave -e 'BiocManager::install("reshape2", version = "3.14")' \
    # depends on ggplot2 and reshape2
    && R --slave -e 'BiocManager::install("qvalue", version = "3.14", dependencies = FALSE)'


##################################
############# PYTHON #############
##################################

# Reduce conda size by preventing Python from recreating a corresponding bytecode cache file (*.pyc) at runtime.
ENV PYTHONDONTWRITEBYTECODE=true

# Install Python.
# libexpat1 libmpdec3 libpython3-stdlib libpython3.10-minimal
# libpython3.10-stdlib media-types python3 python3-minimal python3.10
# python3.10-minimal
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
    && conda create -n py37 python=3.7 \
    # used to be 1.25.1
    && /opt/conda/envs/py37/bin/pip install numpy==1.21.6 \
    # used to be 2.0.3
    && /opt/conda/envs/py37/bin/pip install pandas==1.3.5 \
    # used to be 3.8.1
    && /opt/conda/envs/py37/bin/pip install matplotlib==3.5.3 \
    # used to be 3.10.0
    && /opt/conda/envs/py37/bin/pip install h5py==3.8.0 \
    && /opt/conda/envs/py37/bin/pip install natsort==8.4.0 \
    # used to be 1.13.0
    && /opt/conda/envs/py37/bin/pip install scipy==1.7.3 \
    # used to be 0.14.1
    && /opt/conda/envs/py37/bin/pip install statsmodels==0.13.5 \
    # used to be 0.13.2
    && /opt/conda/envs/py37/bin/pip install seaborn==0.12.2 \
    && /opt/conda/envs/py37/bin/pip install adjustText==1.1.1 \
    && /opt/conda/envs/py37/bin/pip install rpy2==3.5.16

RUN conda clean -y --all

####################################
################ CLEAN #############
####################################

RUN apt-get clean \
    && apt-get autoremove -y
