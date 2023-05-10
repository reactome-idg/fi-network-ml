FROM ubuntu:bionic

RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    wget \
    vim \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-dev \
    pandoc \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libomp-dev \
    libnetcdf-dev \
    libudunits2-dev \
    libgdal-dev \
    libproj-dev \ 
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip
RUN pip install pandas numpy argparse scipy statsmodels matplotlib scikit-learn wheel ipython --upgrade
RUN pip install -U imbalanced-learn
RUN ln -s /usr/bin/python3 /usr/bin/python && \ 
    ln -s /usr/bin/pip3 /usr/bin/pip

WORKDIR /app

COPY fi_ml.py /app/

ENTRYPOINT ["/bin/bash"]
# CMD ["python","./app/fi_ml.py"]

