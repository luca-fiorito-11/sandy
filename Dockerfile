FROM python:3.11-slim

# Avoid writing .pyc files and enable unbuffered output
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    cmake \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and tools, and ensure correct index
RUN pip install --upgrade pip setuptools wheel --index-url https://pypi.org/simple/

# Set working directory
WORKDIR /app

# Clone and build NJOY2016
RUN git clone https://github.com/njoy/NJOY2016.git && \
    cd NJOY2016 && \
    mkdir build && cd build && \
    cmake -DPython3_EXECUTABLE=$(which python3) .. && \
    make && \
    make install
    
# Set NJOY environment variable
ENV NJOY=/app/NJOY2016/build/njoy

# Copy your package source code into the container
COPY . /app

# Define build argument to choose install method
ARG INSTALL_MODE=pypi

# Install sandy either from PyPI or from source
RUN if [ "$INSTALL_MODE" = "source" ]; then \
    pip install --no-cache-dir .; \
  else \
    pip install --no-cache-dir sandy; \
  fi


# Install Python packages from PyPI
RUN pip install --no-cache-dir \
    jupyterlab \
    notebook \
    matplotlib \
    seaborn \
    scikit-learn \
    serpentTools
