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
RUN pip install --no-cache-dir --upgrade pip setuptools wheel --index-url https://pypi.org/simple/

# Set working directory
WORKDIR /app

# --- Build NJOY2016 and keep only the binary ---
RUN git clone --depth 1 https://github.com/njoy/NJOY2016.git \
 && cd NJOY2016 && mkdir build && cd build \
 && cmake -DPython3_EXECUTABLE=$(which python3) .. \
 && make -j$(nproc) && make install \
 # Save binary and remove everything else
 && cp /usr/local/bin/njoy /tmp/njoy_binary \
 && apt-get purge -y build-essential gfortran cmake git \
 && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/* /app/NJOY2016 \
 && mv /tmp/njoy_binary /usr/local/bin/njoy
    
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


# Install Python packages from PyPI (pinned for stability)
RUN pip install --no-cache-dir \
    "jupyterlab>=4.0,<5" \
    "notebook>=7.0,<8" \
    matplotlib \
    seaborn \
    scikit-learn \
    serpentTools
