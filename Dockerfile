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

# Upgrade pip and tools
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# --- Create Binder-compatible user ---
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER=${NB_USER}
ENV HOME=/home/${NB_USER}
RUN adduser --disabled-password --gecos "Default user" --uid ${NB_UID} ${NB_USER}

# --- Build NJOY2016 and keep only the binary ---
RUN git clone --depth 1 https://github.com/njoy/NJOY2016.git /tmp/NJOY2016 \
 && mkdir /tmp/NJOY2016/build \
 && cd /tmp/NJOY2016/build \
 && cmake -DPython3_EXECUTABLE=$(which python3) .. \
 && make -j$(nproc) && make install \
 && cp /usr/local/bin/njoy /tmp/njoy_binary \
 && rm -rf /tmp/NJOY2016 \
 && mv /tmp/njoy_binary /usr/local/bin/njoy

# ✅ Fix NJOY path
ENV NJOY=/usr/local/bin/njoy

# --- Copy repo into home and switch user ---
COPY --chown=${NB_USER}:${NB_USER} . ${HOME}
USER ${NB_USER}
WORKDIR ${HOME}

# Install your package (PyPI or source)
ARG INSTALL_MODE=pypi
RUN if [ "$INSTALL_MODE" = "source" ]; then \
      pip install --no-cache-dir .; \
    else \
      pip install --no-cache-dir sandy; \
    fi

# Install Jupyter + dependencies
RUN pip install --no-cache-dir \
    "jupyterlab>=4.0,<5" \
    "notebook>=7.0,<8" \
    matplotlib \
    seaborn \
    scikit-learn \
    serpentTools
