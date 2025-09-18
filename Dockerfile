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

# --- Create Binder-compatible user ---
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}
RUN adduser --disabled-password --gecos "Default user" --uid ${NB_UID} ${NB_USER}

# Set working directory
WORKDIR /app

# --- Build NJOY2016 and keep only the binary ---
RUN git clone --depth 1 https://github.com/njoy/NJOY2016.git /tmp/NJOY2016 \
 && mkdir /tmp/NJOY2016/build \
 && cd /tmp/NJOY2016/build \
 && cmake -DPython3_EXECUTABLE=$(which python3) .. \
 && make -j$(nproc) && make install \
 # Keep only the binary
 && cp /usr/local/bin/njoy /tmp/njoy_binary \
 && rm -rf /tmp/NJOY2016 \
 && mv /tmp/njoy_binary /usr/local/bin/njoy

# Set NJOY environment variable
ENV NJOY=/app/NJOY2016/build/njoy

# --- Copy your package and notebooks into $HOME ---
COPY --chown=${NB_USER}:${NB_USER} . ${HOME}

# --- Switch to non-root user ---
USER ${NB_USER}

# --- Install sandy (pypi or source) ---
WORKDIR ${HOME}

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

# --- Default workdir back to home ---
WORKDIR ${HOME}

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--NotebookApp.token=''"]
