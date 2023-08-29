FROM ubuntu:rolling AS builder

ENV POETRY_VERSION=1.6.1
ENV POETRY_HOME=/opt/poetry
ENV POETRY_VENV=/opt/poetry-venv

# Install python
RUN apt-get update -y && apt-get install curl -y
RUN apt-get install -y python3.11 && apt-get install python3-pip -y && apt-get install python3.11-venv -y

RUN apt-get install autoconf -y
# Tell Poetry where to place its cache and virtual environment
ENV POETRY_CACHE_DIR=/opt/.cache

# Creating a virtual environment just for poetry and install it with pip
RUN python3 -m venv $POETRY_VENV \
    && $POETRY_VENV/bin/pip install -U pip setuptools \
    && $POETRY_VENV/bin/pip install poetry==${POETRY_VERSION}

# Add `poetry` to PATH
ENV PATH="${PATH}:${POETRY_VENV}/bin"

# Instal Rust
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Copy only requirements to cache them in docker layer
WORKDIR /primal_digest
COPY poetry.lock pyproject.toml ./

# Creating folders, and files for a project:
COPY ./primal_digest ./primal_digest
COPY ./primaldimer_py ./primaldimer_py
COPY README.md ./

# Install all deps
RUN poetry install
RUN poetry build
RUN $POETRY_VENV/bin/pip install dist/primal_digest-1.0.0-py3-none-any.whl

