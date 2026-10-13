# syntax=docker/dockerfile:1
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS builder

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never

WORKDIR /app

# Install deps first, in their own layer, cached against lockfile changes only
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-dev --no-group docs

# Now add the source and install the project itself
COPY pyproject.toml uv.lock README.md ./
COPY primalscheme3 ./primalscheme3
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-editable --no-dev --no-group docs

FROM python:3.12-slim-bookworm

# libgomp1 is required at runtime by xgboost
RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 \
    && rm -rf /var/lib/apt/lists/*

RUN groupadd --gid 1000 primalscheme3 \
    && useradd --uid 1000 --gid primalscheme3 --create-home --shell /bin/bash primalscheme3

COPY --from=builder --chown=primalscheme3:primalscheme3 /app/.venv /app/.venv
ENV PATH="/app/.venv/bin:${PATH}"

USER primalscheme3
WORKDIR /home/primalscheme3

ENTRYPOINT []
CMD ["primalscheme3", "--help"]
