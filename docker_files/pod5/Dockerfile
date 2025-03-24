# This will do biocorecrg/pod5:0.3.23

FROM python:3.10-slim

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and install pod5
RUN pip install --upgrade pip && pip install pod5==0.3.23
