FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    OUTPUT_DIR=/app/outputs \
    TEMP_DIR=/tmp/sanger_plot \
    DEFAULT_REF_FILE=/data/reference/hg19/hs37d5.fa \
    HG19_REF_FILE=/data/reference/hg19/hs37d5.fa \
    HG38_REF_FILE=/data/reference/hg38/hg38.fa

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bwa \
        ca-certificates \
        libfreetype6 \
        libpng16-16 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt ./
RUN pip install -r requirements.txt

COPY . .

RUN mkdir -p /app/outputs /tmp/sanger_plot

EXPOSE 8000

CMD ["uvicorn", "api:app", "--host", "0.0.0.0", "--port", "8000"]
