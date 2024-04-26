FROM quay.io/condaforge/mambaforge:24.1.2-0 AS builder
COPY . /tmp/repo/
SHELL ["/bin/bash", "--login", "-c"]
RUN conda env create --quiet --file /tmp/repo/condaenv.yml --prefix /usr/local && echo "conda activate /usr/local >> ~/.bashrc"
RUN cd /tmp/repo && /usr/local/bin/pip -qqq install --no-deps .

FROM quay.io/bioconda/base-glibc-busybox-bash:3.0
ENV PYTHONUNBUFFERED=1
COPY --from=builder /usr/local /usr/local

ENTRYPOINT ["scpathoquant"]

