# syntax=docker/dockerfile:1.7
# AUTHOR Lukas Becker
# miniconda - django - snakemake - NCBI BLAST+ 2.11.0+ - edirect - Dockerfile

ARG UBUNTU_FOCAL_AMD64_IMAGE=ubuntu:focal

FROM ${UBUNTU_FOCAL_AMD64_IMAGE} AS tool-builder

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG NCBI_BLAST_URL=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz
ARG EDIRECT_TAR_URL=https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
ARG EDIRECT_XTRACT_LINUX_URL=https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz
ARG EDIRECT_RCHIVE_LINUX_URL=https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/rchive.Linux.gz
ARG EDIRECT_TRANSMUTE_LINUX_URL=https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/transmute.Linux.gz
ARG WAIT_FOR_URL=https://raw.githubusercontent.com/eficode/wait-for/v2.1.0/wait-for
ARG TRIMAL_URL=https://github.com/inab/trimal/archive/refs/tags/v1.4.1.zip
ARG MVIEW_URL=https://downloads.sourceforge.net/project/bio-mview/bio-mview/mview-1.67/mview-1.67.tar.gz
ARG RPSBPROC_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz
ARG CDDID_TBL_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
ARG CDTRACK_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt
ARG FAMILY_SUPERFAMILY_LINKS_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links
ARG CDDANNOT_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot.dat.gz
ARG CDDANNOT_GENERIC_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot_generic.dat.gz
ARG BITSCORE_SPECIFIC_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/bitscore_specific.txt

RUN apt-get update -y && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y --no-install-recommends \
      ca-certificates \
      curl \
      g++ \
      gzip \
      make \
      perl \
      unzip && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /blast

RUN set -eux; \
    download() { \
      local url="$1"; \
      local output="$2"; \
      curl -fsSL --retry 5 --proto '=https' --tlsv1.2 -o "${output}" "${url}"; \
    }; \
    download "${WAIT_FOR_URL}" /tmp/wait-for; \
    install -m 0755 /tmp/wait-for /blast/wait-for; \
    download "${NCBI_BLAST_URL}" /tmp/ncbi-blast.tar.gz; \
    tar -xzf /tmp/ncbi-blast.tar.gz -C /blast; \
    download "${EDIRECT_TAR_URL}" /tmp/edirect.tar.gz; \
    tar -xzf /tmp/edirect.tar.gz -C /blast; \
    download "${EDIRECT_XTRACT_LINUX_URL}" /tmp/xtract.Linux.gz; \
    gunzip -c /tmp/xtract.Linux.gz > /blast/edirect/xtract.Linux; \
    chmod 0755 /blast/edirect/xtract.Linux; \
    download "${EDIRECT_RCHIVE_LINUX_URL}" /tmp/rchive.Linux.gz; \
    gunzip -c /tmp/rchive.Linux.gz > /blast/edirect/rchive.Linux; \
    chmod 0755 /blast/edirect/rchive.Linux; \
    download "${EDIRECT_TRANSMUTE_LINUX_URL}" /tmp/transmute.Linux.gz; \
    gunzip -c /tmp/transmute.Linux.gz > /blast/edirect/transmute.Linux; \
    chmod 0755 /blast/edirect/transmute.Linux; \
    mkdir -p /blast/utilities; \
    download "${TRIMAL_URL}" /tmp/trimal-v1.4.1.zip; \
    unzip -q -d /blast/utilities /tmp/trimal-v1.4.1.zip; \
    make -C /blast/utilities/trimal-1.4.1/source; \
    download "${MVIEW_URL}" /tmp/mview-1.67.tar.gz; \
    tar -xzf /tmp/mview-1.67.tar.gz -C /blast/utilities; \
    mkdir -p /blast/utilities/mview; \
    cd /blast/utilities/mview-1.67 && perl install.pl /blast/utilities/mview -y; \
    download "${RPSBPROC_URL}" /tmp/RpsbProc-x64-linux.tar.gz; \
    tar -xzf /tmp/RpsbProc-x64-linux.tar.gz -C /blast/utilities; \
    mkdir -p /blast/utilities/RpsbProc-x64-linux/data; \
    download "${CDDID_TBL_URL}" /blast/utilities/RpsbProc-x64-linux/data/cddid.tbl.gz; \
    gunzip /blast/utilities/RpsbProc-x64-linux/data/cddid.tbl.gz; \
    download "${CDTRACK_URL}" /blast/utilities/RpsbProc-x64-linux/data/cdtrack.txt; \
    download "${FAMILY_SUPERFAMILY_LINKS_URL}" /blast/utilities/RpsbProc-x64-linux/data/family_superfamily_links; \
    download "${CDDANNOT_URL}" /blast/utilities/RpsbProc-x64-linux/data/cddannot.dat.gz; \
    gunzip /blast/utilities/RpsbProc-x64-linux/data/cddannot.dat.gz; \
    download "${CDDANNOT_GENERIC_URL}" /blast/utilities/RpsbProc-x64-linux/data/cddannot_generic.dat.gz; \
    gunzip /blast/utilities/RpsbProc-x64-linux/data/cddannot_generic.dat.gz; \
    download "${BITSCORE_SPECIFIC_URL}" /blast/utilities/RpsbProc-x64-linux/data/bitscore_specific.txt; \
    rm -rf /tmp/*

FROM ${UBUNTU_FOCAL_AMD64_IMAGE}

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG CATHI_DEPENDENCY_LOCK_REVISION=unrecorded
ARG MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh
ARG MICROMAMBA_URL=https://conda.anaconda.org/conda-forge/linux-64/micromamba-1.5.6-0.tar.bz2

LABEL org.opencontainers.image.cathi.dependency-lock-revision="${CATHI_DEPENDENCY_LOCK_REVISION}"
LABEL org.opencontainers.image.base.name="ubuntu:focal"

# netcat for wait-for / libdw1 for rpsbproc / wget and tar for runtime database setup.
RUN apt-get update -y && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y --no-install-recommends \
      bzip2 \
      ca-certificates \
      curl \
      fasttree \
      gzip \
      libdw1 \
      libnet-perl \
      libwww-perl \
      libxml-simple-perl \
      mafft \
      netcat \
      perl \
      tar \
      tzdata \
      wget && \
    rm -rf /var/lib/apt/lists/*

ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/${TZ} /etc/localtime && echo "${TZ}" > /etc/timezone

RUN groupadd --system --gid 10001 cathi && \
    useradd --system --uid 10001 --gid cathi --home-dir /blast --shell /usr/sbin/nologin cathi && \
    mkdir -p /blast/reciprocal_blast /blast/utilities

WORKDIR /blast

RUN set -eux; \
    curl -fsSL --retry 5 --proto '=https' --tlsv1.2 -o /tmp/miniconda.sh "${MINICONDA_URL}"; \
    bash /tmp/miniconda.sh -b -p /blast/miniconda3; \
    curl -fsSL --retry 5 --proto '=https' --tlsv1.2 -o /tmp/micromamba.tar.bz2 "${MICROMAMBA_URL}"; \
    tar -xjf /tmp/micromamba.tar.bz2 -C /tmp bin/micromamba; \
    install -m 0755 /tmp/bin/micromamba /usr/local/bin/micromamba; \
    rm -rf /tmp/miniconda.sh /tmp/micromamba.tar.bz2 /tmp/bin

ENV CATHI_CONDA_ENV=/blast/conda-envs/cathi-runtime
ENV PATH="/blast/edirect:/blast/ncbi-blast-2.11.0+/bin:/blast/miniconda3/bin:${CATHI_CONDA_ENV}/bin:/blast/utilities:/blast/utilities/trimal-1.4.1/source:/blast/utilities/mview:/blast/utilities/RpsbProc-x64-linux:${PATH}"

COPY --from=tool-builder /blast/ncbi-blast-2.11.0+ /blast/ncbi-blast-2.11.0+
COPY --from=tool-builder /blast/edirect /blast/edirect
COPY --from=tool-builder /blast/wait-for /blast/utilities/wait-for
COPY --from=tool-builder /blast/utilities/trimal-1.4.1 /blast/utilities/trimal-1.4.1
COPY --from=tool-builder /blast/utilities/mview /blast/utilities/mview
COPY --from=tool-builder /blast/utilities/mview-1.67 /blast/utilities/mview-1.67
COPY --from=tool-builder /blast/utilities/RpsbProc-x64-linux /blast/utilities/RpsbProc-x64-linux

RUN mkdir -p /blast/reciprocal_blast
COPY environment.yml environment-linux-64.lock.yml /blast/reciprocal_blast/
RUN micromamba create -y -p "${CATHI_CONDA_ENV}" -f /blast/reciprocal_blast/environment.yml --strict-channel-priority && \
    "${CATHI_CONDA_ENV}/bin/python" -m pip check && \
    "${CATHI_CONDA_ENV}/bin/python" -c "from Bio import Entrez; import pandas, matplotlib, seaborn, bokeh, yaml; print('Workflow Python imports ok')" && \
    "${CATHI_CONDA_ENV}/bin/snakemake" --version && \
    micromamba clean -afy && \
    conda clean -afy

COPY /celery_blast /blast/reciprocal_blast
COPY requirements.in requirements.txt requirements.lock.txt /blast/reciprocal_blast/
RUN mkdir -p /blast/reciprocal_blast/build_metadata

ENV BLASTDB=/blast/reciprocal_blast/media/databases
WORKDIR /blast/reciprocal_blast

RUN python3 -m pip install --no-cache-dir -r requirements.lock.txt && python3 -m pip check

RUN mkdir -p /blast/reciprocal_blast/media /blast/reciprocal_blast/tmp /blast/reciprocal_blast/build_metadata && \
    { \
      echo "CATHI_DEPENDENCY_LOCK_REVISION=${CATHI_DEPENDENCY_LOCK_REVISION}"; \
      echo "Dependency input sizes:"; \
      wc -c requirements.in \
            requirements.lock.txt \
            environment.yml \
            environment-linux-64.lock.yml; \
      echo ""; \
      echo "Python:"; \
      python3 --version; \
      echo ""; \
      echo "pip freeze:"; \
      python3 -m pip freeze; \
      echo ""; \
      echo "Conda base export:"; \
      conda list --export; \
      echo ""; \
      echo "Conda workflow environment export:"; \
      conda list -p "${CATHI_CONDA_ENV}" --export; \
      echo ""; \
      echo "Principal tool versions:"; \
      blastp -version | head -n 1; \
      snakemake --version; \
      { mafft --version 2>&1 || true; } | head -n 1; \
      { FastTree 2>&1 || true; } | head -n 1; \
    } > build_metadata/dependency_versions.txt && \
    chown -R cathi:cathi /blast

USER 10001:10001

CMD ["bash"]
