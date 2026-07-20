# syntax=docker/dockerfile:1.7
# AUTHOR Lukas Becker
# miniconda - django - snakemake - NCBI BLAST+ 2.11.0+ - edirect - Dockerfile

ARG UBUNTU_FOCAL_AMD64_IMAGE=ubuntu:focal@sha256:c664f8f86ed5a386b0a340d981b8f81714e21a8b9c73f658c4bea56aa179d54a

FROM ${UBUNTU_FOCAL_AMD64_IMAGE} AS tool-builder

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG NCBI_BLAST_URL=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz
ARG NCBI_BLAST_SHA256=93454cbdf5ba6f541745f31155efd9ba48bc6249fe3659b0aeaea4af62e62b58
ARG EDIRECT_TAR_URL=https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
ARG EDIRECT_TAR_SHA256=987db3477989d2b9c562ded4521e2554c1c673e4c2e3401e9596e08ece32778c
ARG EDIRECT_XTRACT_LINUX_URL=https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz
ARG EDIRECT_XTRACT_LINUX_SHA256=74492d23df44661e32d504cc6d893835c3afd872d6b78139c5e0d3aee1bfdf1d
ARG EDIRECT_RCHIVE_LINUX_URL=https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/rchive.Linux.gz
ARG EDIRECT_RCHIVE_LINUX_SHA256=727949bdce39ef265893300e2c7d9247410c21c980d32644e3aff452b86d0a01
ARG EDIRECT_TRANSMUTE_LINUX_URL=https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/transmute.Linux.gz
ARG EDIRECT_TRANSMUTE_LINUX_SHA256=9ef6654b7664a21f1d20cd248c874688286e9f239372641f1abf06d8c378fb80
ARG WAIT_FOR_URL=https://raw.githubusercontent.com/eficode/wait-for/v2.1.0/wait-for
ARG WAIT_FOR_SHA256=46b07615b92400d67ad73cab1984021c6eccd9b62aec47b26c93c95d0b6049f5
ARG TRIMAL_URL=https://github.com/inab/trimal/archive/refs/tags/v1.4.1.zip
ARG TRIMAL_SHA256=3747c3ba3bc930a97a4bcc095fe81238433ffda27a3d09dde6331a86f9b49a40
ARG MVIEW_URL=https://downloads.sourceforge.net/project/bio-mview/bio-mview/mview-1.67/mview-1.67.tar.gz
ARG MVIEW_SHA256=e5bac78960f8f6c091b2f7ea8a3c6075e9bea5a062391fd3e1e44fca14025e46
ARG RPSBPROC_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz
ARG RPSBPROC_SHA256=eebcffcb629704241f0d8787f7eaf52998d4edc3fc646a2c4544d60892a265b8
ARG CDDID_TBL_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
ARG CDDID_TBL_SHA256=852cb4c22f604b38f34c399c2fdad3d737fd15deaa5736c1a632249f0a6f6264
ARG CDTRACK_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt
ARG CDTRACK_SHA256=263046b8bc6bf4ce19db46913479b99f83914c2b7db0e6db6fb23139084cf828
ARG FAMILY_SUPERFAMILY_LINKS_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links
ARG FAMILY_SUPERFAMILY_LINKS_SHA256=780a29f2ae698131c719b9a8cade4bfe54390380b6a388240695299c09916e9b
ARG CDDANNOT_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot.dat.gz
ARG CDDANNOT_SHA256=9b028a6847322188dad39f0ba854de9c38ca2e057a16705ee47ca0fae42258bd
ARG CDDANNOT_GENERIC_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot_generic.dat.gz
ARG CDDANNOT_GENERIC_SHA256=97b15f39741293e155e0e95b49c70243cedb2b0b2ed59440c4f8f7287f134a13
ARG BITSCORE_SPECIFIC_URL=https://ftp.ncbi.nih.gov/pub/mmdb/cdd/bitscore_specific.txt
ARG BITSCORE_SPECIFIC_SHA256=bf24d622f8b127651ba1c4338965ff8fbe730281c70936f85cddeb13432a0939

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
    download_verify() { \
      local url="$1"; \
      local sha256="$2"; \
      local output="$3"; \
      curl -fsSL --retry 5 --proto '=https' --tlsv1.2 -o "${output}" "${url}"; \
      echo "${sha256}  ${output}" | sha256sum -c -; \
    }; \
    download_verify "${WAIT_FOR_URL}" "${WAIT_FOR_SHA256}" /tmp/wait-for; \
    install -m 0755 /tmp/wait-for /blast/wait-for; \
    download_verify "${NCBI_BLAST_URL}" "${NCBI_BLAST_SHA256}" /tmp/ncbi-blast.tar.gz; \
    tar -xzf /tmp/ncbi-blast.tar.gz -C /blast; \
    download_verify "${EDIRECT_TAR_URL}" "${EDIRECT_TAR_SHA256}" /tmp/edirect.tar.gz; \
    tar -xzf /tmp/edirect.tar.gz -C /blast; \
    download_verify "${EDIRECT_XTRACT_LINUX_URL}" "${EDIRECT_XTRACT_LINUX_SHA256}" /tmp/xtract.Linux.gz; \
    gunzip -c /tmp/xtract.Linux.gz > /blast/edirect/xtract.Linux; \
    chmod 0755 /blast/edirect/xtract.Linux; \
    download_verify "${EDIRECT_RCHIVE_LINUX_URL}" "${EDIRECT_RCHIVE_LINUX_SHA256}" /tmp/rchive.Linux.gz; \
    gunzip -c /tmp/rchive.Linux.gz > /blast/edirect/rchive.Linux; \
    chmod 0755 /blast/edirect/rchive.Linux; \
    download_verify "${EDIRECT_TRANSMUTE_LINUX_URL}" "${EDIRECT_TRANSMUTE_LINUX_SHA256}" /tmp/transmute.Linux.gz; \
    gunzip -c /tmp/transmute.Linux.gz > /blast/edirect/transmute.Linux; \
    chmod 0755 /blast/edirect/transmute.Linux; \
    mkdir -p /blast/utilities; \
    download_verify "${TRIMAL_URL}" "${TRIMAL_SHA256}" /tmp/trimal-v1.4.1.zip; \
    unzip -q -d /blast/utilities /tmp/trimal-v1.4.1.zip; \
    make -C /blast/utilities/trimal-1.4.1/source; \
    download_verify "${MVIEW_URL}" "${MVIEW_SHA256}" /tmp/mview-1.67.tar.gz; \
    tar -xzf /tmp/mview-1.67.tar.gz -C /blast/utilities; \
    mkdir -p /blast/utilities/mview; \
    cd /blast/utilities/mview-1.67 && perl install.pl /blast/utilities/mview -y; \
    download_verify "${RPSBPROC_URL}" "${RPSBPROC_SHA256}" /tmp/RpsbProc-x64-linux.tar.gz; \
    tar -xzf /tmp/RpsbProc-x64-linux.tar.gz -C /blast/utilities; \
    mkdir -p /blast/utilities/RpsbProc-x64-linux/data; \
    download_verify "${CDDID_TBL_URL}" "${CDDID_TBL_SHA256}" /blast/utilities/RpsbProc-x64-linux/data/cddid.tbl.gz; \
    gunzip /blast/utilities/RpsbProc-x64-linux/data/cddid.tbl.gz; \
    download_verify "${CDTRACK_URL}" "${CDTRACK_SHA256}" /blast/utilities/RpsbProc-x64-linux/data/cdtrack.txt; \
    download_verify "${FAMILY_SUPERFAMILY_LINKS_URL}" "${FAMILY_SUPERFAMILY_LINKS_SHA256}" /blast/utilities/RpsbProc-x64-linux/data/family_superfamily_links; \
    download_verify "${CDDANNOT_URL}" "${CDDANNOT_SHA256}" /blast/utilities/RpsbProc-x64-linux/data/cddannot.dat.gz; \
    gunzip /blast/utilities/RpsbProc-x64-linux/data/cddannot.dat.gz; \
    download_verify "${CDDANNOT_GENERIC_URL}" "${CDDANNOT_GENERIC_SHA256}" /blast/utilities/RpsbProc-x64-linux/data/cddannot_generic.dat.gz; \
    gunzip /blast/utilities/RpsbProc-x64-linux/data/cddannot_generic.dat.gz; \
    download_verify "${BITSCORE_SPECIFIC_URL}" "${BITSCORE_SPECIFIC_SHA256}" /blast/utilities/RpsbProc-x64-linux/data/bitscore_specific.txt; \
    rm -rf /tmp/*

FROM ${UBUNTU_FOCAL_AMD64_IMAGE}

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG CATHI_DEPENDENCY_LOCK_REVISION=unrecorded
ARG MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh
ARG MINICONDA_SHA256=1314b90489f154602fd794accfc90446111514a5a72fe1f71ab83e07de9504a7
ARG MICROMAMBA_URL=https://conda.anaconda.org/conda-forge/linux-64/micromamba-1.5.6-0.tar.bz2
ARG MICROMAMBA_SHA256=efe462c7ffcae8b338c7dd7b168ce8d48cfc60b48ab991d02a035c3b8d73633c

LABEL org.opencontainers.image.cathi.dependency-lock-revision="${CATHI_DEPENDENCY_LOCK_REVISION}"
LABEL org.opencontainers.image.base.name="ubuntu:focal@sha256:c664f8f86ed5a386b0a340d981b8f81714e21a8b9c73f658c4bea56aa179d54a"

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
    echo "${MINICONDA_SHA256}  /tmp/miniconda.sh" | sha256sum -c -; \
    bash /tmp/miniconda.sh -b -p /blast/miniconda3; \
    curl -fsSL --retry 5 --proto '=https' --tlsv1.2 -o /tmp/micromamba.tar.bz2 "${MICROMAMBA_URL}"; \
    echo "${MICROMAMBA_SHA256}  /tmp/micromamba.tar.bz2" | sha256sum -c -; \
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
COPY build/remote-artifacts.json /blast/reciprocal_blast/build_metadata/remote-artifacts.json
COPY scripts/generate_release_sbom.py /blast/reciprocal_blast/scripts/generate_release_sbom.py

ENV BLASTDB=/blast/reciprocal_blast/media/databases
WORKDIR /blast/reciprocal_blast

RUN python3 -m pip install --no-cache-dir -r requirements.lock.txt && python3 -m pip check

RUN mkdir -p /blast/reciprocal_blast/media /blast/reciprocal_blast/tmp /blast/reciprocal_blast/build_metadata && \
    python3 scripts/generate_release_sbom.py \
      --artifact-manifest build_metadata/remote-artifacts.json \
      --output build_metadata/scientific-tools.sbom.cdx.json && \
    { \
      echo "CATHI_DEPENDENCY_LOCK_REVISION=${CATHI_DEPENDENCY_LOCK_REVISION}"; \
      echo "Lock file SHA256:"; \
      sha256sum requirements.in \
                requirements.lock.txt \
                environment.yml \
                environment-linux-64.lock.yml \
                build_metadata/remote-artifacts.json; \
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
