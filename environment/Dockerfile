# hash:sha256:c5801b906004e7bdf352793662d77a7d10e810991daf11b34cd474aeb5b851fa
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/omix-r-statistics:r4.4.3-bioconductor3.20-v1

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
