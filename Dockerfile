##
## Build
##
# NOTE: the platform is pinned to linux/amd64 because builds often happen on
# Apple Silicon, but the resulting image needs to run in x86 Linux batch/HPC
# environments.
ARG TARGETPLATFORM=linux/amd64
FROM --platform=$TARGETPLATFORM golang:1.24-alpine AS build

RUN apk add --no-cache gcc musl-dev

WORKDIR /app
COPY go.mod go.sum ./
RUN go mod download

# overwrite this at build time
# put this here so previous layers do not get invalidated
# https://stackoverflow.com/questions/60450479/using-arg-and-env-in-dockerfile
ARG Version=foo-docker-version

# Copy the root library package, tests, CLI package, and implementation packages.
COPY *.go ./
COPY cmd ./cmd
COPY fastq ./fastq
COPY fastqio ./fastqio
COPY sort ./sort

RUN go test ./...
RUN go build -trimpath -ldflags="-X 'squish.Version=$Version'" -o /squish ./cmd/squish

##
## Deploy
##

# NOTE: had issues with alpine on AWS Batch, so use Ubuntu for runtime.
FROM --platform=$TARGETPLATFORM ubuntu:22.04

# Runtime packages:
# - pigz: parallel gzip tooling for downstream shell workflows.
# - graphviz: provides dot, required by `go tool pprof -pdf`.
# - ca-certificates/procps: small operational utilities useful in batch systems.
# - golang-go: provides `go tool pprof` for PDF profile generation. Nextflow is
#   intentionally not installed here because Nextflow runs the container.
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        ca-certificates \
        graphviz \
        golang-go \
        pigz \
        procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /squish /usr/local/bin/squish
RUN which squish
RUN squish -h || true
