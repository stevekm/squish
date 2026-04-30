SHELL:=/bin/bash
.ONESHELL:

# apply code formatting
clean:
	go mod tidy
	gofmt -l -w .

test:
	tmp_cache="$$(mktemp -d)" ; \
	trap 'rm -rf "$$tmp_cache"' EXIT ; \
	GOCACHE="$$tmp_cache" go test -count=1 ./...
.PHONY:test

CMD:=./cmd/squish
BIN:=squish
FASTQGZ:=data/test1.fastq.gz
test-run:
	tmp_cache="$$(mktemp -d)" ; \
	trap 'rm -rf "$$tmp_cache"' EXIT ; \
	GOCACHE="$$tmp_cache" go run "$(CMD)" "$(FASTQGZ)" "output.gz"

# brew install graphviz
PROF:=mem.prof
PDF:=memprofile.pdf
pdf:
	go tool pprof -pdf -output $(PDF) $(PROF)
.PHONY:pdf
# && open -a Preview $(PDF)

# NOTE: you can just ignore this error message;
# fatal: No names found, cannot describe anything.
GIT_TAG:=$(shell git describe --tags)

build:
	tmp_cache="$$(mktemp -d)" ; \
	trap 'rm -rf "$$tmp_cache"' EXIT ; \
	GOCACHE="$$tmp_cache" go build -trimpath -ldflags="-X 'squish.Version=$(GIT_TAG)'" -o ./$(BIN) $(CMD)
.PHONY:build

build-all:
	mkdir -p build ; \
	tmp_cache="$$(mktemp -d)" ; \
	trap 'rm -rf "$$tmp_cache"' EXIT ; \
	for os in darwin linux windows; do \
	for arch in amd64 arm64; do \
	output="build/$(BIN)-v$(GIT_TAG)-$$os-$$arch" ; \
	if [ "$${os}" == "windows" ]; then output="$${output}.exe"; fi ; \
	echo "building: $$output" ; \
	GOOS=$$os GOARCH=$$arch GOCACHE="$$tmp_cache" go build -trimpath -ldflags="-X 'squish.Version=$(GIT_TAG)'" -o "$${output}" $(CMD) ; \
	done ; \
	done


#
# USAGE:
# make test-run-all FASTQIN=data/test1.fastq.gz
#
# input/SRR14575325.gz9.fastq.gz
FASTQIN:=data/SRR6357076_1.fastq.gz
FASTQOUT:=squish.SRR6357076_1
# FASTQIN:=data/SRX1603629_T1_1.fastq.gz
# FASTQOUT:=squish.SRX1603629_T1_1
SAMPLESHEET:=samples.csv
OUTDIR:=output
ENGINE:=memory
BUCKET:=auto
BUCKETS:=512
CLUMPK:=16
$(OUTDIR):
	mkdir -p "$(OUTDIR)"
test-run-all: build $(BIN) $(OUTDIR)
	set -x ; \
	for i in alpha gc qual clump; do \
	echo ">>> RUNNING: method=$$i engine=$(ENGINE) bucket=$(BUCKET) buckets=$(BUCKETS) clumpK=$(CLUMPK)" ; \
	./$(BIN) -outdir "$(OUTDIR)" -engine "$(ENGINE)" -bucket "$(BUCKET)" -buckets "$(BUCKETS)" -clumpK "$(CLUMPK)" -m "$$i" -orderFile "order.$$i.txt" -reportFile "report.$$i.json" -memProf "mem.$$i.prof" -cpuProf "cpu.$$i.prof" "$(FASTQIN)" "$(FASTQOUT)".$$i.fastq.gz ; \
	$(MAKE) pdf PROF="$(OUTDIR)/profile.$$i/mem.$$i.prof" PDF="$(OUTDIR)/profile.$$i/memprofile.$$i.pdf" ; \
	done

test-run-all-external:
	$(MAKE) test-run-all ENGINE=external

nextflow-test-run-all:
	nextflow run main.nf --engine memory --samplesheet "$(SAMPLESHEET)"

nextflow-test-run-all-external:
	nextflow run main.nf --engine external --samplesheet "$(SAMPLESHEET)"

# docker build -t stevekm/squish:latest .
DOCKER_TAG:=stevekm/squish:$(GIT_TAG)
docker-build:
	docker build --build-arg "Version=$(GIT_TAG)" -t $(DOCKER_TAG) .

# docker push stevekm/squish:latest
docker-push:
	docker push $(DOCKER_TAG)

# docker-test-run:
# 	docker run --platform linux/amd64 --rm -ti -v ${PWD}:${PWD} --workdir ${PWD} $(DOCKER_TAG) $(BIN) $(FASTQ)
