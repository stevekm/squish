SHELL:=/bin/bash
.ONESHELL:

# apply code formatting
clean:
	go mod tidy
	gofmt -l -w .

SRC:=main.go
BIN:=squish
FASTQGZ:=data/test1.fastq.gz
test-run:
	go run "$(SRC)" "$(FASTQGZ)" "output.gz"

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
	go build -ldflags="-X 'main.Version=$(GIT_TAG)'" -o ./$(BIN) ./$(SRC)
.PHONY:build

build-all:
	mkdir -p build ; \
	for os in darwin linux windows; do \
	for arch in amd64 arm64; do \
	output="build/$(BIN)-v$(GIT_TAG)-$$os-$$arch" ; \
	if [ "$${os}" == "windows" ]; then output="$${output}.exe"; fi ; \
	echo "building: $$output" ; \
	GOOS=$$os GOARCH=$$arch go build -ldflags="-X 'main.Version=$(GIT_TAG)'" -o "$${output}" $(SRC) ; \
	done ; \
	done

# input/SRR14575325.gz9.fastq.gz
FASTQIN:=input/mid.fastq.gz
FASTQOUT:=output/mid
test-run-all: $(BIN)
	set -x
	for i in alpha gc qual alpha-heap; do
	echo ">>> RUNNING: $$i"
	./$(BIN) -m "$$i" -orderFile "order.$$i.txt" -memProf "mem.$$i.prof" -cpuProf "cpu.$$i.prof" "$(FASTQIN)" "$(FASTQOUT)".$$i.fastq.gz
	$(MAKE) pdf PROF="mem.$$i.prof" PDF="memprofile.$$i.pdf"
	done

# docker build -t stevekm/squish:latest .
DOCKER_TAG:=stevekm/squish:$(GIT_TAG)
docker-build:
	docker build --build-arg "Version=$(GIT_TAG)" -t $(DOCKER_TAG) .

# docker push stevekm/squish:latest
docker-push:
	docker push $(DOCKER_TAG)

# docker-test-run:
# 	docker run --platform linux/amd64 --rm -ti -v ${PWD}:${PWD} --workdir ${PWD} $(DOCKER_TAG) $(BIN) $(FASTQ)