## use BASH
SHELL:=/bin/bash

## name of binary
BINARY?=gobcly

## used for compiling with debug enabled / disabled
DEBUGFLAG?=TRUE

## hyperscan can be installed using homebrew on macOS or Linux 
## eg: "/usr/local/Cellar/hyperscan/5.2.1/include/hs"
## eg: "/home/linuxbrew/.linuxbrew/Cellar/hyperscan/5.4.2/include/hs"
#HSLIB?=$(shell brew info hyperscan | grep "Cellar/hyperscan" | cut -f 1 -d " ")/include/hs
#HSLIB:="/home/linuxbrew/.linuxbrew/Cellar/hyperscan/5.4.2/include/lib"

## OS
GOOS?=$(shell uname -s | perl -ne "chomp; print lc($$_);")
## ARCH
GOARCH?=$(shell uname -m | perl -ne "chomp; s/^x86_64$$/amd64/; print;")

## arch-specific
BINARY:=$(if $(GOOS) && $(GOARCH),$(BINARY).$(GOOS)_$(GOARCH),"")

## default
LDFLAGS:="-X main.Binary=${BINARY} -X main.DebugFlag=${DEBUGFLAG}"

## DIST
DIST:=./dist
DISTBIN:=${DIST}/gobcly

## will write this to a file 'VERSION' in the repo, used by 'govvv'
## expecting that we are using git flow and tagging releases
VERSION=`(git flow release list 2>/dev/null | cut -f 2 -d " ") | (git describe --tags)`

# Builds the project
build:
## this shouldn't be necessary but check anyway
	command -v go >/dev/null 2>&1 || ( echo "Make requires \'go\' in your path" 1>&2 ; exit 1 )
## requires govvv is installed
#	command -v govvv >/dev/null 2>&1 || ( go get -u github.com/ahmetb/govvv )
#	command -v govvv >/dev/null 2>&1 || ( echo "Tried to install \'govvv\' but doesn\'t appear to be in your path!" 1>&2 ; exit 1 )
## write the version file
	echo ${VERSION} >VERSION
## do the build with govvv
	env GOOS=${GOOS} GOARCH=${GOARCH} go build -o ${BINARY} -ldflags ${LDFLAGS} .

# Installs our project: copies binaries
install: export BINARY=gobcly
install:
	#env GOOS=${GOOS} GOARCH=${GOARCH} govvv install -ldflags ${LDFLAGS} .
	env GOOS=${GOOS} GOARCH=${GOARCH} go install -ldflags ${LDFLAGS} .

# Cleans our project: deletes binaries
clean:
	if [ -f ${BINARY} ] ; then rm ${BINARY} ; fi
	if [ -f ${DISTBIN} ] ; then rm ${DISTBIN} && rmdir ${DIST} ; fi

dist: export CGO_CFLAGS="-I/usr/local/lib"
dist: export GOOS=linux
dist: export GOARCH=amd64
dist: export BINARY=gobcly.linux_amd64
dist: build
	if [ -x upx ] ; then upx ${BINARY} ; fi
	mkdir -p ${DIST}
	mv ${BINARY} ${DISTBIN}

.PHONY: clean install dist
