#!/bin/bash

platform="$(uname | tr '[:upper:]' '[:lower:]')"
arch="$(arch | tr '[:upper:]' '[:lower:]')"

if [[ "$platform" == "darwin" ]]; then
    if [[ "$arch" == "arm64" ]]; then
        platform="darwin_arm64"
    else
        platform="darwin_i386"
    fi
else
    platform='linux'
fi

echo "detected platform: $platform"

if [[ ! -f "cmsearch" ]]; then
    ln -s infernal-1.0.2.binaries/cmsearch.$platform cmsearch
fi

# 'realpath' is in the package 'coreutils' https://anaconda.org/conda-forge/coreutils

a=$(realpath cmsearch)
b=$(basename $a)
c=`echo $b | sed 's#cmsearch.##'`

if [[ "$platform" != "$c" ]]; then
    rm cmsearch
    ln -s infernal-1.0.2.binaries/cmsearch.$platform cmsearch
fi
