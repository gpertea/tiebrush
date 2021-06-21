#!/usr/bin/env bash
## 
#directory name must match $BWLIB in Makefile

if [[ ! -d bigwig ]]; then
 git clone https://github.com/alevar/libBigWig.git
 mv libBigWig bigwig
fi

cd bigwig

if [[ "$1" == "clean" ]]; then
  make clean
  exit
fi

if [[ ! -f bigwig/libBigWig.a ]]; then
  make -j 2 lib-static || exit 1
fi


