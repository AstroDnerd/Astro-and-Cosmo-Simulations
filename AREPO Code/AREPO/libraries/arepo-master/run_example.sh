#!/bin/bash            # this line only there to enable syntax highlighting in this file
#Bash file to run examples
eg_name=$1
export SYSTYPE="Ubuntu"
cp Template-Makefile.systype Makefile.systype
mkdir -p ./run/examples/$eg_name
cp -r ./examples/$eg_name/* ./run/examples/$eg_name/
python ./run/examples/$eg_name/create.py ./run/examples/$eg_name
make CONFIG=./run/examples/$eg_name/Config.sh \
BUILD_DIR=./run/examples/$eg_name/build \
EXEC=./run/examples/$eg_name/Arepo
cd ./run/examples/$eg_name/
mpiexec -np 1 ./Arepo param.txt
python ./check.py ./
