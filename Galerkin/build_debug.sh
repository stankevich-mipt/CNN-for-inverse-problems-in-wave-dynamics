#!/bin/sh

rm -r debug
mkdir debug
cd debug

cmake -DCMAKE_BUILD_TYPE=Debug ../
make

cd ../.bin
echo "Done"