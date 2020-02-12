#!/bin/sh

rm -r release
mkdir release
cd release

cmake -DCMAKE_BUILD_TYPE=Release ..
make

cd ../.bin
echo "Done"