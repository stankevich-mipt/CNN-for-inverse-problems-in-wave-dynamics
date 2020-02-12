#!/bin/bash

FILES=*

for file in $(find $PWD -name '*.h' -or -name '*.cpp' -or -name '*.inl');
do
    echo $file
    if [[ -f $file ]]; then
	sed -i -e '$a\' $file
    fi
done
