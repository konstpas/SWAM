#!/bin/bash

old_format=F90
new_format=f90

for f in *.${old_format}; do

    mv -- "$f" "${f%.$old_format}.$new_format"

done
