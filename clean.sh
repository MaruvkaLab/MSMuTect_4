#!/bin/bash

# be very careful about running this script

rm pyinstaller_build/dist/msmutect
find . -iname '*.so' -delete
find . -iname '*.c' -delete
bash reverse_rename.sh