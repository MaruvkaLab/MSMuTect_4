#!/bin/bash

MSMUTECT_PATH=$(readlink -f "$0")
MSMUTECT_DIR=$(dirname $MSMUTECT_PATH)
export PYTHONPATH=$MSMUTECT_DIR
python3 $MSMUTECT_DIR/src/Entry/main.py "$@"
