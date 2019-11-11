#!/bin/bash
# Define resolution
resolution='lowres'
# Perform a smoke test to see all components and modules are installed well.
python3 ./rbm_train_main.py -d ups1 -r $resolution # ups1
