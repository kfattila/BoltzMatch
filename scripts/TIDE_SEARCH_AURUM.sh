#!/bin/bash

# Update path if necessary

CRUX=/home/attila/crux-toolkit-tailor/src/crux
PATH=./CRUX/

DATA_FILE_ORIGINAL=ms_data/AurumR151.ms2
DATA_FILE_BOLTZMATCH_LR=ms_data/AurumR151_boltzmatch_lowres.ms2 


FIX_PARAM="--mz-bin-offset 0.4 --precursor-window 2 --precursor-window-type mass --top-match 1 --concat T --overwrite T --num-threads 1 --use-neutral-loss-peaks F --use-flanking-peaks F --min-peaks 10 --max-precursor-charge 9"

INDEX=$PATH/HUMAN_IDX
TAXON=AURUM

## RUN Tide-search to obtain baseline results
FOLDER_POST=_BASELINE_XCORR_LR
PARAM="--exact-p-value F --mz-bin-width 1.0005079 --use-tailor-calibration T"
OUTPUT_DIR=$PATH$TAXON$FOLDER_POST
$CRUX tide-search $FIX_PARAM $PARAM --output-dir $OUTPUT_DIR $DATA_FILE_ORIGINAL $INDEX
$CRUX assign-confidence --score "xcorr score" --overwrite T --output-dir $OUTPUT_DIR $OUTPUT_DIR/tide-search.txt
$CRUX assign-confidence --score "tailor score" --overwrite T --output-dir $OUTPUT_DIR $OUTPUT_DIR/tide-search.txt

FOLDER_POST=_BASELINE_XPV
PARAM="--exact-p-value T --mz-bin-width 1.0005079 "
OUTPUT_DIR=$PATH$TAXON$FOLDER_POST
$CRUX tide-search $FIX_PARAM $PARAM --output-dir $OUTPUT_DIR $DATA_FILE_ORIGINAL $INDEX
$CRUX assign-confidence --overwrite T --output-dir $OUTPUT_DIR $OUTPUT_DIR/tide-search.txt

## RUN tide-search on spectrum data obtained with BoltzMatch, use XPV or Tailor for score calibration

FOLDER_POST=_BOLZMATCH_XPV
PARAM="--exact-p-value T --mz-bin-width 1.0005079 --cross-corr-penalty F --skip-preprocessing T "
OUTPUT_DIR=$PATH$TAXON$FOLDER_POST
$CRUX tide-search $FIX_PARAM $PARAM --output-dir $OUTPUT_DIR $DATA_FILE_BOLTZMATCH_LR $INDEX
$CRUX assign-confidence --overwrite T --output-dir $OUTPUT_DIR $OUTPUT_DIR/tide-search.txt

FOLDER_POST=_BOLZMATCH_TAILOR_LR
PARAM="--exact-p-value F --mz-bin-width 1.0005079 --use-tailor-calibration T --cross-corr-penalty F --skip-preprocessing T "
OUTPUT_DIR=$PATH$TAXON$FOLDER_POST
$CRUX tide-search $FIX_PARAM $PARAM --output-dir $OUTPUT_DIR $DATA_FILE_BOLTZMATCH_LR $INDEX
$CRUX assign-confidence --score "xcorr score" --overwrite T --output-dir $OUTPUT_DIR $OUTPUT_DIR/tide-search.txt
$CRUX assign-confidence --score "tailor score" --overwrite T --output-dir $OUTPUT_DIR $OUTPUT_DIR/tide-search.txt
