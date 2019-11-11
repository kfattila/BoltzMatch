#!/bin/bash

CRUX=/home/attila/crux-3.2.8.original/src/crux

#FIX_PARAM="--decoy-format none --enzyme custom-enzyme --custom-enzyme [Z]|{Z} --missed-cleavages 0 --overwrite T --peptide-list T"
FIX_PARAM="--decoy-format none --missed-cleavages 1 --overwrite T --peptide-list T"
FOLDER_POST=_IDX

FASTA=/home/data/Fasta/yeast.target-protrev.fasta 
TAXON=./CRUX/YEAST
PARAM="--max-mods 0"
$CRUX tide-index $FIX_PARAM $PARAM --output-dir $TAXON$FOLDER_POST $FASTA $TAXON$FOLDER_POST

FASTA=/home/data/Fasta/ABRF_iPRG_2012.target-protrev.fasta
TAXON=./CRUX/IPRG
PARAM="--max-mods 1 --mods-spec 1M+15.9949"
$CRUX tide-index $FIX_PARAM $PARAM --output-dir $TAXON$FOLDER_POST $FASTA $TAXON$FOLDER_POST

FASTA=/home/data/Fasta/plasmo_Pfalciparum3D7_NCBI.target-protrev.fasta
TAXON=./CRUX/MALARIA
PARAM="--max-mods 1 --enzyme lys-c --mods-spec K+229.16293,1M+15.9949 --nterm-peptide-mods-spec X+229.16293"
$CRUX tide-index $FIX_PARAM $PARAM --output-dir $TAXON$FOLDER_POST $FASTA $TAXON$FOLDER_POST

FASTA=/home/data/Fasta/ipi.HUMAN.v3.87.target-protrev.fasta
TAXON=./CRUX/HUMVAR
PARAM="--max-mods 2 --mods-spec 2K+229.16293,2M+15.9949 --nterm-peptide-mods-spec 1X+229.16293"
$CRUX tide-index $FIX_PARAM $PARAM --output-dir $TAXON$FOLDER_POST $FASTA $TAXON$FOLDER_POST

FASTA=/home/data/Fasta/ipi.HUMAN.v3.87.target-protrev.fasta
TAXON=./CRUX/HUMAN
PARAM="--max-mods 1 --mods-spec 1M+15.9949"
$CRUX tide-index $FIX_PARAM $PARAM --output-dir $TAXON$FOLDER_POST $FASTA $TAXON$FOLDER_POST


