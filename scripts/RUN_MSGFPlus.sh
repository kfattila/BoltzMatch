#!/bin/bash

MSGFPlus=../../msgfplus2019/MSGFPlus.jar
FIXED_PARAM="-maxMissedCleavages 1 -tda 0 -thread 1 -addFeatures 1 -ignoreMetCleavage 1 -n 1 -ti 0,0 -minLength 6 -maxLength 50 -minCharge 1 -maxCharge 9 -ntt 2"

mkdir MSGFPlus

# HumVar
# MSFILE=./ms_data/Linfeng_012511_HapMap39_3.mzML
# FASTA=./fasta/ipi.HUMAN.v3.87.target-protrev.fasta
# PARAM="-t 50ppm -e 1 "

# OUTFILE=./MSGFPlus/humvar_lr.mzid
# java -Xmx64G -jar $MSGFPlus -s $MSFILE -d $FASTA $FIXED_PARAM $PARAM -mod ./oxi_tmt_var_mods.txt -o $OUTFILE
# java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -showDecoy 1 -showQValue 1 -i $OUTFILE

# OUTFILE=./MSGFPlus/humvar_hr.mzid
# java -Xmx64G -jar $MSGFPlus -m 3 -inst 1 -s $MSFILE -d $FASTA $FIXED_PARAM $PARAM -mod ./oxi_tmt_var_mods.txt -o $OUTFILE
# java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -showDecoy 1 -showQValue 1 -i $OUTFILE

# IPRG
MSFILE=./ms_data/iPRG_continuous_scan.mzML
FASTA=./fasta/ABRF_iPRG_2012.target-protrev.fasta
PARAM="-t 10ppm -e 1 "

OUTFILE=./MSGFPlus/iprg_lr.mzid
java -Xmx64G -jar $MSGFPlus -s $MSFILE -d $FASTA $FIXED_PARAM $PARAM -mod ./oxi_mods.txt -o $OUTFILE
java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -showDecoy 1 -showQValue 1 -i $OUTFILE

OUTFILE=./MSGFPlus/iprg_hr.mzid
java -Xmx64G -jar $MSGFPlus -m 0 -inst 2 -s $MSFILE -d $FASTA $FIXED_PARAM $PARAM -mod ./oxi_mods.txt -o $OUTFILE
java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -showDecoy 1 -showQValue 1 -i $OUTFILE

# Malaria
MSFILE=./ms_data/v07548_UofF_malaria_TMT_10.mzML
FASTA=./fasta/plasmo_Pfalciparum3D7_NCBI.target-protrev.fasta
PARAM="-t 50ppm -e 3"
OUTFILE=./MSGFPlus/malaria_lr.mzid
java -Xmx64G -jar $MSGFPlus -inst 0 -s $MSFILE -d $FASTA $FIXED_PARAM $PARAM -mod ./oxi_tmt_mods.txt -o $OUTFILE
java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -showDecoy 1 -showQValue 1 -i $OUTFILE

OUTFILE=./MSGFPlus/malaria_hr.mzid
java -Xmx64G -jar $MSGFPlus -m 0 -inst 1 -s $MSFILE -d $FASTA $FIXED_PARAM $PARAM -mod ./oxi_tmt_mods.txt -o $OUTFILE
java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -showDecoy 1 -showQValue 1 -i $OUTFILE

# Aurum
MSFILE=./ms_data/AurumR151.mzML
FASTA=./fasta/ipi.HUMAN.v3.87.target-protrev.fasta
PARAM="-t 2Da -e 1"

OUTFILE=./MSGFPlus/aurum_lr.mzid
java -Xmx64G -jar $MSGFPlus -s $MSFILE -d $FASTA $FIXED_PARAM $PARAM -mod ./oxi_mods.txt -o $OUTFILE
java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -showDecoy 1 -showQValue 1 -i $OUTFILE

# HSPP2A
MSFILE=./ms_data/HSPP2A.mzML
FASTA=./fasta/ipi.HUMAN.v3.87.target-protrev.fasta
PARAM="-t 50ppm -e 1"

OUTFILE=./MSGFPlus/hspp2a_lr.mzid
java -Xmx64G -jar $MSGFPlus -s $MSFILE -d $FASTA $FIXED_PARAM $PARAM -mod ./oxi_mods.txt -o $OUTFILE
java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -showDecoy 1 -showQValue 1 -i $OUTFILE

# Yeast
MSFILE=./ms_data/yeast-01_full.mzML
FASTA=./fasta/yeast.target-protrev.fasta 
PARAM="-t 3Da -e 1"

OUTFILE=./MSGFPlus/yeast_lr.mzid
java -Xmx64G -jar $MSGFPlus -s $MSFILE -d $FASTA $FIXED_PARAM $PARAM -o $OUTFILE
java -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -showDecoy 1 -showQValue 1 -i $OUTFILE

