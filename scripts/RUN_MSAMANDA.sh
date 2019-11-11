#!/bin/bash

MSAMANDA=~/MSAmanda/MSAmanda.exe

mkdir MSAmanda

mono $MSAMANDA -s ./ms_data/AurumR151.mgf -d ./fasta/ipi.HUMAN.v3.87.target-protrev.fasta -e ./scripts/msamanda_settings_aurum_lr.xml -f 1 -o msamanda/msamanda_out_aurum_lr.csv
mono $MSAMANDA -s ./ms_data/HSPP2A.mgf -d ./fasta/ipi.HUMAN.v3.87.target-protrev.fasta -e ./scripts/msamanda_settings_hspp2a_lr.xml -f 1 -o msamanda/msamanda_out_hspp2a_lr.csv
mono $MSAMANDA -s ./ms_data/yeast-01_full.mgf -d ./fasta/yeast.target-protrev.fasta -e ./scripts/msamanda_settings_yeast_lr.xml -f 1 -o msamanda/msamanda_out_yeast_lr.csv
mono $MSAMANDA -s ./ms_data/iPRG_continuous_scan.mgf -d ./fasta/ABRF_iPRG_2012.target-protrev.fasta -e ./scripts/msamanda_settings_iprg_lr.xml -f 1 -o msamanda/msamanda_out_iprg_lr.csv
mono $MSAMANDA -s ./ms_data/iPRG_continuous_scan.mgf -d ./fasta/ABRF_iPRG_2012.target-protrev.fasta -e ./scripts/msamanda_settings_iprg_hr.xml -f 1 -o msamanda/msamanda_out_iprg_hr.csv
mono $MSAMANDA -s ./ms_data/Linfeng_012511_HapMap39_3.mgf -d ./fasta/ipi.HUMAN.v3.87.target-protrev.fasta -e ./scripts/msamanda_settings_humvar_lr.xml -f 1 -o msamanda/msamanda_out_humvar_lr.csv
mono $MSAMANDA -s ./ms_data/Linfeng_012511_HapMap39_3.mgf -d ./fasta/ipi.HUMAN.v3.87.target-protrev.fasta -e ./scripts/msamanda_settings_humvar_hr.xml -f 1 -o msamanda/msamanda_out_humvar_hr.csv
mono $MSAMANDA -s ./ms_data/v07548_UofF_malaria_TMT_10.mgf -d ./fasta/plasmo_Pfalciparum3D7_NCBI.target-protrev.fasta -e ./scripts/msamanda_settings_malaria_hr.xml -f 1 -o msamanda/msamanda_out_malaria_hr.csv
mono $MSAMANDA -s ./ms_data/v07548_UofF_malaria_TMT_10.mgf -d ./fasta/plasmo_Pfalciparum3D7_NCBI.target-protrev.fasta -e ./scripts/msamanda_settings_malaria_lr.xml -f 1 -o msamanda/msamanda_out_malaria_lr.csv
