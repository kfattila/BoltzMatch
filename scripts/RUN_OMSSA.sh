#!/bin/bash

OMSSA=~/omssa-2.1.9.linux/omssacl
BLAST=~/ncbi-blast-2.9.0+/bin/makeblastdb

mkdir OMSSA

$BLAST -in fasta/ipi.HUMAN.v3.87.target-protrev.fasta -out OMSSA/human -dbtype prot >/dev/null 
$BLAST -in fasta/ABRF_iPRG_2012.target-protrev.fasta -out OMSSA/iprg -dbtype prot >/dev/null 
$BLAST -in fasta/plasmo_Pfalciparum3D7_NCBI.target-protrev.fasta -out OMSSA/malaria -dbtype prot >/dev/null 
$BLAST -in fasta/yeast.target-protrev.fasta -out OMSSA/yeast -dbtype prot >/dev/null 

OUTPUT=./OMSSA/humvar_lr.csv
$OMSSA -d OMSSA/human -fm /home/data/mass_spec_data/Linfeng_012511_HapMap39_3.mgf -te 50 -teppm -to 1.0005079 -i 1,4 -oc $OUTPUT -nt 10 -mnm -mf 3 -mv 1,198,199 -e 0 -he 99999999 -hl 1 -zcc 1 -tem 0 -tom 0 -tez 1 -zh 9 -v 0 -hm 1

OUTPUT=./OMSSA/humvar_hr.csv
$OMSSA -d OMSSA/human -fm /home/data/mass_spec_data/Linfeng_012511_HapMap39_3.mgf -te 50 -teppm -to 0.05 -i 1,4 -oc $OUTPUT -nt 10 -mnm -mf 3 -mv 1,198,199 -e 0 -he 99999999 -hl 1 -zcc 1 -tem 0 -tom 0 -tez 1 -zh 9 -v 0 -hm 1

OUTPUT=./OMSSA/malaria_lr.csv
$OMSSA -d OMSSA/malaria -fm /home/data/mass_spec_data/v07548_UofF_malaria_TMT_10.mgf -te 50 -teppm -to 1.0005079 -i 1,4 -oc $OUTPUT -nt 10 -mnm -mf 3,198,199 -mv 1 -e 5 -he 99999999 -hl 1 -zcc 1 -tem 0 -tom 0 -tez 1 -zh 9 -v 0 -hm 1

OUTPUT=./OMSSA/malaria_hr.csv
$OMSSA -d OMSSA/malaria -fm /home/data/mass_spec_data/v07548_UofF_malaria_TMT_10.mgf -te 50 -teppm -to 0.05 -i 1,4 -oc $OUTPUT -nt 10 -mnm -mf 3,198,199 -mv 1 -e 5 -he 99999999 -hl 1 -zcc 1 -tem 0 -tom 0 -tez 1 -zh 9 -v 0 -hm 1

OUTPUT=./OMSSA/iprg_lr.csv
$OMSSA -d OMSSA/iprg -fm /home/data/mass_spec_data/iPRG_continuous_scan.mgf -te 10 -teppm -to 1.0005079 -i 1,4 -oc $OUTPUT -nt 10 -mnm -mf 3 -mv 1 -e 0 -he 99999999 -hl 1 -zcc 1 -tem 0 -tom 0 -tez 1 -zh 9 -v 0 -hm 1 

OUTPUT=./OMSSA/iprg_hr.csv
$OMSSA -d OMSSA/iprg -fm /home/data/mass_spec_data/iPRG_continuous_scan.mgf -te 10 -teppm -to 0.05 -i 1,4 -oc $OUTPUT -nt 10 -mnm -mf 3 -mv 1 -e 0 -he 99999999 -hl 1 -zcc 1 -tem 0 -tom 0 -tez 1 -zh 9 -v 0 -hm 1

OUTPUT=./OMSSA/yeast_lr.csv
$OMSSA -d OMSSA/yeast -fm /home/data/mass_spec_data/yeast-01_full.mgf -te 3  -to 1.0005079 -i 1,4 -oc $OUTPUT -nt 10 -mnm -mf 3 -e 0 -he 99999999 -hl 1 -zcc 1 -tem 0 -tom 0 -tez 1 -zh 9 -v 0

OUTPUT=./OMSSA/aurum_lr.csv
$OMSSA -d OMSSA/human -fm /home/data/mass_spec_data/AurumR151.mgf -te 2 -to 1.0005079 -i 1,4 -oc $OUTPUT -nt 10 -mnm -mf 3 -mv 1 -e 0 -he 99999999 -hl 1 -zcc 1 -tem 0 -tom 0 -tez 1 -zh 9 -v 0 -hm 1

OUTPUT=./OMSSA/hspp2a_lr.csv
$OMSSA -d OMSSA/human -fm /home/data/mass_spec_data/HSPP2A.mgf -te 50 -teppm -to 1.0005079 -i 1,4 -oc $OUTPUT -nt 10 -mnm -mf 3 -mv 1 -e 0 -he 99999999 -hl 1 -zcc 1 -tem 0 -tom 0 -tez 1 -zh 9 -v 0 -hm 1
