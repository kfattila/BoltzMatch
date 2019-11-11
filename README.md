BoltzMatch is machine-learning-based peptide-spectrum-matching scoring method for spectrum identification in tandem mass spectrometry, shotgun proteomics.

Authors:
Anastasia Voronkova, Pavel Sulimov, and Attila Kertesz-Farkas
National Research University Higher School of Economics (HSE), Moscow, Russia


## INSTALLATION
Download the code. 


## PREREQUISITES
BoltzMatch was tested with 
Python (3.6 )
Pytorch (1.1.0)
Matplotlib (2.0)
pyteomics (3.4)

and other standard python toolboxes such as:
argparse, time, numpy, seaborn, pandas, etc.



## SMOKE TEST
To see whether all components are installed properly, run:

$./RUN_TEST.sh

This script runs BoltzMatch on a small data. It should run not longer than 1-2 minutes.
If this terminates without error, you are all set.



## REPRODUCE OUR EXPERIMENTS
1. 
A modified version of Tide-search/CRUX is needed, which can (a) load observed spectrum data files containing peak intensities of negative values, 
and (b) the application of the cross-correlation penalty can be turned off. Such a modified version of CRUX is available at our repository:
https://github.com/kfattila/crux-toolkit

Download and install this copy of CRUX. Installation guide can be found at http://crux.ms/
 
2.
Unzip the data in folders ms_data and fasta, specify the paths to the datbase search programs in the RUN_ALL.sh script.

3. 
and run the RUN_ALL script:

$./RUN_ALL.sh

This script will run BoltzMatch on all of our benchmark data sets and lanches the gold standard database search programs (MS Amanda, OMMSA, X!Tandem, MS-GF+, etc).
Some programs require .mfg file format and you will need to use a program to convert between different file formats such as msconvert from ProteoWizard ( http://proteowizard.sourceforge.net/ ).


## TRAINING OF BOLTZMATCH 
To train the BoltzMatch on any data sets you can use the rbm_trian_main.py python program. To see the command line options type:
$ python rbm_trian_main.py -h

First, you need to create a dictionary containing the parameters corresponding to your datasets and add these parameters to the source/parameters.py file. In this file, go to the dictionary, called spectrum, and create a new key-value pair for your own data sets. The examples provided there can help.
BoltzMatch works with only mzML file format. Use ProteoWizard to convert between different file mass spec formats.

You may modifiy other parameters in this data files related to training etc., parameter explanations and examples are given in the parameter.py file.

Run BoltzMatch:

$ python rbm_trian_main.py -d NAME_OF_YOUR_DATASET -r RES

where RES can be 'lowres' or 'highres'. This controls the resolution (bin_width of the spectrum discretization step) of the MS-2 data. 

Note that, the parameters of BotlzMatch in 'highres' setting do not fit in standard GPUs. This usually needs to run on CPU.
For bin_width=1.00, BoltzMatch requires adound 40Mb of RAM.
For bin_width=0.05, BoltzMatch requires around 16Gb of RAM.
For bin_width=0.02, BoltzMatch requires adound 90Gb of RAM.


## Contact:

It will be added later.






