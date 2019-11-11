"""
This file contains source code of BoltzMatch.
Copyright 2019 by Anastasia Voronkova, Pavel Sulimov, and Attila Kertesz-Farkas.

All information contained herein is, and remains
the property of authors. The intellectual and technical concepts contained
herein are proprietary to authors.
Distributed under the Apache License, Version 2.0.
"""
import torch

# Specify global biological parameters
dbsearch = {
    'highres' : {'bin_width':0.05, 'bin_offset':0.4},  # Specify the bin_width for the spectrum discretization step
    'lowres' : {'bin_width':1.0005079, 'bin_offset':0.4},
    'max_peak' : 2000,                    # The maximum m/z of peaks considered. All peaks at location higher than this value are discared
    'remove_precursor_peak' : True,       # Remove peak around the precursor mass or not. 
    'remove_precursor_tolerance' : 1.5,   # This is a tolerance window for the previous option.
    'max_theo_pept_peak_charge' : 1,      # This controls the maximal charge states of the theoretical peaks. 
    'min_pept_len' : 7,                   # Minimal peptide length considered.
    'decoy_format' : 0,                   # To generate decoy peptides: 0 for peptide reverse, 1 for peptide shuffle
    'max_intensity' : 1.0,                # Peak intensity scaling value. This is the maximal peak intensity. 
}

# Specify types of variables
types = {
    'float_type': torch.float32,
    'dtype': torch.FloatTensor, # use torch.cuda.FloatTensor in case GPU is available and torch.FloatTensor otherwise
    # 'torch_device': torch.device("cpu") # BoltzMatch training on CPU 
   'torch_device': torch.device("cuda:0") # BoltzMatch training on GPU 
}

# Define training parameters
training = {
    'epochs' :  5,
    'batch_size':  128,
    'learning_rate' : 0.001,
    'cand_from_search_qval' : 0.005, # q-value threshold for candidate peptides for generating spectrum-peptide pairs for BoltzMatch
    'topN_exp' : 100 # value for N most intense peaks to keep from raw spectrum
}

# Define datasets and their specific biological parameters
spectrum = {
    # 'data_set_name' : {
        # 'ms_data_file': path_to_the_mzML_file,
        # 'fasta_file': path_to_the_fasta_file,
        # 'max_mods' : <INT>,          # The number of the modifications allowed by peptides
        # 'tolarence_type' : <STRING>, # Possible values are: "MASS" or "PPM"
        # 'tolarence_window' : <INT>,  # The tolerance window in the units specified above
        # 'missed_cleavages' : <INT>,  # Number of missed cleavages specified for peptide generation
        # 'enzyme': <STRING>,          # the name of the digestion enzyme: possible values: 'trypsin', 'trypsin/p', 'no-digestion' (use this for predigested peptides), use enzymes from parser.expasy_rules[self.enzyme]
        # 'static_mods': <LIST of STRING>, #Static modifications, e.g. {'C+57.02146'}
        # 'modifications' : <LIST OF DICTIONARIES>, # Variable modifications in the format of pyteomics. See: https://pyteomics.readthedocs.io/en/latest/  and https://pyteomics.readthedocs.io/en/latest/api/parser.html#pyteomics.parser.isoforms 
        # 'scan_id': <CODE>,           # code to read scan id from *.mzML file, e.g.: "int(spectrum['id'][5:])"
        # 'mz_array': <CODE>,          # code to read m/z array from *.mzML file, e.g.: "np.float_(spectrum['m/z array'][:])"
        # 'intensity_array': <CODE>,   # code to read intensity array from *.mzML file, e.g.: "np.float_(spectrum['m/z array'][:])",
        # 'charge': <CODE>,            # code to read charge from *.mzML file, e.g. "int(spectrum['ms2 file charge state'][0:1])", 
        # 'precursor_mass': <CODE>,    # code to read precursor mass from *.mzML file, e.g. "float(spectrum['ms2 file charge state'][2:])",
    # },
    'humvar' : {
        'ms_data_file': './ms_data/Linfeng_012511_HapMap39_3.mzML',
        'fasta_file': './fasta/ipi.HUMAN.v3.87.fasta',
        'max_mods' : 2,
        'tolarence_type' : "PPM",
        'tolarence_window' : 50,
        'missed_cleavages' : 0,
        'enzyme':'trypsin',
        'static_mods': {'C+57.02146'},
        'modifications' : {'[15.9949]':['M'],'[229.16293]':['K'],'[229.16293]-':True},
        'scan_id': "int(spectrum['id'][5:])", # code to read scan id from *.mzML file
        'mz_array': "np.float_(spectrum['m/z array'][:])", # code to read m/z array from *.mzML file
        'intensity_array': "np.float_(spectrum['intensity array'][:])", # code to read intensity array from *.mzML file
        'charge': "int(spectrum['ms2 file charge state'][0:1])", # code to read charge from *.mzML file
        'precursor_mass': "float(spectrum['ms2 file charge state'][2:])", # code to read precursor mass from *.mzML file
    },
    'ups1' : {
        'ms_data_file' : './ms_data/UPS1.recalre.mzML',
        'fasta_file' : './fasta/ups1.fasta',
        'tolarence_type' : "PPM",
        'tolarence_window' : 100,
        'max_mods' : 1,
        'enzyme':'trypsin',
        'missed_cleavages' : 1,
        'static_mods': {'C+57.02146'},
        'modifications' : {'[15.9949]':['M']},
        'scan_id': "spectrum_id", # code to read scan id from *.mzML file
        'mz_array': "np.float_(spectrum['m/z array'][:])", # code to read m/z array from *.mzML file
        'intensity_array': "np.float_(spectrum['intensity array'][:])", # code to read intensity array from *.mzML file
        'charge': "int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])", # code to read charge from *.mzML file
        'precursor_mass': "float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])", # code to read precursor mass from *.mzML file
    },
    'malaria' : {
        'ms_data_file' : './ms_data/v07548_UofF_malaria_TMT_10.mzML',
        'fasta_file' : './fasta/plasmo_Pfalciparum3D7_NCBI.fasta',
        'max_mods' : 1,
        'tolarence_type' : "PPM",
        'tolarence_window' : 20,
        'missed_cleavages' : 2,
        'enzyme' : 'lysc',
        'static_mods': {'C+57.02146', 'K+229.16293', 'Nt+229.16293'},
        'modifications' : {'[15.9949]':['M']},
        'scan_id': "int(spectrum['id'][5:])", # code to read scan id from *.mzML file
        'mz_array': "np.float_(spectrum['m/z array'][:])", # code to read m/z array from *.mzML file
        'intensity_array': "np.float_(spectrum['intensity array'][:])", # code to read intensity array from *.mzML file
        'charge': "int(spectrum['ms2 file charge state'][0:1])", # code to read charge from *.mzML file
        'precursor_mass': "float(spectrum['ms2 file charge state'][2:])", # code to read precursor mass from *.mzML file
     },
    'aurum' : {
        'ms_data_file' : './ms_data/AurumR151.mzML',
        'fasta_file' : './fasta/ipi.HUMAN.v3.87.fasta',
        'max_mods' : 1,
        'tolarence_type' : "MASS",
        'tolarence_window' : 2,
        'missed_cleavages' : 0,
        'enzyme':'trypsin',
        'static_mods': {'C+57.02146'},
        'modifications' : {'[15.9949]':['M']},
        'scan_id': "spectrum_id", # code to read scan id from *.mzML file
        'mz_array': "np.float_(spectrum['m/z array'][:])", # code to read m/z array from *.mzML file
        'intensity_array': "np.float_(spectrum['intensity array'][:])", # code to read intensity array from *.mzML file
        'charge': "int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])", # code to read charge from *.mzML file
        'precursor_mass': "float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])", # code to read precursor mass from *.mzML file
    },
    'hspp2a' : {
        'ms_data_file' : './ms_data/HSPP2A.mzML',
        'fasta_file' : './fasta/ipi.HUMAN.v3.87.fasta',    
        'max_mods' : 1,
        'tolarence_type' : "PPM",
        'tolarence_window' : 20,
        'missed_cleavages' : 1,
        'enzyme':'trypsin',
        'static_mods': {'C+57.02146'},
        'modifications' : {'[15.9949]':['M']},
        'scan_id': "spectrum_id", # code to read scan id from *.mzML file
        'mz_array': "np.float_(spectrum['m/z array'][:])", # code to read m/z array from *.mzML file
        'intensity_array': "np.float_(spectrum['intensity array'][:])", # code to read intensity array from *.mzML file
        'charge': "int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])", # code to read charge from *.mzML file
        'precursor_mass': "float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])", # code to read precursor mass from *.mzML file
    },
    'iprg' : {
        'ms_data_file' : './ms_data/iPRG_continuous_scan.mzML',
        'fasta_file' : './fasta/ABRF_iPRG_2012_target.fasta',
        'tolarence_type' : "PPM",
        'max_mods' : 1,
        'tolarence_window' : 10,
        'missed_cleavages' : 1,
        'enzyme':'trypsin',
        'static_mods': {'C+57.02146'},
        'modifications' : {'[15.9949]':['M']},
        'scan_id': "int(spectrum['id'][5:])", # code to read scan id from *.mzML file
        'mz_array': "np.float_(spectrum['m/z array'][:])", # code to read m/z array from *.mzML file
        'intensity_array': "np.float_(spectrum['intensity array'][:])", # code to read intensity array from *.mzML file
        'charge': "int(spectrum['ms2 file charge state'][0:1])", # code to read charge from *.mzML file
        'precursor_mass': "float(spectrum['ms2 file charge state'][2:])", # code to read precursor mass from *.mzML file
    },
    'yeast' : {
        'ms_data_file' : './ms_data/yeast-01_full.mzML',
        'fasta_file' : './fasta/yeast.fasta',
        'tolarence_type' : "MASS",
        'max_mods' : 1,
        'tolarence_window' : 3,
        'enzyme':'trypsin',
        'missed_cleavages' : 0,
        'static_mods': {'C+57.02146'},
        'modifications' : {},
        'scan_id': "int(spectrum['id'][5:])", # code to read scan id from *.mzML file
        'mz_array': "np.float_(spectrum['m/z array'][:])", # code to read m/z array from *.mzML file
        'intensity_array': "np.float_(spectrum['intensity array'][:])", # code to read intensity array from *.mzML file
        'charge': "int(spectrum['ms2 file charge state'][0:1])", # code to read charge from *.mzML file
        'precursor_mass': "float(spectrum['ms2 file charge state'][2:])", # code to read precursor mass from *.mzML file
    }
}

def params_to_string(params):
    """Packing key parameters to string for making proper filenames"""
    string = ''
    for key in params.keys():
        string += '-{}-{}'.format(key, params[key])
    return string

def string_to_params(string):
    """Unpacking key parameters from string"""
    model_tokens = string.split('-')
    DATASET = model_tokens[0]
    RES = model_tokens[1]
    model_name = model_tokens[2]
    
    for i,token in enumerate(model_tokens):
        if token in training.keys():
            try:
                training[token] = int(model_tokens[i+1])
            except ValueError:
                training[token] = float(model_tokens[i+1])
    return (DATASET, RES, model_name)    
