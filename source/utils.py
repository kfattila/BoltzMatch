"""
This file contains source code of BoltzMatch.
Copyright 2019 by Anastasia Voronkova, Pavel Sulimov and Attila Kertesz-Farkas.

All information contained herein is, and remains
the property of authors. The intellectual and technical concepts contained
herein are proprietary to authors.
Distributed under the Apache License, Version 2.0.
"""
import time
import torch
import os.path
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from math import trunc
from source.parameters import types 

# Define types of variables from parameters
dtype = types['dtype']
torch_device = types['torch_device']
float_type = types['float_type']

def prepare_batch_data(dbs, batch_idx, cand_from_search_qval=-1, min_candidates=0, topN_exp = 0):
    """
    Preparing single batch for training RBM
    
    Parameters
    ----------
    dbs: Object of class DBSearch() with predefined biological parameters
    batch_idx: subsample of spectra indices from experimental data
    cand_from_search_qval: q-value threshold for candidate peptides
    min_candidates: minimum number of theoretical spectra candidates for each experimental spectrum
    topN_exp: value for N most intense peaks to keep from raw spectrum
    
    Returns
    -------
    experimental_spectra: array of experimental spectra
    candidate_spectra: array of candidate spectra
    theoretical_spectra:list of theoretical spectra joined into one matrix corresponding to one experimental spectrum
    """
    batch_size = len(batch_idx)
    dim = int(dbs.max_bin)

    experimental_spectra = np.zeros((batch_size, dim))
    candidate_spectra = np.zeros((batch_size, dim))
    theoretical_spectra = []
    total_theo_pep = 0

    for spect_cnt, spectrum_id in enumerate(batch_idx):
        spectrum = dbs.spectrum_collection[spectrum_id]
        num_theo_pep = spectrum.end_pept - spectrum.start_pept
        if num_theo_pep == 0 or (cand_from_search_qval > -1 and spectrum.qvalue >= cand_from_search_qval):
            theoretical_spectra.append([])
            continue

        # Preprocess experimental spectrum
        dbs.discretize_spectrum(spectrum_id)
        dbs.normalize_regions(spectrum_id)

        # Prepare experimental spectrum
        if topN_exp > 0:
            dbs.normalize_regions(spectrum_id, N=topN_exp, initial=False)
        experimental_spectra[spect_cnt, :] = spectrum.spectrum_array

        # Use the best theoretical spectra from the previous search if its qvalue is better than the threshold            
        if spectrum.qvalue < cand_from_search_qval:
            peak_list = spectrum.peptide.unique_peaks
            candidate_spectra[spect_cnt, peak_list] = 1

        # Create theoretical spectra
        theo_pept_id_array = range(spectrum.start_pept, spectrum.end_pept)

        if num_theo_pep < min_candidates:
            start_id = spectrum.start_pept
            num_theo_pep = min_candidates
            end_id = start_id + min_candidates

            if end_id >= len(dbs.peptide_collection):
                end_id = len(dbs.peptide_collection)
                num_theo_pep = end_id - start_id   
            theo_pept_id_array = range(start_id, end_id)
                     
        theo_peaks = torch.zeros((num_theo_pep, dim), dtype=float_type, device=torch_device)
        for pept_cnt, peptide_id in enumerate(theo_pept_id_array):
            theo_peaks[pept_cnt, dbs.peptide_collection[peptide_id].unique_peaks] = dbs.max_intensity

        theoretical_spectra.append(theo_peaks)
        total_theo_pep += num_theo_pep

        # Delete array to free up space
        spectrum.spectrum_array = []

    if total_theo_pep == 0:
        return None, None, None
    return experimental_spectra, candidate_spectra, theoretical_spectra

### Biological data trasformation block ###

def spectrum_normalization(dbs):
    """Basic spectrum normalization for all spectra in collection"""
    for spectrum_id in range(len(dbs.spectrum_collection)):
        dbs.discretize_spectrum(spectrum_id)
        dbs.normalize_regions(spectrum_id)

def spectrum_discretization(dbs):
    """Basic spectrum discretization for all spectra in collection"""
    for spectrum_id in range(len(dbs.spectrum_collection)):
        dbs.discretize_spectrum(spectrum_id)

def spectrum_substract_background(dbs):
    """Background subtraction operation for all spectra in collection"""
    for spectrum_id in range(len(dbs.spectrum_collection)):
        dbs.XCORR_substract_background(spectrum_id)

def transofrm_spectra(dbs, rbm_model):
    """BoltzMatch transformation for all spectra in collection"""
    for spectrum_id in range(len(dbs.spectrum_collection)):
        dbs.discretize_spectrum(spectrum_id)
        dbs.normalize_regions(spectrum_id)

        visible_data = torch.from_numpy(dbs.spectrum_collection[spectrum_id].spectrum_array.reshape(1,-1)).type(dtype).to(torch_device)
        transformed = torch.mm(visible_data, rbm_model.rbm_weights)
        spectrum_bias = torch.mv(visible_data, rbm_model.rbm_vb)[0]

        dbs.spectrum_collection[spectrum_id].spectrum_array = transformed[0].data.cpu().numpy()
        dbs.spectrum_collection[spectrum_id].bias = spectrum_bias.data.cpu().numpy()
    
    peptide_bias = rbm_model.rbm_hb.data.cpu().numpy()

    for peptide_id in range(len(dbs.peptide_collection)):
        dbs.peptide_collection[peptide_id].bias = np.sum(peptide_bias[dbs.peptide_collection[peptide_id].unique_peaks])

### Searching block ###

def run_boltzmatch_search(model, dbs):
    """BoltzMatch search for final results testing"""
    dbs.reset_search_results()
    # dbs.sort_spectra()
    transofrm_spectra(dbs, model)
    dbs.boltzmatch_tailor_scoring()
    dbs.compute_qvalues_tdc()

def run_tide(dbs, filename):
    """Tide search"""
    dbs.reset_search_results()
    dbs.sort_spectra()
    dbs.tide_search()
    dbs.compute_qvalues_tdc()
    dbs.print_results(filename)

def run_tailor(dbs, filename):
    """BoltzMatch search for basic Tailor scoring"""
    dbs.reset_search_results()
    dbs.sort_spectra()
    dbs.boltzmatch_tailor_scoring()
    dbs.compute_qvalues_tdc()
    dbs.print_results(filename)

### Biological data preprocessing block###

def load_peptides_and_spectra(dbs, fasta, spectra, data_type):
    """Loading data from *.mzML and *.fasta files"""
    # print("Loading and preprocessing spectra...")
    dbs.load_data(spectra, data_type=data_type)
    dbs.sort_spectra()

    # print("Loading protein sequences and generating peptides...")
    dbs.load_fasta(fasta)    
    for i in range(len(dbs.protein_collection)):
        dbs.generate_peptides(i)
    dbs.sort_peptides()

    # set up the range of the candidate peptides for every experimental spectrum
    dbs.set_candidate_peptides()

def generate_unique_theoretical_peaks(dbs):
    """Generating theoretical spectra"""
    for peptide_id in range(len(dbs.peptide_collection)):
        dbs.calculate_peptide_fragmentation(peptide_id)
        peaks_b_0 = np.asarray(dbs.peptide_collection[peptide_id].peaks[0]['b'], dtype=np.int16)
        peaks_y_0 = np.asarray(dbs.peptide_collection[peptide_id].peaks[0]['y'], dtype=np.int16)
        if dbs.max_theo_pept_peak_charge == 2:
            peaks_b_1 = np.asarray(dbs.peptide_collection[peptide_id].peaks[1]['b'], dtype=np.int16)
            peaks_y_1 = np.asarray(dbs.peptide_collection[peptide_id].peaks[1]['y'], dtype=np.int16)
            peak_list = np.union1d(np.union1d(peaks_b_0, peaks_y_0), np.union1d(peaks_b_1, peaks_y_1))  
        else:
            peak_list = np.union1d(peaks_b_0, peaks_y_0)
        dbs.peptide_collection[peptide_id].unique_peaks = peak_list

### Results postprocessing block ###

def save_learning_curve(learning_curves, filename_stem):
    """Saving learning curves to *.csv"""
    learning_curves.to_csv(path_or_buf=filename_stem+".learning_curve")

def plot_weights(model_path, output_file_name, scale=0.001):
    """Plotting BoltzMatch weights matrix"""
    from mpl_toolkits.axes_grid1 import make_axes_locatable    
    FONT_SIZE = 24

    # Loading model and extracting weights matrix
    rbm_model = torch.load(model_path)
    weight = rbm_model['weight']

    fig=plt.figure(figsize=(10,10))
    fig.patch.set_facecolor('xkcd:white')

    plt.rc('font', size=FONT_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=FONT_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=FONT_SIZE)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=FONT_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=FONT_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=FONT_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=FONT_SIZE)   # fontsize of the figure title

    w = weight.cpu().detach().numpy()
    im = plt.imshow(w, vmin=-scale, vmax=scale)
    ax = plt.gca()
    ticks = np.arange(0,2500,step=500)
    plt.xticks(ticks)
    plt.yticks(ticks)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.5)

    plt.colorbar(im, cax=cax)
    
    fig.savefig(output_file_name+'_weight.png', bbox_inches='tight')
    fig.savefig(output_file_name+'_weight.pdf', bbox_inches='tight')
    np.savetxt(output_file_name+'_weight.txt', w)

def calculate_fdr_curve(qvalue, target):
    """Calculations for FDR curve"""
    x = []
    y = []
    qval_lim = 0.3
    prev_qval = qvalue[0]
    accepted_psm = 0
    for i in range(len(qvalue)):
        qval = qvalue[i]
        if qval > qval_lim:
            break
        if qval != prev_qval:
            x.append(prev_qval)
            y.append(accepted_psm)
            prev_qval = qval
        try:
            if target[i] == 1:
                accepted_psm += 1
        except:                    
            pass
    x.append(qval)
    y.append(accepted_psm)
    return x,y    

# def plot_fdr_curves(filename, ylim_min, ylim_max, line_width=4, font_size=24):
    # """Plotting FDR curves"""
    # plt.rc('font', size=font_size)          # controls default text sizes
    # plt.rc('axes', titlesize=font_size)     # fontsize of the axes title
    # plt.rc('axes', labelsize=font_size)     # fontsize of the x and y labels
    # plt.rc('xtick', labelsize=font_size-4)    # fontsize of the tick labels
    # plt.rc('ytick', labelsize=font_size-4)    # fontsize of the tick labels
    # plt.rc('legend', fontsize=font_size)    # legend fontsize
    # plt.rc('figure', titlesize=font_size)   # fontsize of the figure title
    # fig = plt.figure(figsize=(10,10))
    # fig.patch.set_facecolor('xkcd:white')
 
    # df = pd.read_csv(filename, sep='\t', index_col=False)
    # scores = pd.DataFrame(df, columns=['qvalue', 'target']).values
    # x,y = calculate_fdr_curve(scores[:,0], scores[:,1])
    # plt.plot(x,y, linewidth=line_width)
    
    # plt.xlim((0,0.1))
    # plt.xticks(np.arange(0.0,0.11,step=0.01))
    # plt.ylim((ylim_min, ylim_max))
    # plt.xlabel('FDR (q-values threshold)')
    # plt.ylabel('Number of accepted PSMs')
    # plt.legend(('BoltzMatch', 'XCorr scoring', 'Scalar  scoring', 'HyperScore'), loc='lower right')
    # plt.grid(True)
    # plt.show()    
    # fig.savefig(filename+".pdf", bbox_inches='tight')


# def plot_learning_curves(filename, ylim_min, ylim_max, line_width=4, font_size=24):
    # """Plotting learning curves"""
    # plt.rc('font', size=font_size)          # controls default text sizes
    # plt.rc('axes', titlesize=font_size)     # fontsize of the axes title
    # plt.rc('axes', labelsize=font_size)     # fontsize of the x and y labels
    # plt.rc('xtick', labelsize=font_size)    # fontsize of the tick labels
    # plt.rc('ytick', labelsize=font_size)    # fontsize of the tick labels
    # plt.rc('legend', fontsize=font_size)    # legend fontsize
    # plt.rc('figure', titlesize=font_size)   # fontsize of the figure title

    # df = pd.read_csv(filename).values

    # hyper_scoring = df[0,2:]
    # scalar_scoring = df[1,2:]
    # xcorr_scoring  = df[2,2:]
    # data_num = len(df)-4
    # psms_1 = df[3:,2]
    # psms_5 = df[3:,3]
    # psms_10 = df[3:,4]
    # fig = plt.figure(figsize=(30,10))
    # fig.patch.set_facecolor('xkcd:white')
    # ax1 = plt.subplot(131)
    # col_idx = 0
    # ax1.plot(psms_1, linewidth=line_width)
    # ax1.plot([0, data_num], [xcorr_scoring[col_idx], xcorr_scoring[col_idx]], linewidth=line_width)
    # ax1.plot([0, data_num], [scalar_scoring[col_idx], scalar_scoring[col_idx]], linewidth=line_width)
    # ax1.plot([0, data_num], [hyper_scoring[col_idx], hyper_scoring[col_idx]], linewidth=line_width)
    # ax1.set_title('FDR = 1 %')
    # ax1.set_xlabel("Epoch")
    # ax1.set_ylabel("Number of accepted PSMs")
    # ax1.set_ylim(ylim_min, ylim_max)
    # ax1.set_xlim(0, data_num)
    # ax1.grid(True)

    # col_idx = 1
    # ax2 = plt.subplot(132)
    # ax2.plot(psms_5, linewidth=line_width)
    # ax2.plot([0, data_num], [xcorr_scoring[col_idx], xcorr_scoring[col_idx]], linewidth=line_width)
    # ax2.plot([0, data_num], [scalar_scoring[col_idx], scalar_scoring[col_idx]], linewidth=line_width)
    # ax2.plot([0, data_num], [hyper_scoring[col_idx], hyper_scoring[col_idx]], linewidth=line_width)
    # ax2.set_title('FDR = 5 %')
    # ax2.set_xlabel("Epoch")
    # ax2.set_ylim(ylim_min, ylim_max)
    # ax2.set_xlim(0, data_num)
    # ax2.grid(True)

    # col_idx = 2
    # ax3 = plt.subplot(133)
    # ax3.plot(psms_10, linewidth=line_width)
    # ax3.plot([0, data_num], [xcorr_scoring[col_idx], xcorr_scoring[col_idx]], linewidth=line_width)
    # ax3.plot([0, data_num], [scalar_scoring[col_idx], scalar_scoring[col_idx]], linewidth=line_width)
    # ax3.plot([0, data_num], [hyper_scoring[col_idx], hyper_scoring[col_idx]], linewidth=line_width)
    # ax3.set_title('FDR = 10 %')
    # ax3.set_xlabel("Epoch")
    # ax3.set_ylim(ylim_min, ylim_max)
    # ax3.set_xlim(0, data_num)
    # ax3.legend(('BoltzMatch', 'XCorr scoring', 'Scalar  scoring', 'HyperScore'), loc='lower right')
    # ax3.grid(True)
    # plt.show()
    # fig.savefig(filename+".pdf", bbox_inches='tight')

### Supplementary functions block ###

def process_timedelta(a):
    """Processing timedelta to be shown in hours, minutes and seconds separately"""
    if a.seconds<60:
        return str(a.seconds)+' sec'
    elif a.seconds<3600:
        return str(trunc(a.seconds/60))+' min '+str(a.seconds%60)+' sec'
    else:
        return str(trunc(a.seconds/3600))+' h '+str(trunc((a.seconds%3600)/60))+' min '+str(a.seconds%60)+' sec'

def Log(string, time, logfile):
    """Writing logs into text file"""
    f = open(logfile, 'a')
    print(string)
    print(time)
    f.write(string+', at time: '+str(time)+'\n')
    f.close()

def choose_indices(batch_size, num_reg_pairs=None):
    """Choosing indices for data pairs for diversifying regularization"""
    Xp_indeces = []
    Xq_indeces = []
    if num_reg_pairs == None:
        num_reg_pairs = batch_size
    for item in combinations(range(num_reg_pairs), 2):
        Xp_indeces.append(item[0])
        Xq_indeces.append(item[1])
    chosen_indices1 = np.array(Xp_indeces, dtype='int32')
    chosen_indices2 = np.array(Xq_indeces, dtype='int32')

    return chosen_indices1, chosen_indices2

def calculate_div_reg(rbm_model, experimental_spectra, chosen_indices1, chosen_indices2, div_reg_rate=0.01):
    """Calculating diversifying regularization"""
    visible_data = torch.from_numpy(experimental_spectra).type(dtype).to(torch_device)
    hidden_representation = torch.matmul(visible_data, rbm_model.rbm_weights) + rbm_model.rbm_hb
    hidden_representation1 = hidden_representation[chosen_indices1, ]
    hidden_representation2 = hidden_representation[chosen_indices2, ]

    return abs(div_reg_rate*torch.mean(hidden_representation1*hidden_representation2))
