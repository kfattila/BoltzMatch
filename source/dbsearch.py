"""
This file contains source code of BoltzMatch.
Copyright 2019 by Anastasia Voronkova, Pavel Sulimov and Attila Kertesz-Farkas.

All information contained herein is, and remains
the property of authors. The intellectual and technical concepts contained
herein are proprietary to authors.
Distributed under the Apache License, Version 2.0.
"""
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import random
import re
import numpy as np
import matplotlib.pyplot as plt
from math import factorial, sqrt
from datetime import datetime
from pyteomics import mzml, parser, mass
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from source.parameters import spectrum as mzml_params
plt.switch_backend('agg')

proton_mass = 1.00727646688 # define the proton mass

class DBSearch:
    """
    Initialize the biological parameters before processing exact dataset
    """
    def __init__(self,
                 bin_width=1.0005079,
                 bin_offset=0.4,
                 remove_precursor_peak=True,
                 remove_precursor_tolerance=1.5,
                 max_peak=2000,
                 skip_preprocessing=False,
                 enzyme='trypsin',
                 missed_cleavages=0,
                 min_pept_len=7,
                 max_pept_len=30,
                 min_pept_mass=200,
                 max_pept_mass=7200,
                 max_mods=1,
                 min_mods=0,
                 decoy_format=0, #0 for reverse, 1 for shuffle, anything else means no decoy generation.
                 semi_cleavage=0, #0 for full tryptic, 1 for semi tryptic, (2 for non-tryptic not supported yet)
                 decoys_only=False, # Generate only decoy peptides, without generating target peptides.
                 modifications={},  # set()
                 static_mods={'C+57.02146'},
                 theo_pept_peaks='by',   #possible peaks: 'abcxyz'
                 max_theo_pept_peak_charge=2,
                 unique_peptides=1,     #0 no unique peptides, much more memory and redundancy,
                                        #1 for unique peptides, less memory, less peptides, but slower peptide generation
                 tolarence_type="PPM",
                 tolarence_window=10,
                 intensity_cutoff_coefficient=0.05
                 ):

        # constants
        self.elts_mono = {
            'H': 1.007825035,
            'C': 12.0,
            'N': 14.003074,
            'O': 15.99491463,
            'P': 30.973762,
            'S': 31.9720707
        }
        self.B = 0.0
        self.mono_h2o = 2 * self.elts_mono['H'] + self.elts_mono['O']
        self.Y = self.mono_h2o
        self.h2o = self.Y
        self.nh3 =3 * self.elts_mono['H'] + self.elts_mono['N']

        # spectrum processing related parameters
        self.skip_preprocessing = skip_preprocessing
        self.remove_precursor_peak = remove_precursor_peak
        self.remove_precursor_tolerance = remove_precursor_tolerance
        self.intensity_cutoff_coefficient = intensity_cutoff_coefficient
        self.max_xcorr_offset = 75
        self.max_intensity = 50.0
        self.num_spectrum_regions = 10

        # spectrum data related parameters
        self.max_peak = max_peak  # Highest peak considered. Peaks having higher m/z than this will be discarded
        self.bin_offset = bin_offset
        self.bin_width = bin_width
        self.max_bin = self.mass2bin(self.max_peak)  # maximum bin
        self.tolarence_type = tolarence_type
        self.tolarence_window = tolarence_window  #precursor_window

        # peptide generation related parameters
        self.missed_cleavages = missed_cleavages
        self.min_pept_len = min_pept_len
        self.max_pept_len = max_pept_len
        self.min_pept_mass = min_pept_mass
        self.max_pept_mass = max_pept_mass
        self.enzyme = enzyme
        self.max_mods = max_mods
        self.min_mods = min_mods
        self.modifications = modifications
        self.theo_pept_peaks = theo_pept_peaks
        self.max_theo_pept_peak_charge = max_theo_pept_peak_charge
        self.decoy_format = decoy_format
        self.decoys_only = decoys_only
        self.semi_cleavage = semi_cleavage
        self.unique_peptides = unique_peptides

        self.spectrum_collection = []
        self.peptide_collection = []
        self.protein_collection = [] # protein sequence related data
        self.pattern = re.compile(r'-?\[([^\]]+)\]-?')
        self.peptide_set = set()

        self.aa_mass = {
            'G': 57.02146,
            'A': 71.03711,
            'S': 87.03203,
            'P': 97.05276,
            'V': 99.06841,
            'T': 101.04768,
            'C': 103.00919, #103.00919+57.02146 = 160.03065
            'L': 113.08406,
            'I': 113.08406,
            'N': 114.04293,
            'D': 115.02694,
            'Q': 128.05858,
            'K': 128.09496,
            'E': 129.04259,
            'M': 131.04049,
            'H': 137.05891,
            'F': 147.06841,
            'U': 150.95364,
            'R': 156.10111,
            'Y': 163.06333,
            'W': 186.07931,
            'O': 255.15829,
            }
        
        # process static modifications
        for stat_mod in static_mods:
            if stat_mod[1] == 't': # for terminal mods: e.g. : Nt+229.99
                if stat_mod[0] == 'C':
                    self.Y += float(stat_mod[2:])
                if stat_mod[0] == 'N':
                    self.B += float(stat_mod[2:])
                continue
            self.aa_mass[stat_mod[0]] += float(stat_mod[1:])

        self.w1 = 0
        self.b1 = 0
        self.w2 = 0
        self.b2 = 0
        
        self.max_n_candidates = 0

    def load_data(self, path_to_file, min_peak_th=10, data_type='ups1'):
        """Loading experimental data from *.mzML file"""
        self.spectrum_collection = []
        # print(eval(mzml_params[data_type]['scan_id']))
        with mzml.read(path_to_file, dtype=dict) as spectra:
            for spectrum_id, spectrum in enumerate(spectra):
                spectrum_record = Spectrum(
                path_to_file,  # path to file
                eval(mzml_params[data_type]['scan_id']), #scan id
                eval(mzml_params[data_type]['mz_array']),  # mz array
                eval(mzml_params[data_type]['intensity_array']),  # intensity array
                eval(mzml_params[data_type]['charge']),  # charge
                eval(mzml_params[data_type]['precursor_mass']),  # precursor mass
                self.max_peak,
                self.remove_precursor_peak,
                self.remove_precursor_tolerance)
                if len(spectrum_record.intensity_array) >= min_peak_th:
                    self.spectrum_collection.append(spectrum_record)
        self.set_spectrum_idx()
       
    def sort_spectra(self, key = "neutral_mass", reverse=False):   #reverse=False => increasing
        """Sorting spectra with respect to mass"""
        self.spectrum_collection = sorted(self.spectrum_collection, key=lambda d: getattr(d, key), reverse=reverse)

    def sort_peptides(self, key = "neutral_mass", reverse=False):
        """Sorting peptides with respect to mass"""
        self.peptide_collection = sorted(self.peptide_collection, key=lambda d: getattr(d, key), reverse=reverse)

    def mass2bin_vec(self, mass, charge=1):
        """Convert mass to bin vector"""
        return ((mass + (charge - 1) * proton_mass) / (charge*self.bin_width) + 1.0 - self.bin_offset).astype(int)

    def mass2bin(self, mass, charge=1):
        """Convert mass to bin"""
        bin = int((mass + (charge - 1) * proton_mass) / (charge*self.bin_width) + 1.0 - self.bin_offset)
        return bin

    def bin2mass(self, bin, charge=1):
        """Convert bin to mass"""
        return (bin - 1.0 + self.bin_offset) * (charge*self.bin_width) - (charge - 1) * proton_mass

    def set_bin_size(self, bin_width, bin_offset, max_peak = 2000):
        """Setting bin size"""
        self.max_peak = max_peak  # Highest peak considered. Peaks having higher m/z than this will be discarded
        self.bin_offset = bin_offset
        self.bin_width = bin_width
        self.max_bin = self.mass2bin(self.max_peak)  # maximum bin        
        
    def discretize_spectrum(self, spectrum_id):
        """Spectrum discretization"""
        spectrum = self.spectrum_collection[spectrum_id]
        spectrum.spectrum_array = np.zeros(self.max_bin)
        spectrum.peak_bins = list(map(self.mass2bin, spectrum.mz_array))

        if not self.skip_preprocessing:
            sqrt_intensity_array = list(map(sqrt,spectrum.intensity_array))
        else:
            sqrt_intensity_array = spectrum.intensity_array

        # keep the maximum of the intensities in every bin
        for i, peak_bin in enumerate(spectrum.peak_bins):
            if peak_bin < self.max_bin and spectrum.spectrum_array[peak_bin] < sqrt_intensity_array[i]:
                spectrum.spectrum_array[peak_bin] = sqrt_intensity_array[i]

    def normalize_regions(self, spectrum_id, N=np.inf, initial=True):
        """
        Performing region normalization
        ------
        Notice: initial=True is used for initial call, in case of selecting topN peaks
        with help of this function, usw initial=false
        -------
        Warning: discretize_spectrum() must be called beforehand!
        """
        # Fill peaks
        spectrum = self.spectrum_collection[spectrum_id]
        largest_bin = max(spectrum.peak_bins)
        highest_intensity = max(spectrum.spectrum_array)
        intensity_cutoff = highest_intensity * self.intensity_cutoff_coefficient  # lower intensity cutoff
        region_size = int(largest_bin / self.num_spectrum_regions) + 1  # size of any region
        if initial == True:
            spectrum.spectrum_array[np.where(spectrum.spectrum_array <= intensity_cutoff)] = 0
        highest_intensity = [max(spectrum.spectrum_array[i*region_size : (i+1)*region_size]) for i in range(self.num_spectrum_regions)]
        
        N = N // self.num_spectrum_regions
        def normalize(id_intens):
            region_start = id_intens[0] * region_size
            if initial == True:
                spectrum.spectrum_array[region_start:region_start + region_size] *= (self.max_intensity / id_intens[1])
            if N < region_size:
                idx = spectrum.spectrum_array[region_start:region_start + region_size].argsort()[:(region_size-N)]
                spectrum.spectrum_array[idx + region_start] = 0

        list(map(normalize, filter(lambda id_intens: id_intens[1], enumerate(highest_intensity))))


    def XCORR_substract_background(self, spectrum_id):
        """
        Operation is as follows: new_observed = observed -
        average_within_window, but average is computed as if the array
        extended infinitely: denominator is same throughout array, even
        near edges (where fewer elements have been summed)
        ------
        Notice: discretize_spectrum must be called beforehand
        """
        spectrum = self.spectrum_collection[spectrum_id]

        multiplier = 1.0 / (self.max_xcorr_offset * 2)
        end = len(spectrum.spectrum_array)
        partial_sums = np.zeros(end + 1)

        partial_sums[0:end] = np.add.accumulate(spectrum.spectrum_array, axis=0)
        partial_sums[end] = partial_sums[end-1]

        l_border = self.max_xcorr_offset
        r_border = end - self.max_xcorr_offset

        partial_sums_left = np.zeros(end)
        partial_sums_left[l_border + 1:end] = partial_sums[0:r_border - 1]

        partial_sums_right = np.repeat(partial_sums[end], end)
        partial_sums_right[0:r_border] = partial_sums[l_border:end]

        spectrum.spectrum_array[:] -= multiplier * (partial_sums_right[:] - partial_sums_left[:] -
                                                                          spectrum.spectrum_array[:])

    def topN(self, spectrum_id, N=50):
        """Keep N most intensive peaks in spectrum"""
        spectrum = self.spectrum_collection[spectrum_id]
        N = np.minimum(np.count_nonzero(spectrum.spectrum_array), N)        
        ind = spectrum.spectrum_array.argsort()[:len(spectrum.spectrum_array)-N]
        spectrum.spectrum_array[ind] = 0

    def spectrum_preprocess(self, spectrum_id):
        """Applies all basic operations, common to be called before Tide search"""
        self.discretize_spectrum(spectrum_id)
        self.normalize_regions(spectrum_id)
        self.XCORR_substract_background(spectrum_id)
    
    def preprocess_all_spectra(self):
        """Applies spectrum_preprocess() to all spectra in collection"""
        print("Preprocessing all spectra (discretization, region normalization, and prepartion for XCORR scoring)...")
        start_time = datetime.now()
        for spect_id in self.spectrum_collection:
            self.spectrum_preprocess(spect_id)
        print("Spectrum preprocess done. Time (h:m:s):\t"+str(datetime.now() - start_time))
    
    def compute_window(self, mass, charge=1):
        """Computing tolerance window"""
        if self.tolarence_type == "MASS":
            out_min = mass - self.tolarence_window
            out_max = mass + self.tolarence_window
        elif self.tolarence_type == "MZ":
            mz_minus_proton = mass - proton_mass   #mass must be precursor_mass
            out_min = (mz_minus_proton - self.tolarence_window) * charge
            out_max = (mz_minus_proton + self.tolarence_window) * charge
        elif self.tolarence_type == "PPM":
            tiny_precursor = self.tolarence_window * 1e-6
            out_min = mass * (1.0 - tiny_precursor)
            out_max = mass * (1.0 + tiny_precursor)
        else:
            out_min = out_max = mass
            print("Uncorrect type of tolerance!")
        return out_min, out_max

    def reset_search_results(self):
        """Resets searh results after spectrum identification to default"""
        for spectrum in self.spectrum_collection:
            spectrum.score = -1000000  # some matching score, like hyperscore or XCORR, etc. bigger the better
            spectrum.confidence = 1000000  # some statistical confidence value, such as E-value, or (exact) p-value, smaller the better
            spectrum.qvalue = 1000000  # statistical q-value value from TDC or exact-pvalue
            # spectrum.n_candidates = 0
            spectrum.peptide = None
            spectrum.isotope = 0
       
    def generate_peptides(self, protein_id, decoy=False):
        """Generating peptides from given protein"""
        if self.enzyme == 'trypsin':
            tide_trypsin = r'([KR](?=[^P]))'
            full_peptides = set(parser.cleave(self.protein_collection[protein_id].seq, tide_trypsin,
                                          self.missed_cleavages, self.min_pept_len))
        elif self.enzyme == 'trypsin/p': # ignore proline rule for trypsin
            tide_trypsin = r'([KR])'
            full_peptides = set(parser.cleave(self.protein_collection[protein_id].seq, tide_trypsin,
                                          self.missed_cleavages, self.min_pept_len))
        elif self.enzyme == 'no-digestion':
            full_peptides = set()
            full_peptides.add(self.protein_collection[protein_id].seq)
        else:
            full_peptides = set(parser.cleave(self.protein_collection[protein_id].seq, parser.expasy_rules[self.enzyme],
                                              self.missed_cleavages, self.min_pept_len))

        peptides = []

        for peptide in full_peptides:  # check peptides
            if peptide.find('X') != -1:
                continue
            if peptide.find('Z') != -1:
                continue
            if peptide.find('B') != -1:
                continue
            pept_len = len(peptide)
            if pept_len < self.min_pept_len:
                continue
            
            if pept_len > self.max_pept_len:
                continue
                
            peptides.append(peptide)
            if self.semi_cleavage == 1:
                for j in range(1, len(peptide)):  # choose the part of peptide
                    rev_j = len(peptide) - j

                    if self.min_pept_len <= rev_j:
                        peptides.append(peptide[j:])

                    if self.min_pept_len <= j:
                        peptides.append(peptide[:rev_j])
        peptides = set(peptides)                
        for peptide in peptides:
            if decoy:
                self.add_peptide_collection(peptide, protein_id, 0)
                continue

            if not self.decoys_only:
                self.add_peptide_collection(peptide, protein_id, 1)

            decoy_peptide = self.get_decoy(peptide, self.decoy_format)
            if decoy_peptide != None:
                self.add_peptide_collection(decoy_peptide, protein_id, 0)
                
    def get_decoy(self, peptide, format=0):
        """Provides decoy peptide from given target"""
        if format == 0: #reverse
            middle = peptide[1:-1][::-1]
            return peptide[0] + middle + peptide[-1]
        if format == 1: #shuffle
            middle = list(peptide[1:-1])
            random.shuffle(middle)
            return peptide[0] + "".join(middle) + peptide[-1]
        return None

    def add_peptide_collection(self, peptide, protein_id, target):
        """Adding PeptideObj() to peptide collection"""
        if self.unique_peptides == 1 and peptide in self.peptide_set:
            return
        #if target == 1:
        self.peptide_set.add(peptide)        

        modified_peptides = set(parser.isoforms(peptide, variable_mods=self.modifications, max_mods=self.max_mods))

        peptide_aa_mass = np.array([self.aa_mass[aa] for aa in peptide])  #Make this faster using AA as ubytes
        for mod_pept in modified_peptides:
        
            mod_cnt = len(re.findall(r"\]", mod_pept))
            if mod_cnt < self.min_mods or mod_cnt > self.max_mods:
                continue
            
            mod_peptide_aa_mass = list(peptide_aa_mass)   #peptide_aa_mass[:] does not work
            offset = 0
            for mod in re.finditer(self.pattern, mod_pept):  # calc modifications' mass
                location = mod.start()

                if mod_pept[location] == "-":
                    mod_peptide_aa_mass[-1] += float(mod.group(1))
                    break

                mod_peptide_aa_mass[location - offset] += float(mod.group(1))
                offset += mod.end() - mod.start()

            # mod_pep_mass = np.sum(mod_peptide_aa_mass) + self.mono_h2o #H and OH 
            mod_pep_mass = np.sum(mod_peptide_aa_mass) + self.Y + self.B #H and OH and static modes
            if mod_pep_mass < self.min_pept_mass or mod_pep_mass >= self.max_pept_mass:
                continue

            pept_obj = PeptideObj(mod_pep_mass, mod_pept, protein_id, target, peptide, "full",
                                                    self.missed_cleavages)
            pept_obj.aa_mass = mod_peptide_aa_mass
            #print(mod_pep_mass, mod_pept, target)
            self.peptide_collection.append(pept_obj)
        
    def calculate_peptide_fragmentation(self, peptide_id):
        """Calculating masses of peptide fragment ions"""
        peptide = self.peptide_collection[peptide_id]
        peptide.peaks = [dict() for x in range(self.max_theo_pept_peak_charge)]
        
        for ion_series in self.theo_pept_peaks:

            if ion_series == 'b':  # generate B ions.
                fragment_ions = np.cumsum(peptide.aa_mass[:-1]) + (self.B + proton_mass)
            if ion_series == 'y':  # generate Y ions.
                fragment_ions = np.cumsum(peptide.aa_mass[1:][::-1]) + (self.Y + proton_mass)
                
            for peak_charge in range(self.max_theo_pept_peak_charge):
                fragment_idx = self.mass2bin_vec(fragment_ions, peak_charge + 1)
                peak_list = list(filter(lambda peak: peak < self.max_bin, fragment_idx ))
                peptide.peaks[peak_charge][ion_series] = peak_list

    def set_candidate_peptides(self):
        """
        Setting candidate peptides
        ------
        Notice:
        1) Spectrum collection must be sorted
        2) Peptide collection must be sorted
        """
        start_pept_id = 0
        end_pept_id = 0        
        pept_num = len(self.peptide_collection)
        for spect_id in range(len(self.spectrum_collection)):
            spectrum = self.spectrum_collection[spect_id]
            min_mass, max_mass = self.compute_window(spectrum.neutral_mass, spectrum.charge)

            #finding candidate peptides using rolling window approach. peptides are assumed to be sorted by neutral mass
            if start_pept_id >= pept_num:
                break
            while self.peptide_collection[start_pept_id].neutral_mass < min_mass:
                start_pept_id += 1
                if start_pept_id >= pept_num:
                    break

            if end_pept_id < start_pept_id:
                end_pept_id = start_pept_id

            if end_pept_id >= pept_num:
                break
            while self.peptide_collection[end_pept_id].neutral_mass < max_mass:
                end_pept_id += 1
                if end_pept_id >= pept_num:
                    break
            self.spectrum_collection[spect_id].n_candidates = end_pept_id-start_pept_id
            self.spectrum_collection[spect_id].start_pept = start_pept_id
            self.spectrum_collection[spect_id].end_pept   = end_pept_id

            if self.max_n_candidates < self.spectrum_collection[spect_id].n_candidates:
                self.max_n_candidates = self.spectrum_collection[spect_id].n_candidates

    def load_fasta(self, path_to_fasta):
        """Loading data from *.fasta file"""
        cnt = 0
        for record in SeqIO.parse(path_to_fasta, "fasta"):
            record.seq.alphabet = generic_protein
            self.protein_collection.append(ProteinObj(cnt, str(record.id), str(record.seq)))
            cnt += 1
            
    def calculate_xcorr_score(self, spect_id, pept_id):
        """Calculating XCORR score"""
        spectrum = self.spectrum_collection[spect_id]

        if not self.peptide_collection[pept_id].peaks:
            self.calculate_peptide_fragmentation(pept_id)

        peaks = []
        for key, peak_list in self.peptide_collection[pept_id].peaks[0].items():   # Match sinlge charged theoretical peaks
            peaks += peak_list

        if spectrum.charge > 2 and self.max_theo_pept_peak_charge > 1:
            for key, peak_list in self.peptide_collection[pept_id].peaks[1].items():  # Match double charged theoretical peaks.
                peaks += peak_list
        peaks = np.unique(peaks)
        score = np.sum(spectrum.spectrum_array[peaks])

        return score/200

 
    def tide_search(self):
        """
        Performing Tide search, python reimplementation of Tide-Search from CRUX. See: crux.ms and https://noble.gs.washington.edu/papers/diamen2011faster.pdf
        ------
        Notice:
        1) Assume that the protein fasta  and the spectrum data files are loaded
        2) Assume that all spectra are discretized, normalized and preprocessed with XCORR_substract_background
        3) Spectrum_collection must be sorted by neutral mass in increasing order
        """
       
        for spect_id in range(len(self.spectrum_collection)):
            spectrum = self.spectrum_collection[spect_id]

            for pept_id in range(spectrum.start_pept, spectrum.end_pept):
                
                xcorr = self.calculate_xcorr_score(spect_id, pept_id)*50
                if xcorr > spectrum.score:
                    spectrum.score = xcorr
                    self.peptide_collection[pept_id].peptide_id = pept_id
                    spectrum.peptide = self.peptide_collection[pept_id]

    def boltzmatch_tailor_scoring(self):
        """
        Performing BoltzMatch scoring using Tailor score calibration
        ------
        Notice:
        1) Assume that the protein fasta  and the spectrum data files are loaded
        2) Assume that all spectra are discretized and normalized
        """
        min_candidates = 20 # set minimum number of candidate peptides for spectrum
               
        for spect_id in range(len(self.spectrum_collection)):
            spectrum = self.spectrum_collection[spect_id]

            if spectrum.n_candidates == 0:
                continue

            start_id = spectrum.start_pept
            end_id = spectrum.end_pept

            if spectrum.n_candidates < min_candidates:
                end_id = start_id + min_candidates
            if end_id >= len(self.peptide_collection):
                end_id = len(self.peptide_collection)
            if (end_id - start_id) < min_candidates:
                start_id = end_id - min_candidates
            if start_id < 0:
                start_id = 0

            scores = np.zeros(end_id - start_id)

            for i, pept_id in enumerate(range(start_id, end_id)):
                scores[i] = self.calculate_xcorr_score(spect_id, pept_id)*50 + spectrum.bias + self.peptide_collection[pept_id].bias
            scores += 10 # need scores to be positive in order to apply division operation in Tailor scoring
            top_hits = max(int(len(scores)*0.05), 5)
            scores_norm = np.sort(scores)[-top_hits]
            norm_scores = scores / scores_norm

            for i, pept_id in enumerate(range(spectrum.start_pept, spectrum.end_pept)):
                if norm_scores[i] > spectrum.score:
                    spectrum.score = norm_scores[i]
                    self.peptide_collection[pept_id].peptide_id = pept_id
                    spectrum.peptide = self.peptide_collection[pept_id]
              
    def print_results(self, filename):
        """Saving spectrum idenitification search-results to tsv file"""
        fout = open(filename, "w")

        names_of_columns = ["file", "scan", "spectrum_id", "charge", "spectrum precursor m/z", "spectrum neutral mass",
                            "peptide mass",
            "score", "confidence", "qvalue", "number of candidates", "target", "protein id", "peptide_id", "peptide sequence", "peptide length",
            "modifications", "cleavage type", "missed cleavages", "original sequence", "quintile"]

        header = "\t".join(names_of_columns) + "\n"
        fout.writelines(header)

        for spectrum_id, spectrum in enumerate(self.spectrum_collection):
            result_string = spectrum.print_spectrum(spectrum_id)
            if result_string:
                fout.writelines(list( "{:d}\t".format(item) if
                type(item) == int else "{:f}\t".format(item) if type(item) == float else "{:s}\t".format(item) for item in result_string))
                fout.write("\n")

        fout.close()

    def compute_qvalues_tdc(self):
        """Computing q-values after scoring applied along the search using concatenated target-decoy approach"""
        self.sort_spectra(key="score", reverse=True)
        
        target_cnt = 0
        decoy_cnt = 1
        for spectrum in self.spectrum_collection:
            if not spectrum.peptide:
                break
            if spectrum.peptide.target == True:
                target_cnt += 1
            else:
                decoy_cnt += 1
            fdr = decoy_cnt /  (target_cnt + 1)
            if fdr > 1.0:
                fdr = 1.0
            spectrum.qvalue = fdr

        #convert fdrs to qvalues:         
        for i in range(len(self.spectrum_collection)-2,-1, -1):
            if self.spectrum_collection[i+1].qvalue <  self.spectrum_collection[i].qvalue:
                self.spectrum_collection[i].qvalue = self.spectrum_collection[i+1].qvalue
    
    def plot_qvalues(self, mode='show', filename=None, qval_lim=0.1):
        """Plotting FDR curve"""
        fig = plt.figure(figsize=(10, 6))

        x = []
        y = []
        qval = 0
        accepted_psm = 0
        for spectrum in self.spectrum_collection:
            if qval > qval_lim:
                break
            if qval != spectrum.qvalue:
                x.append(qval)
                y.append(accepted_psm)
                qval = spectrum.qvalue
            try:
                if spectrum.peptide.target == True:
                    accepted_psm += 1
            except:                    
                pass

        x.append(qval)
        y.append(accepted_psm)

        plt.plot(x,y)
        plt.xlim([0, qval_lim])
        # plt.ylim([0, 2500])
        plt.ylabel('Number of accepted spectra')
        plt.xlabel('Estimated False Discovery Rate (Q-value)')
        fig.patch.set_facecolor('xkcd:white')
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            fig.savefig(filename)
            plt.clf()  
        plt.close()

    def print_accepted_psms(self):
        """Printing PSMs accepted at different level"""
        count_1_percent = 0
        count_5_percent = 0
        count_10_percent = 0
        for spectrum in self.spectrum_collection:
            if spectrum.peptide is None:
                break
            if spectrum.peptide.target == False:
                continue
            if spectrum.qvalue < 0.01:
                count_1_percent += 1
            if spectrum.qvalue < 0.05:
                count_5_percent += 1
            if spectrum.qvalue < 0.10:
                count_10_percent += 1
        counts =  (count_1_percent, count_5_percent, count_10_percent)
        print('Number of accepted PSMs at\t1%% FDR = %d,\t5%% FDR = %d,\t10%% FDR = %d'%counts)
        return counts      
    
    def export_spectra_ms2(self, filename):
        """Exporting results to *.ms2 format"""       
        fp = open(filename, 'w')
        for cnt in range(len(self.spectrum_collection)):
        
            #Print header
            spectrum = self.spectrum_collection[cnt]
            fp.write('S\t%d\t%d\t%lf\n'%(spectrum.scan_id, spectrum.scan_id, spectrum.precursor_mass))
            fp.write('Z\t%d\t%lf\n'%(spectrum.charge, (spectrum.precursor_mass - proton_mass)*spectrum.charge + proton_mass))
            peak_bin_idx = np.nonzero(spectrum.spectrum_array)[0]
            bin_to_peak = [self.bin2mass(peak_bin)+self.bin_width/2 for peak_bin in peak_bin_idx]
            for i,peak in enumerate(bin_to_peak):
                peak_intensity = spectrum.spectrum_array[peak_bin_idx[i]]
                if peak > 0:   
                    fp.write('%f %f\n'%(peak,peak_intensity))
        fp.close()   

    def peptide_batch_generator(self, batch_size, shuffle=True):
        """Generating batch with peptides"""
        protein_id = 0
        peptide_id = 0
        self.peptide_set = set()
        self.peptide_collection = []

        if shuffle:
            random.shuffle(self.protein_collection)
        
        while True:
            if len(self.peptide_collection[peptide_id:]) < batch_size:
                self.peptide_collection = self.peptide_collection[peptide_id:]
            while len(self.peptide_collection[peptide_id:]) < batch_size:
                if protein_id < len(self.protein_collection):
                    self.generate_peptides(protein_id)
                    protein_id += 1
                else:
                    yield self.peptide_collection[peptide_id:]
                    return
            yield self.peptide_collection[peptide_id : peptide_id + batch_size]
            peptide_id += batch_size

    def spectrum_batch_generator(self, batch_size, shuffle=True):
        """Generating batch with spectra"""
        spectrum_num = len(self.spectrum_collection)
        init_indicies = np.arange(spectrum_num)
        spectrum_id = 0

        if shuffle:
            indicies = np.array([], dtype=np.int16)
            half_len = len(init_indicies)//2
            for i in range(half_len):
                indicies = np.append(indicies, [init_indicies[i], init_indicies[i + half_len]])
        else:
            indicies = init_indicies

        while spectrum_id + batch_size < spectrum_num:
            yield indicies[spectrum_id:spectrum_id+batch_size]
            spectrum_id += batch_size      
        else:
            indicies[spectrum_id:]

    def set_spectrum_idx(self):
        """Setting spectrum index"""
        id = 0
        for spectrum in self.spectrum_collection:
            spectrum.id = id
            id += 1

    def get_spectrum_idx(self, id):
        """Getting spectrum by index"""
        for cnt, spectrum in enumerate(self.spectrum_collection):
            if spectrum.id == id:
                return cnt

    def get_spectrum_by_scan(self, scan):
        """Getting spectrum by scan_id"""
        for cnt, spectrum in enumerate(self.spectrum_collection):
            if spectrum.scan_id == scan:
                return cnt

class Spectrum(object):
    """
    Creating experimental spectrum object
    """
    def __init__(self, path_to_file, scan_id, mz_array, intensity_array, charge, precursor_mass, max_peak,
                 remove_precursor_peak, remove_precursor_tolerance):

        # spectrum info
        name_start = path_to_file.rfind('/') + 1
        self.id = 0
        self.path_to_file = path_to_file[name_start:]
        self.scan_id = scan_id
        self.charge = charge
        self.precursor_mass = precursor_mass
        self.neutral_mass = (self.precursor_mass - proton_mass) * self.charge

        mask = np.all([intensity_array > 1e-10], axis=0)
        mask3 = np.all([mz_array < max_peak], axis=0)
        mask2 = np.all([mz_array < self.neutral_mass + 50], axis=0)
        if remove_precursor_peak:
            mask1 = np.all([abs(mz_array - precursor_mass) > remove_precursor_tolerance], axis=0)
        else:
            mask1 = True

        self.mz_array = mz_array[mask & mask1 & mask2 & mask3]  # List of mz peaks
        self.intensity_array = intensity_array[mask & mask1 & mask2 & mask3]  # List of the intensities of the mz peaks

        self.spectrum_array = []  # vector for the discretized spectrum, depending on bin_width and max_mz parameters
        self.peak_bins = []  # list of peak bins


        # peptide-spectrum-match info
        self.threshold = 0 # some metrics for sorting
        self.peak_list = [] # candidate peaks for self supervision
        self.quintile = -1000000 # value of 20%-percentile score
        self.raw_score = -1000000 # raw score for any scoring method
        self.score = -1000000  # some matching score, like hyperscore or XCORR, etc.
        self.confidence = 1000000  # some statistical confidence value, such as E-value, or (exact) p-value
        self.qvalue = 1000000  # statistical q-value value from TDC or exact-pvalue
        self.n_candidates = 0
        self.peptide = None
        self.isotope = 0
        self.start_pept = 0
        self.end_pept = 0
        self.bias = 0

    def print_spectrum(self, spectrum_id):
        """Printing spectrum key info and features"""
        if self.score < -100:
            return 

        modifications = ""
        result_list = [self.path_to_file, self.scan_id, spectrum_id, self.charge, self.precursor_mass, self.neutral_mass,
            float(self.peptide.neutral_mass), float(self.score), self.confidence, self.qvalue, self.n_candidates,
            int(self.peptide.target), self.peptide.protein_id, self.peptide.peptide_id, self.peptide.peptide_seq, len(self.peptide.peptide_seq),
                       modifications, self.peptide.cleavage_info, self.peptide.missed_cleavages, self.peptide.original_peptide_seq, float(self.quintile)]

        return result_list

class PeptideObj(object):
    """
    Creating peptide object right from protein
    """
    def __init__(self, neutral_mass, peptide_seq, protein_id, target, original_peptide_seq, cleavage_info,
                 missed_cleavages):

        self.neutral_mass = neutral_mass
        self.peptide_seq = peptide_seq
        self.protein_id = protein_id
        self.peptide_id = 0
        self.aa_mass = []
        self.target = target  # indicates if the peptide is target: True/False
        self.original_peptide_seq = original_peptide_seq
        self.cleavage_info = cleavage_info  # indicates the type of the cleavage which generated the peptide: full/semi
        self.missed_cleavages = missed_cleavages
        self.peaks = []
        self.weight = -1
        self.bias = 0

class ProteinObj(object):
    """
    Creating protein object right from *.fasta file information
    """
    def __init__(self, protein_id, protein_header, protein_seq, protein_flag = 0):

        self.id = protein_id
        self.header = protein_header
        self.seq = protein_seq
        self.flag = protein_flag  # can be used to indicate something, which can be used to filter protein sequences
