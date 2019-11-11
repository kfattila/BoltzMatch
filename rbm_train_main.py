"""
This file contains source code of BoltzMatch.
Copyright 2019 by Anastasia Voronkova, Pavel Sulimov and Attila Kertesz-Farkas.

All information contained herein is, and remains
the property of authors. The intellectual and technical concepts contained
herein are proprietary to authors.
Distributed under the Apache License, Version 2.0.
"""
import time
import argparse
import numpy as np
import pandas as pd
import source.dbsearch as dbsearch
import source.parameters as params
import source.utils as utils
import source.models as models
import torch
import torch.optim as optim

# Define types of variables from parameters
dtype = params.types['dtype']
torch_device = params.types['torch_device']
float_type = params.types['float_type']
# Training parameters
training_parameters = params.training

# Define log-softmax function in advance
log_softmax = torch.nn.LogSoftmax(dim=0)

epoch_printing_tick = 20 # print learning update after this many epochs
batch_printing_tick = 1000 # print learning update after this many batches


def fit(DATASET, RES):
    """
    Fits the RBM parameters for BoltzMatch transformation (rbm_weights, rbm_hb, rbm_vb).
    
    Parameters
    ----------
    DATASET: The dataset to be used for training BoltzMatch matrix.
    RES: Resolution of the data
    
    Returns
    -------
    filename: Stem of filename. Could be used for further optional learning curves/FDR curves/weights matrix plotting
    """
    print('Dataset: {}, resolution: {}'.format(DATASET, RES))
    # Creating filename at first for saving learning curves, models etc.
    filename_stem = './{}-{}'.format(DATASET, RES)
    param_string = params.params_to_string(training_parameters)
    filename = './results/{}-{}'.format(filename_stem, param_string)
    
    # Import DBSearch parameters
    dbs = dbsearch.DBSearch(
        bin_width=params.dbsearch[RES]['bin_width'],
        bin_offset=params.dbsearch[RES]['bin_offset'],
        tolarence_window=params.spectrum[DATASET]['tolarence_window'],
        tolarence_type=params.spectrum[DATASET]['tolarence_type'],
        remove_precursor_peak=params.dbsearch['remove_precursor_peak'],
        enzyme=params.spectrum[DATASET]['enzyme'],
        missed_cleavages=params.spectrum[DATASET]['missed_cleavages'],
        decoy_format=params.dbsearch['decoy_format'],
        max_theo_pept_peak_charge=params.dbsearch['max_theo_pept_peak_charge'],
        min_pept_len=params.dbsearch['min_pept_len'],        
        max_mods=params.spectrum[DATASET]['max_mods'],
        modifications=params.spectrum[DATASET]['modifications'],
        max_peak=params.dbsearch['max_peak'],
        static_mods=params.spectrum[DATASET]['static_mods'])
    dbs.max_intensity = params.dbsearch['max_intensity']

    # Define RBM model class object
    rbm_model = models.RBMModel(dbs)
    # Creating opimizer for PyTorch
    optimizer = optim.SGD((rbm_model.rbm_weights, rbm_model.rbm_vb, rbm_model.rbm_hb),
        lr=training_parameters['learning_rate'], weight_decay=0,
        momentum=0.9, nesterov=True)
    
    # Template to store results
    learning_curves = pd.DataFrame(columns=['Cost','psms_1','psms_5','psms_10'])

    def append_learninig_progress(pddf, cost, accepted_PSMs):
        """Updating list of cost, PSMs at 1%, 5% and 10% FDR for different methods"""
        return pddf.append({'Cost':cost, 
                            'psms_1':accepted_PSMs[0], 
                            'psms_5':accepted_PSMs[1], 
                            'psms_10':accepted_PSMs[2]}, 
                            ignore_index=True)

    # Load experimental spectra and generate peptides
    print('Loading data...', end='', flush=True)
    start_time = time.time()
    utils.load_peptides_and_spectra(dbs=dbs, 
                            fasta=params.spectrum[DATASET]['fasta_file'],
                            spectra=params.spectrum[DATASET]['ms_data_file'], 
                            data_type=DATASET)
    
    utils.generate_unique_theoretical_peaks(dbs)
    print("Done. Time: {} sec.".format(round(time.time()-start_time, 2)))

    print('Number of spectra: {}'.format(len(dbs.spectrum_collection)))
    print('Number of proteins: {}'.format(len(dbs.protein_collection)))
    print('Number of peptides: {}'.format(len(dbs.peptide_collection)))

    # Run baseline
    # Run BoltzMatch training
    batch_size = training_parameters['batch_size']
    print("Training BoltzMatch parameters...")
    start_time = time.time()
    best_psms = 0 # empty variable for best PSMs calculation

    for epoch in range(training_parameters['epochs']):
        epoch_time = time.time()
        batch_cnt = 0
        epoch_loglikelihood = 0
        spectra_processed = 0
        # Each epoch needs to run weak supervision in order to calculate q-values and create candidate peptides for spectra
        utils.spectrum_normalization(dbs)
        utils.run_tailor(dbs, filename=filename+'-tailor_baseline.dbs')

        for batch_idx in dbs.spectrum_batch_generator(batch_size, shuffle=True):
            batch_cnt += 1
            chosen_indices1, chosen_indices2 = utils.choose_indices(batch_size)
            # Prepare batch data for training
            experimental_spectra, candidate_spectra, theoretical_spectra = \
                                utils.prepare_batch_data(dbs, 
                                                         batch_idx=batch_idx, 
                                                         cand_from_search_qval=training_parameters['cand_from_search_qval'],
                                                         topN_exp=training_parameters['topN_exp'])
            # Process exceptions with empty experimental spectra
            if experimental_spectra is None:
                continue
            # Train RBM here
            # Scoring with the current model
            boltzmatch_scores = rbm_model.scoring(experimental_spectra=experimental_spectra, 
                                                  candidate_spectra=candidate_spectra,
                                                  theoretical_spectra=theoretical_spectra)
            # Diversifying regularization summand
            div_reg = utils.calculate_div_reg(rbm_model=rbm_model, experimental_spectra=experimental_spectra,
                                              chosen_indices1=chosen_indices1, chosen_indices2=chosen_indices2,
                                              div_reg_rate=1000)
            # Build the cost function
            batch_loglikelihood = torch.zeros((1), dtype=float_type).to(torch_device)
            for scores in boltzmatch_scores: 
                if len(scores) == 0:
                    continue
                scores = log_softmax(scores)
                batch_loglikelihood += -scores[0]
                spectra_processed += 1
            batch_loglikelihood += div_reg

            epoch_loglikelihood += batch_loglikelihood.data.cpu().numpy()[0]
            error = (epoch_loglikelihood/spectra_processed)

            if np.isnan(epoch_loglikelihood):
                print("Nan occured at batch: {}".format(batch_cnt+1))
                break

            # Update weights
            rbm_model.weight_update(loss=batch_loglikelihood, optimizer=optimizer)

            # Print information
            if ((batch_cnt) % batch_printing_tick == 0):
                print("Epoch: {}/{}. Batch: {} Epoch time: {} s, partial loss:{}".format(epoch+1, training_parameters['epochs'], batch_cnt, round(time.time()-epoch_time, 2), error))

        print('Diversifying regularization: {}'.format(div_reg.data.cpu().numpy()))
        print("Epoch finished: {}/{}. Total batches: {}, Total spectra: {}, Total time: {} s, loss:{}".format(epoch+1, training_parameters['epochs'], batch_cnt,spectra_processed,
                                                                                                              round(time.time()-start_time, 2), error))
        
        # Test BoltzMatch Scoring
        utils.run_boltzmatch_search(rbm_model, dbs)
        PSMs = dbs.print_accepted_psms()
        learning_curves = append_learninig_progress(pddf=learning_curves, cost=error, accepted_PSMs=PSMs)

        utils.save_learning_curve(learning_curves=learning_curves, filename_stem=filename)
        if best_psms < np.sum(PSMs): # use the sum of PSMs at 1%, 5% and 10% FDR as a criterion for choosing the best fitted model
            best_psms = np.sum(PSMs) # update best PSMs counter
            print('best_psms {} at epoch: {}'.format(best_psms, epoch+1))
            ms2file = output_filename = '.'.join(params.spectrum[DATASET]['ms_data_file'].split('.')[:-1]) + '_boltzmatch_'+RES+'.ms2'
            dbs.export_spectra_ms2(ms2file) # export transformed spectra
            dbs.print_results(filename+'-best.dbs')
            torch.save(rbm_model.state_dict(), filename+'-best.model') # saving the best model

        if PSMs[0] == 0:
            print('Overfit. Zero psms accepted at 1\% FDR :/')
            break

    print("Training Finished with best psms {}. Total training time: {} sec".format(best_psms, round(time.time()-start_time, 2)))
    return filename

# Run the whole code
if __name__ == "__main__":
    # Use arguments passing to command line in Linux OS
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dataset', type=str, default='ups1', help='The dataset to be used for training BoltzMatch matrix. Possible variants: aurum, iprg, humvar, hspp2a, ups1, yeast, malaria. Default: ups1')
    parser.add_argument('-r', '--resolution', type=str, default='lowres', help='Resolution of the data. Possible variants: "lowres" (for low resolution) and "highres" (for high resolution). Default: lowres')
    args = parser.parse_args()
    # Fit the BoltzMatch parameters and save best results
    filename = fit(args.dataset, args.resolution) # define filename stem for further possible plottings
    # Optional results plotting
    if args.resolution is 'lowres':
        utils.plot_weights(filename + "-best.model", filename, scale=0.003)
    