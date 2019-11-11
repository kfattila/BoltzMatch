#!/bin/bash
# Define resolution
resolution='lowres'
# Run all available datasets. These take around 1-2 days.
python ./rbm_train_main.py -d aurum -r $resolution    # aurum
python ./rbm_train_main.py -d iprg -r $resolution     # iprg
python ./rbm_train_main.py -d hspp2a -r $resolution   # hspp2a
python ./rbm_train_main.py -d humvar -r $resolution  # humvar
python ./rbm_train_main.py -d malaria -r $resolution  # malaria
python ./rbm_train_main.py -d yeast -r $resolution    # yeast

resolution='highres'
# Run BoltzMatch using high-res settings. These tests require large amount of memory.
# It might be needed to run on CPU instead of GPU and specify 'torch_device': torch.device("cpu") in the source/parameters.py file. These take around 2 days with CPU
python ./rbm_train_main.py -d iprg -r $resolution    # iprg
python ./rbm_train_main.py -d humvar -r $resolution  # humvar
python ./rbm_train_main.py -d malaria -r $resolution # malaria

# Run Tide-Search/CRUX on the new datasets.

./scripts/TIDE_INDEX.sh

./scripts/TIDE_SEARCH_HUMVAR.sh
./scripts/TIDE_SEARCH_AURUM.sh
./scripts/TIDE_SEARCH_HSPP2A.sh
./scripts/TIDE_SEARCH_IPRG.sh
./scripts/TIDE_SEARCH_MALARIA.sh
./scripts/TIDE_SEARCH_YEAST.sh


# Run gold standard methods

./scripts/RUN_MSAMANDA.sh
./scripts/RUN_MSGFPlus.sh
./scripts/RUN_OMSSA.sh
./scripts/RUN_XTandem.sh



