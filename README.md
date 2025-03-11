# Baseline LDPC Decoder Simulation

This repository contains a Python script to simulate the baseline LDPC (Low-Density Parity-Check) decoder with a BCJR-based estimator for concatenated codes over IDS channels. The script evaluates the performance of the decoder for various deletion, insertion, and substitution probability test points and records the results for further analysis.

## Requirements

- Python 3.x
- Required libraries:
  - `numpy`
  - `scipy`
  - `argparse`
  
You can install the required libraries by running:

pip install numpy scipy

## Usage

python baseline_decoder_simulation.py --main-dir ./Results --test_points_Pd [0.1] --test_points_Ps [sys.float_info.min] --test_points_Pi [sys.float_info.min] --fer_errors 1 --marker_sequence [0,1,0] --Nc 10 --iter_num 100 --matrix small --print_every 10 --CSI unk --seed 1000

## Example Usage

python baseline_decoder_simulation.py --main-dir ./Results --test_points_Pd [0.05] --test_points_Ps [sys.float_info.min] --test_points_Pi [sys.float_info.min] --fer_errors 2 --marker_sequence [1, 0, 1] --Nc 20 --iter_num 200 --matrix medium --print_every 5 --CSI unk --seed 42

## Output

./Results/
    └── matrix/
        └── Nc_<Nc>_Nm_<Nm>_k_<k>_n_<n>/
            └── test_results/
                └── test_results.npy
