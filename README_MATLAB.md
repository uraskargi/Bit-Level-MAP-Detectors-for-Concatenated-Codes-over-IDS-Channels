# LDPC (Low-Density Parity-Check) Code Simulation

This repository contains MATLAB code for simulating a Low-Density Parity-Check (LDPC) code using the Sum-Product Algorithm (SPA) decoder. The simulation involves various types of error models such as insertion, deletion, and substitution errors, with parameters set for simulation flexibility. The goal of this simulation is to evaluate the Bit Error Rate (BER) and Frame Error Rate (FER) for different combinations of insertion, deletion, and substitution probabilities.

## Table of Contents
- [Overview](#overview)
- [Dependencies](#dependencies)
- [Parameters](#parameters)
- [How to Run](#how-to-run)
- [Simulation Flow](#simulation-flow)
- [Output Results](#output-results)

## Overview
This MATLAB code simulates the performance of LDPC codes under the influence of different types of errors. Specifically, it incorporates:
- **Insertion errors**: Adding extra bits to the transmitted message.
- **Deletion errors**: Removing bits from the transmitted message.
- **Substitution errors**: Flipping the bits in the transmitted message.

The LDPC code is decoded using the Sum-Product Algorithm (SPA) and is evaluated for different error scenarios with configurable parameters.

## Dependencies
- MATLAB R2020b or later.
- The required LDPC matrix files should be present:
  - `Matrices/H_new.mat`: The parity-check matrix for the LDPC code.
  
To download these files, visit [David Mackay's LDPC code page](https://www.inference.org.uk/mackay/CodesFiles.html).

## Parameters
The following parameters are configurable:
- **iter_LDPC**: Number of iterations for the LDPC decoder.
- **early_finish**: Flag to enable early termination of the decoding when a valid codeword is found.
- **iter_block**: Number of iterations for the block-level decoder.
- **seed**: Random seed for reproducibility of results.
- **Nc**: Length of each block.
- **marker_code**: Array that defines the marker sequence for the outer encoder.
- **Pi_test**: Array of insertion error probabilities for testing.
- **Pd_test**: Array of deletion error probabilities for testing.
- **Ps_test**: Array of substitution error probabilities for testing.
- **max_fer**: Maximum number of frame errors before terminating the simulation.
- **max_codes**: Maximum number of codewords simulated.
- **start_with_marker**: Boolean flag to start the simulation with a marker.
- **estimation**: Whether the channel probabilities are known (0 for known, 1 for unknown).
- **print_every**: Print the results for every `print_every` codewords.

## How to Run
1. Ensure you have MATLAB installed (version R2020b or later recommended).
2. Download the required LDPC matrix file (e.g., `H_matrix_large.mat` or `H_matrix_small.mat`) and save it in the `Matrices/` directory.
3. Run the main script. The simulation will start and print progress at regular intervals. You can modify the parameters in the script to customize the simulation.
   - Example for running the simulation:

   ```matlab
   load('Matrices/H_new.mat');  % Load the LDPC parity-check matrix
4. The code will simulate the LDPC decoder for different error types (insertion, deletion, and substitution errors) and print out the Bit Error Rate (BER) and Frame Error Rate (FER) at regular intervals.

## Simulation Flow

1. Initial Setup: The simulation begins by initializing parameters such as error probabilities, message bits, and simulation settings.
2. Simulation Loop: For each test point (combination of insertion, deletion, and substitution probabilities), the following steps are executed:
   - A random message is generated.
   - The message is encoded using the LDPC encoder with marker bits.
   - The encoded message is sent through a channel with insertion, deletion, and substitution errors.
   - The received message is decoded using the Sum-Product Algorithm.
   - Bit Error Rate (BER) and Frame Error Rate (FER) are calculated and updated.
 3. Results Output: The simulation prints the BER and FER for each test point and the overall performance after completing the specified number of simulations.

## Output Results

The simulation outputs two key performance metrics: **Bit Error Rate (BER)** and **Frame Error Rate (FER)**. These metrics are printed during the simulation process and at the end of the simulation for each test point.

### Key Metrics:
- **Bit Error Rate (BER)**: This is the average number of bit errors divided by the total number of bits transmitted. It gives an indication of the overall error performance of the system.
  
  \[
  \text{BER} = \frac{\text{Number of bit errors}}{\text{Total number of transmitted bits}}
  \]

- **Frame Error Rate (FER)**: This is the percentage of frames (codewords) that are decoded incorrectly (i.e., they contain one or more bit errors). It gives an indication of the likelihood of receiving an erroneous codeword.

  \[
  \text{FER} = \frac{\text{Number of frames with errors}}{\text{Total number of frames}}
  \]

### Example Output:

During the simulation, the results will be printed at regular intervals (defined by the `print_every` parameter). Here’s an example of how the output will look:

