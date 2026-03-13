# Description
"""
Simulate the Erasure probability after the inner code decoding, and compare with the theoretical values.
"""
from Model.Model import *
from Model.config import DEFAULT_PASSER

from Encode.Helper_Functions import bin_to_dna, int_to_binary_array, inner_redundancy
import sys
from pathlib import Path

_PY_DIR = Path(__file__).resolve().parent
_BUILD_LIB = _PY_DIR / 'build' / f'lib.win-amd64-cpython-{sys.version_info.major}{sys.version_info.minor}'
if _BUILD_LIB.exists():
    sys.path.insert(0, str(_BUILD_LIB))

import fftqspa

import matplotlib.pyplot as plt
import csv

import numpy as np
import logging
import random

import warnings
import gc
import os
import sys
warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.CRITICAL)
plt.rcParams['figure.dpi'] = 300
np.set_printoptions(threshold=np.inf)

def DNAChannel(CodeWrdsTx, Pe, sequencingDepth, innerRedundancy):
    """
    The concatenation of inner code and the DNA storage channel simulation.
    Input:
    - CodeWrdsTx: The transmitted codewords after outer code encoding. Shape: (n_0, k_2)
    - Pe: The base error probability of the DNA storage channel.
    - sequencingDepth: The sequencing depth for the DNA storage channel.
    - innerRedundancy: The total redundancy allocated for the inner code.
    Output:
    - v_score: The voting scores after inner code decoding. Shape: (n_0, k_2)
                for a loss sequence, the voting score for each bit is set to 0.5
    """
    loss_sequence_num = 0  # Number of lost sequences during the DNA storage channel
    ##=============Default Channel parameter=======
    arg = DEFAULT_PASSER
    arg.syn_number = 30  # Synthesis number
    arg.syn_yield = 1
    arg.syn_sub_prob = 0.57 * Pe / 6
    arg.syn_ins_prob = 0
    arg.syn_del_prob = 0

    arg.decay_er = 0
    arg.decay_loss_rate = 0.43 * Pe

    arg.pcrc = 2  # PCR cycle number
    arg.pcrp = 0.8  # PCR efficiency
    arg.pcrBias = 0.05

    ps_seq = 0.57 * Pe / 6  # 测序阶段单向替换概率
    arg.seq_TM = genTm(ps_seq)
    arg.seq_depth = sequencingDepth  # Sequencing Depth
    N_c = np.ceil((2 * arg.pcrp) ** (arg.pcrc))  # The copy times during PCR

    ##===============Coding argument======================
    # BCH code parameter
    k_1 = 14  # length of index
    k_2 = 320  # length of data
    K = [k_1, k_2]
    r_i = innerRedundancy  # Total redundancy for inner code
    # r_1, r_2 = inner_redundancy(Pe, arg.seq_depth, arg.syn_number, N_c, r_i, K)
    r_1, r_2 = 15, 99  # Pre-defined redundancy allocation
    print(f"Redundancy allocation for the inner code: ({r_1}, {r_2}).")
    # Parameters for BCH code lengths
    n_1 = k_1 + r_1  # Codeword length for index BCH code
    n_2 = k_2 + r_2  # Codeword length for data BCH code
    dna_length = int(np.ceil((n_1 + n_2) / 2))
    # LDPC code parameter
    n_0 = len(CodeWrdsTx)  # Outter code length

    ##================Generating index===================
    ids = np.empty((n_0, k_1), dtype=np.uint8)
    for idx in range(n_0):
        ids[idx] = int_to_binary_array(idx, k_1)

    ##================Two layer BCH encode in C++ module=================
    print('BCH encoding...')
    cwr_ids = fftqspa.bch_encode(n_1, k_1, ids.astype(np.int32))
    CodeWrdsTx = np.array(CodeWrdsTx, dtype=np.uint8)
    cwr_data = fftqspa.bch_encode(n_2, k_2, CodeWrdsTx.astype(np.int32))
    cwr2 = np.concatenate((cwr_ids, cwr_data), axis=1)
    del ids, CodeWrdsTx, cwr_ids, cwr_data
    gc.collect()

    ##===================Convert binary data to DNA=======================
    print('Binary to DNA...')
    in_dnas = np.empty(n_0, dtype=f'<U{dna_length}')
    total_bits = n_1 + n_2
    for idx in range(n_0):
        cw = cwr2[idx]
        bit_str = ''.join(cw.astype(int).astype(str))[:total_bits]
        in_dnas[idx] = bin_to_dna(bit_str)
    del cwr2
    gc.collect()

    ##======================DNA Storage Channel Simulation===========================
    print(f'>-------------------------------<')
    print(f'Channel: Pe = {Pe}, d_seq = {arg.seq_depth}.')
    print('Synthesis is proceeding...')
    SYN = Synthesizer(arg)
    dnas_syn = SYN(in_dnas)
    random.shuffle(dnas_syn)
    print('Decay is proceeding...')
    DEC = Decayer(arg)
    dnas_dec = DEC(dnas_syn)
    print('PCR is proceeding...')
    PCR = PCRer(arg=arg)
    dnas_pcr = PCR(dnas_dec)
    print('Sequencing is proceeding...')
    SEQ = Sequencer(arg)
    dnas_seq = SEQ(dnas_pcr)
    del in_dnas, dnas_syn, dnas_dec, dnas_pcr
    gc.collect()

    dnas_sim_result = []
    for dna_set in dnas_seq:
        for dna_error_profile in dna_set['re']:
            for i in range(dna_error_profile[0]):
                dnas_sim_result.append(dna_error_profile[2])

    normized_dnas_result = []
    for i, dna in enumerate(dnas_sim_result):
        if len(dna) != dna_length:
            if len(dna) < dna_length:
                raise Warning(f"Received DNA sequence {i} is shorter than expected length {dna_length}. Padding with 'A's.")
                dna = dna.ljust(dna_length, 'A')
            else:
                dna = dna[:dna_length]
        normized_dnas_result.append(dna)
    print(f'>-------------------------------<')

    ##==================DNA to Binary===================
    print('DNA to Binary...')
    QUANT2BIN = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    index_bin_arrays = []
    inf_bin_arrays = []
    for i, dna in enumerate(normized_dnas_result):
        bin_str = ''
        if ((n_1 + n_2) % 2 == 0):
            for b in dna:
                bin_str += QUANT2BIN[b]
        else:
            for j, b in enumerate(dna):
                if j != len(dna) - 1:
                    bin_str += QUANT2BIN[b]
                else:
                    bin_str += QUANT2BIN[b][1]
        indices_bin_str = bin_str[:n_1]
        inf_bin_str = bin_str[n_1:n_1+n_2]
        index_bin_array = np.array([float(s) for s in indices_bin_str])
        inf_bin_array = np.array([float(s) for s in inf_bin_str])
        index_bin_arrays.append(index_bin_array)
        inf_bin_arrays.append(inf_bin_array)
    simu_indices_arr = np.array(index_bin_arrays)
    simu_inf_arr = np.array(inf_bin_arrays)

    ##===============Two-layer BCH decode + voting=============================
    print('BCH decoding...')
    v_score = fftqspa.bch_decode_and_vote(
        n_1,
        k_1,
        n_2,
        k_2,
        n_0,
        simu_indices_arr.astype(np.int32),
        simu_inf_arr.astype(np.int32),
    )
    loss_sequence_num = int(np.sum(np.all(v_score == 0.5, axis=1)))
    del simu_indices_arr, simu_inf_arr
    gc.collect()
    return v_score

