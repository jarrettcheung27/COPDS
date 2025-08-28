import streamlit as st
import plotly.io as pio

from Model.Model import * 
import plotly.express as px
from Helper_Functions import *
pio.templates.default = "plotly_white"
import Model.config as config
import os
import streamlit as st
import subprocess
import numpy as np
import json
from PIL import Image
import uuid
import pdb # python debugger
# ignore warnings
import warnings
warnings.filterwarnings("ignore")

CONFIG_PATH = "./config/config.json"
DNA_LIB_DIR = "./DNA_Library/"
DEBUG_MODE = True

if DEBUG_MODE:
    coding_config = json.load(open(CONFIG_PATH, "r"))
    file_ids = list(coding_config['files'].keys())
# =========================== assigning parameters ========================= #
print('Running Webapp...')
# ------------ Choosing parameters -------------- # 
st.sidebar.subheader('Coding Parameter')
chunk_size = st.sidebar.number_input('DNA chunk size', min_value = 100, max_value = 400, value = 320)
k_0 = st.sidebar.number_input('Block size $k_0$', min_value = 200, max_value = 8000, value = 1000)
n_0 = st.sidebar.number_input('Outer code length $n_0$', min_value = 200, max_value = 10000, value = 1334)
# k_0 = 1000 # Default value for k_0, LDPC code parameter
k_1 = st.sidebar.number_input('Index length $k_1$', min_value = 1, max_value = 24, value = 11)
r_1 = st.sidebar.number_input('Redundancy $r_1$ of the first layer of TL-BCH', min_value = 5, max_value = 50, value = 15)
k_2 = chunk_size
r_2 = st.sidebar.number_input('Redundancy $r_2$ of the second layer of TL-BCH', min_value = 5, max_value = 200, value =99)

# --------------------Channel parameter----------------
arg = config.DEFAULT_PASSER
st.sidebar.subheader('Parameters of DNA data storage channel')
# synthesis stage
arg.syn_number = st.sidebar.slider('Syn number', min_value = 10, max_value = 50, value = 50)
# arg.syn_sub_prob = st.sidebar.number_input('Syn Error rate', min_value = 0.0, max_value = 0.5, value = 0.00001) / 3 # 3 kinds of errors
arg.syn_yield = st.sidebar.slider('Syn Yield', min_value = 0.98, max_value = 0.995, value = 0.99)

# PCR stage
arg.pcrc = st.sidebar.slider('PCR cycle',min_value = 0, max_value =20,value =2)
arg.pcrp = st.sidebar.number_input('PCR prob',min_value = 0.5, max_value = 1.0,value = 0.9)

# decay stage
arg.decay_er = 0
arg.decay_loss_rate = 0.01

# sequencing stage
seq_platform = st.sidebar.selectbox('Sequencing Platform',['Illumina Sequencing','Nanopore'])
if seq_platform == 'Illumina Sequencing':
    arg.seq_TM = config.TM_NGS
else:
    arg.seq_TM = config.TM_NNP
arg.seq_TM = genTm(0.001) # sequencing Transform Matrix
arg.sam_ratio = st.sidebar.number_input('Sampling ratio',min_value = 0.0, max_value =1.0,value = 0.01)
arg.seq_depth = st.sidebar.slider('Seq Depth', min_value = 1, max_value = 100, value = 10)

# inspect index
index = st.sidebar.slider('inspect index', max_value = 600, value = 0)

LDPC = LDPC_Codec(coding_config)
BCH = BCH_Codec(CONFIG_PATH)  # Initialize BCH codec with config path

#===============================Binary to DNA===============================#
st.subheader('Binary to DNA')
print('Binary to DNA...')
total_bits = chunk_size + k_1 + r_1 + r_2
dna_length = int(np.ceil(total_bits / 2))  # Each DNA nucleotide encodes 2 bits
if (total_bits % 2) != 0:
    is_padding = True
else:
    is_padding = False
coding_config['DNA2Binary'] = {
    "is_padding": is_padding,
    "padding": 0,
}
DNA_Library = []
# Read BCH codeword form the txt file and convert to DNA
for i in range(coding_config['file_num']):
    pool = []
    for block_id in range(coding_config['files'][file_ids[i]]['segmentation']['block_num']):
        input_codeword_name = f"{file_ids[i]}_BCH_encode_out_{block_id}.txt"
        input_codeword_path = os.path.join(LDPC.mid_data_folder, input_codeword_name)        
        block = []
        with open(input_codeword_path, "r") as f:  # bit is separated by ','
            bch_codeword = f.read().strip()
            lines = bch_codeword.split('\n')
            for line in lines:
                block.append([int(x) if x in ['0', '1'] else ValueError("Invalid character in line") for x in line.split(',')])
        pool.append(np.array(block).T)
    dna_pool_filename = os.path.join(DNA_LIB_DIR, f"{file_ids[i]}_in.dna")
    temp = binary_to_dna_pools(pool, dna_length, is_padding, dna_pool_filename)
    DNA_Library.append(temp)
print('Binary to DNA conversion completed.')
print('DNA pools saved to:', DNA_LIB_DIR)

