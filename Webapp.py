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

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

config_path = "./config/config.json"

# =========================== assigning parameters ========================= #

# ------------ Choosing parameters -------------- # 
st.sidebar.subheader('Coding Parameter')
chunk_size = st.sidebar.number_input('DNA chunk size', min_value = 100, max_value = 400, value = 320)
k_0 = st.sidebar.number_input('Block size $k_0$', min_value = 200, max_value = 2000, value = 1000)
n_0 = st.sidebar.number_input('Outer code length $n_0$', min_value = 200, max_value = 10000, value = 1334)
# k_0 = 1000 # Default value for k_0, LDPC code parameter
k_1 = st.sidebar.number_input('Index length $k_1$', min_value = 1, max_value = 24, value = 11)
r_1 = st.sidebar.number_input('Redundancy $r_1$ of the first layer of TL-BCH', min_value = 5, max_value = 50, value = 15)
k_2 = chunk_size
r_2 = st.sidebar.number_input('Redundancy $r_2$ of the second layer of TL-BCH', min_value = 5, max_value = 200, value =99)
coding_config = {"Coding_param" : {"outer": {"n": n_0, "k": k_0}, "inner1": {"n": k_1 + r_1, "k": k_1}, "inner2": {"n": k_2 + r_2, "k": k_2}}}
save_coding_config(config_path, outer_code = (n_0, k_0), inner1 = (k_1 + r_1, k_1), inner2 = (chunk_size + r_2, chunk_size)) # The same coding parameters are used for all images.

# --------------------Channel parameter----------------
arg = config.DEFAULT_PASSER
st.sidebar.subheader('Parameters of DNA data storage channel')
# synthesis stage
arg.syn_number = st.sidebar.slider('Syn number', min_value = 10, max_value = 60, value = 60)
# arg.syn_sub_prob = st.sidebar.number_input('Syn Error rate', min_value = 0.0, max_value = 0.5, value = 0.00001) / 3 # 3 kinds of errors
arg.syn_yield = st.sidebar.slider('Syn Yield', min_value = 0.98, max_value = 0.995, value = 0.99)

# PCR stage
arg.pcrc = st.sidebar.slider('PCR cycle',min_value = 0, max_value =20,value =6)
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
# ------------- file path--------------- #
abs_dir  = "D:/COPDS-main/"
dna_lib_dir = abs_dir + 'DNA_Library/'
in_file_path = "D:/COPDS-main/DNA_Library/Jnu_test.jpg"
in_file_name= "Jnu_test.jpg"
file_name = "Jnu_test"
suffix = "jpg"
out_file_path = "D:/COPDS-main/DNA_Library/" + file_name + "_re." + suffix
dna_pools_dir = os.path.join(dna_lib_dir, 'DNA_pools-' + file_name)


print('\n ==================== DNA 存储仿真平台 ==================== ')
# ======================Select image to store ====================== #
# Create the DNA_Library directory if it doesn't exist
os.makedirs(dna_lib_dir, exist_ok=True)

# =====================Segmentation==================== #
st.subheader('Segmentation')

data,pad = preprocess(in_file_path,int(chunk_size/8))
'Images loaded and split into ', len(data), ' data chunks.'

# =====================Binary to bits===================#
# bytes to bits
data = [bytes2bits(chunk) for chunk in data]

# ========================Splitting=======================#
# Split data into pools of size k_0 (1000), pad the last pool if needed

data_pools = split_data(data, pool_size=k_0, chunk_size=chunk_size)
pools_num = len(data_pools)

save_config_segmentation(in_file_name, pools_num, len(data), pad)
'Data is divided into', len(data_pools), ' pools of size ', len(data_pools[0]), ' chunks.'

# Save input message for debugging
for pool_idx, pool in enumerate(data_pools):
    message_input_file = f"D:\DeSP-main\App_Simulation_Platform\Mid_data\ldpc_encode_input_{pool_idx}.txt"
    with open(message_input_file, "w") as f:
        for chunk in pool:
            f.write(chunk + "\n")

# =================================Encoding===================================#
# LDPC encode
print('LDPC encoding...')
LDPC = LDPC_Codec()
# Transpose matrix for inter-oligos encoding
data_pools = transpose_h2v(data_pools)
'Data is transposed for inter-oligos encoding, each pool is a list of strings, each string is the i-th bit of all chunks.'

# Path to the LDPC encoder executable and H matrix
ldpc_encoder_exe = ".\Encode\LDPC_PEG-v2.0.exe"
h_matrix_path = ".\config\Hmatrix.txt"

data_pools = LDPC.encode(data_pools)

# Convert each pool to a 2-D numpy array and transpose for BCH encoding
data_pools = transpose_v2h(data_pools)

# -----------------------Indexing-----------------------
ids  = np.array([int_to_binary_array(id, k_1) for id in range(n_0)])

# ----------------------BCH encode----------------------
encode_config = {}
encode_config['n1'] = k_1 + r_1
encode_config['k1'] = k_1
encode_config['n2'] = chunk_size + r_2
encode_config['k2'] = chunk_size
encode_config['n0'] = n_0
BCH = BCH_Codec(encode_config = encode_config)
data_pools = BCH.encode(ids, data_pools)

#===============================Binary to DNA===============================#
st.subheader('Binary to DNA')
total_bits = chunk_size + k_1 + r_1 + r_2
dna_length = int(np.ceil(total_bits / 2))  # Each DNA nucleotide encodes 2 bits
if (total_bits % 2) != 0:
    is_padding = True
else:
    is_padding = False
save_padding_info_bin2DNA(in_file_name, is_padding, 0) # Save padding info for binary to DNA conversion

dna_pools = binary_to_dna_pools(data_pools, dna_length, is_padding, dna_pools_dir)

# ============================= DNA simulation Channel ============================= #
st.header('Error simulation of the DNA data storage channel')

st.subheader('Load Data')
Channel = DNA_Channel_Model(Modules = 0, arg = arg) # No model provided

# Run Simulation for each pool
out_dna_pools = []
for dnas in dna_pools:
    out_dnas = Channel(dnas)
    out_dna_pools.append(out_dnas)
'Simulation completed. '


# --------------Save Simulation results ----------------
pools_out_file_path = os.path.join(dna_pools_dir, 'out.dna')
temp = []
for idx, out_dnas in enumerate(out_dna_pools):
    # Extract simulated DNA sequences
    dnas = extract_dnas(out_dnas)
    # Open the file in write mode
    with open(pools_out_file_path, "w") as file:
        # Write DNA sequences as JSON: key is pool index, value is list of DNA sequences
        json.dump({f"pool{idx}": dnas}, file, ensure_ascii=False, indent=2)
    temp.append(dnas)
out_dna_pools = temp
print('Simulated DNAs saved to ', pools_out_file_path)


# ===================== Read the simulated DNA ===================== #
'''
st.header('Decoding')
st.subheader('Load Data')
# Folder selection for DNA pools

uploaded_dna_files = st.file_uploader(
    "Select one or more simulated DNA pool to decode", 
    type=["dna"], 
    accept_multiple_files=True
)
# Obtain the original image file name from the selected .dna file
file_name = ''
if uploaded_dna_files:
    # Use the first uploaded file as reference
    dna_file_name = uploaded_dna_files[0].name
    # Remove 'simu_' prefix if present and split at the last underscore
    base_name = dna_file_name
    if base_name.startswith('simu_'):
        base_name = base_name[5:]
    file_name = base_name.rsplit('_', 1)[0]
out_file_name = file_name + '_reconstructed.jpg'    

if uploaded_dna_files:
    out_dnas_pools = []
    strands_num = 0;
    for uploaded_file in uploaded_dna_files:
        dnas = uploaded_file.read().decode("utf-8").splitlines()
        out_dnas = [dna.strip() for dna in dnas]
        strands_num += len(out_dnas)
        out_dnas_pools.append(out_dnas)
    st.success(f"Loaded {strands_num} strands of length {len(out_dnas[0])} nts")
else:
    st.stop()
    st.warning("No DNA files selected.")
'''

config_path = "./config/config.json"
with open(config_path, 'r', encoding='utf-8') as f:
    coding_config = json.load(f)

#================= DNA to Binary ====================#
st.subheader('DNA to Binary')
out_data_pools = DNA_pools_to_binary(coding_config,in_file_name,out_dna_pools)

# ============================== BCH Decoding ============================= #

BCH = BCH_Codec(decode_config = coding_config)
out_data_pools = BCH.decode(out_data_pools)
print('voting...')
# Cluster the sequences with the same index in each pool
out_prob_pools = [] # Store the voting results for each pool
for pool in out_data_pools:
    voting_result = voting(pool, coding_config["Coding_param"]["outer"]["n"], coding_config["Coding_param"]["inner1"]["k"], coding_config["Coding_param"]["inner2"]["k"])
    out_prob_pools.append(voting_result) # Probability of each bit in the voting result, for LDPC decoding
print(f'Complete!')


# ===================== LDPC Decoding ===================== #
LDPC = LDPC_Codec()
LDPC.decode(out_prob_pools)

# ===================== testing =====================
# calculate the BER
message_input_file = f"./Mid_data/ldpc_encode_input_0.txt"
message_output_file = f"./Mid_data/ldpc_decode_output_0.txt"

# Calculate Bit Error Rate (BER)
def calculate_ber(message_input_file, message_output_file):
    with open(message_input_file, "r") as fin, open(message_output_file, "r") as fout:
        input_lines = [line.strip() for line in fin]
        output_lines = [line.strip() for line in fout]
    total_bits = 0
    error_bits = 0
    for in_bits, out_bits in zip(input_lines, output_lines):
        total_bits += len(in_bits)
        error_bits += sum(b1 != b2 for b1, b2 in zip(in_bits, out_bits))
    ber = error_bits / total_bits if total_bits > 0 else 0
    return ber

ber = calculate_ber(message_input_file, message_output_file)
print(f"Bit Error Rate (BER): {ber:.8f}")


st.header('Image Reconstruction')
# Read decoded output to reconstruct the image

# Read the config file to obtain the configuration.
pool_num = coding_config[in_file_name]["segmentation"]["pools_num"]
chunk_num = coding_config[in_file_name]["segmentation"]["chunk_num"]
pad = coding_config[in_file_name]["segmentation"]["pad"]
decoded_chunks = []
for pool_idx in range(pool_num):
    ldpc_output_file = f"./Mid_data/ldpc_decode_output_{pool_idx}.txt"
    with open(ldpc_output_file, "r") as f:
        pool_bits = [line.strip() for line in f]
        # Transpose: each line is a bit position, columns are chunks
        pool_bits_T = list(zip(*pool_bits))
        for chunk_bits_tuple in pool_bits_T:
            chunk_bits = ''.join(chunk_bits_tuple)
            chunk_bytes = bits2bytes(chunk_bits)
            decoded_chunks.append(chunk_bytes)

# Remove padding chunks if pad decoded_chunks number > original chunk number
decoded_chunks = decoded_chunks[:chunk_num]

# Remove padding bytes in the last chunk if pad > 0
if pad > 0:
    last_chunk = decoded_chunks[-1]
    if len(last_chunk) > pad:
        decoded_chunks[-1] = last_chunk[:-pad]  # Remove padding bytes from the last chunk 

# Save reconstructed image
with open(out_file_path, "wb") as f:
    for chunk in decoded_chunks:
        f.write(chunk)

st.success(f"Image reconstructed and saved as {out_file_path}")

st.subheader('Quality Evaluation')
st.image(out_file_path, width = 300)

# ------------------------ optimizing ----------------------------- #
'''
st.header('Encoding Optimization')
'The encoding-simulation-decoding process should have been run several times to estimate the distributions.'
'To save computation here, we only run the process once and use the observed values and pre-computed values to perform the estimation. Some bias might be introduced because of this.'

loc = 500
scale = 7.5
N = len(data)
Ld = len(data[0])
FA = FT_Analyzer_Simplified(N,Ld,alpha,loc,scale,dnas_seq)

st.subheader('Choosing RS length: ')
fig, rs = FA.choose_rs()
st.write(fig)
'According to D(k), rs can be set to ', rs

st.subheader('Choosing Alpha ')
fig = FA.choose_alpha()
st.write(fig)
'Alpha can be selected from the second graph to meet a specific success possibility requirement.'
'''