import streamlit as st
import plotly.io as pio

from Model.Model import * 
import plotly.express as px
from Helper_Functions import *
pio.templates.default = "plotly_white"
import Model.config as config
import os
import random
import glob
import streamlit as stq
import subprocess
import numpy as np
import json
import matlab.engine


config_path = "./config/config.json"


# ---------------------------- Choosing parameters ---------------------------- # 

# =========================== assigning parameters ========================= #
st.sidebar.subheader('Coding Parameter')
chunk_size = st.sidebar.number_input('DNA chunk size', min_value = 100, max_value = 400, value = 320)
k_0 = st.sidebar.number_input('Block size $k_0$', min_value = 200, max_value = 2000, value = 1000)
n_0 = st.sidebar.number_input('Outer code length $n_0$', min_value = 200, max_value = 10000, value = 1334)
# k_0 = 1000 # Default value for k_0, LDPC code parameter

k_1 = st.sidebar.number_input('Index length $k_1$', min_value = 1, max_value = 24, value = 11)
r_1 = st.sidebar.number_input('Redundancy $r_1$ of the first layer of TL-BCH', min_value = 5, max_value = 50, value = 15)
r_2 = st.sidebar.number_input('Redundancy $r_2$ of the second layer of TL-BCH', min_value = 5, max_value = 200, value =99)

save_coding_config(config_path, outer_code = (n_0, k_0), inner1 = (k_1 + r_1, k_1), inner2 = (chunk_size + r_2, chunk_size)) # The same coding parameters are used for all images.
# --------------------Channel parameter----------------
arg = config.DEFAULT_PASSER
st.sidebar.subheader('Parameters of DNA data storage channel')
arg.syn_number = st.sidebar.slider('Syn number', min_value = 10, max_value = 50, value = 30)
arg.syn_sub_prob = st.sidebar.number_input('Syn Error rate', min_value = 0.0, max_value = 0.5, value = 0.02) / 3 # 3 kinds of substitutions
arg.syn_yield = st.sidebar.slider('Syn Yield', min_value = 0.98, max_value = 0.995, value = 0.99)

arg.pcrc = st.sidebar.slider('PCR cycle',min_value = 0, max_value =20,value =6)
arg.pcrp = st.sidebar.number_input('PCR prob',min_value = 0.5, max_value = 1.0,value = 0.8)

arg.sam_ratio = st.sidebar.number_input('Sampling ratio',min_value = 0.0, max_value =1.0,value = 0.005)
arg.seq_depth = st.sidebar.slider('Seq Depth', min_value = 1, max_value = 100, value = 10)
seq_platform = st.sidebar.selectbox('Sequencing Platform',['Illumina Sequencing','Nanopore'])

index = st.sidebar.slider('inspect index', max_value = 600, value = 0)

if seq_platform == 'Illumina Sequencing':
    arg.seq_TM = config.TM_NGS
else:
    arg.seq_TM = config.TM_NNP

# ------------- file path--------------- #
abs_dir  = "D:/COPDS-main/"
dna_lib_dir = abs_dir + 'DNA_Library/'
in_file_path = "D:/COPDS-main/DNA_Library/Jnu_test.jpg"
in_file_name= "Jnu_test.jpg"
file_name = "Jnu_test"
suffix = "jpg"
out_file_path = "D:/COPDS-main/DNA_Library/" + file_name + "_re." + suffix


print('\n ==================== DNA 存储仿真平台 ==================== ')
# ======================Select image to store ====================== #
# Create the DNA_Library directory if it doesn't exist
os.makedirs(dna_lib_dir, exist_ok=True)

# =====================Segmentation==================== #
st.subheader('Segmentation')

data,pad = preprocess(in_file_path,int(chunk_size/8))
'Images loaded and split into ', len(data), ' data chunks.'

save_config_segmentation(in_file_name, len(data), pad)

# ==================Binary to bits==================#
# bytes to bits
data_temp = []
for chunk in data:
    data_temp.append(bytes2bits(chunk))
data = data_temp

# ==================Spliting==================#
# Split data into pools of size k_0 (1000), pad the last pool if needed

def split_data(data, pool_size):
    """
    Split data into pools of size pool_size.
    If the last pool is smaller than pool_size, pad it with random bits.
    Input:
    - data: List of bit strings (each string is a chunk of bits)
    - pool_size: Size of each pool
    Output:
    - data_pools: List of pools, each pool is a list of bit strings
    """
    data_pools = []
    for i in range(0, len(data), pool_size):
        pool = data[i:i+pool_size]
        # If last pool is smaller than pool_size, pad with random bits
        if len(pool) < pool_size:
            chunk_len = len(pool[0]) if pool else chunk_size
            for _ in range(pool_size - len(pool)):
                random_bits = ''.join(random.choice('01') for _ in range(chunk_len))
                pool.append(random_bits)
        data_pools.append(pool)
    return data_pools
data_pools = split_data(data, k_0)
'Data is divided into', len(data_pools), ' pools of size ', len(data_pools[0]), ' chunks.'

# 'Datas is devided into', len(data_pools), ' pools of size ', len(data_pools[0]), ' chunks.'    
# Prepare input file for the encoder
# for pool_idx, pool in enumerate(data_pools):
#     input_file = f"D:\DeSP-main\App_Simulation_Platform\Mid_data\ldpc_input_{pool_idx}.txt"
#     with open(input_file, "w") as f:
#         for chunk in pool:
#             f.write(chunk + "\n")


# =================================Encoding===================================#
# LDPC encode
print('LDPC encoding...')
Codec = Codec(h_matrix_path="./config/Hmatrix.txt")
# Transpose matrix for inter-oligos encoding
def transpose_h2v(data_pools):
    """
    Transpose the data pools from horizontal to vertical for inter-oligos encoding.
    Each pool is a list of bit strings (chunks).
    Transpose to get a list of strings, each string is the i-th bit of all chunks.
    """
    transposed_data_pools = []
    for pool in data_pools:
        # Each pool is a list of bit strings (chunks)
        # Transpose: list of strings -> list of strings, each string is the i-th bit of all chunks
        transposed_pool = [''.join(chunk[i] for chunk in pool) for i in range(len(pool[0]))]
        transposed_data_pools.append(transposed_pool)
    return transposed_data_pools
data_pools = transpose_h2v(data_pools)
'Data is transposed for inter-oligos encoding, each pool is a list of strings, each string is the i-th bit of all chunks.'

# Path to the LDPC encoder executable and H matrix
ldpc_encoder_exe = ".\Encode\LDPC_PEG-v2.0.exe"
h_matrix_path = ".\config\Hmatrix.txt"

mode = "encode"
data_pools = Codec.LDPC_encode(data_pools, ldpc_encoder_exe)

# Convert each pool to a 2-D numpy array and transpose for BCH encoding
data_pools = transpose_v2h(data_pools)

# -----------------------Indexing-----------------------
print("Generating indices...")
ids = []
n_0 = len(data_pools[0])  # Number of chunks in the first pool
def gen_dix(chunk_num_encoded):
    """
    Generate indices for the data pools.
    Each index is a binary representation of the chunk number.
    Input:
    - chunk_num_encoded: Number of chunks in the encoded data
    Output:
    - ids: 2-D numpy array of shape (n_0, k_1), where n_0 is the number of chunks and k_1 is the index length
    """
    ids = []
    k_1 = int(np.ceil(np.log2(chunk_num_encoded)))
    # Generate indices
    for id in range(n_0):
        ids.append(int_to_binary_array(id, k_1))
    ids = np.array(ids)
    return ids
ids = gen_dix(n_0)

# ----------------------BCH encode----------------------
print("BCH encoding...")
encoded_data_pools = []
matlab_engine = matlab.engine.start_matlab()
matlab_engine.cd(r'D:\DeSP-main', nargout=0)  # Set MATLAB working directory
eng = matlab_engine
# Encode index by BCH encoder.
st.text('BCH encoding...')
cwr_ids = eng.BCH_Encoder(k_1 + r_1, k_1, n_0, ids)

for pool_idx, pool in enumerate(data_pools):
    # Encode information bit by BCH encoder.
    cwr_data = eng.BCH_Encoder(chunk_size + r_2, chunk_size, n_0, pool)
    data_pool = np.concatenate((cwr_ids, cwr_data), axis=1)
    encoded_data_pools.append(data_pool)
data_pools = encoded_data_pools

#================================Binary to DNA==================================#
st.subheader('Binary to DNA')
dna_pools = []
total_bits = chunk_size + k_1 + r_1 + r_2
dna_length = int(np.ceil(total_bits / 2))  # Each DNA nucleotide encodes 2 bits
if (total_bits % 2) != 0:
    is_padding = True
else:
    is_padding = False
save_padding_info_bin2DNA(in_file_name, is_padding, 0) # Save padding info for binary to DNA conversion
for data_pool in data_pools:
    dnas = []
    for binary_codeword in data_pool:
        if is_padding: # check length of each codeword, if less than dna_length*2, pad with 0s.
            binary_codeword = np.pad(binary_codeword, (0, dna_length*2 - len(binary_codeword)), 'constant')
        bit_str = ''.join(binary_codeword.astype(int).astype(str))
        dna = bin_to_dna(bit_str)
        dnas.append(dna)
    dna_pools.append(dnas)

# --------------- Save DNAs to files ------------------ #
# Create a subfolder in DNA_Library named as "<image_name> DNA pools"
dna_pools_dir = os.path.join(dna_lib_dir, 'DNA_pools-' + file_name)
os.makedirs(dna_pools_dir, exist_ok=True)
# Save DNA pools to a .dna file, with different key for each pool
pools_in_file_name = os.path.join(dna_pools_dir, 'in.dna')
for idx, dnas in enumerate(dna_pools):
    # Open the file in write mode
    with open(pools_in_file_name, "w") as file:
        # Write DNA sequences as JSON: key is pool index, value is list of DNA sequences
        json.dump({f"pool{idx}": dnas}, file, ensure_ascii=False, indent=2)
print('DNAs saved to ', pools_in_file_name)

# ================================ DNA simulation Channel ================================ #
st.header('Error simulation of the DNA data storage channel')

st.subheader('Load Data')
Channel = DNA_Channel_Model(Modules = 0, arg = arg) # No model provided

# Run Simulation for each pool
out_dna_pools = []
for dnas in dna_pools:
    out_dnas = Channel(dnas)
    out_dna_pools.append(out_dnas)
'Simulation completed. '


# -------------------Save Simulation results -------------------
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
    config = json.load(f)

#================= DNA to Binary ====================#
st.subheader('DNA to Binary')
QUANT2BIN = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
out_binary_pools = []

out_data_pools = []

# separate index and information part
for out_dnas in out_dna_pools:
    ids = [] # Store the index bits arrays
    infs = [] # Store the information bits arrays
    temp_dnas = [] # Temporary dictionary to store DNA {id, information part, counts}.
    for dna, num in out_dnas.items():
        dna = dna.strip()
        binary = dna_to_bin_array(dna) # Convert DNA to binary array
        id = binary[:config["Coding_param"]["inner1"]["n"]] # Extract index part
        if config[in_file_name]["Bin2DNA"]["is_padding"]: # Extract information part, remove padding bit if exists
            inf = binary[config["Coding_param"]["inner1"]["n"]:-1] # Remove the padding bit
        else:
            inf = binary[config["Coding_param"]["inner1"]["n"]:]
        temp_dnas.append([id, inf, num]) # Store the id, information part and counts
    out_data_pools.append(temp_dnas)

print('BCH Decoding...')
# def BCH_Decoder(k, n, data):
#     readsnum = len(data[""]) # Number of reads in each pool
#     # BCH decode index part
#     eng = matlab.engine.start_matlab()
#     eng.cd(r'D:\DeSP-main', nargout=0)  # Set MATLAB working directory
#     cwr_ids = eng.BCH_Decoder(n, k, readsnum, data[0])
#     # BCH decode information part
#     cwr_data = eng.BCH_Decoder(chunk_size + r_2, chunk_size, readsnum, infs)
#     return cwr_ids, cwr_data

# Decode each pool
eng = matlab.engine.start_matlab()
eng.cd(r'D:\DeSP-main', nargout=0)  # Set MATLAB working directory
temp_pools = []
for pool in out_data_pools:
    ids = np.array([segment[0] for segment in pool])
    infs = np.array([segment[1] for segment in pool])
    ids = eng.BCH_Decoder(config["Coding_param"]["inner1"]["n"], config["Coding_param"]["inner1"]["k"], len(ids), ids)
    data = eng.BCH_Decoder(config["Coding_param"]["inner2"]["n"], config["Coding_param"]["inner2"]["k"], len(infs), infs)
    temp = []
    for i, id in enumerate(ids):
        temp_id = ''.join(str(int(s)) for s in id[1:]) # Convert index bits to string, remove the first bit which is the flag
        temp_data = ''.join(str(int(s)) for s in data[i][1:]) # Convert information bits to string, remove the first bit which is the flag
        temp.append((temp_id, temp_data, pool[i][2])) # Store the id, information part and counts
    temp_pools.append(temp)
out_data_pools  = temp_pools
print('BCH Decoding completed.')

print('voting...')
# Cluster the sequences with the same index in each pool
def voting(pool, n_0, index_length, chunk_size):
    '''
    将相同index的序列聚类投票
    输入：
    ids: list of index bit strings
    infs: list of information bit strings
    n_0: number of chunks in the pool
    chunk_size: size of each chunk in bits
    '''
    # prepare the segments for voting
    result = {} # Store the voting result in a dictionary

    for segment in pool:
        ids = segment[0]  # Index bits
        inf = segment[1]
        num = segment[2]  # Number of sequences with the same index
        
        if int(ids,2) <= n_0: 
            if ids not in result:
                result[ids] = {inf: num}
            else:
                if inf in result[ids]:
                    result[ids][inf] += num
                else:
                    result[ids][inf] = num
    # if some ids are missing, fill them with empty dictionaries, key are the index bits of length index_length in binary
    for i in range(n_0):
        if f"{i:0{index_length}b}" not in result:
            result[f"{i:0{index_length}b}"] = {}
    # sort result by index bits
    result = dict(sorted(result.items(), key=lambda x: int(x[0], 2)))  # Sort by index bits
    # voting
    voting_result = [] # Store the voting result
    for ids, inf_dict in result.items():
        # Sort the information bits by their counts
        sorted_infs = sorted(inf_dict.items(), key=lambda x: x[1], reverse=True)
        # Calculate the average of the most frequent information bits
        avg_inf = np.zeros(chunk_size)
        total_count = 0
        for inf, count in sorted_infs:
            avg_inf += np.array(list(map(int, inf))) * count
            total_count += count
        if total_count > 0:
            avg_inf /= total_count
        voting_result.append(avg_inf)
    # Convert the voting result to a numpy array and transpose it for LDPC decoding
    if len(voting_result) == 0:
        # If no valid segments, return an empty array
        voting_result = np.zeros((chunk_size, n_0))
    voting_result = np.array(voting_result)
    return voting_result.T

out_prob_pools = [] # Store the voting results for each pool
for pool in out_data_pools:
    voting_result = voting(pool, config["Coding_param"]["outer"]["n"], config["Coding_param"]["inner1"]["k"], config["Coding_param"]["inner2"]["k"])
    out_prob_pools.append(voting_result) # Probability of each bit in the voting result, for LDPC decoding
print('Voting completed.')

# Prepare and save the voting results for LDPC decoding
# v_score shape: (chunk_size, n_0)
# Save as text file for LDPC decoder input (each line is a bit position, values are probabilities for each chunk)
for pool_idx, v_score in enumerate(out_prob_pools):
    ldpc_input_file = f"D:\\DeSP-main\\App_Simulation_Platform\\Mid_data\\ldpc_decode_input_{pool_idx}.txt"
    np.savetxt(ldpc_input_file, v_score, fmt="%.3f", delimiter=" ")
print('Voting results saved to LDPC input files.')
# ===================== LDPC Decoding ===================== #
print('LDPC Decoding...')
for pool_idx, v_score in enumerate(out_data_pools):
    LDPC_decoder = LDPC_codec(mode="decode")
    # Call LDPC decoder executable
    ldpc_decoder_exe = "D:\\DeSP-main\\App_Simulation_Platform\\LDPC_PEG-v2.0.exe"
    h_matrix_path = "D:\\DeSP-main\\App_Simulation_Platform\\config\\Hmatrix.txt"
    ldpc_output_file = f"D:\\DeSP-main\\App_Simulation_Platform\\Mid_data\\ldpc_decode_output_{pool_idx}.txt"
    mode = "decode"
    result = subprocess.run([ldpc_decoder_exe, mode, ldpc_input_file, ldpc_output_file, h_matrix_path], capture_output=True, text=True)
    print(result.stdout)  # Print the output from the LDPC decoder
    print(result.stderr)  # Print any error messages from the LDPC decoder
print('LDPC Decoding completed.')

# calculate the BER
input_msg_file = f"D:\\DeSP-main\\App_Simulation_Platform\\Mid_data\\ldpc_input_0.txt"
output_msg_file = f"D:\\DeSP-main\\App_Simulation_Platform\\Mid_data\\ldpc_decode_output_0.txt"

# Calculate Bit Error Rate (BER)
def calculate_ber(input_file, output_file):
    with open(input_file, "r") as fin, open(output_file, "r") as fout:
        input_lines = [line.strip() for line in fin]
        output_lines = [line.strip() for line in fout]
    total_bits = 0
    error_bits = 0
    for in_bits, out_bits in zip(input_lines, output_lines):
        total_bits += len(in_bits)
        error_bits += sum(b1 != b2 for b1, b2 in zip(in_bits, out_bits))
    ber = error_bits / total_bits if total_bits > 0 else 0
    return ber

ber = calculate_ber(input_msg_file, output_msg_file)
print(f"Bit Error Rate (BER): {ber:.6f}")


st.header('Image Reconstruction')
# Read decoded output to reconstruct the image
# Read all LDPC decoded output files and reconstruct the original data
decoded_chunks = []
for pool_idx in range(len(out_prob_pools)):
    ldpc_output_file = f"D:\\DeSP-main\\App_Simulation_Platform\\Mid_data\\ldpc_decode_output_{pool_idx}.txt"
    with open(ldpc_output_file, "r") as f:
        pool_bits = [line.strip() for line in f]
        # Transpose: each line is a bit position, columns are chunks
        pool_bits_T = list(zip(*pool_bits))
        for chunk_bits_tuple in pool_bits_T:
            chunk_bits = ''.join(chunk_bits_tuple)
            chunk_bytes = bits2bytes(chunk_bits)
            decoded_chunks.append(chunk_bytes)

# Read 'pad' and 'chunk_num' from the segmentation config file
segmentation_config_path = f"D:\\DeSP-main\\App_Simulation_Platform\\Mid_data\\Segmentation_config_{file_name}.txt"
if os.path.exists(segmentation_config_path):
    with open(segmentation_config_path, "r") as f:
        for line in f:
            if line.startswith("pad="):
                pad = int(line.strip().split("=")[1])
            elif line.startswith("chunk_num="):
                chunk_num = int(line.strip().split("=")[1])
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
# st.image(out_file_path, width = 300)

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