import streamlit as st
import plotly.io as pio
import matplotlib.pyplot as plt
from Model.Model import * 
import plotly.express as px
from Analysis.Analysis import dna_chunk, plot_oligo_number_distribution, plot_error_distribution, save_simu_result
from Encode.Helper_Functions import *
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
import matplotlib.pyplot as plt

# ============================ Helper Functions ============================ #
def select_image():
    # Prompt user to select an image to store in DNA
    uploaded_file = st.file_uploader("Select images to store in DNA", type=["jpg", "png"])
    if uploaded_file is None:
        st.stop()  # Stop the process until an image is uploaded.

    # Extract file name and suffix
    file_name, suffix = uploaded_file.name.split('.')

    # Get absolute path of the uploaded file
    abs_file_path = os.path.abspath(uploaded_file.name)

    # Get absolute directory for DNA library
    abs_dir = os.path.dirname(abs_file_path)
    dna_lib_dir = os.path.join(abs_dir, 'DNA_Library')
    os.makedirs(dna_lib_dir, exist_ok=True)

    # Prepare file paths
    in_dna_name = os.path.join(dna_lib_dir, file_name + '.dna')
    out_dna_name = os.path.join(dna_lib_dir, 'simu_' + file_name + '.dna')
    out_file_name = os.path.join(dna_lib_dir, file_name + '_re.' + suffix)

    print('DNA 库路径为:', dna_lib_dir)
    print('输出文件路径为：', out_file_name)
    return abs_file_path, file_name, suffix, in_dna_name, out_dna_name, out_file_name

# ---------------------------- Choosing parameters ---------------------------- # 

def inspect(dnas, num_th = 1, inspect_index = 20):
    fig = plt.figure(figsize = (12,6))
    plt.subplot(1,2,1)
    plot_oligo_number_distribution(dnas)
    plt.subplot(1,2,2)
    plot_error_distribution(dnas,th = num_th)
    st.write(fig)
    dc = dna_chunk(dnas[inspect_index],'html')
    table = dc.plot_re_dnas()
    st.markdown(table, unsafe_allow_html = True)
    st.write(dc.plot_voting_result())

# =========================== assigning parameters ========================= #
st.sidebar.subheader('Coding Parameter')
chunk_size = st.sidebar.number_input('DNA chunk size', min_value = 100, max_value = 400, value = 320)
k_0 = st.sidebar.number_input('Block size $k_0$', min_value = 200, max_value = 2000, value = 1000)
n_0 = st.sidebar.number_input('Outer code length $n_0$', min_value = 200, max_value = 10000, value = 1334)
# k_0 = 1000 # Default value for k_0, LDPC code parameter

k_1 = st.sidebar.number_input('Index length $k_1$', min_value = 1, max_value = 24, value = 11)
r_1 = st.sidebar.number_input('Redundancy $r_1$ of the first layer of TL-BCH', min_value = 5, max_value = 50, value = 15)
r_2 = st.sidebar.number_input('Redundancy $r_2$ of the second layer of TL-BCH', min_value = 5, max_value = 200, value =99)

# Channel parameter
arg = config.DEFAULT_PASSER
st.sidebar.subheader('Parameters of DNA data storage channel')
arg.syn_number = st.sidebar.slider('Syn number', min_value = 10, max_value = 50, value = 20)
arg.syn_sub_prob = st.sidebar.number_input('Syn Error rate', min_value = 0.0, max_value = 0.5, value = 0.2) / 3 # 3 kinds of substitutions
arg.syn_yield = st.sidebar.slider('Syn Yield', min_value = 0.98, max_value = 0.995, value = 0.99)

arg.pcrc = st.sidebar.slider('PCR cycle',min_value = 0, max_value =20,value = 5)
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
print('\n ==================== DNA 存储仿真平台 ==================== ')
# ======================Select image to store ====================== #
# Create the DNA_Library directory if it doesn't exist
os.makedirs(dna_lib_dir, exist_ok=True)
# File uploader for selecting images
# in_file_path, file_name, suffix, in_dna_name, out_dna_name, out_file_name = select_image()


# =====================Segmentation==================== #
st.subheader('Segmentation')

data,pad = preprocess(in_file_path,int(chunk_size/8))
'Images loaded and split into ', len(data), ' data chunks.'

# Save segmentation configuration for later reconstruction
def save_config(in_file_name, chunk_num, pad):
    '''
    Save the segmentation configuration to a file in JSON format.
    The configuration includes the number of chunks and padding information.
    input_file_name: Name of the stored image file
    chunk_num: Number of chunks after segmentation
    pad: Padding size (in bytes) for the last chunk
    '''
    config_path = f"./config/config.json"
    # Read existing config if it exists
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            try:
                config_data = json.load(f)
            except json.JSONDecodeError:
                config_data = {}
    else:
        config_data = {}
    # Update config with current image info
    config_data[in_file_name] = {"chunk_num": chunk_num, "pad": pad}

    # Write back to config file in JSON format
    with open(config_path, "w") as f:
        json.dump(config_data, f, indent=4, ensure_ascii=False)

save_config (in_file_name, len(data), pad)

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
'LDPC encoding...'
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
def LDPC_encode(data_pools, ldpc_encoder_exe, mode, h_matrix_path):
    """
    Encode data pools using LDPC encoder.
    Each pool is processed separately.
    Input: 
    - data_pools: List of pools, each pool is a list of bit strings (chunks)
    - ldpc_encoder_exe: Path to the LDPC encoder executable
    - mode: Mode for the LDPC encoder ('encode' or 'decode')
    - h_matrix_path: Path to the H matrix file for LDPC encoding
    Output:
    - encoded_data_pools: List of encoded pools, each pool is a list of encoded bit strings
    """
    encoded_data_pools = []
    for pool_idx, pool in enumerate(data_pools):
        # Prepare input file path
        input_file = f"D:/COPDS-main/Mid_data/ldpc_input_{pool_idx}.txt"
        with open(input_file, "w") as f:
            for chunk in pool:
                f.write(chunk + "\n")
        
        # Call LDPC encoder executable
        output_file = f"D:/COPDS-main/Mid_data/ldpc_output_{pool_idx}.txt"
        process = subprocess.Popen(
            [
                ldpc_encoder_exe,
                mode,
                input_file,
                output_file,
                h_matrix_path
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        status_placeholder = st.empty()
        for line in iter(process.stdout.readline, ''):
            if not line:
                break
            # Show only the latest output line in the app, in one line
            status_placeholder.text(line.strip())
        process.stdout.close()
        process.wait()

        # Read encoded data back
        encoded_pool = []
        with open(output_file, "r") as f:
            for line in f:
                encoded_pool.append(line.strip())
        encoded_data_pools.append(encoded_pool)
    return encoded_data_pools

data_pools = LDPC_encode(data_pools, ldpc_encoder_exe, mode, h_matrix_path)

# Convert each pool to a 2-D numpy array and transpose for BCH encoding
def transpose_v2h(data_pools):
    """
    Transpose the data pools from vertical to horizontal for BCH encoding.
    Each pool is a list of strings, each string is the i-th bit of all chunks.
    Transpose to get a 2-D numpy array of shape (num_chunks, chunk_length).
    """
    transposed_data_pools = []
    for pool in data_pools:
        # Each pool is a list of bit strings (chunks)
        # Convert to 2-D numpy array of shape (num_chunks, chunk_length)
        arr = np.array([list(chunk) for chunk in pool], dtype=np.uint8)
        # Transpose to shape (chunk_length, num_chunks)
        transposed_pool = arr.T
        transposed_data_pools.append(transposed_pool)
    return transposed_data_pools
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
input('Debug')
'''
#==================Binary to DNA==================#
st.subheader('Binary to DNA')
dna_pools = []
total_bits = chunk_size + k_1 + r_1 + r_2
dna_lenght = int(np.ceil(total_bits / 2))  # Each DNA nucleotide encodes 2 bits
for data_pool in data_pools:
    dnas = np.empty(n_0, dtype=f'<U{dna_lenght}') # Create an empty array for DNA sequences
    for idx in range(n_0):
        cw = data_pool[idx]
        bit_str = ''.join(cw.astype(int).astype(str))[:total_bits]
        dnas[idx] = bin_to_dna(bit_str)
    dna_pools.append(dnas)


# Create a subfolder in DNA_Library named as "<image_name> DNA pools"
dna_pools_dir = os.path.join(dna_lib_dir, 'DNA_pools-' + file_name)
os.makedirs(dna_pools_dir, exist_ok=True)

# Save each DNA pool to a separate .dna file in the new subfolder with suffix '_n'
for idx, dnas in enumerate(dna_pools):
    pool_file_name = os.path.join(dna_pools_dir, f'{file_name}_{idx}.dna')
    save_dna_files(dnas, pool_file_name)
'Data converted to DNA sequences, each sequence is of length ', len(dnas[0]), ' nucleotides.'
print('DNAs saved to ', dna_pools_dir)


# --------------------------- error simulation ---------------------------------- #
st.header('Error simulation of the DNA data storage channel')

st.subheader('Load Data')
# Folder selection for DNA pools

# uploaded_dna_files = st.file_uploader(
#     "Select one or more DNA pool to simulate", 
#     type=["dna"], 
#     accept_multiple_files=True
# )

if uploaded_dna_files:
    in_dnas_pools = []
    for uploaded_file in uploaded_dna_files:
        dnas = uploaded_file.read().decode("utf-8").splitlines()
        in_dnas = [dna.strip() for dna in dnas]
        in_dnas_pools.append(in_dnas)
    st.success(f"Loaded {len(in_dnas_pools)*len(in_dnas)} strands of length {len(in_dnas[0])} nts")
    # Update file name for later processing
    # Remove the suffix after the last "_" (including the "_") in the uploaded file name
    file_name = uploaded_file.name.rsplit('_', 1)[0]

    # 'Sequence ', index, ' will be inspected in detail to show how errors are formed in one sequence.\
    #      You can choose another sequence to inspect by altering the **inspect index** in the sidebar.' 
    # 'Three figures will be depicted for each stage:'
    # '1. Oligo copy number distribution and voting error number distribution.'
    # '2. Error types of sequence ', index
    # '3. Voting results of current copies of sequence ', index
    # 'The simulation process now begins:'


    # st.subheader('Synthesis')
    'Synthesis...'
    dnas_syn_pools = []
    for in_dnas in in_dnas_pools:
        SYN = Synthesizer(arg)
        dnas_syn = SYN(in_dnas)
        dnas_syn_pools.append(dnas_syn)
    # inspect(dnas_syn,inspect_index = index)

    # st.subheader('Decay')
    'Decay...'
    dnas_dec_pools = []
    for dnas_syn in dnas_syn_pools:
        DEC = Decayer(arg)
        dnas_dec = DEC(dnas_syn)
        dnas_dec_pools.append(dnas_dec)
    # inspect(dnas_dec,inspect_index = index)

    # st.subheader('PCR')
    'PCR...'
    dnas_pcr_pools = []
    for dnas_dec in dnas_dec_pools:
        PCR = PCRer(arg = arg)
        dnas_pcr = PCR(dnas_dec)
        dnas_pcr_pools.append(dnas_pcr)
    # inspect(dnas_pcr,inspect_index = index)

    # st.subheader('Sampling')
    'Sampling...'
    dnas_sam_pools = []
    for dnas_pcr in dnas_pcr_pools:
        SAM = Sampler(arg = arg)
        dnas_sam = SAM(dnas_pcr)
        dnas_sam_pools.append(dnas_sam)
    # inspect(dnas_sam,inspect_index = index)

    # st.subheader('Sequencing')
    'Sequencing...'
    dnas_seq_pools = []
    for dnas_sam in dnas_sam_pools:
        SEQ = Sequencer(arg)
        dnas_seq = SEQ(dnas_sam)
        dnas_seq_pools.append(dnas_seq)
    # inspect(dnas_seq,inspect_index = index)
    'Simulation completed. '

    # Save all dnas_seq_pools in a subfolder named "<image_name>_simu_DNA_pools" under DNA_Library
    simu_dna_pools_dir = os.path.join(dna_lib_dir, f'simu_DNA_pools-{file_name}')
    os.makedirs(simu_dna_pools_dir, exist_ok=True)
    for idx, dnas_seq in enumerate(dnas_seq_pools):
        pool_file_name = os.path.join(simu_dna_pools_dir, f'simu_{file_name}_{idx}.dna')
        # Extract every DNA after the simulation  pipeline
        dnas_sim_result = []
        for dna_set in dnas_seq:
            for dna_error_profile in dna_set['re']:
                for i in range(dna_error_profile[0]):
                    dnas_sim_result.append(dna_error_profile[2])

        # Open the file in write mode
        with open(pool_file_name, "w") as file:
            # Write each string in the list to the file
            for line in dnas_sim_result:
                file.write(line + "\n")
    'Simulated DNAs saved to ', simu_dna_pools_dir

# else:
#     st.stop()
#     st.warning("No DNA pool files selected.")



# --------------------------- decoding ---------------------------- #
# ===================== Read the simulated DNA ===================== #
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

#===============DNA to Binary====================#
st.subheader('DNA to Binary')
QUANT2BIN = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
out_data_pools = []
for out_dnas in out_dnas_pools:
    IDX = [] # Store the index bits arrays
    INF = [] # Store the information bits arrays
    for i, dna in enumerate(out_dnas):
        bin_str = ''
        if ((chunk_size + k_1 + r_1 + r_2) % 2 == 0):
            for b in dna:
                bin_str += QUANT2BIN[b]
        else:
            for j, b in enumerate(dna):
                if j != len(dna) - 1:
                    bin_str += QUANT2BIN[b]
                else:
                    bin_str += QUANT2BIN[b][1]
        
        ids = bin_str[:k_1 + r_1]
        inf = bin_str[k_1 + r_1:chunk_size + k_1 + r_1 + r_2]
        ids = np.array([float(s) for s in ids])
        inf = np.array([float(s) for s in inf])
        IDX.append(ids)
        INF.append(inf)
    IDX = np.array(IDX)
    INF = np.array(INF)
    # 'show data', IDX, INF
    out_data_pools.append((IDX, INF))

'Intra-oligos Decoding...'
temp = []
for pool_idx, (IDX, INF) in enumerate(out_data_pools):
    readsnum = len(out_data_pools[pool_idx][0]) # Number of reads in each pool
    # BCH decode index part
    eng = matlab.engine.start_matlab()
    eng.cd(r'D:\DeSP-main', nargout=0)  # Set MATLAB working directory
    cwr_ids = eng.BCH_Decoder(k_1 + r_1, k_1, readsnum, IDX)
    # BCH decode information part
    cwr_data = eng.BCH_Decoder(chunk_size + r_2, chunk_size, readsnum, INF)
    temp.append((cwr_ids, cwr_data))
# Store the decoded index and information parts
out_data_pools = temp
# 'show out_data_pools', out_data_pools

'voting...'
# Cluster the sequences with the same index in each pool
out_prob_pools = [] # Store the voting results for each pool
for pool_idx, (IDX, INF) in enumerate(out_data_pools):
    segments_temp = []
    for id in IDX:
        segment_temp = ''.join(str(int(s)) for s in id[1:])
        segments_temp.append(segment_temp)
    indices_dec_str = segments_temp

    segments_temp = []    
    for data in INF:
        segment_temp = ''.join(str(int(s)) for s in data[1:])
        segments_temp.append(segment_temp)
    inf_bit_dec_str = segments_temp

    segments_temp = []
    for i, inf_bit in enumerate(inf_bit_dec_str):
        segment_temp = dict(index=0, num=0, data=' ')
        segment_temp['index'] = int(indices_dec_str[i], 2)
        segment_temp['data'] = inf_bit
        segments_temp.append(segment_temp)
    segments = segments_temp
    segments_temp = []

    for i in range(n_0):
        segment_temp = dict(index=0, num=0, data=[])
        segment_temp['index'] = i
        for segment in segments:
            if segment['index'] == i and len(segment['data']) == chunk_size:
                segment_temp['num'] += 1
                segment_temp['data'].append(segment['data'])
        segments_temp.append(segment_temp)
    segments = segments_temp

    voting_result = []
    for i, segment in enumerate(segments):
        data = []
        if segment['num'] > 1:
            for j in range(chunk_size):
                bit_sum = 0
                for bit_string in segment['data']:
                    bit_sum += int(bit_string[j])
                data.append(float(bit_sum / segment['num']))
        elif segment['num'] == 1:
            data = [float(bit) for bit in segment['data'][0]]
        else:
            data = [0.5 for _ in range(chunk_size)]
        data = np.array(data)
        voting_result.append(data)
    voting_result = np.array(voting_result)
    v_score = voting_result.T
    out_prob_pools.append(v_score) # Probability of each bit in the voting result, for LDPC decoding

# Prepare and save the voting results for LDPC decoding
# v_score shape: (chunk_size, n_0)
# Save as text file for LDPC decoder input (each line is a bit position, values are probabilities for each chunk)
for pool_idx, v_score in enumerate(out_prob_pools):
    ldpc_input_file = f"D:\\DeSP-main\\App_Simulation_Platform\\Mid_data\\ldpc_decode_input_{pool_idx}.txt"
    np.savetxt(ldpc_input_file, v_score, fmt="%.3f", delimiter=" ")
'Voiting results saved to LDPC input files.'   

# ===================== LDPC Decoding ===================== #
'Inter-oligos Decoding...'
for pool_idx, v_score in enumerate(out_data_pools):
    # Call LDPC decoder executable
    ldpc_decoder_exe = "D:\\DeSP-main\\App_Simulation_Platform\\LDPC_PEG-v2.0.exe"
    h_matrix_path = "D:\\DeSP-main\\App_Simulation_Platform\\config\\Hmatrix.txt"
    ldpc_output_file = f"D:\\DeSP-main\\App_Simulation_Platform\\Mid_data\\ldpc_decode_output_{pool_idx}.txt"
    mode = "decode"
    with st.spinner(f"LDPC decoding pool {pool_idx}..."):
        process = subprocess.Popen(
            [
                ldpc_decoder_exe,
                mode,
                ldpc_input_file,
                ldpc_output_file,
                h_matrix_path
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        status_placeholder = st.empty()
        for line in iter(process.stdout.readline, ''):
            if not line:
                break
            status_placeholder.text(line.strip())
        process.stdout.close()
        process.wait()

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
st.write(f"Bit Error Rate (BER): {ber:.6f}")


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

out_file_path = "D:\\DeSP-main\\Data\\" + out_file_name
# Save reconstructed image
with open(out_file_path, "wb") as f:
    for chunk in decoded_chunks:
        f.write(chunk)

st.success(f"Image reconstructed and saved as {out_file_path}")

st.subheader('Quality Evaluation')
st.image(out_file_path, width = 300)

# ------------------------ optimizing ----------------------------- #
'''
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