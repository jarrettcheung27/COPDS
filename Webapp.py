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

# ignore warnings
import warnings
warnings.filterwarnings("ignore")
config_path = "./config/config.json"
dna_lib_dir = "./DNA_Library/"

# =========================== assigning parameters ========================= #

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
coding_config = {"Coding_param" : {"outer": {"n": n_0, "k": k_0}, "inner1": {"n": k_1 + r_1, "k": k_1}, "inner2": {"n": k_2 + r_2, "k": k_2}}}
save_coding_config(config_path, outer_code = (n_0, k_0), inner1 = (k_1 + r_1, k_1), inner2 = (chunk_size + r_2, chunk_size)) # The same coding parameters are used for all images.

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

st.header('Upload file')
# 使用file_uploader 上传多张图片
uploaded_files = st.file_uploader(
    "Select images to store in DNA",
    type=["jpg", "jpeg", "png"],
    accept_multiple_files=True,
    key="file_uploader"
)
if uploaded_files is not None:
    for uploaded_file in uploaded_files:
        # 使用 PIL 打开图片
        image = Image.open(uploaded_file)
        st.image(image, caption="上传的图片", width=100)

# If no files are uploaded, use a default image
if not st.session_state.get("file_uploader"):
    # Use a default image if no files are uploaded
    default_image_path = "D:/COPDS-main/DNA_Library/Jnu_test.jpg"
    st.session_state["file_uploader"] = [default_image_path]

# Initialize coding configuration
coding_config = {'files': {
                    'num': 0},
                 'ECC': {},
                 'DNA2Binary': {}
                 }
image_ids = [str(uuid.uuid4()) for _ in range(len(uploaded_files))]  # Generate unique IDs for each image
image_names = [file.name for file in uploaded_files]  # Get the names of the uploaded files
name_to_id_mapping = {image_names[i]: image_ids[i] for i in range(len(uploaded_files))}

for i, file in enumerate(uploaded_files):
    coding_config['files'][i] = {
        "id": image_ids[i],
        "file_name": file.name,
        "file_type": file.type,
        "suffix": file.name.split('.')[-1],
        "file_size": file.size,
        "segmentation": {
            "chunk_num": 0,
            "block_num": 0,
            "pad_size": 0
        }
    }
coding_config['files']["num"] = len(uploaded_files)
print(coding_config)


print('\n ==================== DNA 存储仿真平台 ==================== ')
# ======================Select image to store ====================== #
# Create the DNA_Library directory if it doesn't exist
os.makedirs(dna_lib_dir, exist_ok=True)

# =====================Segmentation==================== #
# 数据分层结构 Library->pool(image)->block(coding unit)->chunk(dna)
st.subheader('Segmentation')
total_chunk = 0
Library = []
for i, file in enumerate(uploaded_files):
    block, pad_size = split_image(file, chunk_size)
    total_chunk += len(block)
    # bytes to bits
    block = [bytes2bits(chunk) for chunk in block]
    coding_config['files'][i]['segmentation']['chunk_num'] = len(block)
    # Split data into pools of size k_0 (1000), pad the last pool if needed
    pool = split_data(block, block_size=k_0, chunk_size=chunk_size)
    coding_config['files'][i]['segmentation']['block_num'] = len(pool)
    coding_config['files'][i]['segmentation']['pad_size'] = pad_size

    Library.append(pool)
print('Images loaded and split into ', total_chunk, ' data chunks.')
print('config:', coding_config)
# Save input message for debugging
'''
for pool_idx, pool in enumerate(data_pools):
    message_input_file = f"D:\DeSP-main\App_Simulation_Platform\Mid_data\ldpc_encode_input_{pool_idx}.txt"
    with open(message_input_file, "w") as f:
        for chunk in pool:
            f.write(chunk + "\n")
'''

# -----------------------Indexing-----------------------
ids  = np.array([int_to_binary_array(id, k_1) for id in range(n_0)])

# =================================Encoding===================================#
print('Processing Encode...')
# ---------------------Save configuration---------------------
coding_config['ECC'] = {
    "outer": {
        "n": n_0,
        "k": k_0,
        "h_matrix_path": ".\config\Hmatrix.txt",
        "ldpc_encoder_exe": ".\Encode\LDPC_PEG-v2.0.exe"
    },
    "inner1": {
        "n": k_1 + r_1,
        "k": k_1
    },
    "inner2": {
        "n": chunk_size + r_2,
        "k": chunk_size
    }
}
LDPC = LDPC_Codec()
BCH = BCH_Codec(coding_config = coding_config)
temp = []
for pool in Library:
    # LDPC encode
    # Transpose matrix for inter-oligos encoding
    pool = transpose_h2v(pool)
    # 'Data is transposed for inter-oligos encoding, each pool is a list of strings, each string is the i-th bit of all chunks.'

    # Path to the LDPC encoder executable and H matrix
    pool = LDPC.encode(pool)
    # Convert each pool to a 2-D numpy array and transpose for BCH encoding
    pool = transpose_v2h(pool)
    # ----------------------BCH encode----------------------
    pool = BCH.encode(ids, pool)
    temp.append(pool)
Library = temp
print('Encoding completed.')
print('config:', coding_config)
print(len(Library), 'pools generated.')
for i, pool in enumerate(Library):
    print(f'Pool {i}: {len(pool)} blocks')

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
for i, pool in enumerate(Library):
    dna_pool_filename = os.path.join(dna_lib_dir, f"{coding_config['files'][i]['id']}_in.dna")
    temp = binary_to_dna_pools(pool, dna_length, is_padding, dna_pool_filename)
    DNA_Library.append(temp)
print('Binary to DNA conversion completed.')
print('DNA pools saved to:', dna_lib_dir)

# ==================Save coding configuration================= #
with open(config_path, "w") as f:
    json.dump(coding_config, f, indent=4, ensure_ascii=False)
print('Coding configuration saved to:', config_path)
# ============================= DNA simulation Channel ============================= #
st.header('Error simulation of the DNA data storage channel')
st.subheader('Load Data')
Channel = DNA_Channel_Model(Modules = 0, arg = arg) # No model provided

for i, dna_pools in enumerate(DNA_Library):
    # Run Simulation for each pool
    temp = []
    for dnas in dna_pools:
        out_dnas = Channel(dnas)
        temp.append(out_dnas)
    DNA_Library[i] = temp
print('Simulation completed.')


# --------------Save Simulation results ----------------
for i, pool in enumerate(DNA_Library):
    dna_pool_filename = os.path.join(dna_lib_dir, f"{coding_config['files'][i]['id']}_out.dna")
    temp = []
    for idx, out_dnas in enumerate(pool):
        # Extract simulated DNA sequences
        dnas = extract_dnas(out_dnas)
        # Open the file in write mode
        if idx == 0:
            with open(dna_pool_filename, "w") as file:
                # Write DNA sequences as JSON: key is pool index, value is list of DNA sequences
                json.dump({f"chunk_{idx}": dnas}, file, ensure_ascii=False, indent=2)
        else:   
            # Append DNA sequences to the existing file
            with open(dna_pool_filename, "a") as file:
                json.dump({f"chunk_{idx}": dnas}, file, ensure_ascii=False, indent=2)
    temp.append(dnas)
    DNA_Library[i] = temp
print('Simulated DNAs saved to ', dna_lib_dir)


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
pool_num = coding_config[in_file_name]["segmentation"]["block_num"]
chunk_num = coding_config[in_file_name]["segmentation"]["chunk_num"]
pad_size = coding_config[in_file_name]["segmentation"]["pad_size"]
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

# Remove padding bytes in the last chunk if pad_size > 0
if pad_size > 0:
    last_chunk = decoded_chunks[-1]
    if len(last_chunk) > pad_size:
        decoded_chunks[-1] = last_chunk[:-pad_size]  # Remove padding bytes from the last chunk 

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