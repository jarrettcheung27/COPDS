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
coding_config = {'file_num':0,
                 'files': {},
                 'ECC': {},
                 'DNA2Binary': {}
                 }
file_ids = [str(uuid.uuid4()) for _ in range(len(uploaded_files))]  # Generate unique IDs for each image
file_names = [file.name for file in uploaded_files]  # Get the names of the uploaded files
name_to_id_mapping = {file_names[i]: file_ids[i] for i in range(len(uploaded_files))}

# ---------------------seg configuration---------------------
for i, file in enumerate(uploaded_files):
    coding_config['files'][file_ids[i]] = {
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
coding_config['file_num'] = len(uploaded_files)

# ---------------------coding configuration---------------------
coding_config['ECC'] = {
    "outer": {
        "n": n_0,
        "k": k_0,
        "h_matrix_path": "./config/Hmatrix.txt",
        "ldpc_encoder_exe": "./Encode/LDPC_PEG-v2.0.exe"
    },
    "inner1": {
        "n": k_1 + r_1,
        "k": k_1
    },
    "inner2": {
        "n": chunk_size + r_2,
        "k": chunk_size
    },
    "BCH_codec_exe": "./Encode/BCH_Codec.exe"
}

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
    coding_config['files'][file_ids[i]]['segmentation']['chunk_num'] = len(block)
    # Split data into pools of size k_0 (1000), pad the last pool if needed
    pool = split_data(block, block_size=k_0, chunk_size=chunk_size)
    coding_config['files'][file_ids[i]]['segmentation']['block_num'] = len(pool)
    coding_config['files'][file_ids[i]]['segmentation']['pad_size'] = pad_size

    Library.append(pool)
print('Images loaded and split into ', total_chunk, ' data chunks.')
# ==================Save coding configuration================= #
with open(config_path, "w") as f:
    json.dump(coding_config, f, indent=4, ensure_ascii=False)
print('Coding configuration saved to:', config_path)

# Save input message for debugging
'''
for pool_idx, pool in enumerate(data_pools):
    message_input_file = f"D:\DeSP-main\App_Simulation_Platform\Mid_data\ldpc_encode_input_{pool_idx}.txt"
    with open(message_input_file, "w") as f:
        for chunk in pool:
            f.write(chunk + "\n")
'''

# =================================Encoding===================================#
print('Processing Encode...')
LDPC = LDPC_Codec(coding_config)
BCH = BCH_Codec(config_path)  # Initialize BCH codec with config path
temp = []
for i, pool in enumerate(Library):

    # Transpose matrix for inter-oligos encoding
    pool = transpose_h2v(pool)

    # Process LDPC encode for each pool, and save the encoded output
    LDPC.encode(pool = pool,file_id = file_ids[i])

    # Process BCH encode for each pool
    print('BCH Encoding...')
    BCH.encode(file_id = file_ids[i])
print('Encoding completed.')


