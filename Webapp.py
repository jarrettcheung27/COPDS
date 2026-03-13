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
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
warnings.filterwarnings("ignore")

CONFIG_PATH = "./config/config.json"
DNA_LIB_DIR = "./DNA_Library/"
DEBUG_MODE = True

if DEBUG_MODE:
    coding_config = json.load(open(CONFIG_PATH, "r"))
    file_ids = list(coding_config['files'].keys())
# =========================== assigning parameters ========================= #
print('Running Webapp...')
with st.form(key="selectMode"): # key 是这个表单的标识符
    st.session_state['mode'] = st.selectbox('Mode', ['Encode & Store in DNA', 'Restore from DNA',])
    but_submitted  = st.form_submit_button("Submit")

if st.session_state['mode'] == 'Encode & Store in DNA':
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
    st.subheader('Parameters of DNA data storage channel')
    with st.form(key="channelParams"): # 创建一个表单
        # synthesis stage
        arg.syn_number = st.slider('Syn number', min_value = 10, max_value = 50, value = 20)
        # arg.syn_sub_prob = st.number_input('Syn Error rate', min_value = 0.0, max_value = 0.5, value = 0.00001) / 3 # 3 kinds of errors
        arg.syn_yield = st.slider('Syn Yield', min_value = 0.98, max_value = 0.995, value = 0.9999)

        # PCR stage
        arg.pcrc = st.slider('PCR cycle',min_value = 0, max_value =20,value =1)
        arg.pcrp = st.number_input('PCR prob',min_value = 0.5, max_value = 1.0,value = 0.9)

        # decay stage
        arg.decay_er = 0
        arg.decay_loss_rate = 0.01

        # sequencing stage
        seq_platform = st.selectbox('Sequencing Platform',['Illumina Sequencing','Nanopore'])
        if seq_platform == 'Illumina Sequencing':
            arg.seq_TM = config.TM_NGS
        else:
            arg.seq_TM = config.TM_NNP
        arg.seq_TM = genTm(0.00001) # sequencing Transform Matrix
        arg.sam_ratio = st.slider('Sampling ratio',min_value = 0.0, max_value =1.0,value = 1.0)
        arg.seq_depth = st.slider('Seq Depth', min_value = 1, max_value = 100, value = 10)

        # inspect index
        index = st.slider('inspect index', max_value = 600, value = 0)

        but_submitted  = st.form_submit_button("Submit") # 提交表单

    st.header('Upload file')
    with st.form(key='select_file'):
        # 使用file_uploader 上传多张图片
        uploaded_files = st.file_uploader(
            "Select images to store in DNA",
            type=["jpg", "jpeg", "png"],
            accept_multiple_files=True,
            key="file_uploader"
        )

        but_submitted = st.form_submit_button("Submit")

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
    uploaded_file_num = len(uploaded_files)
    coding_config['file_num'] = uploaded_file_num

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
    # ======================Select image to store ====================== #
    # Create the DNA_Library directory if it doesn't exist
    os.makedirs(DNA_LIB_DIR, exist_ok=True)

    # =====================Segmentation==================== #
    # 数据分层结构 Library->pool(image)->block(coding unit)->chunk(dna)
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
    f'{coding_config['file_num']} image split into ', total_chunk, ' data chunks.'


    # ==================Save coding configuration================= #
    with open(CONFIG_PATH, "w") as f:
        json.dump(coding_config, f, indent=4, ensure_ascii=False)
    print('Coding configuration saved to:', CONFIG_PATH)


    # =================================Encoding===================================#
    'Processing encode...'
    LDPC = LDPC_Codec(coding_config)
    BCH = BCH_Codec(coding_config)  # Initialize BCH codec with config path
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
    #===============================Binary to DNA===============================#
    'Converting Binary data to DNA...'
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

    # ==================Update coding configuration================= #
    coding_config_exist = json.load(open(CONFIG_PATH, "r"))
    # Merge existing and new file info instead of replacing
    existing_files = coding_config_exist.get('files', {})
    new_files = coding_config.get('files', {})
    # Add/overwrite new entries
    existing_files.update(new_files)
    coding_config_exist['files'] = existing_files
    coding_config_exist['file_num'] = len(existing_files)
    # Keep latest ECC and DNA2Binary sections from current run (if present)
    for _k in ['ECC', 'DNA2Binary']:
        if _k in coding_config:
            coding_config_exist[_k] = coding_config[_k]
    # Use merged config for saving
    coding_config = coding_config_exist
    with open(CONFIG_PATH, "w") as f:
        json.dump(coding_config, f, indent=4, ensure_ascii=False)
    print('Coding configuration saved to:', CONFIG_PATH)
    
    # ============================= DNA simulation Channel ============================= #
    'Error simulation of the DNA data storage channel...'
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
        dna_pool_filename = os.path.join(DNA_LIB_DIR, f"{file_ids[i]}_out.dna")
        temp = []
        out_pool = {}
        for idx, out_dnas in enumerate(pool):
            # Extract simulated DNA sequences
            out_pool[f"chunk_{idx}"] = extract_dnas(out_dnas)
            # Open the file in write mode

        with open(dna_pool_filename, "w") as file:
            # Write DNA sequences as JSON: key is pool index, value is list of DNA sequences
            json.dump(out_pool, file, ensure_ascii=False, indent=2)
        temp.append(dnas)
        DNA_Library[i] = temp
    print('Simulated DNAs saved to ', DNA_LIB_DIR)
    st.success('Files are stored into DNA pools.')

elif st.session_state['mode'] == 'Restore from DNA':
        # ===================== Load the configuration file ===================== #
    config_path = "./config/config.json"
    with open(config_path, 'r', encoding='utf-8') as f:
        coding_config = json.load(f)
    file_ids = list(coding_config['files'].keys())
    restore_file_folder = "./Restored_files/"
    chunk_size = coding_config['ECC']['inner2']['k']
    # ===================== Select file to decode and restore===================== #
    store_file_num = coding_config['file_num']

    if store_file_num == 0:
        st.warning("No files are stored in the DNA Library. Please upload files first.")
        st.stop()
    else:
        file_names = [coding_config['files'][id]['file_name'] for id in file_ids]
        file_ids = list(coding_config['files'].keys())
        with st.form(key="restoreFiles"): # 创建表单
            out_file_names = st.multiselect(
                "Select files to restore",
                options=file_names,
                default=file_names[0]  # Default to the first file
            )
            but_submitted  = st.form_submit_button("Submit")
    # -------------- obtain the IDs of the selected files --------------
    out_file_ids = [file_ids[i] for i, name in enumerate(file_names) if name in out_file_names]
    print(f'Selected files: {out_file_names}, IDs: {out_file_ids}')

    # -------------- Load the DNA pools from the selected files --------------
    DNA_Library = []
    strands_num = 0;
    for file_id in out_file_ids:
        # f'Loading the DNA of the selected files {coding_config["files"][file_id]["file_name"]}...'
        print(f'Loading DNA pool for file ID: {file_id}')
        out_dna_pool = []
        DNA_pool = json.load(open(os.path.join(DNA_LIB_DIR, f"{file_id}_out.dna"), "r"))
        block_num = coding_config['files'][file_id]['segmentation']['block_num']
        for i in range(block_num):
            out_dna_pool.append(DNA_pool[f'chunk_{i}'])
            strands_num += len(DNA_pool[f'chunk_{i}'])
        DNA_Library.append(out_dna_pool)
    f"Loaded {strands_num} strands of length {chunk_size} nts"

    #================= DNA to Binary ====================#
    'Converting DNA to Binary data...'
    temp = []
    for pool in DNA_Library: # 提取出已选择的文件对应的DNA池，并转换为二进制数据
        out_data_pool = DNA_pool_to_binary_pool(coding_config,pool)
        temp.append(out_data_pool)
        print(out_data_pool[0][0].shape)
        print(out_data_pool[0][1].shape)
        print(out_data_pool[0][2].shape)
    Library = temp

    print('DNA to binary conversion completed.')
    print('Number of pools:', len(Library))
    print('Number of blocks in each pool:', [len(pool) for pool in Library])

    'Processing decode...'
    print('Voting and decoding...')
    BCH = BCH_Codec(coding_config = coding_config)
    Prob_Library = [] # Store the voting results for each pool
    for i, pool in enumerate(Library):
        # ============================== BCH Decoding ============================= #
        nums = [pool[j][2] for j in range(len(pool))] # record the reads number for each sequence.
        BCH.pool_to_txt(pool = pool, file_id = file_ids[i])
        pool = BCH.decode(file_ids[i])
        temp = []
        for j, block in enumerate(pool):
            block.append(nums[j])
            temp.append(block)
        pool = temp
        # Cluster the sequences with the same index in each pool
        prob_pool = [] # Store the voting results for each pool
        for j, block in enumerate(pool):
            print(f'Voting for block {j} in pool {i}...')
            temp = voting(block, coding_config)
            prob_pool.append(temp) # Probability of each bit in the voting result, for LDPC decoding
            # pdb.set_trace()
        Prob_Library.append(prob_pool)

    # ===================== LDPC Decoding ===================== #

    # 并行 LDPC 解码（按文件并行）
    print('LDPC Decoding (multiprocessing)...')
    # 每个元素: (prob_pool, file_id, coding_config)
    ldpc_tasks = [(Prob_Library[idx], out_file_ids[idx])
                for idx in range(len(Prob_Library))]

    def _ldpc_decode(task):
        prob_pool, file_id = task
        # 每个进程各自实例
        ldpc = LDPC_Codec(coding_config)
        ldpc.decode(prob_pool, file_id=file_id)
        return file_id

    # workers = min(len(ldpc_tasks), os.cpu_count()- 1)
    workers = 1
    print(f'The number of worker is {workers}')

    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_ldpc_decode, t): t[1] for t in ldpc_tasks}
            for fut in as_completed(futures):
                fid = fut.result()
                print(f'LDPC decode finished for file {fid}')
    else:
        # 回退到串行
        for t in ldpc_tasks:
            _ldpc_decode(t)

    print('LDPC decoding completed.')

    'Reconstructing the files...'
    # Read decoded output to reconstruct the image

    # Read the config file to obtain the configuration.
    decoded_chunks = []
    for i in range(len(Library)):
        file_id = out_file_ids[i]
        block_num = coding_config['files'][file_id]["segmentation"]["block_num"]
        chunk_num = coding_config['files'][file_id]["segmentation"]["chunk_num"]
        pad_size = coding_config['files'][file_id]["segmentation"]["pad_size"]
        for block_idx in range(block_num):
            ldpc_output_file = "./Mid_data/" + f"{file_id}_LDPC_decode_out_{block_idx}.txt"
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
        file_name = coding_config['files'][file_id]["file_name"]

        if not os.path.exists(restore_file_folder):
            os.makedirs(restore_file_folder)
        out_file_path = os.path.join(restore_file_folder, file_name)
        if os.path.exists(out_file_path):
            os.remove(out_file_path)  # Remove existing file to avoid appending to it
        with open(out_file_path, "wb") as f:
            for chunk in decoded_chunks:
                f.write(chunk)
        
        try:
            img = Image.open(out_file_path)
            st.image(img, caption="Restored Image", width=200)
            st.success(f"图像 '{file_name}' 提取成功，已保存到 {out_file_path}")
        except Exception as e:
            st.error(f"图像 '{file_name}' 提取失败: {e}")
