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
import sys
from PIL import Image
import uuid
import pdb # python debugger
# ignore warnings
import warnings
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
warnings.filterwarnings("ignore")

CONFIG_PATH = "./config/config.json"
DNA_LIB_DIR = "./DNA_Library/"
MID_DATA_DIR = "./Mid_data/"
DEBUG_MODE = False
def load_config_or_default():
    if os.path.exists(CONFIG_PATH):
        try:
            with open(CONFIG_PATH, "r", encoding="utf-8") as f:
                return json.load(f)
        except json.JSONDecodeError:
            pass
    return {"file_num": 0, "files": {}, "ECC": {}, "DNA2Binary": {}}

def save_config(config_data):
    with open(CONFIG_PATH, "w", encoding="utf-8") as f:
        json.dump(config_data, f, indent=4, ensure_ascii=False)

def list_stored_files_ui():
    coding_config = load_config_or_default()
    files = coding_config.get("files", {})
    if not files:
        st.info("DNA库中没有已存储的文件。")
        return
    rows = []
    total_blocks = 0
    for file_id, info in files.items():
        block_num = info.get("segmentation", {}).get("block_num", 0)
        total_blocks += block_num
        rows.append({
            "DNA 池 id": file_id,
            "文件名": info.get("file_name", ""),
            "DNA 块数": block_num,
        })
    st.subheader("DNA库已存文件")
    st.dataframe(rows, use_container_width=True)
    st.info(f"总DNA块数（block_num）: {total_blocks}")

def clear_dna_library():
    os.makedirs(DNA_LIB_DIR, exist_ok=True)
    removed = 0
    for name in os.listdir(DNA_LIB_DIR):
        if name.lower().endswith(".dna"):
            file_path = os.path.join(DNA_LIB_DIR, name)
            try:
                os.remove(file_path)
                removed += 1
            except OSError:
                continue
    coding_config = load_config_or_default()
    coding_config["file_num"] = 0
    coding_config["files"] = {}
    save_config(coding_config)
    st.success(f"已清空DNA库，共删除 {removed} 个 .dna 文件，并更新配置。")

def clear_mid_data():
    os.makedirs(MID_DATA_DIR, exist_ok=True)
    removed = 0
    for name in os.listdir(MID_DATA_DIR):
        file_path = os.path.join(MID_DATA_DIR, name)
        if os.path.isfile(file_path):
            try:
                os.remove(file_path)
                removed += 1
            except OSError:
                continue
    return removed

if DEBUG_MODE:
    coding_config = json.load(open(CONFIG_PATH, "r"))
    file_ids = list(coding_config['files'].keys())
# =========================== assigning parameters ========================= #
print('Running Webapp...')

if st.button("查询已存文件"):
    list_stored_files_ui()

if st.button("清空 DNA 库"):
    st.warning("已执行清空DNA库操作，恢复功能将不再可用。")
    clear_dna_library()

st.subheader("选择模式")
with st.form(key="selectMode"): # key 是这个表单的标识符
    st.session_state['mode'] = st.selectbox('模式', ['编码并存入 DNA 库', '从 DNA 库中读取并译码',])
    mode_submitted  = st.form_submit_button("提交")
# ------------ Choosing parameters -------------- # 
st.sidebar.subheader('编码参数')
chunk_size = st.sidebar.number_input('数据分块大小', min_value = 320, max_value = 320, value = 320)
k_0 = st.sidebar.number_input('外码信息位长度 $k_0$', min_value = 1024, max_value = 1024, value = 1024)
n_0 = st.sidebar.number_input('外码码长 $n_0$', min_value = 2048, max_value = 2048, value = 2048)
k_1 = st.sidebar.number_input('索引长度 $k_1$', min_value = 14, max_value = 14, value = 14)
r_1 = st.sidebar.number_input('第一层TL-BCH冗余 $r_1$', min_value = 15, max_value = 15, value = 15)
k_2 = chunk_size
r_2 = st.sidebar.number_input('第二层TL-BCH冗余 $r_2$', min_value = 99, max_value = 99, value = 99)
# --------------------Channel parameter----------------
arg = config.DEFAULT_PASSER
st.subheader('DNA数据存储信道参数')
with st.form(key="channelParams"): # 创建一个表单
    # synthesis stage
    arg.syn_sub_prob = st.slider('信道噪声水平', min_value = 0.01, max_value = 0.4, value = 0.01) / 3 # 3 kinds of errors
    arg.syn_number = st.slider('合成数量', min_value = 10, max_value = 50, value = 20)
    arg.syn_yield = st.slider('合成产率', min_value = 0.98, max_value = 0.9999, value = 0.9999)

    # PCR stage
    arg.pcrc = st.slider('PCR循环次数',min_value = 0, max_value =20,value =1)
    arg.pcrp = st.number_input('PCR概率',min_value = 0.5, max_value = 1.0,value = 0.9)

    # decay stage
    arg.decay_er = 0
    arg.decay_loss_rate = 0.01

    # sequencing stage
    seq_platform = st.selectbox('测序平台',['Illumina测序','纳米孔'])
    if seq_platform == 'Illumina测序':
        arg.seq_TM = config.TM_NGS
    else:
        arg.seq_TM = config.TM_NNP
    arg.sam_ratio = st.slider('采样比例',min_value = 0.0, max_value =1.0,value = 1.0)
    arg.seq_depth = st.slider('测序深度', min_value = 1, max_value = 100, value = 10)

    # inspect index
    index = st.slider('查看索引', max_value = 600, value = 0)

    channel_submitted  = st.form_submit_button("提交") # 提交表单
if st.session_state['mode'] == '编码并存入 DNA 库':
    st.header('上传文件')
    with st.form(key='select_file'):
        # 使用file_uploader 上传多张图片
        uploaded_files = st.file_uploader(
            "选择要存入DNA的图片",
            type=["jpg", "jpeg", "png"],
            accept_multiple_files=True,
            key="file_uploader"
        )

        encode_submitted = st.form_submit_button("提交")

    if uploaded_files is not None:
        for uploaded_file in uploaded_files:
            # 使用 PIL 打开图片
            image = Image.open(uploaded_file)
            st.image(image, caption="上传的图片", width=100)

    if not uploaded_files:
        st.info("请至少上传一张图片，然后点击提交。")
        st.stop()

    if not encode_submitted:
        st.info("上传图片并点击提交开始编码与存储。")
        st.stop()

    if os.path.exists(CONFIG_PATH):
        try:
            coding_config_exist = json.load(open(CONFIG_PATH, "r", encoding="utf-8"))
        except json.JSONDecodeError:
            coding_config_exist = {'file_num':0, 'files': {}, 'ECC': {}, 'DNA2Binary': {}}
    else:
        coding_config_exist = {'file_num':0, 'files': {}, 'ECC': {}, 'DNA2Binary': {}}

    # Initialize coding configuration
    coding_config = {'file_num':0,
                    'files': {},
                    'ECC': {},
                    'DNA2Binary': {}
                    }
    file_ids = [str(uuid.uuid4()) for _ in range(len(uploaded_files))]  # Generate unique IDs for each image

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

    fftqspa_module = f"./Encode/Nonbinary/fftqspa.cp{sys.version_info.major}{sys.version_info.minor}-win_amd64.pyd"

    # ---------------------coding configuration---------------------
    coding_config['ECC'] = {
        "outer": {
            "n": n_0,
            "k": k_0,
            "h_matrix_path": "./Encode/Nonbinary/Parity_files_2048/512_256_16aryCode17.dat",
            "ldpc_encoder_exe": fftqspa_module,
            "parity_path": "./Encode/Nonbinary/Parity_files_2048/512_256_16aryCode17.dat",
            "mapping_path": "./Encode/Nonbinary/Mapping_files/SignalSet_BPSK-4.txt",
            "max_iter": 50,
            "code_ary": 16
        },
        "inner1": {
            "n": k_1 + r_1,
            "k": k_1
        },
        "inner2": {
            "n": chunk_size + r_2,
            "k": chunk_size
        },
        "BCH_codec_exe": fftqspa_module
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
    print(f"{len(uploaded_files)} image(s) split into {total_chunk} data chunks.")


    # =================================Encoding===================================#
    st.info('正在编码...')
    LDPC = LDPC_Codec(coding_config)
    BCH = BCH_Codec(coding_config)  # Initialize BCH codec with config path
    bch_encoded_library = []
    for i, pool in enumerate(Library):

        # Transpose matrix for inter-oligos encoding
        pool = transpose_h2v(pool)

        # Process LDPC encode for each pool, and save the encoded output
        LDPC.encode(pool = pool,file_id = file_ids[i])

        # Process BCH encode for each pool
        print('BCH Encoding...')
        bch_encoded_pool = BCH.encode(file_id = file_ids[i])
        bch_encoded_library.append(bch_encoded_pool)
    print('Encoding completed.')
    #===============================Binary to DNA===============================#
    st.info('正在转换二进制数据为DNA...')
    if not bch_encoded_library or not bch_encoded_library[0]:
        st.error("BCH输出为空，编码已中止。")
        st.stop()

    total_bits = int(bch_encoded_library[0][0].shape[1])
    dna_length = int(np.ceil(total_bits / 2))  # Each DNA nucleotide encodes 2 bits
    is_padding = (total_bits % 2) != 0
    coding_config['DNA2Binary'] = {
        "is_padding": is_padding,
        "padding": int(dna_length * 2 - total_bits),
    }
    DNA_Library = []
    # Convert full BCH codewords (id + data) to DNA
    for i in range(coding_config['file_num']):
        pool = bch_encoded_library[i]
        dna_pool_filename = os.path.join(DNA_LIB_DIR, f"{file_ids[i]}_in.dna")
        temp = binary_to_dna_pools(pool, dna_length, is_padding, dna_pool_filename)
        DNA_Library.append(temp)
    print('Binary to DNA conversion completed.')
    print('DNA pools saved to:', DNA_LIB_DIR)

    # ==================Update coding configuration================= #
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
    with open(CONFIG_PATH, "w", encoding="utf-8") as f:
        json.dump(coding_config, f, indent=4, ensure_ascii=False)
    print('Coding configuration saved to:', CONFIG_PATH)
    
    # ============================= DNA simulation Channel ============================= #
    st.info('正在模拟DNA数据存储信道...')
    Channel = DNA_Channel_Model(Modules = 0, arg = arg) # No model provided

    for i, dna_pools in enumerate(DNA_Library):
        # Run Simulation for each pool
        temp = []
        for dnas in dna_pools:
            out_dnas = Channel(dnas)
            temp.append(out_dnas)
        DNA_Library[i] = temp
    print('仿真完成。')

    # --------------Save Simulation results ----------------
    for i, pool in enumerate(DNA_Library):
        dna_pool_filename = os.path.join(DNA_LIB_DIR, f"{file_ids[i]}_out.dna")
        out_pool = {}
        for idx, out_dnas in enumerate(pool):
            # Extract simulated DNA sequences
            out_pool[f"chunk_{idx}"] = extract_dnas(out_dnas)
            # Open the file in write mode

        with open(dna_pool_filename, "w") as file:
            # Write DNA sequences as JSON: key is pool index, value is list of DNA sequences
            json.dump(out_pool, file, ensure_ascii=False, indent=2)
    print('Simulated DNAs saved to ', DNA_LIB_DIR)
    st.success('文件已存入 DNA 库。')
    removed_mid_files = clear_mid_data()
    # 给出DNA库的文件可跳转的文件夹绝对路径链接
    st.markdown(f"DNA库路径: `{os.path.abspath(DNA_LIB_DIR)}`")

elif st.session_state['mode'] == '从 DNA 库中读取并译码':
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
        st.warning("DNA库中没有已存储的文件，请先上传文件。")
        st.stop()
    else:
        file_names = [coding_config['files'][id]['file_name'] for id in file_ids]
        file_ids = list(coding_config['files'].keys())
        with st.form(key="restoreFiles"): # 创建表单
            out_file_names = st.multiselect(
                "选择要恢复的文件",
                options=file_names,
                default=[file_names[0]]  # Default to the first file
            )
            restore_submitted  = st.form_submit_button("提交")

    if not restore_submitted:
        st.info("选择文件并点击提交开始恢复。")
        st.stop()
    # -------------- obtain the IDs of the selected files --------------
    out_file_ids = [file_ids[i] for i, name in enumerate(file_names) if name in out_file_names]
    if not out_file_ids:
        st.warning("未选择任何要恢复的文件。")
        st.stop()
    print(f'Selected files: {out_file_names}, IDs: {out_file_ids}')

    # -------------- Load the DNA pools from the selected files --------------
    DNA_Library = []
    loaded_file_ids = []
    strands_num = 0;
    for file_id in out_file_ids:
        # f'Loading the DNA of the selected files {coding_config["files"][file_id]["file_name"]}...'
        print(f'Loading DNA pool for file ID: {file_id}')
        out_dna_pool = []
        out_dna_path = os.path.join(DNA_LIB_DIR, f"{file_id}_out.dna")
        if not os.path.exists(out_dna_path):
            st.error(f"未找到 DNA 库文件: {out_dna_path}")
            continue
        DNA_pool = json.load(open(out_dna_path, "r", encoding="utf-8"))
        block_num = coding_config['files'][file_id]['segmentation']['block_num']
        for i in range(block_num):
            chunk_key = f'chunk_{i}'
            if chunk_key not in DNA_pool:
                st.error(f"{out_dna_path} 中缺少 {chunk_key}")
                out_dna_pool = []
                break
            out_dna_pool.append(DNA_pool[chunk_key])
            strands_num += len(DNA_pool[chunk_key])
        if not out_dna_pool:
            continue
        DNA_Library.append(out_dna_pool)
        loaded_file_ids.append(file_id)
    out_file_ids = loaded_file_ids
    if not out_file_ids:
        st.warning("没有可用于恢复的有效 DNA 库数据。")
        st.stop()
    st.info(f"已加载 {strands_num} 条长度为 {chunk_size} nts 的DNA序列用于恢复。")

    def infer_inner2_n_from_pool(dna_pool, cfg):
        inner1_n = cfg['ECC']['inner1']['n']
        is_padding = cfg.get('DNA2Binary', {}).get('is_padding', False)
        for block in dna_pool:
            if isinstance(block, dict) and block:
                first_dna = next(iter(block.keys()))
                total_bits = len(first_dna) * 2
                if is_padding:
                    total_bits -= 1
                inferred_n2 = total_bits - inner1_n
                if inferred_n2 > 0:
                    return inferred_n2
        return cfg['ECC']['inner2']['n']

    def observed_total_bits_from_pool(dna_pool):
        for block in dna_pool:
            if isinstance(block, dict) and block:
                first_dna = next(iter(block.keys()))
                return len(first_dna) * 2
        return None

    restore_configs = []
    legacy_data_only_files = []
    for idx, dna_pool in enumerate(DNA_Library):
        file_id = out_file_ids[idx]
        file_cfg = copy.deepcopy(coding_config)
        observed_bits = observed_total_bits_from_pool(dna_pool)
        if observed_bits is not None:
            padding_bits = 1 if file_cfg.get('DNA2Binary', {}).get('is_padding', False) else 0
            expected_full_bits = file_cfg['ECC']['inner1']['n'] + file_cfg['ECC']['inner2']['n'] + padding_bits
            expected_data_only_bits = file_cfg['ECC']['inner2']['n'] + padding_bits
            legacy_data_only = (
                observed_bits in {expected_data_only_bits, expected_data_only_bits + 1}
                and observed_bits < expected_full_bits
            )
            if legacy_data_only:
                legacy_data_only_files.append(file_cfg['files'][file_id]['file_name'])
                restore_configs.append(file_cfg)
                continue

        inferred_n2 = infer_inner2_n_from_pool(dna_pool, file_cfg)
        configured_n2 = file_cfg['ECC']['inner2']['n']
        if inferred_n2 != configured_n2:
            st.warning(
                f"{file_cfg['files'][file_id]['file_name']}: 从 DNA 库检测到 inner2.n={inferred_n2}，"
                f"为恢复兼容性已覆盖配置值 {configured_n2}。"
            )
            file_cfg['ECC']['inner2']['n'] = inferred_n2
        restore_configs.append(file_cfg)

    if legacy_data_only_files:
        st.error(
            "无法恢复缺少BCH索引位的旧版 DNA 库格式: "
            + ", ".join(legacy_data_only_files)
            + "。请使用更新后的流程重新编码这些文件，然后再恢复。"
        )
        st.stop()

    #================= DNA to Binary ====================#
    st.info('正在转换DNA为二进制数据...')
    temp = []
    for i, pool in enumerate(DNA_Library): # 提取出已选择的文件对应的 DNA 库，并转换为二进制数据
        out_data_pool = DNA_pool_to_binary_pool(restore_configs[i], pool)
        temp.append(out_data_pool)
        print(out_data_pool[0][0].shape)
        print(out_data_pool[0][1].shape)
        print(out_data_pool[0][2].shape)
    Library = temp

    print('DNA to binary conversion completed.')
    print('Number of pools:', len(Library))
    print('Number of blocks in each pool:', [len(pool) for pool in Library])

    st.info('正在译码..')
    print('Voting and decoding...')
    Prob_Library = [] # Store the voting results for each pool
    for i, pool in enumerate(Library):
        prob_pool = [] # Store the voting results for each pool
        for j, block in enumerate(pool):
            print(f'Voting for block {j} in pool {i}...')
            temp = bch_decode_and_vote_cpp(block, restore_configs[i])
            prob_pool.append(temp) # Probability of each bit in the voting result, for LDPC decoding
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

    st.info('正在重建文件...')
    # Read decoded output to reconstruct the image

    # Read the config file to obtain the configuration.
    for i in range(len(Library)):
        decoded_chunks = []
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
            st.success(f"图像 '{file_name}' 提取成功，已保存到 {out_file_path}")
            img = Image.open(out_file_path)
            st.image(img, caption="恢复的图像", width=400)
        except Exception as e:
            st.error(f"图像 '{file_name}' 提取失败: {e}")

    removed_mid_files = clear_mid_data()
