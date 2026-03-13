import copy
import json
import os
import uuid
import hashlib
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Helper_Functions import (
    split_image,
    bytes2bits,
    split_data,
    transpose_h2v,
    LDPC_Codec,
    BCH_Codec,
    binary_to_dna_pools,
    extract_dnas,
    DNA_pool_to_binary_pool,
    bch_decode_and_vote_cpp,
    bits2bytes,
)
from Model.Model import DNA_Channel_Model
import Model.config as model_config


def sha256_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as file_obj:
        for chunk in iter(lambda: file_obj.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    input_path = project_root / "Test_Files" / "default.jpg"
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}")
        return 1

    config_path = project_root / "config" / "config.json"
    dna_library_dir = project_root / "DNA_Library"
    mid_data_dir = project_root / "Mid_data"
    restored_dir = project_root / "Restored_files"

    dna_library_dir.mkdir(parents=True, exist_ok=True)
    mid_data_dir.mkdir(parents=True, exist_ok=True)
    restored_dir.mkdir(parents=True, exist_ok=True)

    chunk_size = 320
    k_0 = 1024
    n_0 = 2048
    k_1 = 14
    r_1 = 15
    r_2 = 99

    file_id = str(uuid.uuid4())
    file_name = input_path.name
    file_bytes = input_path.read_bytes()

    coding_config = {
        "file_num": 1,
        "files": {
            file_id: {
                "file_name": file_name,
                "file_type": "image/jpeg",
                "suffix": input_path.suffix.lstrip("."),
                "file_size": len(file_bytes),
                "segmentation": {
                    "chunk_num": 0,
                    "block_num": 0,
                    "pad_size": 0,
                },
            }
        },
        "ECC": {
            "outer": {
                "n": n_0,
                "k": k_0,
                "h_matrix_path": "./Encode/Nonbinary/Parity_files_2048/512_256_16aryCode17.dat",
                "ldpc_encoder_exe": "./Encode/Nonbinary/fftqspa.cp312-win_amd64.pyd",
                "parity_path": "./Encode/Nonbinary/Parity_files_2048/512_256_16aryCode17.dat",
                "mapping_path": "./Encode/Nonbinary/Mapping_files/SignalSet_BPSK-4.txt",
                "max_iter": 50,
                "code_ary": 16,
            },
            "inner1": {"n": k_1 + r_1, "k": k_1},
            "inner2": {"n": chunk_size + r_2, "k": chunk_size},
            "BCH_codec_exe": "./Encode/Nonbinary/fftqspa.cp312-win_amd64.pyd",
        },
        "DNA2Binary": {"is_padding": False, "padding": 0},
    }

    with input_path.open("rb") as file_obj:
        blocks, pad_size = split_image(file_obj, chunk_size)
    bit_chunks = [bytes2bits(chunk) for chunk in blocks]
    pools = split_data(bit_chunks, block_size=k_0, chunk_size=chunk_size)

    coding_config["files"][file_id]["segmentation"]["chunk_num"] = len(bit_chunks)
    coding_config["files"][file_id]["segmentation"]["block_num"] = len(pools)
    coding_config["files"][file_id]["segmentation"]["pad_size"] = pad_size

    ldpc = LDPC_Codec(coding_config)
    bch = BCH_Codec(coding_config)

    vertical_pools = transpose_h2v(pools)
    ldpc.encode(pool=vertical_pools, file_id=file_id)
    bch_encoded_pool = bch.encode(file_id=file_id)

    total_bits = int(bch_encoded_pool[0].shape[1])
    dna_length = int(np.ceil(total_bits / 2))
    is_padding = (total_bits % 2) != 0
    coding_config["DNA2Binary"] = {
        "is_padding": is_padding,
        "padding": int(dna_length * 2 - total_bits),
    }

    in_dna_path = dna_library_dir / f"{file_id}_in.dna"
    dna_pools = binary_to_dna_pools(bch_encoded_pool, dna_length, is_padding, str(in_dna_path))

    channel_arg = copy.deepcopy(model_config.DEFAULT_PASSER)
    channel_arg.syn_yield = 1.0
    channel_arg.syn_number = 30
    channel_arg.syn_sub_prob = 0.0
    channel_arg.syn_ins_prob = 0.0
    channel_arg.syn_del_prob = 0.0
    channel_arg.decay_er = 0.0
    channel_arg.decay_loss_rate = 0.0
    channel_arg.pcrc = 0
    channel_arg.pcrp = 0.9
    channel_arg.sam_ratio = 1.0
    channel_arg.seq_depth = 30
    channel_arg.seq_TM = 0

    channel = DNA_Channel_Model(Modules=0, arg=channel_arg)
    simulated_pool = {}
    for idx, dnas in enumerate(dna_pools):
        out_dnas = channel(dnas, print_state=False)
        simulated_pool[f"chunk_{idx}"] = extract_dnas(out_dnas)

    out_dna_path = dna_library_dir / f"{file_id}_out.dna"
    out_dna_path.write_text(json.dumps(simulated_pool, ensure_ascii=False, indent=2), encoding="utf-8")

    prob_pool = []
    out_dna_dict = json.loads(out_dna_path.read_text(encoding="utf-8"))
    decoded_input_pool = [out_dna_dict[f"chunk_{i}"] for i in range(coding_config["files"][file_id]["segmentation"]["block_num"])]
    binary_pool = DNA_pool_to_binary_pool(coding_config, decoded_input_pool)
    for block in binary_pool:
        prob_pool.append(bch_decode_and_vote_cpp(block, coding_config))

    ldpc.decode(prob_pool, file_id=file_id)

    decoded_chunks = []
    block_num = coding_config["files"][file_id]["segmentation"]["block_num"]
    chunk_num = coding_config["files"][file_id]["segmentation"]["chunk_num"]
    for block_idx in range(block_num):
        ldpc_out_file = mid_data_dir / f"{file_id}_LDPC_decode_out_{block_idx}.txt"
        with ldpc_out_file.open("r", encoding="utf-8") as file_obj:
            pool_bits = [line.strip() for line in file_obj if line.strip()]
        pool_bits_t = list(zip(*pool_bits))
        for chunk_bits_tuple in pool_bits_t:
            decoded_chunks.append(bits2bytes("".join(chunk_bits_tuple)))

    decoded_chunks = decoded_chunks[:chunk_num]
    if pad_size > 0 and decoded_chunks:
        decoded_chunks[-1] = decoded_chunks[-1][:-pad_size]

    restored_path = restored_dir / f"e2e_{file_name}"
    with restored_path.open("wb") as file_obj:
        for chunk in decoded_chunks:
            file_obj.write(chunk)

    source_sha = sha256_file(input_path)
    restored_sha = sha256_file(restored_path)
    identical = source_sha == restored_sha

    if config_path.exists():
        try:
            existing_cfg = json.loads(config_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            existing_cfg = {"file_num": 0, "files": {}, "ECC": {}, "DNA2Binary": {}}
    else:
        existing_cfg = {"file_num": 0, "files": {}, "ECC": {}, "DNA2Binary": {}}

    existing_cfg.setdefault("files", {})
    existing_cfg["files"][file_id] = coding_config["files"][file_id]
    existing_cfg["file_num"] = len(existing_cfg["files"])
    existing_cfg["ECC"] = coding_config["ECC"]
    existing_cfg["DNA2Binary"] = coding_config["DNA2Binary"]
    config_path.write_text(json.dumps(existing_cfg, ensure_ascii=False, indent=4), encoding="utf-8")

    print("[E2E] file_id:", file_id)
    print("[E2E] input:", input_path)
    print("[E2E] restored:", restored_path)
    print("[E2E] source_sha256:", source_sha)
    print("[E2E] restored_sha256:", restored_sha)
    print("[E2E] exact_match:", identical)

    return 0 if identical else 2


if __name__ == "__main__":
    raise SystemExit(main())
