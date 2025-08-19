import json
from Helper_Functions import *

# ------------- file path--------------- #
abs_dir  = "D:/COPDS-main/"
dna_lib_dir = abs_dir + 'DNA_Library/'
in_file_path = "D:/COPDS-main/DNA_Library/Jnu_test.jpg"
in_file_name= "Jnu_test.jpg"
file_name = "Jnu_test"
suffix = "jpg"
out_file_path = "D:/COPDS-main/DNA_Library/" + file_name + "_re." + suffix
dna_pools_dir = os.path.join(dna_lib_dir, 'DNA_pools-' + file_name)

# Read decoded output to reconstruct the image
in_file_name= "Jnu_test.jpg"
config_path = "./config/config.json"
with open(config_path, 'r', encoding='utf-8') as f:
    coding_config = json.load(f)
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