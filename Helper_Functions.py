import logging
import random
import sys
from math import sqrt
from reedsolo import RSCodec
import numpy as np
from ordered_set import OrderedSet
from scipy.spatial.distance import kulczynski1
import os
import subprocess
import json
import streamlit as st
import matplotlib.pyplot as plt
from Analysis.Analysis import dna_chunk, plot_oligo_number_distribution, plot_error_distribution, save_simu_result
import matplotlib.pyplot as plt
import matlab.engine
# from Sequencing_Cost_Optimization_Analy.Crossover_Prob import Optimal_Allocation_InnerCode


#-----------------segmentation----------------#
def file_to_indexed_dnas(file_name, chunk_size, index_length = None):
    data = preprocess(file_name, chunk_size)
    index_l = index_length
    if index_l == None:
        index_l = index_len(len(data))
    return data_to_dnas(data,index_l), index_l

def preprocess(file_name, chunk_size, is_text = False):
    data = data_from_file(file_name,is_text)
    return segments(data,chunk_size,is_text)

def data_from_file(file_name, is_text = False):
    try:
      if is_text == False:
        f = open(file_name, 'rb')
      else:
        f = open(file_name,'r')
    except: 
      logging.error("%s file not found", file_name)
      sys.exit(0)
    data = f.read()
    f.close()    
    return data

def segments(data, chunk_size,is_text = False):

    pad = -len(data) % chunk_size
    if pad > 0:
      logging.debug("Padded the file with %d zero to have a round number of blocks of data", pad)    
    if is_text == False:
        data += b'\0' * pad #zero padding.
    else:
        data += ' ' * pad 
    size = len(data)

    chunk_num = int(size/chunk_size)
    data_array = [None]*chunk_num
    for num in range(chunk_num):
        start = chunk_size * num
        end = chunk_size * (num+1)
        chunk_binary = data[start:end]
        data_array[num] = chunk_binary

    return data_array, pad

def lines_from_file(file_name):
    lines = []
    with open(file_name,'r') as f:
        while True:
            l = f.readline().split('\n')[0]
            if(l == ''):
                 break
            lines.append(l)
        f.close()
        return lines


def parse_int(f):
    line = f.readline()
    return int(line.split('\n')[0].split(' ')[1])

def load_dna(file_name):
    with open(file_name) as f:
        dnas = f.readlines()
    in_dnas = [dna.split('\n')[0] for dna in dnas]
    return in_dnas

def split_data(data, pool_size, chunk_size):
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

# Save segmentation configuration for later reconstruction
def save_config_segmentation(in_file_name, pools_num, chunk_num, pad):
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
    # Update segmentation config with current image info

    config_data[in_file_name]["segmentation"] = {"pools_num": pools_num, "chunk_num": chunk_num, "pad": pad}

    # Write back to config file in JSON format
    with open(config_path, "w") as f:
        json.dump(config_data, f, indent=4, ensure_ascii=False)

#-----------------scanner-------------------#
class Scanner:
    def __init__(self,max_repeat = 3, gc_interval = [0.45,0.55]):#maximum number of repeated nucleotide & allowable GC content
        self.max_repeat = max_repeat
        self.gc_interval = gc_interval
    #repetion information of nucleotide    
    def scan_repeats(self,dna,record_position = False):
        repeats = []
        prv = dna[0]
        r_num = 1
        for i,c in enumerate(dna[1:]):
            if prv == c:
                r_num += 1
            else:
                if(r_num > self.max_repeat):
                    if(record_position):
                        repeats.append([prv,r_num,i-r_num+1])
                    else:
                        repeats.append([prv,r_num])
                r_num = 1
                prv = c

        if(r_num > self.max_repeat): 
                    if(record_position):
                        repeats.append([prv,r_num,i-r_num+1])
                    else:
                        repeats.append([prv,r_num])
        return repeats
    
    #obtain maximum repeation number
    def max_repeats(self,dna):
        rs = self.scan_repeats(dna)
        if rs == []: return 0
        else: return max([r[1] for r in rs])
    
    #repeated segments relative to the maximum allowed repeats
    def repeats_point(self,dna):
        rs = self.scan_repeats(dna)
        if rs == []: return 0
        else:
            return sum([r[1] / (self.max_repeat + 1) for r in rs])
    
    #return GC content
    def Gc(self,dna):
        gc = dna.count('G') + dna.count('C')
        l = len(dna)
        return float(gc) / l
    
    #test if gc content satisfy requirement
    def gc_pass(self,dna):
        if self.gc_interval[0]  < self.Gc(dna) < self.gc_interval[1]:
            return True
        else:
            return False
    
    # If GC content and repetition rate satisfy requirement, pass = true
    def Pass(self,dna,with_primer = False):
        if self.gc_pass(dna) and self.repeats_point(dna)  == 0:
            return True
        else:
            return False
        
    #average gc content fo all dna strands
    def ave_gc(self,dnas):
        return sum([self.Gc(dna) for dna in dnas])/len(dnas)
    
    #total repeated nucleotide of all dnas strands
    def rp_total(self,dnas):
        return(sum([self.repeats_point(dna) for dna in dnas]))
    
    # select the best dna strand form copies
    def select_best(self,dnas):
        min_rp = 10000
        best_dna = dnas[0]
        for dna in dnas:
            if(self.gc_pass(dna)):
                if self.repeats_point(dna) < min_rp:
                    min_rp = self.repeats_point(dna)
                    best_dna = dna
        return best_dna,min_rp

    def analyze(self,dnas):
        gcs = [self.Gc(dna) for dna in dnas]
        gc_out_num = sum([not self.gc_pass(dna) for dna in dnas])#total number of nda that couldn't pass the requirement of gc content
        ave_Gc = self.ave_gc(dnas)
        
        rps = [self.repeats_point(dna) for dna in dnas]
        rp_too_long = sum([rp > 0 for rp in rps])
        dic = {
            'gc_list': gcs,
            'gc_out': gc_out_num,
            'homo_list': rps,
            'average_gc': ave_Gc,
            'homo_too_long': rp_too_long
        }
        return dic

#------------------xor---------------------#

xor_map = {
    'A':{
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T'
    },
    'C':{
        'A': 'C',
        'C': 'A',
        'G': 'T',
        'T': 'G'
    },
    'G':{
        'A': 'G',
        'C': 'T',
        'G': 'A',
        'T': 'C'
    },
    'T':{
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
}

def xor_dna(d1, d2):
    if(len(d1) != len(d1)):
        logging.error("length not equal")
    dr = ''.join([xor_map[c1][c2] for c1,c2 in zip(d1,d2)])
    return dr

def xor_ord(ord_array1, ord_array2):
       return [ord1 ^ ord2 for (ord1,ord2) in zip(ord_array1,ord_array2)]

def xor(byte_array1, byte_array2):
    return bytes([b1 ^ b2 for (b1,b2) in zip(byte_array1,byte_array2)])

#------------------random------------------#

def happen(prob):
    i = random.uniform(0,1)
    if i<=prob:
        return 1
    else:
        return 0
    
def random_base():
    r = random.random()
    if(r < 0.25): 
        return 'A'
    elif(r < 0.5):
        return 'C'
    elif(r < 0.75):
        return 'G'
    else:
        return 'T'
    
def random_dna(num):
    return ''.join([random_base() for i in range(num)])

def Random_Bit_Stream(n):
    '''
    Desccription: Generate random bit stream of length n.
    Input:
        n: length of bit stream to be generated.
    Output:
        bit_stream: Random bit string of length n.
    '''
    # Generate a random bit stream of length n
    bit_stream = ''.join(str(random.randint(0, 1)) for _ in range(n))
    return bit_stream

#----------------transfromation functions---------------------#
BASE = ['A','C','G','T']
QUANT = {'A': 0, 'C':1, 'G':2, 'T':3}

def dna_to_int_array(dna_str):
    #convert a string like ACTCA to an array of ints like [10, 2, 4]
    s = ''.join('{0:02b}'.format(QUANT[dna_str[t]]) for t in range(0, len(dna_str),1))
    return [int(s[t:t+8],2) for t in range(0,len(s), 8)]


#transform a number to qua_len bases
#example: num_to_dna(6,3) returns 'ACC'
def num_to_dna(num, qua_len):
    arr = []
    while True:
        lef = num % 4
        arr.append(lef)
        num = int(num / 4)
        if 0 == num:
            break

    outString = ''
    for n in arr:
        outString = BASE[n] + outString
    
    dl = qua_len - len(arr)
    if(dl < 0):
        logging.error('space not enough for encoding num')
        return -1
    elif(dl == 0):
        return outString
    else:
        return 'A' * dl + outString
    return outString

# dna <-> bytes
def byte_to_dna(s):
    #convert byte data (\x01 \x02) to DNA data: ACTC
    bin_data = ''.join('{0:08b}'.format(s[t]) for t in range(0,len(s)))
    return bin_to_dna(bin_data)

def dna_to_byte(dna):
    #convert a string like ACTCA to a string of bytes like \x01 \x02
    num = [QUANT[b] for b in dna]
    s = ''.join('{0:02b}'.format(num[t]) for t in range(0, len(num),1))
    data = b''.join(bytes([int(s[t:t+8],2)]) for t in range(0, len(s), 8))
    return data

def bin_to_dna(bin_str):
    s = ''.join(BASE[int(bin_str[t:t+2],2)] for t in range(0, len(bin_str),2))
    return s

def dna_to_num(dna):
    return sum([num * 4**i for i,num in enumerate([QUANT[b] for b in dna][::-1])])

def bytes2bits(byte_string):
    # Convert each byte to its binary representation and join them
    bit_sequence = ''.join(f'{byte:08b}' for byte in byte_string)
    return bit_sequence

def bits2bytes(bit_string):
    # Ensure the bit string length is a multiple of 8
    if len(bit_string) % 8 != 0:
        raise ValueError("Bit string length must be a multiple of 8")
    
    # Convert each group of 8 bits to a byte
    byte_array = bytearray(int(bit_string[i:i+8], 2) for i in range(0, len(bit_string), 8))
    return bytes(byte_array)

def dna_to_bin_str(dna):
    # Description
    '''
    Convert dna strings to binary list like '01000100010'.
    input: dna string
    output: binary list: binary_data, whose element is a dna string.
    '''
    num = [QUANT[b] for b in dna]
    bin_str = ''.join('{0:02b}'.format(num[t]) for t in range(0, len(num),1))
    return bin_str

def dna_to_bin_array(dna):
    # Convert dna strings to binary array like [0,1,0,0,0,1,0,0,0,1,0].
    # input: dna string to be converted
    # output: bin_str: binary ndarray
    num = [QUANT[b] for b in dna]
    bin_str = ''.join('{0:02b}'.format(num[t]) for t in range(0, len(num),1))
    bin_array =  [int(s) for s in bin_str]
    return np.array(bin_array)
    
def int_to_binary_array(number, length):
    # Step 1 & 2: Convert the integer to a binary string and remove the '0b' prefix
    binary_str = bin(number)[2:]

    # Step 3: Pad the binary string with zeros at the beginning to ensure the desired length
    padded_binary_str = binary_str.zfill(length)

    # Step 4: Convert the binary string to a list of integers
    binary_array = [int(bit) for bit in padded_binary_str]

    return np.array(binary_array)
    
def int_to_binary_array(number, length):
    # Step 1 & 2: Convert the integer to a binary string and remove the '0b' prefix
    binary_str = bin(number)[2:]
    
    # Step 3: Pad the binary string with zeros at the beginning to ensure the desired length
    padded_binary_str = binary_str.zfill(length)
    
    # Step 4: Convert the binary string to a list of integers
    binary_array = [float(int(bit)) for bit in padded_binary_str]
    
    return np.array(binary_array)
def binary_to_dna_pools(data_pools, dna_length, is_padding,dna_pools_dir):
    """
    Convert binary codewords to DNA sequences pools.
    Each codeword is a binary array, and each pool is a list of codewords.
    Args:
        data_pools (list of np.ndarray): List of binary codewords pools.
        dna_length (int): Length of each DNA sequence in nucleotides.
        is_padding (bool): Whether to pad the binary codewords to fit the DNA length.
    """
    # Convert binary codewords to DNA sequences
    dna_pools = []
    for data_pool in data_pools:
        dnas = []
        for binary_codeword in data_pool:
            if is_padding: # check length of each codeword, if less than dna_length*2, pad with 0s.
                binary_codeword = np.pad(binary_codeword, (0, dna_length*2 - len(binary_codeword)), 'constant')
            bit_str = ''.join(binary_codeword.astype(int).astype(str))
            dna = bin_to_dna(bit_str)
            dnas.append(dna)
        dna_pools.append(dnas)
    
    # ------------ Save DNAs to files --------------- #
    # Create a subfolder in DNA_Library named as "<image_name> DNA pools"
    os.makedirs(dna_pools_dir, exist_ok=True)
    # Save DNA pools to a .dna file, with different key for each pool
    pools_in_file_name = os.path.join(dna_pools_dir, 'in.dna')
    for idx, dnas in enumerate(dna_pools):
        # Open the file in write mode
        with open(pools_in_file_name, "w") as file:
            # Write DNA sequences as JSON: key is pool index, value is list of DNA sequences
            json.dump({f"pool{idx}": dnas}, file, ensure_ascii=False, indent=2)
    print('Synthesised DNAs saved to ', pools_in_file_name)
    return dna_pools

def DNA_pools_to_binary(coding_config, in_file_name, out_dna_pools):
    print('DNA to binary...')
    out_data_pools = []
    # separate index and information part
    for out_dnas in out_dna_pools:
        temp_dnas = [] # Temporary dictionary to store DNA {id, information part, counts}.
        for dna, num in out_dnas.items():
            dna = dna.strip()
            binary = dna_to_bin_array(dna) # Convert DNA to binary array
            id = binary[:coding_config["Coding_param"]["inner1"]["n"]] # Extract index part
            if coding_config[in_file_name]["Bin2DNA"]["is_padding"]: # Extract information part, remove padding bit if exists
                inf = binary[coding_config["Coding_param"]["inner1"]["n"]:-1] # Remove the padding bit
            else:
                inf = binary[coding_config["Coding_param"]["inner1"]["n"]:]
            temp_dnas.append([id, inf, num]) # Store the id, information part and counts
        out_data_pools.append(temp_dnas)
    return out_data_pools
# Save padding info of binary to DNA conversion
def save_padding_info_bin2DNA(file_name, is_padding, padding):
    '''
    Save the padding information for binary to DNA conversion.
    Input:
    - file_name: Name of the stored image file
    - is_padding: Boolean indicating if padding was added
    - padding: padding bits added (0 or 1)
    '''
    config_path = f"./config/config.json"
    # Read existing padding info if it exists
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            try:
                config_data = json.load(f)
            except json.JSONDecodeError:
                config_data = {}
    else:
        config_data = {}
    # Update padding info
    config_data[file_name]["Bin2DNA"] = {"is_padding": is_padding, "padding": padding}
    # Write back to padding info file in JSON format
    with open(config_path, "w") as f:
        json.dump(config_data, f, indent=4, ensure_ascii=False)

#-----------------------indexing dna chunks------------------------#
def data_to_dnas(data,index_length = 8):
    dnas = []
    for i,d in enumerate(data):
        d = num_to_dna(i,index_length) + byte_to_dna(d)
        dnas.append(d)
    return dnas

def dnas_to_data(dnas,chunk_num,index_length = 8):
    data_chunks = [b'' for b in range(chunk_num)]
    for dna in dnas:
        index = dna_to_num(dna[:index_length])
        payload = dna_to_byte(dna[index_length:])
        data_chunks[index] = payload
    return b''.join(data_chunks)

def index_len(chunk_num):
    return int(sqrt(sqrt(chunk_num))) + 1

#------------------------RS------------------------------------------#
def rs_decode(data, rs_obj = None, rs = None, max_hamming = 1):
    if not rs_obj:
        if rs: rs_obj = RSCodec(rs)
        else: 
            print('rs decoder not assigned.')
            return
    try:
        data_corrected = list(rs_obj.decode(data)[0])
    except:
        logging.debug('can not correct ori data')
        return -1, None 
    #we will encode the data again to evaluate the correctness of the decoding
    data_again = list(rs_obj.encode(data_corrected))
    if np.count_nonzero(data != list(data_again)) > max_hamming: #measuring hamming distance between raw input and expected raw input
        #too many errors to correct in decoding
        logging.debug('too many errors!')
        return -1, None
    return 0, data_corrected

#------------------save .dna files--------------------------------#
def save_dna_files(dnas,file_name = 'out.dna'):
    with open(file_name, 'w') as f:
        # f.write('Fountain code\n')
        # f.write('CN: ' + str(self.num_chunks) +'\n')
        # f.write('CL: ' + str(self.chunk_size) + '\n')
        # f.write('RS: ' + str(self.rs) + '\n')
        f.writelines('\n'.join([dna for dna in dnas]))
        f.close()
        print(f'DNAs have been saved to {file_name}')

#------------------read .dna files--------------------------------#
def read_dna_file(file_name):
    """
    Read a DNA file and return a list of DNA sequences.
    
    Args:
        encode_filepath (str): The path to the DNA file.
        
    Returns:
        list: A list of DNA sequences read from the file.
    """
    # Open the file in read mode
    with open(file_name, 'r') as file:
    # Read all lines from the file and add them to a list
        lines = file.readlines()
    # Strip newline characters from each line
    dnas = [line.strip() for line in lines]
    return dnas

#------------------extract dnas from output-----------------------#
def extract_dnas(out_dnas):
    '''
    Extracts simulated DNA sequences from the output of the DNA channel model.
    Input:
    - out_dnas: List of dictionaries, containing error profiles and corresponding DNA sequences.
    Output:
    - result: Dictionary of simulated DNA sequences, key:DNA, value:number.
    '''
    result = {}
    for dna_set in out_dnas:
        for dna_error_profile in dna_set['re']:
            result[dna_error_profile[2]] = dna_error_profile[0]
    return result


#--------------------------BCH---------------------------------#

def effective_reduandancy(r_theory, k):

    """
    Description: Return the effective redundancy for a given optimal reduandancy in theory and information bit k.
    Return: The the closest effective redundancy r.
    """
    eff_r = dict()
    eff_r[5] = [5, 10, 15]
    eff_r[6] = [24, 27, 33, 39, 45, 47]
    eff_r[7] = [77, 84, 91, 98, 99]
    eff_r[9] = [9, 18, 27, 36, 45, 54, 63, 72, 81, 90, 99, 108, 117, 126, 135, 144]
    # r = r_theory
    n = r_theory + k
    if (n >= 2**4 and n < 2**5):
        order = 5
        r = min(eff_r[5], key=lambda x: abs(x - r_theory))  # The closest effective redundancy to r_theory
    elif (n >= 2**5 and n < 2**6):
        order = 6
        r = min(eff_r[6], key=lambda x: abs(x - r_theory))  # The closest effective redundancy to r_theory
    elif (n >= 2**6 and n < 2**7):
        order = 7
        r = min(eff_r[7], key=lambda x: abs(x - r_theory))  # The closest effective redundancy to r_theory
    elif (n >= 2**8 and n < 2**9):
        order = 9
        r = min(eff_r[9], key=lambda x: abs(x - r_theory))  # The closest effective redundancy to r_theory
    if ((r + k) >= 2**order): # no effective
        r = effective_reduandancy(r, k) # Increase the order and refind.
    return r

class BCH_Codec:
    def __init__(self, encode_config = None, decode_config = None):
        if encode_config is None:
            self.n1 = decode_config["Coding_param"]["inner1"]["n"]
            self.k1 = decode_config["Coding_param"]["inner1"]["k"]
            self.k2 = decode_config["Coding_param"]["inner2"]["k"]
            self.n2 = decode_config["Coding_param"]["inner2"]["n"]
            self.n0 = decode_config["Coding_param"]["outer"]["n"]
        else:
            self.n1 = encode_config["n1"]
            self.k1 = encode_config["k1"]
            self.k2 = encode_config["k2"]
            self.n2 = encode_config["n2"]
            self.n0 = encode_config["n0"]


    def encode(self, ids, data_pools):
        print("BCH encoding...")
        matlab_engine = matlab.engine.start_matlab()
        matlab_engine.cd(r'D:\DeSP-main', nargout=0)  # Set MATLAB working directory
        eng = matlab_engine
        # Encode index by BCH encoder.
        st.text('BCH encoding...')
        cwr_ids = eng.BCH_Encoder(self.n1, self.k1, self.n0, ids)
        for pool_idx, pool in enumerate(data_pools):
            # Encode information bit by BCH encoder.
            cwr_data = eng.BCH_Encoder(self.n2, self.k2, self.n0, pool)
            data_pools[pool_idx] = np.concatenate((cwr_ids, cwr_data), axis=1)
        print("BCH encode completed.")
        return data_pools
    def decode(self, out_data_pools):
        print('BCH Decoding...')
        # Decode each pool
        eng = matlab.engine.start_matlab()
        eng.cd(r'D:\DeSP-main', nargout=0)  # Set MATLAB working directory
        temp_pools = []
        for pool in out_data_pools:
            ids = np.array([segment[0] for segment in pool])
            infs = np.array([segment[1] for segment in pool])
            ids = eng.BCH_Decoder(self.n1, self.k1, len(ids), ids)
            data = eng.BCH_Decoder(self.n2, self.k2, len(infs), infs)
            temp = []
            for i, id in enumerate(ids):
                temp_id = ''.join(str(int(s)) for s in id[1:]) # Convert index bits to string, remove the first bit which is the flag
                temp_data = ''.join(str(int(s)) for s in data[i][1:]) # Convert information bits to string, remove the first bit which is the flag
                temp.append((temp_id, temp_data, pool[i][2])) # Store the id, information part and counts
            temp_pools.append(temp)
        out_data_pools  = temp_pools
        print('BCH decode completed.')
        return out_data_pools
#----------------------LDPC----------------------------#
class LDPC_Codec:
    def __init__(self, h_matrix_path='D:/COPDS-main/config/Hmatrix.txt',LDPC_codec_path = "D:/COPDS-main/Encode/LDPC_PEG-v2.0.exe", mid_data_folder = "D:/COPDS-main/Mid_data/"):
        """
        Initialize the LDPC codec with the specified mode and H matrix path.
        :param mode: 'encode' or 'decode'
        :param h_matrix_path: Path to the H matrix file for encoding/decoding
        """
        self.h_matrix = h_matrix_path
        self.LDPC_codec_path = LDPC_codec_path
        self.mid_data_folder = mid_data_folder # folder for mid data
        if not os.path.exists(self.mid_data_folder):
            os.makedirs(self.mid_data_folder)
    def encode(self, data_pools):
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
        mode = "encode"
        encoded_data_pools = []
        for pool_idx, pool in enumerate(data_pools):
            # Prepare input file path
            input_file = self.mid_data_folder + f"ldpc_encode_input_{pool_idx}.txt"
            with open(input_file, "w") as f:
                for chunk in pool:
                    f.write(chunk + "\n")
            
            # Call LDPC encoder executable
            output_file = self.mid_data_folder + f"ldpc_encode_output_{pool_idx}.txt"
            # Call LDPC decoder executable
            mode = "encode"
            result = subprocess.run([self.LDPC_codec_path, mode, input_file, output_file, self.h_matrix], capture_output=True, text=True)
            # print(result.stdout)  # Print the output from the LDPC decoder
            print(result.stderr)  # Print any error messages from the LDPC decoder
            # Read encoded data back
            encoded_pool = []
            with open(output_file, "r") as f:
                for line in f:
                    encoded_pool.append(line.strip())
            encoded_data_pools.append(encoded_pool)
        return encoded_data_pools

    def decode(self, out_prob_pools):
        """
        Decode the input data using LDPC decoding.
        Input:
        - segment_num: Number of segments (pools) to decode 
        """
        print('LDPC Decoding...')
        # Save as text file for LDPC decoder input (each line is a bit position, values are probabilities for each chunk)
        for pool_idx, v_score in enumerate(out_prob_pools):
            input_file = self.mid_data_folder + f"ldpc_decode_input_{pool_idx}.txt"
            output_file = self.mid_data_folder + f"ldpc_decode_output_{pool_idx}.txt"
            
            np.savetxt(input_file, v_score, fmt="%.3f", delimiter=" ")
            # Call LDPC decoder executable
            mode = "decode"
            result = subprocess.run([self.LDPC_codec_path, mode, input_file, output_file, self.h_matrix], capture_output=True, text=True)
            # print(result.stdout)  # Print the output from the LDPC decoder
            # print(result.stderr)  # Print any error messages from the LDPC decoder
        print('Completed!')

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

def save_coding_config(config_path, outer_code = (1000, 800), inner1 = (20,10), inner2 = (80,20)):
    '''
    Save coding parameters to config file in JSON format.
    outer_code: (n_0, k_0)
    inner1: (n_1, k_1)
    inner2: (n_2, k_2)
    '''
    # Read existing padding info if it exists
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            try:
                config_data = json.load(f)
            except json.JSONDecodeError:
                config_data = {}
    else:
        config_data = {}
    # Update padding info
    config_data["Coding_param"] = {"outer": {"n": outer_code[0], "k": outer_code[1]}, "inner1": {"n": inner1[0], "k": inner1[1]}, "inner2": {"n": inner2[0], "k": inner2[1]}}
    # Write back to padding info file in JSON format
    with open(config_path, "w") as f:
        json.dump(config_data, f, indent=4, ensure_ascii=False)

# ----------------------voting----------------------------#
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


# ----------------------Webapp----------------------------#
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


# ----------------------Data analysis----------------------------#
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

