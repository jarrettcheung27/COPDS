from Encode.Helper_Functions import dna_to_int_array, rs_decode, preprocess, load_dna, bytes2bits
from Analysis.Analysis import dna_chunk, save_simu_result,  error_distribution
from Encode.DNAFountain import DNAFountain, Glass
import numpy as np
from reedsolo import RSCodec, ReedSolomonError
from pyldpc import make_ldpc, encode, decode, get_message
from scipy.stats import gumbel_r, poisson
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import clear_output

# def error_profile(out_dnas, rs = 2):
#     lost_num = 0
#     fail_num = 0
#     mis_judge_num = 0
#     ignore_index = []
#     for i,dna in enumerate(out_dnas):
        
#         if dna['num'] == 0: 
#             lost_num += 1
#             continue
        
#         dc = dna_chunk(dna)
#         re_dna = dc.voting_result()
        
#         re_data = dna_to_int_array(re_dna)
#         flag, data_corrected = rs_decode(re_data, rs = rs)

#         if flag == -1:
#             fail_num += 1
#         else:
#             ori_data = dna_to_int_array(dna['ori'])
#             if ori_data[:-2] != data_corrected:
#                 mis_judge_num += 1
#                 fail_num += 1
#                 ignore_index.append(i)
#     return lost_num, fail_num, mis_judge_num, ignore_index

class InterLdpcEncoder:
    def __init__(self, file_name, chunk_size = 20, encode = False):
        # Segment input data into datachunks of size chunk_size(byte)
        # input:
            # filename
            # chunk_size: size of data segment after dividing.
            # encode: whether proceed encode process.
        #output: 
            # self.data: data after segmentation.
        
        # declare files path
        file_name, file_type = file_name.split('.')
        self.ori_path = 'files/' + file_name + '.' + file_type
        self.encode_path = 'files/' + file_name + '.dna'
        self.simu_path = 'files/' + 'simu_' + file_name + '.dna'
        self.simu_decode_path = 'files/' + 'simu_re_' + file_name + '.' + file_type
        
        # data segmentation
        self.chunk_size = chunk_size
        self.data, self.pad = preprocess(self.ori_path, chunk_size = chunk_size)#preprocess data and split it into chunks
        self.chunk_num = len(self.data)
        print(f'Data split into {len(self.data)} data chunks.')
        self.bit_sequences = [] # Storage bit stream
        # channel model
        # self.model = model
        # distribution
        # self.loss_nums = []
        # self.fail_nums = []
        # self.decode_lines = [] #success number?
        if encode: self.encode()

    def bytes2bits(self):
        for byte_string in self.data:
            # Convert each byte to its binary representation and join them
            binary_sequence = ''.join(f'{byte:08b}' for byte in byte_string)
            self.bit_sequences.append(binary_sequence)

    # def Inter_Oligos_LDPC_Encode(self):
        
                

def Change_Direction(sequences): 
    '''
    Change Coding direction horizontal/vertical
    input: Sequencs to be converted
    output: 
        Converted_sequences: Sequences in another direction.
    '''
    bit_lists = [list(bit_stream) for bit_stream in sequences]
    # Use zip to aggregate the i-th elements from each bit list
    aggregated = zip(*bit_lists)
    # Join the elements of each aggregated group into new bit streams
    Converted_sequences = [''.join(group) for group in aggregated]
    return Converted_sequences


def LDPC_encode(bit_sequences, n, k, G):
    '''
    LDPC encode for bit sequences
    input:
        bit_sequences: original sequences in bit string.
        n: code length after encoding
        k: message length.
        G: Generate matrix for LDPC.
    output:
        Encoded_sequences in bit array.
    '''
    #LDPC encode
    # print(f'Length of message: {k}, LDPC code rate: {k/n}' )
    encoded_sequences = []
    for i, bit_sequence in enumerate(bit_sequences):
        # padding sequence length to k
        if k !=  len(bit_sequence): 
            s = sequences_padding(bit_sequence, k)
        else:
            bit_list = list(map(int, bit_sequence))
            s = np.array(bit_list) #Convert bit string to bit array.
            # print("type of G: ", type(G))
            # print("size of G: ", G.shape)
            # print("type of s: ", type(s))
            # print("size of s: ", s.size)
            # display(s)
        encoded_sequence = np.dot(G,s) % 2
        encoded_sequences.append(encoded_sequence)
        # Clear history outputs, and display recent state.
        clear_output(wait=True)
        print('Execute LDPC encoding:')
        print(f'{i + 1}/{len(bit_sequences)} messages have been successfully encoded by LDPC encoder.')
        
    # print(f'Encoded symbol: {encoded_sequences[0]}')
    # print(f'Oligos number: {len(encoded_sequences[0])}')
    # print(f'Oligos length: {len(encoded_sequences)}')
    return encoded_sequences

def RS_encode1(sequences, rs_length, m = 8):
    '''
    One layer Rs encode
    input: 
        sequences:input sequences in bits
        re_length: ecc symbols
    output: encoded sequences in byte  
    '''
    rsc = RSCodec(rs_length, c_exp = m)  # ecc symbols number
    encoded_sequences = []
    for sequence in sequences:
        encoded_sequence = rsc.encode(bits_to_bytearray(sequence))
        encoded_sequences.append(encoded_sequence)
    return encoded_sequences

def RS_decode1(sequences, rs_length, m = 8):
    '''
    One layer Rs decode
    input: 
        sequences:input sequences in byte
        re_length: ecc symbols
    output: decoded sequences in bits   
    '''
    rsc = RSCodec(rs_length, c_exp = m)  # ecc symbols number
    fail_num = 0 # number of sequences that can not be corrected
    success_num = 0 # number of sequences that can be corrected
    decoded_sequences = []
    for i, sequence in enumerate(sequences):
        print(f'RS decoder is decoding... {(i + 1)/len(sequences) * 100: .3f}', '%')
        try:
            rsc.check(sequence)
            decoded_sequence = rsc.decode(sequence)
            decoded_sequence = bytes2bits(decoded_sequence[0])
            success_num += 1
        except ReedSolomonError as e:
            decoded_sequence = sequence[:-rs_length]
            decoded_sequence = bytes2bits(decoded_sequence)
            fail_num += 1
        decoded_sequences.append(decoded_sequence)
        # print(f'{i} sequences have been decoded.')
    print(f'{success_num} sequences have been successfully decoded by RS decoder, while {fail_num} fails. Success rate: {success_num/len(sequences)}.')
    return decoded_sequences

def RS_encode2(sequences, index_length, rs_length_index, rs_length_data):
    '''
    Two layer Rs encode
    Input: 
        sequences: input sequences in bits strings.
        index_length: length of index part in the sequences in bits.
        rs_length_index: ecc symbols length for index in bytes.
        rs_length_data: ecc symbols length for data in bytes.
    Output: 
        encoded_sequences: encoded sequences in byte.
    '''
    rsc_index = RSCodec(rs_length_index)  # encoder for index.
    rsc_data = RSCodec(rs_length_data) # encoder for data.
    encoded_sequences = []
    for sequence in sequences:
        bytes_sequence = bits_to_bytearray(sequence)# Convert a bit string to byte array
        encoded_index = rsc_index.encode(bytes_sequence[: int(index_length/8)]) # Encode index
        encoded_data = rsc_data.encode(bytes_sequence[int(index_length/8) :]) #Encode data
        encoded_sequence = encoded_index + encoded_data
        encoded_sequences.append(encoded_sequence)
    return encoded_sequences

def RS_decode2(sequences, index_length, rs_length_index, rs_length_data):
    '''
    Two layer Rs decode
    Input: 
        sequences: input sequences in byte.
        index_length: length of index part in the sequences in bits.
        rs_length_index: ecc symbols length for index in bytes.
        rs_length_data: ecc symbols length for data in bytes.
    Output: 
        decoded_sequences: decoded sequences in bits.
    '''
    rsc_index = RSCodec(rs_length_index)  # encoder for index.
    rsc_data = RSCodec(rs_length_data) # encoder for data.

    index_fail_num = 0 # number of sequences whose index can not be corrected.
    index_success_num = 0 # number of sequences whose index can be corrected.
    
    data_fail_num = 0 # number of sequences whose data can not be corrected.
    data_success_num = 0 # number of sequences whose data can be corrected.
    decoded_sequences = []
    for i, sequence in enumerate(sequences):
        print(f'RS decoder is decoding... {(i + 1) / len(sequences) * 100: .3f}', '%')
        index = sequence[:(int(index_length/8) + rs_length_index)]
        data = sequence[(int(index_length/8) + rs_length_index):]
        try:
            rsc_index.check(index)
            decoded_index = rsc_index.decode(index)
            decoded_index = bytes2bits(decoded_index[0])
            index_success_num += 1
        except ReedSolomonError as e:
            decoded_index = index[:-rs_length_index]
            decoded_index = bytes2bits(decoded_index)
            index_fail_num += 1

        try:
            rsc_data.check(data)
            decoded_data = rsc_data.decode(data)
            decoded_data = bytes2bits(decoded_data[0])
            data_success_num += 1
        except ReedSolomonError as e:
            decoded_data = data[:-rs_length_data]
            decoded_data = bytes2bits(decoded_data)
            data_fail_num += 1    
        decoded_sequences.append(decoded_index + decoded_data)
        # print(f'{i} sequences have been decoded.')
    print(f'{index_success_num} indices have been successfully decoded by RS decoder, while {index_fail_num} fails. Success rate: {index_success_num/len(sequences)}.')
    print(f'{data_success_num} data chunk have been successfully decoded by RS decoder, while {data_fail_num} fails. Success rate: {data_success_num/len(sequences)}.')
    return decoded_sequences
    
def sequences_padding(ori_sequences, length, padding_value=0):
    '''
    Pad a sequence to a certain length for ldpc encode,return the result in array format.
    input:
        ori_sequences: sequences to be padded.
        length: sequences' length after padding.
        padding_value: value in padding position.
    '''
    # Convert the bit stream string to a list of binary number
    bit_list = list(map(int, ori_sequences))
    result = bit_list + [padding_value] * (length - len(bit_list))
    
    # Convert the list of binary number back to an array
    padded_bit_stream = np.array(result)
    
    return padded_bit_stream
    
def bits_to_bytearray(bits):
    '''
    Convert bit strings to byte array
    Input:
        bits: a bit string
    Output:
        byte_array: a byte array
    '''
    byte_array = bytearray()
    for i in range(0, len(bits), 8):
        byte = bits[i:i+8]
        # Ensure the byte has exactly 8 bits by padding with zeros if necessary
        # byte += [0] * (8 - len(byte))
        # Convert the list of bits to a single byte
        byte_value = 0
        for bit in byte:
            byte_value = (byte_value << 1) | bit
        byte_array.append(byte_value)
    return byte_array
    