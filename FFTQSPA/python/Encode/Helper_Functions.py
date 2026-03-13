import numpy as np
from Sequencing_Cost_Optimization_Analy.Crossover_Prob import Optimal_Allocation_InnerCode

#----------------transfromation functions---------------------#
BASE = ['A','C','G','T']

def bin_to_dna(bin_str):
    s = ''.join(BASE[int(bin_str[t:t+2],2)] for t in range(0, len(bin_str),2))
    return s

def int_to_binary_array(number, length):
    # Step 1 & 2: Convert the integer to a binary string and remove the '0b' prefix
    binary_str = bin(number)[2:]
    
    # Step 3: Pad the binary string with zeros at the beginning to ensure the desired length
    padded_binary_str = binary_str.zfill(length)
    
    # Step 4: Convert the binary string to a list of integers
    binary_array = [float(int(bit)) for bit in padded_binary_str]
    
    return np.array(binary_array)

#--------------------------BCH---------------------------------#
def inner_redundancy(P_e, d_seq, n_syn, N_c, r_i, K):
    """
    Return the closest effective optimal redundancy allocation r1 and r2 for a given R_i and P_e.
    Args:
        P_e: P_e: Overall error probability.
        d_seq: sequenncing depth.
        n_syn: synthesis number.
        N_c: the copy time during PCR.
        r_i: total redundancy for inner code.
        K: length of index bits and information bits [index,inf].
    Returns:
        r1, r2: Redundancy for the index bits and redundancy for the data bits.
    """
    delta_1, r, T = Optimal_Allocation_InnerCode(P_e, d_seq, n_syn, N_c, r_i, K)  # optimal inner code radundancy allocation for the given condition.

    # 从BCH的有效冗余中选取距离理论最优值最近的冗余。
    r1 = effective_reduandancy(r[0], K[0])
    r2 = effective_reduandancy(r[1], K[1])
    return r1, r2

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
    order = 0
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
