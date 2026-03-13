import math
import numpy as np
import csv
from scipy.ndimage import maximum

from BER_NA import Rate_blocklength_NA, Calc_Capacity, Calc_Dispersion
from Erasure_Prob import Prob_loss_sequence
from Crossover_Prob import Optimal_Allocation_InnerCode
#==================Methods====================
#==================Parameter==================
# Coding parameter
# Msg_size =
# block_length_inner =  320
block_length_outer = 7000 # Block length of outer code

# Parameter
#*************************************************#
d_seq = 6 # Sequencing depth
r_total_max = 120 # maximum redundancy for inner code
r_total_min = 10 # minmuim redundancy for inner code
r_step = 3 # redundancy increment in this simulation
# if d_seq == 5:
#     PE = [0.1, 0.12, 0.14, 0.16] # Overall error probability (For d_seq = 5)
# elif d_seq == 10:
#     PE = [0.2] # Overall error probability (For d_seq = 10)
# elif d_seq == 15:
#     PE = [0.12, 0.14, 0.16, 0.18, 0.20] # Overall error probability (For d_seq = 15)
PE = [0.1]
#*************************************************#
n_syn = 30 # Synthesis number
N_pcr = 2 # PCR cycle number~
P_pcr = 0.9 # PCR efficiency
N_c = np.ceil((2*P_pcr)**(N_pcr)) # The copy times during PCR
# Coding parameter
K = (14,320)
ber = 1e-6 # Maximal block error rate

# Data file path
#========================Result data: reading cost vs R_i, R_o for different sequencing depth curve========================
numeric_result_path = f"D:\DeSP-main\Data\Cost_Optimization_result\Fix_indexing_cost\Data\Theory\BER_1e-6\Analysis-d{d_seq}.csv"


##=====================Creat file for simulation result=====================
# Creat a CSV file for simulation result
header = [['Seq_depth', 'PE', 'R_i', 'R_o', 'r_1', 'r_2', 't_1', 't_2', 'Reading cost']]  # Header for the .csv file.
'''
with open(numeric_result_path, 'w', newline='') as file:
    writer = csv.writer(file)
    # Write the list of lists to the CSV file
    writer.writerows(header)
'''
for P_e in PE:
    P_l = 0.43 * P_e  # Probability of sequence loss
    delta_2 = Prob_loss_sequence(n_syn, d_seq, N_pcr, P_pcr, P_l)  # Equivalent erasure probability for BEC
    print("#**************************************#",
        f"\nSequencing depth: {d_seq}, Pe: {P_e}, delta_2: {delta_2: .10f}.",
        "\n#**************************************#")
    for i, r_i in enumerate(range(r_total_min, r_total_max, r_step)):
        delta_1, r, T = Optimal_Allocation_InnerCode(P_e, d_seq, n_syn, N_c, r_i, K) # Equivalent crossover probability for BSC
        C = Calc_Capacity(delta_1, delta_2)  # Calculate capacity
        V = Calc_Dispersion(delta_1, delta_2, C)  # Calculate channel dispersion
        R_o = Rate_blocklength_NA(block_length_outer, C, V, ber)
        R_i = (K[0]+K[1])/(r_i + K[0]+K[1]) # Calculate the inner code rate according to the amount of redundancy for inner code.
        print(f"========================Exp {i+1}/{math.ceil((r_total_max - r_total_min) / r_step)}=========================")
        print(f"Inner code rate: {R_i: .6f}; ",f"Outer code rate: {R_o: .6f}.")
        print(f"delta_1: {delta_1: .6f}", f". Optimal r: ({int(r[0])}, {int(r[1])}); ", f"T: ({int(T[0])}, {int(T[1])}).")
        ## Calculating reading cost
        # Method 1
        reading_cost = d_seq/(R_o*R_i)*0.5
        print('Reading cost: ',reading_cost,' Bases/bit.')
        data = [d_seq, P_e, R_i, R_o, r[0], r[1], T[0], T[1], reading_cost]
        '''
        with open(numeric_result_path, mode='a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(data)
        print("Simulation Data has been written to .csv file in ", numeric_result_path)
        '''