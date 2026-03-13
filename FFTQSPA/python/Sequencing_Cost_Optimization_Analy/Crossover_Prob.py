import numpy as np
from scipy.stats import norm
from scipy.integrate import quad
from scipy.stats import binom
import math
import scipy.special as sp
import matplotlib.pyplot as plt
import pandas as pd

##=======pre define function===============
# Function to calculate t_1
def calculate_t(r, k):
    # Apply the equation for t_1
    t = np.floor(r / np.ceil(np.log2(k + r + 1)))
    return t
def calc_FER_of_index(code_length, P_e, t_1):
    '''
    Calculate the FER of index
    n: code length of index / bits.
    P_e: base edit prabability.
    t_1: error-correcting capability / bits

    prob: FER of index
    '''
    # Tranform code length from bit to base.
    M = np.ceil(code_length/2).astype(int)
    # Step 1: Calculate the threshold for the summation
    n_threshold = math.ceil((3/4) * t_1)

    # Step 2: Initialize the probability to 0
    prob = 0

    # Step 3: Perform the summation from n_threshold to r_n1 + M_1
    for n in range(n_threshold, M + 1):
        # Calculate the term P_e^n * (1 - P_e)^(r_n1 + M_1 - n)
        term = sp.comb(M, n) * (P_e ** n) * ((1 - P_e) ** (M - n))
        # Add this term to the total probability
        prob += term
    return prob if prob < 1 else 1 # In case, Prob > 1 due to precision issue.
# Define the function to be integrated
def integrand(x, n, P_e):

    norm_x = (x - (4 / 3) * n * P_e) / ((4 / 3) * np.sqrt(n * P_e * (1 - P_e)))
    # Standard normal PDF (phi(x))
    phi_x = norm.pdf(norm_x)

    return phi_x * x

# Function to calculate error bit in information bit
def Q_function(z):
    '''
    Q-function in probability theory.
    '''
    # Define the integrand for the Q function
    def integrand(y):
        return (1 / np.sqrt(2 * np.pi)) * np.exp(-y**2 / 2)

    # Perform the integration from z to infinity
    result, _ = quad(integrand, z, np.inf)
    return result

def h(delta):
    """
    Calculate the binary entropy function h(delta) for BSC.
    """
    if delta == 0 or delta == 1:
        return 0
    return -delta * np.log2(delta) - (1 - delta) * np.log2(1 - delta)
def func_cal_z(delta, n, k):
    """Calculate z based on the given formula."""
    # Calculate the entropy
    entropy = h(delta)

    # Calculate the term O(1) (here we will assume it's small enough to be ignored, or you can define a specific value)
    O_1 = 0  # You can replace this with an appropriate constant if needed

    # Calculate z
    numerator = n * (1 - entropy) - k + (1 / 2) * np.log2(n) + O_1
    denominator = np.sqrt(n * delta * (1 - delta)) * np.log2((1 - delta) / delta)

    return numerator / denominator
def pmf_lambda(lambda_, gamma, d_seq, n_syn, P_e, N_c):
    """
    Calculate the PMF of copies number lambda in each cluster after sampling and PCRing.
    Args:
        lambda_: copies number in each cluster
        gamma: remaining copies of reach sequence after decaying
        d_seq: sequencing depth
        n_syn: synthesis number
        c
    Return:
        The PMF of copies number lambda
    """
    sample_rate = d_seq / (n_syn * (1 - 0.43 * P_e)*N_c)
    return binom.pmf(lambda_, gamma*N_c, sample_rate)

def cal_ps_tilde(n_syn, d_seq, N_c, P_e, E_1, epsilon):
    """Calculate the total probability P(X > M_k / 2) considering K and Q as random variables."""
    ps_tilde = 0.0
    for gamma in range(0, n_syn + 1):
        # Calculate the probability of remaining copies gamma of each sequence due to decay.
        prob_gamma = binom.pmf(gamma, n_syn, 1 - 0.43 * P_e)
        for lambda_ in range(0, gamma + 1):
            prob_lambda = pmf_lambda(lambda_, gamma, d_seq, n_syn, P_e, N_c)
            for lambda_c in range(0, lambda_ + 1):
                # Calculate the number of sequences within a cluster with correct indices.
                prob_lambda_c = binom.pmf(lambda_c, lambda_, 1 - E_1)
                # Calculate the probability that less than λ_c/2 bit are correct in the original λ_c reads.
                prob_error_bit = binom.cdf(math.ceil(lambda_c/3), math.ceil(lambda_c), 1 - epsilon)*0.1
                # The bit error rate in the consensus sequence after voting, assumeing that half of the bits are different in the read with worng index in a cluster.(A bit error occurs when the number of error bits exceeds half of the total bits in the same column.)

                ps_tilde += prob_gamma * prob_lambda * prob_lambda_c * prob_error_bit
                #===========================Debuging===========================#
                """
                if prob_gamma * prob_lambda * prob_lambda_c * prob_error_bit > 0.01:
                    print("========================",
                          "\ngamma: ", gamma, ", prob_gamma: ", prob_gamma,
                          "\nlambda_: ", lambda_, ", prob_lambda_: ", prob_lambda,
                          "\nlambda_c: ", lambda_c,", prob_lambda_c: ", lambda_c,
                          "\nerror bit probability: ", prob_error_bit,
                          "\nps_tilde increment: ", prob_gamma * prob_lambda * prob_lambda_c * prob_error_bit)
                """
                #==============================================================#
    return ps_tilde

def Prob_crossover(P_e, d_seq, n_syn,N_c, K, r, T):
    """
    Args:
        P_e: Overall error probability
        n_syn: Synthesis number
        d_seq: Sequence depth
        K: length of index bits and information bits [index,inf].
        r: Redundancy distribution for inner code, [first layer, second layer].
        T: Error coorection capability for inner code, [first layer, second layer].
    return:
        delta_1: Probability of crossover probability after voting.
    """

    # Calculate the probability of error index.
    Ps = 0.57 * P_e
    FER_index = calc_FER_of_index(K[0]+r[0], Ps, T[0])

    # Calculate ber in information bits.
    delta = 0.38 * P_e # Equivalent coressover probability within a sequence.
    z = func_cal_z(delta, K[1]+r[1], K[1])
    block_error_rate_inf = Q_function(z) # ber in the informarion bits.
    # print("block error rate: ", block_error_rate_inf)
    Bit_error_rate = np.float64(1 - (1-block_error_rate_inf)**(1/(K[1]+r[1])))
    # Calculate the crossover probability after voting.
    # print("Bit error rate:", Bit_error_rate)
    delta_1 = cal_ps_tilde(n_syn, d_seq, N_c, P_e, FER_index, Bit_error_rate)
    return delta_1

def Optimal_Allocation_InnerCode(P_e, d_seq, n_syn, N_c, r_i, K):
    """
    Return the optimal redundancy allocation r for a given R_i and P_e.
    Args:
        R_i: inner code rate.
        P_e: P_e: Overall error probability.。
        K: length of index bits and information bits [index,inf].
        N_c: the copy time during PCR.
        r_i: total redundancy for inner code.
    Returns:
        delta_1_min: minimal crossover probability for inner coderatre = R_i.
    """
    # if r_i < 10:
    #     R_i_max = 1 / (10 / (K[0] + K[1]) + 1)
    #     raise ValueError(f"Code rate of inner code must be less than {R_i_max}")
    delta_1_min = 1 # Initial the value o f delta.
    for r_1 in range(r_i + 1):
        T = []  # Error coorection capability for inner code, [first layer, second layer].
        r = []  # The amount redundancy for inner code, [first layer, second layer].
        T.append(Correction_Capability(r_1, K[0])) # Error coorection capability of the first layer coding.
        r.append(r_1)
        r.append(r_i - r_1) # Amount of redundancy for the second layer.
        T.append(Correction_Capability(r[1], K[1])) # Error coorection capability of the second layer coding.
        delta_1 = Prob_crossover(P_e, d_seq, n_syn, N_c, K, r, T) # Calculate the corresponding crossover probability for every redundancy allocation.
        # print("delta_1 = ", delta_1,"\n")
        if delta_1 < delta_1_min:
            delta_1_min = delta_1 # Obtian the minimal crossover probability.
            r_opt = r # Obtian the optimal redundancy allocation.
            T_opt = T # Error-correction capability corresponding to the optimal redundancy allocation.
    return delta_1_min, r_opt, T_opt

def Correction_Capability(r, k):
    """
    Calculate the error-correction capability of BCH code for a given amount of redundancy and length of information bits.
    Args:
        r: amount of redundancy
        k:length of information bits
    Returns:
        t:error-correction capability
    """
    t = np.floor(r / np.ceil(np.log2(k + r + 1)))
    return t

# Simulation
'''
#==================Channel Parameter==================
n_syn = 3 # Synthesis number
d_seq = [4, 5, 6, 7, 8] # Sequence depth
N_pcr = 2 # PCR cycle number
P_pcr = 0.9 # PCR efficiency
N_c = np.ceil((2 * P_pcr) ** (N_pcr))  # The copy times during PCR
K = [14, 320] # Length of index bits and information bits [index,inf].  
r_i = 99 # Total redundancy for inner code.
PE = np.linspace(0.02, 0.2, 10) # Overall error probability
Crossover_probs = []
for d in d_seq:
    row = []
    for P_e in PE:
        P_l = 0.43 * P_e  # Probability of sequence loss
        delta_2, _, _ = Optimal_Allocation_InnerCode(P_e, d, n_syn, N_c, r_i, K) # Lowest Equivalent crossover probability for BSC
        row.append(delta_2)
    Crossover_probs.append(row)

#===================Plot Figure======================
Graphical_result_path = "D:\DeSP-main\Data\Cost_Optimization_result\Fix_indexing_cost\Figure\Delta2_vs_PE\\3_Delta2.eps"
simu_file_name_base = "D:\DeSP-main\Data\Cost_Optimization_result\Fix_indexing_cost\Figure\Delta1_vs_PE\Ave_Erasure_Prob_simu_d_"
markers = ['^', 'o', 's', 'p', 'h']
colors = ['#EE5940', '#F2AF30', '#579B85', '#2093AE', '#A06C7D']
plt.figure(figsize=(8, 5))
for i, d in enumerate(d_seq):
    df = pd.read_csv(f"{simu_file_name_base}{d}.csv")
    plt.plot(df['Pe'], df['erasure_probability'], linewidth = 2,linestyle = '-',color = colors[i],marker = markers[i], markersize = 6, markerfacecolor='none', label = '$d_{seq}$ = '+f'{d}, simulated')  # simulated results 
    plt.plot(PE, Crossover_probs[i], linewidth = 1,linestyle = '--',color = colors[i],marker = markers[i], markersize = 6, markerfacecolor='none', label='$d_{seq}$ = '+f'{d}, theoretical')  # theoretical results
plt.legend()
plt.legend(fontsize=10)
# plt.title('Influence of $𝑃_{𝑠𝑢𝑏}$ to Synthesis error')
plt.yscale('log')
# plt.ylim([4.5, 7])
plt.xlabel('$P_e$', fontsize = 14)
plt.ylabel('Crossover probability $\delta_2$ of BSC', fontsize = 14)
plt.tick_params(axis="x", which = "both", direction="in")
plt.tick_params(axis="y", which = "both", direction="in")
plt.xticks(fontsize=14)  # Set font size for x-axis ticks
plt.yticks(fontsize=14)  # Set font size for y-axis ticks
plt.grid(True, which='both', linestyle=':', linewidth=0.1)
plt.savefig(Graphical_result_path, bbox_inches='tight', format='eps')
plt.savefig(Graphical_result_path.replace('.eps', '.png'))
print("Figure saved successfully.")
'''




