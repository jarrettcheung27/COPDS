import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import stats
from scipy.stats import norm
from scipy.integrate import quad
import math
import pandas as pd
from scipy.stats import binom
import math
import scipy.special as sp

# from Inbox.Sequencing_Cost_Optimization_Simu2 import delta_2


#==================Probability of loss sequence after voting==================
def Prob_loss_sequence(n_syn, d_seq, N_pcr,P_pcr, P_l):
    """
    Args:
        n_syn: Synthesis number
        d_seq: Sequence depth
        N_pcr: PCR cycle number
        P_pcr: PCR efficiency
        P_l: Probability of sequence loss
    return:
        Pr_erasure: Probability of loss sequence after voting.
    """
    Pr_erasure = 0 # Probability of loss sequence after voting.
    for reamin_num in range(0, n_syn + 1):
        ## The distribution of remaining copies for a given sequence in the oligo pool after Decay.
        Pr_reamin_num = binom.pmf(reamin_num, n_syn, 1 - P_l)

        ## Probability of zero sequencing reads in a cluster.
        N_c = math.ceil((2 * P_pcr) ** N_pcr) # Number of copies by PCR.
        P = d_seq / (n_syn * (1 - P_l) * N_c) # Probability of having a sequenceing read (sample rate).
        N = reamin_num * N_c # Total copies number of a particular sequence in the pool.
        # For a certain number N of total sequence, return the proability that sampled copies equal to 0.
        Pr_lambda = binom.pmf(0, N, P)
        Pr_erasure +=  Pr_lambda * Pr_reamin_num # Summer over all remaining copies value
    return Pr_erasure


#==================Channel Parameter==================
n_syn = 3 # Synthesis number
d_seq = [4, 5, 6, 7, 8] # Sequence depth
N_pcr = 2 # PCR cycle number
P_pcr = 0.9 # PCR efficiency
PE = np.linspace(0.02, 0.2, 10) # Overall error probability
delta_1 = []
for d in d_seq:
    row = []
    for P_e in PE:
        P_l = 0.43 * P_e  # Probability of sequence loss
        row.append(Prob_loss_sequence(n_syn, d, N_pcr,P_pcr, P_l))
    delta_1.append(row)

#===================Plot Figure======================
Graphical_result_path = "D:\DeSP-main\Data\Cost_Optimization_result\Fix_indexing_cost\Figure\Delta1_vs_PE\\3_Delta1.eps"
simu_file_name_base = "D:\DeSP-main\Data\Cost_Optimization_result\Fix_indexing_cost\Figure\Delta1_vs_PE\Ave_Erasure_Prob_simu_d_"
markers = ['^', 'o', 's', 'p', 'h']
colors = ['#EE5940', '#F2AF30', '#579B85', '#2093AE', '#A06C7D']
plt.figure(figsize=(8, 5))
for i, d in enumerate(d_seq):
    df = pd.read_csv(f"{simu_file_name_base}{d}.csv")
    plt.plot(df['Pe'], df['erasure_probability'], linewidth = 2,linestyle = '-',color = colors[i],marker = markers[i], markersize = 6, markerfacecolor='none', label = '$d_{seq}$ = '+f'{d}, simulated')  # simulated results 
    plt.plot(PE, delta_1[i], linewidth = 1,linestyle = '--',color = colors[i],marker = markers[i], markersize = 6, markerfacecolor='none', label='$d_{seq}$ = '+f'{d}, theoretical')  # theoretical results
plt.legend()
plt.legend(fontsize=10)
# plt.title('Influence of $𝑃_{𝑠𝑢𝑏}$ to Synthesis error')
plt.yscale('log')
# plt.ylim([4.5, 7])
plt.xlabel('$P_e$', fontsize = 14)
plt.ylabel('Erasure probability $\delta_1$ of BEC', fontsize = 14)
plt.tick_params(axis="x", which = "both", direction="in")
plt.tick_params(axis="y", which = "both", direction="in")
plt.xticks(fontsize=14)  # Set font size for x-axis ticks
plt.yticks(fontsize=14)  # Set font size for y-axis ticks
plt.grid(True, which='both', linestyle=':', linewidth=0.1)
plt.savefig(Graphical_result_path, bbox_inches='tight', format='eps')
plt.savefig(Graphical_result_path.replace('.eps', '.png'))
print("Figure saved successfully.")