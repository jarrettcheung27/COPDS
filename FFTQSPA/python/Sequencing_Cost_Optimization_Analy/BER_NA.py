import math
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
# import pandas as pd
import logging

import warnings

from Crossover_Prob import Optimal_Allocation_InnerCode
from Erasure_Prob import Prob_loss_sequence

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.CRITICAL)
# plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300
# np.set_printoptions(threshold=np.inf)
# %load_ext autoreload
# %autoreload
def Rate_blocklength_NA(block_length, capacity, dispersion, ber):
    """
    Calculate the rate for a given block length and bit error rate (BER) based on the given equation.

    Args:
        block_length (int): The block length (n).
        ber (float): block error rate (\epsilon).

    Returns:
        float: The calculated rate.
    """
    if block_length <= 0:
        raise ValueError("Block length must be greater than zero.")
    if not (0 < ber < 1):
        raise ValueError("Block error rate must be in the range (0, 1).")

    # Compute the Q-inverse of the bit error rate
    Q_inv = norm.ppf(1 - ber)

    # Calculate the main terms of the rate equation
    rate = capacity - Q_inv * math.sqrt(dispersion / block_length)

    # Add the O(log n) term (approximated here as negligible)
    # If a specific form for the O(log n) term is known, it should be added here.

    return rate

def Calc_Capacity(delta_1, delta_2):
    """
    Calculate the channel capacity based on the given equation.

    Args:
        delta_1 (float): Crossover probability for BSC.
        delta_2 (float): Erasure probability for BEC.

    Returns:
        float: The calculated channel capacity.
    """
    if not (0 <= delta_1 <= 1/2):
        raise ValueError(f"delta_1 = {delta_1} exceed the proper range for Normal Approximation, delta_1 must be in the range (0, 0.5).")
    if not (0 <= delta_2 <= 1):
        raise ValueError(f"delta_2 = {delta_2} exceed the proper range for Normal Approximation, delta_2 must be in the range [0, 1].")

    capacity = (1 - delta_2) * (1 - h(delta_1))
    return capacity

def Calc_Dispersion(delta_1, delta_2,C):
    """
    Calculate the channel dispersion based on the given equation.

    Args:
        delta_1 (float): Crossover probability for BSC.
        delta_2 (float): Erasure probability for BEC.
        C: Channel capacity

    Returns:
        float: The calculated channel dispersion.
    """
    if not (0 <= delta_1 <= 1/2):
        raise ValueError("delta_1 must be in the range (0, 0.5).")
    if not (0 <= delta_2 <= 1):
        raise ValueError("delta_2 must be in the range [0, 1].")

    term1 = (1 - delta_2) * (1 - 2 * h(delta_1))
    term2 = (1 - delta_2) * (1 - delta_1) * math.log2(1 - delta_1) ** 2
    term3 = (1 - delta_2) * delta_1 * math.log2(delta_1) ** 2

    dispersion = term1 + term2 + term3 - C ** 2
    return dispersion


def h(p):
    """
    Binary entropy function.

    Args:
        p (float): Probability value (0 <= p <= 1).

    Returns:
        float: Binary entropy of p.
    """
    if p <= 0 or p >= 1:
        return 0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)

#========================Plot the Rate-codelength curve========================
Graphical_result_path = "D:\DeSP-main\Data\Cost_Optimization_result\Fix_indexing_cost\Figure\Composite_Rate-Codelength\\3_Rate-blocklength_NA_BSC-BEC.eps"
bers = [1e-3,1e-4,1e-5,1e-6] # Target block error rates
delta_1 = 0.1
delta_2 = 0.1
C = Calc_Capacity(delta_1, delta_2)
V = Calc_Dispersion(delta_1, delta_2,C)
block_length = np.linspace(1,8000,100)
colors = ['#8590d5',"#8a8c94","#DD8F81BC","#F55235"]
plt.figure(figsize=(8, 5))
for i, ber in enumerate(bers):
    rate = [Rate_blocklength_NA(block_length[i], C, V, ber) for i in range(len(block_length))]
    plt.plot(block_length, rate , color = colors[i], linewidth = 1, label=f'Target block error rate = {ber:.0e}') # Plot the curve of rate-blocklength trade-off.
plt.legend(fontsize=10)
plt.ylim([0, 0.5])
plt.xlabel('Blocklength, n', fontsize = 14)
plt.ylabel('Rate, bit/ch.use', fontsize = 14)
plt.tick_params(axis="x", which = "both", direction="in")
plt.tick_params(axis="y", which = "both", direction="in")
plt.xticks(fontsize=12)  # Set font size for x-axis ticks
plt.yticks(fontsize=12)  # Set font size for y-axis ticks
plt.grid(True, which='both', linestyle=':', linewidth=0.1)
plt.savefig(Graphical_result_path, bbox_inches='tight', format='eps')
plt.savefig(Graphical_result_path.replace('.eps', '.png'))
print(f"Figure saved to {Graphical_result_path} successfully.")