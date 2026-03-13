import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob
import pandas as pd
import os
import logging
import warnings

from numpy.f2py.crackfortran import endifs

class Data2Figure:
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.CRITICAL)
    # plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['figure.dpi'] = 300
    # np.set_printoptions(threshold=np.inf)
    # %load_ext autoreload
    # %autoreload
    d_seq = 5
    ##======================Reading cost vs inner code rate======================##
    def CostVsInnerRate(numeric_result_path, Graphical_result_path):
        d_seq = [4, 5, 6, 7, 8]
        markers = ['^','o','s','p','h']
        colors = ['#EE5940','#F2AF30','#579B85','#2093AE','#A06C7D']
        plt.figure(figsize=(8, 6))
        for i, d in enumerate(d_seq):
            # Read data from .csv files
            df = pd.read_csv(numeric_result_path + f"\simulation-d{d}.csv") # Read
            plt.plot(df['R_i'], df['Reading cost'], label = '$d_{seq}$ = '+ str(df['Seq_depth'][0]), color = colors[i], marker=markers[i], markersize=3) # Plot BERs in one layer scheme.
        plt.legend()
        # plt.title('Influence of $𝑃_{𝑠𝑢𝑏}$ to Synthesis error')
        # plt.yscale('log')
        plt.ylim([4.5,7])
        plt.xlabel('$R_i$')
        plt.ylabel('Reading Cost / $bases/bit$')
        plt.grid(True, which='both', linestyle='-', linewidth=0.5)
        plt.savefig(Graphical_result_path)
    ##======================R_o vs R_i======================##
    def plotOuterRate_InnerRate(numeric_result_path, Graphical_result_path):
        markers = ['^','o','s','p','h']
        colors = ['#EE5940','#F2AF30','#579B85','#2093AE','#A06C7D']
        i = 0
        plt.figure(figsize=(8, 6))
        # Read data from .csv files
        df = pd.read_csv(numeric_result_path) # Read
        plt.plot(df['R_i'], df['R_o'], linewidth=2, linestyle='-', color='#EE5940', marker='o', markersize=6,
                     markerfacecolor='none') # Plot BERs in one layer scheme.
        # plt.title('Influence of $𝑃_{𝑠𝑢𝑏}$ to Synthesis error')
        # plt.yscale('log')
        plt.xlabel('$R_i$')
        plt.ylabel('$R_o$')
        plt.grid(True, which='both', linestyle=':', linewidth=0.5)
        plt.savefig(Graphical_result_path)


    ##======================R_o vs inner code redundancy allocation and error-correction capability======================##
    """
    markers = ['^','o','s','p','h']
    colors = ['#EE5940','#F2AF30','#579B85','#2093AE','#A06C7D']
    numeric_result_path = f".\Data\Cost_Optimization_result\Fix_indexing_cost\Data\Simulation-d{d_seq}\Simulation-2.csv"
    # plt.figure(figsize=(8, 6))
    # Read data from .csv files
    df = pd.read_csv(numeric_result_path) # Read
    
    # Create the figure and axis for redundancy.
    fig, ax1 = plt.subplots()
    
    ## Plot redundancy on ax1
    # Plot redundancy for the first layer
    ax1.plot(df['R_i'], df['r_1'], color=colors[0], label='$r_1$')
    ax1.set_xlabel('$R_i$')
    ax1.set_ylabel('Redundancy')
    ax1.tick_params(axis='y')
    
    # Plot redundancy for the second layer
    ax1.plot(df['R_i'], df['r_2'], color=colors[3], label='$r_2$')
    ax1.set_xlabel('$R_i$')
    ax1.set_ylabel('Redundancy')
    ax1.tick_params(axis='y')
    
    plt.grid(True, which='both', linestyle='-', linewidth=0.5)
    plt.legend(loc = 'upper left')
    
    
    # Create a second y-axis for error-correction capability
    ax2 = ax1.twinx()
    
    ## Plot error-correction capability on ax2
    # Plot error-correction capability for the first layer
    ax2.plot(df['R_i'], df['t_1'], color=colors[0], linestyle = '--', linewidth = 1, label='$t_1$')
    ax1.set_xlabel('$R_i$')
    ax2.set_ylabel('Error-Correction Capability')
    ax2.tick_params(axis='y')
    
    # Plot error-correction capability for the first layer
    ax2.plot(df['R_i'], df['t_2'], color=colors[3], linestyle = '--', linewidth = 1, label='$t_2$')
    ax1.set_xlabel('$R_i$')
    ax2.set_ylabel('Error-Correction Capability')
    ax2.tick_params(axis='y')
    
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend(loc = 'upper right')
    
    # plt.show()
    Graphical_result_path = f".\Data\Cost_Optimization_result\Fix_indexing_cost\Figure\Sequence-Depth_{d_seq}\Redundancy_t_vs_Ri"
    plt.savefig(Graphical_result_path)
    """

    ##======================Reading cost vs inner code rate for different sequencing depth (surface)======================##
    def plotReadCostSurface(numeric_result_path, Graphical_result_path):
        # Load the data from CSV.

        data = np.genfromtxt(numeric_result_path, delimiter=',')

        # Extract X, Y, and Z from the data
        # X: first column (excluding the first row)
        x = data[1:, 1]
        # Y: first row (excluding the first column)
        y = data[0, 2:]
        # Z: the remaining grid of values
        Z = data[1:, 2:]

        # Create a meshgrid.
        # Note: The shape of Z is (len(x), len(y)), so we use indexing='ij' to match the dimensions.
        X, Y = np.meshgrid(x, y, indexing='ij')

        # Find the minimum value and its coordinates, ignoring NaNs
        min_index = np.unravel_index(np.nanargmin(Z), Z.shape)
        x_min = X[min_index]
        y_min = Y[min_index]
        z_min = Z[min_index]

        # Plotting the surface
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot the transparent surface (alpha=0.5 for 50% transparency)
        surface = ax.plot_surface(X, Y, Z, cmap='coolwarm',alpha=0.6)

        # Highlight the optimal point with an "x" marker
        ax.scatter(x_min, y_min, z_min, color='black', marker='x', s=20, label='Optimal Point')

        # Add annotation
        annotation = f'Minimal Cost = {z_min:.2f}\n'+rf'($R_i$={x_min:.2f},$R_o$ = 0.72, $d_{{seq}}$={int(y_min)})'
        ax.text(x_min, y_min, z_min + 0.05 * (np.nanmax(Z) - np.nanmin(Z)), annotation,
                fontsize=14, color='black', ha='center', va='bottom')

        # Add labels and a color bar
        ax.set_xlabel('$R_i$', fontsize=12)
        ax.set_ylabel('Sequencing Depth', fontsize=12)
        ax.set_zlabel('Read Cost, $bases/bit$', fontsize=12)

        ax.view_init(elev=30, azim=150)  # adjust these values to get the desired perspective
        # fig.colorbar(surface, shrink=0.5, aspect=8)

        plt.savefig(Graphical_result_path, bbox_inches='tight', pad_inches=0.3)
        # plt.savefig("Graphical_result_path", bbox_inches='tight')

    ##======================Reading cost vs inner code rate for different sequencing depth (color plane)======================##
    """
    # Load the data from CSV.
    # We use np.genfromtxt which is robust to missing values.
    data_filepath = ".\Data\Cost_Optimization_result\Fix_indexing_cost\Data\Plot_Surface\simulation-summary.csv"
    data = np.genfromtxt(data_filepath, delimiter=',')
    
    # Extract X, Y, and Z from the data
    # X: first column (excluding the first row)
    x = data[0, 2:-2]
    # Y: first row (excluding the first column)
    y = data[5:-10:, 1]
    # Z: the remaining grid of values
    Z = data[5:-10, 2:-2]
    
    # Create a meshgrid.
    # Note: The shape of Z is (len(x), len(y)), so we use indexing='ij' to match the dimensions.
    X, Y = np.meshgrid(x, y)
    
    # Plotting the color plane
    plt.figure(figsize=(8, 6))
    plt.imshow(Z, aspect='auto', extent=[x.min(), x.max(), y.min(), y.max()], cmap='bone_r', origin='lower')
    
    # Add labels and a color bar
    plt.colorbar(label='Reading Cost')  # Add a color bar to the plot
    plt.xlabel('$R_i$')
    plt.ylabel('Sequencing Depth')
    # plt.show()
    
    Graphical_result_path = f".\Data\Cost_Optimization_result\Fix_indexing_cost\Figure\Plot_ColorPlane\simulation-summary"
    plt.savefig(Graphical_result_path)
    """
    """
    # Load the data from CSV.
    # We use np.genfromtxt which is robust to missing values.
    data_filepath = ".\Data\Cost_Optimization_result\Fix_indexing_cost\Data\Plot_Surface\simulation-summary.csv"
    data = np.genfromtxt(data_filepath, delimiter=',')
    
    # Extract X, Y, and Z from the data
    # X: first column (excluding the first row)
    x = data[0, 2:-2]
    # Y: first row (excluding the first column)
    y = data[5:-10:, 1]
    # Z: the remaining grid of values
    Z = data[5:-10, 2:-2]
    
    # Create a meshgrid.
    # Note: The shape of Z is (len(x), len(y)), so we use indexing='ij' to match the dimensions.
    X, Y = np.meshgrid(x, y)
    
    # Plotting the color plane
    plt.figure(figsize=(8, 6))
    plt.imshow(Z, aspect='auto', extent=[x.min(), x.max(), y.min(), y.max()], cmap='bone_r', origin='lower')
    
    # Add labels and a color bar
    plt.colorbar(label='Reading Cost')  # Add a color bar to the plot
    plt.xlabel('$R_i$')
    plt.ylabel('Sequencing Depth')
    # plt.show()
    
    Graphical_result_path = f".\Data\Cost_Optimization_result\Fix_indexing_cost\Figure\Plot_ColorPlane\simulation-summary"
    plt.savefig(Graphical_result_path)
    """

    ##======================Reading cost vs outer code rate to achieve 10^{-3} FER ======================##
    def plotReadCost_Simu(numeric_result_path, Graphical_result_path):
        # Define your prefix (e.g., all files starting with "data_")
        prefix = "Cost_Optimization_"
        Ro = [0.62, 0.67, 0.72, 0.77, 0.82]
        ReadingCost = []
        # Create the search pattern
        search_pattern = os.path.join(numeric_result_path, f"{prefix}*.csv")
        # Use glob to find matching files
        matching_files = glob.glob(search_pattern)
        # Check if any files matched

        if matching_files:
            for i in range(5):
                file_to_read = matching_files[i]  # Choose the first match (or loop over them)
                df = pd.read_csv(file_to_read)
                averageCost = df['Reading Cost'].mean()
                ReadingCost.append(averageCost)
        plt.figure(figsize=(8, 5))
        plt.plot(Ro, ReadingCost, linewidth=2, linestyle='-', color='#EE5940', marker='o', markersize=6,
                     markerfacecolor='none')  # Plot BERs in one layer scheme.
        plt.xlim([0.6, 0.84])
        plt.xlabel('$R_o$', fontsize=14)
        plt.ylabel('Read Cost, $bases/bit$', fontsize=14)
        plt.tick_params(axis="x", which="both", direction="in")
        plt.tick_params(axis="y", which="both", direction="in")
        plt.xticks(fontsize=14)  # Set font size for x-axis ticks
        plt.yticks(fontsize=14)  # Set font size for y-axis ticks
        plt.grid(True, which='both', linestyle=':', linewidth=0.1)
        plt.savefig(Graphical_result_path, bbox_inches='tight', format='eps')

    def plot_FER_vs_PE(numeric_result_path, Graphical_result_path = None):
        """
        Plot FER vs Pe figure from the simulation results stored in a .csv file.
        
        :param numeric_result_path: simulation results .csv file path
        :param Graphical_result_path: output figure file path
        """
        df = pd.read_csv(numeric_result_path)
        # ascending order sort by Pe
        df = df.sort_values(by='Pe')
        Pe_values = df['Pe']
        FER_values = df['FER']
        
        plt.figure(figsize=(8, 6))
        plt.semilogy(Pe_values, FER_values, marker='o', linestyle='-', color='#EE5940', markersize=6,
                     markerfacecolor='none', label = f'{CODE_ARY}-ary LDPC + TL-BCH, $d_{{seq}}=15$')
        plt.xlabel('Overall error probability, $P_e$', fontsize=14)
        plt.ylabel('Frame Error Rate (FER)', fontsize=14)
        plt.xticks(fontsize=14)  # Set font size for x-axis ticks
        plt.yticks(fontsize=14)  # Set font size for y-axis ticks
        plt.grid(True, which='both', linestyle=':', linewidth=0.5)
        plt.ylim([1e-5, 1e-1])
        plt.legend(fontsize=12)
        plt.savefig(Graphical_result_path, bbox_inches='tight', format='png')
        # plt.show()
    
    def plot_FER_vs_PE_multi_ary(result_dir, Graphical_result_path, sequencingDepth=10, innerRedundancy=114):
        """
        Plot FER vs Pe curves of multiple code aries (2/4/16) in one figure.
        """
        code_ary_list = [2, 4, 16]
        colors = ['#EE5940', '#2093AE', '#579B85']
        markers = ['o', 's', '^']

        plt.figure(figsize=(8, 6))

        for i, ary in enumerate(code_ary_list):
            file_path = os.path.join(
                result_dir,
                f"FFTQSPA_DNA_Channel_{ary}-ary_SequencingDepth{sequencingDepth}_InnerRedundancy{innerRedundancy}.csv"
            )
            if not os.path.exists(file_path):
                continue

            df = pd.read_csv(file_path).sort_values(by='Pe')
            plt.semilogy(
                df['Pe'],
                df['FER'],
                marker=markers[i],
                linestyle='-',
                color=colors[i],
                markersize=6,
                markerfacecolor='none',
                label=f'{ary}-ary LDPC'
            )

        plt.xlabel('Overall error probability, $P_e$', fontsize=14)
        plt.ylabel('Frame Error Rate (FER)', fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.grid(True, which='both', linestyle=':', linewidth=0.5)
        plt.ylim([1e-6, 0])
        plt.legend(fontsize=12)
        plt.savefig(Graphical_result_path, bbox_inches='tight', format='png')
#============================main====================================#
sequencingDepth = 10
innerRedundancy = 114 
CODE_LEN = 2048
Data2Figure.plot_FER_vs_PE_multi_ary(
    result_dir = f"results\codelen_{CODE_LEN}",
    Graphical_result_path = f"results\codelen_{CODE_LEN}\FER_vs_PE_multi_ary_SequencingDepth{sequencingDepth}_InnerRedundancy{innerRedundancy}.png",
    sequencingDepth = sequencingDepth, 
    innerRedundancy = innerRedundancy
)
