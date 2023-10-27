# Script to run glv on ts data with 10+1 taxa (+1 is called "other" and chosen as reference in alr)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.gridspec as gridspec
import math
import os
from datetime import datetime

from clv.generalized_lotka_volterra import *


# function to import the datafile and put it into the right format
def create_data(filepath):
    data = pd.read_csv(f'{filepath}', sep=",", header=0)
    names = list(data.columns)[1:]
    usol = data.to_numpy()
    ts = usol[:, 0]
    data_y = usol[:, 1:]
    # set dimensions of the dataset
    global n_samples, n_taxa
    n_samples, n_taxa = data_y.shape
    print("The time coodinates have shape {}".format(ts.shape))
    print("The data has shape {}".format(data_y.shape)) 
    return [ts], [data_y], names

def load_data(filepath):

    global T, Y, Names, U, P, Y_pc, log_Y

    ### Read data
    T, Y, Names = create_data(filepath)

    # create empty matrix for external effects
    U = [ np.zeros((x.shape[0], 1)) for x in Y ]

    ### Prepare data
    # P then does not contain zero values
    P = []
    Y_pc = []
    log_Y = []
    for y in Y:
        p = (y + 1e-5)
        P.append(p)
        # mass = y.sum(axis=1)
        # p = y / y.sum(axis=1,keepdims=True)
        # p = (p + 1e-5) / (p + 1e-5).sum(axis=1,keepdims=True)
        # P.append(p)
        # Y_pc.append((mass.T*p.T).T)
        # log_Y.append(np.log(mass.T*p.T).T)

    # plot time series over all taxa
    fig, ax = plt.subplots()
    for i in np.arange(n_taxa):
        ax.plot(T[0], P[0][:,i], label = f"x{i} ({Names[i]})")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f'{folderpath_out}plot_dataset.png', dpi=300, bbox_inches='tight')
    plt.close()


# Generalized Lotka Volterra
def apply_glv():

    global glv_elastic_net, glv_ridge
    
    ### glv
    print("train glv elastic net:")
    glv_elastic_net = GeneralizedLotkaVolterra(P, T, pseudo_count=1e-5)
    glv_elastic_net.train()
    print(f"final regularizers: {glv_elastic_net.get_regularizers()}")

    ### glv with ridge regularization
    print("train glv ridge:")
    glv_ridge = GeneralizedLotkaVolterra(P, T, pseudo_count=1e-5)
    glv_ridge.train_ridge()
    print(f"final regularizers: {glv_ridge.get_regularizers()}")


## Save and plot estimated g and A
def save_results(glv, regularization = "elastic_net"):
    
    # output variables glv
    A_glv, g_glv, B_glv = glv.get_params()

    g = round(pd.DataFrame(g_glv), 3)
    A = round(pd.DataFrame(A_glv), 3)

    # save g and A as csv file
    g.to_csv(f'{folderpath_out}glv_{regularization}_g.csv')
    A.to_csv(f'{folderpath_out}glv_{regularization}_A.csv')


    # prepare plotting
    max_value = max(A_glv.max(), g_glv.max(),
                    abs(A_glv.min()), abs(g_glv.min()))

    ### Plot g and A in one graph

    # Creating figure and gridspec
    fig = plt.figure(figsize=(10, 7))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, n_taxa])  

    # Creating subplots
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    # Plotting the vector g as a heatmap
    sns.heatmap(g, cmap="RdBu", vmin=-max_value, vmax=max_value, annot=True, yticklabels=False, 
                cbar=False, ax=ax0)  # Set cbar=False to remove colorbar
    ax0.set_title('growth vector g')
    ax0.xaxis.tick_top()
    ax0.tick_params(left=False, top=False)

    # Plotting the matrix A as a heatmap
    sns.heatmap(A, cmap="RdBu", vmin=-max_value, vmax=max_value, annot=True, yticklabels=False, 
                cbar=False, ax=ax1) 
    ax1.set_title('interaction matrix A')
    ax1.xaxis.tick_top()
    ax1.tick_params(left=False, top=False)

    # Adding a colorbar for the entire figure
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # Adjust the position and size of the colorbar
    fig.colorbar(ax1.collections[0], cax=cbar_ax)

    # save the plot
    plt.savefig(f'{folderpath_out}glv_{regularization}_g_A_heatmaps.png', dpi=300, bbox_inches='tight')
    plt.close()


    ### Plot only heatmap of A
    # make heatmap
    plt.figure(figsize=(n_taxa, 0.75*n_taxa))
    sns.heatmap(A, fmt="", cmap="RdBu", center=0, vmin = -max_value, vmax = max_value, annot=True, yticklabels=False)
    plt.gca().xaxis.tick_top()
    plt.gca().tick_params(left=False, top=False)
    # save the plot
    plt.savefig(f'{folderpath_out}glv_{regularization}_A_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()


def predict_glv(glv, regularization = "elastic_net"):
    # load standard color palette
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    ### Predict time series using estimates from glv
    # c0 = [P[0][0, :]] * (n_taxa-1)
    print(Y[0][0, :])
    pred = glv.predict(c0 = np.array([Y[0][0, :]]), times = T[0])
    print(pred[0,:])

    n_row = n_taxa

    fig, axs = plt.subplots(n_row, 1)
    fig.set_figwidth(3)
    fig.set_figheight(3 * n_taxa)
    fig.tight_layout(pad=3.0)

    for i in np.arange(n_taxa):
        # plot each taxon timeline separately
        axs[i].plot(T[0], P[0][:,i], linewidth = 0.6, color = colors[1])
        axs[i].plot(T[0], pred[:,i], label = "Prediction")
        # axs[math.floor(i/n_col), (i%n_col)].plot(T[0], P[0][:,i], linewidth = 0.8, linestyle='--')
        # axs[math.floor(i/n_col), (i%n_col)].scatter(T[0], P[0][:,i], linewidth = 0.8, s = 1, color = colors[1])
        axs[i].set_title(f"{Names[i]} (x{i})")
        axs[i].title.set_size(10)
                        
    # plt.setp(axs, ylim=(min(pred[1:,].min(), P[0].min()) - 0.01, max(pred[1:,].max(), P[0].max()) + 0.05))

    # save the plot
    plt.savefig(f'{folderpath_out}glv_{regularization}_fits.png', dpi=300, bbox_inches='tight')
    plt.close()


def run_glv(data_name, filename, run):

    global folderpath_out

    filepath = folderpath_in + filename

    ### generate output directory
    
    # folderpath for output
    folderpath_out = f"C:/Users/Maria/Documents/Masterstudium/Masterarbeit/glv_output/output_{data_name}_run_{run}/"

    # create output folder
    if not os.path.exists(folderpath_out):
        os.makedirs(folderpath_out)
    
    # log output and errors
    log_file = open(f'{folderpath_out}log_file_{data_name}.log', 'w')
    # Redirect stdout and stderr
    sys.stdout = log_file
    sys.stderr = log_file
    
    ### load data
    load_data(filepath)

    apply_glv()

    save_results(glv=glv_elastic_net)
    save_results(glv=glv_ridge, regularization="ridge")

    predict_glv(glv_elastic_net, regularization = "elastic_net")
    predict_glv(glv_ridge, regularization = "ridge")


if __name__ == "__main__":

    folderpath_in = "C:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/explore/data/final_datasets/"

    ### example input
    data_name = "miaSim"
    # path of data file
    filename = "miaSim_GLV_4species_new.csv"
    # specify name of run
    run = datetime.now().strftime('%m-%d_%H-%M')
    # start run
    run_glv(data_name, filename, run)

    # datasets = [
    #     ['miaSim', 'miaSim_GLV_4species_new.csv'],
    #     [],
    #     [],
    #     [],
    #     ['3DLV', 'ts_3DLV.csv'],
    #     ['VanderPol', 'ts_VanderPol.csv'],
    #     ['VanderPol_noise_0-1', 'ts_VanderPol_noise_0-1.csv'],
    #     ['VanderPol_noise_0-2', 'ts_VanderPol_noise_0-2.csv'],
    #     ['VanderPol_noise_0-5', 'ts_VanderPol_noise_0-5.csv'],
    #     ['VanderPol_noise_1', 'ts_VanderPol_noise_1.csv']
    # ]
    
    # for run in ["00"]:
    #     print(f"run {run} starts")
    #     for data_name, filename in datasets:
    #         print(data_name, "started")
    #         run_glv(data_name, filename, run)
    #         print(data_name, "finished")
    