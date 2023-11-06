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

    # plot time series over all taxa
    fig, ax = plt.subplots()
    for i in np.arange(n_taxa):
        ax.plot(T[0], Y[0][:,i], label = f"x{i+1} ({Names[i]})")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f'{folderpath_out}plot_dataset_raw.png', dpi=300, bbox_inches='tight')
    plt.close()

    ### Prepare data: normalize and add pseudo counts to Y    
    Y_norm = []
    for y in Y:
        y_norm = (y - y.min()) / (y.max() - y.min() + 1e-8)
        y_norm = [y_norm + 1e-6]
        Y_norm.extend(y_norm)

    # P = []
    # Y_pc = []
    # log_Y = []
    # for y in Y_norm:
    #     mass = y.sum(axis=1)
    #     p = y / y.sum(axis=1,keepdims=True)
    #     p = (p + 1e-5) / (p + 1e-5).sum(axis=1,keepdims=True)
    #     P.append(p)
    #     Y_pc.append((mass.T*p.T).T)
    #     log_Y.append(np.log(mass.T*p.T).T)

    P = Y_norm

    # plot time series over all taxa
    fig, ax = plt.subplots()
    for i in np.arange(n_taxa):
        ax.plot(T[0], P[0][:,i], label = f"x{i+1} ({Names[i]})")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f'{folderpath_out}plot_dataset_prep.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # save P as csv
    df_P = pd.DataFrame(data=P[0], columns=Names)
    df_P.insert(0, "Time", T[0])
    df_P.to_csv(f'{folderpath_out}ts_glv_input_data_{data_name}.csv', index=False)


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
    g.to_csv(f'{folderpath_out}glv_{data_name}_{regularization}_g.csv')
    A.to_csv(f'{folderpath_out}glv_{data_name}_{regularization}_A.csv')

    # prepare plotting
    max_value = max(A_glv.max(), g_glv.max(),
                    abs(A_glv.min()), abs(g_glv.min()))

    ### Plot g and A in one graph

    # Creating figure and gridspec
    fig = plt.figure(figsize=(1.3*n_taxa, n_taxa))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, n_taxa])  

    # Creating subplots
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    # Plotting the vector g as a heatmap
    sns.heatmap(g, cmap="RdBu", vmin=-max_value, vmax=max_value, annot=True,
                cbar=False, ax=ax0)  # Set cbar=False to remove colorbar
    ax0.xaxis.tick_top()
    ax0.set_yticklabels([f'x{i+1}' for i in range(g.shape[0])], rotation=0)
    ax0.tick_params(left=False, top=False)

    # Plotting the matrix A as a heatmap
    sns.heatmap(A, cmap="RdBu", vmin=-max_value, vmax=max_value, annot=True, yticklabels=False, 
                cbar=False, ax=ax1) 
    ax1.xaxis.tick_top()
    ax1.tick_params(left=False, top=False)

    # set x labels
    ax0.set_xticklabels(['1'], fontsize = 10)
    ax1.set_xticklabels([f'x{i+1}' for i in range(A.shape[1])], fontsize = 10)

    # set title
    title_text = ax1.set_title(f"Estimated Interaction Effects for {data_name}", fontsize = 12, pad = 15)
    title_text.set_position((0.4, 1.1))
    
    # fig.tight_layout(pad=1.0)

    # Adding a colorbar for the entire figure
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # Adjust the position and size of the colorbar
    fig.colorbar(ax1.collections[0], cax=cbar_ax)

    # save the plot
    plt.savefig(f'{folderpath_out}glv_{data_name}_{regularization}_g_A_heatmaps.png', dpi=300, bbox_inches='tight')
    plt.close()


    ### Plot only heatmap of A
    # make heatmap
    plt.figure(figsize=(n_taxa, 0.75*n_taxa))
    sns.heatmap(A, fmt="", cmap="RdBu", center=0, vmin = -max_value, vmax = max_value, annot=True, yticklabels=False)
    plt.gca().xaxis.tick_top()
    plt.gca().tick_params(left=False, top=False)
    # save the plot
    plt.savefig(f'{folderpath_out}glv_{data_name}_{regularization}_A_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()


def predict_glv(glv, regularization = "elastic_net"):
    # load standard color palette
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    ### Predict time series using estimates from glv
    pred = glv.predict(c0 = np.array([P[0][0, :]]), times = T[0])

    # save pred as csv
    df_pred = pd.DataFrame(data=pred, columns=Names)
    df_pred.insert(0, "Time", T[0])
    df_pred.to_csv(f'{folderpath_out}ts_glv_prediction_{data_name}.csv', index=False)

    # plot prediction
    n_col = 2
    n_row = math.ceil(n_taxa/n_col)

    fig, axs = plt.subplots(n_row, n_col)
    fig.suptitle(f"Fits for {data_name} with gLV ({regularization})", fontsize = 16)
    fig.set_figwidth(4*n_col)
    fig.set_figheight(3*n_row)

    for i, ax in enumerate(axs.flat):
        if i >= n_taxa:
            ax.axis('off')  # Turn off the axis for empty plots
        else:
            # plot each taxon timeline separately
            ax.plot(T[0], P[0][:,i], linewidth = 0.7, label = "data", color = colors[1])
            ax.plot(T[0][1:], pred[1:,i], linewidth = 2, label = "fit")
            ax.set_title(f"x{i+1}", fontsize = 12)
            ax.title.set_size(10)

    # add overall legend below the plots
    handles, labels = [], []
    for ax in axs.ravel():
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h)
                labels.append(l)
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol = 2)

    fig.tight_layout(pad=1.0)

    # save the plot
    plt.savefig(f'{folderpath_out}glv_{data_name}_{regularization}_fits.png', dpi=300, bbox_inches='tight')
    plt.close()


def run_glv(data_name, filename, run):

    global folderpath_out

    filepath = folderpath_in + filename

    ### generate output directory
    
    # folderpath for output
    folderpath_out = f"C:/Users/Maria/Documents/Masterstudium/Masterarbeit/glv_output/output_{data_name}/output_{data_name}_run_{run}/"

    # create output folder
    if not os.path.exists(folderpath_out):
        os.makedirs(folderpath_out)
        
    # log output and errors
    log_file = open(f'{folderpath_out}log_{data_name}.log', 'w')
    # Redirect stdout to logfile
    org_stdout = sys.stdout
    sys.stdout = log_file
    
    ### load data
    load_data(filepath)

    apply_glv()

    save_results(glv=glv_elastic_net)
    save_results(glv=glv_ridge, regularization="ridge")

    predict_glv(glv_elastic_net, regularization = "elastic_net")
    predict_glv(glv_ridge, regularization = "ridge")

    # redirect output to console
    sys.stdout = org_stdout


if __name__ == "__main__":

    folderpath_in = "C:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/explore/data/final_datasets/"

    # ### example input
    # data_name = "miaSim"
    # # path of data file
    # filename = "miaSim_GLV_4species_new.csv"
    # # specify name of run
    # run = datetime.now().strftime('%m-%d_%H-%M')
    # # start run
    # run_glv(data_name, filename, run)

    datasets = [
        ['miaSim', 'miaSim_GLV_4species_new.csv'],
        ['miaSim_noise_0-005', 'ts_miaSim_GLV_4species_new_noise_0-005.csv'],
        ['miaSim_noise_0-01', 'ts_miaSim_GLV_4species_new_noise_0-01.csv'],
        ['miaSim_noise_0-02', 'ts_miaSim_GLV_4species_new_noise_0-02.csv'],
        ['miaSim_noise_0-04', 'ts_miaSim_GLV_4species_new_noise_0-04.csv'],
        ['3DLV', 'ts_3DLV.csv'],
        ['VanderPol', 'ts_VanderPol.csv'],
        ['VanderPol_noise_0-1', 'ts_VanderPol_noise_0-1.csv'],
        ['VanderPol_noise_0-2', 'ts_VanderPol_noise_0-2.csv'],
        ['VanderPol_noise_0-5', 'ts_VanderPol_noise_0-5.csv'],
        ['VanderPol_noise_1', 'ts_VanderPol_noise_1.csv']
    ]
    
    for run in ["{:02}".format(i) for i in range(1)]: # change if multiple runs should be started
        print(f"run {run} starts")
        for data_name, filename in datasets:
            print(data_name, "started")
            run_glv(data_name, filename, run)
            print(data_name, "finished")
    