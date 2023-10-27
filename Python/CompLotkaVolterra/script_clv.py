# Script to run clv on ts data with 10+1 taxa (+1 is called "other" and chosen as reference in alr)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.gridspec as gridspec
import math
import os
from datetime import datetime

from clv.compositional_lotka_volterra import *


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

    global T, Y, Names, U, P, Y_pc, log_Y, denom, denom_name, Names_alr

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
        mass = y.sum(axis=1)
        p = y / y.sum(axis=1,keepdims=True)
        p = (p + 1e-5) / (p + 1e-5).sum(axis=1,keepdims=True)
        P.append(p)
        Y_pc.append((mass.T*p.T).T)
        log_Y.append(np.log(mass.T*p.T).T)

    # plot time series over all taxa
    fig, ax = plt.subplots()
    for i in np.arange(n_taxa):
        ax.plot(T[0], P[0][:,i], label = f"x{i} ({Names[i]})")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f'{folderpath_out}plot_dataset.png', dpi=300, bbox_inches='tight')
    plt.close()


    ### Select the denominator for alr transformation
    if "other" in Names:
        # use the genus "other" as deonominator
        denom = Names.index("other") # for 10_most_abundant_rel_counts datasets
    else:
        denom = choose_denom(P)
    
    Names_alr = Names.copy()
    denom_name = Names_alr.pop(denom)

    ### plot alr transformed dataset
    # plot alr transformed time series
    ALR = construct_alr(P, denom)
    fig, ax = plt.subplots()
    for i in np.arange(n_taxa-1):
        ax.plot(T[0], ALR[0][:,i], label = Names_alr[i])
    ax.set_title(f'ALR transformed dataset - chosen denominator is \"{denom_name}\" (x{denom})')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f'{folderpath_out}plot_dataset_ALR.png', dpi=300, bbox_inches='tight')
    plt.close()

    # save ALR transformed data csv
    df_ALR = pd.DataFrame(data=ALR[0], columns=Names_alr)
    df_ALR.insert(0, "Time", T[0])
    df_ALR.to_csv(f'{folderpath_out}ts_ALR_denom-{denom}-{denom_name}_{data_name}.csv', index=False)

    # save P as csv
    df_P = pd.DataFrame(data=P[0], columns=Names)
    df_P.insert(0, "Time", T[0])
    df_P.to_csv(f'{folderpath_out}ts_clv_input_data_{data_name}.csv', index=False)

# Compositional Lotka Volterra
def apply_clv():

    global clv_elastic_net, clv_ridge
    
    ### clv
    print("train clv elastic net:")
    clv_elastic_net = CompositionalLotkaVolterra(P, T, denom=denom, pseudo_count=1e-5)
    clv_elastic_net.train()
    print(f"final regularizers: {clv_elastic_net.get_regularizers()}")

    ### clv with ridge regularization
    print("train clv ridge:")
    clv_ridge = CompositionalLotkaVolterra(P, T, denom=denom, pseudo_count=1e-5)
    clv_ridge.train_ridge()
    print(f"final regularizers: {clv_ridge.get_regularizers()}")


## Save and plot estimated g and A
def save_results(clv, regularization = "elastic_net"):
    
    # output variables clv
    A_clv, g_clv, B_clv = clv.get_params()

    g = round(pd.DataFrame(g_clv), 3)
    A = round(pd.DataFrame(A_clv), 3)

    # save g and A as csv file
    g.to_csv(f'{folderpath_out}clv_{regularization}_g.csv')
    A.to_csv(f'{folderpath_out}clv_{regularization}_A.csv')


    # prepare plotting
    max_value = max(A_clv.max(), g_clv.max(),
                    abs(A_clv.min()), abs(g_clv.min()))

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
    plt.savefig(f'{folderpath_out}clv_{regularization}_g_A_heatmaps.png', dpi=300, bbox_inches='tight')
    plt.close()


    ### Plot only heatmap of A
    # make heatmap
    plt.figure(figsize=(n_taxa, 0.75*n_taxa))
    sns.heatmap(A, fmt="", cmap="RdBu", center=0, vmin = -max_value, vmax = max_value, annot=True, yticklabels=False)
    plt.gca().xaxis.tick_top()
    plt.gca().tick_params(left=False, top=False)
    # save the plot
    plt.savefig(f'{folderpath_out}clv_{regularization}_A_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()


def predict_clv(clv, regularization = "elastic_net"):
    # load standard color palette
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    ### Predict time series using estimates from clv
    pred = clv.predict(p0 = np.array([P[0][0, :]]), times = T[0])
    # save pred as csv
    df_pred = pd.DataFrame(data=pred, columns=Names)
    df_pred.insert(0, "Time", T[0])
    df_pred.to_csv(f'{folderpath_out}ts_prediction_{data_name}.csv', index=False)

    n_col = 3
    n_row = math.ceil(n_taxa/3)

    fig, axs = plt.subplots(n_row, n_col)
    fig.suptitle(f"Fits for {data_name} with cLV ({regularization})")
    fig.set_figwidth(3*n_row)
    fig.set_figheight(4*n_col)
    fig.tight_layout(pad=2.0)

    for i in np.arange(n_taxa):
        # plot each taxon timeline separately
        axs[math.floor(i/n_col), (i%n_col)].plot(T[0], P[0][:,i], linewidth = 0.7, label = "data", color = colors[1])
        axs[math.floor(i/n_col), (i%n_col)].plot(T[0][1:], pred[1:,i], linewidth = 2, label = "fit")
        # axs[math.floor(i/n_col), (i%n_col)].plot(T[0], P[0][:,i], linewidth = 0.8, linestyle='--')
        # axs[math.floor(i/n_col), (i%n_col)].scatter(T[0], P[0][:,i], linewidth = 0.8, s = 1, color = colors[1])
        axs[math.floor(i/n_col), (i%n_col)].set_title(f"{Names[i]} (x{i})")
        axs[math.floor(i/n_col), (i%n_col)].title.set_size(10)
                        
    # plt.setp(axs, ylim=(min(pred[1:,].min(), P[0].min()) - 0.01, max(pred[1:,].max(), P[0].max()) + 0.05))
    # add overall legend below the plots
    handles, labels = [], []
    for ax in axs.ravel():
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h)
                labels.append(l)
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.95, 0.99), ncol = 2)

    # save the plot
    plt.savefig(f'{folderpath_out}clv_{regularization}_fits.png', dpi=300, bbox_inches='tight')
    plt.close()


def run_clv(data_name, filename, run):

    global folderpath_out

    filepath = folderpath_in + filename

    ### generate output directory
    
    # folderpath for output
    folderpath_out = f"C:/Users/Maria/Documents/Masterstudium/Masterarbeit/clv_output/output_{data_name}/output_{data_name}_run_{run}/"

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

    apply_clv()

    save_results(clv=clv_elastic_net)
    save_results(clv=clv_ridge, regularization="ridge")

    predict_clv(clv_elastic_net, regularization = "elastic_net")
    predict_clv(clv_ridge, regularization = "ridge")

    # redirect output to console
    sys.stdout = org_stdout


if __name__ == "__main__":

    folderpath_in = "C:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/explore/data/final_datasets/"

    # ### example input
    # data_name = "female"
    # # path of data file
    # subject = "female"
    # filename = f"ts_{subject}_Genus_10_most_abundant_rel_counts.csv"
    # # specify name of run
    # run = datetime.now().strftime('%m-%d_%H-%M')
    # # start run
    # run_clv(data_name, filename, run)

    compositional_data = [
        ['donorA', 'ts_donorA_Genus_10_most_abundant_rel_counts.csv'],
        ['donorB', 'ts_donorB_Genus_10_most_abundant_rel_counts.csv'],
        ['male', 'ts_male_Genus_10_most_abundant_rel_counts.csv'],
        ['female', 'ts_female_Genus_10_most_abundant_rel_counts.csv'],
        ['Silverman_all', 'ts_Silverman_Vall_all_Genus_10_most_abundant_rel_counts.csv'],
        # ['Silverman_daily', 'ts_Silverman_Vall_daily_Genus_10_most_abundant_rel_counts.csv'],
        ['Silverman_hourly', 'ts_Silverman_Vall_hourly_Genus_10_most_abundant_rel_counts.csv'],
        ['Bucci', 'ts_bucci_subject_all_rel_counts_denoised.csv']
    ]
    
    for run in ["{:02}".format(i) for i in range(1)]:
        print(f"run {run} starts")
        for data_name, filename in compositional_data:
            print(data_name, "started")
            run_clv(data_name, filename, run)
            print(data_name, "finished")
    