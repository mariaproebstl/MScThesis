# script to run deepmod for ODEs

# General imports
import numpy as np
import torch
import matplotlib.pylab as plt
import pandas as pd
import os
import sys
import logging
import tensorflow as tf
from tensorflow.core.util import event_pb2
import shutil
from datetime import datetime
import argparse

# DeepMoD functions
from deepymod import DeepMoD
from deepymod.data import Dataset, get_train_test_loader
from deepymod.model.func_approx import NN
from deepymod.model.constraint import LeastSquares
from deepymod.model.sparse_estimators import Threshold, PDEFIND
from deepymod.training import train
from deepymod.training.sparsity_scheduler import TrainTestPeriodic
from deepymod.model.libraryODE import LibraryODE

# clv functions
from compositional_lotka_volterra import choose_denom
from compositional_lotka_volterra import construct_alr

# # Settings for reproducibility
# np.random.seed(30)
# torch.manual_seed(0)

# Configuring GPU or CPU
if False: #torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"

################################## Commands #######################################
# C:/Users/Maria/anaconda3/envs/DeePyMoD/python.exe c:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/Python/DeePyMoD/script_deepmod_ODE.py 
# +
# # human ts 
# -data_name "humanTS_donorA" -filename 'ts_donorA_Phylumlevel_rel_count.csv' -hl_number 15 -hl_size 60 -max_iterations 50000
# -data_name "humanTS_donorB" -filename 'ts_donorB_Phylumlevel_rel_count.csv' -hl_number 15 -hl_size 60 -max_iterations 50000
# -data_name "humanTS_male" -filename 'ts_male_Phylumlevel_rel_count.csv' -hl_number 15 -hl_size 60 -max_iterations 50000
# -data_name "humanTS_female" -filename 'ts_female_Phylumlevel_rel_count.csv' -hl_number 15 -hl_size 60 -max_iterations 50000
# # bucci
# -data_name "ts_bucci_subject_1_run$i" -filename "ts_bucci_subject_1_rel_count.csv" -hl_number 10 -hl_size 60 -max_iterations 50000
# -data_name "ts_bucci_subject_2_run$i" -filename "ts_bucci_subject_2_rel_count.csv" -hl_number 10 -hl_size 60 -max_iterations 50000
# -data_name "ts_bucci_subject_3_run$i" -filename "ts_bucci_subject_3_rel_count.csv" -hl_number 10 -hl_size 60 -max_iterations 50000
# -data_name "ts_bucci_subject_4_run$i" -filename "ts_bucci_subject_4_rel_count.csv" -hl_number 10 -hl_size 60 -max_iterations 50000
# -data_name "ts_bucci_subject_5_run$i" -filename "ts_bucci_subject_5_rel_count.csv" -hl_number 10 -hl_size 60 -max_iterations 50000
# # miaSim S4
# -data_name "ts_miaSim_S4_run$i" -filename "miaSim_GLV_4species_oscillating_zero.csv" -hl_number 10 -hl_size 60 -max_iterations 50000


# C:/Users/Maria/anaconda3/envs/DeePyMoD/python.exe c:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/Python/DeePyMoD/script_deepmod_ODE.py -data_name "ts_miaSim_S4_test" -filename "miaSim_GLV_4species_oscillating_zero.csv" -hl_number 5 -hl_size 30 -max_iterations 100


################################### Variables ######################################

# help functions to parse arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-data_name', required = True)
    parser.add_argument('-filename', required = True)
    parser.add_argument('-max_samples', type = int)
    parser.add_argument('-int_order', type = int, default = 2)
    parser.add_argument('-hl_number', type = int, default = 5)
    parser.add_argument('-hl_size', type = int, default = 40)
    parser.add_argument('-add_alr', default = False)
    parser.add_argument('-threshold', type = float, default = 0.1)
    parser.add_argument('-max_iterations', type = int, default = 100000)
    parser.add_argument('-n_runs', type = int, default = 1)
    
    args = parser.parse_args()
    return args

def initialize_args(args):
    global data_name, filename, max_samples, int_order, hl_number, hl_size, add_alr, threshold, max_iterations, n_runs
    
    data_name = f"{args.data_name}_{datetime.now().strftime('%m-%d_%H-%M')}"
    filename = args.filename
    max_samples = args.max_samples
    int_order = args.int_order
    hl_number = args.hl_number
    hl_size = args.hl_size
    add_alr = args.add_alr
    threshold = args.threshold
    max_iterations = args.max_iterations
    n_runs = args.n_runs

    # # some example files/names
    # name = "ts_bucci_subject_1"
    # filename = "ts_bucci_subject_1_rel_count.csv"
    # name = "ts_miaSim_S4"
    # filename = "miaSim_GLV_4species_oscillating_zero.csv"
    # name = "humanTS_donorA"
    # filename = "ts_donorA_Phylumlevel_rel_count.csv"

    # name = "ts_bucci_subject_1"
    # data_name = f"{name}_{datetime.now().strftime('%m-%d_%H-%M')}"
    # filename = "ts_bucci_subject_1_rel_count.csv"
    # # max_samples = 100
    # int_order = 2
    # hl_number = 10
    # hl_size = 60
    # add_alr = True
    # threshold = 0.1
    # max_iterations = 50000
    # n_runs = 1
    
    
# function to initialize all variables
def set_variables():
    global filepath, folderpath_out, folderpath_plots, folderpath_data, write_iterations

    args = parse_args()
    initialize_args(args) # = None

    # specify how often data is written to tensorboard and checks train loss , by default 25.
    write_iterations = 25
    
    # folderpaths for output
    folderpath_out = f"../output/output_{data_name}"
    folderpath_plots = f'{folderpath_out}/Plots'
    folderpath_data = f'{folderpath_out}/Data'

    # create output folder
    os.makedirs(folderpath_out)
    os.makedirs(folderpath_plots)
    os.makedirs(folderpath_data)
    
    # initialize log file
    log_filename = f'{folderpath_out}/log_{data_name}.log'
    logging.basicConfig(filename=log_filename,
                        filemode='a',
                        format='%(asctime)s, %(message)s',
                        level=logging.INFO)
    
    # path of data file (input)
    filepath = "../../data/" + filename
    
    
    logging.info(f"""the parameters are initialized for {data_name}:\n
                ALR transformation: {add_alr}\n
                threshold: {threshold}\n
                hidden layers: number={hl_number}, size={hl_size}\n
                order of interactions: {int_order} \n
                max. iterations: {max_iterations}\n
                number of runs: {n_runs}\n
                \n
                device = {device}""")


# help function to add alr transformation
def add_alr_transformation(T, P):
    denom = choose_denom(P)

    ALR = construct_alr(P, denom)
    ALR

    # plot alr
    fig, ax = plt.subplots()
    for i in np.arange(n_otu-1):
        ax.plot(T, ALR[0][:, i])  # , label = f"x{i+1}"
    ax.set_title(f"chosen denominator is otu {denom+1}")
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f'{folderpath_plots}/plot_dataset_alr.png', dpi = 200) #, bbox_inches='tight'
    plt.close()

    return ALR[0]

# function to import the datafile and put it into the right format
def create_data():
    data = np.genfromtxt(filepath, delimiter=",")
    usol = data[1:, :]  # removes header
    if "max_samples" in globals():
        ts = usol[0:max_samples, 0]
        data_y = usol[0:max_samples, 1:]
    else:
        ts = usol[:, 0]
        data_y = usol[:, 1:]

    # set dimensions of the dataset
    global n_samples, n_otu
    n_samples, n_otu = data_y.shape

    # plot the raw data
    fig, ax = plt.subplots()
    for i in np.arange(n_otu):
        ax.plot(ts, data_y[:, i], label=f"otu {i+1}")
    ax.set_xlabel("Time")
    ax.set_ylabel("Abundance")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f'{folderpath_plots}/plot_dataset.png', 
                bbox_inches='tight', dpi = 200)
    plt.close()

    # include alr transformation (and save plot)
    if add_alr:
        data_y = add_alr_transformation(ts, [data_y])
        n_otu = n_otu-1

    T = torch.from_numpy(ts.reshape(-1, 1)).float()
    Y = torch.from_numpy(data_y).float()

    return T, Y


################################ Configuring DeepMoD ##########################################

# help function
def access_TFRecordDataset(out_var, log_path):

    out_var_dir = log_path + "/" + out_var + "/"
    tmp_file = os.listdir(out_var_dir)[-1]
    file_dir = out_var_dir + tmp_file

    out = np.array([])
    index = np.array([])
    i = 0
    for serialized_example in tf.data.TFRecordDataset(file_dir):
        event = event_pb2.Event.FromString(serialized_example.numpy())
        for value in event.summary.value:
            # Extract relevant information from the event
            val = value.simple_value
            out = np.append(out, val)
            index = np.append(index, (i+1)*write_iterations)
            i += 1
    
    # save values
    df_tmp = pd.DataFrame({'Iteration': index, 'Value': out})
    df_tmp.to_csv(f"{log_path}/data/{out_var}.csv", index=False)
    
    return [index, out]


def run_deepmod_and_save_results(dataset, network_shape):

    # Split dataset
    train_dataloader, test_dataloader = get_train_test_loader(
        dataset, train_test_split=0.8)

    # network
    hidden_layer = list(np.repeat(network_shape[0], network_shape[1]))
    network = NN(1, hidden_layer, n_otu)

    # library function
    library = LibraryODE(int_order=int_order, intercept=False)

    # Configuration of the sparsity estimator and sparsity scheduler
    estimator = Threshold(threshold)
    sparsity_scheduler = TrainTestPeriodic(
        periodicity=100, patience=200, delta=1e-5)

    constraint = LeastSquares()
    # instantiate the model
    model = DeepMoD(network, library, estimator, constraint)  # .to(device)
    # define optimizer
    optimizer = torch.optim.Adam(
        model.parameters(), betas=(0.99, 0.99), amsgrad=True, lr=5e-3
    )

    # Run DeepMoD

    # create directory for train output
    log_path = f'{folderpath_out}/train_log'
    if not os.path.exists(log_path):
        os.makedirs(log_path)
        os.makedirs(f'{log_path}/Plots/')
        os.makedirs(f'{log_path}/Data/')

    # log print output of train()
    old_stdout = sys.stdout
    log_file = open(f"{folderpath_out}/log_iterations.log", "w")
    sys.stdout = log_file

    train(
        model,
        train_dataloader,
        test_dataloader,
        optimizer,
        sparsity_scheduler,
        log_dir=log_path,
        max_iterations=max_iterations,
        write_iterations=write_iterations,
        delta=1e-2,
        sparsity_update = False
    )

    # close log file again
    sys.stdout = old_stdout
    log_file.close()

    ####################### save results ########################

    # library
    # save structure of the library (list of coefficients contained in the library)
    library_values = model.library.get_content(dataset.data)
    df_library_values = pd.DataFrame()
    idx = 0
    for ls in library_values:
        df_tmp = pd.DataFrame(ls)
        df_library_values[f"x{idx+1}"] = df_tmp
        idx += 1

    df_library_values.to_csv(
        f"{folderpath_data}/model_library_values.csv")

    # number of coefficients per otu
    n_coefs = len(library_values[0])

    # save sparsity mask and estimated coefficients
    # sparsity masks
    df_sparsity_masks = pd.DataFrame()
    idx = 0
    for ls in model.sparsity_masks:
        np_tmp = ls.numpy()
        df_tmp = pd.DataFrame(np_tmp)
        df_sparsity_masks[f"x{idx+1}"] = df_tmp
        idx += 1

    df_sparsity_masks.to_csv(
        f"{folderpath_data}/model_sparsity_masks.csv")

    # estimation coefficients
    df_estimated_coeffs = pd.DataFrame()
    idx = 0
    for ls in model.estimator_coeffs():
        df_tmp = pd.DataFrame(ls)
        df_estimated_coeffs[f"x{idx+1}"] = df_tmp
        idx += 1

    df_estimated_coeffs.to_csv(
        f"{folderpath_data}/model_estimated_coeffs.csv")

    # Analysis/Visualization of the loss

    # # get list of all output values that were calculated during train()
    # os.listdir(log_path)
    # loss_vars = ["loss_mse_output_0", "remaining_MSE_test_val_0",
    #               "loss_l1_output_0", "loss_reg_output_0"]

    for otu_tmp in np.arange(n_otu):

        loss_mse = access_TFRecordDataset(
            f"loss_mse_output_{otu_tmp}", log_path)
        loss_reg = access_TFRecordDataset(
            f"loss_reg_output_{otu_tmp}", log_path)
        MSE_test = access_TFRecordDataset(
            f"remaining_MSE_test_val_{otu_tmp}", log_path)
        loss_l1 = access_TFRecordDataset(f"loss_l1_output_{otu_tmp}", log_path)

        # plot mse and reg loss
        fig, ax = plt.subplots()
        ax.plot(loss_mse[0], loss_mse[1],
                c='#002635', marker='o', label='MSE loss')
        ax.plot(loss_reg[0], loss_reg[1],
                c='gray', marker='o', ls='--', alpha=0.6, label='Reg loss')
        ax.set_yscale('log')
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Cost')
        ax.set_title(f'loss for x{otu_tmp+1}')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(f'{folderpath_plots}/loss_plot_x{otu_tmp+1}.png',
            bbox_inches='tight', dpi = 200)
        plt.close()

        # plot and save estimated coefs by iteration
        output = []

        for coef in np.arange(n_coefs):
            output_coef = access_TFRecordDataset(
                f"estimator_coeffs_output_{otu_tmp}_coeff_{coef}", log_path)
            output.append(output_coef)

        fig, ax = plt.subplots()
        for coef in np.arange(n_coefs):
            ax.scatter(output[coef][0], output[coef][1], 
                       label=f'{library_values[otu_tmp][coef]}', s=1)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Coefficient')
        ax.set_title(f'Coefficients for x{otu_tmp+1}')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(
            f'{folderpath_plots}/estimated_coeffs_x{otu_tmp+1}.png',
            bbox_inches='tight', dpi = 200)
        plt.close()

    
    # # move prediction plots and values to parent folder
    # shutil.move(f"{log_path}/Plots/",
    #             folderpath_out)
    # shutil.move(f"{log_path}/Data/",
    #             folderpath_out)
    # fetch all Plots
    for file_name in os.listdir(f"{log_path}/Plots/"):
        # construct full file path
        source = f"{log_path}/Plots/" + file_name
        destination = folderpath_plots + "/" + file_name
        shutil.move(source, destination)
    # fetch all data files
    for file_name in os.listdir(f"{log_path}/Data/"):
        # construct full file path
        source = f"{log_path}/Data/" + file_name
        destination = folderpath_data + "/" + file_name
        shutil.move(source, destination)
    
    # remove log folder with all training output
    shutil.rmtree(log_path)


if __name__ == "__main__":
        
    set_variables()
    
    # load dataset
    dataset = Dataset(
        create_data,
        device=device,
    )
    logging.info("preparation of the dataset is done")

    model = run_deepmod_and_save_results(dataset, network_shape=[hl_size, hl_number])

    logging.info("run finished")
