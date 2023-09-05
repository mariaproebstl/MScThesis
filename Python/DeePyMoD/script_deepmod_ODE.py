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
# # human ts 
# python3 script_deepmod_ODE.py -data_name "humanTS_donorA" -filename 'ts_donorA_Phylumlevel_rel_count.csv' -hl_number 15 -hl_size 60 -max_iterations $iterations -threshold 0.1
# python3 script_deepmod_ODE.py -data_name "humanTS_donorB" -filename 'ts_donorB_Phylumlevel_rel_count.csv' -hl_number 15 -hl_size 60 -max_iterations $iterations -threshold 0.1
# python3 script_deepmod_ODE.py -data_name "humanTS_male" -filename 'ts_male_Phylumlevel_rel_count.csv' -hl_number 15 -hl_size 60 -max_iterations $iterations -threshold 0.1
# python3 script_deepmod_ODE.py -data_name "humanTS_female" -filename 'ts_female_Phylumlevel_rel_count.csv' -hl_number 15 -hl_size 60 -max_iterations 50000
# # bucci
# python3 script_deepmod_ODE.py -data_name "ts_bucci_subject_1_run$i" -filename "ts_bucci_subject_1_rel_count.csv" -hl_number 10 -hl_size 60 -max_iterations 50000
# python3 script_deepmod_ODE.py -data_name "ts_bucci_subject_2_run$i" -filename "ts_bucci_subject_2_rel_count.csv" -hl_number 10 -hl_size 60 -max_iterations 50000
# python3 script_deepmod_ODE.py -data_name "ts_bucci_subject_3_run$i" -filename "ts_bucci_subject_3_rel_count.csv" -hl_number 10 -hl_size 60 -max_iterations 50000
# python3 script_deepmod_ODE.py -data_name "ts_bucci_subject_4_run$i" -filename "ts_bucci_subject_4_rel_count.csv" -hl_number 10 -hl_size 60 -max_iterations 50000
# python3 script_deepmod_ODE.py -data_name "ts_bucci_subject_5_run$i" -filename "ts_bucci_subject_5_rel_count.csv" -hl_number 10 -hl_size 60 -max_iterations 50000
# # miaSim S4
# python3 script_deepmod_ODE.py -data_name "ts_miaSim_S4_run$i" -filename "miaSim_GLV_4species_oscillating_zero.csv" -hl_number 10 -hl_size 60 -max_iterations 50000


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
    
    
# function to initialize all variables
def set_variables():
    global filepath, folderpath_out
    
    # # name of the dataset including current date_time
    # data_name = f"humanTS_female_{datetime.now().strftime('%m-%d_%H-%M')}"
    # 
    # # filename = "ts_bucci_subject_5_rel_count.csv"
    # # filename = "miaSim_GLV_4species_oscillating_zero.csv"
    # filename = "ts_female_Phylumlevel_rel_count.csv"
    # 
    # # add alr transformation to the data
    # add_alr = False
    # 
    # ### specify the model parameters for training
    # 
    # # order of interactions included in the library (2 or 3)
    # int_order = 2
    # 
    # # threshold
    # threshold = 0.5
    # 
    # # network
    # # number of hidden layers
    # hl_number = 15
    # # size of hidden layers
    # hl_size = 50
    # 
    # # maximal iterations
    # max_iterations = 100000
    # 
    # # number of runs
    # n_runs = 2
    
    
    # python3 script_deepmod_ODE.py 
    # -data_name 'humanTS_female' -filename 'ts_female_Phylumlevel_rel_count.csv'
    # -hl_number 5 -hl_size 50 -n_runs 3

    args = parse_args()
    initialize_args(args)
    
    
    # folderpath output
    folderpath_out = f"../output/output_{data_name}"
    # create output folder
    if not os.path.exists(folderpath_out):
        os.makedirs(folderpath_out)
    
    # initialize log file
    log_filename = f'{folderpath_out}/log_{data_name}.log'
    logging.basicConfig(filename=log_filename,
                        filemode='a',
                        format='%(asctime)s, %(message)s',
                        level=logging.INFO)
    
    # path of data file (input)
    filepath = "../data/" + filename
    
    
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
    plt.savefig(f'{folderpath_out}/plot_dataset_alr.png')
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
    ax.legend()
    plt.savefig(f'{folderpath_out}/plot_dataset.png')
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

    res = np.array([])
    i = 0
    for serialized_example in tf.data.TFRecordDataset(file_dir):
        event = event_pb2.Event.FromString(serialized_example.numpy())
        for value in event.summary.value:
            # Extract relevant information from the event
            val = value.simple_value
            res = np.append(res, val)
            i += 1

    return res


def run_deepmod_and_save_results(iter, dataset, folderpath_out, network_shape):

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
    log_path = f'{folderpath_out}/train_log_{iter}'
    if not os.path.exists(log_path):
        os.makedirs(log_path)

    # log print output of train()
    old_stdout = sys.stdout
    log_file = open(f"{folderpath_out}/log_iteration{iter}.log", "w")
    sys.stdout = log_file

    train(
        model,
        train_dataloader,
        test_dataloader,
        optimizer,
        sparsity_scheduler,
        log_dir=log_path,
        max_iterations=max_iterations,
        delta=1e-3,
    )

    # close log file again
    sys.stdout = old_stdout
    log_file.close()

    # move prediction plots to parent folder
    for i in np.arange(n_otu):
        shutil.move(f"{log_path}/prediction_x{i+1}.png",
                    f"{folderpath_out}/prediction_x{i+1}.png")

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
        f"{folderpath_out}/model_library_values_{iter}.csv")

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
        f"{folderpath_out}/model_sparsity_masks_{iter}.csv")

    # estimation coefficients
    df_estimated_coeffs = pd.DataFrame()
    idx = 0
    for ls in model.estimator_coeffs():
        df_tmp = pd.DataFrame(ls)
        df_estimated_coeffs[f"x{idx+1}"] = df_tmp
        idx += 1

    df_estimated_coeffs.to_csv(
        f"{folderpath_out}/model_estimated_coeffs_{iter}.csv")

    # Analysis/Visualization of the loss

    # # get list of all output values that were calculated during train()
    # os.listdir(log_path)
    # loss_vars = ["loss_mse_output_0", "remaining_MSE_test_val_0",
    #               "loss_l1_output_0", "loss_reg_output_0"]

    for otu_tmp in np.arange(n_otu):

        train_loss = access_TFRecordDataset(
            f"loss_mse_output_{otu_tmp}", log_path)
        test_loss = access_TFRecordDataset(
            f"remaining_MSE_test_val_{otu_tmp}", log_path)
        l1_loss = access_TFRecordDataset(f"loss_l1_output_{otu_tmp}", log_path)
        reg_loss = access_TFRecordDataset(
            f"loss_reg_output_{otu_tmp}", log_path)

        # plot and save loss
        fig, ax = plt.subplots()
        ax.plot(train_loss, c='#002635', marker='o', label='Train loss')
        ax.plot(test_loss, c='red', marker='o',
                ls='--', alpha=0.6, label='Test loss')
        ax.plot(reg_loss, c='gray', marker='o',
                ls='--', alpha=0.6, label='Reg loss')
        ax.plot(l1_loss, c='darkgray', marker='o',
                ls='--', alpha=0.6, label='L1 loss')
        ax.set_yscale('log')
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Cost')
        ax.set_title(f'loss for x{otu_tmp+1}')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(f'{folderpath_out}/loss_plot_x{otu_tmp+1}_iter{iter}.png',
            bbox_inches='tight')
        plt.close()

        # plot and save estimated coefs by iteration
        output = []

        for coef in np.arange(n_coefs):
            output_coef = access_TFRecordDataset(
                f"estimator_coeffs_output_{otu_tmp}_coeff_{coef}", log_path)
            output.append(output_coef)

        n_iteration = len(output[0])

        fig, ax = plt.subplots()
        for coef in np.arange(n_coefs):
            ax.scatter(np.arange(25*n_iteration, step=25),
                       output[coef], label=f'{library_values[otu_tmp][coef]}', s=1)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Coefficient')
        ax.set_title(f'Coefficients for x{otu_tmp+1}')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(
            f'{folderpath_out}/estimated_coeffs_x{otu_tmp+1}_{iter}.png',
            bbox_inches='tight')
        plt.close()

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

    # run deepmod several times - start loop
    for iter in np.arange(n_runs):
        logging.info(f"iteration {iter} started")

        output_iter = f'{folderpath_out}/output_iteration_{iter}'
        if not os.path.exists(output_iter):
            os.makedirs(output_iter)

        model = run_deepmod_and_save_results(
            iter, dataset, output_iter, network_shape=[hl_size, hl_number])

        logging.info(f"iteration {iter} finished")
