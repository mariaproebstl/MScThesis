# script to run deepmod for ODEs

# General imports
import numpy as np
import torch
import matplotlib.pylab as plt
import pandas as pd
import os
import logging
import tensorflow as tf
from tensorflow.core.util import event_pb2

# DeepMoD functions
from deepymod import DeepMoD
from deepymod.data import Dataset, get_train_test_loader
from deepymod.model.func_approx import NN
from deepymod.model.constraint import LeastSquares
from deepymod.model.sparse_estimators import Threshold, PDEFIND
from deepymod.training import train
from deepymod.training.sparsity_scheduler import TrainTestPeriodic
from deepymod.model.libraryODE import LibraryODE

# # Settings for reproducibility
# np.random.seed(30)
# torch.manual_seed(0)

# Configuring GPU or CPU
if False:  # torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"


################################### Variables ######################################

# name of the dataset
data_name = "miaSim_S4"
# data_name = "ODE_bucci"

# folderpath output
folderpath_out = f"output_{data_name}"
# create output folder
if not os.path.exists(folderpath_out):
    os.makedirs(folderpath_out)

# initialize log file
logging.basicConfig(filename = f'{folderpath_out}/log_{data_name}.log',
                    filemode='w', 
                    format='%(asctime)s, %(message)s', 
                    level=logging.INFO)

# path of data file (input)
# filename = "ts_bucci_subject_5_rel_count.csv"
# filepath = "C:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/explore/data/01b-timeseries-CLVpaper/" + filename
filename = "miaSim_GLV_4species_oscillating_zero.csv"
filepath = "C:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/explore/data/01e-timeseries-miaSim/" + filename

# # specify the number of otu (n_otu) and
# n_otu = 10
# # number of timepoints (n_samples) in the given dataset:
# n_samples = 26
# specify, if not the whole dataset should be used, but only the time points up to max_samples
max_samples = 100

# order of interactions included in the model (2 or 3)
int_order = 2

# specify the network
# number of hidden layers
hl_number = 10
# size of hidden layers
hl_size = 60

logging.info("values are initialized")

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

# function to import the datafile and put it into the right format
def create_data():
    data = np.genfromtxt(filepath, delimiter=",")
    usol = data[1:, :]  # removes header
    if max_samples in globals():
        ts = usol[0:max_samples, 0]
        data_y = usol[0:max_samples, 1:]
    else:
        ts = usol[:, 0]
        data_y = usol[:, 1:]
    T = torch.from_numpy(ts.reshape(-1, 1)).float()
    Y = torch.from_numpy(data_y).float()
    # set dimensions of the dataset
    global n_samples, n_otu
    n_samples, n_otu = data_y.shape
    return T, Y


def prepare_data(folderpath_out):
    # load data
    data = create_data()

    # plot the raw data
    fig, ax = plt.subplots()
    for i in np.arange(n_otu):
        ax.plot(data[0], data[1][:, i], label=f"otu {i+1}")
    ax.set_xlabel("Time")
    ax.set_ylabel("Abundance")
    ax.legend()
    plt.savefig(f'{folderpath_out}/plot_dataset.png')
    plt.close()

    # create the dataset
    dataset = Dataset(
        create_data,
        device=device,
    )

    return dataset


################################ Configuring DeepMoD ##########################################

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
    estimator = Threshold(0.1)
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

    log_path = f'{folderpath_out}/train_log_{iter}'
    if not os.path.exists(log_path):
        os.makedirs(log_path)

    train(
        model,
        train_dataloader,
        test_dataloader,
        optimizer,
        sparsity_scheduler,
        log_dir = log_path,
        max_iterations=100000,
        delta=1e-3,
    )

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

    otu = 0

    train_loss = access_TFRecordDataset(f"loss_mse_output_{otu}", log_path)
    test_loss = access_TFRecordDataset(f"remaining_MSE_test_val_{otu}", log_path)
    l1_loss = access_TFRecordDataset(f"loss_l1_output_{otu}", log_path)
    reg_loss = access_TFRecordDataset(f"loss_reg_output_{otu}", log_path)

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
    ax.set_title(f'loss for otu {otu}')
    ax.legend()
    plt.savefig(f'{folderpath_out}/loss_plot_otu{otu}_iter{iter}.png')
    plt.close()

    # plot and save estimated coefs by iteration
    for otu_tmp in np.arange(n_otu):
        output = []

        for coef in np.arange(n_coefs):
            output_coef = access_TFRecordDataset(f"estimator_coeffs_output_{otu_tmp}_coeff_{coef}", log_path)
            output.append(output_coef)

        n_iteration = len(output[0])

        fig, ax = plt.subplots()
        for coef in np.arange(n_coefs):
            ax.scatter(np.arange(25*n_iteration, step=25),
                       output[coef], label=f'{library_values[otu_tmp][coef]}', s=1)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Coefficient')
        ax.set_title(f'Coefficients for x{otu_tmp+1}')
        ax.legend()
        plt.savefig(f'{folderpath_out}/estimated_coeffs_x{otu_tmp+1}_{iter}.png')
        plt.close()


if __name__ == "__main__":
    
    # load dataset
    data = prepare_data(folderpath_out)
    logging.info("preparation of the dataset is done")

    # run deepmod several times - start loop
    for iter in np.arange(5):
        logging.info(f"iteration {iter} started")

        output_iter = f'{folderpath_out}/output_iteration_{iter}'
        if not os.path.exists(output_iter):
            os.makedirs(output_iter)

        model = run_deepmod_and_save_results(iter, data, output_iter, network_shape=[hl_size, hl_number])
        print(f"finished run_deepmod_and_save_results iteration {iter}")

        logging.info(f"iteration {iter} finished")
