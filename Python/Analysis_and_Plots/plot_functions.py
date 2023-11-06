# imports
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec


#### Heatmaps

# function to plot a heatmap of A with rectangles around the largest 20% abs(values)
def plot_heatmap(matrix_A, ax, fig, title = "", colnames = None, Mat = "", rec_per = 20, add_caption = True, add_outline = False):

    n_taxa = matrix_A.shape[0]

    fig.set_figwidth(n_taxa)
    fig.set_figheight(0.75*n_taxa)

    # make heatmap
    if Mat == "weightsMat":
        cmap = LinearSegmentedColormap.from_list("white_to_green", ["white", "darkgreen"])
        sns.heatmap(matrix_A, annot=True, fmt="", cmap=cmap, vmin=0, ax=ax)
    else:
        sns.heatmap(matrix_A, annot=True, fmt="", cmap="RdBu", center=0, ax=ax)
    
    # Setting x and y axis labels to the column names of the data
    if colnames is None:
        colnames = [f"x{i}" for i in range(1, n_taxa + 1)]
    
    tick_positions = np.arange(0.5, n_taxa, 1)  # Position at the center of each cell
    ax.set_xticks(tick_positions)
    ax.set_yticks(tick_positions)
    
    ax.set_xticklabels(colnames, rotation=90) if max(len(name) for name in colnames) >= 6 else ax.set_xticklabels(colnames)
    ax.set_yticklabels(colnames, rotation=0)
        
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.tick_params(left=False, top=False)
    
    # Ass border around the largest 20% of the values in A
    threshold = np.percentile(abs(matrix_A), (100-rec_per)) # threshold for the top 20% of values
    positions = np.argwhere(abs(matrix_A) >= threshold) # positions of the top 20% of values
    # Add rectangles around the top 20% of values
    for pos in positions:
        rect = Rectangle((pos[1], pos[0]), 1, 1, linewidth=1, edgecolor='black', facecolor='none')
        ax.add_patch(rect)
        
    # Add a border around the entire plot
    if add_outline:
        rect = Rectangle((0, 0), 1, 1, linewidth=2, edgecolor='black', facecolor='none', transform=ax.transAxes)
        ax.add_patch(rect)

    # add title
    ax.set_title(title)
    # add caption
    if add_caption:
        fig.text(0.8, 0.05, f'threshold = {round(threshold,2)} (top {rec_per}%)', ha='right', va='bottom', fontsize=12, transform=fig.transFigure)

    
## add function to plot vector g and interaction matrix A as heatmaps in one plot:
def plot_heatmaps(A, g, n_taxa, title = "", colnames = None):
    
    if g is None:
        max_value = max(A.max(), abs(A.min()))
    else:
        max_value = max(A.max(), g.max(), abs(A.min()), abs(g.min()))
        g = round(pd.DataFrame(g), 3)
    
    A = round(pd.DataFrame(A), 3)

    # Creating figure and gridspec
    fig = plt.figure(figsize=(1 * n_taxa, 0.8 * n_taxa))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, n_taxa], wspace=0.1)  

    # Creating subplots
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    # Plotting the vector g as a heatmap
    if g is None:
        sns.heatmap(np.zeros((n_taxa, 1)), cmap="RdBu", center=0, vmin=-max_value, vmax=max_value, 
                    cbar=False, ax=ax0)  #  linewidths=0.5, linecolor='black'
        # ax0.set_title('growth vector g')
        ax0.xaxis.tick_top()
        ax0.set_yticklabels([f'x{i+1}' for i in range(n_taxa)], rotation=0)
        ax0.tick_params(left=False, top=False)
    else:
        sns.heatmap(g, cmap="RdBu", vmin=-max_value, vmax=max_value, annot=True, 
                    cbar=False, ax=ax0)  #  linewidths=0.5, linecolor='black'
        # ax0.set_title('growth vector g')
        ax0.xaxis.tick_top()
        ax0.set_yticklabels([f'x{i+1}' for i in range(g.shape[0])], rotation=0)
        ax0.tick_params(left=False, top=False)

    # Plotting the matrix A as a heatmap
    sns.heatmap(A, cmap="RdBu", vmin=-max_value, vmax=max_value, annot=True, yticklabels=False, 
                cbar=False, ax=ax1)  #  linewidths=0.5, linecolor='black'
    # ax1.set_title('interaction matrix A')
    ax1.xaxis.tick_top()
    ax1.tick_params(left=False, top=False)
    if colnames is not None:
        ax0.set_xticklabels(["1 (growth rates)"], rotation=90, fontsize = 10)
        ax1.set_xticklabels([f'x{i+1} ({colnames[i]})' for i in range(A.shape[1])], rotation=90, fontsize = 10) if max(len(name) for name in colnames) >= 6 else ax1.set_xticklabels(colnames, fontsize = 10)
    else:
        ax0.set_xticklabels(['1'], fontsize = 10)
        ax1.set_xticklabels([f'x{i+1}' for i in range(A.shape[1])], fontsize = 10)


    # Adding a colorbar for the entire figure
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(ax1.collections[0], cax=cbar_ax)

    # outline the largest 20% of the values in A
    threshold = np.percentile(abs(A), 80) # threshold for the top 20% of values
    positions = np.argwhere(abs(A) >= threshold) # positions of the top 20% of values
    # Add rectangles around the top 20% of values
    for pos in positions:
        rect = Rectangle((pos[1], pos[0]), 1, 1, linewidth=1, edgecolor='black', facecolor='none')
        ax1.add_patch(rect)
        
    # Add a border around the entire plot
    for axis in [ax0, ax1]:
        rect = Rectangle((0, 0), 1, 1, linewidth=2, edgecolor='black', facecolor='none', transform=axis.transAxes)
        axis.add_patch(rect)
    
    title_text = ax1.set_title(title, fontsize = 16)
    title_text.set_position((0.4, 1.1))


#### Prediction Plots

# plot fits or prediction vs. observed data for one taxon / one plot
def plot_pred(T_pred, data_pred, T_org, data_org, ax, title = "", legend = True):
    # load colorpalette for later plots
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # plot each taxon timeline separately
    ax.plot(T_org, data_org, linewidth = 0.6, label = "data", color = colors[1]) # light color
    ax.plot(T_pred, data_pred, label = "pred", color = colors[0]) # dark color
    # xlabel
    ax.set_xlabel("time")
    ax.set_title(title, pad = 10)
    if legend:
        ax.legend()

# plot fits or prediction vs. observed data for all taxa in one line
def plot_pred_line(T_pred, data_pred, T_org, data_org, n_taxa, title = ""):
    # load colorpalette for later plots
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    n_col = 1
    n_row = n_taxa

    fig, axs = plt.subplots(n_row, n_col)
    fig.suptitle(title)
    fig.set_figwidth(4)
    fig.set_figheight(10)
    fig.tight_layout(pad=1.2)

    for i in np.arange(n_taxa):
        # plot each taxon timeline separately
        # axs[i].plot(T_org, data_org[:,i], linewidth = 0.7, label = "data", color = colors[1])
        # axs[i].plot(T_pred, data_pred[:,i], linewidth = 2, label = "fit")
        plot_pred(T_org, data_org[:,i], T_pred, data_pred[:,i], ax = axs[i])
        # y label (x1, x2, ...)
        axs[i].set_ylabel(f"x{i+1}", rotation = 0, fontsize = 14)
        axs[i].yaxis.set_label_coords(-0.2, 0.5)
        # xlabel
    axs[i].set_xlabel("time")
                        
    # plt.setp(axs, ylim=(min(pred[1:,].min(), P[0].min()) - 0.01, max(pred[1:,].max(), P[0].max()) + 0.05))
    # add overall legend below the plots
    handles, labels = [], []
    for ax in axs.ravel():
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h)
                labels.append(l)
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol = 2)

