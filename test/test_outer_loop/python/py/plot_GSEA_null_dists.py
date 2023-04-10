from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt

def plot_GSEA_null_dists(
    ds_null,
    gene_set_list,
    gene_sets,
    scores_type = 'positive', # 'positive', 'negative' or 'both'
    plot_sizes = True,
    n_cols = 3,
    plot_x_size = 10,
    x_label = 'x',
    y_label = 'y',
    title = None,
    density_plot_color = 'orange',
    filename = None
    ):

    n_plots = len(gene_set_list) + 1
    n_rows = int(np.ceil(n_plots/n_cols))
    x_size = plot_x_size
    y_size = plot_x_size*0.9 + n_cols*(n_rows - n_cols) 

    fig, axs = plt.subplots(n_rows, n_cols, figsize = (x_size, y_size), dpi = 1000, tight_layout = True)
    
    if title is not None:
        st = fig.suptitle(title, fontsize="x-large")

    # Plot individual histograms and density plots
    r = 0
    for gs in gene_set_list:
        subplot_x = int(np.floor(r/n_cols))
        subplot_y = int(r % n_cols)
        #print('x = {} y = {}'.format(subplot_x, subplot_y))
        max_hist = 0
        if scores_type == 'positive':
            ES_scores = ds_null.loc[gs, ds_null.loc[gs,:] > 0].values
        elif scores_type == 'negative':
            ES_scores = ds_null.loc[gs, ds_null.loc[gs,:] < 0].values
        elif scores_type == 'both':
            ES_scores = ds_null.loc[gs, :].values
        else:
            raise('ERROR: Unknown scores_type: {}'.format(scores_type))
            
        pos_n_bins = 20 #np.int(np.sqrt(len(ES_scores)))
        pos_hist, pos_bins = np.histogram(ES_scores.flatten(), bins = pos_n_bins, density = False)
        pos_bin_width = np.abs(pos_bins[0] - pos_bins[1])/2
        pos_bins = np.append(np.insert(pos_bins, 0, pos_bins[0] - pos_bin_width), pos_bins[-1] + pos_bin_width)
        pos_mid_bins = [(x + y)/2 for x, y in zip(pos_bins[:-1], pos_bins[1:])]
        pos_hist = np.append(np.insert(pos_hist, 0, 0), 0)
        pos_hist = pos_hist/pos_hist.max()
        if np.max(pos_hist) > max_hist:
            max_hist = np.max(pos_hist) 

        axs[subplot_x][subplot_y].set_xlabel(x_label, labelpad = 3, fontsize = 10)
        axs[subplot_x][subplot_y].set_ylabel(y_label, labelpad = 3, fontsize = 10)
        if len(gs) < 35:
            title_font = 9
        else:
            title_font = 9 * 35/len(gs)
        axs[subplot_x][subplot_y].set_title('{}'.format(gs), fontsize = title_font)
        axs[subplot_x][subplot_y].bar(pos_mid_bins, pos_hist, width = 2*pos_bin_width)

        if scores_type == 'positive':
            x_d = np.linspace(-0.1, 1.15*ES_scores.max().max(), 100)
        elif scores_type == 'negative':
            x_d = np.linspace(1.15*ES_scores.min().min(), 0.1, 100)
        elif scores_type == 'both':
            x_d = np.linspace(1.15*ES_scores.min().min(), 1.15*ES_scores.max().max(), 100)
        
        scale = (x_d.max() - x_d.min())/20
        density = sum(norm(xi, scale = scale).pdf(x_d) for xi in ES_scores)
        density = density/density.max()
        axs[subplot_x][subplot_y].plot(x_d, density, linewidth = 2, color = density_plot_color) 
        axs[subplot_x][subplot_y].plot(ES_scores, np.full_like(ES_scores, -0.05*max_hist), '|k', markeredgewidth=0.5)
        axs[subplot_x][subplot_y].set_yticklabels(axs[subplot_x][subplot_y].get_yticklabels(), 
                                              rotation = 0, minor = False, fontsize = 7) 
        axs[subplot_x][subplot_y].set_xticklabels(axs[subplot_x][subplot_y].get_xticklabels(), 
                                              rotation = 0, minor = False, fontsize = 7) 

        r += 1
    
        # Plot combined density plots
      
    j = 0
    for gs in gene_set_list:
        subplot_x = int(np.floor(r/3))
        subplot_y = int(r % 3)
        #print('x = {} y = {}'.format(subplot_x, subplot_y))
        max_hist = 0
        if scores_type == 'positive':
            ES_scores = ds_null.loc[gs, ds_null.loc[gs,:] > 0].values
        elif scores_type == 'negative':
            ES_scores = ds_null.loc[gs, ds_null.loc[gs,:] < 0].values
        elif scores_type == 'both':
            ES_scores = ds_null.loc[gs, :].values
        else:
            raise('ERROR: Unknown scores_type: {}'.format(scores_type))
        
        pos_n_bins = 20 #np.int(np.sqrt(len(ES_scores)))
        pos_hist, pos_bins = np.histogram(ES_scores.flatten(), bins = pos_n_bins, density = False)
        pos_bin_width = np.abs(pos_bins[0] - pos_bins[1])/2
        pos_bins = np.append(np.insert(pos_bins, 0, pos_bins[0] - pos_bin_width), pos_bins[-1] + pos_bin_width)
        pos_mid_bins = [(x + y)/2 for x, y in zip(pos_bins[:-1], pos_bins[1:])]
        pos_hist = np.append(np.insert(pos_hist, 0, 0), 0)
        pos_hist = pos_hist/pos_hist.max()
        if np.max(pos_hist) > max_hist:
            max_hist = np.max(pos_hist) 

        axs[subplot_x][subplot_y].set_xlabel(x_label, labelpad = 3, fontsize = 10)
        axs[subplot_x][subplot_y].set_ylabel(y_label, labelpad = 3, fontsize = 10)
        axs[subplot_x][subplot_y].set_title('Combined Null Dists.', fontsize = 9)
        if len(gs) < 35:
            title_font = 9
        else:
            title_font = 9 * 35/len(gs)

        if scores_type == 'positive':
            x_d = np.linspace(-0.1, 1.15*ES_scores.max().max(), 100)
        elif scores_type == 'negative':
            x_d = np.linspace(1.15*ES_scores.min().min(), 0.1, 100)
        elif scores_type == 'both':
            x_d = np.linspace(1.15*ES_scores.min().min(), 1.15*ES_scores.max().max(), 100)

        scale = (x_d.max() - x_d.min())/20
        density = sum(norm(xi, scale = scale).pdf(x_d) for xi in ES_scores)
        density = density/density.max()
        axs[subplot_x][subplot_y].plot(x_d, density, linewidth = 2, color = density_plot_color) 
        axs[subplot_x][subplot_y].plot(ES_scores, np.full_like(ES_scores, -0.05*max_hist), '|k', markeredgewidth=0.5)
        axs[subplot_x][subplot_y].set_yticklabels(axs[subplot_x][subplot_y].get_yticklabels(), 
                                              rotation = 0, minor = False, fontsize = 7) 
        axs[subplot_x][subplot_y].set_xticklabels(axs[subplot_x][subplot_y].get_xticklabels(), 
                                              rotation = 0, minor = False, fontsize = 7) 

        if plot_sizes is True:
            gs_size = gene_sets[gs]['size']

            if scores_type == 'positive':
                percentile = 95
            elif scores_type == 'negative':
                percentile = 5
            elif scores_type == 'both':
                percentile = 97
            
            axs[subplot_x][subplot_y].text(np.percentile(ES_scores, percentile), j/len(gene_set_list), fontsize = 7,
                                   s = gs_size, ha='center')
        j += 1
    
    # remove empty plots 
    
    for i in range(r+1, n_cols*n_rows):
        subplot_x = int(np.floor(i/3))
        subplot_y = int(i % 3)
        axs[subplot_x][subplot_y].set_visible(False)

        fig.tight_layout()

        # shift subplots down:
        st.set_y(0.95)
        fig.subplots_adjust(top=0.90)

    if filename is not None:
        plt.savefig('{}.png'.format(filename), dpi = 600, facecolor="white", edgecolor="white")

    plt.show()
        
    return
