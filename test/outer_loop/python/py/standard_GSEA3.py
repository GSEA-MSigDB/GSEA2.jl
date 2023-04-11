from statsmodels.sandbox.stats.multicomp import multipletests
import os
import numpy as np
import pandas as pd
from scipy.stats import norm
import random
import matplotlib.pyplot as plt
from .match_target_vs_features import match_target_vs_features
from .KL_GSEA_core3 import KL_GSEA_core3
from .information_coefficient import information_coefficient
from .write_dataset import write_dataset
from .read_dataset import read_dataset
from time import process_time
from .plot_GSEA_null_dists import plot_GSEA_null_dists

def standard_GSEA3(
                target,
                expression_dataset,
                gene_sets,
                results_directory,
                dict_to_map_categorical_to_numeric_vals = None,
                GSEA_enrichment_statistic = 'KS',
                GSEA_alpha = 1.0,
                GSEA_gene_list_norm_type = 'zscore',
                GSEA_NES_normalization = 'zscore', # 'mean_scaling' 'mean_scaling_skew' #'zscore_skew' #'mean_scaling' # 'zscore' #
                n_rand_perm = 100,
                perm_type = 'gene', #'phenotype'  # 'phenotype' or 'gene'
                target_type = 'categorical',
                dataset_type = 'continuous',
                gene_selection_metric = information_coefficient,
                plot_gene_scores_heatmap = True,
                plot_enrichment_plots_for_top_gs = True,
                plot_gs_null_distributions = True,
                save_gs_null_distributions = False,
                max_gs_dist_to_plot = 100,
                n_top_gs_to_plot = 5,
                gene_CI_n_bootstraps = 0,
                #gene_p_val_n_permutations = 0,
                random_num_gene_seed = 1729,
                max_gene_set_size = 500,
                min_gene_set_size = 15,
                print_timings = False,
                use_histogram_to_compute_p_vals = True,
                plot_GSEA_scores_vs_pvals = False,
                results_files_prefix = 'GSEA'):

    start = process_time()
    
    np.random.seed(random_num_gene_seed)

    if not os.path.exists(results_directory):
        os.mkdir(results_directory)

    # Redefine gene sets with the intersection of the gene sets and the gene list

    gene_universe = set(expression_dataset.index)
    
    for gs in list(gene_sets.keys()):
        gene_sets[gs]['genes'] = set(gene_sets[gs]['genes']).intersection(gene_universe)  
        gene_sets[gs]['size'] = len(gene_sets[gs]['genes'])
        
    # Filter gene sets

    print('Number of gene sets as input: {}'.format(len(list(gene_sets.keys()))))

    gene_sets_names = []
    for gs in list(gene_sets.keys()):
        size = gene_sets[gs]['size'] 
        if size >= min_gene_set_size and size <= max_gene_set_size:
            gene_sets_names.append(gs)

    gene_sets = {k: gene_sets[k] for k in gene_sets_names}

    print('Number of gene sets after filtering: {}'.format(len(gene_sets_names)))

    if print_timings is True:        
        print('time filtering gene sets: {} secs.'.format(round(process_time() - start, 5)))

    start = process_time()
        
    #gene_scores_rand = pd.DataFrame(0, index = expression_dataset.index, columns = range(n_rand_perm))
    #gene_scores = pd.DataFrame(0, index = expression_dataset.index, columns = ['Observed'] + list(range(n_rand_perm)))

    # Compute observed gene scores

    if plot_gene_scores_heatmap == True:
        filepath_prefix = '{}/{}_gene_selection'.format(results_directory, results_files_prefix)
    else:
        filepath_prefix = None

    gene_scores_obs =  match_target_vs_features(
                                    target = target,
                                    target_type = target_type,
                                    target_ascending = False,
                                    dict_to_map_categorical_to_numeric_vals = dict_to_map_categorical_to_numeric_vals,
                                    features = expression_dataset,
                                    features_type = dataset_type,
                                    normalize_features = True,
                                    metric = gene_selection_metric,
                                    CI_n_bootstraps = gene_CI_n_bootstraps,
                                    p_val_n_permutations = n_rand_perm,
                                    title = 'Gene Scores',
                                    title_font_size = 12,
                                    n_features_plot = 50,
                                    font_scale = 0.65,
                                    row_compression = 1.0,
                                    figsize = 'auto',
                                    cluster_within_category = False,
                                    plot_heatmap = plot_gene_scores_heatmap,
                                    plot_histograms_of_p_vals = False,
                                    plot_dpi = 500,
                                    random_num_gene_seed = random_num_gene_seed,
                                    fill_na_with_zeroes = True,
                                    save_random_perm_data = True,
                                    save_random_perm_data_file = '{}/{}_rand_perm_gene_scores.txt'.format(results_directory, results_files_prefix),
                                    print_timings = False,
                                    filepath_prefix = filepath_prefix)

    write_dataset(gene_scores_obs, '{}/{}_gene_selection_scores.txt'.format(results_directory, results_files_prefix))

    # Compute random permutation gene scores

    if perm_type == 'phenotype':

        gene_scores_rand = read_dataset('{}/{}_rand_perm_gene_scores.txt'.format(results_directory, results_files_prefix))
            
    elif perm_type == 'gene':

        gene_scores_rand = pd.DataFrame(0, index = gene_scores_obs.index, columns = range(n_rand_perm))
                                            
        for perm in range(n_rand_perm):
            np.random.seed(random_num_gene_seed + perm)
            rand_gene_order = list(np.random.permutation(range(expression_dataset.shape[0])))
            gene_scores_obs2 = gene_scores_obs.loc[:, 'IC']
            gene_scores_rand.loc[:, perm] = gene_scores_obs2.iloc[rand_gene_order].values
            
    else:
        raise TypeError('ERROR: Unknown perm_type: {}'.format(perm_type))

    if print_timings is True:        
        print('time computing gene scores: {} secs.'.format(round(process_time() - start, 5)))

    #gene_scores = pd.concat([gene_scores_obs, gene_scores_rand], axis = 1)
        
    # Compute GSEA scores

    start = process_time()

#    gs_scores_obs = pd.Series(0, index = gene_sets_names)
#    gs_scores_rand = pd.DataFrame(0, index = gene_sets_names, columns = range(n_rand_perm))

    gs_scores_obs = KL_GSEA_core3(
                                    gene_score = gene_scores_obs.loc[:, 'IC'],
                                    GSEA_enrichment_statistic = GSEA_enrichment_statistic,
                                    plot_gene_names = False,
                                    alpha = GSEA_alpha,
                                    IT_statistics_weights = 'half', #proportional', #'half',
                                    KL_area_norm = False,
                                    KL_extended_metric = False,
                                    KL_norm_metric = False,
                                    IT_combine_metrics = True,
                                    produce_intermediate_plots = False,
                                    produce_final_plot = False,
                                    annotation_text_font_size = 14,
                                    gene_sets = gene_sets,
                                    max_gene_set_size = max_gene_set_size,
                                    min_gene_set_size = min_gene_set_size,
                                    sample_norm_type = GSEA_gene_list_norm_type)

    gs_scores_obs = gs_scores_obs.iloc[:, 0]

    gs_scores_rand = KL_GSEA_core3(
    #gs_scores = KL_GSEA_core2(
                                    gene_score = gene_scores_rand,
                                    GSEA_enrichment_statistic = GSEA_enrichment_statistic,
                                    plot_gene_names = False,
                                    alpha = GSEA_alpha,
                                    IT_statistics_weights = 'half', #proportional', #'half',
                                    KL_area_norm = False,
                                    KL_extended_metric = False,
                                    KL_norm_metric = False,
                                    IT_combine_metrics = True,
                                    produce_intermediate_plots = False,
                                    produce_final_plot = False,
                                    annotation_text_font_size = 14,
                                    gene_sets = gene_sets,
                                    max_gene_set_size = max_gene_set_size,
                                    min_gene_set_size = min_gene_set_size,                                 
                                    sample_norm_type = GSEA_gene_list_norm_type)

    #gs_scores_obs = gs_scores.iloc[:, 0]
    #gs_scores_rand = gs_scores.iloc[:, 1:]
    
    # print('gs_scores_obs: {}'.format(np.round(gs_scores_obs, 3)))
    # print('gs_scores_rand: {}'.format(np.round(gs_scores_rand, 3)))

    if print_timings is True:        
        print('time computing GSEA scores: {} secs.'.format(round(process_time() - start, 5)))
    
    # Normalize enrichment scores

    start = process_time()

    gs_pos_mean = pd.Series(0, index = gs_scores_rand.index)
    gs_neg_mean = pd.Series(0, index = gs_scores_rand.index)
    gs_pos_std = pd.Series(0, index = gs_scores_rand.index)
    gs_neg_std = pd.Series(0, index = gs_scores_rand.index)
    gs_pos_mean_right = pd.Series(0, index = gs_scores_rand.index)
    gs_neg_mean_left = pd.Series(0, index = gs_scores_rand.index)
    gs_pos_std_right = pd.Series(0, index = gs_scores_rand.index)
    gs_neg_std_left = pd.Series(0, index = gs_scores_rand.index)
    
    for gs in gs_scores_rand.index: 
        gs_pos_mean.loc[gs] = gs_scores_rand.loc[gs, gs_scores_rand.loc[gs, :] >= 0].mean()
        gs_neg_mean.loc[gs] = gs_scores_rand.loc[gs, gs_scores_rand.loc[gs, :] < 0].mean()
        gs_pos_std.loc[gs] = gs_scores_rand.loc[gs, gs_scores_rand.loc[gs, :] >= 0].std()
        gs_neg_std.loc[gs] = gs_scores_rand.loc[gs, gs_scores_rand.loc[gs, :] < 0].std()
        gs_pos_mean_right.loc[gs] = gs_scores_rand.loc[gs, gs_scores_rand.loc[gs, :] >= gs_pos_mean.loc[gs]].mean()
        gs_neg_mean_left.loc[gs] = gs_scores_rand.loc[gs, gs_scores_rand.loc[gs, :] < gs_neg_mean.loc[gs]].mean()    
        gs_pos_std_right.loc[gs] = gs_scores_rand.loc[gs, gs_scores_rand.loc[gs, :] >= gs_pos_mean.loc[gs]].std()
        gs_neg_std_left.loc[gs] = gs_scores_rand.loc[gs, gs_scores_rand.loc[gs, :] < gs_neg_mean.loc[gs]].std()
        
    gs_scores_rand_norm = pd.DataFrame(0, index = gs_scores_rand.index, columns = gs_scores_rand.columns)

    pos_null_dist = []
    neg_null_dist = []

    for gs in gs_scores_rand_norm.index: 
        for col in gs_scores_rand_norm.columns:
            val = gs_scores_rand.loc[gs, col]
            if val >= 0:
                if GSEA_NES_normalization == 'mean_scaling':
                    val_norm = val/gs_pos_mean.loc[gs]
                elif GSEA_NES_normalization == 'mean_scaling_skew':
                    val_norm = val/gs_pos_mean_right.loc[gs]
                elif GSEA_NES_normalization == 'zscore':
                    val_norm = 1 + (val - gs_pos_mean.loc[gs])/(3 * gs_pos_std.loc[gs])         
                elif GSEA_NES_normalization == 'zscore_skew':
                    val_norm = 1 + (val - gs_pos_mean_right.loc[gs])/(3 * gs_pos_std_right.loc[gs])                      
                elif GSEA_NES_normalization is None:
                    val_norm = val
                else:
                    raise TypeError('ERROR - Unknown GSEA_NES_normalization: {}'.format(GSEA_NES_normalization))
                    
                pos_null_dist.append(val_norm)
                gs_scores_rand_norm.loc[gs, col] = val_norm
                
            else:

                if GSEA_NES_normalization == 'mean_scaling':
                    val_norm =val/np.abs(gs_neg_mean.loc[gs])
                elif GSEA_NES_normalization == 'mean_scaling_skew':
                    val_norm = val/np.abs(gs_neg_mean_left.loc[gs])     
                elif GSEA_NES_normalization == 'zscore':
                    val_norm = -1 + (val - gs_neg_mean.loc[gs])/(3 * gs_neg_std.loc[gs])     
                elif GSEA_NES_normalization == 'zscore_skew':
                    val_norm = -1 + (val - gs_neg_mean_left.loc[gs])/(3 * gs_neg_std_left.loc[gs]) 
                elif GSEA_NES_normalization is None:
                    val_norm = val
                else:
                    raise TypeError('ERROR - Unknown GSEA_NES_normalization: {}'.format(GSEA_NES_normalization))
                    
                neg_null_dist.append(val_norm)
                gs_scores_rand_norm.loc[gs, col] = val_norm

    gs_scores_obs_norm = pd.Series(0, index = gs_scores_obs.index)
    
    for gs in gs_scores_obs.index:
        val = gs_scores_obs.loc[gs]
        
        if val >= 0:
            if GSEA_NES_normalization == 'mean_scaling':
                val_norm = val/gs_pos_mean.loc[gs]
            elif GSEA_NES_normalization == 'mean_scaling_skew':
                val_norm = val/gs_pos_mean_right.loc[gs]
            elif GSEA_NES_normalization == 'zscore':
                val_norm = 1 + (val - gs_pos_mean.loc[gs])/(3 * gs_pos_std.loc[gs])         
            elif GSEA_NES_normalization == 'zscore_skew':
                val_norm = 1 + (val - gs_pos_mean_right.loc[gs])/(3 * gs_pos_std_right.loc[gs])                      
            elif GSEA_NES_normalization is None:
                val_norm = val
            else:
                raise TypeError('ERROR - Unknown GSEA_NES_normalization: {}'.format(GSEA_NES_normalization))
                    
            gs_scores_obs_norm.loc[gs] = val_norm

        else:

            if GSEA_NES_normalization == 'mean_scaling':
                val_norm =val/np.abs(gs_neg_mean.loc[gs])
            elif GSEA_NES_normalization == 'mean_scaling_skew':
                val_norm = val/np.abs(gs_neg_mean_left.loc[gs])     
            elif GSEA_NES_normalization == 'zscore':
                val_norm = -1 + (val - gs_neg_mean.loc[gs])/(3 * gs_neg_std.loc[gs])     
            elif GSEA_NES_normalization == 'zscore_skew':
                val_norm = -1 + (val - gs_neg_mean_left.loc[gs])/(3 * gs_neg_std_left.loc[gs]) 
            elif GSEA_NES_normalization is None:
                val_norm = val
            else:
                raise TypeError('ERROR - Unknown GSEA_NES_normalization: {}'.format(GSEA_NES_normalization))

            gs_scores_obs_norm.loc[gs] = val_norm

    gs_scores_obs_norm.sort_values(ascending = False, inplace = True)

    for i in range(gs_scores_obs_norm.size):
        if gs_scores_obs_norm.iloc[i] < 0:
            first_neg_val = i
            break

    pos_null_size = len(pos_null_dist)
    neg_null_size = len(neg_null_dist)
    pos_null_dist = pd.Series(pos_null_dist)
    neg_null_dist = pd.Series(neg_null_dist)

    gs_p_vals = pd.Series(0, index = gs_scores_obs_norm.index)
    gs_p_vals_FDRs_symbol = pd.Series('', index = gs_scores_obs_norm.index)

    pos_n_bins = np.int(np.sqrt(len(pos_null_dist)))
    if pos_n_bins > 25:
        pos_n_bins = 25            
    pos_hist, pos_bins = np.histogram(pos_null_dist.values.flatten(), bins = pos_n_bins, density = False)
    pos_bin_width = np.abs(pos_bins[0] - pos_bins[1])/2
    pos_bins = np.append(np.insert(pos_bins, 0, pos_bins[0] - pos_bin_width), pos_bins[-1] + pos_bin_width)
    global_pos_mid_bins = [(x + y)/2 for x, y in zip(pos_bins[:-1], pos_bins[1:])]
    pos_hist = np.append(np.insert(pos_hist, 0, 0), 0)
    global_pos_hist = pos_hist/pos_hist.sum()

    neg_n_bins = np.int(np.sqrt(len(neg_null_dist)))
    if neg_n_bins > 25:
        neg_n_bins = 25            
    neg_hist, neg_bins = np.histogram(neg_null_dist.values.flatten(), bins = neg_n_bins, density = False)
    neg_bin_width = np.abs(neg_bins[0] - neg_bins[1])/2
    neg_bins = np.append(np.insert(neg_bins, 0, neg_bins[0] - neg_bin_width), neg_bins[-1] + neg_bin_width)
    global_neg_mid_bins = [(x + y)/2 for x, y in zip(neg_bins[:-1], neg_bins[1:])]
    neg_hist = np.append(np.insert(neg_hist, 0, 0), 0)
    global_neg_hist = neg_hist/neg_hist.sum()

    if use_histogram_to_compute_p_vals == False:

            # Compute p-vals directly from null distribution
            
        for i in range(first_neg_val):
            val = gs_scores_obs_norm.iloc[i]
            gs_p_vals.iloc[i] = (pos_null_dist >= val).sum()/pos_null_size
            if gs_p_vals.iloc[i] < 1.0/pos_null_size:
                gs_p_vals.iloc[i] = 1.0/pos_null_size
                gs_p_vals_FDRs_symbol.iloc[i] = '\u2264'

        for i in range(first_neg_val, gs_scores_obs_norm.size):
            val = gs_scores_obs_norm.iloc[i]
            gs_p_vals.iloc[i] = (neg_null_dist <= val).sum()/neg_null_size
            if gs_p_vals.iloc[i] < 1.0/neg_null_size:
                gs_p_vals.iloc[i] = 1.0/neg_null_size
                gs_p_vals_FDRs_symbol.iloc[i] = '\u2264'

    else:
          # Compute p-vals from null distribution histogram
            
        pos_hist_cumsum = np.round(pos_hist.cumsum()[::-1], 10)
        for i in range(first_neg_val):
            pos_val = gs_scores_obs_norm.iloc[i]
            pos_diffs = list(np.abs(pos_val - pos_mid_bins))
            pos_min_index = pos_diffs.index(min(pos_diffs))
            gs_p_vals.iloc[i] = pos_hist_cumsum[pos_min_index]
            if gs_p_vals.iloc[i] < 1.0/pos_null_size:                
                gs_p_vals.iloc[i] = 1.0/pos_null_size
                gs_p_vals_FDRs_symbol.iloc[i] = '\u2264'
            #print('val: {} pos_min_index: {} p-val: {}'.format(pos_val,  pos_min_index, gs_p_vals.iloc[i]))
            
        neg_hist_cumsum = np.round(neg_hist.cumsum(), 10)
        for i in range(first_neg_val, gs_scores_obs_norm.size):       
            neg_val = gs_scores_obs_norm.iloc[i]
            neg_diffs = list(np.abs(neg_val - neg_mid_bins))
            neg_min_index = neg_diffs.index(min(neg_diffs))
            gs_p_vals.iloc[i] = neg_hist_cumsum[neg_min_index]
            if gs_p_vals.iloc[i] < 1.0/neg_null_size:                
                gs_p_vals.iloc[i] = 1.0/neg_null_size
                gs_p_vals_FDRs_symbol.iloc[i] = '\u2264'
            #print('val: {} neg_min_index: {} p-val: {}'.format(neg_val,  neg_min_index, gs_p_vals.iloc[i]))            

    if plot_GSEA_scores_vs_pvals == True:
        plt.xlabel('GSEA Scores')
        plt.ylabel('Log(p-value)')
        plt.title('Plot of p-values vs. GSEA scores')
        #plt.plot(gs_scores_obs_norm, np.log(gs_p_vals))
        plt.plot(gs_scores_obs_norm, np.log(gs_p_vals))
        # add dotted lines at 1.0/neg_null_size 1.0/pos_null_size to show minimum p-vals

    # Compute FDRs
            
    pos_FDRs = multipletests(gs_p_vals.iloc[0:first_neg_val], method = 'fdr_bh')[1]
    neg_FDRs = multipletests(gs_p_vals.iloc[first_neg_val:gs_scores_obs_norm.size], method = 'fdr_bh')[1]
    gs_FDRs = pd.Series(list(pos_FDRs) + list(neg_FDRs), index = gs_p_vals.index)
        
    # Round p-vals and FDRs
    for gs in gs_p_vals.index:
        if gs_p_vals.loc[gs] < 0.001:
            gs_p_vals.loc[gs] = '{:.3E}'.format(gs_p_vals.loc[gs])
            gs_FDRs.loc[gs] = '{:.3E}'.format(gs_FDRs.loc[gs])
        else:
            gs_p_vals.loc[gs] = '{:.3f}'.format(gs_p_vals.loc[gs])
            gs_FDRs.loc[gs] = '{:.3f}'.format(gs_FDRs.loc[gs])           

    gs_scores_obs_norm = gs_scores_obs_norm.round(decimals = 3)
    gs_scores_obs = gs_scores_obs.round(decimals = 3)

    # Make results table
      
    top_g_sets = [x for x in list(gs_scores_obs_norm.index[range(n_top_gs_to_plot)]) + 
                                 list(gs_scores_obs_norm.index[range(gs_scores_obs_norm.size - 
                                        n_top_gs_to_plot, gs_scores_obs_norm.size)])]

    gs_p_vals_labels = pd.Series(['{}{}'.format(x,y) for x,y in zip(gs_p_vals_FDRs_symbol, gs_p_vals)], index = gs_p_vals.index)
    gs_FDRs_labels = pd.Series(['{}{}'.format(x,y) for x,y in zip(gs_p_vals_FDRs_symbol, gs_FDRs)], index = gs_FDRs.index)
    rank = pd.Series(range(1, gs_scores_obs_norm.size + 1), index = gs_scores_obs_norm.index)
        
    res_table = pd.concat([rank,
                           pd.Series([gene_sets[x]['size'] for x in gene_sets_names], index = gene_sets_names),
                           gs_scores_obs_norm, 
                           gs_scores_obs, 
                           gs_p_vals_labels,
                           gs_FDRs_labels], axis = 1)

    res_table.columns = ['Rank', 'Size', 'NES', 'ES', 'p-val', 'FDR']

    # Display top results
    
    top_res_table = res_table.loc[top_g_sets, :]
    
    #display(top_res_table)   
    
    # Save results table
    
    write_dataset(res_table, '{}/{}_GSEA_results_table.txt'.format(results_directory, results_files_prefix))

    if print_timings is True:        
        print('time computing p-vals and FDRs: {} secs.'.format(round(process_time() - start, 5)))
    
    # Plot enrichment plots for top gene sets
        
    if plot_enrichment_plots_for_top_gs == True:
        
        for g_set in top_g_sets:
            
            KL_GSEA_core3(
                gene_score = gene_scores_obs.loc[:, 'IC'],
                GSEA_enrichment_statistic = GSEA_enrichment_statistic,
                title = '{}   ({}, statistic: {})'.format(g_set, results_files_prefix, GSEA_enrichment_statistic),
                plot_gene_names = False,
                alpha = 1,
                IT_statistics_weights = 'half', # proportional', #'half',
                KL_area_norm = False,
                KL_extended_metric = False,
                KL_norm_metric = False,
                IT_combine_metrics = True,
                produce_intermediate_plots = False,
                produce_final_plot = True,
                annotation_text_font_size = 14,
                gene_sets = gene_sets[g_set]['genes'],
                enrichment_plot_subtitle = 'NES: {}     ES: {}     p-val: {}     FDR: {}'.format(
                                                                                 gs_scores_obs_norm.loc[g_set],
                                                                                 gs_scores_obs.loc[g_set],
                                                                                 gs_p_vals.loc[g_set],
                                                                                 gs_FDRs.loc[g_set]),
                html_file_path = '{}/{}_top_gene_set_{}.html'.format(results_directory, results_files_prefix, g_set),
                image_filename = '{}/{}_top_gene_set_{}.png'.format(results_directory, results_files_prefix, g_set))

    # Plot per gene set and global null distributions          
                
    if plot_gs_null_distributions == True:
 
        if gs_scores_rand.shape[0] < max_gs_dist_to_plot:
            g_sets = gs_scores_rand.index
        else:
            g_sets = random.choices(gs_scores_rand.index, k = max_gs_dist_to_plot)

        # Unnormalized per gene set distributions
        
        plot_GSEA_null_dists(
            ds_null = gs_scores_rand,
            gene_set_list = g_sets,
            gene_sets = gene_sets,
            n_cols = 3,
            plot_x_size = 10,
            scores_type = 'positive',
            x_label = 'Unnormalized GSEA Scores (ES)',
            y_label = 'Density',
            title = 'Unnormalized GSEA Null Distributions (Positive)',
            density_plot_color = 'orange',
            filename = ('{}/{}_unnorm_GSEA_pos_null_dists.png'.format(results_directory, results_files_prefix)))

        plot_GSEA_null_dists(
            ds_null = gs_scores_rand,
            gene_set_list = g_sets,
            gene_sets = gene_sets,
            n_cols = 3,
            plot_x_size = 10,
            scores_type = 'negative',
            x_label = 'Unnormalized GSEA Scores (ES)',
            y_label = 'Density',
            title = 'Unnormalized GSEA Null Distributions (Negative)',
            density_plot_color = 'orange',
            filename = ('{}/{}_unnorm_GSEA_neg_null_dists.png'.format(results_directory, results_files_prefix)))

        plot_GSEA_null_dists(
            ds_null = gs_scores_rand,
            gene_set_list = g_sets,
            gene_sets = gene_sets,
            n_cols = 3,
            plot_x_size = 10,
            scores_type = 'both',
            x_label = 'Unnormalized GSEA Scores (ES)',
            y_label = 'Density',
            title = 'Unnormalized GSEA Null Distributions (Positive and Negative)',
            density_plot_color = 'orange',
            filename = ('{}/{}_unnorm_GSEA_null_dists.png'.format(results_directory, results_files_prefix)))

        # Normalized per gene set distributions

        plot_GSEA_null_dists(
            ds_null = gs_scores_rand_norm,
            gene_set_list = g_sets,
            gene_sets = gene_sets,
            n_cols = 3,
            plot_x_size = 10,
            scores_type = 'positive',
            x_label = 'Normalized GSEA Scores (NES)',
            y_label = 'Density',
            title = 'Normalized GSEA Null Distributions (Positive)',
            density_plot_color = 'orange',
            filename = ('{}/{}_norm_GSEA_pos_null_dists.png'.format(results_directory, results_files_prefix)))

        plot_GSEA_null_dists(
            ds_null = gs_scores_rand_norm,
            gene_set_list = g_sets,
            gene_sets = gene_sets,
            n_cols = 3,
            plot_x_size = 10,
            scores_type = 'negative',
            x_label = 'Normalized GSEA Scores (NES)',
            y_label = 'Density',
            title = 'Normalized GSEA Null Distributions (Negative)',
            density_plot_color = 'orange',
            filename = ('{}/{}_norm_GSEA_neg_null_dists.png'.format(results_directory, results_files_prefix)))

        plot_GSEA_null_dists(
            ds_null = gs_scores_rand_norm,
            gene_set_list = g_sets,
            gene_sets = gene_sets,
            n_cols = 3,
            plot_x_size = 10,
            scores_type = 'both',
            x_label = 'Normalized GSEA Scores (NES)',
            y_label = 'Density',
            title = 'Normalized GSEA Null Distributions (Positive and Negative)',
            density_plot_color = 'orange',
            filename = ('{}/{}_norm_GSEA_null_dists.png'.format(results_directory, results_files_prefix)))

        # Global distributions
        
        fig, axs = plt.subplots(1, 2, figsize = (10, 5), dpi = 500, tight_layout = True)

        st = fig.suptitle('Global null and observed positive and negative distributions of NES scores', fontsize="x-large")

        # Positive densities
        
        x_d = np.linspace(-0.1, 1.15*np.max(pos_null_dist), 100)
        #scale = (x_d.max() - x_d.min())/20
        density = sum(norm(xi, scale = 0.1).pdf(x_d) for xi in pos_null_dist)
        density = density/np.max(density)
        gs_scores_obs_norm_pos = gs_scores_obs_norm.iloc[0: first_neg_val]
        density_obs = sum(norm(xi, scale = 0.1).pdf(x_d) for xi in gs_scores_obs_norm_pos)
        density_obs = density_obs/np.max(density_obs)
        axs[0].set_xlabel('Normalized Enrichment Score (NES)', labelpad = 3, fontsize = 13)
        axs[0].set_ylabel('Probability Density', labelpad = 3, fontsize = 13)
        axs[0].plot(x_d, density, linewidth = 3, color = 'green')
        axs[0].plot(x_d, density_obs, linewidth = 3, color = 'blue')
        axs[0].plot(pos_null_dist, np.full_like(pos_null_dist, -0.05), '|k', markeredgewidth=0.5)
        #axs[0].axis([-0.1, 1.15*np.max(pos_null_dist), -0.1, 1.2*np.max(density)])
        axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation = 0, minor = False, fontsize = 12) 
        axs[0].set_yticklabels(axs[0].get_yticklabels(), rotation = 0, minor = False, fontsize = 12) 
       
        #axs[0].title('{} Global positive null and observed distribution of NES scores'.format(results_files_prefix))

        # Negative densities
        
        x_d = np.linspace(1.15*np.min(neg_null_dist), 0.1, 100)
        #scale = (x_d.max() - x_d.min())/20
        density = sum(norm(xi, scale = 0.1).pdf(x_d) for xi in neg_null_dist)
        density = density/np.max(density)
        gs_scores_obs_norm_neg = gs_scores_obs_norm.iloc[first_neg_val: gs_scores_obs_norm.size]
        density_obs = sum(norm(xi, scale = 0.1).pdf(x_d) for xi in gs_scores_obs_norm_neg)
        density_obs = density_obs/np.max(density_obs)
        axs[1].set_xlabel('Normalized Enrichment Score (NES)', labelpad = 3, fontsize = 13)
        axs[1].set_ylabel('Probability Density', labelpad = 3, fontsize = 13)
        axs[1].plot(x_d, density, linewidth = 3, color = 'green')
        axs[1].plot(x_d, density_obs, linewidth = 3, color = 'blue')        
        axs[1].plot(neg_null_dist, np.full_like(neg_null_dist, -0.05), '|k', markeredgewidth=0.5)
        #axs[1].axis([1.15*np.min(neg_null_dist), 0.1, -0.1, 1.15*np.max(density)])
        axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation = 0, minor = False, fontsize = 12) 
        axs[1].set_yticklabels(axs[1].get_yticklabels(), rotation = 0, minor = False, fontsize = 12) 
        
        #axs[1].title('{} Global negative null and observed distribution of NES scores'.format(results_files_prefix))

        fig.savefig('{}/{}_global_null_dist.png'.format(results_directory, results_files_prefix), dpi = 600, 
            facecolor="white",edgecolor="white")

        fig.show()
        
        # Save null distributions

        if save_gs_null_distributions is True:

            write_dataset(gs_scores_rand, '{}/{}_unnorm_gene_set_null_dists.txt'.format(results_directory, results_files_prefix))
            write_dataset(gs_scores_rand_norm, '{}/{}_norm_gene_set_null_dists.txt'.format(results_directory, results_files_prefix))
        
    return(res_table)
