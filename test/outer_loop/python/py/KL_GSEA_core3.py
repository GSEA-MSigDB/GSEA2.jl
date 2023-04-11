from warnings import warn
import pandas as pd
from numpy import absolute, in1d, asarray, round, log, exp, zeros
from .make_enrichment_plot import make_enrichment_plot
from time import process_time

import pyximport
import numpy as np
pyximport.install(setup_args={"include_dirs":np.get_include(), "extra_compile_args": ['-O3', '-march=native']}, inplace=True)
from .compute_KL_PQ_enrichment import compute_KL_PQ_enrichment     
from .compute_KL_PRURm_enrichment import compute_KL_PRURm_enrichment  
from .compute_KL_PRURp_enrichment import compute_KL_PRURp_enrichment  
from .compute_KL_PR_enrichment import compute_KL_PR_enrichment
from .compute_KL_PU_enrichment import compute_KL_PU_enrichment     
from .compute_KS_AUC_enrichment import compute_KS_AUC_enrichment    
from .compute_KS_SUP_enrichment import compute_KS_SUP_enrichment     

import sys

def KL_GSEA_core3(
    gene_score, # a series or dataframe of gene scores
    gene_sets, # dictionary containing one of more gene sets
    filter_gene_sets = False,
    min_gene_set_size = 15,
    max_gene_set_size = 500,
    gene_set_weights = None,
    GSEA_enrichment_statistic = 'KS_SUP',
        # Standard (KS) GSEA:
        #     KS_SUP: Kolmogorov-Smirnov Supremum 
        #     KS_AUC: Kolmogorov-Smirnov Area Under the Curve
        # Information Theoretical (KL) GSEA:
        #     KL_PR: Kullback-Leibler Divergence of P and R (KL-Li)
        #     KL_PQ: Kullback-Leibler Divergence of P and Q
        #     KL_PRURm: Kullback-Leibler Divergence of P and R minus Divergence of U and R (KL-Lio-)
        #     KL_PRURp: Kullback-Leibler Divergence of P and R plus Divergence of U and R (KL-Lio+)
    IT_statistics_weights = 'half', # half and half as in the original JS metric
                          # 'proportional', # proportional so that R is the mixture of P and U (R = pi_P * P + pi_U * U)
    KL_area_norm = False,
    KL_extended_metric = True,
    KL_norm_metric = False,
    IT_combine_metrics = False,
    alpha = 1.0,
    produce_intermediate_plots = False,
    produce_final_plot = False,
    plot_gene_names = False,
    title = None,
    gene_score_name = None,
    annotation_text_font_size = 20,
    annotation_text_width = 100,
    annotation_text_yshift = 50,
    html_file_path = None,
    sample_norm_type = None, # 'zscore', 'rank' or None
    sample_norm_clims = [-3, 3],
    plotly_html_file_path = None,
    print_timings = False,
    enrichment_plot_subtitle = None,
    image_filename = None):

    # Normalize gene scores

    start = process_time()

    # if gene_score is a series convert it into a dataframe
    if isinstance(gene_score, pd.Series):
        gene_score = gene_score.to_frame()

    if sample_norm_type == 'rank':
        gene_score = gene_score.rank(axis = 0, method = 'dense', numeric_only = None, na_option = 'keep',
                                         ascending = True, pct = False)
        gene_score = (gene_score - gene_score.min())/(gene_score.max() - gene_score.min())
    elif sample_norm_type == 'zscore':
        D = gene_score.T            
        D1 = ((D - D.mean(axis=0))/D.std(axis=0)).clip(lower = sample_norm_clims[0], upper = sample_norm_clims[1])
        gene_score = D1.T
    elif sample_norm_type is not None:
        sys.exit('ERROR: unknown sample_norm_type: {}'.format(sample_norm_type))

    if print_timings is True:
        print('time normalizing gene scores: {} secs.'.format(round(process_time() - start, 5)))
        
    # Define gene set indicator matrix and gene_score sorting index matrix

    start = process_time()
    gene_list = list(gene_score.index)
    n_rows = gene_score.shape[0]
    n_cols = gene_score.shape[1]

    print('Number of gene sets as input: {}'.format(len(list(gene_sets.keys()))))
    
    for gs in list(gene_sets.keys()):
        removed_count = 0
        genes = gene_sets[gs]['genes'].intersection(gene_list)
        if len(genes) == 0:
            gene_sets.pop(gs)
            removed_count += 1
        else:
            gene_sets[gs]['genes'] = genes

    print('Number of gene sets removed for having zero overlap with gene list: {}'.format(removed_count))

    # Filter gene sets

    if filter_gene_sets is True:
    
        gene_sets_names = []
        for gs in list(gene_sets.keys()):
            size = gene_sets[gs]['size'] 
            if size >= min_gene_set_size and size <= max_gene_set_size:
                gene_sets_names.append(gs)

        gene_sets = {k: gene_sets[k] for k in gene_sets_names}

        print('Number of gene sets after filtering: {}'.format(len(gene_sets)))

    I_h_matrix = pd.DataFrame(0, index = gene_score.index, columns = list(gene_sets.keys()))
    
    for gs in list(gene_sets.keys()):
        I_h_matrix.loc[gene_sets[gs]['genes'], gs] = 1

    I_h_matrix = I_h_matrix.values
    #I_h_matrix = np.ascontiguousarray(I_h_matrix.values)
    
#    gene_score_sorting_matrix = pd.DataFrame(0, index = gene_score.index, columns = gene_score.columns)

#    for p in gene_score.columns:
#        gene_score_instance = gene_score.loc[:, p].values
#        # kind{‘quicksort’, ‘mergesort’, ‘heapsort’, ‘stable’},
#        gene_score_sorting_matrix.loc[:, p] = gene_score_instance.argsort()
#    gene_score_sorting_matrix = gene_score.values.argsort(axis = 1)

#    gene_score_sorting_matrix = pd.DataFrame(gene_score.values.argsort(axis = 0)[::-1], index = gene_score.index, columns = gene_score.columns)
    gene_score_columns = gene_score.columns
    gene_score_rows =  gene_score.index
    gene_score = gene_score.values
#    gene_score_sorting_matrix = gene_score.values.argsort(axis = 0)[::-1]
    gene_score_sorting_matrix = gene_score.argsort(axis = 0)[::-1]

    #sorted_gene_score_matrix = np.empty([n_rows, n_cols], order='C')
    sorted_gene_score_matrix = np.empty([n_rows, n_cols])
    for p in range(n_cols):
        sorted_gene_score_matrix[:, p] = gene_score[gene_score_sorting_matrix[:, p], p]

    #sorted_gene_score_matrix = np.ascontiguousarray(sorted_gene_score_matrix)
        
    if print_timings is True:        
        print('time computing gene set indicator and gene_score sorting matrices: {} secs.'.format(round(process_time() - start, 5)))
        
    # loop to compute GSEA Enrichment Scores (ES) over gene_score and gene sets instances

    start = process_time()
        
    ES_matrix = pd.DataFrame(0, index = list(gene_sets.keys()), columns = gene_score_columns)

    #####
    #print('sorted_gene_score_matrix: {}'.format(sorted_gene_score_matrix.flags))
    #print('I_h_matrix: {}'.format(I_h_matrix.flags))    
    #####
    #for gs in list(gene_sets.keys()):

    for gs in range(len(gene_sets)):
        #print('computing GSEA scores for gene set: {} ({}/{})'.format(list(gene_sets.keys())[gs], gs, len(gene_sets)))
        #for p in gene_score.columns:
        for p in range(n_cols):
            #gene_score_instance = gene_score.loc[:, p].values
            #sorted_gene_score = gene_score_instance[gene_score_sorting_matrix.loc[:, p].values]
            #I_h_instance = I_h_matrix.loc[:, gs].values
            #sorted_I_h = I_h_instance[gene_score_sorting_matrix.loc[:, p].values]

            #gene_score_instance = gene_score.iloc[:, p].values
            #sorted_gene_score = gene_score_instance[gene_score_sorting_matrix[:, p]]
            ###sorted_gene_score = gene_score[gene_score_sorting_matrix[:, p], p]
            #I_h_instance = I_h_matrix.loc[:, gs].values
            #sorted_I_h = I_h_instance[gene_score_sorting_matrix[:, p]]
            ###sorted_I_h = I_h_matrix[gene_score_sorting_matrix[:, p], gs]

            # Compute KS
        
            ###Delta, ES = compute_KS_or_AUC_enrichment(sorted_gene_score, sorted_I_h, alpha, 1)
            #Delta, ES = compute_KS_or_AUC_enrichment(#gene_score[gene_score_sorting_matrix[:, p], p],
            #                                          sorted_gene_score_matrix[:, p],
            #                                         I_h_matrix[gene_score_sorting_matrix[:, p], gs],
            #                                         alpha,
            #                                         1)
            #print(gene_score_columns[p])
            A = np.ascontiguousarray(sorted_gene_score_matrix[:, p])
            B = np.ascontiguousarray(I_h_matrix[gene_score_sorting_matrix[:, p], gs])

            # Call cython functions implementing each GSEA enrichment statistic    
            
            if GSEA_enrichment_statistic == 'KS_SUP': # Kolmogorov-Smirnov Supremum 
                Delta, ES = compute_KS_SUP_enrichment(A, B, alpha)
            elif GSEA_enrichment_statistic == 'KS_AUC': # Kolmogorov-Smirnov Area Under the Curve
                Delta, ES = compute_KS_AUC_enrichment(A, B, alpha)                
            elif GSEA_enrichment_statistic == 'KL_PR': # Kullback-Leibler Divergence of P and R
                Delta, ES = compute_KL_PR_enrichment(A, B, alpha)
            elif GSEA_enrichment_statistic == 'KL_PQ': # Kullback-Leibler Divergence of P and Q
                Delta, ES = compute_KL_PQ_enrichment(A, B, alpha)
            elif GSEA_enrichment_statistic == 'KL_PU': # Kullback-Leibler Divergence of P and U
                Delta, ES = compute_KL_PU_enrichment(A, B, alpha)
            elif GSEA_enrichment_statistic == 'KL_PRURm': # Kullback-Leibler Divergence of P and R minus Divergence of U and R 
                Delta, ES = compute_KL_PRURm_enrichment(A, B, alpha)
            elif GSEA_enrichment_statistic == 'KL_PRURp': # Kullback-Leibler Divergence of P and R plus Divergence of U and R
                Delta, ES = compute_KL_PRURp_enrichment(A, B, alpha)
            else:
                raise('ERROR: unknown GSEA_enrichment_statistic: {}'.format(GSEA_enrichment_statistic))

            ES_matrix.iloc[gs, p] = ES

            
            if produce_final_plot is True:

                if enrichment_plot_subtitle is None:
                    enrichment_plot_subtitle = 'ES: {}'.format(np.round(ES, 3))


                gene_score_series = pd.Series(A, index = gene_score_rows[gene_score_sorting_matrix[:, p]])
                make_enrichment_plot(Delta, B, gene_score_series, ES, None, None, title, gene_score_name,
                                     annotation_text_font_size, annotation_text_width,
                                     annotation_text_yshift, plot_gene_names, html_file_path, plotly_html_file_path,
                                     subtitle = enrichment_plot_subtitle, image_filename = image_filename)
        
    if print_timings is True:        
        print('time computing GSEA Enrichment Scores (ES) over gene_score and gene sets instances: {} secs.'.format(round(process_time() - start, 5)))

    ES_matrix = ES_matrix.round(decimals = 3)
    
    return ES_matrix
