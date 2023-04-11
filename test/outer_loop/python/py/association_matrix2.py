from numpy import asarray
from pandas import DataFrame

from .information_coefficient import information_coefficient
from .plot_heat_map import plot_heat_map

import pyximport
import numpy as np
pyximport.install(setup_args={"include_dirs":np.get_include()}, inplace=True)
from .signal_to_noise_GSEA4 import signal_to_noise_GSEA4

def association_matrix2(df_A, df_B, assoc_function = information_coefficient, axis = 0, title = None, name_A = None, name_B = None, make_plot = False):

    """ 
    Computes an association matrix between two dataframes using a given association function. 

        Parameters: 
            x (series): first vector.
            y (series): second vector.
            n_grid (int): the size of the grid where the kernel probability 
                 densities are computed.
            df_A (DataFrame): first input dataframe.
            df_B (DataFrame): second input dataframe.
            assoc_function (function): function to compute the association.
            axis (0 or 1): axis to use as internal dimension to compute association across.
            title (str): title of the plot.
            name_A (str): name of the first dataframe for the plot.
            name_B (str): name of the second dataframe for the plot.
            make_plot (logical): if the plot of the association matrix will be produced.

        Returns: 
            assoc_matrix (dataframe): association matrix.
        
    """
    
    if axis == 1:

        if assoc_function == signal_to_noise_GSEA4:

            scores = signal_to_noise_GSEA4(df_A.iloc[0,:].values, df_B.values)
            assoc_matrix = DataFrame(scores, index = df_B.index, columns = df_A.index).T
            
        else:
            assoc_matrix = DataFrame(0, index = df_A.index, columns = df_B.index)
            for a in df_A.index:
                for b in df_B.index: # print('a: {}   b: {}'.format(df_A.loc[a, :], df_B.loc[b, :]))
                    #print('-------------------------------------------------------')
                    print('\na: {}   \nb: {}'.format(df_A.loc[a, :], df_B.loc[b, :]))
                    assoc_matrix.loc[a, b] = assoc_function(df_A.loc[a, :], df_B.loc[b, :])
                    #print('IC= {}'.format(assoc_matrix.loc[a, b]))
                    #print('IC2= {}'.format(information_coefficient(df_A.loc[a, :], df_B.loc[b, :])))
                
    elif axis == 0:
        assoc_matrix = DataFrame(0, index = df_A.columns, columns = df_B.columns)
        #print(assoc_matrix)
        for a in df_A.columns:
            for b in df_B.columns:
                #print('a: {}   b: {}'.format(df_A.loc[:, a], df_B.loc[:, b]))
                x = df_A.loc[:, a]
                y = df_B.loc[:, b]
                assoc_matrix.loc[a, b] = assoc_function(x, y)

    if make_plot == True:
        plot_heat_map(
            assoc_matrix,
            cluster_axis = "01",
            title = title,
            xaxis_title = name_B,
            yaxis_title = name_A,
            html_file_path = file_path)

    return assoc_matrix
