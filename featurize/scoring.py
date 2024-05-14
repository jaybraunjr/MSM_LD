from deeptime.decomposition import VAMP
"""

This will score the features with VAMP2 among others using deeptime and pyemma

"""
import numpy as np
import deeptime as dt
import pyemma
from pyemma.util.contexts import settings
import matplotlib.pyplot as plt




def return_vamp2_score(features):
    """
    Returns the VAMP2 score of the features using deeptime

    Parameters
    ----------
    features : ndarray
        The features to be scored

    Returns
    -------
    score : float
        The VAMP2 score of the features
    """
    score = dt.decomposition.vamp_score(features, lag=1, dim=1, scaling="kinetic_map", score_method="VAMP2")
    return score

    


def plot_vamp2_score(return_vamp2_score, dim, lag, scaling, score_method):
    
    """
    Plots the VAMP2 score of the features using deeptime using the input features of multiple datasets

    Parameters
    ----------
    features : ndarray
        The features to be plotted

    Returns
    -------
    bar plot of the VAMP2 score of the features
    """
    plt.bar(features, return_vamp2_score)
    plt.xlabel("Features")
    plt.ylabel("VAMP2 Score")
    plt.title("VAMP2 Score of Features")
    plt.show()
    return  plt.show()


def plot_vamp_score_clustering():
    '''
    This will plot the VAMP score of the features using pyemma. It will 
    determine the optimal amount of k-means clusters to use for the features based on the VAMP-2 score

    '''
    with settings(show_progress_bars=False):
        score = pyemma.coordinates.vamp(features, lag=1, dim=1, scaling="kinetic_map", score_method="VAMP2")
        print(score)
        return score