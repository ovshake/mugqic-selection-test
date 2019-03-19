import pandas as pd 
import numpy as np 
import sklearn
from sklearn.decomposition import PCA
import umap
import time 
import seaborn as sns 
import matplotlib.pyplot as plt 
import matplotlib
import logging 
import sys


SEED = 123
ID = 1

BLANK = './.'

def setup_log(filename):
    """
    @Params: filename: path to the log file
    """
    logging.StreamHandler(sys.stdout)
    if filename is None:
        logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S' , level=logging.INFO) 
    else:
        logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S' , filename = filename , level=logging.INFO) 
    
def sample_label_mapping(filename):
    """
    @Params: filename: Path to file which contains Sample Label Mappings 
    @Returns: df: A Dict Object with each sample mapped with the population it belongs 
    """
    df = pd.read_csv(filename , sep = '\t') 
    sample_label_dict = {} 
    for i in range(df.shape[0]):
        sample_label_dict[df.loc[i]['sample']] = df.loc[i]['pop']
    # print(df.head())
    return sample_label_dict 

def label_name_mapping(filename):
    """
    @Params: filename: Path to file which contains Population Abbreviation , Pop. Full name pair
    @Returns: label_name_dict: A Dict Object with each abbreviation mapped with the full name. 
    """
    df = pd.read_csv(filename , sep = '\t')  
    label_name_dict = {} 
    for i in range(df.shape[0]):
        label_name_dict[df.loc[i]['Population Code']] = df.loc[i]['Population Description']
    return label_name_dict 
