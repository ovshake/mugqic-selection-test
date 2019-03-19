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

def read_gene_data(filename):
    """
    @Params: filename: Path to the gene data
    @Returns: df: a DataFrame object with the data
    """
    df = pd.read_csv(filename, sep = '\t') 
    return df 

def process_column(col):
    """
    @Params: col: A column with the gene data of the format x/y
    @Returns: col: Modified column with a numeric representation of x/y
                   with 2*x + y. 
                   ./. is replaced with -1
    """
    for i,x in enumerate(col):
        if  x == BLANK:
            col[i] = -1 
        else:
            activations = list(map(int, x.split('/')))
            col[i] = 2 * activations[0] + activations[1]    
    return col 


def preprocess_data(raw_data, sample_label_dict):
    """
    @Params: raw_data: A DataFrame containing the raw gene data as given in the .vcf file
             sample_label_dict: Dict having mapping between Sample and the label
    @Returns: X: A NumPy array with each row representing a multidimensional data-point.
              y: List containing labels
    """
    # print(raw_data.columns)
    samples = list(raw_data.columns)[9:]  #list(sample_label_df.index.values)  # sample_label_df['sample'].values 
    X , y = [] , [] 
    for sample in samples:
        # print('Cur Sample is ' , sample)
        raw_col_data = raw_data[sample].values
        processed_col_data = process_column(raw_col_data) 
        X.append(processed_col_data) 
        y.append(sample_label_dict[sample])
    
    return np.asarray(X) , y

def plot(X , Y , label_name_map, plot_dir , type_):
    """
    @Params: X: NumPy array containing data
            y: list containing labels
            label_name_map: mapping between abbreviation and full name of the populations
            plot_dir: path of directory storing the plots
            type_: str representing the type of algorithm (PCA / UMAP)
    """
    title = '[{}]{} Dimensionality Reduction Scatterplot'
    dim_red_x = [x[0] for x in X] 
    dim_red_y = [x[1] for x in X] 
    labels = [label_name_map[y] for y in Y] 
    del X 
    X_col_name = '{} - X'
    Y_col_name = '{} - Y'
    data = {X_col_name.format(type_): dim_red_x , Y_col_name.format(type_) : dim_red_y , 'Label': labels}
    df = pd.DataFrame.from_dict(data)
    a4_dims = (14, 14)
    fig,ax = plt.subplots(figsize=a4_dims)
    ax = sns.scatterplot(x = X_col_name.format(type_), y = Y_col_name.format(type_) , data = df , hue = 'Label',  ax = ax)  
    ax.set_title(title.format(ID,type_))
    ax.legend(loc = (0.80, 0))
    plt.savefig(plot_dir + '/' + title.format(ID , type_) + '.jpeg')
    plt.close() 

def write_reduced_coordinates(X , filename):
    """
    @Params: X: NumPy array containing data
            filename: path to the .txt file which will store the 2-D data 
    """
    with open(filename , 'w') as f:
        for x_,y_ in X:
            # print(x_, y_)
            f.write(str(x_) + ',' + str(y_) + '\n') 