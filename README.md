# mugqic-selection-test
This repository is my attempt at the selection test for the project **Human history and data visualization** for **GSoC 2019**. 
## Directory Layout
### src
It contains ```main.py``` and ```utils.py```. ```utils.py``` has all the functions including those which are responsible for applying the PCA and UMAP
algorithm to the dataset. ```main.py``` contains the driver code to call the functions and take parameters using command line arguments. More
information in **How to Run** section.

**Note**: The number in the filename of the plots, the title of the plots and the ```.txt``` files having the 2-d co-ordinates, enclosed in ```[.]``` 
is the unique ID and can be directly mapped with log. 

### plots
This contains the PCA and UMAP plots. 

### reduced_dim_coordinates
This contain the ```.txt``` files which has the 2-D Co-ordinates. 

### data
This contains the data files provided. Since the gene data is very big and my laptop is unable to handle it, I have applied the dataset on a small  sub-part
of the dataset ```data_sample.csv```. I will be uploading the complete results when I manage to run the complete dataset. 

## How to run 
To run the code the libraries required (and the versions which I used are):
* Pandas - 0.23.4
* NumPy - 1.14.2
* sklearn - 0.20.0
* umap - 0.3.5
* Seaborn - 0.9.0
* Matplotlib - 3.0.0

### Steps to run the code
1. clone this repo 
2. run the following command ```python main.py --plot_dir ../plots --crds_dir ../reduced_dim_coordinates --genedata_src ../data/data_sample.csv --label_src ../data/affy_samples.20141118.panel --full_pop_name ../data/20131219.populations.tsv --logfile log.txt```


```main.py``` takes in multiple arguments which are as follows
* ```--logfile```: filepath of the log file, if not given, the default path is stdout 
* ```--plot_dir```: Directory path for the plots (mandatory)
* ```--crds_dir```: Directory path for the reduced dimensional coordinates (mandatory) 
* ```--genedata_src```: path to the gene data .csv file (mandatory)
* ```--label_src```: path to the sample-label .txt file (mandatory)
* ```--full_pop_name```: path to name-verbose name .txt file (mandatory)

