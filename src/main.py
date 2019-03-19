import utils 
import argparse
import logging
parser = argparse.ArgumentParser(description='Applying Dimensionality Reduction Algorithms on Genes')
parser.add_argument('--logfile' , help = 'filepath of the log file' , default=None) 
parser.add_argument('--plot_dir' , help = 'Directory path for the plots')
parser.add_argument('--crds_dir' , help = 'Directory path for the reduced dimensional coordinates') 
parser.add_argument('--genedata_src' , help = 'path to the gene data .csv file')
parser.add_argument('--label_src' , help = 'path to the sample label .txt file')
parser.add_argument('--full_pop_name' , help= 'path to the full population - label mapping')

args = parser.parse_args() 

if __name__ == '__main__':
    utils.setup_log(args.logfile)
    sample_label_df = utils.sample_label_mapping(args.label_src)
    label_name_df = utils.label_name_mapping(args.full_pop_name) 
    gene_data = utils.read_gene_data(args.genedata_src) 
    X , y = utils.preprocess_data(gene_data , sample_label_df) 
    utils.apply_PCA(X, y, label_name_df, args.plot_dir, args.crds_dir)
    utils.apply_umap(X, y, label_name_df, args.plot_dir, args.crds_dir)
    doc_strings = utils.function_docs()
    logging.info('The functions used are')
    logging.info(doc_strings)
    