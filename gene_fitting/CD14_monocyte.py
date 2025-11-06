import sys
sys.path.insert(0, '/home/cat/monod/src/monod')

import extract_data, cme_toolbox, inference, analysis
import importlib
import numpy as np

import pickle

import logging, sys
logging.basicConfig(stream=sys.stdout)
log = logging.getLogger()
log.setLevel(logging.DEBUG)
import warnings
warnings.filterwarnings("ignore") #warning suppression within script is not respected by colab
warnings.simplefilter('ignore')

# Reload modules to get latest changes
importlib.reload(inference)
importlib.reload(extract_data)
importlib.reload(analysis)
importlib.reload(cme_toolbox)

gridsize = [2,2,2]
samp_lb, samp_ub = [-7, -2, -2.5],[-5, 0, -0.] 
# poisson average log length is 5.
# default sample bounds are: "Poisson":{'samp_lb':[-8, -3], 'samp_ub':[-5, 0],'gridsize':[6, 7]}}

# Default physical bounds
phys_ub, phys_lb = None, None

# filt_param = {'min_means':[0.01]*3, 'max_maxes':[350, 350, 10000], 'min_maxes':[1,1, 1]}
filt_param = {'min_means':[0]+ [0.01]*2, 'max_maxes':[10000]+[ 350, 10000], 'min_maxes':[0]+ [1,1]}

h5ad_filepath = '../combined_adata_monocytes_all.h5ad'  # Path to save the .h5ad file

# Define model.
fitmodel = cme_toolbox.CMEModel('ProteinBursty','Poisson', fit_unspliced = False,
        min_fudge = 1)


fit_10x_10k = inference.perform_inference(h5ad_filepath, fitmodel, use_lengths=False,
                                         transcriptome_filepath=None,
                                         num_cores=48, n_genes=1, genes_to_fit=['ENSG00000170458'],
                                         filt_param = filt_param, phys_lb=phys_lb, phys_ub=phys_ub,
                                        gradient_params = {'max_iterations':10,'init_pattern':'moments','num_restarts':1},
                                        samp_lb=samp_lb,samp_ub=samp_ub,
                                         gridsize=gridsize, exclude_sigma=True, dataset_string='CD14_monocyte_1_10_moments')

with open('/home/cat/proMonod/2024_10x_v4/individual genes/CD14_monocyte_1_10_moments/ProteinBursty_Poisson_2x2x2/monod_adata.pkl', 'rb') as file:
    fit = pickle.load(file)

analysis.run_qc(fit)#, marg=0)
