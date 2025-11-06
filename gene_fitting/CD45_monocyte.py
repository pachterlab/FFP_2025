import sys
sys.path.insert(0, '/home/user/monod_source') # Put path to Monod source code here.

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

gridsize = [1,3,3]
# NB -4.5 is too large for unspliced, but since there's only one grid point for unspliced, it doesn't matter.
samp_lb, samp_ub = [-5, -2, -2.5],[-4.5, 0, -0.] 
# poisson average log length is 5.

# Default physical bounds
phys_ub, phys_lb = None, None

filt_param = {'min_means':[0]+ [0.01]*2, 'max_maxes':[10000]+[ 350, 10000], 'min_maxes':[0]+ [1,1]}

h5ad_filepath = '../combined_adata_monocytes_all.h5ad'  # Path to save the .h5ad file

# Define model.
fitmodel = cme_toolbox.CMEModel('ProteinBursty','Poisson', fit_unspliced = False,
        min_fudge = 1)

fit_10x_10k = inference.perform_inference(h5ad_filepath, fitmodel, use_lengths=False,
                                         transcriptome_filepath=None,
                                         num_cores=48, n_genes=1, genes_to_fit=['ENSG00000081237'],
                                         filt_param = filt_param, phys_lb=phys_lb, phys_ub=phys_ub,
                                        gradient_params = {'max_iterations':10,'init_pattern':None,'num_restarts':10},
                                        samp_lb=samp_lb,samp_ub=samp_ub,
                                         gridsize=gridsize, exclude_sigma=True, dataset_string='CD45_moncyte_10_10')

with open('/home/cat/proMonod/2024_10x_v4/individual genes/CD45_moncyte_10_10/ProteinBursty_Poisson_1x3x3/monod_adata.pkl', 'rb') as file:
    fit = pickle.load(file)

analysis.run_qc(fit)#, marg=0)
