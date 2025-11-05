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

gridsize = [1,3,3]
samp_lb, samp_ub = [-5, -2, -2.5],[-4.5, 0, -0.] 
# poisson average log length is 5.
# default sample bounds are: "Poisson":{'samp_lb':[-8, -3], 'samp_ub':[-5, 0],'gridsize':[6, 7]}}

# Default physical bounds
phys_ub, phys_lb = None, None

# filt_param = {'min_means':[0.01]*3, 'max_maxes':[350, 350, 10000], 'min_maxes':[1,1, 1]}
filt_param = {'min_means':[0]+ [0.01]*2, 'max_maxes':[10000]+[ 350, 10000], 'min_maxes':[0]+ [1,1]}

h5ad_filepath = './combined_adata_tcells_0_with_ambiguous.h5ad'  # Path to save the .h5ad file

# Define model.
fitmodel = cme_toolbox.CMEModel('ProteinBursty','Poisson', fit_unspliced = False,
        min_fudge = 1)


fit_10x_10k = inference.perform_inference(h5ad_filepath, fitmodel, use_lengths=False,
                                         transcriptome_filepath=None,
                                         num_cores=48, n_genes=1, genes_to_fit=['ENSG00000167286'],
                                         filt_param = filt_param, phys_lb=phys_lb, phys_ub=phys_ub,
                                        gradient_params = {'max_iterations':10,'init_pattern':None,'num_restarts':10},
                                        samp_lb=samp_lb,samp_ub=samp_ub,
                                         gridsize=gridsize, exclude_sigma=True, dataset_string='CD3D_Tcell0_10_10')

with open('./CD3D_Tcell0_10_10/ProteinBursty_Poisson_1x3x3/monod_adata.pkl', 'rb') as file:
    fit = pickle.load(file)

analysis.run_qc(fit)#, marg=0)
