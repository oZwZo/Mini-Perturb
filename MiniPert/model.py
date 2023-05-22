import os, sys
import pickle as pkl
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats.distributions import chi2
from itertools import combinations

from multiprocessing import Manager
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning, ValueWarning
warnings.simplefilter('ignore', ValueWarning)

class Perturb_NBGLM(object):
	def __init__(self,
                adata : sc.AnnData,
                perturb_key : str,
                plate_key : str,
                celltype_key : str = None,
                other_covariates : list = None,
                interaction_terms : list = None,
                genes : list = None,
                n_jobs : int = 1,
                alphas : list = None,
                L1_weights : list = None,
					):
		r"""
		Generalized Linear Model following Negative Binomial family for mini-bulk RNA-seq perturbation data.

		Arguments
		-------
		adata : scanpy anndata, required. In the dimension of Sample * Genes.
				The obs dataframe must contain all the covariates needed for model fitting.

		perturb_key : str, the obs column name that indicates the perturbation given to each sample

		plate_key : str, the obs column name stating the Plate/Batch each sample belongs

		celltype_key : str, the obs column name stating the cell type of each sample

		other_covariates : list, default None, a list of obs column name that should be 
							included in the model as independent variables.

		interaction_terms : list of tuple default None. A pair of variables will be added to the model as independent vairiables.
							Each element in the tuple can be the obs name or a category from any design key.
							e.g.1 `interaction_terms=[("TF", "Cell_type")]` if you want to study every interaction pair between TF and Cell type.
							e.g.2 `interaction_terms=[("Ebf1", "HoxB8"), ("Elf1" , "HoxB8")]` if you only want to study the HoxB8 specific effect of Ebf1 and Elf1 perturbation

		genes : list default None. If specified, we only fit models for the given genes.

		n_jobs : int, default 1. The number of processors to use when fitting models.

		alphas : list, default 0.1. The range of regularization strength to search within.
					If alphas = [0], there is no regularization.

		L1_weights : list, default 0. The range of L1_regularization to search. 
						If L1_weights = [0], the model turns into a ridge regression. 
						If L1_weights = [1], the model becomes a LASSO regression.

		"""

		self.n_jobs = n_jobs
		self.alphas = alphas
		self.L1_weights = L1_weights

		if genes is not None:
			assert np.all([g in adata.var_names for g in genes]), "Invalid argument `genes`. Found items that are not in adata"
			self.adata = adata[:,genes].copy()
			self.genes = genes
		else:
			self.genes = adata.var_names
			self.adata = adata.copy()

		self.obs = self.adata.obs
		self.obs_key = list(self.adata.obs_keys())

		# the name of independent variables for glm
		self.independent_variables = []
		self.N_perturb_variables = 0
		self.N_plate_variables = 0
		self.N_interaction_variables = 0

        # a lot of sanity check 
		self.perturb_key = perturb_key
		self.plate_key = plate_key
		self.celltype_key = celltype_key
		self.other_covariates = other_covariates
		self.key_check()

		if interaction_terms is not None:
			self.add_interaction(interaction_terms)

		# construct input df
		Gene_df = pd.DataFrame(self.adata.X, index=self.adata.obs_names,  columns=self.adata.var_names)
		Gene_df['Size_factor'] = Gene_df.sum(axis=1) / 1e4
		self.input_data = Gene_df.merge(self.design_df, left_index=True, right_index=True)

		# GLM formula
		self.formula = 'Gene ~ ' + " + ".join(self.independent_variables)
		print(f"Fitting Negative Binomial regression \n\tfor {len(self.genes)} genes \n\twith {len(self.independent_variables)} independent variables + 1 intercept")
		print("Gene ~ " +  " + ".join(self.independent_variables[:5]) + ' +  ...  + ' + " + ".join(self.independent_variables[-5:]))
		self.pvalue_df = None # initialize

	def add_variable_name(self, key):
		orig_len = len(self.independent_variables)
		self.independent_variables += self.adata.obs[key].unique().tolist()
		assert len(set(self.independent_variables)) == orig_len + self.adata.obs[key].nunique(), "duplicated ind var name detected when adding "

	def add_variable_to_design(self, key):
		series_type = self.obs[key].dtype

		if series_type == 'object':
			df = pd.get_dummies(self.obs[key]).astype(int)
		elif series_type in [int, np.int32, np.in64]:
			if self.obs[key].nunique() == 2:
				df = pd.get_dummies(self.obs[key]).astype(int)
				df = df.rename({col:f"{key}_{col}" for col in df.columns})
			else:
				df = self.obs[[key]]
		else:
			df = self.obs[[key]]

		self.design_df = self.design_df.merge(df, left_index=True, right_index=True)

    
	def key_check(self):

		# required arguments
		assert self.perturb_key in self.obs_key , "Invalid perturb_key"
		assert self.plate_key in self.obs_key , "Invalid plate_key"
		self.perturb_variable = self.adata.obs[self.perturb_key].unique().tolist()
		self.plate_variable = self.adata.obs[self.plate_key].unique().tolist()
		self.add_variable_name(self.perturb_key)
		self.add_variable_name(self.plate_key)

		perturb_matrix = pd.get_dummies(self.obs[self.perturb_key]).astype(int)
		plate_matrix = pd.get_dummies(self.obs[self.plate_key]).astype(int)
		self.design_df = perturb_matrix.merge(plate_matrix, left_index=True, right_index=True)

		# optional arguments
		if self.celltype_key is not None:
			assert self.celltype_key in self.obs_key, "Invalid celltype_key"
			self.add_variable_name(self.celltype_key)
			self.add_variable_to_design(self.celltype_key)

		if self.other_covariates is not None:
			for covariate in self.other_covariates:
				assert covariate in self.obs_key, f"Invalid covariate `{covariate}`"
				self.add_variable_name(covariate)
				self.add_variable_to_design(covariate)

	def add_interaction(self,interaction_list):
		"""
		Add interaction automatically according to the type of the variables
		Input
		------
		interaction_list: list of tuple, [('key1', 'key2')]
		"""
		# for each pair of interaction terms
		for ind_var1, ind_var2 in interaction_list:
			
			if ind_var1 in self.obs_key and self.obs_key.dtype == "object":
				IndV1_ls = self.obs[ind_var1].unique()
			elif ind_var1 in self.independent_variables:
				IndV1_ls = [ind_var1]
			else:
				raise ValueError("invalid interaction term `%s`"%ind_var1)
			
			if ind_var2 in self.obs_key and self.obs_key.dtype == "object":
				IndV2_ls = self.obs[ind_var2].unique()
			elif ind_var2 in self.independent_variables:
				IndV2_ls = [ind_var2]
			else:
				raise ValueError("invalid interaction term `%s`"%ind_var2)
			
			inter_var_names = []
			inter_var_dfs = []
			for IndV1 in IndV1_ls:
				if IndV1 not in self.independent_variables:
					self.independent_variables.append(IndV1)
					self.add_variable_to_design(IndV1)

				for IndV2 in IndV2_ls:
					if IndV2 not in self.independent_variables:
						self.independent_variables.append(IndV2)
						self.add_variable_to_design(IndV2)

					v1_v2 = f"{IndV1} : {IndV2}"
					inter_var_names.append(v1_v2)
					inter_var_dfs.append(self.design_df[IndV1].multiply(self.design_df[IndV2]).to_frame(name=v1_v2))
			
			self.independent_variables += inter_var_names
			IntPairs_df = pd.concat(inter_var_dfs, axis=1)
			self.design_df = self.design_df.merge(IntPairs_df, left_index=True, right_index=True)
					

	def fit_gene(self, gene):
		"""
		fit for each gene
		"""
		formula_gene = self.formula.replace("Gene", gene)
		glm = smf.glm(formula = formula_gene,
						offset = np.log(self.input_data['Size_factor']), 
						data = self.input_data,
						family = sm.families.NegativeBinomial())
		try:
			glm_result = glm.fit()
			beta = glm_result.params.values
		except:
			glm_result = None
			beta = np.zeros((len(self.independent_variables) + 1,))

		# glm_result
		# for alpha in self.alphas:
		# 	for L1_w in self.L1_weights:
		# 		glm_result = glm.fit_regularized(L1_wt=0, alpha=0.1) 

		return glm, glm_result, beta

	def fit(self, n_jobs=None):
		"""
		fit for all genes
		"""
		print("start fittings...")
		if n_jobs is None:
			n_jobs = self.n_jobs
			chunksize = 10
		else:
			chunksize = len(self.genes) // (n_jobs *10)

		res = process_map(self.fit_gene, self.genes, max_workers=n_jobs, chunksize=chunksize)
		self.model_dict = {g:res[i][0] for i,g in enumerate(self.genes)}
		self.result_dict = {g:res[i][1] for i,g in enumerate(self.genes)}
		
		print("\n Done.")
		print("Adding fitted coeffcients to adata.varm with key `Betas`")
		Betas = np.stack([r[2] for r in res])
		self.adata.varm['Betas'] = Betas

		self.Differerntial_analysis([]);
		self.beta_df = pd.DataFrame(Betas.T,  
									index=self.pvalue_df.index, 
									columns=self.pvalue_df.columns)

	def save_to(self, save_path):
		
		checkpoint = {
			"Glms" : self.model_dict,
			"GlmResults" : self.result_dict,
			"Betas" : self.adata.varm['Betas']
		}
		with open(save_path, 'wb') as f:
			pkl.dump(checkpoint, f)
			f.close()
		
		print('MiniPert checkpoint file saved to %s'%save_path)

	def load_from(self, save_path):
		with open(save_path, 'rb') as f:
			checkpoint = pkl.load(f)
			f.close()
		
		print('Checkpoint file loaded from %s'%save_path)
		
		self.model_dict = checkpoint['Glms']
		self.result_dict = checkpoint['GlmResults']
		self.adata.varm['Betas'] = checkpoint['Betas']

		_ = self.Differerntial_analysis(condition=[]) # no return if passing []
		self.beta_df = pd.DataFrame(checkpoint['Betas'].T,  
									index=self.pvalue_df.index, 
									columns=self.pvalue_df.columns)

	def Differerntial_analysis(self, condition=None, threshold=0.05):
		"""
		Perform Gene differential analysis with the coefficients.
		The p-values are obtained by F-test with the null-hypothesis of beta = 0

		Input
		------
		gene: list of str. 
		threshold : float, default 0.05, the threshold used to claim statistical significance

		Return
		-------
		DEGs: dict of list, genes, p-values
		"""
		if self.pvalue_df is None:
			print("\n Testing DEGs of all conditions")
			get_p = lambda g : self.result_dict[g].pvalues.values if self.result_dict[g] is not None else np.ones(1+len(self.independent_variables))
			pmatrix = np.stack([get_p(gene) for gene in self.genes]).T
			self.pvalue_df = pd.DataFrame(pmatrix, 
				columns = self.genes, index = ['intercept'] + self.independent_variables)


		DEGs = {}
		if condition is None:
			condition = self.independent_variables

		for c in condition:
			assert c in self.independent_variables, f"invalid gene {g}"
			idx = np.where(self.pvalue_df.loc[c] < threshold)[0]
			DEGs[c] = [list(self.genes[idx]), self.pvalue_df.loc[c, self.genes[idx]].values]
		return DEGs
	
	def regress_out_plate_effect(self):
        
		plate_betas = self.beta_df.loc[self.plate_variable].values
		plate_design = self.input_data.loc[:,self.plate_variable].values

		plate_effect = plate_design @ plate_betas
		logX = sc.pp.log1p(self.adata.X, copy=True)
		offset = self.input_data['Size_factor'].apply(np.log).values
		correct_X = logX - offset.reshape(-1,1) - plate_effect
		correct_X = np.exp(correct_X)

		self.adata.layers['Plate_Correct_Count'] = correct_X
		return correct_X
